"""

This module defines pipelines - sets of tasks in nemo that we sometimes want to do on different inputs
(e.g., real data or simulated data).

"""

import os
import sys
import glob
import shutil
import time
import astropy.io.fits as pyfits
import astropy.table as atpy
from astLib import astWCS
import numpy as np
from . import startUp
from . import filters
from . import photometry
from . import catalogs
from . import maps
from . import signals
from . import completeness
from . import MockSurvey
#import IPython

#------------------------------------------------------------------------------------------------------------
def filterMapsAndMakeCatalogs(config, rootOutDir = None, copyFilters = False, measureFluxes = True, 
                              invertMap = False, verbose = True, useCachedMaps = True):
    """Runs the map filtering and catalog construction steps according to the given configuration.
    
    Args:
        config (:obj: 'startup.NemoConfig'): Nemo configuration object.
        rootOutDir (str): If None, use the default given by config. Otherwise, use this to override where the
            output filtered maps and catalogs are written.
        copyFilters (bool, optional): If True, and rootOutDir is given (not None), then filters will be
            copied from the default output location (from a pre-existing nemo run) to the appropriate
            directory under rootOutDir. This is used by, e.g., contamination tests based on sky sims, where
            the same kernels as used on the real data are applied to simulated maps. If rootOutDir = None,
            setting copyKernels = True has no effect.
        measureFluxes (bool, optional): If True, measure fluxes. If False, just extract S/N values for 
            detected objects.
        invertMap (bool, optional): If True, multiply all maps by -1; needed by 
            :meth:maps.estimateContaminationFromInvertedMaps).
    
    Returns:
        Optimal catalog (keeps the highest S/N detection when filtering at multiple scales).
    
    Note:
        See bin/nemo for how this pipeline is applied to real data, and maps.sourceInjectionTest
        for how this is applied to source-free sims that are generated on the fly.
        
    """
    
    if config.parDict['twoPass'] == False:
        catalog=_filterMapsAndMakeCatalogs(config, rootOutDir = None, copyFilters = False, 
                                           measureFluxes = True, invertMap = False, verbose = True, 
                                           useCachedMaps = True)
    
    else:
        
        # Two pass pipeline
        # On 1st pass, find sources (and maybe clusters) with canned settings, masking nothing.
        # On 2nd pass, the 1st pass catalog will be used to mask or subtract sources from maps used for 
        # noise estimation only.
        
        # Sanity checks first
        # No point doing this if point source masks or catalogs are used
        if 'maskPointSourcesFromCatalog' in config.parDict.keys():
            raise Exception("There is no point running in two-pass mode if maskPointSourcesFromCatalog is set.")
        # No point doing this if we're not using the map itself for the noise term in the filter
        for f in config.parDict['mapFilters']:
            for key in f.keys():
                if key == 'noiseParams' and f['noiseParams']['method'] != 'dataMap':
                    raise Exception("There is no point running if filter noise method != 'dataMap'.")

        # Pass 1 - find point sources, save nothing
        # NOTE: We need to do this for each map in the list, if we have a multi-frequency filter
        pass1PtSrcSettings={'label': "Beam",
                            'class': "BeamMatchedFilter",
                            'params': {'noiseParams': {'method': "model",
                                                       'noiseGridArcmin': 40.0,
                                                       'numNoiseBins': 2},
                            'saveFilteredMaps': False,
                            'outputUnits': 'uK',
                            'edgeTrimArcmin': 0.0}}
        config.parDict['mapFilters']=[pass1PtSrcSettings]
        config.parDict['photFilter']=None
        config.parDict['measureShapes']=True    # Double-lobed extended source at f090 causes havoc in one tile
        orig_unfilteredMapsDictList=list(config.unfilteredMapsDictList)
        pass1CatalogsList=[]
        surveyMasksList=[] # ok, these should all be the same, otherwise we have problems...
        for mapDict in orig_unfilteredMapsDictList:
            # We use whole tile area (i.e., don't trim overlaps) so that we get everything if under MPI
            # Otherwise, powerful sources in overlap regions mess things up under MPI
            # Serial mode doesn't have this issue as it can see the whole catalog over all tiles
            # But since we now use full area, we may double subtract ovelap sources when in serial mode
            # So the removeDuplicates call fixes that, and doesn't impact anything else here
            surveyMasksList.append(mapDict['surveyMask'])
            mapDict['surveyMask']=None
            config.unfilteredMapsDictList=[mapDict]
            catalog=_filterMapsAndMakeCatalogs(config, verbose = False, writeAreaMasks = False)
            catalog, numDuplicatesFound, names=catalogs.removeDuplicates(catalog)       
            pass1CatalogsList.append(catalog)

        # Pass 2 - subtract point sources in the maps used for noise term in filter only
        # To avoid ringing in the pass 2, we siphon off the super bright things found in pass 1
        # We subtract those from the maps used in pass 2 - we then need to add them back at the end
        config.restoreConfig()
        siphonSNR=50
        for mapDict, catalog, surveyMask in zip(orig_unfilteredMapsDictList, pass1CatalogsList, surveyMasksList):
            mapDict['noiseMaskCatalog']=catalog[catalog['SNR'] < siphonSNR]
            mapDict['subtractPointSourcesFromCatalog']=[catalog[catalog['SNR'] > siphonSNR]]
            mapDict['maskSubtractedPointSources']=True
            mapDict['surveyMask']=surveyMask
        config.unfilteredMapsDictList=orig_unfilteredMapsDictList
        catalog=_filterMapsAndMakeCatalogs(config, verbose = False)
        
        # Merge back in the bright sources that were subtracted in pass 1
        mergeList=[catalog]
        for pass1Catalog in pass1CatalogsList:
            mergeList.append(pass1Catalog[pass1Catalog['SNR'] > siphonSNR])
        catalog=atpy.vstack(mergeList)
    
    return catalog
    
#------------------------------------------------------------------------------------------------------------
def _filterMapsAndMakeCatalogs(config, rootOutDir = None, copyFilters = False, measureFluxes = True, 
                               invertMap = False, verbose = True, useCachedMaps = True,
                               writeAreaMasks = True):
    """Runs the map filtering and catalog construction steps according to the given configuration.
    
    Args:
        config (:obj: 'startup.NemoConfig'): Nemo configuration object.
        rootOutDir (str): If None, use the default given by config. Otherwise, use this to override where the
            output filtered maps and catalogs are written.
        copyFilters (bool, optional): If True, and rootOutDir is given (not None), then filters will be
            copied from the default output location (from a pre-existing nemo run) to the appropriate
            directory under rootOutDir. This is used by, e.g., contamination tests based on sky sims, where
            the same kernels as used on the real data are applied to simulated maps. If rootOutDir = None,
            setting copyKernels = True has no effect.
        measureFluxes (bool, optional): If True, measure fluxes. If False, just extract S/N values for 
            detected objects.
        invertMap (bool, optional): If True, multiply all maps by -1; needed by 
            :meth:maps.estimateContaminationFromInvertedMaps).
    
    Returns:
        Optimal catalog (keeps the highest S/N detection when filtering at multiple scales).
    
    Note:
        See bin/nemo for how this pipeline is applied to real data, and maps.sourceInjectionTest
        for how this is applied to source-free sims that are generated on the fly.
        
    """

    # If running on sims (source-free or with injected sources), this ensures we use the same kernels for 
    # filtering the sim maps as was used on the real data, by copying kernels to the sims dir. The kernels 
    # will then be loaded automatically when filterMaps is called. Yes, this is a bit clunky...
    if rootOutDir is not None:
        filteredMapsDir=rootOutDir+os.path.sep+"filteredMaps"
        diagnosticsDir=rootOutDir+os.path.sep+"diagnostics"
        dirList=[rootOutDir, filteredMapsDir, diagnosticsDir]
        for d in dirList:
            if os.path.exists(d) == False:
                os.makedirs(d, exist_ok = True)
        if copyFilters == True:
            for tileName in config.tileNames:
                fileNames=glob.glob(config.diagnosticsDir+os.path.sep+tileName+os.path.sep+"filter*#%s*.fits" % (tileName))
                kernelCopyDestDir=diagnosticsDir+os.path.sep+tileName
                if os.path.exists(kernelCopyDestDir) == False:
                    os.makedirs(kernelCopyDestDir, exist_ok = True)
                for f in fileNames:
                    shutil.copyfile(f, kernelCopyDestDir+os.path.sep+os.path.split(f)[-1]) 
    else:
        rootOutDir=config.rootOutDir
        filteredMapsDir=config.filteredMapsDir
        diagnosticsDir=config.diagnosticsDir
    
    # We re-sort the filters list here - in case we have photFilter defined
    photFilter=config.parDict['photFilter']
    filtersList=[]
    if photFilter is not None:
        for f in config.parDict['mapFilters']:
            if f['label'] == photFilter:
                filtersList.append(f)
    for f in config.parDict['mapFilters']:
        if photFilter is not None:
            if f['label'] == photFilter:
                continue
        filtersList.append(f)
    if photFilter is not None:
        assert(filtersList[0]['label'] == photFilter)
    photFilteredMapDict=None
    
    # Make filtered maps for each filter and tile
    catalogDict={}
    for tileName in config.tileNames:
        # Now have per-tile directories (friendlier for Lustre)
        tileFilteredMapsDir=filteredMapsDir+os.path.sep+tileName
        tileDiagnosticsDir=diagnosticsDir+os.path.sep+tileName
        for d in [tileFilteredMapsDir, tileDiagnosticsDir]:
            os.makedirs(d, exist_ok = True)
        if verbose == True: print(">>> Making filtered maps - tileName = %s ..." % (tileName))
        # We could load the unfiltered map only once here?
        # We could also cache 'dataMap' noise as it will always be the same
        for f in filtersList:

            label=f['label']+"#"+tileName            
            catalogDict[label]={}
            if 'saveDS9Regions' in f['params'] and f['params']['saveDS9Regions'] == True:
                DS9RegionsPath=config.filteredMapsDir+os.path.sep+tileName+os.path.sep+"%s_filteredMap.reg"  % (label)
            else:
                DS9RegionsPath=None
                
            filteredMapDict=filters.filterMaps(config.unfilteredMapsDictList, f, tileName, 
                                               filteredMapsDir = tileFilteredMapsDir,
                                               diagnosticsDir = tileDiagnosticsDir, selFnDir = config.selFnDir, 
                                               verbose = True, undoPixelWindow = True,
                                               useCachedMaps = useCachedMaps)
            
            if f['label'] == photFilter:
                photFilteredMapDict={}
                photFilteredMapDict['SNMap']=filteredMapDict['SNMap']
                photFilteredMapDict['data']=filteredMapDict['data']

            # Forced photometry on user-supplied list of objects, or detect sources
            if 'forcedPhotometryCatalog' in config.parDict.keys() and config.parDict['forcedPhotometryCatalog'] is not None:
                catalog=photometry.makeForcedPhotometryCatalog(filteredMapDict, 
                                                               config.parDict['forcedPhotometryCatalog'],
                                                               useInterpolator = config.parDict['useInterpolator'],
                                                               DS9RegionsPath = DS9RegionsPath)
            else:
                # Normal mode
                catalog=photometry.findObjects(filteredMapDict, threshold = config.parDict['thresholdSigma'], 
                                               minObjPix = config.parDict['minObjPix'], 
                                               findCenterOfMass = config.parDict['findCenterOfMass'], 
                                               removeRings = True, 
                                               ringThresholdSigma = config.parDict['ringThresholdSigma'],
                                               rejectBorder = config.parDict['rejectBorder'], 
                                               objIdent = config.parDict['objIdent'], 
                                               longNames = config.parDict['longNames'],
                                               useInterpolator = config.parDict['useInterpolator'], 
                                               measureShapes = config.parDict['measureShapes'],
                                               invertMap = invertMap,
                                               DS9RegionsPath = DS9RegionsPath)
            
            # We write area mask here, because it gets modified by findObjects if removing rings
            # NOTE: condition added to stop writing tile maps again when running nemoMass in forced photometry mode
            maskFileName=config.selFnDir+os.path.sep+"areaMask#%s.fits" % (tileName)
            surveyMask=np.array(filteredMapDict['surveyMask'], dtype = int)
            if writeAreaMasks == True:
                if os.path.exists(maskFileName) == False and os.path.exists(config.selFnDir+os.path.sep+"areaMask.fits") == False:
                    maps.saveFITS(maskFileName, surveyMask, filteredMapDict['wcs'], compressed = True)
            
            if measureFluxes == True:
                photometry.measureFluxes(catalog, filteredMapDict, config.diagnosticsDir,
                                         photFilteredMapDict = photFilteredMapDict,
                                         useInterpolator = config.parDict['useInterpolator'])

            else:
                # Get S/N only - if the reference (fixed) filter scale has been given
                # This is (probably) only used by maps.estimateContaminationFromInvertedMaps
                if photFilter is not None:
                    photometry.getSNRValues(catalog, photFilteredMapDict['SNMap'], 
                                            filteredMapDict['wcs'], prefix = 'fixed_', 
                                            useInterpolator = config.parDict['useInterpolator'],
                                            invertMap = invertMap)
            
            catalogDict[label]['catalog']=catalog

    # Merged/optimal catalogs
    optimalCatalog=catalogs.makeOptimalCatalog(catalogDict, constraintsList = config.parDict['catalogCuts'])
    
    return optimalCatalog

#------------------------------------------------------------------------------------------------------------
def makeSelFnCollection(config, mockSurvey):
    """Makes a collection of selection function dictionaries (one per footprint specified in selFnFootprints
    in the config file, plus the full survey mask), that contain information on noise levels, area covered, 
    and completeness. 
    
    Returns a dictionary (keys: 'full' - corresponding to whole survey, plus other keys named by footprint).
    
    """
    
    # Q varies across tiles
    tckQFitDict=signals.loadQ(config)
        
    # We only care about the filter used for fixed_ columns
    photFilterLabel=config.parDict['photFilter']
    for filterDict in config.parDict['mapFilters']:
        if filterDict['label'] == photFilterLabel:
            break

    # We'll only calculate completeness for this given selection
    SNRCut=config.parDict['selFnOptions']['fixedSNRCut']

    # Handle any missing options for calcCompleteness (these aren't used by the default fast method anyway)
    if 'numDraws' not in config.parDict['selFnOptions'].keys():
        config.parDict['selFnOptions']['numDraws']=2000000
    if 'numIterations' not in config.parDict['selFnOptions'].keys():
        config.parDict['selFnOptions']['numIterations']=100
    
    # We can calculate stats in different extra areas (e.g., inside optical survey footprints)
    footprintsList=[]
    if 'selFnFootprints' in config.parDict.keys():
        footprintsList=footprintsList+config.parDict['selFnFootprints']
        
    # Run the selection function calculation on each tile in turn
    selFnCollection={'full': []}
    for footprintDict in footprintsList:
        if footprintDict['label'] not in selFnCollection.keys():
            selFnCollection[footprintDict['label']]=[]
            
    for tileName in config.tileNames:
        RMSTab=completeness.getRMSTab(tileName, photFilterLabel, config.selFnDir, diagnosticsDir = config.diagnosticsDir)
        compMz=completeness.calcCompleteness(RMSTab, SNRCut, tileName, mockSurvey, config.parDict['massOptions'], tckQFitDict, 
                                           numDraws = config.parDict['selFnOptions']['numDraws'],
                                           numIterations = config.parDict['selFnOptions']['numIterations'],
                                           method = config.parDict['selFnOptions']['method'])
        selFnDict={'tileName': tileName,
                   'RMSTab': RMSTab,
                   'tileAreaDeg2': RMSTab['areaDeg2'].sum(),
                   'compMz': compMz}
        selFnCollection['full'].append(selFnDict)
        
        # Generate footprint intersection masks (e.g., with HSC) and RMS tables, which are cached
        # May as well do this bit here (in parallel) and assemble output later
        for footprintDict in footprintsList:
            completeness.makeIntersectionMask(tileName, config.selFnDir, footprintDict['label'], masksList = footprintDict['maskList'])
            tileAreaDeg2=completeness.getTileTotalAreaDeg2(tileName, config.selFnDir, footprintLabel = footprintDict['label'])
            if tileAreaDeg2 > 0:
                RMSTab=completeness.getRMSTab(tileName, photFilterLabel, config.selFnDir, diagnosticsDir = config.diagnosticsDir, 
                                              footprintLabel = footprintDict['label'])
                compMz=completeness.calcCompleteness(RMSTab, SNRCut, tileName, mockSurvey, config.parDict['massOptions'], tckQFitDict,
                                                   numDraws = config.parDict['selFnOptions']['numDraws'],
                                                   numIterations = config.parDict['selFnOptions']['numIterations'],
                                                   method = config.parDict['selFnOptions']['method'])
                selFnDict={'tileName': tileName,
                           'RMSTab': RMSTab,
                           'tileAreaDeg2': RMSTab['areaDeg2'].sum(),
                           'compMz': compMz}
                selFnCollection[footprintDict['label']].append(selFnDict)
            
        # Optional mass-limit maps
        if 'massLimitMaps' in list(config.parDict['selFnOptions'].keys()):
            for massLimitDict in config.parDict['selFnOptions']['massLimitMaps']:
                completeness.makeMassLimitMap(SNRCut, massLimitDict['z'], tileName, photFilterLabel, mockSurvey, 
                                            config.parDict['massOptions'], tckQFitDict, config.diagnosticsDir,
                                            config.selFnDir) 
    
    return selFnCollection
                
#------------------------------------------------------------------------------------------------------------
def makeMockClusterCatalog(config, numMocksToMake = 1, combineMocks = False, writeCatalogs = True, 
                           writeInfo = True, verbose = True):
    """Generate a mock cluster catalog using the given nemo config.
    
    Returns:
        List of catalogs (each is an astropy Table object)
    
    """
    
    # Having changed nemoMock interface, we may need to make output dir
    if os.path.exists(config.mocksDir) == False:
        os.makedirs(config.mocksDir, exist_ok = True)
        
    # Noise sources in mocks
    if 'applyPoissonScatter' in config.parDict.keys():
        applyPoissonScatter=config.parDict['applyPoissonScatter']
    else:
        applyPoissonScatter=True
    if 'applyIntrinsicScatter' in config.parDict.keys():
        applyIntrinsicScatter=config.parDict['applyIntrinsicScatter']
    else:
        applyIntrinsicScatter=True
    if 'applyNoiseScatter' in config.parDict.keys():
        applyNoiseScatter=config.parDict['applyNoiseScatter']
    else:
        applyNoiseScatter=True
    if verbose: print(">>> Mock noise sources (Poisson, intrinsic, measurement noise) = (%s, %s, %s) ..." % (applyPoissonScatter, applyIntrinsicScatter, applyNoiseScatter))
    
    # Q varies across tiles
    tckQFitDict=signals.loadQ(config)
    
    # We only care about the filter used for fixed_ columns
    photFilterLabel=config.parDict['photFilter']
    for filterDict in config.parDict['mapFilters']:
        if filterDict['label'] == photFilterLabel:
            break
        
    # The same as was used for detecting objects
    thresholdSigma=config.parDict['thresholdSigma']

    # We need an assumed scaling relation for mock observations
    scalingRelationDict=config.parDict['massOptions']
    
    if verbose: print(">>> Setting up mock survey ...")
    # NOTE: Sanity check is possible here: area in RMSTab should equal area from areaMask.fits
    # If it isn't, there is a problem...
    # Also, we're skipping the individual tile-loading routines here for speed
    checkAreaConsistency=False
    wcsDict={}
    RMSMap=pyfits.open(config.selFnDir+os.path.sep+"RMSMap_%s.fits" % (photFilterLabel))
    RMSTab=atpy.Table().read(config.selFnDir+os.path.sep+"RMSTab.fits")
    count=0
    totalAreaDeg2=0
    RMSMapDict={}
    areaDeg2Dict={}
    if checkAreaConsistency == True:
        areaMap=pyfits.open(config.selFnDir+os.path.sep+"areaMask.fits")
    t0=time.time()
    for tileName in config.tileNames:
        count=count+1
        if tileName == 'PRIMARY':
            if tileName in RMSMap:
                extName=tileName
                data=RMSMap[extName].data
            else:
                data=None
            if data is None:
                for extName in RMSMap:
                    data=RMSMap[extName].data
                    if data is not None:
                        break
            RMSMapDict[tileName]=RMSMap[extName].data
            wcsDict[tileName]=astWCS.WCS(RMSMap[extName].header, mode = 'pyfits')
        else:
            RMSMapDict[tileName]=RMSMap[tileName].data
            wcsDict[tileName]=astWCS.WCS(RMSMap[tileName].header, mode = 'pyfits')
        # Area from RMS table
        areaDeg2=RMSTab[RMSTab['tileName'] == tileName]['areaDeg2'].sum()
        areaDeg2Dict[tileName]=areaDeg2
        totalAreaDeg2=totalAreaDeg2+areaDeg2
        # Area from map (slower)
        if checkAreaConsistency == True:
            areaMask, wcsDict[tileName]=completeness.loadAreaMask(tileName, config.selFnDir)
            areaMask=areaMap[tileName].data
            map_areaDeg2=(areaMask*maps.getPixelAreaArcmin2Map(areaMask, wcsDict[tileName])).sum()/(60**2)
            if abs(map_areaDeg2-areaDeg2) > 1e-4:
                raise Exception("Area from areaMask.fits doesn't agree with area from RMSTab.fits")
    RMSMap.close()
    if checkAreaConsistency == True:
        areaMap.close()
    t1=time.time()
    if verbose: print("... took %.3f sec ..." % (t1-t0))
    
    # Useful for testing:
    if 'seed' in config.parDict.keys():
        seed=config.parDict['seed']
    else:
        seed=None
        
    if seed is not None:
        np.random.seed(seed)
        
    # We're now using one MockSurvey object for the whole survey
    massOptions=config.parDict['massOptions']
    minMass=5e13
    zMin=0.0
    zMax=2.0
    defCosmo={'H0': 70.0, 'Om0': 0.30, 'Ob0': 0.05, 'sigma8': 0.80, 'ns': 0.95, 'delta': 500, 'rhoType': 'critical'}
    for key in defCosmo:
        if key not in massOptions.keys():
            massOptions[key]=defCosmo[key]
    H0=massOptions['H0']
    Om0=massOptions['Om0']
    Ob0=massOptions['Ob0']
    sigma8=massOptions['sigma8']
    ns=massOptions['ns']
    delta=massOptions['delta']
    rhoType=massOptions['rhoType']
    mockSurvey=MockSurvey.MockSurvey(minMass, totalAreaDeg2, zMin, zMax, H0, Om0, Ob0, sigma8, ns, 
                                     delta = delta, rhoType = rhoType, enableDrawSample = True)
    print("... mock survey parameters:")
    for key in defCosmo.keys():
        print("    %s = %s" % (key, str(massOptions[key])))
    for key in ['tenToA0', 'B0', 'Mpivot', 'sigma_int']:
        print("    %s = %s" % (key, str(scalingRelationDict[key])))
    print("    total area = %.1f square degrees" % (totalAreaDeg2))
    print("    random seed = %s" % (str(seed)))
    
    if verbose: print(">>> Making mock catalogs ...")
    catList=[]
    for i in range(numMocksToMake):       
        mockTabsList=[]
        t0=time.time()
        for tileName in config.tileNames:
            # It's possible (depending on tiling) that blank tiles were included - so skip
            # We may also have some tiles that are almost but not quite blank
            if RMSMapDict[tileName].sum() == 0 or areaDeg2Dict[tileName] < 0.5:
                continue
            mockTab=mockSurvey.drawSample(RMSMapDict[tileName], scalingRelationDict, tckQFitDict, wcs = wcsDict[tileName], 
                                          photFilterLabel = photFilterLabel, tileName = tileName, makeNames = True,
                                          SNRLimit = thresholdSigma, applySNRCut = True, 
                                          areaDeg2 = areaDeg2Dict[tileName],
                                          applyPoissonScatter = applyPoissonScatter, 
                                          applyIntrinsicScatter = applyIntrinsicScatter,
                                          applyNoiseScatter = applyNoiseScatter)
            if mockTab is not None:
                mockTabsList.append(mockTab)
        tab=atpy.vstack(mockTabsList)
        catList.append(tab)
        t1=time.time()
        if verbose: print("... making mock catalog %d took %.3f sec ..." % (i+1, t1-t0))
        
        # Write catalog and .reg file
        if writeCatalogs == True:
            #colNames=['name', 'RADeg', 'decDeg', 'template', 'redshift', 'redshiftErr', 'true_M500', 'true_fixed_y_c', 'fixed_SNR', 'fixed_y_c', 'fixed_err_y_c']
            #colFmts =['%s',   '%.6f',  '%.6f',   '%s',       '%.3f',     '%.3f',        '%.3f',      '%.3f',           '%.1f',      '%.3f',      '%.3f']
            mockCatalogFileName=config.mocksDir+os.path.sep+"mockCatalog_%d.csv" % (i+1)
            catalogs.writeCatalog(tab, mockCatalogFileName)
            catalogs.writeCatalog(tab, mockCatalogFileName.replace(".csv", ".fits"))

            addInfo=[{'key': 'fixed_SNR', 'fmt': '%.1f'}]
            catalogs.catalog2DS9(tab, mockCatalogFileName.replace(".csv", ".reg"), constraintsList = [], 
                                    addInfo = addInfo, color = "cyan") 
            
    if combineMocks == True:
        tab=None
        for i in range(numMocksToMake):
            mockCatalogFileName=config.mocksDir+os.path.sep+"mockCatalog_%d.fits" % (i+1)
            stackTab=atpy.Table().read(mockCatalogFileName)
            if tab == None:
                tab=stackTab
            else:
                tab=atpy.vstack([tab, stackTab])
        outFileName=config.mocksDir+os.path.sep+"mockCatalog_combined.fits"
        if os.path.exists(outFileName) == True:
            os.remove(outFileName)
        tab.write(outFileName)
    
    # Write a small text file with the parameters used to generate the mocks into the mocks dir (easier than using headers)
    if writeInfo == True:
        mockKeys=['massOptions', 'makeMockCatalogs', 'applyPoissonScatter', 'applyIntrinsicScatter', 'applyNoiseScatter']
        with open(config.mocksDir+os.path.sep+"mockParameters.txt", "w") as outFile:
            for m in mockKeys:
                if m in config.parDict.keys():
                    outFile.write("%s: %s\n" % (m, config.parDict[m]))
    
    return catList
