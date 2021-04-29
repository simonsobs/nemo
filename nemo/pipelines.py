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
from scipy import ndimage
import copy
from pixell import enmap
import nemo
from . import startUp
from . import filters
from . import photometry
from . import catalogs
from . import maps
from . import signals
from . import completeness
from . import MockSurvey
import nemoCython
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
        catalog=_filterMapsAndMakeCatalogs(config, rootOutDir = rootOutDir, copyFilters = copyFilters, 
                                           measureFluxes = measureFluxes, invertMap = invertMap, 
                                           verbose = verbose, useCachedMaps = useCachedMaps)
    
    else:
        
        # Two pass pipeline
        # On 1st pass, find sources (and maybe clusters) with canned settings, masking nothing.
        # On 2nd pass, the 1st pass catalog will be used to mask or subtract sources from maps used for 
        # noise estimation only.
        
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
        config.parDict['maskPointSourcesFromCatalog']=[]    # This is only applied on the 2nd pass
        config.parDict['measureShapes']=True    # Double-lobed extended source at f090 causes havoc in one tile
        orig_unfilteredMapsDictList=list(config.unfilteredMapsDictList)
        config.parDict['forcedPhotometryCatalog']=None # If in this mode, only wanted on 2nd pass
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
        config.parDict['measureShapes']=True    # We'll keep this for pass 2 as well
        siphonSNR=50
        for mapDict, catalog, surveyMask in zip(orig_unfilteredMapsDictList, pass1CatalogsList, surveyMasksList):
            #catalogs.catalog2DS9(catalog[catalog['SNR'] > siphonSNR], config.diagnosticsDir+os.path.sep+"pass1_highSNR_siphoned.reg")
            mapDict['noiseMaskCatalog']=catalog[catalog['SNR'] < siphonSNR]
            mapDict['subtractPointSourcesFromCatalog']=[catalog[catalog['SNR'] > siphonSNR]]
            mapDict['maskSubtractedPointSources']=True
            mapDict['surveyMask']=surveyMask
        config.unfilteredMapsDictList=orig_unfilteredMapsDictList
        catalog=_filterMapsAndMakeCatalogs(config, verbose = False)
        
        # Merge back in the bright sources that were subtracted in pass 1
        # (but we don't do that in forced photometry mode)
        mergeList=[catalog]
        if config.parDict['forcedPhotometryCatalog'] is None:
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
                                               removeRings = config.parDict['removeRings'],
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
                    maps.saveFITS(maskFileName, surveyMask, filteredMapDict['wcs'], compressed = True,
                                  compressionType = 'PLIO_1')
            
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
            map_areaDeg2=(areaMask*maps.getPixelAreaArcmin2Map(areaMask.shape, wcsDict[tileName])).sum()/(60**2)
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
        tab.meta['NEMOVER']=nemo.__version__
        tab.write(outFileName, overwrite = True)
    
    # Write a small text file with the parameters used to generate the mocks into the mocks dir (easier than using headers)
    if writeInfo == True:
        mockKeys=['massOptions', 'makeMockCatalogs', 'applyPoissonScatter', 'applyIntrinsicScatter', 'applyNoiseScatter']
        with open(config.mocksDir+os.path.sep+"mockParameters.txt", "w") as outFile:
            for m in mockKeys:
                if m in config.parDict.keys():
                    outFile.write("%s: %s\n" % (m, config.parDict[m]))
    
    return catList

#------------------------------------------------------------------------------------------------------------
def extractSpec(config, tab, method = 'CAP', diskRadiusArcmin = 4.0, highPassFilter = False,
                estimateErrors = True, saveFilteredMaps = False):
    """Returns a table containing the spectral energy distribution, extracted using either compensated
    aperture photometry (CAP) at each object location in the input catalog, or using a matched filter.
    Maps at different frequencies will first be matched to the lowest resolution beam, using a Gaussian
    kernel.
    
    For the CAP method, at each object location, the temperature fluctuation is measured within a disk of
    radius diskRadiusArcmin, after subtracting the background measured in an annulus between
    diskRadiusArcmin < r < sqrt(2) * diskRadiusArcmin (i.e., this should be similar to the method
    described in Schaan et al. 2020).
    
    For the matched filter method, the catalog must contain a `template` column, as produced by the main
    `nemo` script, with template names in the format Arnaud_M2e14_z0p4 (for example). This will be used to
    set the signal scale used for each object. All definitions of filters in the config will be ignored,
    in favour of a filter using a simple CMB + white noise model. Identical filters will be used for all
    maps (i.e., the method of Saro et al. 2014).
    
    Args:
        config (:obj:`startup.NemoConfig`): Nemo configuration object.
        tab (:obj:`astropy.table.Table`): Catalog containing input object positions. Must contain columns
            'name', 'RADeg', 'decDeg'.
        method (str, optional):
        diskRadiusArcmin (float, optional): If using CAP method: disk aperture radius in arcmin, within
            which the signal is measured. The background will be estimated in an annulus between
            diskRadiusArcmin < r < sqrt(2) * diskRadiusArcmin.
        highPassFilter (bool, optional): If using CAP method: if set, subtract the large scale
            background using maps.subtractBackground, with the smoothing scale set to 
            2 * sqrt(2) * diskRadiusArcmin.
        estimateErrors (bool, optional): If used CAP method: if set, estimate uncertainties by placing
            random apertures throughout the map. For now, this is done on a tile-by-tile basis, and
            doesn't take into account inhomogeneous noise within a tile.
        saveFilteredMaps (bool, optional): If using matchedFilter method: save the filtered maps under
            the `nemoSpecCache` directory (which is created in the current working directory, if it
            doesn't already exist).
    
    Returns:
        Catalog containing spectral energy distribution measurements for each object.
        For the CAP method, units of extracted signals are uK arcmin^2.
        For the matchedFilter method, extracted signals are deltaT CMB amplitude in uK.
            
    """
    
    diagnosticsDir=config.diagnosticsDir
        
    # Choose lowest resolution as the reference beam - we match to that
    refBeam=None
    refFWHMArcmin=0
    refIndex=0
    beams=[]
    for i in range(len(config.unfilteredMapsDictList)):
        mapDict=config.unfilteredMapsDictList[i]
        beam=signals.BeamProfile(mapDict['beamFileName'])
        if beam.FWHMArcmin > refFWHMArcmin:
            refBeam=beam
            refFWHMArcmin=beam.FWHMArcmin
            refIndex=i
        beams.append(beam)

    # Sort the list of beams and maps so that the one with the reference beam is in index 0
    config.unfilteredMapsDictList.insert(0, config.unfilteredMapsDictList.pop(refIndex))
    beams.insert(0, beams.pop(refIndex))
    
    # Figure out how much we need to Gaussian blur to match the reference beam
    # NOTE: This was an alternative to proper PSF-matching that wasn't good enough for ACT beams
    #for i in range(1, len(config.unfilteredMapsDictList)):
        #mapDict=config.unfilteredMapsDictList[i]
        #beam=beams[i]
        #degPerPix=np.mean(np.diff(beam.rDeg))
        #assert(abs(np.diff(beam.rDeg).max()-degPerPix) < 0.001)
        #resMin=1e6
        #smoothPix=0
        #attFactor=1.0
        #for j in range(1, 100):
            #smoothProf=ndimage.gaussian_filter1d(beam.profile1d, j)
            #smoothProf=smoothProf/smoothProf.max()
            #res=np.sum(np.power(refBeam.profile1d-smoothProf, 2))
            #if res < resMin:
                #resMin=res
                #smoothPix=j
                #attFactor=1/smoothProf.max()
        #smoothScaleDeg=smoothPix*degPerPix
        #mapDict['smoothScaleDeg']=smoothScaleDeg
        #mapDict['smoothAttenuationFactor']=1/ndimage.gaussian_filter1d(beam.profile1d, smoothPix).max()
    
    # For testing on CMB maps here
    refMapDict=config.unfilteredMapsDictList[0]
            
    # PSF matching via a convolution kernel
    kernelDict={}   # keys: obsFreqGHz
    for i in range(1, len(config.unfilteredMapsDictList)):
        mapDict=config.unfilteredMapsDictList[i]
        beam=beams[i]
        degPerPix=np.mean(np.diff(beam.rDeg))
        assert(abs(np.diff(beam.rDeg).max()-degPerPix) < 0.001)
        
        # Calculate convolution kernel
        symRefProf=np.zeros((refBeam.profile1d.shape[0]*2)-1)
        symRefProf[:refBeam.profile1d.shape[0]]=refBeam.profile1d[::-1]
        symRefProf[refBeam.profile1d.shape[0]-1:]=refBeam.profile1d

        symProf=np.zeros((beam.profile1d.shape[0]*2)-1)
        symProf[:beam.profile1d.shape[0]]=beam.profile1d[::-1]
        symProf[beam.profile1d.shape[0]-1:]=beam.profile1d        
        
        symRefProf=np.fft.fftshift(symRefProf)
        symProf=np.fft.fftshift(symProf)
        fSymRef=np.fft.fft(symRefProf)
        fSymBeam=np.fft.fft(symProf)
        fSymConv=fSymRef/fSymBeam
        fSymConv[fSymBeam < 1e-1]=0 # Was 1e-2; this value avoids ringing, smaller values do not
        symMatched=np.fft.ifft(fSymBeam*fSymConv).real
        symConv=np.fft.fftshift(np.fft.ifft(fSymConv).real)

        # This allows normalization in same way as Gaussian smooth method
        symConv=symConv/symConv.sum()
        convedProf=ndimage.convolve(symProf, np.fft.fftshift(symConv))[beam.profile1d.shape[0]-1:]
        attenuationFactor=1/convedProf.max() # norm

        conv=symConv[beam.profile1d.shape[0]-1:]
        convKernel=copy.deepcopy(refBeam)
        convKernel.profile1d=conv
                        
        # Additional numerical fudge factor, calculated from comparing the kernel+beam in 2d with ref beam
        # NOTE: We could properly fix 1d -> 2d projection on chunky pixels - this should be good enough for now
        tileName=config.tileNames[0]
        shape=(config.tileCoordsDict[tileName]['header']['NAXIS2'], 
               config.tileCoordsDict[tileName]['header']['NAXIS1'])
        wcs=astWCS.WCS(config.tileCoordsDict[tileName]['header'], mode = 'pyfits')
        #with pyfits.open(mapDict['mapFileName']) as img:
            #data=img[0].data
            #wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
            #shape=data.shape
        degreesMap=np.ones([shape[0], shape[1]], dtype = float)*1e6
        RADeg, decDeg=wcs.pix2wcs(degreesMap.shape[1]/2, degreesMap.shape[0]/2)
        degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, 1.0)
        beamMap=signals.makeBeamModelSignalMap(degreesMap, wcs, beam, amplitude = None)
        refBeamMap=signals.makeBeamModelSignalMap(degreesMap, wcs, refBeam, amplitude = None)
        matchedBeamMap=maps.convolveMapWithBeam(beamMap*attenuationFactor, wcs, convKernel, maxDistDegrees = 1.0)
        
        refBeamMap=refBeamMap/refBeamMap.max()
        matchedBeamMap=matchedBeamMap/matchedBeamMap.max() 
        fudge=matchedBeamMap.sum()/refBeamMap.sum()
        attenuationFactor=attenuationFactor*fudge
        
        # NOTE: If we're NOT passing in 2d kernels, don't need to organise by tile
        kernelDict[mapDict['obsFreqGHz']]={'smoothKernel': convKernel, 
                                            'smoothAttenuationFactor': attenuationFactor}
        
    if method == 'CAP':
        catalog=_extractSpecCAP(config, tab, kernelDict, diskRadiusArcmin = 4.0, highPassFilter = False,
                                estimateErrors = True)
    elif method == 'matchedFilter':
        catalog=_extractSpecMatchedFilter(config, tab, kernelDict, saveFilteredMaps = saveFilteredMaps)
    else:
        raise Exception("'method' should be 'CAP' or 'matchedFilter'")
    
    return catalog

#------------------------------------------------------------------------------------------------------------
def _extractSpecMatchedFilter(config, tab, kernelDict, saveFilteredMaps = False, noiseMethod = 'dataMap'):
    """See extractSpec.
    
    """
        
    cacheDir="nemoSpecCache"+os.path.sep+os.path.basename(config.rootOutDir)
    os.makedirs(cacheDir, exist_ok = True)

    # Build filter configs
    allFilters={'class': 'ArnaudModelMatchedFilter',
                'params': {'noiseParams': {'method': noiseMethod, 'noiseGridArcmin': 40.0},
                           'saveFilteredMaps': False,
                           'saveRMSMap': False,
                           'savePlots': False,
                           'saveDS9Regions': False,
                           'saveFilter': False,
                           'outputUnits': 'yc',
                           'edgeTrimArcmin': 0.0,
                           'GNFWParams': 'default'}}
    
    filtersList=[]
    templatesUsed=np.unique(tab['template']).tolist()
    for t in templatesUsed:
        newDict=copy.deepcopy(allFilters)
        M500MSun=float(t.split("_M")[-1].split("_")[0])
        z=float(t.split("_z")[-1].replace("p", "."))
        newDict['params']['M500MSun']=M500MSun
        newDict['params']['z']=z
        newDict['label']=t
        filtersList.append(newDict)
    
    # Filter and extract
    # NOTE: We assume index 0 of the unfiltered maps list is the reference for which the filter is made
    catalogList=[]
    for tileName in config.tileNames:
        print("... rank %d: tileName = %s ..." % (config.rank, tileName))
        diagnosticsDir=cacheDir+os.path.sep+tileName
        os.makedirs(diagnosticsDir, exist_ok = True)
        for f in filtersList:
            tempTileTab=None # catalogs are organised by tile and template
            filterObj=None
            for mapDict in config.unfilteredMapsDictList:
                if tempTileTab is None:
                    shape=(config.tileCoordsDict[tileName]['header']['NAXIS2'], 
                           config.tileCoordsDict[tileName]['header']['NAXIS1'])
                    wcs=astWCS.WCS(config.tileCoordsDict[tileName]['header'], mode = 'pyfits')
                    tempTileTab=catalogs.getCatalogWithinImage(tab, shape, wcs)
                    tempTileTab=tempTileTab[tempTileTab['template'] == f['label']]
                if tempTileTab is None or len(tempTileTab) == 0:
                    continue
                if mapDict['obsFreqGHz'] == config.unfilteredMapsDictList[0]['obsFreqGHz']:
                    filteredMapDict, filterObj=filters.filterMaps([mapDict], f, tileName, 
                                                                  filteredMapsDir = cacheDir,
                                                                  diagnosticsDir = diagnosticsDir,
                                                                  selFnDir = cacheDir, 
                                                                  verbose = True, 
                                                                  undoPixelWindow = True,
                                                                  returnFilter = True)
                else:
                    mapDict['smoothKernel']=kernelDict[mapDict['obsFreqGHz']]['smoothKernel']
                    mapDict['smoothAttenuationFactor']=kernelDict[mapDict['obsFreqGHz']]['smoothAttenuationFactor']
                    mapDictToFilter=maps.preprocessMapDict(mapDict.copy(), tileName = tileName)
                    filteredMapDict['data']=filterObj.applyFilter(mapDictToFilter['data'])
                    RMSMap=filterObj.makeNoiseMap(filteredMapDict['data'])
                    filteredMapDict['SNMap']=np.zeros(filterObj.shape)
                    mask=np.greater(filteredMapDict['surveyMask'], 0)
                    filteredMapDict['SNMap'][mask]=filteredMapDict['data'][mask]/RMSMap[mask]
                    filteredMapDict['data']=enmap.apply_window(filteredMapDict['data'], pow=-1.0)
                if saveFilteredMaps == True:
                    outFileName=cacheDir+os.path.sep+'%d_' % (mapDict['obsFreqGHz'])+f['label']+'#'+tileName+'.fits'
                    # Add conversion to delta T in here?
                    maps.saveFITS(outFileName, filteredMapDict['data'], filteredMapDict['wcs'])
                freqTileTab=photometry.makeForcedPhotometryCatalog(filteredMapDict, 
                                                                   tempTileTab,
                                                                   useInterpolator = config.parDict['useInterpolator'])
                photometry.measureFluxes(freqTileTab, filteredMapDict, cacheDir,
                                         useInterpolator = config.parDict['useInterpolator'],
                                         ycObsFreqGHz = mapDict['obsFreqGHz'])
                # We don't take tileName from the catalog, some objects in overlap areas may only get cut here
                if len(freqTileTab) == 0:
                    tempTileTab=None
                    continue
                tempTileTab, freqTileTab, rDeg=catalogs.crossMatch(tempTileTab, freqTileTab, radiusArcmin = 2.5)
                colNames=['deltaT_c', 'y_c', 'SNR']
                suff='_%d' % (mapDict['obsFreqGHz'])
                for colName in colNames:
                    tempTileTab[colName+suff]=freqTileTab[colName]
                    if 'err_'+colName in freqTileTab.keys():
                        tempTileTab['err_'+colName+suff]=freqTileTab['err_'+colName]
            if tempTileTab is not None and len(tempTileTab) > 0:
                catalogList.append(tempTileTab)
    
    if len(catalogList) > 0:
        catalog=atpy.vstack(catalogList)
    else:
        catalog=[]
    
    return catalog
    
#------------------------------------------------------------------------------------------------------------
def _extractSpecCAP(config, tab, kernelDict, method = 'CAP', diskRadiusArcmin = 4.0, highPassFilter = False, 
                    estimateErrors = True):
    """See extractSpec.
    
    """
        
    # Define apertures like Schaan et al. style compensated aperture photometry filter
    innerRadiusArcmin=diskRadiusArcmin
    outerRadiusArcmin=diskRadiusArcmin*np.sqrt(2)
    
    catalogList=[]
    for tileName in config.tileNames:
        
        # This loads the maps, applies any masks, and smooths to approx. same scale
        mapDictList=[]
        freqLabels=[]
        for mapDict in config.unfilteredMapsDictList:           
            mapDict=maps.preprocessMapDict(mapDict.copy(), tileName = tileName)
            if highPassFilter == True:
                mapDict['data']=maps.subtractBackground(mapDict['data'], mapDict['wcs'], 
                                                        smoothScaleDeg = (2*outerRadiusArcmin)/60)
            freqLabels.append(int(round(mapDict['obsFreqGHz'])))
            mapDictList.append(mapDict)
        wcs=mapDict['wcs']
        shape=mapDict['data'].shape
                
        # Extract spectra
        pixAreaMap=maps.getPixelAreaArcmin2Map(shape, wcs)
        maxSizeDeg=(outerRadiusArcmin*1.2)/60
        tileTab=catalogs.getCatalogWithinImage(tab, shape, wcs)
        for label in freqLabels:
            tileTab['diskT_uKArcmin2_%s' % (label)]=np.zeros(len(tileTab))
            tileTab['err_diskT_uKArcmin2_%s' % (label)]=np.zeros(len(tileTab))
            tileTab['diskSNR_%s' % (label)]=np.zeros(len(tileTab))
        for row in tileTab:
            degreesMap=np.ones(shape, dtype = float)*1e6 # NOTE: never move this
            degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, 
                                                                           row['RADeg'], row['decDeg'],
                                                                           maxSizeDeg)
            innerMask=degreesMap < innerRadiusArcmin/60
            outerMask=np.logical_and(degreesMap >= innerRadiusArcmin/60, degreesMap < outerRadiusArcmin/60)
            for mapDict, label in zip(mapDictList, freqLabels):
                d=mapDict['data'] 
                diskFlux=(d[innerMask]*pixAreaMap[innerMask]).sum()-(d[outerMask]*pixAreaMap[outerMask]).sum()
                row['diskT_uKArcmin2_%s' % (label)]=diskFlux
            
        # Estimate noise in every measurement (on average) from spatting down on random positions
        # This will break if noise is inhomogeneous though. But at least it's done separately for each tile.
        # We can later add something that scales / fits using the weight map?
        if estimateErrors == True:
            randTab=catalogs.generateRandomSourcesCatalog(mapDict['surveyMask'], wcs, 1000)
            for label in freqLabels:
                randTab['diskT_uKArcmin2_%s' % (label)]=np.zeros(len(randTab))
            for row in randTab:
                degreesMap=np.ones(shape, dtype = float)*1e6 # NOTE: never move this
                degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, 
                                                                                row['RADeg'], row['decDeg'], 
                                                                                maxSizeDeg)
                innerMask=degreesMap < innerRadiusArcmin/60
                outerMask=np.logical_and(degreesMap >= innerRadiusArcmin/60, degreesMap < outerRadiusArcmin/60)
                for mapDict, label in zip(mapDictList, freqLabels):
                    d=mapDict['data'] 
                    diskFlux=(d[innerMask]*pixAreaMap[innerMask]).sum()-(d[outerMask]*pixAreaMap[outerMask]).sum()
                    row['diskT_uKArcmin2_%s' % (label)]=diskFlux
            noiseLevels={}
            for label in freqLabels:
                if signals.fSZ(float(label)) < 0:
                    SNRSign=-1
                else:
                    SNRSign=1
                noiseLevels[label]=np.percentile(abs(randTab['diskT_uKArcmin2_%s' % (label)]), 68.3)
                tileTab['err_diskT_uKArcmin2_%s' % (label)]=noiseLevels[label]
                tileTab['diskSNR_%s' % (label)]=SNRSign*(tileTab['diskT_uKArcmin2_%s' % (label)]/noiseLevels[label])
                
        catalogList.append(tileTab)
    
    catalog=atpy.vstack(catalogList)
        
    return catalog
