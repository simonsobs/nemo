"""

This module defines pipelines - sets of tasks in nemo that we sometimes want to do on different inputs
(e.g., real data or simulated data).

"""

import os
import sys
import glob
import shutil
import time
import astropy
import astropy.io.fits as pyfits
import astropy.table as atpy
from astLib import astWCS
import numpy as np
from scipy import ndimage, interpolate
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

#------------------------------------------------------------------------------------------------------------
def filterMapsAndMakeCatalogs(config, rootOutDir = None, useCachedFilters = False, useCachedRMSMap = False,\
                              useCachedFilteredMaps = False, measureFluxes = True, invertMap = False, \
                              verbose = True, writeAreaMask = False, writeFlagMask = False):
    """Runs the map filtering and catalog construction steps according to the given configuration.
    
    Args:
        config (:obj: 'startup.NemoConfig'): Nemo configuration object.
        rootOutDir (str): If None, use the default given by config. Otherwise, use this to override where the
            output filtered maps and catalogs are written.
        useCachedFilters (:obj:`bool`, optional): If True, and previously made filters are found, they will be
            read from disk, rather than re-calculated (used by source injection simulations).
        useCachedRMSMap (:obj:`bool`, optional): If True, use the previously estimated noise map, which has
            been saved to disk. This is only useful for source injection simulations.
        useCachedFilteredMaps (:obj:`bool`, optional): If True, and previously written maps are found, these
            will be read from disk, rather than re-made. This is useful for performing forced photometry using
            external catalogs without having to run nemo again. Otherwise, this should be set to False.
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
    
    # Multi-pass pipeline, enabled with use of filterSets parameter in config file
    if config.filterSets != [] and useCachedFilters == False and useCachedFilteredMaps == False:
        # If we wanted to save results from each step, could set-up filterSet specific diagnostics dir here
        if rootOutDir is None:
            rootOutDir=config.rootOutDir
        for setNum in config.filterSets:
            print(">>> Filter set: %d" % (setNum))
            config.setFilterSet(setNum)
            if setNum == config.filterSets[-1]:
                writeAreaMask=True
                writeFlagMask=True
            config.filterSetOptions[setNum]['catalog']=_filterMapsAndMakeCatalogs(config, verbose = True,
                                                                                  useCachedFilters = False,
                                                                                  useCachedFilteredMaps = False,
                                                                                  writeAreaMask = writeAreaMask,
                                                                                  writeFlagMask = writeFlagMask)

            if config.filterSetOptions[setNum]['addSiphonedFromSets'] is not None:
                toStack=[config.filterSetOptions[setNum]['catalog']]
                for siphonSetNum in config.filterSetOptions[setNum]['addSiphonedFromSets']:
                    toStack.append(config.filterSetOptions[siphonSetNum]['catalog'])
                config.filterSetOptions[setNum]['catalog']=atpy.vstack(toStack)

            if config.rank == 0:
                if config.filterSetOptions[setNum]['saveCatalog'] == True:
                    if 'label' not in config.filterSetOptions[setNum].keys():
                        label="filterSet%d" % (setNum)
                    else:
                        label=config.filterSetOptions[setNum]['label']
                    outFileName=rootOutDir+os.path.sep+label+"_catalog.fits"
                    catalogs.writeCatalog(config.filterSetOptions[setNum]['catalog'], outFileName)
                    catalogs.catalog2DS9(config.filterSetOptions[setNum]['catalog'], outFileName.replace(".fits", ".reg"))

        catalog=config.filterSetOptions[config.filterSets[-1]]['catalog']

    else:
        # Default single pass behaviour (also used by source injection tests)
        catalog=_filterMapsAndMakeCatalogs(config, rootOutDir = rootOutDir, useCachedFilters = useCachedFilters,
                                           useCachedFilteredMaps = useCachedFilteredMaps,
                                           useCachedRMSMap = useCachedRMSMap, measureFluxes = measureFluxes,
                                           invertMap = invertMap, verbose = verbose,
                                           writeAreaMask = writeAreaMask, writeFlagMask = writeFlagMask)

    if verbose == True and config.rank == 0:
        print("... after map filtering and making catalogs: time since start = %.3f sec" % (time.time()-config._timeStarted))

    return catalog

#------------------------------------------------------------------------------------------------------------
def _filterMapsAndMakeCatalogs(config, rootOutDir = None, useCachedFilters = False, useCachedRMSMap = False,\
                               useCachedFilteredMaps = False, measureFluxes = True, invertMap = False, \
                               verbose = True, writeAreaMask = False, writeFlagMask = False):
    """Runs the map filtering and catalog construction steps according to the given configuration.
    
    Args:
        config (:obj: 'startup.NemoConfig'): Nemo configuration object.
        rootOutDir (str): If None, use the default given by config. Otherwise, use this to override where the
            output filtered maps and catalogs are written.
        useCachedFilters (:obj:`bool`, optional): If True, and previously made filters are found, they will be
            read from disk, rather than re-calculated (used by source injection simulations).
        useCachedFilteredMaps (:obj:`bool`, optional): If True, and previously written maps are found, these
            will be read from disk, rather than re-made. This is useful for performing forced photometry using
            external catalogs without having to run nemo again. Otherwise, this should be set to False.
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
            os.makedirs(d, exist_ok = True)
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
    if photFilter is not None and len(config.parDict['mapFilters']) > 1:
        assert(filtersList[0]['label'] == photFilter)
    photFilteredMapDict=None

    # For source injection stuff (yes, this is messy - see below for why we do this)
    undoPixelWindow=True
    if useCachedRMSMap == True:
        undoPixelWindow=False

    # Make filtered maps for each filter and tile
    catalogDict={}
    areaMaskDict=maps.TileDict({}, tileCoordsDict = config.tileCoordsDict)
    flagMaskDict=maps.TileDict({}, tileCoordsDict = config.tileCoordsDict)
    stitchedFilteredMapDict=maps.TileDict({}, tileCoordsDict = config.tileCoordsDict)
    stitchedSNMapDict=maps.TileDict({}, tileCoordsDict = config.tileCoordsDict)
    stitchedRMSMapDict=maps.TileDict({}, tileCoordsDict = config.tileCoordsDict)
    for tileName in config.tileNames:
        if verbose == True: print(">>> [rank = %d] Making filtered maps - tileName = %s " % (config.rank, tileName))
        # Operations that only need to be done once go here
        if 'findAndMaskExtended' in config.parDict.keys():
            maps.makeExtendedSourceMask(config, tileName)
        for f in filtersList:

            label=f['label']+"#"+tileName            
            catalogDict[label]={}
            if 'saveDS9Regions' in f['params'] and f['params']['saveDS9Regions'] == True:
                DS9RegionsPath=config.filteredMapsDir+os.path.sep+tileName+os.path.sep+"%s_filteredMap.reg"  % (label)
            else:
                DS9RegionsPath=None

            filteredMapFileName=filteredMapsDir+os.path.sep+tileName+os.path.sep+"%s_filteredMap.fits"  % (label)
            SNMapFileName=filteredMapsDir+os.path.sep+tileName+os.path.sep+"%s_SNMap.fits" % (label)
            if useCachedFilteredMaps == True and os.path.exists(filteredMapFileName):
                print("... loading cached filtered map %s ..." % (filteredMapFileName))
                filteredMapDict={}
                with pyfits.open(filteredMapFileName) as img:
                    filteredMapDict['data']=img[0].data
                    filteredMapDict['wcs']=astWCS.WCS(img[0].header, mode = 'pyfits')
                    filteredMapDict['mapUnits']=filteredMapDict['wcs'].header['BUNIT']
                    if 'BEAMNSR' in filteredMapDict['wcs'].header.keys():
                        filteredMapDict['beamSolidAngle_nsr']=filteredMapDict['wcs'].header['BEAMNSR']
                        filteredMapDict['obsFreqGHz']=filteredMapDict['wcs'].header['FREQGHZ']
                with pyfits.open(SNMapFileName) as img:
                    filteredMapDict['SNMap']=img[0].data
                filteredMapDict['surveyMask'], wcs=completeness.loadAreaMask(tileName, config.selFnDir)
                filteredMapDict['label']=f['label']
                filteredMapDict['tileName']=tileName
            else:
                # We now keep the reference filter in memory to save some I/O
                if f['label'] == photFilter:
                    returnFilter=True
                else:
                    returnFilter=False
                filterResults=filters.filterMaps(config.unfilteredMapsDictList, f, tileName,
                                                 diagnosticsDir = config.diagnosticsDir, selFnDir = config.selFnDir,
                                                 verbose = True, undoPixelWindow = undoPixelWindow,
                                                 useCachedFilter = useCachedFilters, returnFilter = returnFilter)
                if returnFilter == False:
                    filteredMapDict=filterResults
                else:
                    filteredMapDict, config.cachedFilters[tileName]=filterResults[0], filterResults[1]

            if useCachedRMSMap == True and photFilter is not None: # i.e., only an option for cluster insertion sims
                # This is messy:
                # 1. the saved RMS map doesn't have pixel window correction undone
                # 2. we make the S/N map before doing the pixel window correction (it would cancel)
                # 3. but if we're running a source injection sim using the cached RMS map, we'd make a new SN map
                #    here that has the pixel window correction undone for S, but not for N
                # 4. so, if we're doing source injection sims, we DON'T want to undo pixel window in call to filterMaps
                #    above
                # 5. but then we DO want to undo the pixel window after we've made our new S/N map here
                RMSMap, wcs=completeness.loadRMSMap(tileName, config.selFnDir, photFilter)
                validMask=np.greater(RMSMap, 0)
                SNMap=np.zeros(filteredMapDict['data'].shape)+filteredMapDict['data']
                SNMap[validMask]=SNMap[validMask]/RMSMap[validMask]
                filteredMapDict['SNMap']=SNMap
                mask=np.equal(filteredMapDict['data'], 0)
                filteredMapDict['data']=enmap.apply_window(filteredMapDict['data'], pow=-1.0)
                filteredMapDict['data'][mask]=0 # just in case we rely elsewhere on zero == no data

            # New behavior - only save stitched tiles [and make stitched maps even if not tiled]
            if 'saveFilteredMaps' in f['params'] and f['params']['saveFilteredMaps'] == True:
                stitchedFilteredMapDict[tileName]=filteredMapDict['data'].astype(np.float32)
                stitchedSNMapDict[tileName]=filteredMapDict['SNMap'].astype(np.float32)
                stitchedRMSMapDict[tileName]=filteredMapDict['RMSMap'].astype(np.float32)
                # If needed for debugging new behaviour
                # maps.saveFITS(filteredMapFileName, filteredMapDict['data'], filteredMapDict['wcs'])
                # maps.saveFITS(SNMapFileName, filteredMapDict['SNMap'], filteredMapDict['wcs'])

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
            
            # We collect the area mask here, because it gets modified by findObjects if removing rings
            # NOTE: area mask should always be the same across all filters, could add a consistency check here
            if writeAreaMask == True and tileName not in areaMaskDict.keys():
                areaMaskDict[tileName]=np.array(filteredMapDict['surveyMask'], dtype = np.uint8)

            if writeFlagMask == True and tileName not in flagMaskDict.keys():
                flagMaskDict[tileName]=filteredMapDict['flagMask']
            
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
            del filteredMapDict
            catalogDict[label]['catalog']=catalog

    # Merged/optimal catalogs
    optimalCatalog=catalogs.makeOptimalCatalog(catalogDict, constraintsList = config.parDict['catalogCuts'])
    
    # Gathering catalogs
    if config.MPIEnabled == True:
        # Every process needs the whole catalog, for running in multipass mode
        optimalCatalogList=config.comm.allgather(optimalCatalog)
        if config.rank == 0: print("... gathered catalogs")
        toStack=[]  # We sometimes return [] if no objects found - we can't vstack those
        for collectedTab in optimalCatalogList:
            if type(collectedTab) == astropy.table.table.Table and len(collectedTab) > 0:
                toStack.append(collectedTab)
        optimalCatalog=atpy.vstack(toStack)
        # Strip out duplicates (this is necessary when run in tileDir mode under MPI)
        if len(optimalCatalog) > 0:
            optimalCatalog, numDuplicatesFound, names=catalogs.removeDuplicates(optimalCatalog)

    # Gathering and stitching tiles
    # NOTE: Order matters, bumped up RMS map to make sure written before next stage
    # It doesn't matter if rank 0 then takes its time with the other maps that aren't needed again
    labelsList=['area mask', 'flag mask', 'RMS map', 'filtered map', 'S/N map']
    tileDictsList=[areaMaskDict, flagMaskDict, stitchedRMSMapDict, stitchedFilteredMapDict, stitchedSNMapDict]
    writeMEFList=[writeAreaMask, writeFlagMask, True, False, False]
    writeStitchedList=[writeAreaMask, writeFlagMask, True, True, True]
    MEFPaths=[config.selFnDir+os.path.sep+"areaMask.fits",
              config.selFnDir+os.path.sep+"flagMask.fits",
              config.selFnDir+os.path.sep+"RMSMap_%s.fits" % (photFilter),
              None, None]
    stitchedPaths=[config.selFnDir+os.path.sep+"stitched_areaMask.fits",
                   config.selFnDir+os.path.sep+"stitched_flagMask.fits",
                   config.selFnDir+os.path.sep+"stitched_RMSMap_%s.fits" % (photFilter),
                   filteredMapsDir+os.path.sep+"stitched_%s_filteredMap.fits" % (photFilter),
                   filteredMapsDir+os.path.sep+"stitched_%s_SNMap.fits" % (photFilter)]
    compressionTypeList=['PLIO_1', 'PLIO_1', None, None, None]
    for i in range(len(tileDictsList)):
        label=labelsList[i]
        tileDict=tileDictsList[i]
        writeMEF=writeMEFList[i]
        writeStitched=writeStitchedList[i]
        MEFPath=MEFPaths[i]
        stitchedPath=stitchedPaths[i]
        compressionType=compressionTypeList[i]
        # MPI stuff
        if (writeMEF == True or writeStitched == True) and config.MPIEnabled == True:
            # Hmm. This isn't reliable on wits-core at least
            # gathered_tileDicts=config.comm.gather(tileDict, root = 0)
            # if config.rank == 0:
            #     print("... gathered %s" % (label))
            #     for rankTileDict in gathered_tileDicts:
            #         for key in rankTileDict:
            #             if key not in tileDict:
            #                 tileDict[key]=rankTileDict[key]
            # else:
            #     del tileDict
            if config.rank > 0:
                config.comm.send(tileDict, dest = 0)
                del tileDict
            elif config.rank == 0:
                gathered_tileDicts=[]
                gathered_tileDicts.append(tileDict)
                for source in range(1, config.size):
                    gathered_tileDicts.append(config.comm.recv(source = source))
                print("... gathered %s" % (label))
                for t in gathered_tileDicts:
                    for tileName in t.keys():
                        tileDict[tileName]=t[tileName]
        # Write MEFs and stitched versions
        if config.rank == 0:
            if writeMEF == True and MEFPath is not None and len(tileDict) > 0:
                tileDict.saveMEF(MEFPath, compressionType = compressionType)
            if writeStitched == True and stitchedPath is not None and len(tileDict) > 0:
                tileDict.saveStitchedFITS(stitchedPath, config.origWCS, compressionType = compressionType)
            del tileDict

    return optimalCatalog

#------------------------------------------------------------------------------------------------------------
def makeRMSTables(config):
    """Makes a collection of selection function dictionaries (one per footprint specified in selFnFootprints
    in the config file, plus the full survey mask), that contain information on noise levels and area covered,
    and completeness. Adds footprint columns to object catalog.

    """

    # We only care about the filter used for fixed_ columns
    if config.parDict['photFilter'] is None:
        return None
    photFilterLabel=config.parDict['photFilter']
    for filterDict in config.parDict['mapFilters']:
        if filterDict['label'] == photFilterLabel:
            break
    
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
        RMSTab=completeness.getRMSTab(tileName, photFilterLabel, config.selFnDir)
        selFnDict={'tileName': tileName,
                   'RMSTab': RMSTab,
                   'tileAreaDeg2': RMSTab['areaDeg2'].sum()}
        selFnCollection['full'].append(selFnDict)

        # Generate footprint intersection masks (e.g., with HSC) and RMS tables, which are cached
        # May as well do this bit here (in parallel) and assemble output later
        for footprintDict in footprintsList:
            completeness.makeIntersectionMask(tileName, config.selFnDir, footprintDict['label'], masksList = footprintDict['maskList'])
            tileAreaDeg2=completeness.getTileTotalAreaDeg2(tileName, config.selFnDir, footprintLabel = footprintDict['label'])
            if tileAreaDeg2 > 0:
                RMSTab=completeness.getRMSTab(tileName, photFilterLabel, config.selFnDir,
                                              footprintLabel = footprintDict['label'])
                selFnDict={'tileName': tileName,
                           'RMSTab': RMSTab,
                           'tileAreaDeg2': RMSTab['areaDeg2'].sum()}
                selFnCollection[footprintDict['label']].append(selFnDict)

    if config.MPIEnabled == True:
        gathered_selFnCollections=config.comm.gather(selFnCollection, root = 0)
        if config.rank == 0:
            print("... gathered RMS tables")
            all_selFnCollection={'full': []}
            for key in selFnCollection.keys():
                if key not in all_selFnCollection.keys():
                    all_selFnCollection[key]=[]
            for selFnCollection in gathered_selFnCollections:
                for key in all_selFnCollection.keys():
                    all_selFnCollection[key]=all_selFnCollection[key]+selFnCollection[key]
            selFnCollection=all_selFnCollection

    # Combine RMSTab files (we can downsample further later if needed)
    # We add a column for the tileName just in case want to select on this later
    if config.rank == 0:
        strLen=0
        for tileName in config.allTileNames:
            if len(tileName) > strLen:
                strLen=len(tileName)
        for footprint in selFnCollection.keys():
            if footprint == "full":
                label=""
            else:
                label="_"+footprint
            outFileName=config.selFnDir+os.path.sep+"RMSTab"+label+".fits"
            tabList=[]
            for selFnDict in selFnCollection[footprint]:
                tileName=selFnDict['tileName']
                tileTab=selFnDict['RMSTab']
                tileTab['tileName']=atpy.Column(np.array([tileName]*len(tileTab),
                                                dtype = '<U%d' % (strLen)), "tileName")
                tabList.append(tileTab)
            if len(tabList) > 0:
                tab=atpy.vstack(tabList)
                tab.sort('y0RMS')
                tab.meta['NEMOVER']=nemo.__version__
                tab.write(outFileName, overwrite = True)

        # Add footprint columns to object catalog
        catFileName=config.rootOutDir+os.path.sep+"%s_optimalCatalog.fits" % (os.path.split(config.rootOutDir)[-1])
        tab=atpy.Table().read(catFileName)
        for footprintDict in footprintsList:
            for maskPath in footprintDict['maskList']:
                m, wcs=maps.chunkLoadMask(maskPath)
                tab=catalogs.addFootprintColumnToCatalog(tab, footprintDict['label'], m, wcs)
        catalogs.writeCatalog(tab, catFileName)
        catalogs.writeCatalog(tab, catFileName.replace(".fits", ".csv"))

#------------------------------------------------------------------------------------------------------------
def makeMockClusterCatalog(config, numMocksToMake = 1, combineMocks = False, writeCatalogs = True,\
                           writeInfo = True, verbose = True, QSource = 'fit'):
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
    
    Q=signals.QFit(QSource = QSource, selFnDir = config.selFnDir, tileNames = config.allTileNames)

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
                                     delta = delta, rhoType = rhoType)
    print("... mock survey parameters:")
    print("    QSource = %s" % (QSource))
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
            mockTab=mockSurvey.drawSample(RMSMapDict[tileName], scalingRelationDict, Q, wcs = wcsDict[tileName], 
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
            tab.meta['QSOURCE']=QSource
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
        tab.meta['QSOURCE']=QSource
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
    kernelDict={}   # keys: tile, obsFreqGHz
    for tileName in config.tileNames:
        if tileName not in kernelDict.keys():
            kernelDict[tileName]={}
        for i in range(1, len(config.unfilteredMapsDictList)):
            mapDict=config.unfilteredMapsDictList[i]
            beam=beams[i]
            degPerPix=np.mean(np.diff(beam.rDeg))
            assert(abs(np.diff(beam.rDeg).max()-degPerPix) < 0.001)
            
            # Calculate convolution kernel
            sizePix=beam.profile1d.shape[0]*2
            if sizePix % 2 == 0:
                sizePix=sizePix+1
            symRDeg=np.linspace(-0.5, 0.5, sizePix)
            assert((symRDeg == 0).sum())
            
            symProf=interpolate.splev(abs(symRDeg), beam.tck)
            symRefProf=interpolate.splev(abs(symRDeg), refBeam.tck)

            fSymRef=np.fft.fft(np.fft.fftshift(symRefProf))
            fSymBeam=np.fft.fft(np.fft.fftshift(symProf))
            fSymConv=fSymRef/fSymBeam
            fSymConv[fSymBeam < 1e-1]=0 # Was 1e-2; this value avoids ringing, smaller values do not
            symMatched=np.fft.ifft(fSymBeam*fSymConv).real
            symConv=np.fft.ifft(fSymConv).real

            # This allows normalization in same way as Gaussian smooth method
            symConv=symConv/symConv.sum()
            convedProf=ndimage.convolve(symProf, np.fft.fftshift(symConv))
            attenuationFactor=1/convedProf.max() # norm

            # Make profile object
            peakIndex=np.argmax(np.fft.fftshift(symConv))       
            convKernel=signals.BeamProfile(profile1d = np.fft.fftshift(symConv)[peakIndex:], rDeg = symRDeg[peakIndex:])
            
            ## Check plots
            #import pylab as plt
            #plt.figure(figsize=(10,8))
            #plt.plot(abs(symRDeg*60), symRefProf, label = 'ref', lw = 3)
            #plt.plot(abs(symRDeg*60), convedProf*attenuationFactor, label = 'kernel convolved')
            #integralRatio=np.trapz(symRefProf)/np.trapz(convedProf*attenuationFactor)
            #plt.title("%.3f" % (integralRatio))
            #plt.semilogy()
            #plt.legend()
            #ratio=(convedProf*attenuationFactor)/symRefProf
            #plt.figure(figsize=(10,8))
            #plt.plot(abs(symRDeg*60), ratio, label = 'ratio')
            #plt.plot(abs(symRDeg*60), [1.0]*len(symRDeg), 'r-')
            #plt.legend()
            
            # Fudging 2d kernel to match (fix properly later)
            # NOTE: Now done at higher res but doesn't make much difference
            # (but DOES blow up in some tiles if you use e.g. have the resolution)
            wcs=astWCS.WCS(config.tileCoordsDict[tileName]['header'], mode = 'pyfits').copy()
            wcs.header['CDELT1']=np.diff(refBeam.rDeg)[0]
            wcs.header['CDELT2']=np.diff(refBeam.rDeg)[0]
            wcs.header['NAXIS1']=int(np.ceil(2*refBeam.rDeg.max()/wcs.header['CDELT1'])) 
            wcs.header['NAXIS2']=int(np.ceil(2*refBeam.rDeg.max()/wcs.header['CDELT2']))
            wcs.updateFromHeader()
            shape=(wcs.header['NAXIS2'], wcs.header['NAXIS1'])
            degreesMap=np.ones([shape[0], shape[1]], dtype = float)*1e6
            RADeg, decDeg=wcs.pix2wcs(int(degreesMap.shape[1]/2), int(degreesMap.shape[0]/2))
            degreesMap, xBounds, yBounds=maps.makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, 1.0)
            beamMap=signals.makeBeamModelSignalMap(degreesMap, wcs, beam, amplitude = None)
            refBeamMap=signals.makeBeamModelSignalMap(degreesMap, wcs, refBeam, amplitude = None)
            matchedBeamMap=maps.convolveMapWithBeam(beamMap*attenuationFactor, wcs, convKernel, maxDistDegrees = 1.0)
            
            # Find and apply radial fudge factor
            yRow=np.where(refBeamMap == refBeamMap.max())[0][0]
            rowValid=np.logical_and(degreesMap[yRow] < refBeam.rDeg.max(), matchedBeamMap[yRow] != 0)
            ratio=refBeamMap[yRow][rowValid]/matchedBeamMap[yRow][rowValid]
            zeroIndex=np.argmin(degreesMap[yRow][rowValid])
            assert(degreesMap[yRow][rowValid][zeroIndex] == 0)
            tck=interpolate.splrep(degreesMap[yRow][rowValid][zeroIndex:], ratio[zeroIndex:])
            fudge=interpolate.splev(convKernel.rDeg, tck)
            #fudge[fudge < 0.5]=1.0
            #fudge[fudge > 1.5]=1.0
            fudgeKernel=signals.BeamProfile(profile1d = convKernel.profile1d*fudge, rDeg = convKernel.rDeg)
            
            ## Check plot
            #import pylab as plt
            #plt.figure(figsize=(10,8))
            #plt.plot(convKernel.rDeg, fudge, lw = 3, label = 'fudge')
            #plt.plot(convKernel.rDeg, [1.0]*len(fudge), 'r-')
            #plt.title("fudge")
            ##plt.ylim(0, 2)
            #plt.legend()
            #plt.show()
            
            # 2nd fudge factor - match integrals of 2d kernels
            fudgeMatchedBeamMap=maps.convolveMapWithBeam(beamMap*attenuationFactor, wcs, fudgeKernel, maxDistDegrees = 1.0)
            attenuationFactor=refBeamMap.sum()/fudgeMatchedBeamMap.sum()
                        
            # Check at map pixelization that is actually used
            #shape=(config.tileCoordsDict[tileName]['header']['NAXIS2'], 
                   #config.tileCoordsDict[tileName]['header']['NAXIS1'])
            #wcs=astWCS.WCS(config.tileCoordsDict[tileName]['header'], mode = 'pyfits').copy()
            #degreesMap=np.ones([shape[0], shape[1]], dtype = float)*1e6
            #RADeg, decDeg=wcs.pix2wcs(int(degreesMap.shape[1]/2), int(degreesMap.shape[0]/2))
            #degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, 1.0)
            #beamMap=signals.makeBeamModelSignalMap(degreesMap, wcs, beam, amplitude = None)
            #refBeamMap=signals.makeBeamModelSignalMap(degreesMap, wcs, refBeam, amplitude = None)
            #fudgeMatchedBeamMap=maps.convolveMapWithBeam(beamMap*attenuationFactor, wcs, fudgeKernel, maxDistDegrees = 1.0)            
            ## Check plot
            #import pylab as plt
            #yRow=np.where(refBeamMap == refBeamMap.max())[0][0]
            #rowValid=np.logical_and(degreesMap[yRow] < refBeam.rDeg.max(), fudgeMatchedBeamMap[yRow] != 0)
            #plt.figure(figsize=(10,8))
            #plt.plot(degreesMap[yRow][rowValid]*60, refBeamMap[yRow][rowValid], lw = 3, label = 'ref')
            #plt.plot(degreesMap[yRow][rowValid]*60, fudgeMatchedBeamMap[yRow][rowValid], label = 'fudged')
            #integralRatio=np.trapz(fudgeMatchedBeamMap[yRow][rowValid])/np.trapz(refBeamMap[yRow][rowValid])
            #plt.title("native map res - %.3f" % (integralRatio))
            #plt.semilogy()
            #plt.ylim(1e-5)
            #plt.legend()
            #plt.show()
            #from astLib import astImages
            #astImages.saveFITS("ref.fits", refBeamMap, wcs)
            #astImages.saveFITS("fudgematched.fits", fudgeMatchedBeamMap, wcs)
            #astImages.saveFITS("diff.fits", refBeamMap-fudgeMatchedBeamMap, wcs)
            #import IPython
            #IPython.embed()
            #sys.exit()
            
            # NOTE: If we're NOT passing in 2d kernels, don't need to organise by tile
            kernelDict[tileName][mapDict['obsFreqGHz']]={'smoothKernel': fudgeKernel, 
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
                    mapDict['smoothKernel']=kernelDict[tileName][mapDict['obsFreqGHz']]['smoothKernel']
                    mapDict['smoothAttenuationFactor']=kernelDict[tileName][mapDict['obsFreqGHz']]['smoothAttenuationFactor']
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
            degreesMap, xBounds, yBounds=maps.makeDegreesDistanceMap(degreesMap, wcs,
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
                degreesMap, xBounds, yBounds=maps.makeDegreesDistanceMap(degreesMap, wcs,
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
