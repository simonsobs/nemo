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
from . import startUp
from . import filters
from . import photometry
from . import catalogs
from . import maps
from . import signals
from . import MockSurvey
import IPython

#------------------------------------------------------------------------------------------------------------
def filterMapsAndMakeCatalogs(config, rootOutDir = None, copyFilters = False, measureFluxes = True, 
                              invertMap = False):
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
        A dictionary containing filtered maps and catalogs.
    
    Note:
        See bin/nemo for how this pipeline is applied to real data, and maps.estimateContaminationFromSkySim
        for how this is applied to source-free sims that are generated on the fly.
        
    """

    # If running on sims (source-free or with injected sources), this ensures we use the same kernels for 
    # filtering the sim maps as was used on the real data, by copying kernels to the sims dir. The kernels 
    # will then be loaded automatically when filterMaps is called. Yes, this is a bit clunky...
    if rootOutDir != None:
        dirList=[rootOutDir]
        if copyFilters == True:
            kernelCopyDestDir=rootOutDir+os.path.sep+"diagnostics"
            dirList.append(kernelCopyDestDir)
        for d in dirList:
            if os.path.exists(d) == False:
                os.makedirs(d)
        if copyFilters == True:
            for tileName in config.tileNames:
                fileNames=glob.glob(config.diagnosticsDir+os.path.sep+"filter*#%s*.fits" % (tileName))
                for f in fileNames:
                    shutil.copyfile(f, kernelCopyDestDir+os.path.sep+os.path.split(f)[-1]) 
    else:
        rootOutDir=config.rootOutDir
            
    imageDict=filters.filterMaps(config.unfilteredMapsDictList, config.parDict['mapFilters'], 
                                 tileNames = config.tileNames, rootOutDir = rootOutDir,
                                 undoPixelWindow = config.parDict['undoPixelWindow'])
    
    # Find objects in filtered maps
    photometry.findObjects(imageDict, threshold = config.parDict['thresholdSigma'], 
                           minObjPix = config.parDict['minObjPix'], 
                           findCenterOfMass = config.parDict['findCenterOfMass'], 
                           rejectBorder = config.parDict['rejectBorder'], 
                           selFnDir = config.selFnDir, objIdent = config.parDict['objIdent'], 
                           longNames = config.parDict['longNames'],
                           useInterpolator = config.parDict['useInterpolator'], 
                           measureShapes = config.parDict['measureShapes'],
                           invertMap = invertMap)
    
    # Measure fluxes
    if measureFluxes == True:
        photometry.measureFluxes(imageDict, config.parDict['photFilter'], config.diagnosticsDir, 
                                 unfilteredMapsDict = config.parDict['unfilteredMaps'],
                                 useInterpolator = config.parDict['useInterpolator'])
    else:
        # Get S/N only - if the reference (fixed) filter scale has been given
        # This is (probably) only used by maps.estimateContaminationFromInvertedMaps
        if 'photFilter' in list(config.parDict.keys()):
            photFilter=config.parDict['photFilter']
        else:
            photFilter=None
        if photFilter != None:
            photometry.getSNValues(imageDict, SNMap = 'file', prefix = 'fixed_', template = photFilter, 
                                   invertMap = invertMap)
                    
    # Merged/optimal catalogs
    catalogs.makeOptimalCatalog(imageDict, constraintsList = config.parDict['catalogCuts'])
    
    # This is useful for multi-freq, e.g., relativistic SZ corrections; tracking which objects are in 148 GHz only parts of the map
    photometry.addFreqWeightsToCatalog(imageDict, config.parDict['photFilter'], 
                                       config.diagnosticsDir)
    
    return imageDict
    
#------------------------------------------------------------------------------------------------------------
def makeMockClusterCatalog(config, numMocksToMake = 1, writeCatalogs = True, writeInfo = True, 
                           verbose = True):
    """Generate a mock cluster catalog using the given nemo config.
    
    Returns:
        List of catalogs (each is an astropy Table object)
    
    """
    
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
    tckQFitDict=signals.fitQ(config)
    
    # We only care about the filter used for fixed_ columns
    photFilterLabel=config.parDict['photFilter']
    for filterDict in config.parDict['mapFilters']:
        if filterDict['label'] == photFilterLabel:
            break
        
    # The same as was used for detecting objects
    thresholdSigma=config.parDict['thresholdSigma']

    # We need an assumed scaling relation for mock observations
    scalingRelationDict=config.parDict['massOptions']
    
    # Set-up MockSurvey objects to be used to generate multiple mocks here
    # The reason for doing this is that enableDrawSample = True makes this very slow for 3000 mass bins...
    if verbose: print(">>> Setting up mock surveys dictionary ...")
    t0=time.time()
    mockSurveyDict={}
    wcsDict={}
    RMSMapDict={}
    for tileName in config.tileNames:
        # Need area covered 
        areaImg=pyfits.open(config.selFnDir+os.path.sep+"areaMask#%s.fits.gz" % (tileName))
        areaMask=areaImg[0].data
        wcs=astWCS.WCS(areaImg[0].header, mode = 'pyfits')
        areaDeg2=(areaMask*maps.getPixelAreaArcmin2Map(areaMask, wcs)).sum()/(60**2)

        # Need RMS map to apply selection function
        RMSImg=pyfits.open(config.selFnDir+os.path.sep+"RMSMap_%s#%s.fits.gz" % (photFilterLabel, tileName))
        RMSMap=RMSImg[0].data
    
        # For a mock, we could vary the input cosmology...
        minMass=5e13
        zMin=0.0
        zMax=2.0
        H0=70.
        Om0=0.30
        Ob0=0.05
        sigma_8=0.8
        mockSurvey=MockSurvey.MockSurvey(minMass, areaDeg2, zMin, zMax, H0, Om0, Ob0, sigma_8, enableDrawSample = True)
        
        mockSurveyDict[tileName]=mockSurvey
        RMSMapDict[tileName]=RMSMap
        wcsDict[tileName]=wcs
    
    t1=time.time()
    if verbose: print("... took %.3f sec ..." % (t1-t0))
    
    if verbose: print(">>> Making mock catalogs ...")
    catList=[]
    for i in range(numMocksToMake):       
        mockTabsList=[]
        for tileName in config.tileNames:
            t0=time.time()
            mockTab=mockSurveyDict[tileName].drawSample(RMSMapDict[tileName], scalingRelationDict, tckQFitDict, wcs = wcsDict[tileName], 
                                                       photFilterLabel = photFilterLabel, tileName = tileName, makeNames = True,
                                                       SNRLimit = thresholdSigma, applySNRCut = True, 
                                                       applyPoissonScatter = applyPoissonScatter, 
                                                       applyIntrinsicScatter = applyIntrinsicScatter,
                                                       applyNoiseScatter = applyNoiseScatter)
            t1=time.time()
            mockTabsList.append(mockTab)
            if verbose: print("... making mock catalog %d for tileName = %s took %.3f sec ..." % (i+1, tileName, t1-t0))
            
        tab=atpy.vstack(mockTabsList)
        catList.append(tab)
        
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
            
    if 'combineMocks' in config.parDict.keys() and config.parDict['combineMocks'] == True:
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
