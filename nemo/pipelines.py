"""

This module defines pipelines - sets of tasks in nemo that we sometimes want to do on different inputs
(e.g., real data or simulated data).

"""

import os
import sys
import glob
import shutil
from . import startUp
from . import filters
from . import photometry
from . import catalogs
import IPython

#------------------------------------------------------------------------------------------------------------
def filterMapsAndMakeCatalogs(config, rootOutDir = None, copyKernels = False, measureFluxes = True, 
                              invertMap = False):
    """Runs the map filtering and catalog construction steps according to the given configuration.
    
    Args:
        config (:obj: 'startup.NemoConfig'): Nemo configuration object.
        rootOutDir (str): If None, use the default given by config. Otherwise, use this to override where the
            output filtered maps and catalogs are written.
        copyKernels (bool, optional): If True, and rootOutDir is given (not None), then kernels will be
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
        if copyKernels == True:
            kernelCopyDestDir=rootOutDir+os.path.sep+"diagnostics"
            dirList.append(kernelCopyDestDir)
        for d in dirList:
            if os.path.exists(d) == False:
                os.makedirs(d)
        if copyKernels == True:
            for extName in config.extNames:
                fileNames=glob.glob(config.diagnosticsDir+os.path.sep+"kern2d*#%s*.fits" % (extName))
                for f in fileNames:
                    shutil.copyfile(f, kernelCopyDestDir+os.path.sep+os.path.split(f)[-1]) 
    else:
        rootOutDir=config.rootOutDir
            
    imageDict=filters.filterMaps(config.unfilteredMapsDictList, config.parDict['mapFilters'], 
                                 extNames = config.extNames, rootOutDir = rootOutDir,
                                 undoPixelWindow = config.parDict['undoPixelWindow'])
    
    # Find objects in filtered maps
    photometry.findObjects(imageDict, threshold = config.parDict['thresholdSigma'], 
                           minObjPix = config.parDict['minObjPix'], 
                           findCenterOfMass = config.parDict['findCenterOfMass'], 
                           rejectBorder = config.parDict['rejectBorder'], 
                           diagnosticsDir = config.diagnosticsDir, objIdent = config.parDict['objIdent'], 
                           longNames = config.parDict['longNames'],
                           useInterpolator = config.parDict['useInterpolator'], 
                           measureShapes = config.parDict['measureShapes'],
                           invertMap = invertMap)
    
    # Measure fluxes
    if measureFluxes == True:
        photometry.measureFluxes(imageDict, config.parDict['photometryOptions'], config.diagnosticsDir, 
                                 unfilteredMapsDict = config.parDict['unfilteredMaps'],
                                 useInterpolator = config.parDict['useInterpolator'])
    else:
        # Get S/N only - if the reference (fixed) filter scale has been given
        # This is (probably) only used by maps.estimateContaminationFromInvertedMaps
        if 'photFilter' in list(config.parDict['photometryOptions'].keys()):
            photFilter=config.parDict['photometryOptions']['photFilter']
        else:
            photFilter=None
        if photFilter != None:
            photometry.getSNValues(imageDict, SNMap = 'file', prefix = 'fixed_', template = photFilter, 
                                   invertMap = invertMap)
                    
    # Merged/optimal catalogs
    catalogs.makeOptimalCatalog(imageDict, constraintsList = config.parDict['catalogCuts'])
    
    # This is useful for multi-freq, e.g., relativistic SZ corrections; tracking which objects are in 148 GHz only parts of the map
    photometry.addFreqWeightsToCatalog(imageDict, config.parDict['photometryOptions'], 
                                       config.diagnosticsDir)
    
    return imageDict
    
            
            
