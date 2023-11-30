"""

This module contains basic set-up stuff (making directories, parsing config etc.) used by all the scripts 
in bin/ (nemo, nemoMass, nemoSelFn etc.).

"""

import os
import sys
import yaml
import copy
import astropy.io.fits as pyfits
import astropy.table as atpy
from astLib import astWCS, astImages
from nemo import signals
import numpy as np
import pickle
import time
from . import maps

#------------------------------------------------------------------------------------------------------------
def parseConfigFile(parDictFileName, verbose = False):
    """Parse a Nemo .yml config file.
    
    Args:
        parDictFileName (:obj:`str`): Path to a nemo .yml configuration file.
        verbose (:obj:`bool`, optional): If True, warning messages may be printed to the console, if there
            are any.
    
    Returns:
        A dictionary of parameters.
    
    """
    
    if verbose:
        print(">>> Parsing config file %s" % (parDictFileName))
    with open(parDictFileName, "r") as stream:
        parDict=yaml.safe_load(stream)
        # We've moved masks out of the individual map definitions in the config file
        # (makes config files simpler as we would never have different masks across maps)
        # To save re-jigging how masks are treated inside filter code, add them back to map definitions here
        maskKeys=['pointSourceMask', 'surveyMask', 'flagMask', 'maskPointSourcesFromCatalog', 'apodizeUsingSurveyMask',
                  'maskSubtractedPointSources', 'RADecSection', 'maskHoleDilationFactor', 'reprojectToTan']
        for mapDict in parDict['unfilteredMaps']:
            for k in maskKeys:
                if k in parDict.keys():
                    mapDict[k]=parDict[k]
                else:
                    mapDict[k]=None
            # Also add key for type of weight map (inverse variance is default for enki maps)
            if 'weightsType' not in mapDict.keys():
                mapDict['weightsType']='invVar'
        # Apply global filter options (defined in allFilters) to mapFilters
        # Note that anything defined in mapFilters has priority
        # Bit ugly... we only support up to three levels of nested dictionaries...
        if 'allFilters' in parDict.keys():
            mapFiltersList=[]
            for filterDict in parDict['mapFilters']:
                newDict=copy.deepcopy(parDict['allFilters'])
                for key in filterDict.keys():
                    if type(filterDict[key]) == dict: 
                        if key not in newDict.keys():
                            newDict[key]={}
                        for subkey in filterDict[key].keys():
                            if type(filterDict[key][subkey]) == dict:
                                if subkey not in filterDict[key].keys():
                                    newDict[key][subkey]={}
                                for subsubkey in filterDict[key][subkey].keys():
                                    if type(filterDict[key][subkey][subsubkey]) == dict:
                                        if subsubkey not in filterDict[key][subkey].keys():
                                            newDict[key][subkey][subsubkey]={}                                    
                                    # No more levels please...
                                    newDict[key][subkey][subsubkey]=filterDict[key][subkey][subsubkey]                                    
                            else:
                                newDict[key][subkey]=filterDict[key][subkey]
                    else:
                        newDict[key]=filterDict[key]
                mapFiltersList.append(newDict)
            parDict['mapFilters']=mapFiltersList
        # We always need RMSMap and freqWeightsMap to do any photometry
        # So we may as well force inclusion if they have not been explicitly given
        if 'photFilter' not in parDict.keys():
            # This is to allow source finding folks to skip this option in .yml
            # (and avoid having 'fixed_' keywords in output (they have only one filter scale)
            parDict['photFilter']=None
        else:
            photFilter=parDict['photFilter']
            for filtDict in parDict['mapFilters']:
                if filtDict['label'] == photFilter:
                    filtDict['params']['saveRMSMap']=True
                    filtDict['params']['saveFreqWeightMap']=True
                    filtDict['params']['saveFilter']=True
        # Global noise mask catalog, if given, now goes into filter settings where it belongs
        if 'noiseMaskCatalog' in parDict.keys() and parDict['noiseMaskCatalog'] is not None:
            for filtDict in parDict['mapFilters']:
                filtDict['params']['noiseMaskCatalog']=parDict['noiseMaskCatalog']
        # tileNames must be case insensitive in .yml file 
        # we force upper case here (because FITS will anyway)
        if 'tileDefinitions' in parDict.keys():
            if type(parDict['tileDefinitions']) == list:
                for tileDef in parDict['tileDefinitions']:
                    tileDef['tileName']=tileDef['tileName'].upper()
        if 'tileNameList' in parDict.keys():
            newList=[]
            for entry in parDict['tileNameList']:
                newList.append(entry.upper())
            parDict['tileNameList']=newList
        if 'reprojectToTan' not in parDict.keys():
            parDict['reprojectToTan']=False
        # We shouldn't have to give this unless we're using it
        if 'catalogCuts' not in parDict.keys():
            parDict['catalogCuts']=[]
        # Don't measure object shapes by default
        if 'measureShapes' not in parDict.keys():
            parDict['measureShapes']=False
        # Don't reject objects in map border areas by default
        if 'rejectBorder' not in parDict.keys():
            parDict['rejectBorder']=0
        # By default, undo the pixel window function
        if 'undoPixelWindow' not in parDict.keys():
            parDict['undoPixelWindow']=True
        if 'fitQ' not in parDict.keys():
            parDict['fitQ']=False
        if 'calcSelFn' not in parDict.keys():
            parDict['calcSelFn']=False
        # We need a better way of giving defaults than this...
        if 'selFnOptions' in parDict.keys() and 'method' not in parDict['selFnOptions'].keys():
            parDict['selFnOptions']['method']='fast'
        if 'selFnOptions' in parDict.keys() and parDict['selFnOptions']['method'] not in ['fast', 'injection']:
                raise Exception("Valid completeness estimation methods are 'fast' or 'injection' - edit selFnOptions['method'] in config.")
        if 'selFnOptions' in parDict.keys() and 'QSource' not in parDict['selFnOptions'].keys():
            if parDict['fitQ'] == True:
                parDict['selFnOptions']['QSource']='fit'
            else:
                parDict['selFnOptions']['QSource']='injection'
        # Check of tile definitions
        if 'useTiling' not in list(parDict.keys()):
            parDict['useTiling']=False
        if 'tileDefinitions' in parDict.keys() and type(parDict['tileDefinitions']) == list:
            checkList=[]
            for entry in parDict['tileDefinitions']:
                if entry['tileName'] in checkList:
                    raise Exception("Duplicate tileName '%s' in tileDefinitions - fix in config file" % (entry['tileName']))
                checkList.append(entry['tileName'])
        if 'stitchTiles' not in list(parDict.keys()):
            if parDict['useTiling'] == True:
                parDict['stitchTiles']=True
            else:
                parDict['stitchTiles']=False
        # Optional override of default GNFW parameters (used by Arnaud model), if used in filters given
        if 'GNFWParams' not in list(parDict.keys()):
            parDict['GNFWParams']='default'
        for filtDict in parDict['mapFilters']:
            filtDict['params']['GNFWParams']=parDict['GNFWParams']
        # Optional forced photometry
        if 'forcedPhotometryCatalog' not in parDict.keys():
            parDict['forcedPhotometryCatalog']=None
        # Used for finding and removing rings around bright sources
        if 'removeRings' not in parDict.keys():
            parDict['removeRings']=True
        if 'ringThresholdSigma' not in parDict.keys():
            parDict['ringThresholdSigma']=3
        # Applies to source injection recover sims only (whether print message or trigger exception)
        if 'haltOnPositionRecoveryProblem' not in parDict.keys():
            parDict['haltOnPositionRecoveryProblem']=False
        # Mass/scaling relation/cosmology options - set fiducial values here if not chosen in config
        # NOTE: We SHOULD use M200c not M500c here (to avoid CCL Tinker08 problem)
        # But we don't, currently, as old runs/tests used M500c and Arnaud-like scaling relation
        if 'massOptions' not in parDict.keys():
            parDict['massOptions']={}
        defaults={'tenToA0': 4.95e-5, 'B0': 0.08, 'Mpivot': 3.0e+14, 'sigma_int': 0.2,
                  'relativisticCorrection': True, 'rhoType': 'critical', 'delta': 500,
                  'H0': 70.0, 'Om0': 0.3, 'Ob0': 0.05, 'sigma8': 0.80, 'ns': 0.95,
                  'concMassRelation': 'Bhattacharya13'}
        for key in defaults:
            if key not in parDict['massOptions'].keys():
                parDict['massOptions'][key]=defaults[key]

    # Stuff which is now mandatory... left here until we update everywhere in code + docs
    parDict['stitchTiles']=True

    # This isn't actually being used, but has been left in for now
    parDict['_file_last_modified_ctime']=os.path.getctime(parDictFileName)
    
    # To aid user friendliness - spot any out-of-date / removed / renamed parameters here
    # Use None for those that are totally removed
    oldKeyMap={'makeTileDir': 'useTiling', 'tileDefLabel': None, 'twoPass': None,
               'clusterInjectionModels': 'sourceInjectionModels'}
    for k in oldKeyMap.keys():
        if k in list(parDict.keys()) and oldKeyMap[k] is None:
            del parDict[k]
            if verbose:
                print("... WARNING: config parameter '%s' is no longer used by Nemo and will be ignored." % (k))
        if k in list(parDict.keys()) and type(oldKeyMap[k]) == str:
            if verbose:
                print("... WARNING: config parameter '%s' (old usage) has been renamed to '%s' (current usage) - you may wish to update your config file." % (k, oldKeyMap[k]))
            parDict[oldKeyMap[k]]=parDict[k]
            del parDict[k]

    if verbose:
        print("... config loaded successfully")

    return parDict

#------------------------------------------------------------------------------------------------------------
class NemoConfig(object):
    """An object that manages Nemo's configuration (paths to maps, output directories, filter settings etc.).
    
    Attributes:
        parDict (:obj:`dict`): Dictionary containing the contents of the config file.
        rootOutDir (:obj:`str`): Path to the directory where all output will be written.
        filteredMapsDir (:obj:`str`): Name of the directory where filtered maps will be written.
        diagnosticsDir (:obj:`str`): Path to the directory where miscellaneous diagnostic data (e.g., filter 
            kernel plots) will be written.
        unfilteredMapsDictList (:obj:`list`): List of dictionaries corresponding to maps needed.
        tileNames (:obj:`list`): List of map tiles (extension names) to operate on.
        MPIEnabled (:obj:`bool`): If True, use MPI to divide `tileNames` list among processes.
        comm (:obj:`MPI.COMM_WORLD`): Used by MPI.
        rank (:obj:`int`): Used by MPI.
        size (:obj:`int`): Used by MPI.
    
    """
    
    def __init__(self, config, makeOutputDirs = True, setUpMaps = True, writeTileInfo = False,
                 selFnDir = None, calcSelFn = False, sourceInjectionTest = False, MPIEnabled = False,
                 divideTilesByProcesses = True, verbose = True, strictMPIExceptions = True):
        """Creates an object that manages Nemo's configuration (paths to maps, output directories,
        filter settings etc.).
        
        Args:
            config (:obj:`str` or :obj:`dict`): Either the path to a nemo .yml configuration
                file, or a dictionary containing nemo configuration parameters.
            makeOutputDirs (:obj:`bool`, optional): If True, create output directories (where maps, 
                catalogs are stored).
            setUpMaps (:obj:`bool`, optional): If True, set-up data structures for handling maps 
                (inc. breaking into tiles if wanted).
            writeTileInfo (:obj:`bool`, optional): If True and set-up to break maps into tiles, write
                info on the tile geometry (such as a DS9 region file) into `selFnDir`.
            selFnDir (:obj:`str`, optional): Path to the selFn directory (use to override the 
                default location).
            calcSelFn (:obj:`bool`, optional): Overrides the value given in the config file if True.
            sourceInjectionTest (:obj:`bool`, optional): Overrides the value given in the config file
                if True.
            MPIEnabled (:obj:`bool`, optional): If True, use MPI to divide the map into tiles, 
                distributed among processes. This requires `tileDefinitions` and `tileNoiseRegions` 
                to be given in the .yml config file.
            divideTilesByProcesses (:obj:`bool`, optional): If True, divides up the list of tiles
                optimally among the available MPI processes.
            strictMPIExceptions (:obj:`bool`): If True, MPI will abort if an Exception is encountered
                (the downside is that you will not get the full traceback, but at least you will not waste
                CPU cycles). If False, MPI probably will not abort if an Exception is encountered, but you 
                will get the full traceback (the downside is your MPI job may never complete). These options 
                are a compromise due to how mpi4py handles MPI errors (the default handling for mpi4py 
                corresponds to strictMPIExceptions = False).
            verbose (:obj:`bool`): If True, print some info to the terminal while we set-up the config file.
    
        """
        self.MPIEnabled=MPIEnabled
        if self.MPIEnabled == True:
            from mpi4py import MPI
            # This is needed to get MPI to abort if one process crashes (due to mpi4py error handling)
            # If this part is disabled, we get nice exceptions, but the program will never finish if a process dies
            # Here we get the error message at least but not the traceback before MPI Aborts
            if strictMPIExceptions == True:
                sys_excepthook=sys.excepthook
                def mpi_excepthook(v, t, tb):
                    sys_excepthook(v, t, tb)
                    print("Exception: %s" % (t.args[0]))
                    MPI.COMM_WORLD.Abort(1)
                sys.excepthook=mpi_excepthook
            self.comm=MPI.COMM_WORLD
            self.size=self.comm.Get_size()
            self.rank=self.comm.Get_rank()
            if self.size == 1:
                raise Exception("If you want to use MPI, run with e.g., mpiexec -np 4 nemo configFile.yml -M")
        else:
            self.rank=0
            self.comm=None
            self.size=1

        if self.rank != 0:
            self.verbose=False
        else:
            self.verbose=verbose

        # Timekeeping for benchmarking
        if self.rank == 0:
            self._timeStarted=time.time()

        if type(config) == str:
            self.parDict=parseConfigFile(config, verbose = self.verbose)
            self.configFileName=config
        elif type(config) == dict:
            self.parDict=config
            self.configFileName=''
        else:
            raise Exception("'config' must be either a path to a .yml file, or a dictionary of parameters.")
        
        # Handle a couple of optional command-line args. These only override if set to True, otherwise ignored
        if calcSelFn == True:
            self.parDict['calcSelFn']=True
        if sourceInjectionTest == True:
            self.parDict['sourceInjectionTest']=True

        # We want the original map WCS and shape (for using stitchMaps later)
        try:
            with pyfits.open(self.parDict['unfilteredMaps'][0]['mapFileName']) as img:
                # NOTE: Zapping keywords here that appear in old ACT maps but which confuse astropy.wcs
                # Also handling compressed maps
                for ext in img:
                    if img[ext].data is not None:
                        break
                self.origWCS=astWCS.WCS(img[ext].header, mode = 'pyfits', zapKeywords = ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'])
                self.origShape=(img[ext].header['NAXIS2'], img[ext].header['NAXIS1'])
        except:
            # We don't always need or want this... should we warn by default if not found?
            self.origWCS=None
            self.origShape=None
                
        # We keep a copy of the original parameters dictionary in case they are overridden later and we want to
        # restore them (e.g., if running source-free sims).
        self._origParDict=copy.deepcopy(self.parDict)
                                
        # Output dirs
        self.rootOutDir=None
        if 'outputDir' in list(self.parDict.keys()):
            self.rootOutDir=os.path.abspath(self.parDict['outputDir'])
            # if os.path.exists(self.rootOutDir) == False: This isn't helpful really
            #     self.rootOutDir=None
        if self.rootOutDir is None:
            if self.configFileName.find(".yml") == -1 and makeOutputDirs == True:
                raise Exception("File must have .yml extension")
            self.rootOutDir=os.getcwd()+os.path.sep+os.path.split(self.configFileName.replace(".yml", ""))[-1]
        self.filteredMapsDir=self.rootOutDir+os.path.sep+"filteredMaps"
        self.diagnosticsDir=self.rootOutDir+os.path.sep+"diagnostics"
        self.selFnDir=self.rootOutDir+os.path.sep+"selFn"
        dirList=[self.rootOutDir, self.diagnosticsDir, self.filteredMapsDir, self.selFnDir]
        madeOutputDirs=None
        if self.rank == 0 and makeOutputDirs == True:
            for d in dirList:
                os.makedirs(d, exist_ok = True)
            madeOutputDirs=True

        # Optional override of selFn directory location
        if selFnDir is not None:
            self.selFnDir=selFnDir

        if setUpMaps == True:
            if self.verbose == True:
                print(">>> Setting up maps")
            self._setUpMaps(writeTileInfo = writeTileInfo)
        else:
            pickleFileName=self.selFnDir+os.path.sep+"tileCoordsDict.pkl"
            if os.path.exists(pickleFileName) == False:
                raise Exception("You can only use setUpMaps = False if a previous run has created the file %s" % (pickleFileName))
            with open(pickleFileName, "rb") as pickleFile:
                unpickler=pickle.Unpickler(pickleFile)
                self.tileCoordsDict=unpickler.load()
                self.tileNames=self.tileCoordsDict.keys()

        # For when we want to test on only a subset of tiles
        if 'tileNameList' in list(self.parDict.keys()):
            newList=[]
            for name in self.tileNames:
                if name in self.parDict['tileNameList']:
                    newList.append(name)
            if newList == []:
                raise Exception("tileNameList given in nemo config file but no extensions in images match")
            self.tileNames=newList
        else:
            # If we don't have / didn't set-up maps, we would still want the list of tile names
            if os.path.exists(self.selFnDir+os.path.sep+"tileCoordsDict.pkl") == True:
                with open(self.selFnDir+os.path.sep+"tileCoordsDict.pkl", "rb") as pickleFile:
                    unpickler=pickle.Unpickler(pickleFile)
                    tileCoordsDict=unpickler.load()
                assert(tileCoordsDict != {})
                self.tileCoordsDict=tileCoordsDict
                self.tileNames=list(tileCoordsDict.keys())

        # For convenience, keep the full list of tile names
        # (for when we don't need to be running in parallel - see, e.g., signals.getFRelWeights)
        self.allTileNames=self.tileNames.copy()

        # MPI: Divide up number of pixels to process as evenly as we can
        # For stuff with weight < 0.5: assign two tiles per process, otherwise one
        if self.MPIEnabled == True and divideTilesByProcesses == True:
            totalPix=0
            tilePix=np.zeros(len(self.tileNames))
            index=0
            for key in self.tileNames:
                tilePix[index]=self.tileCoordsDict[key]['numPix']
                totalPix=totalPix+tilePix[index]
                index=index+1
            targetPixPerProcess=totalPix/self.size
            pixTab=atpy.Table()
            pixTab['tileName']=self.tileNames
            pixTab['numPix']=tilePix
            pixTab['weight']=pixTab['numPix']/max(pixTab['numPix'])
            pixTab.sort('numPix')
            # Optimal number - if scaling isn't an issue [+1 because rank 0 doesn't process tiles]
            balancedNumProcesses=1+int(round((pixTab['weight'] <= 0.5).sum()/2))+(pixTab['weight'] > 0.5).sum()
            rankExtNames={}
            rankCounter=1
            weightCounter=0
            weightLimit=0.51
            for e in self.tileNames:
                if rankCounter not in rankExtNames:
                    rankExtNames[rankCounter]=[]
                rankExtNames[rankCounter].append(e)
                weightCounter=weightCounter+pixTab['weight'][pixTab['tileName'] == e]
                if weightCounter > weightLimit:
                    rankCounter=rankCounter+1
                    weightCounter=0
                if rankCounter > self.size-1:
                    rankCounter=1
            if self.rank in rankExtNames.keys():
                self.tileNames=rankExtNames[self.rank]
            else:
                self.tileNames=[]
            if self.rank == 0 and verbose == True:
                print(">>> Total tiles = %d ; total processes = %d ; balanced number of processes = %d" % (len(self.allTileNames), self.size, balancedNumProcesses))

        # # MPI: just divide up tiles pointed at by tileNames among processes
        # if self.MPIEnabled == True and divideTilesByProcesses == True:
        #     # New - bit clunky but distributes more evenly
        #     rankExtNames={}
        #     rankCounter=1
        #     for e in self.tileNames:
        #         if rankCounter not in rankExtNames:
        #             rankExtNames[rankCounter]=[]
        #         rankExtNames[rankCounter].append(e)
        #         rankCounter=rankCounter+1
        #         if rankCounter > self.size-1:
        #             rankCounter=1
        #     if self.rank in rankExtNames.keys():
        #         self.tileNames=rankExtNames[self.rank]
        #     else:
        #         self.tileNames=[]
        # if self.rank == 0:
        #     print(">>> Total tiles = %d ; total processes = %d" % (len(self.allTileNames), self.size))

        # We're now writing items per tile into their own dir (friendlier for Lustre)
        # NOTE: No longer writing individual tile filtered maps - only stitched versions
        if makeOutputDirs == True:
            for tileName in self.tileNames:
                for d in [self.diagnosticsDir, self.selFnDir]:
                    os.makedirs(d+os.path.sep+tileName, exist_ok = True)

        # Identify filter sets, for enabling new multi-pass filtering and object finding
        self._identifyFilterSets()

        # Cache filters in memory to save some I/O
        self.cachedFilters={}
        
        # For debugging...
        if self.MPIEnabled == True and verbose == True:
            print((">>> rank = %d [PID = %d]: tileNames = %s" % (self.rank, os.getpid(), str(self.tileNames))))
  
  
    def _identifyFilterSets(self):
        """Inspect the config dictionary to identify filter sets, which are used in multi-pass map filtering
        and object detection.

        """
        self.filterSets=[]
        if 'filterSetOptions' in self.parDict.keys():
            self.filterSetOptions=self.parDict['filterSetOptions']
            for filtDict in self.parDict['mapFilters']:
                if 'filterSets' in filtDict.keys():
                    for f in filtDict['filterSets']:
                        if f not in self.filterSets:
                            self.filterSets.append(f)
            self.filterSets.sort()
            self.filterSetLabels={}
            for setNum in self.filterSetOptions.keys():
                if 'label' in self.filterSetOptions[setNum].keys():
                    self.filterSetLabels[setNum]=self.filterSetOptions[setNum]['label']
                else:
                    self.filterSetLabels[setNum]=None


    def addAutoTileDefinitions(self, DS9RegionFileName = None, cacheFileName = None):
        """Runs the autotiler to add automatic tile definitions into the parameters dictionary in-place.

        Args:
            DS9RegionFileName (str, optional): Path to DS9 regions file to be written.
            cacheFileName (str, optional): Path to output a cached .yml file which will can be read instead on
                repeated runs (for speed).

        """

        # This is a bad idea if we're doing re-runs with -T to check tiling
        # if cacheFileName is not None and os.path.exists(cacheFileName):
        #     with open(cacheFileName, "r") as stream:
        #         self.parDict['tileDefinitions']=yaml.safe_load(stream)
        #     return None

        if 'tileDefinitions' in self.parDict.keys() and type(self.parDict['tileDefinitions']) == dict:
            # If we're not given a survey mask, we'll make one up from the map image itself
            if 'mask' in self.parDict['tileDefinitions'].keys() and self.parDict['tileDefinitions']['mask'] is not None:
                surveyMaskPath=self.parDict['tileDefinitions']['mask']
                maskProvided=True
            else:
                surveyMaskPath=self.parDict['unfilteredMaps'][0]['mapFileName']
                maskProvided=False
            if maskProvided == True:
                surveyMask, wcs=maps.chunkLoadMask(surveyMaskPath) # memory efficient
            else:
                with pyfits.open(surveyMaskPath, do_not_scale_image_data = True) as img:
                    for ext in img:
                        if ext.data is not None:
                            break
                    surveyMask=np.array(ext.data, dtype = np.uint8)
                    wcs=astWCS.WCS(ext.header, mode = 'pyfits')
                    # One day we will write a routine to deal with the multi-plane thing sensibly...
                    # But today is not that day
                    if surveyMask.ndim == 3:
                        surveyMask=surveyMask[0, :]
                    assert(surveyMask.ndim == 2)
                    surveyMask[surveyMask != 0]=1
                del img
            self._tileDefinitionsMaskPath=surveyMaskPath
            self.parDict['tileDefinitions']=maps.autotiler(surveyMask, wcs,
                                                           self.parDict['tileDefinitions']['targetTileWidthDeg'],
                                                           self.parDict['tileDefinitions']['targetTileHeightDeg'])
            if self.verbose:
                print("... breaking map into %d tiles" % (len(self.parDict['tileDefinitions'])))

            if DS9RegionFileName is not None:
                maps.saveTilesDS9RegionsFile(self.parDict, DS9RegionFileName)

            if cacheFileName is not None:
                stream=yaml.dump(self.parDict['tileDefinitions'])
                with open(cacheFileName, "w") as outFile:
                    outFile.write(stream)


    def getTileCoordsDict(self):
        """Construct a dictionary that describes how a large map is broken up into smaller tiles
        (see :ref:`Tiling` for information on the relevant configuration file parameters).

        Returns:
            A dictionary indexed by tile name, where each entry is a dictionary containing information
            on pixel coordinates of each tile within the larger map, and the WCS of each tile.

        """
        # Spin through a map, figuring out the actual coords to clip based on the tile definitions
        clipCoordsDict={}

        # We can take any map, because we earlier verified they have consistent WCS and size
        wcs=None
        try:
            wcsPath=self._tileDefinitionsMaskPath # For set-up of tiling with just a mask, nemo-sim-kit related
        except:
            wcsPath=self.parDict['unfilteredMaps'][0]['mapFileName'] # The old behavior
        with pyfits.open(wcsPath) as img:
            for ext in img:
                if ext.data is not None:
                    break
            wcs=astWCS.WCS(ext.header, mode = 'pyfits')
        assert(wcs is not None)

        # Untiled
        if self.parDict['useTiling'] == False:
            clipCoordsDict[ext.name]={'clippedSection': [0, wcs.header['NAXIS1'], 0, wcs.header['NAXIS2']],
                                      'header': wcs.header,
                                      'areaMaskInClipSection': [0, wcs.header['NAXIS1'], 0, wcs.header['NAXIS2']],
                                      'reprojectToTan': self.parDict['reprojectToTan']}

        # Tiled - this takes about 4 sec
        if self.parDict['useTiling'] == True:
            if self.verbose: print(">>> Finding tile coords")
            # Extract tile definitions (may have been inserted by autotiler before calling here)
            tileNames=[]
            coordsList=[]
            for tileDict in self.parDict['tileDefinitions']:
                ra0, ra1, dec0, dec1=tileDict['RADecSection']
                x0, y0=wcs.wcs2pix(ra0, dec0)
                x1, y1=wcs.wcs2pix(ra1, dec1)
                xMin=min([x0, x1])
                xMax=max([x0, x1])
                yMin=min([y0, y1])
                yMax=max([y0, y1])
                coordsList.append([xMin, xMax, yMin, yMax])
                tileNames.append(tileDict['tileName'])
            # Define clip regions in terms of pixels, adding overlap region
            tileOverlapDeg=self.parDict['tileOverlapDeg']
            mapData=np.ones([wcs.header['NAXIS2'], wcs.header['NAXIS1']], dtype = np.uint8)
            for c, name, tileDict in zip(coordsList, tileNames, self.parDict['tileDefinitions']):
                y0=c[2]
                y1=c[3]
                x0=c[0]
                x1=c[1]
                ra0, dec0=wcs.pix2wcs(x0, y0)
                ra1, dec1=wcs.pix2wcs(x1, y1)
                # Be careful with signs here... and we're assuming approx pixel size is ok
                if x0-tileOverlapDeg/wcs.getPixelSizeDeg() > 0:
                    ra0=ra0+tileOverlapDeg
                if x1+tileOverlapDeg/wcs.getPixelSizeDeg() < mapData.shape[1]:
                    ra1=ra1-tileOverlapDeg
                if y0-tileOverlapDeg/wcs.getPixelSizeDeg() > 0:
                    dec0=dec0-tileOverlapDeg
                if y1+tileOverlapDeg/wcs.getPixelSizeDeg() < mapData.shape[0]:
                    dec1=dec1+tileOverlapDeg
                if ra1 > ra0:
                    ra1=-(360-ra1)
                clip=astImages.clipUsingRADecCoords(mapData, wcs, ra1, ra0, dec0, dec1)

                # Storing clip coords etc. so can stitch together later
                # areaMaskSection here is used to define the region that would be kept (takes out overlap)
                ra0, dec0=wcs.pix2wcs(x0, y0)
                ra1, dec1=wcs.pix2wcs(x1, y1)
                clip_x0, clip_y0=clip['wcs'].wcs2pix(ra0, dec0)
                clip_x1, clip_y1=clip['wcs'].wcs2pix(ra1, dec1)
                clip_x0=int(round(clip_x0))
                clip_x1=int(round(clip_x1))
                clip_y0=int(round(clip_y0))
                clip_y1=int(round(clip_y1))
                if name not in clipCoordsDict:
                    clipCoordsDict[name]={'clippedSection': clip['clippedSection'], 'header': clip['wcs'].header,
                                          'areaMaskInClipSection': [clip_x0, clip_x1, clip_y0, clip_y1],
                                          'reprojectToTan': self.parDict['reprojectToTan'],
                                          'numPix': int(abs(clip_x1-clip_x0)*abs(clip_y1-clip_y0))}
                    if self.verbose:
                        print("... adding %s [%d, %d, %d, %d ; %d, %d]" % (name, ra1, ra0, dec0, dec1, ra0-ra1, dec1-dec0))

        return clipCoordsDict


    def _setUpMaps(self, writeTileInfo = False):

        if self.rank == 0:

            # Check any mask files are valid (e.g., -ve values can cause things like -ve area if not caught)
            maskKeys=['surveyMask', 'pointSourceMask']
            for key in maskKeys:
                if key in self.parDict.keys() and self.parDict[key] is not None:
                    maps.checkMask(self.parDict[key])

            # Check all maps have same WCS before we start
            self._checkWCSConsistency()

            if writeTileInfo == True:
                DS9RegionFileName=self.selFnDir+os.path.sep+"tiles.reg"
                cacheFileName=self.selFnDir+os.path.sep+"tileDefinitions.yml"
            else:
                DS9RegionFileName=None
                cacheFileName=None

            self.addAutoTileDefinitions(DS9RegionFileName = DS9RegionFileName, cacheFileName = cacheFileName)
            self.tileCoordsDict=self.getTileCoordsDict()
            assert(self.tileCoordsDict != {})

            bcastParDict=self.parDict
            bcastTileCoordsDict=self.tileCoordsDict
            if writeTileInfo == True:
                with open(self.selFnDir+os.path.sep+"tileCoordsDict.pkl", "wb") as pickleFile:
                    pickler=pickle.Pickler(pickleFile)
                    pickler.dump(self.tileCoordsDict)
        else:
            bcastParDict=None
            bcastTileCoordsDict=None
        if self.MPIEnabled == True:
            bcastParDict=self.comm.bcast(bcastParDict, root = 0)
            bcastTileCoordsDict=self.comm.bcast(bcastTileCoordsDict, root = 0)
            self.comm.barrier()
        self.tileNames=list(bcastTileCoordsDict.keys())
        self.parDict=bcastParDict
        self.tileCoordsDict=bcastTileCoordsDict

        self.unfilteredMapsDictList=maps.MapDictList(self.parDict['unfilteredMaps'],
                                                     tileCoordsDict = self.tileCoordsDict)

        # We keep a copy of the original maps set-up in case we want to override later
        self._origUnfilteredMapsDictList=copy.deepcopy(self.unfilteredMapsDictList)


    def _checkWCSConsistency(self):
        # Check consistency of WCS across maps
        mapKeys=['mapFileName', 'weightsFileName', 'pointSourceMask', 'surveyMask', 'flagMask']
        refWCS=None
        for mapDict in self.parDict['unfilteredMaps']:
            for key in mapKeys:
                if key in mapDict.keys() and mapDict[key] is not None:
                    with pyfits.open(mapDict[key]) as img:
                        wcs=None
                        for ext in img:
                            if ext is not None:
                                wcs=astWCS.WCS(img[ext].header, mode = 'pyfits')
                        if wcs is None:
                            raise Exception("Map %s doesn't have a WCS." % (mapDict[key]))
                        if refWCS is None:
                            refWCS=wcs
                        else:
                            try:
                                #assert(refWCS.getCentreWCSCoords() == wcs.getCentreWCSCoords())
                                #assert(refWCS.getImageMinMaxWCSCoords() == wcs.getImageMinMaxWCSCoords())
                                assert(refWCS.header['CTYPE1'] == wcs.header['CTYPE1'])
                                assert(refWCS.header['CTYPE2'] == wcs.header['CTYPE2'])
                                assert(refWCS.header['NAXIS1'] == wcs.header['NAXIS1'])
                                assert(refWCS.header['NAXIS2'] == wcs.header['NAXIS2'])
                                assert(refWCS.getXPixelSizeDeg() == wcs.getXPixelSizeDeg())
                                assert(refWCS.getYPixelSizeDeg() == wcs.getYPixelSizeDeg())
                            except:
                                raise Exception("WCS of %s is not consistent with other maps (all maps and masks must have the same WCS)." % (mapDict[key]))


    def restoreConfig(self):
        """Restores the parameters dictionary (self.parDict) and unfiltered maps list to the original
        state specified in the config .yml file.
        
        """      
        self.parDict=copy.deepcopy(self._origParDict)
        self.unfilteredMapsDictList=copy.deepcopy(self._origUnfilteredMapsDictList)


    def setFilterSet(self, setNum):
        """For use with multi-pass mode, in conjunction with filterSets and filterSetOptions config
        file parameters.

        """
        self.restoreConfig()

        options=None
        if setNum in self.filterSetOptions.keys():
            options=self.filterSetOptions[setNum]
            if 'saveCatalog' not in options.keys():
                options['saveCatalog']=False
            #if 'maskSubtractedRegions' not in options.keys():
                #options['maskSubtractedRegions']=False
            if 'maskHoleDilationFactor' not in options.keys():
                options['maskHoleDilationFactor']=None
            if 'addSiphonedFromSets' not in options.keys():
                options['addSiphonedFromSets']=None
            if 'ignoreSurveyMask' not in options.keys():
                options['ignoreSurveyMask']=False

        # We could add some checks here for options that don't make sense
        # e.g., if addSiphonedFromSets is present, better have subtractModelFromSets present

        permittedOverrides=['thresholdSigma', 'objIdent', 'findCenterOfMass', 'measureShapes']
        if options is not None:
            for override in permittedOverrides:
                if override in options.keys():
                    self.parDict[override]=options[override]

        # Set the filters and add catalogs used for model subtraction from filter noise term
        saveKeys=['saveFilteredMaps', 'saveFilter', 'saveRMSMap', 'savePlots', 'saveDS9Regions']
        filtersToActivate=[]
        for filtDict in self.parDict['mapFilters']:
            if setNum in filtDict['filterSets']:
                if options is not None:
                    if 'mapToUse' in options.keys():
                        filtDict['params']['mapToUse']=options['mapToUse']
                    if 'noiseModelCatalogFromSets' in options.keys():
                        filtDict['params']['noiseModelCatalog']=[]
                        for noiseSubIndex in options['noiseModelCatalogFromSets']:
                            filtDict['params']['noiseModelCatalog'].append(self.filterSetOptions[noiseSubIndex]['catalog'])
                # NOTE: We prevent any filter set apart from the last one from writing maps, filters to disk
                # Similarly, no point doing forced photometry on any run except the last
                if setNum != self.filterSets[-1]:
                    for saveKey in saveKeys:
                        if saveKey in filtDict['params'].keys():
                            filtDict['params'][saveKey]=False
                    self.parDict['forcedPhotometryCatalog']=None
                # However, we do allow intermediate filter sets to write maps, if explicitly asked for
                # (e.g., useful for debugging a point source run before a cluster run)
                if type(options) == dict and 'saveFilteredMaps' in options.keys():
                    filtDict['params']['saveFilteredMaps']=options['saveFilteredMaps']
                filtersToActivate.append(filtDict)
        self.parDict['mapFilters']=filtersToActivate

        # Add catalogs for model subtraction from maps
        if options is not None and 'subtractModelFromSets' in options.keys():
            for mapDict in self.unfilteredMapsDictList:
                for subtractSetIndex in options['subtractModelFromSets']:
                    if 'mapToUse' in self.filterSetOptions[subtractSetIndex].keys():
                        if 'label' in mapDict.keys() and mapDict['label'] != self.filterSetOptions[subtractSetIndex]['mapToUse']:
                            continue
                    mapDict['subtractModelFromCatalog']=self.filterSetOptions[subtractSetIndex]['catalog']

        # Add catalogs for masking and filling
        if options is not None and 'maskAndFillFromSets' in options.keys():
            for mapDict in self.unfilteredMapsDictList:
                for subtractSetIndex in options['maskAndFillFromSets']:
                    if 'mapToUse' in self.filterSetOptions[subtractSetIndex].keys():
                        if 'label' in mapDict.keys() and mapDict['label'] != self.filterSetOptions[subtractSetIndex]['mapToUse']:
                            continue
                    mapDict['maskAndFillFromCatalog']=self.filterSetOptions[subtractSetIndex]['catalog']

        # Other map-level preprocessing keys
        if options is not None:
            for mapDict in self.unfilteredMapsDictList:
                #mapDict['maskSubtractedRegions']=options['maskSubtractedRegions']
                mapDict['maskHoleDilationFactor']=options['maskHoleDilationFactor']
                if options['ignoreSurveyMask'] == True:
                    mapDict['surveyMask']=None
