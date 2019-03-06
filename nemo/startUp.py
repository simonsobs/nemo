"""

This module contains basic set-up stuff (making directories, parsing config etc.) used by all the scripts 
in bin/ (nemo, nemoMass, nemoSelFn etc.).

"""

import os
import sys
import yaml
import copy
import IPython
from . import maps

#------------------------------------------------------------------------------------------------------------
def parseConfigFile(parDictFileName):
    """Parse a nemo .yml config file.
    
    Args:
        parDictFileName (:obj:`str`): Path to a nemo .yml configuration file.
    
    Returns:
        A dictionary of parameters.
    
    """
    
    with open(parDictFileName, "r") as stream:
        parDict=yaml.safe_load(stream)
        # We've moved masks out of the individual map definitions in the config file
        # (makes config files simpler as we would never have different masks across maps)
        # To save re-jigging how masks are treated inside filter code, add them back to map definitions here
        maskKeys=['pointSourceMask', 'surveyMask']
        for mapDict in parDict['unfilteredMaps']:
            for k in maskKeys:
                if k in parDict.keys():
                    mapDict[k]=parDict[k]
                else:
                    mapDict[k]=None
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
        if 'photometryOptions' in parDict.keys():
            photometryOptions=parDict['photometryOptions']
            if 'photFilter' in list(photometryOptions.keys()):
                photFilter=photometryOptions['photFilter']
                for filtDict in parDict['mapFilters']:
                    if filtDict['label'] == photFilter:
                        filtDict['params']['saveRMSMap']=True
                        filtDict['params']['saveFreqWeightMap']=True
        # Don't measure object shapes by default
        if 'measureShapes' not in parDict.keys():
            parDict['measureShapes']=False
        # Don't reject objects in map border areas by default
        if 'rejectBorder' not in parDict.keys():
            parDict['rejectBorder']=0
        # This is to allow source finding folks to skip this option in .yml
        # (and avoid having 'fixed_' keywords in output (they have only one filter scale)
        if 'photometryOptions' not in parDict.keys():
            parDict['photometryOptions']={}
        # We need a better way of giving defaults than this...
        if 'selFnOptions' in parDict.keys() and 'method' not in parDict['selFnOptions'].keys():
            parDict['selFnOptions']['method']='fast'
    
    return parDict

#------------------------------------------------------------------------------------------------------------
class NemoConfig(object):
    """An object that keeps track of nemo's configuration, maps, and output directories etc..
    
    Attributes:
        parDict (:obj:`dict`): Dictionary containing the contents of the config file.
        rootOutDir (:obj:`str`): Path to the directory where all output will be written.
        filteredMapsDir (:obj:`str`): Name of the directory where filtered maps will be written.
        diagnosticsDir (:obj:`str`): Path to the directory where miscellaneous diagnostic data (e.g., filter 
            kernel plots) will be written.
        unfilteredMapsDictList (:obj:`list`): List of dictionaries corresponding to maps needed.
        extNames (:obj:`list`): List of map tiles (extension names) to operate on.
        MPIEnabled (:obj:`bool`): If True, use MPI to divide `extNames` list among processes.
        comm (:obj:`MPI.COMM_WORLD`): Used by MPI.
        rank (:obj:`int`): Used by MPI.
        size (:obj:`int`): Used by MPI.
    
    """
    
    def __init__(self, configFileName, makeOutputDirs = True, MPIEnabled = False):
        """Creates an object that keeps track of nemo's configuration, maps, output directories etc..
        
        Args:
            configFileName (:obj:`str`): Path to a nemo .yml configuration file.
            makeOutputDirs (:obj:`bool`): If True, create output directories (where maps, catalogs are stored).
            MPIEnabled (:obj:`bool`): If True, use MPI to divide the map into tiles, distributed among processes.
                This requires `tileDefinitions` and `tileNoiseRegions` to be given in the .yml config file.
    
        """

        print(">>> Running .yml config file: %s" % (configFileName))

        self.parDict=parseConfigFile(configFileName)
        self.configFileName=configFileName
        
        # We keep a copy of the original parameters dictionary in case they are overridden later and we want to
        # restore them (e.g., if running source-free sims).
        self._origParDict=copy.deepcopy(self.parDict)
                    
        self.MPIEnabled=MPIEnabled
        if self.MPIEnabled == True:
            from mpi4py import MPI
            self.comm=MPI.COMM_WORLD
            self.size=self.comm.Get_size()
            self.rank=self.comm.Get_rank()
            if self.size == 1:
                raise Exception("if you want to use MPI, run with e.g., mpiexec --np 4 nemo ...")
        else:
            self.rank=0
            self.comm=None
            self.size=1
            
        # Output dirs
        if 'outputDir' in list(self.parDict.keys()):
            self.rootOutDir=parDict['outDir']
        else:
            if configFileName.find(".yml") == -1:
                raise Exception("File must have .yml extension")
            self.rootOutDir=configFileName.replace(".yml", "")
        self.filteredMapsDir=self.rootOutDir+os.path.sep+"filteredMaps"
        self.diagnosticsDir=self.rootOutDir+os.path.sep+"diagnostics"
        self.mocksDir=self.rootOutDir+os.path.sep+"mocks"
        self.selFnDir=self.rootOutDir+os.path.sep+"selFn"
        dirList=[self.rootOutDir, self.filteredMapsDir, self.mocksDir, self.selFnDir]
        if self.rank == 0 and makeOutputDirs == True:
            for d in dirList:
                if os.path.exists(d) == False:
                    os.makedirs(d)

        # Optional override of default GNFW parameters (used by Arnaud model), if used in filters given
        if 'GNFWParams' not in list(self.parDict.keys()):
            self.parDict['GNFWParams']='default'
        for filtDict in self.parDict['mapFilters']:
            filtDict['params']['GNFWParams']=self.parDict['GNFWParams']

        # tileDeck file handling - either make one, or handle loading of one
        # MPI: if the tileDeck doesn't exist, only one process makes it - the others wait until it is done
        if self.rank == 0:
            self.unfilteredMapsDictList, self.extNames=maps.makeTileDeck(self.parDict)
            madeTileDeck=True
        else:
            madeTileDeck=None
        if self.MPIEnabled == True:
            madeTileDeck=self.comm.bcast(madeTileDeck, root = 0)
            if self.rank != 0 and madeTileDeck == True:
                self.unfilteredMapsDictList, self.extNames=maps.makeTileDeck(self.parDict)
                
        # For when we want to test on only a subset of tiles
        if 'extNameList' in list(self.parDict.keys()):
            newList=[]
            for name in self.extNames:
                if name in self.parDict['extNameList']:
                    newList.append(name)
            if newList == []:
                raise Exception("extNameList given in nemo config file but no extensions in images match")
            self.extNames=newList

        # MPI: just divide up tiles pointed at by extNames among processes
        if self.MPIEnabled == True:
            numTilesPerNode=int(len(self.extNames)/self.size)
            startIndex=numTilesPerNode*self.rank
            if self.rank == self.size-1:
                endIndex=len(self.extNames)
            else:
                endIndex=numTilesPerNode*(self.rank+1)
        else:
            startIndex=0
            endIndex=len(self.extNames)
        self.extNames=self.extNames[startIndex:endIndex]
        
        # For debugging...
        print(("... rank = %d [PID = %d]: extNames = %s" % (self.rank, os.getpid(), str(self.extNames))))
  
  
    def restoreConfig(self):
        """Restores the parameters dictionary (self.parDict) to the original state specified in the config 
        .yml file.
        
        """      
        self.parDict=copy.deepcopy(self._origParDict)
