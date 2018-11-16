"""

This module contains basic set-up stuff (making directories, parsing config etc.)
used by all the scripts in bin (nemo, nemoMass, nemoSelFn etc.)

"""

import os
import sys
import yaml
import copy
import IPython
from . import mapTools

#------------------------------------------------------------------------------------------------------------
def parseConfigFile(parDictFileName):
    """Parse a nemo .yml config file.
    
    Returns a dictionary of parameters
    
    """
    
    with open(parDictFileName, "r") as stream:
        parDict=yaml.safe_load(stream)
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
        # This is to allow source finding folks to skip this option in .yml
        # (and avoid having 'fixed_' keywords in output (they have only one filter scale)
        if 'photometryOptions' not in parDict.keys():
            parDict['photometryOptions']={}
        # We need a better way of giving defaults than this...
        if 'selFnOptions' in parDict.keys() and 'method' not in parDict['selFnOptions'].keys():
            parDict['selFnOptions']['method']='fast'
    
    return parDict
                        
#------------------------------------------------------------------------------------------------------------
def startUp(parDictFileName, ignoreMPI = False):
    """Does start-up tasks that are common to the nemo scripts (nemo, nemoMass, nemoSelFn etc.).
    
    Set ignoreMPI = True to disregard useMPI given in config file - this will return a complete set of 
    extNames (raher than just those for a given node)
    
    Returns: 
        * parDict (dictionary containing the contents of the config file)
        * rootOutDir (name of the directory where all output will be written)
        * filteredMapsDir (name of the directory where filtered maps will be written)
        * diagnosticsDir (name of the directory where e.g. filter kernels will be written)
        * unfilteredMapsDictList (list of dictionaries corresponding to maps needed)
        * extNames (list of map tiles to operate on)
        * comm, rank, size (used by MPI)
        
    """

    print(">>> Running .yml config file: %s" % (parDictFileName))

    parDict=parseConfigFile(parDictFileName)
    
    # Useful to throttle this way sometimes
    if ignoreMPI == True:
        parDict['useMPI']=False
        
    MPIEnabled=parDict['useMPI']
    if MPIEnabled == True:
        from mpi4py import MPI
        comm=MPI.COMM_WORLD
        size=comm.Get_size()
        rank=comm.Get_rank()
        if size == 1:
            raise Exception("if you want to use MPI, run with e.g., mpirun --np 4 nemo ...")
    else:
        rank=0
        comm=None
        size=1
        
    # Output dirs
    if 'outputDir' in list(parDict.keys()):
        rootOutDir=parDict['outDir']
    else:
        if parDictFileName.find(".yml") == -1:
            raise Exception("File must have .yml extension")
        rootOutDir=parDictFileName.replace(".yml", "")
    filteredMapsDir=rootOutDir+os.path.sep+"filteredMaps"
    diagnosticsDir=rootOutDir+os.path.sep+"diagnostics"
    mocksDir=rootOutDir+os.path.sep+"mocks"
    dirList=[rootOutDir, filteredMapsDir, mocksDir]
    if rank == 0:
        for d in dirList:
            if os.path.exists(d) == False:
                os.makedirs(d)

    # Optional override of default GNFW parameters (used by Arnaud model), if used in filters given
    if 'GNFWParams' not in list(parDict.keys()):
        parDict['GNFWParams']='default'
    for filtDict in parDict['mapFilters']:
        filtDict['params']['GNFWParams']=parDict['GNFWParams']

    # tileDeck file handling - either make one, or handle loading of one
    # MPI: if the tileDeck doesn't exist, only one process makes it - the others wait until it is done
    if rank == 0:
        unfilteredMapsDictList, extNames=mapTools.makeTileDeck(parDict)
        madeTileDeck=True
    else:
        madeTileDeck=None
    if MPIEnabled == True:
        madeTileDeck=comm.bcast(madeTileDeck, root = 0)
        if rank != 0 and madeTileDeck == True:
            unfilteredMapsDictList, extNames=mapTools.makeTileDeck(parDict)
            
    # For when we want to test on only a subset of tiles
    if 'extNameList' in list(parDict.keys()):
        newList=[]
        for name in extNames:
            if name in parDict['extNameList']:
                newList.append(name)
        if newList == []:
            raise Exception("extNameList given in .par file but no extensions in images match")
        extNames=newList

    # MPI: just divide up tiles pointed at by extNames among processes
    if MPIEnabled == True:
        numTilesPerNode=int(len(extNames)/size)
        startIndex=numTilesPerNode*rank
        if rank == size-1:
            endIndex=len(extNames)
        else:
            endIndex=numTilesPerNode*(rank+1)
    else:
        startIndex=0
        endIndex=len(extNames)
    extNames=extNames[startIndex:endIndex]
    
    # For debugging...
    print(("... rank = %d: extNames = %s" % (rank, str(extNames))))
    
    return parDict, rootOutDir, filteredMapsDir, diagnosticsDir, unfilteredMapsDictList, extNames, comm, rank, size
