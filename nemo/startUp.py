"""

This module contains basic set-up stuff (making directories, parsing config etc.)
used by all the scripts in bin (nemo, nemoMass, nemoSelFn etc.)

"""

import os
import sys
import yaml
from . import mapTools

#------------------------------------------------------------------------------------------------------------
def startUp(parDictFileName):
    """Does start-up tasks that are common to the nemo scripts (nemo, nemoMass, nemoSelFn etc.).
    
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

    with open(parDictFileName, "r") as stream:
        parDict=yaml.safe_load(stream)

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
            raise Exception("File must have .par extension")
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
