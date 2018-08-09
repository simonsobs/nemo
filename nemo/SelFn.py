"""

A class which handles application of the selection function, as calculated by nemoSelFn.

"""

import os
import sys
import numpy as np
import IPython
import astropy.table as atpy
from astLib import *
import pylab as plt
import astropy.io.fits as pyfits
from scipy import interpolate
from scipy import stats
from astLib import *
from nemo import simsTools
from nemo import selFnTools
from nemo import MockSurvey
from nemo import startUp

class SelFn(object):
        
    def __init__(self, parDictFileName, SNRCut):
        """Initialise a SelFn object.
        
        This is a class that uses the output from nemoSelFn to re-calculate the selection function
        (completeness fraction on (M, z) grid) with different cosmological / scaling relation parameters
        (see SelFn.update).
                
        """
        
        self.SNRCut=SNRCut
        
        # ignoreMPI gives us the complete list of extNames, regardless of how this parameter is set in the config file
        parDict, rootOutDir, filteredMapsDir, diagnosticsDir, unfilteredMapsDictList, extNames, comm, rank, size=startUp.startUp(parDictFileName, ignoreMPI = True)
        
        self.diagnosticsDir=diagnosticsDir
        self.extNames=extNames
        self.tckQFitDict=simsTools.fitQ(parDict, diagnosticsDir, filteredMapsDir)
        
        # We only care about the filter used for fixed_ columns
        self.photFilterLabel=parDict['photometryOptions']['photFilter']
        
        self.scalingRelationDict=parDict['massOptions']
        
        # Tile-weighted average noise and area will not change... we'll just re-calculate fitTab and put in place here
        self.selFnDictList=[]
        self.totalAreaDeg2=0.0
        for extName in self.extNames:
            tileAreaDeg2=selFnTools.getTileTotalAreaDeg2(extName, self.diagnosticsDir)
            y0Noise=selFnTools.calcTileWeightedAverageNoise(extName, self.photFilterLabel, self.diagnosticsDir)
            selFnDict={'extName': extName,
                       'y0Noise': y0Noise,
                       'tileAreaDeg2': tileAreaDeg2}
            self.selFnDictList.append(selFnDict)
            self.totalAreaDeg2=self.totalAreaDeg2+tileAreaDeg2

        # Since a fiducial cosmology (OmegaM0 = 0.3, OmegaL0 = 0.7, H0 = 70 km/s/Mpc) was used in the object detection/filtering stage, we use the same one here      
        minMass=5e13
        zMin=0.0
        zMax=2.0
        H0=70.
        Om0=0.30
        Ob0=0.05
        sigma_8=0.8
        self.mockSurvey=MockSurvey.MockSurvey(minMass, self.totalAreaDeg2, zMin, zMax, H0, Om0, Ob0, sigma_8, enableDrawSample = True)
        
        #t1=time.time()
        #print("SelFn set-up time: %.3f sec" % (t1-t0))
        
        # An inital run...
        self.update(H0, Om0, Ob0, sigma_8)


    def update(self, H0, Om0, Ob0, sigma_8, scalingRelationDict = None):
        """Re-calculates survey-average selection function given new set of cosmological / scaling relation parameters.
        
        This updates self.mockSurvey at the same time - i.e., this is the only thing that needs to be updated.
        
        Resulting (M, z) completeness grid stored as self.compMz
        
        To apply the selection function and get the expected number of clusters in the survey do e.g.:
        
        selFn.mockSurvey.calcNumClustersExpected(selFn = selFn.compMz)
        
        (yes, this is a bit circular)
        
        """
        
        if scalingRelationDict != None:
            self.scalingRelationDict=scalingRelationDict
        
        #t0=time.time()
        self.mockSurvey.update(H0, Om0, Ob0, sigma_8)
        tileAreas=[]
        compMzCube=[]
        for selFnDict in self.selFnDictList:
            tileAreas.append(selFnDict['tileAreaDeg2'])
            selFnDict['fitTab']=selFnTools.calcCompleteness(selFnDict['y0Noise'], self.SNRCut, selFnDict['extName'], self.mockSurvey, 
                                                            self.scalingRelationDict, self.tckQFitDict, self.diagnosticsDir)
            compMzCube.append(selFnTools.makeMzCompletenessGrid(selFnDict['fitTab'], self.mockSurvey))
        tileAreas=np.array(tileAreas)
        fracArea=tileAreas/self.totalAreaDeg2
        compMzCube=np.array(compMzCube)
        self.compMz=np.average(compMzCube, axis = 0, weights = fracArea)
        #t1=time.time()

        # If testing
        #selFnTools.makeMzCompletenessPlot(self.compMz, self.mockSurvey.log10M, self.mockSurvey.z, "test", "testMz.pdf")
        #print("time taken = %.3f sec" % (t1-t0))
        #print("update")
        #IPython.embed()
        #sys.exit()
