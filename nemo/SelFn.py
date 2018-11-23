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
        
    def __init__(self, parDictFileName, selFnDir, SNRCut, footprintLabel = None, zStep = 0.01, enableDrawSample = False,
                 downsampleRMS = True, applyMFDebiasCorrection = True):
        """Initialise a SelFn object.
        
        This is a class that uses the output from nemoSelFn to re-calculate the selection function
        (completeness fraction on (M, z) grid) with different cosmological / scaling relation parameters
        (see SelFn.update).
        
        Use footprintLabel to specify a footprint (e.g., 'DES', 'HSC', 'KiDS' etc.) defined in 
        selFnFootprints in the .yml config file. The default None uses the whole survey.
        
        Set downsampleRMS = True to speed up the completeness calculation (called each time update is called)
        considerably.
        
        """
        
        self.SNRCut=SNRCut
        self.footprintLabel=footprintLabel
        self.downsampleRMS=downsampleRMS
        self.applyMFDebiasCorrection=applyMFDebiasCorrection
        self.selFnDir=selFnDir

        self.tckQFitDict=simsTools.loadQ(self.selFnDir+os.path.sep+"QFit.pickle")
        parDict=startUp.parseConfigFile(parDictFileName)
        self.extNames=self.tckQFitDict.keys()
        
        # ignoreMPI gives us the complete list of extNames, regardless of how this parameter is set in the config file
        #config=startUp.NemoConfig(parDictFileName, makeOutputDirs = False, ignoreMPI = True)
        #parDict=config.parDict
        #self.extNames=config.extNames
        
        # Sanity check that any given footprint is defined - if not, give a useful error message
        if footprintLabel != None:
            if 'selFnFootprints' not in parDict.keys():
                raise Exception("No footprints defined in .yml config file")
            else:
                labelsList=[]
                for footprintDict in parDict['selFnFootprints']:
                    labelsList.append(footprintDict['label'])
                if footprintLabel not in labelsList:
                    raise Exception("Footprint '%s' not found in selFnFootprints - check .yml config file" % (footprintLabel))
        
        
        # We only care about the filter used for fixed_ columns
        self.photFilterLabel=parDict['photometryOptions']['photFilter']
        
        self.scalingRelationDict=parDict['massOptions']
        
        # Tile-weighted average noise and area will not change... we'll just re-calculate fitTab and put in place here
        self.selFnDictList=[]
        self.totalAreaDeg2=0.0
        for extName in self.extNames:
            RMSTab=selFnTools.getRMSTab(extName, self.photFilterLabel, self.selFnDir, footprintLabel = self.footprintLabel)
            if type(RMSTab) == atpy.Table:
                tileAreaDeg2=RMSTab['areaDeg2'].sum()
                if downsampleRMS == True:
                    RMSTab=selFnTools.downsampleRMSTab(RMSTab)
                selFnDict={'extName': extName,
                           'RMSTab': RMSTab,
                           'tileAreaDeg2': tileAreaDeg2}
                self.selFnDictList.append(selFnDict)
                self.totalAreaDeg2=self.totalAreaDeg2+tileAreaDeg2

        # Set initial fiducial cosmology - can be overridden using update function     
        minMass=5e13
        zMin=0.0
        zMax=2.0
        H0=70.
        Om0=0.30
        Ob0=0.05
        sigma_8=0.8
        self.mockSurvey=MockSurvey.MockSurvey(minMass, self.totalAreaDeg2, zMin, zMax, H0, Om0, Ob0, sigma_8, zStep = zStep, enableDrawSample = enableDrawSample)
        
        # An initial run...
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
        
        self.mockSurvey.update(H0, Om0, Ob0, sigma_8)
        
        tileAreas=[]
        compMzCube=[]
        for selFnDict in self.selFnDictList:
            tileAreas.append(selFnDict['tileAreaDeg2'])
            selFnDict['compMz']=selFnTools.calcCompleteness(selFnDict['RMSTab'], self.SNRCut, selFnDict['extName'], self.mockSurvey, 
                                                            self.scalingRelationDict, self.tckQFitDict)
            compMzCube.append(selFnDict['compMz'])
        tileAreas=np.array(tileAreas)
        fracArea=tileAreas/self.totalAreaDeg2
        compMzCube=np.array(compMzCube)
        self.compMz=np.average(compMzCube, axis = 0, weights = fracArea)
                    

    def projectCatalogToMz(self, tab):
        """Project a catalog (as astropy Table) into the (log10 M500, z) grid. Takes into account the uncertainties
        on y0, redshift - but if redshift error is non-zero, is a lot slower.
                
        Returns (log10 M500, z) grid
        
        """
        
        catProjectedMz=np.zeros(self.mockSurvey.clusterCount.shape)
        tenToA0, B0, Mpivot, sigma_int=self.scalingRelationDict['tenToA0'], self.scalingRelationDict['B0'], \
                                       self.scalingRelationDict['Mpivot'], self.scalingRelationDict['sigma_int']

        for row in tab:
            extName=row['template'].split("#")[-1]
            z=row['redshift']
            zErr=row['redshiftErr']
            y0=row['fixed_y_c']*1e-4
            y0Err=row['fixed_err_y_c']*1e-4
            P=simsTools.calcPM500(y0, y0Err, z, zErr, self.tckQFitDict[extName], self.mockSurvey, 
                                  tenToA0 = tenToA0, B0 = B0, Mpivot = Mpivot, sigma_int = sigma_int, 
                                  applyMFDebiasCorrection = self.applyMFDebiasCorrection, fRelWeightsDict = {148.0: 1.0},
                                  return2D = True)
            # Paste into (M, z) grid
            catProjectedMz=catProjectedMz+P # For return2D = True, P is normalised such that 2D array sum is 1
        
        return catProjectedMz

    
