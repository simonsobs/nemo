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
from scipy import interpolate
from scipy import stats
from astLib import *
from nemo import simsTools
plt.matplotlib.interactive(False)

# Standard z, M ranges over which selection function is calculated
zRange=[0.05, 0.1, 0.15, 0.2, 0.3, 0.4, 0.5, 0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
MRange=np.linspace(5e13, 1e15, 20)
        
class SelFn(object):
    
    def __init__(self, nemoOutputDir, parDict, SNRCut):
        """Creates a SelFn object, which enables calculation of completeness (or detection probability) for
        given M500, z, SNR etc.
        
        Here, SNRCut would be for fixed_SNR. May be specified in the .par file (but __not__ read/used here).
        
        e.g., selFnOptions = {'fixed_SNR_cut': 5.0}
        
        Uses output from nemoSelFn, which is currently dumped in nemo diagnostics dir.
        
        """
        
        self.fitTab=atpy.Table().read(nemoOutputDir+os.path.sep+"diagnostics"+os.path.sep+"fitSimTab.fits")
        self.wcs=astWCS.WCS(nemoOutputDir+os.path.sep+"diagnostics"+os.path.sep+"areaMask.fits")
        #self.MLimits=self.fitTab['MFitSlope']*SNRCut+self.fitTab['MFitIntercept']
        self.ycLimits=self.fitTab['ycFitSlope']*SNRCut+self.fitTab['ycFitIntercept']
        self.SNRCut=SNRCut
        self.zRange=np.array(np.unique(self.fitTab['z']))
        
        # Calculate survey-averaged completeness (weighted by area of each cell)
        # We also do the same with the yc limit vs SNR fit parameters
        self.ycLimit_surveyAverage=np.zeros(self.zRange.shape)
        self.ycFitSlope_surveyAverage=np.zeros(self.zRange.shape)
        self.ycFitIntercept_surveyAverage=np.zeros(self.zRange.shape)
        for i in range(self.zRange.shape[0]):
            z=self.zRange[i]
            mask=np.equal(self.fitTab['z'], z)
            self.ycLimit_surveyAverage[i]=np.average(self.ycLimits[mask], weights = self.fitTab['fracSurveyArea'][mask])
            self.ycFitSlope_surveyAverage[i]=np.average(self.fitTab['ycFitSlope'][mask], weights = self.fitTab['fracSurveyArea'][mask])
            self.ycFitIntercept_surveyAverage[i]=np.average(self.fitTab['ycFitIntercept'][mask], weights = self.fitTab['fracSurveyArea'][mask])
        
        # Retrieves the survey-averaged mass limit at specified z, using a spline interpolation.
        self.getSurveyAverage_ycLimitAtRedshift=interpolate.InterpolatedUnivariateSpline(self.zRange, self.ycLimit_surveyAverage)

        # Similarly for ycFitSlope, ycFitIntercept (for SNRs of clusters drawn from mass function)
        self.getSurveyAverage_ycFitSlopeAtRedshift=interpolate.InterpolatedUnivariateSpline(self.zRange, self.ycFitSlope_surveyAverage)
        self.getSurveyAverage_ycFitInterceptAtRedshift=interpolate.InterpolatedUnivariateSpline(self.zRange, self.ycFitIntercept_surveyAverage)

        # For mass calculations
        diagnosticsDir=nemoOutputDir+os.path.sep+"diagnostics"
        filteredMapsDir=nemoOutputDir+os.path.sep+"filteredMaps"
        self.tckQFit=simsTools.fitQ(parDict, diagnosticsDir, filteredMapsDir)

        # For caching of MLimits_ tables
        self.diagnosticsDir=diagnosticsDir
        self.MLimits=[]
        
        
    def calcMLimits(self):
        """Makes self.MLimits, corresponding to self.ycLimits at given SNR. Takes a while to run, but result is 
        cached.
        
        You shouldn't need to run this directly - this is only needed by makeM500LimitMap
        
        """
        
        if len(self.MLimits) == 0:
            MLimitsPath=self.diagnosticsDir+os.path.sep+"MLimits_SNR%.2f.npy" % (self.SNRCut)
            if os.path.exists(MLimitsPath) == False:
                print "... cached M500 limits file %s does not exist - making it (this takes ~1 hr for E-D56) ..." % (MLimitsPath)
                count=0
                self.MLimits=np.zeros(self.ycLimits.shape[0])
                for i in range(self.ycLimits.shape[0]):
                    count=count+1
                    print "... %d/%d ..." % (count, len(self.ycLimits))
                    yc=self.ycLimits[i]
                    M500Dict=simsTools.calcM500Fromy0(yc, yc/self.SNRCut, self.fitTab['z'][i], 0.0, 
                                                      tckQFit = self.tckQFit, mockSurvey = None, 
                                                      applyMFDebiasCorrection = False, calcErrors = False)
                    self.MLimits[i]=M500Dict['M500Uncorr']
                np.save(MLimitsPath, self.MLimits)
            else:
                print "... loading cached M500 limits file %s ..." % (MLimitsPath)
                self.MLimits=np.load(MLimitsPath)      
            

    def findCell(self, RADeg, decDeg):
        """Returns the ID in self.fitTab of the cell that contains the given RA, dec coordinates.
        
        """
        x, y=self.wcs.wcs2pix(RADeg, decDeg)
        xMask=np.logical_and(np.greater(x, self.fitTab['x0']), np.less(x, self.fitTab['x1']))
        yMask=np.logical_and(np.greater(y, self.fitTab['y0']), np.less(y, self.fitTab['y1']))
        mask=np.logical_and(xMask, yMask)
        
        return int(self.fitTab['ID'][mask][0])


    def surveyAverage_ycLimitAtCompleteness(self, completenessLim):
        """Returns the survey-averaged yc limit and the corresponding zs at which it is evaluated at the 
        given completeness level.
        
        """
        
        for m in np.linspace(0.0, 2.0, 2001):
            testComplete=np.mean(stats.norm.cdf(self.ycLimit_surveyAverage*m, loc = self.ycLimit_surveyAverage, 
                                                scale = self.ycLimit_surveyAverage/self.SNRCut))
            if testComplete >= completenessLim:
                break
        
        return self.ycLimit_surveyAverage*m, self.zRange
    
    
    def surveyAverage_MLimitAtCompleteness(self, completenessLim):
        """Returns the survey-averaged M500 limit and the corresponding zs at which it is evaluated at the 
        given completeness level.
        
        """
        
        ycLimit, zRange=self.surveyAverage_ycLimitAtCompleteness(completenessLim)
        
        MLimit=np.zeros(len(ycLimit))
        for i in range(len(ycLimit)):
            M500Dict=simsTools.calcM500Fromy0(ycLimit[i], ycLimit[i]/self.SNRCut, 
                                              zRange[i], 0.0, tckQFit = self.tckQFit, mockSurvey = None, 
                                              applyMFDebiasCorrection = False, calcErrors = False)
            MLimit[i]=M500Dict['M500Uncorr']
        
        return MLimit, zRange

    
    def detectionProbability(self, yc, z, ycErr, RADeg = None, decDeg = None):
        """Returns detection probability, yc limit at SNRCut, and corresponding M500 limit at SNRCut for cluster of 
        given yc at given z. This could also be thought of as completeness.
        
        For real clusters, nemoMass gives M500, M500_errPlus, M500_errMinus, which can be used to feed in logM500Err
        (error bars from nemoMass are log normal). For simulated clusters drawn from the mass function (using 
        MockSurvey), set logM500Err according to the SNR.
        
        If RADeg, decDeg are given, the probability and mass limit is evaluated at that position. Otherwise, the
        survey-averaged detection probability is returned.
        
        """
        
        # Given position - NOTE: this can be very noisy as function of z if numPositions is small in nemoSelFn?
        if np.any(RADeg) != None and np.any(decDeg) != None:
            # These are over the range of zs in the tile we're trying
            cellID=self.findCell(RADeg, decDeg)
            cellLimits=self.ycLimits[np.equal(self.fitTab['ID'], cellID)]
            tck=interpolate.splrep(self.zRange, cellLimits)
            ycLimitAtClusterRedshift=interpolate.splev(z, tck)
            # Sanity check plot
            plt.plot(np.unique(self.fitTab['z']), cellLimits, 'r.')
            plotRange=np.linspace(0, 2, 100)
            plt.plot(plotRange, interpolate.splev(plotRange, tck), 'k-')
            
        # Survey-wide average
        else:
            ycLimitAtClusterRedshift=self.getSurveyAverage_ycLimitAtRedshift(z)

        # Defining detP as the fraction of the yc probability distribution that's over the yc selection limit
        detP=stats.norm.sf(ycLimitAtClusterRedshift, loc = yc, scale = ycErr)
        
        # Corresponding M500 limit
        M500Dict=simsTools.calcM500Fromy0(ycLimitAtClusterRedshift, ycLimitAtClusterRedshift/self.SNRCut, 
                                          z, 0.0, tckQFit = self.tckQFit, mockSurvey = None, 
                                          applyMFDebiasCorrection = False, calcErrors = False)
        
        return detP, ycLimitAtClusterRedshift, M500Dict['M500Uncorr']
        

    def ycLimitMap(self, z, wcs, outFileName = None):
        """Returns an image of the yc completeness limit. Writes to a .fits file if outFileName is given. 
        
        By definition, the limit is 50% completeness at whatever SNRCut was initially specified for the
        selFn object.
        
        """
        
        limitMap=np.zeros([wcs.header['NAXIS2'], wcs.header['NAXIS1']])
        mask=np.equal(self.fitTab['z'], z)
        for row, yc in zip(self.fitTab[mask], self.ycLimits[mask]):
            limitMap[int(row['y0']):int(row['y1']), int(row['x0']):int(row['x1'])]=yc
        
        if outFileName != None:
            astImages.saveFITS(outFileName, limitMap, wcs)
        
        return limitMap
    
        
    def makeM500LimitMap(self, z, wcs, outFileName = None):
        """Makes a .fits image map of the M500 completeness limit. 
        
        By definition, the limit is 50% completeness at whatever SNRCut was initially specified for the
        selFn object.
        
        NOTE: this is slow - takes ~500 sec per map on E-D56.
        
        """
        
        self.calcMLimits()
        
        limitMap=np.zeros([wcs.header['NAXIS2'], wcs.header['NAXIS1']])
        mask=np.equal(self.fitTab['z'], z)
        for row, M in zip(self.fitTab[mask], self.MLimits[mask]):
            limitMap[int(row['y0']):int(row['y1']), int(row['x0']):int(row['x1'])]=M
        
        if outFileName != None:
            astImages.saveFITS(outFileName, limitMap, wcs)
        
        return limitMap
        
        #---
        # Or...
        #ycLimitMap=self.ycLimitMap(z, wcs, outFileName = None)
                
        ## Slow... takes ~500 sec on E-D56
        #t0=time.time()
        #ycs=np.unique(ycLimitMap)
        #ycs=ycs[np.where(ycs != 0)]
        #MLimitMap=np.zeros(ycLimitMap.shape)
        #count=0
        #for yc in ycs:
            #t00=time.time()
            #count=count+1
            #print "... %d/%d ..." % (count, len(ycs))
            #M500Dict=simsTools.calcM500Fromy0(yc, yc/self.SNRCut, 
                                              #z, 0.0, tckQFit = self.tckQFit, mockSurvey = None, 
                                              #applyMFDebiasCorrection = False, calcErrors = False)
            #t11=time.time()
            #MLimitMap[np.equal(ycLimitMap, yc)]=M500Dict['M500Uncorr']
            #t22=time.time()
            #print t22-t00, t22-t11, t11-t00
        #t1=time.time()
        #print "total", t1-t0
        
        #if outFileName != None:
            #astImages.saveFITS(outFileName, MLimitMap, wcs)
            
        #return MLimitMap

