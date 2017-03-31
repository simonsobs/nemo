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
plt.matplotlib.interactive(False)

class SelFn(object):
    
    def __init__(self, nemoOutputDir, SNRCut):
        """Creates a SelFn object, which enables calculation of completeness (or detection probability) for
        given M500, z, SNR etc.
        
        Here, SNRCut would be for fixed_SNR. May be specified in the .par file (but not read here).
        
        e.g., selFnOptions = {'fixed_SNR_cut': 5.0}
        
        Uses output from nemoSelFn, which is currently dumped in nemo diagnostics dir.
        
        """
        
        self.fitTab=atpy.Table().read(nemoOutputDir+os.path.sep+"diagnostics"+os.path.sep+"fitSimTab.fits")
        self.wcs=astWCS.WCS(nemoOutputDir+os.path.sep+"diagnostics"+os.path.sep+"areaMask.fits")
        self.MLimits=self.fitTab['MFitSlope']*SNRCut+self.fitTab['MFitIntercept']
        self.ycLimits=self.fitTab['ycFitSlope']*SNRCut+self.fitTab['ycFitIntercept']
        self.SNRCut=SNRCut
        self.zRange=np.array(np.unique(self.fitTab['z']))
        
        # Calculate survey-averaged completeness (weighted by area of each cell)
        self.surveyAverageMLimit=np.zeros(self.zRange.shape)
        for i in range(self.zRange.shape[0]):
            z=self.zRange[i]
            mask=np.equal(self.fitTab['z'], z)
            self.surveyAverageMLimit[i]=np.average(self.MLimits[mask], weights = self.fitTab['fracSurveyArea'][mask])
        
        # Retrieves the survey-averaged mass limit at specified z, using a spline interpolation.
        self.getSurveyAverageMLimitAtRedshift=interpolate.InterpolatedUnivariateSpline(self.zRange, self.surveyAverageMLimit)


    def findCell(self, RADeg, decDeg):
        """Returns the ID in self.fitTab of the cell that contains the given RA, dec coordinates.
        
        """
        x, y=self.wcs.wcs2pix(RADeg, decDeg)
        xMask=np.logical_and(np.greater(x, self.fitTab['x0']), np.less(x, self.fitTab['x1']))
        yMask=np.logical_and(np.greater(y, self.fitTab['y0']), np.less(y, self.fitTab['y1']))
        mask=np.logical_and(xMask, yMask)
        
        return int(self.fitTab['ID'][mask][0])

    
    def surveyAverageMLimitAtCompleteness(self, completenessLim):
        """Returns the survey-averaged M500 limit and the corresponding zs at which it is evaluated at the 
        given completeness level.
        
        """
        
        for m in np.linspace(1.0, 2.0, 1001):
            testComplete=np.mean(stats.norm.cdf(np.log10(self.surveyAverageMLimit*m), loc = np.log10(self.surveyAverageMLimit), 
                                                scale = np.log10(1.+1./self.SNRCut)))
            if testComplete > completenessLim:
                break
        
        return np.power(10, np.log10(self.surveyAverageMLimit*m)), self.zRange
        
        
    def logM500DetectionProbability(self, logM500, z, logM500Err, RADeg = None, decDeg = None):
        """Returns detection probability and mass limit for cluster of given log10(M500 / 10^14 MSun) at given z. 
        This could also be thought of as completeness.
        
        For real clusters, nemoMass gives M500, M500_errPlus, M500_errMinus, which can be used to feed in logM500Err
        (error bars from nemoMass are log normal). For simulated clusters drawn from the mass function (using 
        MockSurvey), set logM500Err according to the SNR.
        
        If RADeg, decDeg are given, the probability and mass limit is evaluated at that position. Otherwise, the
        survey-averaged detection probability is returned.
        
        """
               
        # Given position
        if np.any(RADeg) != None and np.any(decDeg) != None:
            # These are over the range of zs in the tile we're trying
            cellID=self.findCell(RADeg, decDeg)
            cellMLimits=self.MLimits[np.equal(self.fitTab['ID'], cellID)]
            tck=interpolate.splrep(self.zRange, cellMLimits)
            MLimitAtClusterRedshift=interpolate.splev(z, tck)
            # Sanity check plot
            #plt.plot(np.unique(self.fitTab['z']), cellMLimits, 'r.')
            #plotRange=np.linspace(0, 2, 100)
            #plt.plot(plotRange, interpolate.splev(plotRange, tck), 'k-')
            
        # Survey-wide average
        else:
            MLimitAtClusterRedshift=self.getSurveyAverageMLimitAtRedshift(z)

        # Defining detP as the fraction of the M500 probability distribution that's over the mass selection limit
        detP=stats.norm.sf(np.log10(MLimitAtClusterRedshift), loc = logM500, scale = logM500Err)
        
        return detP, MLimitAtClusterRedshift
    
    
    def y0Completeness(self, y0, z):
        """Returns completeness for clusters of given y0~ (fixed_y_c in nemo speak; units 1e-4) at given 
        z. This could also be thought of as detection probability.
        
        Note that z can be an array, but M500 cannot.

        NOTE: This routine is also currently wrong and hasn't been fixed yet...
        
        """
        
        print "Fix y0completeness calc"
        IPython.embed()
        sys.exit()
        
        # We could avoid some duplication here...
        zRange=np.unique(self.fitTab['z'])
        completeness=np.zeros(zRange.shape)
        for i in range(zRange.shape[0]):
            zMask=np.equal(self.fitTab['z'], zRange[i])
            selMask=np.less(self.ycLimits[zMask], y0)
            # Weighted by area of each tile
            completeness[i]=np.sum(self.fitTab['fracSurveyArea'][zMask][selMask])/np.sum(self.fitTab['fracSurveyArea'][zMask])
        tck=interpolate.splrep(zRange, completeness)
        
        return interpolate.splev(z, tck)
            
            
    def makeM500LimitMap(self, z, wcs, outFileName):
        """Makes a .fits image map of the M500 completeness limit. 
        
        By definition, the limit is 50% completeness at whatever SNRCut was initially specified for the
        selFn object.
        
        """
        
        limitMap=np.zeros([wcs.header['NAXIS2'], wcs.header['NAXIS1']])
        mask=np.equal(self.fitTab['z'], z)
        for row, M in zip(self.fitTab[mask], self.MLimits[mask]):
            limitMap[int(row['y0']):int(row['y1']), int(row['x0']):int(row['x1'])]=M
        astImages.saveFITS(outFileName, limitMap, wcs)

