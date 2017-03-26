"""

A class which handles application of the selection function, as calculated by nemoSelFn


"""

import os
import sys
import numpy as np
import IPython
import astropy.table as atpy
import pylab as plt
from scipy import interpolate
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
        self.MLimits=self.fitTab['MFitSlope']*SNRCut+self.fitTab['MFitIntercept']
        self.ycLimits=self.fitTab['ycFitSlope']*SNRCut+self.fitTab['ycFitIntercept']
        self.SNRCut=SNRCut
    
        
    def M500Completeness(self, M500, z):
        """Returns completeness for cluster of given M500 (in 10^14 MSun) at given z. This could also be 
        thought of as detection probability.
        
        Note that z can be an array, but M500 cannot.
        
        """
        
        # This is too slow currently
        zRange=np.unique(self.fitTab['z'])
        completeness=np.zeros(zRange.shape)
        for i in range(zRange.shape[0]):
            zMask=np.equal(self.fitTab['z'], zRange[i])
            selMask=np.less(self.MLimits[zMask], M500)
            # Weighted by area of each tile
            completeness[i]=np.sum(self.fitTab['fracSurveyArea'][zMask][selMask])/np.sum(self.fitTab['fracSurveyArea'][zMask])
        tck=interpolate.splrep(zRange, completeness)
        
        # Sanity check plot
        #plotRange=np.linspace(zRange.min(), zRange.max(), 1000)
        #plt.plot(zRange, completeness, 'r.')
        #plt.plot(plotRange, interpolate.splev(plotRange, tck), 'k-')
        #plt.xlabel("z")
        #plt.ylabel("completeness")
        #plt.title("M500 = %.2f x 10$^{14}$ M$_{\odot}$, z = %.3f" % (M500, z))
        
        return interpolate.splev(z, tck)
    
    
    def y0Completeness(self, y0, z):
        """Returns completeness for clusters of given y0~ (fixed_y_c in nemo speak; units 1e-4) at given 
        z. This could also be thought of as detection probability.
        
        Note that z can be an array, but M500 cannot.

        """
        
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
        
        """
        
        limitMap=np.zeros([wcs.header['NAXIS2'], wcs.header['NAXIS1']])
        mask=np.equal(self.fitTab['z'], z)
        for row, M in zip(self.fitTab[mask], self.MLimits[mask]):
            limitMap[int(row['y0']):int(row['y1']), int(row['x0']):int(row['x1'])]=M
        astImages.saveFITS(outFileName, limitMap, wcs)

