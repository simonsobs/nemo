"""

Playing with the halo mass function... updated for the latest hmf which uses astropy

"""

import os
import sys
import numpy as np
import IPython
import astropy.table as atpy
import pylab as plt
import subprocess
import hmf
from hmf import cosmo
from astropy.cosmology import FlatLambdaCDM
from nemo import simsTools
from nemo import catalogTools
import pickle
from scipy import interpolate
from scipy import stats
from astLib import *
import time

class MockSurvey(object):
    
    def __init__(self, minMass, areaDeg2, zMin, zMax, H0, Om0, Ob0, sigma_8, enableDrawSample = False):
        """Initialise a MockSurvey object. This first calculates the probability of drawing a cluster of 
        given M500, z, assuming the Tinker mass function, and the given (generous) selection limits. 
        An additional selection function can be dialled in later when using drawSample.
        
        NOTE: We've hard coded everything to use M500 wrt critical density at this point.
        
        NOTE: MockSurvey.mf.m has factor of h^-1 in it.
                
        """

        zRange=np.linspace(zMin, zMax, 201)
        areaSr=np.radians(np.sqrt(areaDeg2))**2
                
        # Globally change hmf's cosmology - at least, according to the docs...
        # NOTE: for astropy 2.0+ need to actually set Tcmb0 here as it defaults to zero
        # Here we use the value from Fixsen (2009): http://adsabs.harvard.edu/abs/2009ApJ...707..916F
        cosmo_model=FlatLambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0, Tcmb0 = 2.72548)
        cosmo.Cosmology(cosmo_model = cosmo_model)
        
        self.minMass=minMass
        
        # It's much faster to generate one mass function and then update its parameters (e.g., z)
        # NOTE: Mmin etc. are log10 MSun h^-1; dndm is h^4 MSun^-1 Mpc^-3
        # Internally, it's better to stick with how hmf does this, i.e., use these units
        # Externally, we still give  inputs without h^-1
        self.mf=hmf.MassFunction(z = zRange[0], Mmin = 13., Mmax = 16., delta_wrt = 'crit', delta_h = 500.0,
                            sigma_8 = sigma_8, cosmo_model = cosmo_model)#, force_flat = True, cut_fit = False)
            
        self.log10M=np.log10(self.mf.m/self.mf.cosmo.h)
        self.areaSr=areaSr
        self.zBinEdges=zRange
        self.z=(zRange[:-1]+zRange[1:])/2.

        self._doClusterCount()
        
        # Stuff to enable us to draw mock samples:
        # We work with the mass function itself, and apply selection function at the end
        # Selection function gives detection probability for inclusion in sample for given fixed_SNR cut
        # We make draws = self.numClusters 
        # For each draw, we first draw z from the overall z distribution, then log10M from the mass distribution at the given z
        # We then roll 0...1 : if we roll < the detection probability score for M, z cell then we include it in the mock sample 
        if enableDrawSample == True:
            self.enableDrawSample=True
            # For drawing from overall z distribution
            zSum=self.clusterCount.sum(axis = 1)
            pz=np.cumsum(zSum)/self.numClusters
            self.tck_zRoller=interpolate.splrep(pz, self.z)
            # For drawing from each log10M distribution at each point on z grid
            # This can blow up if there are identical entries in pM... so we truncate spline fit when reach 1.0
            # NOTE: Should add sanity check code that the interpolation used here is accurate (make plots)
            #print "WARNING: should add code to check interpolation used for drawing mock samples is accurate enough"
            self.tck_log10MRoller=[]
            for i in range(len(self.z)):
                MSum=self.clusterCount[i].sum()
                pM=np.cumsum(self.clusterCount[i])/MSum
                # To avoid multiple pM == 1 blowing up the spline fit
                try:
                    maxIndex=np.min(np.where(np.diff(pM) == 0))
                except:
                    maxIndex=len(pM)
                self.tck_log10MRoller.append(interpolate.splrep(pM[:maxIndex], self.log10M[:maxIndex]))
            # Sanity check
            for i in range(len(self.z)):
                if np.any(np.isnan(self.tck_log10MRoller[i][1])) == True:
                    print("nans in self.tck_log10MRoller[%d]" % (i))
                    IPython.embed()
                    sys.exit()
                    
            
    def update(self, H0, Om0, Ob0, sigma_8):
        """Recalculate cluster counts if cosmological parameters updated.
        
        """
        cosmo_model=FlatLambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0)
        self.mf.update(cosmo_model = cosmo_model, sigma_8 = sigma_8)
        self._doClusterCount()
        

    def _doClusterCount(self):
        """Updates cluster count etc. after mass function object is updated.
        
        """
        
        mf=self.mf
        zRange=self.zBinEdges
    
        # Number density by z and total cluster count (in redshift shells)
        # Can use to make P(m, z) plane
        numberDensity=[]
        clusterCount=[]
        totalVolumeMpc3=0.
        for i in range(len(zRange)-1):
            zShellMin=zRange[i]
            zShellMax=zRange[i+1]
            zShellMid=(zShellMax+zShellMin)/2.  
            mf.update(z = zShellMid)
            n=hmf.integrate_hmf.hmf_integral_gtm(mf.m/mf.cosmo.h, mf.dndm*(mf.cosmo.h**4))  # Need to account for h^-1 in mass, h^4 in dndm
            n=abs(np.gradient(n))# Above is cumulative integral (n > m), need this for actual number count 
            numberDensity.append(n)
            shellVolumeMpc3=mf.cosmo.comoving_volume(zShellMax).value-mf.cosmo.comoving_volume(zShellMin).value
            shellVolumeMpc3=shellVolumeMpc3*(self.areaSr/(4*np.pi))
            totalVolumeMpc3=totalVolumeMpc3+shellVolumeMpc3
            clusterCount.append(n*shellVolumeMpc3)
        numberDensity=np.array(numberDensity)
        clusterCount=np.array(clusterCount)  
        self.volumeMpc3=totalVolumeMpc3
        self.numberDensity=numberDensity
        self.clusterCount=clusterCount
        self.numClusters=np.sum(clusterCount)
        self.numClustersByRedshift=np.sum(clusterCount, axis = 1)
        
        
    def addSelFn(self, selFn, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, sigma_int = 0.2):
        """Given SelFn object selFn, calculates completeness over the (self.z, self.mf.M) grid.
        
        Result stored as self.M500Completeness
        
        Can then just multiply by self.clusterCount and sum to get expected number of clusters.
        
        """
        
        self.selFn=selFn
        
        # We may need these elsewhere...
        self.scalingRelationDict={'tenToA0': tenToA0, 'B0': B0, 'Mpivot': Mpivot, 'sigma_int': sigma_int}

        # We should apply the intrinsic scatter in M500 at fixed y0~ somewhere here
        
        # This takes ~95 sec
        print("... calculating (M, z) detection probabilities in each tile (takes ~100 sec on E-D56) ...")
        self.M500Completeness=np.zeros([len(self.selFn.ycLimitTab), self.clusterCount.shape[0], self.clusterCount.shape[1]])
        t0=time.time()        
        ycLimits=self.selFn.ycLimitTab['ycLimit']
        ycErr=ycLimits/self.selFn.SNRCut
        M=(self.mf.m/self.mf.cosmo.h)
        logM=np.log10(M)
        for i in range(len(self.z)):
            z=self.z[i]
            for j in range(M.shape[0]):
                yc, theta500Arcmin, Q=simsTools.y0FromLogM500(logM[j], z, self.selFn.tckQFit, tenToA0 = tenToA0,
                                                              B0 = B0, Mpivot = Mpivot, sigma_int = sigma_int)
                self.M500Completeness[:, i, j]=stats.norm.sf(ycLimits, loc = yc, scale = ycErr)
        t1=time.time()
        
        # This takes ~7.5 sec
        M=(self.mf.m/self.mf.cosmo.h)
        logM=np.log10(M)
        self.M500Completeness_surveyAverage=np.zeros(self.clusterCount.shape)
        for i in range(len(self.z)):
            z=self.z[i]
            ycLimitAtClusterRedshift=selFn.getSurveyAverage_ycLimitAtRedshift(z)
            for j in range(M.shape[0]):
                yc, theta500Arcmin, Q=simsTools.y0FromLogM500(logM[j], z, selFn.tckQFit, tenToA0 = tenToA0,
                                                              B0 = B0, Mpivot = Mpivot, sigma_int = sigma_int)
                ycErr=ycLimitAtClusterRedshift/selFn.SNRCut
                detP=stats.norm.sf(ycLimitAtClusterRedshift, loc = yc, scale = ycErr)
                self.M500Completeness_surveyAverage[i, j]=detP


    def calcNumClustersExpected(self, M500Limit = 0.1, zMin = 0.0, zMax = 2.0, applySelFn = False, 
                                useSurveyAverageSelFn = True):
        """Calculate the number of clusters expected above a given mass limit. If applySelFn = True, apply
        the selection function (in which case M500Limit isn't important, so long as it is low).
        
        NOTE: units of M500Limit are 1e14 MSun.
        
        """
        
        if applySelFn == True:
            if useSurveyAverageSelFn == True:
                numClusters=self.M500Completeness_surveyAverage*self.clusterCount
            else:
                numClusters=0
                for i in range(len(self.selFn.ycLimitTab)):
                    numClusters=numClusters+self.M500Completeness[i]*self.clusterCount*self.selFn.ycLimitTab['fracSurveyArea'][i]
        else:
            numClusters=self.clusterCount
        
        zMask=np.logical_and(np.greater(self.z, zMin), np.less(self.z, zMax))
        mMask=np.greater(self.mf.m/self.mf.cosmo.h, M500Limit*1e14)
        
        return numClusters[:, mMask][zMask].sum()
        

    def getPLog10M(self, z):
        """Returns P(log10M) at given z, which corresponds to self.log10M.
        
        """

        self.mf.update(z = z)
        numberDensity=hmf.integrate_hmf.hmf_integral_gtm(self.mf.m, self.mf.dndm)
        PLog10M=numberDensity/np.trapz(numberDensity, self.mf.m)

        return PLog10M
    
    
    def drawSample(self, SNRLimit = None, wcs = None, areaMask = None, RMSMap = None):
        """Draw a cluster sample from the MockSurvey. We can either apply the survey averaged
        selection function, or use a user-specified area mask and noise map. If the latter, we'll
        also give mock clusters RA and dec coords while we're at it.
        
        NOTE: survey-averaged version has bugs currently (and e.g., it won't make sense to use
        SNRLimit with that, as the SelFn object has some SNRCut applied to it anyway).
        
        Returns an astropy Table object containing the mock catalog. Optional SNRCut can be
        applied to that.
        
        NOTE: units of M500 are 1e14 MSun
        
        """
        
        if SNRLimit == None:
            SNRLimit=0.
                    
        # If we want to draw coords
        if np.any(areaMask) != None and np.any(RMSMap) != None:
            ysInMask, xsInMask=np.where(areaMask != 0)
            useRMSMap=True
        else:
            useRMSMap=False

        # This takes ~16 sec for [200, 300]-shaped z, log10M500 grid using survey-averaged selection function
        #t0=time.time()
        mockCatalog=[]
        for i in range(int(self.numClusters)):
            #t0=time.time()
            #print "... %d/%d ..." % (i, self.numClusters)
            # Draw z
            zRoll=np.random.uniform(0, 1)
            z=interpolate.splev(zRoll, self.tck_zRoller)
            zIndex=np.where(abs(self.z-z) == abs(self.z-z).min())[0][0]
            # Draw M|z
            MRoll=np.random.uniform(0, 1)
            log10M=interpolate.splev(MRoll, self.tck_log10MRoller[zIndex])
            log10MIndex=np.where(abs(self.log10M-log10M) == abs(self.log10M-log10M).min())[0][0]
            
            # Using mask and noise map
            if useRMSMap == True:
                coordIndex=np.random.randint(0, len(xsInMask))
                y0Noise=RMSMap[ysInMask[coordIndex], xsInMask[coordIndex]]
                true_y0, theta500Arcmin, Q=simsTools.y0FromLogM500(log10M, z, self.selFn.tckQFit)
                measured_y0=np.exp(np.random.normal(np.log(true_y0), self.scalingRelationDict['sigma_int']))
                measured_y0=np.random.normal(measured_y0, y0Noise)
                measured_y0=true_y0+np.random.normal(0, y0Noise)
                RADeg, decDeg=wcs.pix2wcs(xsInMask[coordIndex], ysInMask[coordIndex]) # Not the bottleneck
            else:   # Survey-averaged selection function
                detP=self.M500Completeness_surveyAverage[zIndex, log10MIndex]
                PRoll=np.random.uniform(0, 1)
                if PRoll < detP:
                    # y0 from M500... we add scatter (both intrinsic and from error bar) here
                    # We should then feed this back through to get what our inferred mass would be...
                    # (i.e., we go true mass -> true y0 -> "measured" y0 -> inferred mass)
                    # NOTE: we've applied the selection function, so we've already applied the noise...?
                    # Survey-averaged, so here is the 1-sigma noise on y0~
                    y0Noise=self.selFn.ycLimit_surveyAverage.mean()/self.selFn.SNRCut
                    true_y0, theta500Arcmin, Q=simsTools.y0FromLogM500(log10M, z, self.selFn.tckQFit)
                    measured_y0=np.exp(np.random.normal(np.log(true_y0), self.scalingRelationDict['sigma_int']))
                    measured_y0=np.random.normal(measured_y0, y0Noise)
                    measured_y0=true_y0+np.random.normal(0, y0Noise)
                    RADeg=0.
                    decDeg=0.
            
            # inferred mass (full blown M500 UPP, mass function shape de-biased)
            # NOTE: we may get confused about applying de-biasing term or not here...
            if measured_y0 / y0Noise > SNRLimit:
                M500Dict=simsTools.calcM500Fromy0(measured_y0, y0Noise, z, 0.0, 
                                                    tenToA0 = self.scalingRelationDict['tenToA0'],
                                                    B0 = self.scalingRelationDict['B0'],
                                                    Mpivot = self.scalingRelationDict['Mpivot'], 
                                                    sigma_int = self.scalingRelationDict['sigma_int'],
                                                    tckQFit = self.selFn.tckQFit, mockSurvey = self, 
                                                    applyMFDebiasCorrection = True, calcErrors = True)
                # Add to catalog
                if RADeg != 0 and decDeg != 0:
                    name=catalogTools.makeACTName(RADeg, decDeg, prefix = 'MOCK-CL')
                else:
                    name='MOCK-CL %d' % (i+1)
                    
                objDict={'name': name, 
                         'RADeg': RADeg,
                         'decDeg': decDeg,
                         'true_M500': np.power(10, log10M)/1e14,
                         'true_fixed_y_c': true_y0/1e-4,
                         'fixed_y_c': measured_y0/1e-4,
                         'err_fixed_y_c': y0Noise/1e-4,
                         'fixed_SNR': measured_y0/y0Noise,
                         'redshift': float(z),
                         'redshiftErr': 0}
                for key in M500Dict:
                    objDict[key]=M500Dict[key]
                mockCatalog.append(objDict)
            
            #t1=time.time()
            #print "... took %.3f sec ..." % (t1-t0)
            
        # Convert to table
        # NOTE: left out Uncorr confusion for now... i.e., we use 'Uncorr' here...
        tab=atpy.Table()
        keyList=['name', 'redshift', 'redshiftErr', 'true_M500', 'true_fixed_y_c', 'fixed_SNR', 'fixed_y_c', 'err_fixed_y_c', 
                 'M500Uncorr', 'M500Uncorr_errPlus', 'M500Uncorr_errMinus']
        for key in keyList:
            arr=[]
            for objDict in mockCatalog:
                arr.append(objDict[key])
            tab.add_column(atpy.Column(arr, key))
        
        return tab
        
