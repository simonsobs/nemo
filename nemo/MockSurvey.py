"""Playing with the halo mass function... updated for the latest hmf which uses astropy


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
import cPickle
from scipy import interpolate
from scipy import stats
from astLib import *
import time
plt.matplotlib.interactive(False)

class MockSurvey(object):
    
    def __init__(self, minMass, areaDeg2, zMin, zMax, H0, Om0, Ob0, sigma_8):
        """Initialise a MockSurvey object. This first calculates the probability of drawing a cluster of 
        given M500, z, assuming the Tinker mass function, and the given (generous) selection limits. 
        An additional selection function can be dialled in later when using drawSample.
        
        NOTE: We've hard coded everything to use M500 wrt critical density at this point.
        
        NOTE: MockSurvey.mf.m has factor of h^-1 in it.
                
        """

        zRange=np.linspace(zMin, zMax, 201)
        areaSr=np.radians(np.sqrt(areaDeg2))**2
                
        # Globally change hmf's cosmology - at least, according to the docs...
        cosmo_model=FlatLambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0)
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

        # We should apply the intrinsic scatter in M500 at fixed y0~ somewhere here
        
        # This takes ~95 sec
        print "... calculating (M, z) detection probabilities in each tile (takes ~100 sec on E-D56) ..."
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
    
    
    def drawSample(self):
        """Draw a cluster sample from the MockSurvey. Returns it as a table object, with columns M500, z
        
        NOTE: units of M500 are 1e14 MSun
        
        """
        
        print "new draw sample"
        IPython.embed()
        sys.exit()
        ## Drawing a sample from our mass function
        #N=100
        #icdf=interpolate.InterpolatedUnivariateSpline((mf.ngtm / mf.ngtm[0])[::-1], np.log10(mf.m[::-1]), k=3)
        #x=np.random.random(N)
        #m=10**icdf(x)
        ## To do the reverse of the above ^^^
        #cdf=interpolate.InterpolatedUnivariateSpline(np.log10(mf.m), (mf.ngtm / mf.ngtm[0]), k=3)
        #PLog10M=cdf(13.5)

        
            
    def drawSample_old(self, dlog10M = 0.01, dz = 0.001):
        """Draws a cluster sample from the MockSurvey. Returns it as an atpy.Table object, with columns
        name, z, dz, M500, dM500.
        
        NOTE: units of M500 are 1e14 MSun
        
        """
        
        if self.enableDrawSample == False:
            raise Exception, "enableDrawSample == False"
        
        zRolls=np.random.uniform(0, 1, self.numClusters)
        mRolls=np.random.uniform(0, 1, self.numClusters)
        
        zs=interpolate.splev(zRolls, self.tckRedshift)
        ms=[]
        for z, mRoll in zip(zs, mRolls):
            zMidIndex=np.where(abs(self.z-z) == abs(self.z-z).min())[0][0]
            ms.append(interpolate.splev(mRoll, self.tckLog10MByRedshift[zMidIndex]))
        ms=np.array(ms)
        msErr=(np.power(10, ms+dlog10M)-np.power(10, ms))/1e14
        ms=np.power(10, ms)/1e14
        
        tab=atpy.Table()
        tab.table_name="mock sample"
        tab.add_column(atpy.Column(zs, "z"))
        tab.add_column(atpy.Column([dz]*len(tab), "dz"))
        tab.add_column(atpy.Column(ms, "M500"))
        tab.add_column(atpy.Column(msErr, "dM500"))
        
        return tab
        
