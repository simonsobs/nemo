"""

This module defines the MockSurvey class, used for obtaining de-biased cluster mass estimates and
selection function calculations.

"""

import os
import sys
import numpy as np
import IPython
import astropy.table as atpy
import pylab as plt
import subprocess
from astropy.cosmology import FlatLambdaCDM
from colossus.cosmology import cosmology
from colossus.lss import mass_function
from . import signals
from . import catalogs
import pickle
from scipy import interpolate
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy import stats
from astLib import *
import time

#------------------------------------------------------------------------------------------------------------
class MockSurvey(object):
    
    def __init__(self, minMass, areaDeg2, zMin, zMax, H0, Om0, Ob0, sigma_8, zStep = 0.01, enableDrawSample = False):
        """Initialise a MockSurvey object. This first calculates the probability of drawing a cluster of 
        given M500, z, assuming the Tinker mass function, and the given (generous) selection limits. 
        An additional selection function can be dialled in later when using drawSample.
        
        NOTE: We've hard coded everything to use M500 wrt critical density at this point.
        
        NOTE: MockSurvey.mf.m has factor of h^-1 in it.
                
        """

        zRange=np.arange(zMin, zMax+zStep, zStep)
        areaSr=np.radians(np.sqrt(areaDeg2))**2
        self.areaSr=areaSr
        self.areaDeg2=areaDeg2
        self.zBinEdges=zRange
        self.z=(zRange[:-1]+zRange[1:])/2.
        
        params={'flat': True, 'H0': H0, 'Om0': Om0, 'Ob0': Ob0, 'sigma8': sigma_8, 'ns': 0.95}
        self.cosmoModel=cosmology.setCosmology('nemo', params)
        self.cosmoModel.checkForChangedCosmology()
        
        self.log10M=np.arange(13, 16, 0.01)
        self.M=np.power(10, self.log10M)*self.cosmoModel.h
        self.mdef='500c'
        self.model='tinker08'
        
        self.enableDrawSample=enableDrawSample
        self.update(H0, Om0, Ob0, sigma_8)
        
            
    def update(self, H0, Om0, Ob0, sigma_8):
        """Recalculate cluster counts if cosmological parameters updated.
                
        """

        self.H0=H0
        self.Om0=Om0
        self.Ob0=Ob0
        self.sigma_8=sigma_8

        params={'flat': True, 'H0': H0, 'Om0': Om0, 'Ob0': Ob0, 'sigma8': sigma_8, 'ns': 0.95}
        self.cosmoModel=cosmology.setCosmology('nemo', params)
        self.cosmoModel.checkForChangedCosmology()
        
        self._doClusterCount()

        # For quick Q, fRel calc (these are in MockSurvey rather than SelFn as used by drawSample)
        self.theta500Splines=[]
        self.fRelSplines=[]
        self.Ez=self.cosmoModel.Ez(self.z)  
        self.DAz=self.cosmoModel.angularDiameterDistance(self.z)/self.cosmoModel.h 
        self.criticalDensity=(self.cosmoModel.rho_c(self.z)*np.power(1000, 3))*np.power(self.cosmoModel.h, 2)
        for k in range(len(self.z)):
            zk=self.z[k]
            interpLim_minLog10M=self.log10M.min()
            interpLim_maxLog10M=self.log10M.max()
            interpPoints=100
            fitM500s=np.power(10, np.linspace(interpLim_minLog10M, interpLim_maxLog10M, interpPoints))
            fitTheta500s=np.zeros(len(fitM500s))
            fitFRels=np.zeros(len(fitM500s))
            criticalDensity=self.criticalDensity[k]
            DA=self.DAz[k]
            Ez=self.Ez[k]
            R500Mpc=np.power((3*fitM500s)/(4*np.pi*500*criticalDensity), 1.0/3.0)                     
            fitTheta500s=np.degrees(np.arctan(R500Mpc/DA))*60.0
            fitFRels=signals.calcFRel(zk, fitM500s, self.cosmoModel)
            tckLog10MToTheta500=interpolate.splrep(np.log10(fitM500s), fitTheta500s)
            tckLog10MToFRel=interpolate.splrep(np.log10(fitM500s), fitFRels)
            self.theta500Splines.append(tckLog10MToTheta500)
            self.fRelSplines.append(tckLog10MToFRel)
        
        # Stuff to enable us to draw mock samples (see drawSample)
        # Interpolators here need to be updated each time we change cosmology
        if self.enableDrawSample == True:

            # For drawing from overall z distribution
            zSum=self.clusterCount.sum(axis = 1)
            pz=np.cumsum(zSum)/self.numClusters
            self.zRoller=_spline(pz, self.z, k = 3)
            
            # For drawing from each log10M distribution at each point on z grid
            # And quick fRel, Q calc using interpolation
            # And we may as well have E(z), DA on the z grid also
            self.log10MRollers=[]
            for i in range(len(self.z)):
                ngtm=self._cumulativeNumberDensity(self.z[i])
                mask=ngtm > 0
                self.log10MRollers.append(_spline((ngtm[mask] / ngtm[0])[::-1], np.log10(self.M[mask][::-1]/self.cosmoModel.h), k=3))

    
    def _cumulativeNumberDensity(self, z):
        """Returns N > M (per cubic Mpc), using Colossus routines.
        
        """
    
        dndlnM=mass_function.massFunction(self.M/self.cosmoModel.h, z, mdef = self.mdef, 
                                          model = self.model, q_out = 'dndlnM')
        dndM=dndlnM/self.M
        ngtm=integrate.cumtrapz(dndlnM[::-1], np.log(self.M/self.cosmoModel.h), initial = 0)[::-1]
        
        MUpper=np.arange(np.log(self.M[-1]), np.log(10**18), np.log(self.M[1])-np.log(self.M[0]))
        extrapolator=_spline(np.log(self.M), np.log(dndlnM), k=1)
        MF_extr=extrapolator(MUpper)
        intUpper=integrate.simps(np.exp(MF_extr), dx=MUpper[2] - MUpper[1], even='first')
        ngtm=ngtm+intUpper*self.cosmoModel.h
    
        return ngtm
    
    
    def _comovingVolume(self, z):
        """Returns co-moving volume in Mpc^3 (all sky) to some redshift z, using Colossus routines (taking
        care of the fact that Colossus returns all distances in Mpc/h).
        
        NOTE: Assumes flat cosmology
        
        """
        return (4/3)*np.pi*np.power(self.cosmoModel.comovingDistance(0, z)/self.cosmoModel.h, 3)

        
    def _doClusterCount(self):
        """Updates cluster count etc. after mass function object is updated.
        
        """

        zRange=self.zBinEdges
        self.M=np.power(10, self.log10M)*self.cosmoModel.h
        
        # Number density by z and total cluster count (in redshift shells)
        # Can use to make P(m, z) plane
        numberDensity=[]
        clusterCount=[]
        totalVolumeMpc3=0.
        for i in range(len(zRange)-1):
            zShellMin=zRange[i]
            zShellMax=zRange[i+1]
            zShellMid=(zShellMax+zShellMin)/2.  
            dndlnM=mass_function.massFunction(self.M/self.cosmoModel.h, zShellMid, mdef = self.mdef, 
                                              model = self.model, q_out = 'dndlnM')
            dndM=dndlnM/self.M
            # NOTE: this differs from hmf by several % at the high-mass end (binning or interpolation?)
            n=(dndM*self.cosmoModel.h**4)*np.gradient(self.M/self.cosmoModel.h)
            numberDensity.append(n)
            shellVolumeMpc3=self._comovingVolume(zShellMax)-self._comovingVolume(zShellMin)
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


    def calcNumClustersExpected(self, M500Limit = 0.1, zMin = 0.0, zMax = 2.0, selFn = None):
        """Calculate the number of clusters expected above a given mass limit. If selFn is not None, apply
        the selection function (in which case M500Limit isn't important, so long as it is low).
        
        selFn should be an (M, z) grid that corresponds with self.log10M, self.z
        
        NOTE: units of M500Limit are 1e14 MSun.
        
        """
        
        if type(selFn) == np.ndarray:
            numClusters=selFn*self.clusterCount
        else:
            numClusters=self.clusterCount
        
        zMask=np.logical_and(np.greater(self.z, zMin), np.less(self.z, zMax))
        mMask=np.greater(self.M/self.cosmoModel.h, M500Limit*1e14)
        
        return numClusters[:, mMask][zMask].sum()
        

    def getPLog10M(self, z):
        """Returns P(log10M) at given z, which corresponds to self.log10M.
        
        """
        numberDensity=self._cumulativeNumberDensity(z)
        PLog10M=numberDensity/np.trapz(numberDensity, self.M)

        return PLog10M


    def drawSample(self, y0Noise, scalingRelationDict, tckQFitDict, wcs = None, photFilterLabel = None, 
                   tileName = None, SNRLimit = None, makeNames = False, z = None, numDraws = None,
                   areaDeg2 = None, applySNRCut = False, applyPoissonScatter = True, 
                   applyIntrinsicScatter = True, applyNoiseScatter = True):
        """Draw a cluster sample from the mass function, generate mock y0~ values by applying the given 
        scaling relation parameters, and then (optionally, if both SNRLimit is given and applySNRCut 
        is True) apply the survey selection function.
        
        Here, y0Noise can be either a single number (if using e.g. survey average) or an RMSMap
        (2d array). Coords will be generated if an RMSMap and a wcs are supplied.
        
        scalingRelationDict should contain keys 'tenToA0', 'B0', 'Mpivot', 'sigma_int' (this is the
        contents of massOptions in nemo .yml config files).
        
        Most likely you should set SNRLimit to thresholdSigma, as given in the nemo .yml config file
        (the resulting catalog could be cut in z, fixed_SNR afterwards anyway).
        
        photFilterLabel and tileName are only used for adding a 'template' key to each object: useful
        for quickly checking which tile (tileName) an object is in.
        
        If z is given, the sample is drawn from only that redshift. The default (z = None) is to use
        the range given by self.z
        
        If numDraws = None, the number of draws is set by self.numClustersByRedshift. If numDraws is
        given, the number of draws will be divided equally between each z.
        
        If areaDeg2 is given, the cluster counts will be scaled accordingly (otherwise, they 
        correspond to self.areaDeg2). This will be ignored if numDraws is also set.
        
        Use applyPoissonScatter, applyIntrinsicScatter, applyNoiseScatter to control whether Poisson 
        noise (in the expected counts / number of draws from the mass function), intrinsic scatter, 
        and/or measurement noise scatter will be applied (i.e., if all three options are set to False,
        fixed_y_c values will be the same as true_y_c, although each object will still have an error 
        bar in the output catalog, corresponding to where it is found in the RMS map).
        
        This routine is used in the nemoMock script.

        Returns catalog as an astropy Table (and an array of length self.z corresponding to low mass limit
        if numDraws is used).
                
        """
        
        # WARNING: testing only!
        #print("WARNING: np.random.seed set to fixed value in drawSample - you don't want this if not testing!")
        #np.random.seed(100)
                
        if z == None:
            zRange=self.z
        else:
            # Pick the nearest z on the grid
            zIndex=np.argmin(abs(z-self.z))
            zRange=[self.z[zIndex]] 
        
        # Add Poisson noise (we do by z to keep things simple on the z grid later)
        numClustersByRedshift=np.zeros(len(zRange), dtype = int)
        for k in range(len(zRange)):
            zk=zRange[k]
            zIndex=np.argmin(abs(zk-self.z))  
            if applyPoissonScatter == False:
                numClustersByRedshift[k]=int(round(self.numClustersByRedshift[zIndex]))
            else:
                numClustersByRedshift[k]=np.random.poisson(int(round(self.numClustersByRedshift[zIndex])))

        if areaDeg2 != None:
            numClustersByRedshift=int(round(numClustersByRedshift*(areaDeg2/self.areaDeg2)))
        
        numClusters=numClustersByRedshift.sum()
            
        if numDraws != None:
            numClusters=numDraws            

        tenToA0, B0, Mpivot, sigma_int=[scalingRelationDict['tenToA0'], scalingRelationDict['B0'], 
                                        scalingRelationDict['Mpivot'], scalingRelationDict['sigma_int']]
                    
        # If given y0Noise as RMSMap, draw coords (assuming clusters aren't clustered - which they are...)
        # NOTE: switched to using valid part of RMSMap here rather than areaMask - we need to fix the latter to same area
        # It isn't a significant issue though
        if type(y0Noise) == np.ndarray and y0Noise.ndim == 2:
            RMSMap=y0Noise
            ysInMask, xsInMask=np.where(RMSMap != 0)
            coordIndices=np.random.randint(0, len(xsInMask), numClusters)
            ys=ysInMask[coordIndices]
            xs=xsInMask[coordIndices]
            if wcs != None:
                RADecCoords=wcs.pix2wcs(xs, ys)
                RADecCoords=np.array(RADecCoords)
                RAs=RADecCoords[:, 0]
                decs=RADecCoords[:, 1]
            y0Noise=RMSMap[ys, xs]
        else:
            y0Noise=np.ones(numClusters)*y0Noise
            RAs=np.zeros(numClusters)
            decs=np.zeros(numClusters)
        
        ## WARNING: For testing only!
        #y0Noise[:]=1e-6

        # Fancy names or not?
        if makeNames == True:
            names=[]
            for RADeg, decDeg in zip(RAs, decs):
                names.append(catalogs.makeACTName(RADeg, decDeg, prefix = 'MOCK-CL'))
        else:
            names=np.arange(numClusters)+1
        
        # New way - on the redshift grid
        t0=time.time()
        currentIndex=0
        measured_y0s=np.zeros(y0Noise.shape)
        true_y0s=np.zeros(y0Noise.shape)
        log10Ms=np.zeros(y0Noise.shape)
        zs=np.zeros(y0Noise.shape)
        zErrs=np.zeros(y0Noise.shape)
        minTrueMass=np.zeros(len(zRange))   # if using numDraws
        for k in range(len(zRange)):
            
            zk=zRange[k]
            zIndex=np.argmin(abs(zk-self.z))      
            if numDraws != None:
                numClusters_zk=int(round(numDraws/len(zRange)))
            else:
                numClusters_zk=numClustersByRedshift[k]
            if numClusters_zk == 0:
                continue
            
            # Some fiddling here to avoid rounding issues with array sizes (+/-1 draw here shouldn't make a difference)
            nextIndex=currentIndex+numClusters_zk
            if nextIndex >= len(y0Noise):
                nextIndex=len(y0Noise)-1
            mask=np.arange(currentIndex, nextIndex)
            numClusters_zk=len(mask)
            y0Noise_zk=y0Noise[mask]
            currentIndex=nextIndex
                                
            # For some cosmo parameters, splined masses can end up outside of valid range, so catch this
            #log10Ms_zk=self.log10MRollers[k](np.random.uniform(0, maxDraw, numClusters_zk))
            log10Ms_zk=self.log10MRollers[k](np.random.random_sample(numClusters_zk))
            log10Ms_zk[log10Ms_zk < self.log10M.min()]=self.log10M.min()
            log10Ms_zk[log10Ms_zk > self.log10M.max()]=self.log10M.max()
            
            theta500s_zk=interpolate.splev(log10Ms_zk, self.theta500Splines[k], ext = 3)
            Qs_zk=interpolate.splev(theta500s_zk, tckQFitDict[tileName], ext = 3)

            # For some cosmo parameters, fRel can wander outside its range for crazy masses
            # So we just cap it at 0.1 here just to avoid -ve in log
            fRels_zk=interpolate.splev(log10Ms_zk, self.fRelSplines[k], ext = 3)
            fRels_zk[fRels_zk <= 0]=0.1
            fRels_zk[fRels_zk > 1]=1.0
            
            try:
                true_y0s_zk=tenToA0*np.power(self.Ez[k], 2)*np.power(np.power(10, log10Ms_zk)/Mpivot, 1+B0)*Qs_zk*fRels_zk
                if applyIntrinsicScatter == True:
                    scattered_y0s_zk=np.exp(np.random.normal(np.log(true_y0s_zk), sigma_int, len(true_y0s_zk)))
                else:
                    scattered_y0s_zk=true_y0s_zk
                if applyNoiseScatter == True:
                    measured_y0s_zk=np.random.normal(scattered_y0s_zk, y0Noise_zk)
                else:
                    measured_y0s_zk=scattered_y0s_zk
            except:
                raise Exception("Negative y0 values (probably spline related) for H0 = %.6f Om0 = %.6f sigma_8 = %.6f at z = %.3f" % (self.H0, self.Om0, self.sigma_8, zk))
                
            true_y0s[mask]=true_y0s_zk
            measured_y0s[mask]=measured_y0s_zk
            log10Ms[mask]=log10Ms_zk
            zs[mask]=zk
                
        tab=atpy.Table()
        tab.add_column(atpy.Column(names, 'name'))
        tab.add_column(atpy.Column(RAs, 'RADeg'))
        tab.add_column(atpy.Column(decs, 'decDeg'))
        tab.add_column(atpy.Column(np.power(10, log10Ms)/1e14, 'true_M500'))
        tab.add_column(atpy.Column(true_y0s/1e-4, 'true_fixed_y_c'))
        tab.add_column(atpy.Column(measured_y0s/1e-4, 'fixed_y_c'))
        tab.add_column(atpy.Column(y0Noise/1e-4, 'fixed_err_y_c'))
        tab.add_column(atpy.Column(measured_y0s/y0Noise, 'fixed_SNR'))
        tab.add_column(atpy.Column(zs, 'redshift'))
        tab.add_column(atpy.Column(zErrs, 'redshiftErr'))
        if photFilterLabel != None and tileName != None:
            tab.add_column(atpy.Column([photFilterLabel+"#"+tileName]*len(tab), 'template'))
                
        # Apply selection?
        if applySNRCut == True:
            selMask=measured_y0s > y0Noise*SNRLimit
            tab=tab[selMask]

        return tab
        
