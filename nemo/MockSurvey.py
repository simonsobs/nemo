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
import astropy.constants
from . import simsTools
from . import catalogTools
import pickle
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy import stats
from astLib import *
import time

#------------------------------------------------------------------------------------------------------------
# Conversion constants
Mpc_in_cm=astropy.constants.pc.value*100*1e6
MSun_in_g=astropy.constants.M_sun.value*1000

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
                
        # Globally change hmf's cosmology - at least, according to the docs...
        # NOTE: for astropy 2.0+ need to actually set Tcmb0 here as it defaults to zero
        # Here we use the value from Fixsen (2009): http://adsabs.harvard.edu/abs/2009ApJ...707..916F
        self.cosmo_model=FlatLambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0, Tcmb0 = 2.72548)
        cosmo.Cosmology(cosmo_model = self.cosmo_model)
        
        # For drawSample, we use astLib routines (a few times faster than cosmo_model)
        # This is just to make sure we used the same parameters
        # H0, OmegaM0, OmegaL0 used for E(z), theta500 calcs in Q
        astCalc.H0=H0
        astCalc.OMEGA_M0=Om0
        astCalc.OMEGA_L0=1.0-Om0
        
        self.minMass=minMass

        # It's much faster to generate one mass function and then update its parameters (e.g., z)
        # NOTE: Mmin etc. are log10 MSun h^-1; dndm is h^4 MSun^-1 Mpc^-3
        # Internally, it's better to stick with how hmf does this, i.e., use these units
        # Externally, we still give  inputs without h^-1
        self.mf=hmf.MassFunction(z = zRange[0], Mmin = 13., Mmax = 16., delta_wrt = 'crit', delta_h = 500.0,
                                 sigma_8 = sigma_8, cosmo_model = self.cosmo_model)#, force_flat = True, cut_fit = False)
            
        self.log10M=np.log10(self.mf.m/self.mf.cosmo.h)
        self.areaSr=areaSr
        self.areaDeg2=areaDeg2
        self.zBinEdges=zRange
        self.z=(zRange[:-1]+zRange[1:])/2.

        self.enableDrawSample=enableDrawSample
        self.update(H0, Om0, Ob0, sigma_8)
        
            
    def update(self, H0, Om0, Ob0, sigma_8):
        """Recalculate cluster counts if cosmological parameters updated.
                
        """
        # We're using both astLib and astropy... 
        # astLib is used for E(z) etc. in selFnTools where it's quicker
        # We're also keeping track inside MockSurvey itself just for convenience
        # NOTE: Working on switching all cosmology over to astropy for consistency with hmf - remove when done
        self.H0=H0
        self.Om0=Om0
        self.Ob0=Ob0
        self.sigma_8=sigma_8
        astCalc.H0=H0
        astCalc.OMEGA_M0=Om0
        astCalc.OMEGA_L0=1.0-Om0
        try:
            self.cosmo_model=FlatLambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0, Tcmb0 = 2.72548)
            #cosmo.Cosmology(cosmo_model = self.cosmo_model) # Makes no difference...
            self.mf.update(cosmo_model = self.cosmo_model, sigma_8 = sigma_8)
        except:
            raise Exception("failed to update mf when H0 = %.3f Om0 = %.3f Ob0 = %.3f sigma_8 = %.3f" % (H0, Om0, Ob0, sigma_8))
        self._doClusterCount()

        # For quick Q, fRel calc (these are in MockSurvey rather than SelFn as used by drawSample)
        self.theta500Splines=[]
        self.fRelSplines=[]
        self.Ez=self.cosmo_model.efunc(self.z)
        self.DAz=self.cosmo_model.angular_diameter_distance(self.z).value
        self.criticalDensity=self.cosmo_model.critical_density(self.z).value
        self.criticalDensity=(self.criticalDensity*np.power(Mpc_in_cm, 3))/MSun_in_g
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
            fitFRels=simsTools.calcFRel(zk, fitM500s)
            tckLog10MToTheta500=interpolate.splrep(np.log10(fitM500s), fitTheta500s)
            tckLog10MToFRel=interpolate.splrep(np.log10(fitM500s), fitFRels)
            self.theta500Splines.append(tckLog10MToTheta500)
            self.fRelSplines.append(tckLog10MToFRel)
            
        # Stuff to enable us to draw mock samples (see drawSample)
        # Interpolators here need to be updated each time we change cosmology
        if self.enableDrawSample == True:

            # Now using hmf.sample style...
            # For drawing from overall z distribution
            zSum=self.clusterCount.sum(axis = 1)
            pz=np.cumsum(zSum)/self.numClusters
            self.zRoller=_spline(pz, self.z, k = 3)
            
            # For drawing from each log10M distribution at each point on z grid
            # And quick fRel, Q calc using interpolation
            # And we may as well have E(z), DA on the z grid also
            self.log10MRollers=[]
            for i in range(len(self.z)):
                #print("updating z = %.2f" % (self.z[i]))
                self.mf.update(z = self.z[i])
                mask=self.mf.ngtm > 0
                # NOTE: / h here to match what we do in cluster count (also ensures we don't need to worry about little h in draw sample)
                self.log10MRollers.append(_spline((self.mf.ngtm[mask] / self.mf.ngtm[0])[::-1], np.log10(self.mf.m[mask][::-1]/self.mf.cosmo.h), k=3))
        

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
            try:
                n=hmf.integrate_hmf.hmf_integral_gtm(mf.m/mf.cosmo.h, mf.dndm*(mf.cosmo.h**4))  # Need to account for h^-1 in mass, h^4 in dndm
            except:
                raise Exception("Integrating hmf mass function probably failed due to mf.update using cosmo_model without Tcmb0 given?")
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
        mMask=np.greater(self.mf.m/self.mf.cosmo.h, M500Limit*1e14)
        
        return numClusters[:, mMask][zMask].sum()
        

    def getPLog10M(self, z):
        """Returns P(log10M) at given z, which corresponds to self.log10M.
        
        """

        self.mf.update(z = z)
        numberDensity=hmf.integrate_hmf.hmf_integral_gtm(self.mf.m, self.mf.dndm)
        PLog10M=numberDensity/np.trapz(numberDensity, self.mf.m)

        return PLog10M


    def drawSample(self, y0Noise, scalingRelationDict, tckQFitDict, wcs = None, photFilterLabel = None, 
                   extName = None, SNRLimit = None, makeNames = False, z = None, numDraws = None,
                   areaDeg2 = None, applySNRCut = False):
        """Draw a cluster sample from the mass function, generate mock y0~ values by applying the given 
        scaling relation parameters, and then (optionally, if both SNRLimit is given and applySNRCut 
        is True) apply the survey selection function.
        
        Here, y0Noise can be either a single number (if using e.g. survey average) or an RMSMap
        (2d array). Coords will be generated if an RMSMap and a wcs are supplied.
        
        scalingRelationDict should contain keys 'tenToA0', 'B0', 'Mpivot', 'sigma_int' (this is the
        contents of massOptions in nemo .yml config files).
        
        Most likely you should set SNRLimit to thresholdSigma, as given in the nemo .yml config file
        (the resulting catalog could be cut in z, fixed_SNR afterwards anyway).
        
        photFilterLabel and extName are only used for adding a 'template' key to each object: useful
        for quickly checking which tile (extName) an object is in.
        
        If z is given, the sample is drawn from only that redshift. The default (z = None) is to use
        the range given by self.z
        
        If numDraws = None, the number of draws is set by self.numClustersByRedshift. If numDraws is
        given, the number of draws will be divided equally between each z.
        
        If areaDeg2 is given, the cluster counts will be scaled accordingly (otherwise, they 
        correspond to self.areaDeg2). This will be ignored if numDraws is also set.
        
        This routine is used in the nemoMock script.

        Returns catalog as an astropy Table (and an array of length self.z corresponding to low mass limit
        if numDraws is used).
                
        """
                
        if z == None:
            zRange=self.z
            numClusters=int(round(self.numClustersByRedshift.sum()))
        else:
            # Pick the nearest z on the grid
            zIndex=np.argmin(abs(z-self.z))
            zRange=[self.z[zIndex]] 
            numClusters=int(round(self.numClustersByRedshift[zIndex]))
        
        if areaDeg2 != None:
            numClusters=int(round(numClusters*(areaDeg2/self.areaDeg2)))
            
        # Add Poisson noise
        numClusters=np.random.poisson(numClusters)

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
        
        # Fancy names or not?
        if makeNames == True:
            names=[]
            for RADeg, decDeg in zip(RAs, decs):
                names.append(catalogTools.makeACTName(RADeg, decDeg, prefix = 'MOCK-CL'))
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
                numClusters_zk=int(round(self.numClustersByRedshift[zIndex]))
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
            Qs_zk=interpolate.splev(theta500s_zk, tckQFitDict[extName], ext = 3)

            # For some cosmo parameters, fRel can wander outside its range for crazy masses
            # So we just cap it at 0.1 here just to avoid -ve in log
            fRels_zk=interpolate.splev(log10Ms_zk, self.fRelSplines[k], ext = 3)
            fRels_zk[fRels_zk <= 0]=0.1
            fRels_zk[fRels_zk > 1]=1.0
            
            try:
                true_y0s_zk=tenToA0*np.power(self.Ez[k], 2)*np.power(np.power(10, log10Ms_zk)/Mpivot, 1+B0)*Qs_zk*fRels_zk
                true_y0s_zk=true_y0s_zk*fRels_zk
                scattered_y0s_zk=np.exp(np.random.normal(np.log(true_y0s_zk), sigma_int, len(true_y0s_zk)))        
                measured_y0s_zk=np.random.normal(scattered_y0s_zk, y0Noise_zk)
            except:
                raise Exception("Negative y0 values (probably spline related) for H0 = %.6f Om0 = %.6f sigma_8 = %.6f at z = %.3f" % (self.H0, self.Om0, self.sigma_8, zk))
                
            true_y0s[mask]=true_y0s_zk
            measured_y0s[mask]=measured_y0s_zk
            log10Ms[mask]=log10Ms_zk
            zs[mask]=zk
                
        
        # Make table
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
        if photFilterLabel != None and extName != None:
            tab.add_column(atpy.Column([photFilterLabel+"#"+extName]*len(tab), 'template'))

        # Apply selection?
        if applySNRCut == True:
            selMask=measured_y0s > y0Noise*SNRLimit
            tab=tab[selMask]

        return tab
        
