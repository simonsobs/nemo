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
from . import simsTools
from . import catalogTools
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
                                 sigma_8 = sigma_8, cosmo_model = cosmo_model)#, force_flat = True, cut_fit = False)
            
        self.log10M=np.log10(self.mf.m/self.mf.cosmo.h)
        self.areaSr=areaSr
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
        self.H0=H0
        self.Om0=Om0
        self.Ob0=Ob0
        self.sigma_8=sigma_8
        astCalc.H0=H0
        astCalc.OMEGA_M0=Om0
        astCalc.OMEGA_L0=1.0-Om0
        try:
            cosmo_model=FlatLambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0, Tcmb0 = 2.72548)
            self.mf.update(cosmo_model = cosmo_model, sigma_8 = sigma_8)
        except:
            raise Exception("failed to update mf when H0 = %.3f Om0 = %.3f Ob0 = %.3f sigma_8 = %.3f" % (H0, Om0, Ob0, sigma_8))
        self._doClusterCount()
        
        # Stuff to enable us to draw mock samples (see drawSample)
        # Interpolators here need to be updated each time we change cosmology
        if self.enableDrawSample == True:
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
    
    
    def drawSample(self, SNRLimit, areaMask, RMSMap, wcs, scalingRelationDict, tckQFitDict, 
                   photFilterLabel = None, extName = None):
        """Draw a cluster sample from the mass function, generate mock y0~ values by applying the given 
        scaling relation parameters, and then apply the survey selection function using the given
        SNRLimit.
        
        scalingRelationDict should contain keys 'tenToA0', 'B0', 'Mpivot', 'sigma_int' (this is the
        contents of massOptions in nemo .yml config files).
        
        Most likely you should set SNRLimit to thresholdSigma, as given in the nemo .yml config file
        (the resulting catalog could be cut in z, fixed_SNR afterwards anyway).
        
        photFilterLabel and extName are only used for adding a 'template' key to each object: useful
        for quickly checking which tile (extName) an object is in.
        
        This routine is used in the nemoMock script.

        Returns a catalog (list containing object dictionaries).
                
        """
        
        tenToA0, B0, Mpivot, sigma_int=[scalingRelationDict['tenToA0'], scalingRelationDict['B0'], 
                                        scalingRelationDict['Mpivot'], scalingRelationDict['sigma_int']]
            
        numClusters=int(round(self.numClusters))
        
        # Draw coords (assuming clusters aren't clustered - which they are...)
        ysInMask, xsInMask=np.where(areaMask != 0)
        coordIndices=np.random.randint(0, len(xsInMask), numClusters)
        ys=ysInMask[coordIndices]
        xs=xsInMask[coordIndices]
        RADecCoords=wcs.pix2wcs(xs, ys)
        RADecCoords=np.array(RADecCoords)
        RAs=RADecCoords[:, 0]
        decs=RADecCoords[:, 1]

        # Roll for redshifts
        zs=interpolate.splev(np.random.uniform(0, 1, numClusters), self.tck_zRoller)
                
        # Draw from mass function - takes ~6 sec for E-D56
        log10Ms=[]
        t0=time.time()
        for z in zs:
            # Draw M|z
            zIndex=np.where(abs(self.z-z) == abs(self.z-z).min())[0][0]
            MRoll=np.random.uniform(0, 1)
            log10M=interpolate.splev(MRoll, self.tck_log10MRoller[zIndex])
            log10Ms.append(float(log10M))
        log10Ms=np.array(log10Ms)
        t1=time.time()
        
        # Mock observations - takes < 1 min for E-D56
        # NOTE: we need to re-jig where fRel is applied generally for multi-freq analysis
        mockCatalog=[]
        for x, y, RADeg, decDeg, z, log10M in zip(xs, ys, RAs, decs, zs, log10Ms):
            # NOTE: Cosmology changes from astLib defaults handled in __init__
            # We're still using astCalc here as it's several times faster for da than cosmo_model above
            Ez=astCalc.Ez(z)
            Hz=Ez*astCalc.H0
            G=4.301e-9  # in MSun-1 km2 s-2 Mpc
            criticalDensity=(3*np.power(Hz, 2))/(8*np.pi*G)
            M500=np.power(10, log10M)
            R500Mpc=np.power((3*M500)/(4*np.pi*500*criticalDensity), 1.0/3.0)                     
            theta500Arcmin=np.degrees(np.arctan(R500Mpc/astCalc.da(z)))*60.0
            Q=interpolate.splev(theta500Arcmin, tckQFitDict[extName])
            fRel=simsTools.calcFRel(z, M500)
            true_y0=tenToA0*np.power(Ez, 2)*np.power(M500/Mpivot, 1+B0)*Q*fRel   
            if true_y0 < 0:
                continue    # This happens if extremely low mass / wander out of range where Q fit is valid
            # Mock "observations" (apply intrinsic scatter and noise)...
            scattered_y0=np.exp(np.random.normal(np.log(true_y0), sigma_int))        
            y0Noise=RMSMap[y, x]
            measured_y0=np.random.normal(scattered_y0, y0Noise)
            if measured_y0 > y0Noise*SNRLimit and y0Noise > 0:
                name=catalogTools.makeACTName(RADeg, decDeg, prefix = 'MOCK-CL')
                objDict={'name': name, 
                         'RADeg': RADeg,
                         'decDeg': decDeg,
                         'true_M500': M500/1e14,
                         'true_fixed_y_c': true_y0/1e-4,
                         'fixed_y_c': measured_y0/1e-4,
                         'fixed_err_y_c': y0Noise/1e-4,
                         'fixed_SNR': measured_y0/y0Noise,
                         'redshift': z,
                         'redshiftErr': 0}
                # Usefull to quickly check which tile something is in
                if photFilterLabel != None and extName != None:
                    objDict['template']=photFilterLabel+"#"+extName
                mockCatalog.append(objDict)

        return mockCatalog
        
