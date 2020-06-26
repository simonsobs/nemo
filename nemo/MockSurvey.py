"""

This module defines the MockSurvey class, used for mass function calculations, obtaining de-biased cluster
mass estimates, selection function calculations, and generating mock catalogs.

"""

import os
import sys
import numpy as np
#import IPython
import astropy.table as atpy
import pylab as plt
import subprocess
from astropy.cosmology import FlatLambdaCDM
import pyccl as ccl
#from colossus.cosmology import cosmology
#from colossus.lss import mass_function
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
    """An object that provides routines calculating cluster counts (using `CCL <https://github.com/LSSTDESC/CCL>`_) 
    and generating mock catalogs for a given set of cosmological and mass scaling relation parameters.
    The Tinker et al. (2008) halo mass function is used (hardcoded at present, but in principle this can
    easily be swapped for any halo mass function supported by CCL).
    
    Attributes:
        areaDeg2 (:obj:`float`): Survey area in square degrees.
        zBinEdges (:obj:`np.ndarray`): Defines the redshift bins for the cluster counts.
        z (:obj:`np.ndarray`): Centers of the redshift bins.
        log10M (:obj:`np.ndarray`): Centers of the log10 mass bins for the cluster counts
            (in MSun, with mass defined according to `delta` and `rhoType`).
        a (:obj:`np.ndarray`): Scale factor (1/(1+z)).
        delta (:obj:``float`): Overdensity parameter, used for mass definition (e.g., 200, 500).
        rhoType (:obj:`str`): Density definition, either 'matter' or 'critical', used for mass definition.
        mdef (:obj:`pyccl.halos.massdef.MassDef`): CCL mass definition object, defined by `delta` and `rhoType`.
        transferFunction (:obj:`str`): Transfer function to use, as understood by CCL (e.g., 'eisenstein_hu', 
            'boltzmann_camb').
        H0 (:obj:`float`): The Hubble constant at redshift 0, in km/s/Mpc.
        Om0 (:obj:`float`): Dimensionless total (dark + baryonic) matter density parameter at redshift 0.
        Ob0 (:obj:`float`): Dimensionless baryon matter density parameter at redshift 0.
        sigma8 (:obj:`float`): Defines the amplitude of the matter power spectrum.
        ns (:obj:`float`): Scalar spectral index of matter power spectrum.
        volumeMpc3 (:obj:`float`): Co-moving volume in Mpc3 for the given survey area and cosmological
            parameters.
        numberDensity (:obj:`np.ndarray`): Number density of clusters (per cubic Mpc) on the 
            (z, log10M) grid.
        clusterCount (:obj:`np.ndarray`): Cluster counts on the (z, log10M) grid.
        numClusters (:obj:`float`): Total number of clusters in the survey area above the minimum mass
            limit.
        numClustersByRedshift (:obj:`np.ndarray`): Number of clusters in the survey area above the
            minimum mass limit, as a function of redshift.
    
    """
    def __init__(self, minMass, areaDeg2, zMin, zMax, H0, Om0, Ob0, sigma8, ns, zStep = 0.01, 
                 enableDrawSample = False, delta = 500, rhoType = 'critical', 
                 transferFunction = 'boltzmann_camb'):
        """Create a MockSurvey object, for performing calculations of cluster counts or generating mock
        catalogs. The Tinker et al. (2008) halo mass function is used (hardcoded at present, but in 
        principle this can easily be swapped for any halo mass function supported by CCL).
        
        Args:
            minMass (:obj:`float`): The minimum mass, in MSun. This should be set considerably lower than
                the actual survey completeness limit, otherwise completeness calculations will be wrong.
            areaDeg2 (:obj:`float`): Specifies the survey area in square degrees, which scales the
                resulting cluster counts accordingly.
            zMin (:obj:`float`): Minimum redshift for the (z, log10M) grid.
            zMax (:obj:`float`): Maximum redshift for the (z, log10M) grid.
            H0 (:obj:`float`): The Hubble constant at redshift 0, in km/s/Mpc.
            Om0 (:obj:`float`): Dimensionless total (dark + baryonic) matter density parameter at redshift 0.
            Ob0 (:obj:`float`): Dimensionless baryon matter density parameter at redshift 0.
            sigma8 (:obj:`float`): Defines the amplitude of the matter power spectrum.
            ns (:obj:`float`): Scalar spectral index of matter power spectrum.  
            zStep (:obj:`float`, optional): Sets the linear spacing between redshift bins.
            enableDrawSample (:obj:`bool`, optional): This needs to be set to True to enable use of the
                :func:`self.drawSample` function. Setting this to False avoids some overhead.
            delta (:obj:``float`): Overdensity parameter, used for mass definition (e.g., 200, 500).
            rhoType (:obj:`str`): Density definition, either 'matter' or 'critical', used for mass definition.
            transferFunction (:obj:`str`): Transfer function to use, as understood by CCL (e.g., 'eisenstein_hu', 
                'boltzmann_camb').
                
        """

        zRange=np.arange(zMin, zMax+zStep, zStep)
        areaSr=np.radians(np.sqrt(areaDeg2))**2
        self.areaSr=areaSr
        self.areaDeg2=areaDeg2
        self.zBinEdges=zRange
        self.z=(zRange[:-1]+zRange[1:])/2.
        self.a=1./(1+self.z)
        
        self.delta=delta
        self.rhoType=rhoType
        self.mdef=ccl.halos.MassDef(self.delta, self.rhoType)
        self.transferFunction=transferFunction

        self.H0=-1
        self.Om0=-1
        self.Ob0=-1
        self.sigma8=-1
        self.ns=-1
        self._get_new_cosmo(H0, Om0, Ob0, sigma8, ns)

        # NOTE: This is just MSun now (NOT MSun/h)
        self.log10M=np.arange(13, 16, 0.01)
        self.M=np.power(10, self.log10M)

        self.enableDrawSample=enableDrawSample
        self.update(H0, Om0, Ob0, sigma8, ns)


    def _get_new_cosmo(self, H0, Om0, Ob0, sigma8, ns):
        if ((self.H0 != H0) or (self.Om0 != Om0) or
            (self.Ob0 != Ob0) or (self.sigma8 != sigma8)):
            self.H0=H0
            self.Om0=Om0
            self.Ob0=Ob0
            self.sigma8=sigma8
            self.ns=ns
            self.cosmoModel = ccl.Cosmology(Omega_c=Om0-Ob0,
                                            Omega_b=Ob0,
                                            h=0.01*H0,
                                            sigma8=sigma8,
                                            n_s=ns,
                                            transfer_function=self.transferFunction)
            self.mfunc = ccl.halos.MassFuncTinker08(self.cosmoModel,
                                                    self.mdef)

            
    def update(self, H0, Om0, Ob0, sigma8, ns):
        """Recalculate cluster counts for the updated cosmological parameters given.
        
        Args:
            H0 (:obj:`float`): The Hubble constant at redshift 0, in km/s/Mpc.
            Om0 (:obj:`float`): Dimensionless total (dark + baryonic) matter density parameter at redshift 0.
            Ob0 (:obj:`float`): Dimensionless baryon matter density parameter at redshift 0.
            sigma8 (:obj:`float`): Defines the amplitude of the matter power spectrum.
            ns (:obj:`float`): Scalar spectral index of matter power spectrum.  
                
        """

        self._get_new_cosmo(H0, Om0, Ob0, sigma8, ns)

        self._doClusterCount()

        # For quick Q, fRel calc (these are in MockSurvey rather than SelFn as used by drawSample)
        self.theta500Splines=[]
        self.fRelSplines=[]
        self.Ez=ccl.h_over_h0(self.cosmoModel,self.a)
        self.DAz=ccl.angular_diameter_distance(self.cosmoModel,self.a)
        self.criticalDensity=ccl.physical_constants.RHO_CRITICAL*(self.Ez*self.cosmoModel['h'])**2
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
            fitFRels=signals.calcFRel(zk, fitM500s, Ez)
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
                self.log10MRollers.append(_spline((ngtm[mask] / ngtm[0])[::-1], np.log10(self.M[mask][::-1]), k=3))
    
    def _cumulativeNumberDensity(self, z):
        """Returns N > M (per cubic Mpc), using Colossus routines.
        
        """

        h=self.cosmoModel['h']
        dndlnM=self.mfunc.get_mass_function(self.cosmoModel,
                                            self.M, 1/(1+z)) / np.log(10) #/ h**3
        dndM=dndlnM/self.M
        ngtm=integrate.cumtrapz(dndlnM[::-1], np.log(self.M), initial = 0)[::-1]
        
        MUpper=np.arange(np.log(self.M[-1]), np.log(10**18), np.log(self.M[1])-np.log(self.M[0]))
        extrapolator=_spline(np.log(self.M), np.log(dndlnM), k=1)
        MF_extr=extrapolator(MUpper)
        intUpper=integrate.simps(np.exp(MF_extr), dx=MUpper[2] - MUpper[1], even='first')
        ngtm=ngtm+intUpper
    
        return ngtm
    
    
    def _comovingVolume(self, z):
        """Returns co-moving volume in Mpc^3 (all sky) to some redshift z, using Colossus routines (taking
        care of the fact that Colossus returns all distances in Mpc/h).
        
        NOTE: Assumes flat cosmology
        
        """
        return 4.18879020479 * ccl.comoving_radial_distance(self.cosmoModel, 1./(1+z))**3

        
    def _doClusterCount(self):
        """Updates cluster count etc. after mass function object is updated.
        
        """

        zRange=self.zBinEdges
        h = self.cosmoModel['h']
        self.M=np.power(10, self.log10M) # in M_sun
        norm_mfunc=1. / np.log(10)

        # Number density by z and total cluster count (in redshift shells)
        # Can use to make P(m, z) plane
        numberDensity=[]
        clusterCount=[]
        totalVolumeMpc3=0.
        for i in range(len(zRange)-1):
            zShellMin=zRange[i]
            zShellMax=zRange[i+1]
            zShellMid=(zShellMax+zShellMin)/2.
            dndlnM=self.mfunc.get_mass_function(self.cosmoModel, self.M,
                                                1./(1+zShellMid)) * norm_mfunc
            dndM = dndlnM / self.M
            # NOTE: this differs from hmf by several % at the high-mass end (binning or interpolation?)
            n=dndM * np.gradient(self.M)
            numberDensity.append(n)
            shellVolumeMpc3=self._comovingVolume(zShellMax)-self._comovingVolume(zShellMin)
            shellVolumeMpc3=shellVolumeMpc3*(self.areaSr/(4*np.pi))
            totalVolumeMpc3+=shellVolumeMpc3
            clusterCount.append(n*shellVolumeMpc3)
        numberDensity=np.array(numberDensity)
        clusterCount=np.array(clusterCount)
        self.volumeMpc3=totalVolumeMpc3
        self.numberDensity=numberDensity
        self.clusterCount=clusterCount
        self.numClusters=np.sum(clusterCount)
        self.numClustersByRedshift=np.sum(clusterCount, axis = 1)


    def calcNumClustersExpected(self, MLimit = 1e13, zMin = 0.0, zMax = 2.0, compMz = None):
        """Calculate the number of clusters expected above a given mass limit, for the
        mass definition set when the MockSurvey object was constructed.
        
        Args:
            MLimit (:obj:`float`, optional): Mass limit above which to count clusters, in MSun.
            zMin (:obj:`float`, optional): Count clusters above this minimum redshift.
            zMax (:obj:`float`, optional): Count clusters below this maximum redshift.
            compMz (:obj:`np.ndarray`, optional): If given, a 2d array with the same dimensions
                and binning as the (z, log10M) grid, as calculated by the
                :class:`nemo.completeness.SelFn` class, that describes the completeness as
                a function of mass and redshift.
        
        Returns:
            The number of clusters in the survey area, according to the chose sample selection
            cuts.
        
        """
        
        if type(compMz) == np.ndarray:
            numClusters=compMz*self.clusterCount
        else:
            numClusters=self.clusterCount
        
        zMask=np.logical_and(np.greater(self.z, zMin), np.less(self.z, zMax))
        mMask=np.greater(self.M, MLimit)
        
        return numClusters[:, mMask][zMask].sum()
        

    def getPLog10M(self, z):
        """Returns the log10(mass) probability distribution at the given z, for the logarithmic mass
        binning and mass definition set when the MockSurvey object was constructed.
        
        Args:
            z (:obj:`float`): Redshift at which to calculate P(log10(M)).
        
        Returns:
            Array corresponding to the log10(mass) probability distribution.
        
        """
        numberDensity=self._cumulativeNumberDensity(z)
        PLog10M=numberDensity/np.trapz(numberDensity, self.M)

        return PLog10M


    def drawSample(self, y0Noise, scalingRelationDict, tckQFitDict, wcs = None, photFilterLabel = None, 
                   tileName = None, SNRLimit = None, makeNames = False, z = None, numDraws = None,
                   areaDeg2 = None, applySNRCut = False, applyPoissonScatter = True, 
                   applyIntrinsicScatter = True, applyNoiseScatter = True, seed = None):
        """Draw a cluster sample from the mass function, generating mock y0~ values (called `fixed_y_c` in
        Nemo catalogs) by applying the given scaling relation parameters, and then (optionally) applying
        a survey selection function.
        
        Args:
            y0Noise (:obj:`float` or :obj:`np.ndarray`): Either a single number (if using e.g., a survey
                average) or a noise map (2d array). A noise map must be provided here if you want the
                output catalog to contain RA, dec coordinates (in addition, a WCS object must also be
                provided - see below).
            scalingRelationDict (:obj:`dict`): A dictionary containing keys 'tenToA0', 'B0', 'Mpivot',
                'sigma_int' that describes the scaling relation between y0~ and mass (this is the
                format of `massOptions` in Nemo .yml config files).
            tckQFitDict (:obj:`dict`): A dictionary of interpolation spline knots indexed by tileName,
                that can be used to estimate Q, the filter mismatch function (see 
                :func:`nemo.signals.loadQ`).
            wcs (:obj:`astWCS.WCS`, optional): WCS object corresponding to `y0Noise`, if `y0Noise` is
                as noise map (2d image array). Needed if you want the output catalog to contain RA, dec
                coordinates.
            photFilterLabel (:obj:`str`, optional): Name of the reference filter (as defined in the
                Nemo .yml config file) that is used to define y0~ (`fixed_y_c`) and the filter mismatch 
                function, Q.
            tileName (:obj:`str`, optional): Name of the tile for which the sample will be generated.
            SNRLimit (:obj:`float`, optional): Signal-to-noise detection threshold used for the
                output catalog (corresponding to a cut on `fixed_SNR` in Nemo catalogs). Only applied
                if `applySNRCut` is also True (yes, this can be cleaned up).
            makeNames (:obj:`bool`, optional): If True, add names of the form MOCK CL JHHMM.m+/-DDMM
                to the output catalog.
            z (:obj:`float`, optional): If given produce a sample at the nearest z in the MockSurvey
                z grid. The default behaviour is to use the full redshift grid specified by `self.z`.
            numDraws (:obj:`int`, optional): If given, the number of draws to perform from the mass
                function, divided equally among the redshift bins. The default is to use the values
                contained in `self.numClustersByRedshift`.
            areaDeg2 (:obj:`float`, optional): If given, the cluster counts will be scaled to this
                area. Otherwise, they correspond to `self.areaDeg2`. This parameter will be ignored
                if `numDraws` is also given.
            applySNRCut (:obj:`bool`, optional): If True, cut the output catalog according to the
                `fixed_SNR` threshold set by `SNRLimit`.
            applyPoissonScatter (:obj:`bool`, optional): If True, add Poisson noise to the cluster
                counts (implemented by modifiying the number of draws from the mass function).
            applyIntrinsicScatter (:obj:`bool`, optional): If True, apply intrinsic scatter to the
                SZ measurements (`fixed_y_c`), as set by the `sigma_int` parameter in 
                `scalingRelationDict`.
            applyNoiseScatter (:obj:`bool`, optional): If True, apply measurement noise, generated
                from the given noise level or noise map (`y0Noise`), to the output SZ measurements
                (`fixed_y_c`).
            seed (:obj:`int`): If given, use this value for the random seed (may be useful for
                testing).
                
        Returns:
            A catalog as an :obj:`astropy.table.Table` object, in the same format as produced by
            the main `nemo` script.
        
        Notes:
            If both `applyIntrinsicScatter`, `applyNoiseScatter` are set to False, then the output
            catalog `fixed_y_c` values will be exactly the same as `true_y_c`, although each object
            will still have an error bar listed in the output catalog, corresponding to its location
            in the noise map (if given).
                
        """
        
        if seed is not None:
            np.random.seed(seed)
                
        if z is None:
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
                names.append(catalogs.makeName(RADeg, decDeg, prefix = 'MOCK-CL'))
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
                raise Exception("Negative y0 values (probably spline related) for H0 = %.6f Om0 = %.6f sigma8 = %.6f at z = %.3f" % (self.H0, self.Om0, self.sigma8, zk))
                
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
        if photFilterLabel is not None and tileName is not None:
            tab.add_column(atpy.Column([photFilterLabel]*len(tab), 'template'))
            tab.add_column(atpy.Column([tileName]*len(tab), 'tileName'))
                
        # Apply selection?
        if applySNRCut == True:
            selMask=measured_y0s > y0Noise*SNRLimit
            tab=tab[selMask]

        return tab
        
