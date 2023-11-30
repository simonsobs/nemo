"""

This module defines the MockSurvey class, used for mass function calculations, obtaining de-biased cluster
mass estimates, selection function calculations, and generating mock catalogs.

"""

import os
import sys
import numpy as np
import astropy.table as atpy
import pylab as plt
import subprocess
from astropy.cosmology import FlatLambdaCDM
on_rtd=os.environ.get('READTHEDOCS', None)
if on_rtd is None:
    import pyccl as ccl
from . import signals
from . import catalogs
from . import maps
import pickle
from scipy import interpolate
from scipy import integrate
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy import stats
from astLib import *
import time

# CLASS-SZ is quite a big import (tensorflow) if we don't need it
noCLASS=os.environ.get('NEMO_NOCLASSSZ', None)
if noCLASS is None:
    try:
        from classy_sz import Class
    except:
        pass

#------------------------------------------------------------------------------------------------------------
# Global settings for CLASS-SZ, if used
CLASS_SZ_SETTINGS={# 'mass function' : 'T08M200c',
                   'hm_consistency': 0,
                   # 'concentration parameter' : 'B13',
                   'B':1.,
                   'N_ncdm' : 1,
                   'N_ur' : 2.0328,
                   'm_ncdm' : 0.0,
                   'T_ncdm' : 0.71611,
                   'z_min': 1.e-3,
                   'z_max': 3.05,
                   'redshift_epsrel': 1e-6,
                   'redshift_epsabs': 1e-100,
                   'M_min': 1e13,
                   'M_max': 1e17,
                   'mass_epsrel':1e-6,
                   'mass_epsabs':1e-100,
                   'ndim_redshifts' :5000,
                   'ndim_masses' : 500,
                   'n_m_dndlnM' : 1000,
                   'n_z_dndlnM' : 1000,
                   'HMF_prescription_NCDM': 1,
                   'no_spline_in_tinker': 1,
                   'use_m500c_in_ym_relation' : 0,
                   'use_m200c_in_ym_relation' : 1,
                   'y_m_relation' : 1,
                   'output': 'dndlnM,m500c_to_m200c,m200c_to_m500c'}

#------------------------------------------------------------------------------------------------------------
class MockSurvey(object):
    """An object that provides routines calculating cluster counts (using `CCL <https://ccl.readthedocs.io/en/latest/>`_) 
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
    def __init__(self, minMass, areaDeg2, zMin, zMax, H0, Om0, Ob0, sigma8, ns,
                 maxMass = 1e16, zStep = 0.01, numMassBins = 200, delta = 500, rhoType = 'critical',
                 transferFunction = 'boltzmann_camb', massFunction = 'Tinker08',
                 c_m_relation = 'Bhattacharya13', theoryCode = 'CCL'):
        """Create a MockSurvey object, for performing calculations of cluster counts or generating mock
        catalogs.
        
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
            maxMass (:obj:`float`, optional): The maximum mass, in MSun, used in mass function calculations.
            zStep (:obj:`float`, optional): Sets the linear spacing between redshift bins.
            numMassBins (:obj:`int`, optional): Number of mass bins to use.
            delta (:obj:``float`, optional): Overdensity parameter, used for mass definition (e.g., 200, 500).
            rhoType (:obj:`str`, optional): Density definition, either 'matter' or 'critical', used for mass definition.
            transferFunction (:obj:`str`, optional): Transfer function to use, as understood by CCL (e.g., 'eisenstein_hu',
                'boltzmann_camb').
            massFunction (:obj:`str`, optional): Name of the mass function to use, currently either 'Tinker08' or
                'Tinker10'. Mass function calculations are done by CCL.
            c_m_relation (:'obj:`str`, optional): Name of the concentration -- mass relation to assume, as understood by
                CCL (this may be used internally for conversion between mass definitions, as needed).
            theoryCode (:obj:`str`, optional): Either 'CCL' or 'CLASS-SZ'.

        """
        
        if areaDeg2 == 0:
            raise Exception("Cannot create a MockSurvey object with zero area")
        self.areaDeg2=areaDeg2
        self.areaSr=np.radians(np.sqrt(areaDeg2))**2

        zRange=np.arange(zMin, zMax+zStep, zStep)
        self.zBinEdges=zRange
        self.z=(zRange[:-1]+zRange[1:])/2.
        self.a=1./(1+self.z)
        self.numMassBins=numMassBins

        self.delta=delta
        self.rhoType=rhoType
        self.c_m_relation=c_m_relation
        self.mdef=ccl.halos.MassDef(self.delta, self.rhoType)
        self.transferFunction=transferFunction
        self.massFuncName=massFunction
        self.theoryCode=theoryCode
        if self.theoryCode not in ['CCL', 'CLASS-SZ']:
            raise Exception("theoryCode should be either 'CCL' or 'CLASS-SZ' - given %s" % (self.theoryCode))
        if self.theoryCode == 'CLASS_SZ':
            CLASS_SZ_SETTINGS['mass function']=self.massFuncName
        
        # Just for convenience when used elsewhere
        self.mdefLabel="M%d%s" % (self.delta, self.rhoType[0])
        
        self.H0=-1
        self.Om0=-1
        self.Ob0=-1
        self.sigma8=-1
        self.ns=-1
        self._get_new_cosmo(H0, Om0, Ob0, sigma8, ns)

        # NOTE: These are just MSun now (NOT MSun/h); always defined according to mdef
        self.log10M=np.linspace(np.log10(minMass), np.log10(maxMass), numMassBins)
        self.M=np.power(10, self.log10M)
        self.log10MBinEdges=np.linspace(self.log10M.min()-(self.log10M[1]-self.log10M[0])/2,
                                        self.log10M.max()+(self.log10M[1]-self.log10M[0])/2, len(self.log10M)+1)

        # Below is needed for Q calc when not using M500c definition (for now at least)
        if (self.delta == 500 and self.rhoType == 'critical') == False:
            self._M500cDef=ccl.halos.MassDef(500, "critical")
            self._transToM500c=ccl.halos.mass_translator(mass_in = self.mdef, mass_out = self._M500cDef,
                                                         concentration = self.c_m_relation)

        self.update(H0, Om0, Ob0, sigma8, ns)


    def setSurveyArea(self, areaDeg2):
        """Change the area of the survey to a user-specified value, updating the cluster
        counts accordingly.

        Args:
            areaDeg2 (:obj:`float`): Area of the survey in square degrees.

        """

        if areaDeg2 == 0:
            raise Exception("Cannot create a MockSurvey object with zero area")
        areaSr=np.radians(np.sqrt(areaDeg2))**2
        if areaDeg2 != self.areaDeg2:
            self.areaSr=areaSr
            self.areaDeg2=areaDeg2
            self._doClusterCount()


    def _get_new_cosmo(self, H0, Om0, Ob0, sigma8, ns):
        if ((self.H0 != H0) or (self.Om0 != Om0) or
            (self.Ob0 != Ob0) or (self.sigma8 != sigma8)):

            self.H0=H0
            self.Om0=Om0
            self.Ob0=Ob0
            self.sigma8=sigma8
            self.ns=ns

            # CCL cosmology - we may not want to entirely use CLASS-SZ, at least to start
            self.cosmoModel=ccl.Cosmology(Omega_c=Om0-Ob0,
                                          Omega_b=Ob0,
                                          h=0.01*H0,
                                          sigma8=sigma8,
                                          n_s=ns,
                                          transfer_function=self.transferFunction)

            if self.theoryCode == 'CCL':
                if self.massFuncName == 'Tinker10':
                    self.mfunc=ccl.halos.MassFuncTinker10(mass_def = self.mdef)
                elif self.massFuncName == 'Tinker08':
                    self.mfunc=ccl.halos.MassFuncTinker08(mass_def = self.mdef)

            elif self.theoryCode == 'CLASS-SZ':
                classDict={'Omega_b': Ob0, 'Omega_cdm':  Om0-Ob0, 'H0': H0, 'sigma8': sigma8,
                           'tau_reio':  0.0561, 'n_s': ns} # tau_reio not important but needs to be set
                self.cosmoCLASS=Class()
                self.cosmoCLASS.set(CLASS_SZ_SETTINGS) # Not sure if we want to update every time?
                self.cosmoCLASS.set(classDict)
                # self.cosmo.set(class_sz_ym_params)
                self.cosmoCLASS.compute_class_szfast()

            
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
        self.Ez2=np.power(self.Ez, 2)
        self.DAz=ccl.angular_diameter_distance(self.cosmoModel,self.a)
        self.criticalDensity=ccl.physical_constants.RHO_CRITICAL*(self.Ez*self.cosmoModel['h'])**2
        for k in range(len(self.z)):
            # NOTE: Q fit uses theta500, as does fRel (hardcoded M500 - T relation in there)
            # This bit here may not be strictly necessary, since we don't need to map on to binning
            if self.delta == 500 and self.rhoType == "critical":
                interpLim_minLog10M500c=self.log10M.min()
                interpLim_maxLog10M500c=self.log10M.max()
            else:
                interpLim_minLog10M500c=np.log10(self._transToM500c(self.cosmoModel, self.M.min(), self.a[k]))
                interpLim_maxLog10M500c=np.log10(self._transToM500c(self.cosmoModel, self.M.max(), self.a[k]))
            zk=self.z[k]
            interpPoints=100
            fitM500s=np.power(10, np.linspace(interpLim_minLog10M500c, interpLim_maxLog10M500c, interpPoints))
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


    def _cumulativeNumberDensity(self, z):
        """Returns N > M (per cubic Mpc).
        
        """

        h=self.cosmoModel['h']
        dndlnM=self.mfunc(self.cosmoModel, self.M, 1/(1+z)) / np.log(10)
        dndM=dndlnM/self.M
        ngtm=integrate.cumtrapz(dndlnM[::-1], np.log(self.M), initial = 0)[::-1]
        
        MUpper=np.arange(np.log(self.M[-1]), np.log(10**18), np.log(self.M[1])-np.log(self.M[0]))
        extrapolator=_spline(np.log(self.M), np.log(dndlnM), k=1)
        MF_extr=extrapolator(MUpper)
        intUpper=integrate.simps(np.exp(MF_extr), dx=MUpper[2] - MUpper[1], even='first')
        ngtm=ngtm+intUpper
    
        return ngtm
    
    
    def _comovingVolume(self, z):
        """Returns co-moving volume in Mpc^3 (all sky) to some redshift z.
                
        """
        return 4.18879020479 * ccl.comoving_radial_distance(self.cosmoModel, 1./(1+z))**3

        
    def _doClusterCount(self):
        """Updates cluster count etc. after mass function object is updated.
        
        """

        assert(self.areaSr == np.radians(np.sqrt(self.areaDeg2))**2)
        self.fsky=self.areaSr/(4*np.pi) # in case it was updated

        if self.theoryCode == 'CCL':
            # Previous way - this agrees with CLASS-SZ method below to 0.1 % or something
            # zRange=self.zBinEdges
            # h = self.cosmoModel['h']
            # self.M=np.power(10, self.log10M) # in M_sun
            # norm_mfunc=1. / np.log(10)
            # numberDensity=[]
            # clusterCount=[]
            # totalVolumeMpc3=0.
            # for i in range(len(zRange)-1):
            #     zShellMin=zRange[i]
            #     zShellMax=zRange[i+1]
            #     zShellMid=(zShellMax+zShellMin)/2.
            #     dndlnM=self.mfunc(self.cosmoModel, self.M, 1./(1+zShellMid)) * norm_mfunc
            #     dndM = dndlnM / self.M
            #     n=dndM * np.gradient(self.M)
            #     numberDensity.append(n)
            #     shellVolumeMpc3=self._comovingVolume(zShellMax)-self._comovingVolume(zShellMin)
            #     shellVolumeMpc3=shellVolumeMpc3*(self.areaSr/(4*np.pi))
            #     totalVolumeMpc3+=shellVolumeMpc3
            #     clusterCount.append(n*shellVolumeMpc3)
            # numberDensity=np.array(numberDensity)
            # clusterCount=np.array(clusterCount)
            # self.volumeMpc3=totalVolumeMpc3
            # self.numberDensity=numberDensity
            # self.clusterCount=clusterCount
            # self.numClusters=np.sum(clusterCount)
            # self.numClustersByRedshift=np.sum(clusterCount, axis = 1)

            # New way (same as for CLASS-SZ below)
            # Missing: self.volumeMpc3, self.numberDensity
            dndmdz=np.zeros((self.z.shape[0], self.log10M.shape[0]))
            norm_mfunc=1. / np.log(10)
            dVdzdOmega=(((signals.CLight/1000)/self.cosmoModel['H0'])*np.power(1+self.z, 2)*self.cosmoModel.angular_diameter_distance(self.a)**2)/self.cosmoModel.h_over_h0(self.a)
            for i in range(self.a.shape[0]):
                dndlnM=self.mfunc(self.cosmoModel, self.M, self.a[i]) * norm_mfunc
                dndmdz[i]=4.*np.pi*self.fsky*dVdzdOmega[i]*dndlnM
            dndmdz=dndmdz.transpose()
            self.HMFRange=np.array([np.min(dndmdz),np.max(dndmdz)])
            dndz=np.trapz(dndmdz, x = np.log(self.M), axis = 0)
            self.numClusters=np.trapz(dndz, x = self.z)
            self.numClustersByRedshift=dndz*np.gradient(self.z)
            norm=np.sum(dndmdz)/self.numClusters # Should be good to 0.1%
            self.clusterCount=dndmdz.transpose()/norm

        elif self.theoryCode == 'CLASS-SZ':
            # Still need to implement (maybe): self.volumeMpc3, self.numberDensity
            lnms=np.log(np.power(10,  self.log10M)*self.cosmoCLASS.h())
            dndmdz = np.zeros((self.log10M.shape[0], self.z.shape[0]))
            for (im,mm) in enumerate(lnms):
                dndmdz[im,:]=4.*np.pi*self.fsky*np.vectorize(self.cosmoCLASS.get_volume_dVdzdOmega_at_z)(self.z)*np.vectorize(self.cosmoCLASS.get_dndlnM_at_z_and_M)(self.z,np.exp(mm))
            self.HMFRange=np.array([np.min(dndmdz),np.max(dndmdz)])
            dndz=np.trapz(dndmdz,x = lnms,axis = 0)
            self.numClusters=np.trapz(dndz, x = self.z)
            self.numClustersByRedshift=dndz*np.gradient(self.z)
            norm=np.sum(dndmdz)/self.numClusters # Should be good to 0.1%
            self.clusterCount=dndmdz.transpose()/norm


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


    def drawSample(self, y0Noise, scalingRelationDict, QFit = None, wcs = None, photFilterLabel = None,\
                   tileName = None, SNRLimit = None, makeNames = False, z = None, numDraws = None,\
                   areaDeg2 = None, applySNRCut = False, applyPoissonScatter = True,\
                   applyIntrinsicScatter = True, applyNoiseScatter = True,\
                   applyRelativisticCorrection = True, verbose = False, biasModel = None):
        """Draw a cluster sample from the mass function, generating mock y0~ values (called `fixed_y_c` in
        Nemo catalogs) by applying the given scaling relation parameters, and then (optionally) applying
        a survey selection function.
        
        Args:
            y0Noise (:obj:`float` or :obj:`np.ndarray`): Either a single number (if using e.g., a survey
                average), an RMS table (with columns 'areaDeg2' and 'y0RMS'), or a noise map (2d array).
                A noise map must be provided here if you want the output catalog to contain RA, dec
                coordinates (in addition, a WCS object must also be provided - see below).
            scalingRelationDict (:obj:`dict`): A dictionary containing keys 'tenToA0', 'B0', 'Mpivot',
                'sigma_int' that describes the scaling relation between y0~ and mass (this is the
                format of `massOptions` in Nemo .yml config files).
            QFit (:obj:`nemo.signals.QFit`, optional): Object that handles the filter mismatch
                function, *Q*. If not given, the output catalog will not contain `fixed_y_c` columns,
                only `true_y_c` columns.
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
            applyRelativisticCorrection (:obj:`bool`, optional): If True, apply the relativistic
                correction.
                
        Returns:
            A catalog as an :obj:`astropy.table.Table` object, in the same format as produced by
            the main `nemo` script.
        
        Notes:
            If both `applyIntrinsicScatter`, `applyNoiseScatter` are set to False, then the output
            catalog `fixed_y_c` values will be exactly the same as `true_y_c`, although each object
            will still have an error bar listed in the output catalog, corresponding to its location
            in the noise map (if given).
                
        """

        t0=time.time()
        if z is None:
            zRange=self.z
        else:
            # Pick the nearest z on the grid
            zIndex=np.argmin(abs(z-self.z))
            zRange=[self.z[zIndex]]

        # Add Poisson noise
        if applyPoissonScatter == True:
            numClusters=np.random.poisson(lam = self.numClusters)
        else:
            numClusters=int(round(self.numClusters))

        if areaDeg2 is not None:
            numClusters=int(round(numClusters*(areaDeg2/self.areaDeg2)))

        if numDraws is not None:
            numClusters=numDraws

        tenToA0, B0, Mpivot, sigma_int=[scalingRelationDict['tenToA0'], scalingRelationDict['B0'],
                                        scalingRelationDict['Mpivot'], scalingRelationDict['sigma_int']]

        # If given y0Noise as RMSMap, draw coords (assuming clusters aren't clustered - which they are...)
        # NOTE: switched to using valid part of RMSMap here rather than areaMask - we need to fix the latter to same area
        # It isn't a significant issue though
        if type(y0Noise) == np.ndarray and y0Noise.ndim == 2:
            # This generates even density RA, dec coords on the whole sky taking into account the projection
            # Consequently, this is inefficient if fed individual tiles rather than a full sky noise map
            assert(wcs is not None)
            RMSMap=y0Noise
            xsList=[]
            ysList=[]
            maxCount=10000
            count=0
            while(len(xsList) < numClusters):
                count=count+1
                if count > maxCount:
                    raise Exception("Failed to generate enough random coords in %d iterations" % (maxCount))
                theta=np.degrees(np.pi*2*np.random.uniform(0, 1, numClusters))
                phi=np.degrees(np.arccos(2*np.random.uniform(0, 1, numClusters)-1))-90
                xyCoords=np.array(wcs.wcs2pix(theta, phi))
                xs=np.array(np.round(xyCoords[:, 0]), dtype = int)
                ys=np.array(np.round(xyCoords[:, 1]), dtype = int)
                mask=np.logical_and(np.logical_and(xs >= 0, xs < RMSMap.shape[1]), np.logical_and(ys >= 0, ys < RMSMap.shape[0]))
                xs=xs[mask]
                ys=ys[mask]
                mask=RMSMap[ys, xs] > 0
                xsList=xsList+xs[mask].tolist()
                ysList=ysList+ys[mask].tolist()
            xs=np.array(xsList)[:numClusters]
            ys=np.array(ysList)[:numClusters]
            del xsList, ysList
            RADecCoords=wcs.pix2wcs(xs, ys)
            RADecCoords=np.array(RADecCoords)
            RAs=RADecCoords[:, 0]
            decs=RADecCoords[:, 1]
            y0Noise=RMSMap[ys, xs]
        elif type(y0Noise) == atpy.Table:
            noisetck=interpolate.splrep(np.cumsum(y0Noise['areaDeg2']/y0Noise['areaDeg2'].sum()), y0Noise['y0RMS'], k = 1)
            rnd=np.random.uniform(0, 1, numClusters)
            vals=interpolate.splev(rnd, noisetck, ext = 3)
            if np.any(vals < 0) or np.any(vals == np.nan):
                raise Exception("Failed to make interpolating spline for RMSTab in tileName = %s" % (tileName))
            y0Noise=vals
            RAs=np.zeros(numClusters)
            decs=np.zeros(numClusters)
        else:
            y0Noise=np.ones(numClusters)*y0Noise
            RAs=np.zeros(numClusters)
            decs=np.zeros(numClusters)

        # Fancy names or not?
        if makeNames == True:
            names=[]
            for RADeg, decDeg in zip(RAs, decs):
                names.append(catalogs.makeName(RADeg, decDeg, prefix = 'MOCK-CL'))
        else:
            names=np.arange(numClusters)+1

        # Draw zs and masses using rejection sampling (from Inigo Zubeldia)
        # NOTE: We can't parallelize this with e.g. multiprocessing because we can't pickle swig objects from pyccl
        zs=np.zeros(y0Noise.shape)
        zErrs=np.zeros(y0Noise.shape)
        log10Ms=np.zeros(y0Noise.shape)
        log10M500cs=np.zeros(y0Noise.shape)
        Qs=np.zeros(y0Noise.shape)
        fRels=np.zeros(y0Noise.shape)
        t0=time.time()
        for clusterIndex in range(numClusters):
            log10Ms[clusterIndex], zs[clusterIndex], log10M500cs[clusterIndex], Qs[clusterIndex], fRels[clusterIndex]=self._drawSampleRow(QFit = QFit)
        t1=time.time()

        # For some cosmo parameters, fRel can wander outside its range for crazy masses
        # So we just cap it at 0.1 here just to avoid -ve in log
        fRels[fRels <= 0]=0.1
        fRels[fRels > 1]=1.0
        Ez2s=np.power(ccl.h_over_h0(self.cosmoModel, 1/(1+zs)), 2)
        try:
            true_y0s=tenToA0*Ez2s*np.power(np.power(10, log10Ms)/Mpivot, 1+B0)*Qs
        except:
            raise Exception("Negative y0 values (probably spline related) for H0 = %.6f Om0 = %.6f sigma8 = %.6f at z = %.3f" % (self.H0, self.Om0, self.sigma8, zk))
        if applyRelativisticCorrection == True:
            true_y0s=true_y0s*fRels

        # Add noise and intrinsic scatter everywhere
        if applyIntrinsicScatter == True:
            scattered_y0s=np.exp(np.random.normal(np.log(true_y0s), sigma_int, len(true_y0s)))
        else:
            scattered_y0s=true_y0s
        if applyNoiseScatter == True:
            measured_y0s=np.random.normal(scattered_y0s, y0Noise)
        else:
            measured_y0s=scattered_y0s

        # NOTE: We're now allowing user to specify mass definition rather than hardcoding M500c
        # So, label the output true mass column appropriately
        massColLabel="true_M%d%s" % (self.delta, self.rhoType[0])
        tab=atpy.Table()
        tab.add_column(atpy.Column(names, 'name'))
        tab.add_column(atpy.Column(RAs, 'RADeg'))
        tab.add_column(atpy.Column(decs, 'decDeg'))
        tab.add_column(atpy.Column(np.power(10, log10Ms)/1e14, massColLabel))
        if 'true_M500c' not in tab.keys():
            tab.add_column(atpy.Column(np.power(10, log10M500cs)/1e14, 'true_M500c'))
        if QFit is None:
            tab.add_column(atpy.Column(true_y0s/1e-4, 'true_y_c'))
        else:
            tab.add_column(atpy.Column(Qs, 'true_Q'))
            tab.add_column(atpy.Column(true_y0s/1e-4, 'true_fixed_y_c'))
            tab.add_column(atpy.Column(measured_y0s/1e-4, 'fixed_y_c'))
            tab.add_column(atpy.Column(y0Noise/1e-4, 'fixed_err_y_c'))
            tab['true_fixed_SNR']=tab['true_fixed_y_c']/tab['fixed_err_y_c']  # True truth, but pre-intrinsic and measurement scatter
            # tab['true_fixed_SNR']=(scattered_y0s/1e-4)/tab['fixed_err_y_c']     # With intrinsic scatter, no measurement scatter
            # tab['true_fixed_SNR']=tab['fixed_y_c']/tab['fixed_err_y_c']         # Like forced photometry case on a real map at true location
            # Apply optimization bias first, then it'll feed through to SNR automatically
            if biasModel is not None:
                corrFactors=biasModel['func'](tab['true_fixed_SNR'], biasModel['params'][0], biasModel['params'][1], biasModel['params'][2])
                tab['fixed_y_c']=tab['fixed_y_c']*corrFactors
            tab['fixed_SNR']=tab['fixed_y_c']/tab['fixed_err_y_c']

        tab.add_column(atpy.Column(zs, 'redshift'))
        tab.add_column(atpy.Column(zErrs, 'redshiftErr'))
        if photFilterLabel is not None and tileName is not None:
            tab.add_column(atpy.Column([photFilterLabel]*len(tab), 'template'))
            tab.add_column(atpy.Column([tileName]*len(tab), 'tileName'))

        # Apply selection?
        if applySNRCut == True:
            selMask=tab['fixed_SNR'] > tab['fixed_err_y_c']*SNRLimit
            tab=tab[selMask]
        t1=time.time()

        return tab

    def _drawSampleRow(self, QFit = None):
        hmfEval=0.
        hmfSample=1.
        while hmfSample > hmfEval:
            ln10MSample=np.random.rand()*(self.log10M[-1]-self.log10M[0])+self.log10M[0]
            zSample=np.random.rand()*(self.z[-1]-self.z[0])+self.z[0]
            hmfSample=np.random.rand()*(self.HMFRange[-1]-self.HMFRange[0])+self.HMFRange[0]
            if self.theoryCode == 'CCL':
                aSample=1/(1+zSample)
                norm_mfunc=1. / np.log(10)
                dVdzdOmega=(((signals.CLight/1000)/self.cosmoModel['H0'])*np.power(1+zSample, 2)*self.cosmoModel.angular_diameter_distance(aSample)**2)/self.cosmoModel.h_over_h0(aSample)
                dndlnM=self.mfunc(self.cosmoModel, np.power(10, ln10MSample), aSample) * norm_mfunc
                hmfEval=4.*np.pi*self.fsky*dVdzdOmega*dndlnM
            elif self.theoryCode == 'CLASS-SZ':
                lnMhSample=np.log(10**ln10MSample*self.cosmoCLASS.h())
                hmfEval=4.*np.pi*self.fsky*self.cosmoCLASS.get_volume_dVdzdOmega_at_z(zSample)*self.cosmoCLASS.get_dndlnM_at_z_and_M(zSample,np.exp(lnMhSample))
        log10M=ln10MSample
        z=zSample
        # Add other quantities - where these are interpolated, they are on a z grid
        zIndex=np.argmin(abs(z-self.z))
        zk=self.z[zIndex]
        if self.delta == 500 and self.rhoType == "critical":
            log10M500c=log10M
        else:
            log10M500c=np.log10(self._transToM500c(self.cosmoModel, np.power(10, log10M), 1/(1+zk)))
        theta500=interpolate.splev(log10M500c, self.theta500Splines[zIndex], ext = 3)
        if QFit is not None:
            Q=QFit.getQ(theta500, z = zk, tileName = tileName)
        else:
            Q=1.0
        fRel=interpolate.splev(log10M500c, self.fRelSplines[zIndex], ext = 3)
        return log10M, z, log10M500c, Q, fRel
