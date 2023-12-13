"""

This module contains routines for modeling cluster and source signals.

"""

import os
import sys
from pixell import enmap, curvedsky, utils, pointsrcs
import astropy
import astropy.wcs as enwcs
import astropy.io.fits as pyfits
import astropy.constants as constants
from astLib import *
from scipy import ndimage
from scipy import interpolate
from scipy import stats
import time
import astropy.table as atpy
import nemo
from . import maps
from . import catalogs
from . import photometry
from . import filters
from . import gnfw
from . import completeness
from . import plotSettings
import numpy as np
import numpy.fft as fft
import math
import pylab as plt
import pickle
import operator
import nemo
import glob
import shutil
import yaml
import warnings
np.random.seed()

#------------------------------------------------------------------------------------------------------------
# Global constants (we could move others here but then need to give chunky obvious names, not just e.g. h)
TCMB=2.72548
CLight=299792458 # m/s
Mpc_in_cm=constants.pc.value*100*1e6
MSun_in_g=constants.M_sun.value*1000

# Default cosmology (e.g., for fitQ)
#fiducialCosmoModel=FlatLambdaCDM(H0 = 70.0, Om0 = 0.3, Ob0 = 0.05, Tcmb0 = TCMB)

# Default cosmology (e.g., for fitQ) - now based on CCL rather than astropy
Om0=0.3
Ob0=0.05
H0=70
sigma8=0.8
ns=0.95
transferFunction="boltzmann_camb"
on_rtd=os.environ.get('READTHEDOCS', None)
if on_rtd is None:
    import pyccl as ccl
    fiducialCosmoModel=ccl.Cosmology(Omega_c=Om0-Ob0, Omega_b=Ob0, h=0.01*H0, sigma8=sigma8, n_s=ns,
                                    transfer_function=transferFunction)

    # For CCL-based mass conversions
    M200mDef=ccl.halos.MassDef(200, "matter")
    M200cDef=ccl.halos.MassDef(200, "critical")
    M500cDef=ccl.halos.MassDef(500, "critical")
else:
    fiducialCosmoModel=None
    M200mDef=None
    M200cDef=None
    M500cDef=None

#------------------------------------------------------------------------------------------------------------
class BeamProfile(object):
    """Describes the beam profile (i.e., the point spread function for some instrument in real space). This
    can be either read from a white-space delimited text file (with the angle in degrees in the first column
    and the response in the second column; or as a beam transform file with *l* in the first column, and
    *B*\ :sub:`l` in the second column), or can be set directly using arrays.
    
    Args:
        beamFileName(:obj:`str`, optional): Path to text file containing a beam profile (or transform) in
            the ACT format.
        profile1d (:obj:`np.ndarray`, optional): One dimensional beam profile, with index 0 at the centre.
        rDeg (:obj:`np.ndarray`, optional): Corresponding angular distance in degrees from the centre for
            the beam profile.

    Attributes:
        profile1d (:obj:`np.ndarray`): One dimensional beam profile, with index 0 at the centre.
        rDeg (:obj:`np.ndarray`): Corresponding angular distance in degrees from the centre for the 
            beam profile.
        tck (:obj:`tuple`): Spline knots for interpolating the beam onto different angular binning 
            (in degrees), for use with :meth:`scipy.interpolate.splev`.
        FWHMArcmin (float): Estimate of the beam FWHM in arcmin.

    """
    
    def __init__(self, beamFileName = None, profile1d = None, rDeg = None):
        
        if beamFileName is not None:
            beamData=np.loadtxt(beamFileName).transpose()
            # Identify if beam file is a profile or a transform
            if beamData[0][1]-beamData[0][0] >= 1:
                ell=beamData[0]
                Bell=beamData[1]
                if len(np.unique(np.diff(ell))) != 1:
                    raise Exception("If using a beam transform file, need delta ell = 1 between all ell values.")
                Bell=Bell/Bell[0]
                rDeg=np.linspace(0.0, 0.5, 1800) # This may need adjusting
                prof=curvedsky.harm2profile(Bell, np.radians(rDeg))
                prof=prof/prof[0]
                self.profile1d=prof
                self.rDeg=rDeg
                self.Bell=Bell
                self.ell=ell
            else:
                self.profile1d=beamData[1]
                self.rDeg=beamData[0]
                self.Bell=curvedsky.profile2harm(self.profile1d, np.radians(self.rDeg))
                self.Bell=self.Bell/self.Bell[0]
                self.ell=np.arange(len(self.Bell))
                #prof=curvedsky.harm2profile(Bell, np.radians(self.rDeg))
        else:
            self.profile1d=profile1d
            self.rDeg=rDeg

        if self.profile1d is not None and self.rDeg is not None:
            self.tck=interpolate.splrep(self.rDeg, self.profile1d)
        
        # This is really just for sorting a list of beams by resolution
        self.FWHMArcmin=self.rDeg[np.argmin(abs(self.profile1d-0.5))]*60*2

#------------------------------------------------------------------------------------------------------------
class QFit(object):
    """A class for managing the filter mismatch function, referred to as `Q` in the ACT papers from
    `Hasselfield et al. (2013) <http://adsabs.harvard.edu/abs/2013JCAP...07..008H>`_ onwards.
    
    Args:
        QSource (:obj:`str`, optional): The source to use for Q (the filter mismatch function) - either
            'fit' (to use results from the original Q-fitting routine), 'injection' (to use Q derived
            from source injection simulations), or 'hybrid' (to use 'fit' at scales less than the
            reference filter scale, and 'injection' at scales greater than the reference filter scale).
            For the 'injection' and 'hybrid' methods, `selFnDir` must be supplied.
        selFnDir (:obj:`str`, optional): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
        QFitFileName (`str`, optional): Path to a FITStable containing Q fits for all tiles - this is
            normally ``selFn/QFit.fits``). This is only used if QSource is set to ``fit``.
        tileNames (:obj:`list`): If given, the Q-function will be defined only for these tiles.
    
    Attributes:
        fitDict (:obj:`dict`): Dictionary of interpolation objects, indexed by `tileName`. You should not
            need to access this directly - use :meth:`getQ` instead.
    
    """
        
    def __init__(self, QSource = 'fit', selFnDir = None, QFitFileName = None, tileNames = None):
        self._zGrid=np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0])
        self._theta500ArcminGrid=np.logspace(np.log10(0.1), np.log10(55), 10)
        self.zMin=(self._zGrid).min()
        self.zMax=(self._zGrid).max()
        self.zDependent=None
        self.zDepThetaMax=None
        self.selFnDir=selFnDir

        self.fitDict={}

        self.QSource=QSource
        if self.QSource not in ['fit', 'injection', 'hybrid']:
            raise Exception("QSource must be either 'fit', 'injection', or 'hybrid'")

        if self.QSource == 'fit' or self.QSource == 'hybrid':
            if self.selFnDir is not None and QFitFileName is None:
                self.loadQ(self.selFnDir+os.path.sep+"QFit.fits", tileNames = tileNames)
            if QFitFileName is not None:
                self.loadQ(QFitFileName, tileNames = tileNames)

        elif self.QSource == 'injection':
            theta500s, thetaQ=self._loadInjectionData()
            self.fitDict[None]=interpolate.InterpolatedUnivariateSpline(theta500s, thetaQ, ext = 1)


    def _loadInjectionData(self):
        # Stuff from the source injection sims (now required for completeness calculation)
        if self.selFnDir is None:
            raise Exception("selFnDir must be supplied when using 'injection' or 'hybrid' QSource")
        injDataPath=self.selFnDir+os.path.sep+"sourceInjectionData.fits"
        inputDataPath=self.selFnDir+os.path.sep+"sourceInjectionInputCatalog.fits"
        injTab=atpy.Table().read(injDataPath)
        inputTab=atpy.Table().read(inputDataPath)
        SNRCut=5.0
        theta500s, binCentres, compThetaGrid, thetaQ=completeness._parseSourceInjectionData(injTab, inputTab, SNRCut)

        return theta500s, thetaQ

        
    def loadQ(self, QFitFileName, tileNames = None):
        """Load the filter mismatch function Q (see `Hasselfield et al. 2013 
        <https://ui.adsabs.harvard.edu/abs/2013JCAP...07..008H/abstract>`_) as a dictionary of spline fits.
        
        Args:
            QFitFileName(:obj:`str`): The path to a .fits table (containing Q fits for all tiles - this is
                normally ``selFn/QFit.fits``).
            tileNames (optional, list): A list of tiles for which the Q function spline fit coefficients 
                will be extracted. If source is a :obj:`nemo.startUp.NemoConfig` object, this should be set to 
                ``None``.

        Returns:
            None

        """

        if self.QSource == 'hybrid':
            injThetas, injQs=self._loadInjectionData()
            refTheta=None

        # Inspect file and get tile names if MEF
        if tileNames is None:
            tileNames=[]
            with pyfits.open(QFitFileName) as QTab:
                for ext in QTab:
                    if type(ext) == astropy.io.fits.hdu.table.BinTableHDU:
                        tileNames.append(ext.name)
        zMin=self._zGrid.max()
        zMax=self._zGrid.min()
        QStack=[]
        thetaStack=[]
        for tileName in tileNames:
            QTab=atpy.Table().read(QFitFileName, hdu = tileName)
            if QTab['z'].min() < zMin:
                self.zMin=QTab['z'].min()
            if QTab['z'].max() > zMax:
                self.zMax=QTab['z'].max()
            if self.QSource == 'hybrid':
                # Splice on the injection Q estimate at scales larger than reference scale
                # NOTE: Doing things this way to ensure using same thetas across all tiles
                # (so we can average easily below)
                if refTheta is None:
                    refTheta=QTab['theta500Arcmin'][QTab['Q'] > 1].min()
                fitQs=QTab['Q'][QTab['theta500Arcmin'] <= refTheta]
                fitThetas=QTab['theta500Arcmin'][QTab['theta500Arcmin'] <= refTheta]
                hybQTab=atpy.Table()
                hybQTab['theta500Arcmin']=list(fitThetas)+list(injThetas[injThetas > refTheta])
                hybQTab['Q']=list(fitQs)+list(injQs[injThetas > refTheta])
                hybQTab.meta=QTab.meta
                QTab=hybQTab
            QStack.append(QTab['Q'].data)
            thetaStack.append(QTab['theta500Arcmin'].data)
            self.fitDict[tileName]=self._makeInterpolatorFromQTab(QTab)

        # Average Q across all tiles
        medQ=np.median(np.array(QStack), axis = 0)
        medTheta=np.median(np.array(thetaStack), axis = 0)
        medQTab=atpy.Table()
        medQTab['Q']=medQ
        medQTab['theta500Arcmin']=QTab['theta500Arcmin']
        if 'z' in QTab.keys():
            medQTab['z']=QTab['z']
        medQTab.meta=QTab.meta
        self.fitDict[None]=self._makeInterpolatorFromQTab(medQTab)


    def _makeInterpolatorFromQTab(self, QTab):
        """Inspects QTab, and makes an interpolator object - 2d if there is z-dependence, 1d if not.
        
        """
        
        if QTab.meta['ZDEPQ'] == 0:
            QTab.sort('theta500Arcmin')
            spline=interpolate.InterpolatedUnivariateSpline(QTab['theta500Arcmin'], QTab['Q'], ext = 1)
            if self.zDependent == True:
                raise Exception("QFit contains a mixture of z-dependent and z-independent tables")
            self.zDepThetaMax=None
            self.zDependent=False
        elif QTab.meta['ZDEPQ'] == 1:
            with warnings.catch_warnings():
                warnings.simplefilter('ignore', UserWarning)
                spline=interpolate.LSQBivariateSpline(QTab['z'], QTab['theta500Arcmin'], QTab['Q'],
                                                      self._zGrid, self._theta500ArcminGrid)
            zs=np.unique(QTab['z'])
            thetaMaxs=[]
            for z in zs:
                thetaMaxs.append(QTab['theta500Arcmin'][QTab['z'] == z].max())
            self.zDepThetaMax=interpolate.InterpolatedUnivariateSpline(zs, thetaMaxs)
            if self.zDependent == False:
                raise Exception("QFit contains a mixture of z-dependent and z-independent tables")
            self.zDependent=True
        else:
            raise Exception("Valid ZDEPQ values are 0 or 1 only")
        
        return spline
        
        
    def getQ(self, theta500Arcmin, z = None, tileName = None):
        """Return the value of Q (the filter mismatch function) using interpolation.
        
        Args:
            theta500Arcmin (:obj:`float` or :obj:`np.ndarray`): The angular scale at which *Q* will
                be calculated. This can be an array or a single value.
            z (:obj:`float`, optional): Redshift, only used if *Q* is a function of
                redshift, otherwise it is ignored. This must be a single value only, 
                i.e., not an array.
            tileName (:obj:`str`, optional): The name of the tile to use for the *Q* function.
                If None, or if the given tileName is not found, then an average over all tiles
                is used.
            
        Returns:
            The value of *Q* (an array or a single float, depending on the input).
            
        Note:
            In the case where *Q* is a function of redshift, values outside of the range for which
            *Q* has been calculated will be filled with zeros (i.e., there is no extrapolation in
            redshift).
                    
        """

        if tileName not in self.fitDict.keys():
            tileName=None
        
        if z is not None:
            if type(z) == np.ndarray and z.shape == (1,):
                z=float(z)
            if type(z) is not float and type(z) is not np.float64:
                raise Exception("z must be a float, and not, e.g., an array")
        
        if self.zDependent == True:
            Qs=self.fitDict[tileName](z, theta500Arcmin)[0]
            thetaMask=theta500Arcmin > self.zDepThetaMax(z)
            Qs[thetaMask]=0.0
            if z < self.zMin or z > self.zMax:
                Qs=0
        else:
            # Univariate case handles own valid bounds checking
            Qs=self.fitDict[tileName](theta500Arcmin)

        if (Qs < 0).sum() > 0:
            #print("WARNING: negative Q value in tileName = %s" % (tileName))
            Qs[Qs < 0]=0
        
        return Qs

#------------------------------------------------------------------------------------------------------------
def fSZ(obsFrequencyGHz, TCMBAlpha = 0.0, z = None):
    """Returns the frequency dependence of the (non-relativistic) Sunyaev-Zel'dovich effect.
    
    Args:
        obsFrequencyGHz (float): Frequency in GHz at which to calculate fSZ.
        TCMBAlpha (float, optional): This should always be zero unless you really do want to make a model
            where CMB temperature evolves T0*(1+z)^{1-TCMBAlpha}.
        z (float, optional): Redshift - needed only if TCMBAlpha is non-zero.
    
    Returns:
        Value of SZ spectral shape at given frequency (neglecting relativistic corrections).
        
    """

    h=constants.h.value
    kB=constants.k_B.value
    sigmaT=constants.sigma_T.value
    me=constants.m_e.value
    c=constants.c.value
    x=(h*obsFrequencyGHz*1e9)/(kB*TCMB)
    if TCMBAlpha != 0 and z is not None:
        assert(z >= 0)
        x=x*np.power(1+z, TCMBAlpha)
    fSZ=x*((np.exp(x)+1)/(np.exp(x)-1))-4.0
    
    return fSZ

#------------------------------------------------------------------------------------------------------------
def calcRDeltaMpc(z, MDelta, cosmoModel, delta = 500, wrt = 'critical'):
    """Calculate RDelta (e.g., R500c, R200m etc.) in Mpc, for a halo with the given mass and redshift.
    
    Args:
        z (float): Redshift.
        MDelta (float): Halo mass in units of solar masses, using the definition set by `delta` and `wrt`.
        cosmoModel (:obj:`pyccl.Cosmology`): Cosmology object.
        delta (float, optional): Overdensity (e.g., typically 500 or 200).
        wrt (str, optional): Use 'critical' or 'mean' to set the definition of density with respect to the
            critical density or mean density at the given redshift.
    
    Returns:
        RDelta (in Mpc)
    
    """

    if type(MDelta) == str:
        raise Exception("MDelta is a string - use, e.g., 1.0e+14 (not 1e14 or 1e+14)")

    Ez=ccl.h_over_h0(cosmoModel, 1/(1+z))
    if wrt == 'critical':
        wrtDensity=ccl.physical_constants.RHO_CRITICAL*(Ez*cosmoModel['h'])**2
    elif wrt == 'mean':
        wrtDensity=ccl.omega_x(cosmoModel, 1/(1+z), 'matter')*ccl.physical_constants.RHO_CRITICAL*(Ez*cosmoModel['h'])**2
    else:
        raise Exception("wrt should be either 'critical' or 'mean'")
    RDeltaMpc=np.power((3*MDelta)/(4*np.pi*delta*wrtDensity), 1.0/3.0)
        
    return RDeltaMpc

#------------------------------------------------------------------------------------------------------------
def calcR500Mpc(z, M500c, cosmoModel):
    """Calculate R500 (in Mpc), with respect to critical density.

    Args:
        z (float): Redshift.
        M500c (float): Mass within R500c (i.e., with respect to critical density) in units of solar masses.
        cosmoModel (`:obj:`pyccl.Cosmology`): Cosmology object.
    
    Returns:
        R500c (in Mpc)
    
    """
    
    R500Mpc=calcRDeltaMpc(z, M500c, cosmoModel, delta = 500, wrt = 'critical')
    
    return R500Mpc

#------------------------------------------------------------------------------------------------------------
def calcTheta500Arcmin(z, M500, cosmoModel):
    """Given `z`, `M500` (in MSun), returns the angular size equivalent to R:sub:`500c`, with respect to the
    critical density.
    
    Args:
        z (float): Redshift.
        M500 (float): Mass within R500c (i.e., with respect to critical density) in units of solar masses.
        cosmoModel (`:obj:`pyccl.Cosmology`): Cosmology object.
    
    Returns:
        theta500c (in arcmin)
    
    """
    
    R500Mpc=calcR500Mpc(z, M500, cosmoModel)
    #theta500Arcmin=np.degrees(np.arctan(R500Mpc/cosmoModel.angular_diameter_distance(z).value))*60.0
    theta500Arcmin=np.degrees(np.arctan(R500Mpc/ccl.angular_diameter_distance(cosmoModel, 1/(1+z))))*60.0
    
    return theta500Arcmin
    
#------------------------------------------------------------------------------------------------------------
def makeArnaudModelProfile(z, M500, GNFWParams = 'default', cosmoModel = None, binning = 'log'):
    """Given z, M500 (in MSun), returns dictionary containing Arnaud model profile (well, knots from spline 
    fit, 'tckP' - assumes you want to interpolate onto an array with units of degrees) and parameters 
    (particularly 'y0', 'theta500Arcmin').
    
    Use GNFWParams to specify a different shape. If GNFWParams = 'default', then the default parameters as listed
    in gnfw.py are used, i.e., 
    
    GNFWParams = {'P0': 8.403, 'c500': 1.177, 'gamma': 0.3081, 'alpha': 1.0510, 'beta':  5.4905, 'tol': 1e-7,
                  'npts': 100}
    
    Otherwise, give a dictionary that specifies the wanted values. This would usually be specified as
    GNFWParams in the filter params in the nemo .par file (see the example .par files).
    
    If cosmoModel is None, use default (Om0, Ol0, H0) = (0.3, 0.7, 70 km/s/Mpc) cosmology.
    
    Used by ArnaudModelFilter
    
    """

    if cosmoModel is None:
        cosmoModel=fiducialCosmoModel

    if GNFWParams == 'default':
        GNFWParams=gnfw._default_params
    
    # Adjust tol for speed vs. range of b covered
    if binning == 'linear': # Old
        bRange=np.linspace(0, 30, 1000)
    elif binning == 'log':
        bRange=np.logspace(np.log10(1e-6), np.log10(100), 300)
    else:
        raise Exception("'binning' must be 'linear' or 'log' (given '%s')." % (binning))

    # Older, much slower code - can now be removed
    # cylPProfile=[]
    # tol=1e-6
    # for i in range(len(bRange)):
    #     b=bRange[i]
    #     cylPProfile.append(gnfw.integrated(b, params = GNFWParams))
    # cylPProfile=np.array(cylPProfile)
    # cylPProfile_orig=cylPProfile/cylPProfile.max()

    # Much faster, based on routines in pixell but for our A10-style GNFW function
    cylPProfile=gnfw.tsz_profile_los(bRange, c = GNFWParams['c500'], alpha = GNFWParams['alpha'],
                                     beta = GNFWParams['beta'], gamma = GNFWParams['gamma'])
    cylPProfile=cylPProfile/cylPProfile.max()

    # Calculate R500Mpc, theta500Arcmin corresponding to given mass and redshift
    theta500Arcmin=calcTheta500Arcmin(z, M500, cosmoModel)
    
    # Map between b and angular coordinates
    # NOTE: c500 now taken into account in gnfw.py
    thetaDegRange=bRange*(theta500Arcmin/60.)
    tckP=interpolate.splrep(thetaDegRange, cylPProfile)
    
    return {'tckP': tckP, 'theta500Arcmin': theta500Arcmin, 'rDeg': thetaDegRange}

#------------------------------------------------------------------------------------------------------------
def makeBattagliaModelProfile(z, M500c, GNFWParams = 'default', cosmoModel = None):
    """Given z, M500 (in MSun), returns dictionary containing Battaglia+2012 model profile (well, knots from
    spline fit, 'tckP' - assumes you want to interpolate onto an array with units of degrees) and parameters 
    (particularly 'y0', 'theta500Arcmin').
    
    Use GNFWParams to specify a different shape. If GNFWParams = 'default', then the default parameters as
    listed in Battaglia et al. 2012 are used, i.e., GNFWParams = {'gamma': 0.3, 'alpha': 1.0, 'beta': 4.49,
    'c500': 1.408, 'tol': 1e-7, 'npts': 100}. Note that the definitions/sign convention is slightly
    different in Battaglia+2012 compared to Arnaud+2010 (we follow the latter). 
    
    Otherwise, give a dictionary that specifies the wanted values. This would usually be specified as
    GNFWParams in the filter params in the nemo .par file (see the example .par files).
    
    If cosmoModel is None, use default (Om0, Ol0, H0) = (0.3, 0.7, 70 km/s/Mpc) cosmology.
    
    Used by BattagliaModelFilter
    
    """

    if cosmoModel is None:
        cosmoModel=fiducialCosmoModel
    
    if GNFWParams == 'default':
        # NOTE: These are Table 1 values from Battaglia+2012 for M500c
        GNFWParams={'P0': 7.49, 'gamma': 0.3, 'alpha': 1.0, 'beta': 4.49, 'c500': 1.408, 'tol': 1e-7, 'npts': 100}
    
    # Redshift dependence
    # (we do P0 here anyway but since we have arbitrary normalization that seems pointless)
    # These are all defined for M200c in Battaglia+2012
    # Parameters for shape are for M500c in Table 1 of Battaglia+2012
    # NOTE: Some transforming between A10 <-> B12 conventions here
    P0=GNFWParams['P0']
    P0_alpha_m=0.226
    P0_alpha_z=-0.957
    xc=1/GNFWParams['c500']
    xc_alpha_m=-0.0833
    xc_alpha_z=0.853
    beta=GNFWParams['beta']-0.3
    beta_alpha_m=0.0480
    beta_alpha_z=0.615

    M200c=M500cToMdef(M500c, z, M200cDef, cosmoModel)

    P0z=P0*np.power(M200c/1e14, P0_alpha_m)*np.power(1+z, P0_alpha_z)
    xcz=xc*np.power(M200c/1e14, xc_alpha_m)*np.power(1+z, xc_alpha_z)
    betaz=beta*np.power(M200c/1e14, beta_alpha_m)*np.power(1+z, beta_alpha_z)
    
    # Some more B12 -> A10 notation conversion
    GNFWParams['P0']=P0z
    GNFWParams['beta']=betaz+0.3
    GNFWParams['c500']=1/xcz
    GNFWParams['gamma']=0.3
    GNFWParams['alpha']=1.0    
    
    # Slower, original code
    # bRange=np.logspace(np.log10(1e-6), np.log10(100), 300)
    # cylPProfile=[]
    # tol=1e-6
    # for i in range(len(bRange)):
    #     b=bRange[i]
    #     cylPProfile.append(gnfw.integrated(b, params = GNFWParams))
    #     if i > 0 and abs(cylPProfile[i] - cylPProfile[i-1]) < tol:
    #         break
    # cylPProfile=np.array(cylPProfile)
    # bRange=bRange[:i+1]
    # cylPProfile=cylPProfile/cylPProfile.max()

    # Much faster, based on routines in pixell but for our A10-style GNFW function
    bRange=np.logspace(np.log10(1e-6), np.log10(100), 300)
    cylPProfile=gnfw.tsz_profile_los(bRange, c = GNFWParams['c500'], alpha = GNFWParams['alpha'],
                                     beta = GNFWParams['beta'], gamma = GNFWParams['gamma'])
    cylPProfile=cylPProfile/cylPProfile.max()

    # Calculate R500Mpc, theta500Arcmin corresponding to given mass and redshift
    theta500Arcmin=calcTheta500Arcmin(z, M500c, cosmoModel)
    
    # Map between b and angular coordinates
    # NOTE: c500 now taken into account in gnfw.py
    thetaDegRange=bRange*(theta500Arcmin/60.)
    tckP=interpolate.splrep(thetaDegRange, cylPProfile)
    
    return {'tckP': tckP, 'theta500Arcmin': theta500Arcmin, 'rDeg': thetaDegRange}


#------------------------------------------------------------------------------------------------------------
def makeBeamModelSignalMap(degreesMap, wcs, beam, amplitude = None):
    """Makes a 2d signal only map containing the given beam.
    
    Args:
        degreesMap (:obj:`np.ndarray`): Map of angular distance from the object position.
        wcs (:obj:`astWCS.WCS`): WCS corresponding to degreesMap.
        beam (:obj:`BeamProfile` or str): Either a BeamProfile object, or a string that gives the path to a 
            text file that describes the beam profile.
        amplitude (:obj: float, optional): Specifies the amplitude of the input signal (in map units, 
            e.g., uK), before beam convolution. This is only needed if this routine is being used to inject
            sources into maps. It is not needed for making filter kernels.
        
    Returns:
        signalMap (:obj:`np.ndarray`)

    Note:
        The pixel window function is not applied here; use pixell.enmap.apply_window to do that (see 
        nemo.filters.filterMaps).
        
    """
    
    if amplitude is None:
        amplitude=1.0
    
    if type(beam) == str:
        beam=BeamProfile(beamFileName = beam)        
    profile1d=amplitude*beam.profile1d

    # Turn 1d profile into 2d
    r2p=interpolate.interp1d(beam.rDeg, profile1d, bounds_error=False, fill_value=0.0)
    signalMap=r2p(degreesMap)
    
    return signalMap

#------------------------------------------------------------------------------------------------------------
def _paintSignalMap(shape, wcs, tckP, beam = None, RADeg = None, decDeg = None, amplitude = None,
                    maxSizeDeg = 10.0, convolveWithBeam = True, vmin = 1e-12, omap = None,
                    obsFrequencyGHz = None, TCMBAlpha = 0, z = None):
    """Use Sigurd's fast object painter to paint given signal into map.

    Notes:
        Setting vmin arbitrarily small seems the best way to avoid unwanted behavior (e.g., cutting off
        clusters before they are painted in properly).

    """

    # Now allowing arrays of sky coords and amplitudes (but one profile shape)
    if RADeg is None and decDeg is None:
        RADeg, decDeg=wcs.getCentreWCSCoords()
    dtype=np.float32 # only float32 supported by fast srcsim
    amp=1.0
    if convolveWithBeam == True:
        if beam is None:
            raise Exception("No beam supplied.")
        if type(beam) == str:
            beam=BeamProfile(beamFileName = beam)
        rht=utils.RadialFourierTransform()
        rprof=interpolate.splev(np.degrees(rht.r), tckP, ext = 1)
        lbeam=np.interp(rht.l, beam.ell, beam.Bell)
        lprof=rht.real2harm(rprof)
        lprof*=lbeam
        rprof=rht.harm2real(lprof)
        r, rprof=rht.unpad(rht.r, rprof)
    else:
        rDeg=np.logspace(np.log10(1e-6), np.log10(maxSizeDeg), 5000)
        rprof=interpolate.splev(rDeg, tckP, ext = 1)
        r=np.radians(rDeg)
    if amplitude is not None:
        amp=rprof[0]*amplitude
        rprof=rprof/rprof[0]

    if type(RADeg) == np.ndarray and type(decDeg) == np.ndarray:
        poss=np.array([np.radians(decDeg), np.radians(RADeg)]).astype(dtype)
    else:
        poss=np.array([[np.radians(decDeg)], [np.radians(RADeg)]]).astype(dtype)
    if type(amp) == np.ndarray:
        amps=np.array(amp, dtype = dtype)
    else:
        amps=np.array([amp], dtype = dtype)

    if omap is not None:
        ra1=RADeg-(maxSizeDeg/2)/np.cos(np.radians(decDeg))
        ra0=RADeg+(maxSizeDeg/2)/np.cos(np.radians(decDeg))
        dec1=decDeg-maxSizeDeg/2
        dec0=decDeg+maxSizeDeg/2
        clip=astImages.clipUsingRADecCoords(omap, wcs, ra1, ra0, dec0, dec1)
        modelClip=pointsrcs.sim_objects(clip['data'].shape, clip['wcs'].AWCS, poss, amps, (r, abs(rprof)), vmin = vmin,
                                        rmax = np.radians(maxSizeDeg), #prof_equi = False,
                                        pixwin = False)
        if obsFrequencyGHz is not None:
            modelClip=maps.convertToDeltaT(modelClip, obsFrequencyGHz = obsFrequencyGHz,
                                           TCMBAlpha = TCMBAlpha, z = z)
        xMin, xMax, yMin, yMax=clip['clippedSection']
        omap[yMin:yMax, xMin:xMax]=omap[yMin:yMax, xMin:xMax]+modelClip
        signalMap=omap
    else:
        signalMap=pointsrcs.sim_objects(shape, wcs.AWCS, poss, amps, (r, abs(rprof)), vmin = vmin,
                                        rmax = np.radians(maxSizeDeg), #prof_equi = False,
                                        pixwin = False)
        if obsFrequencyGHz is not None:
            signalMap=maps.convertToDeltaT(signalMap, obsFrequencyGHz = obsFrequencyGHz,
                                           TCMBAlpha = TCMBAlpha, z = z)

    if rprof[0] < 0:
        signalMap=signalMap*-1

    return signalMap

#------------------------------------------------------------------------------------------------------------
def makeArnaudModelSignalMap(z, M500, shape, wcs, beam = None, RADeg = None, decDeg = None,\
                             GNFWParams = 'default', amplitude = None, maxSizeDeg = 15.0,\
                             convolveWithBeam = True, cosmoModel = None, painter = 'pixell',
                             omap = None, obsFrequencyGHz = None, TCMBAlpha = 0):
    """Makes a 2d signal only map containing an Arnaud model cluster.
    
    Args:
        z (float): Redshift; used for setting angular size.
        M500 (float): Mass within R500, defined with respect to critical density; units are solar masses.
        shape (:obj:`tuple`): A tuple describing the map (numpy array) shape in pixels (height, width).
        wcs (:obj:`astWCS.WCS`): An astWCS object.
        beam (:obj:`str` or :obj:`signals.BeamProfile`): Either the file name of the text file that describes
            the beam with which the map will be convolved, or a :obj:`signals.BeamProfile` object. If None,
            no beam convolution is applied.
        RADeg (float, or array, optional): If None, the signal will be inserted at the center of the generated map.
        decDeg (float, or array, optional): If None, the signal will be inserted at the center of the generated map.
        GNFWParams (dict, optional): Used to specify a different profile shape to the default (which follows 
            Arnaud et al. 2010). If GNFWParams = 'default', then the default parameters as listed in 
            gnfw.py are used, i.e., GNFWParams = {'gamma': 0.3081, 'alpha': 1.0510, 'beta': 5.4905, 
            'tol': 1e-7, 'npts': 100}. Otherwise, give a dictionary that specifies the wanted values. This 
            would usually be specified using the GNFWParams key in the .yml config used when running nemo
            (see the examples/ directory).
        amplitude (float, or array, optional): Amplitude of the cluster, i.e., the central decrement (in map units,
            e.g., uK), or the central Comptonization parameter (dimensionless), before beam convolution. 
            Not needed for generating filter kernels.
        maxSizeDeg (float, optional): Use to limit the region over which the beam convolution is done, 
            for optimization purposes.
        convolveWithBeam (bool, optional): If False, no beam convolution is done.
    
    Returns:
        signalMap (:obj:`np.ndarray`).

    Note:
        The pixel window function is not applied here; use pixell.enmap.apply_window to do that (see 
        nemo.filters.filterMaps).    

    Note:
        If arrays of sky coordinates and amplitudes are given, multiple clusters (with the same profile) will
        be painted in the map. This is useful for running cluster injection simulations for estimating
        completeness, but only works with the 'pixell' object painter.
        
    """

    if RADeg is None or decDeg is None:
        RADeg, decDeg=wcs.getCentreWCSCoords()

    # Making the 1d profile itself is the slowest part (~1 sec)
    signalDict=makeArnaudModelProfile(z, M500, GNFWParams = GNFWParams, cosmoModel = cosmoModel)
    tckP=signalDict['tckP']

    if painter == 'legacy': # Old method
        degreesMap=np.ones(shape, dtype = np.float64)*1e6
        degreesMap, xBounds, yBounds=maps.makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg,
                                                                 maxSizeDeg)
        rDeg=np.linspace(0.0, maxSizeDeg, 5000)
        profile1d=interpolate.splev(rDeg, tckP, ext = 1)
        if amplitude is not None:
            profile1d=profile1d*amplitude
        r2p=interpolate.interp1d(rDeg, profile1d, bounds_error=False, fill_value=0.0)
        signalMap=r2p(degreesMap)
        if convolveWithBeam == True:
            signalMap=maps.convolveMapWithBeam(signalMap, wcs, beam, maxDistDegrees = maxSizeDeg)
    elif painter == 'pixell': # New method - using Sigurd's object painter
        signalMap=_paintSignalMap(shape, wcs, tckP, beam = beam, RADeg = RADeg, decDeg = decDeg,
                                  amplitude = amplitude, maxSizeDeg = maxSizeDeg,
                                  convolveWithBeam = convolveWithBeam, omap = omap,
                                  obsFrequencyGHz = obsFrequencyGHz, TCMBAlpha = TCMBAlpha, z = z)
    else:
        raise Exception("'painter' must be 'legacy' or 'pixell' (given '%s')." % (painter))

    return signalMap

#------------------------------------------------------------------------------------------------------------
def makeBattagliaModelSignalMap(z, M500, shape, wcs, beam = None, RADeg = None, decDeg = None,\
                                GNFWParams = 'default', amplitude = None, maxSizeDeg = 15.0,\
                                convolveWithBeam = True, cosmoModel = None, omap = None,\
                                obsFrequencyGHz = None, TCMBAlpha = 0):
    """Makes a 2d signal only map containing a Battaglia+2012 model cluster (taking into account the redshift
    evolution described in Table 1 and equation 11 there).

    Args:
        z (float): Redshift; used for setting angular size.
        M500 (float): Mass within R500, defined with respect to critical density; units are solar masses.
        shape (:obj:`tuple`): A tuple describing the map (numpy array) shape in pixels (height, width).
        wcs (:obj:`astWCS.WCS`): An astWCS object.
        beam (:obj:`str` or :obj:`signals.BeamProfile`): Either the file name of the text file that describes
            the beam with which the map will be convolved, or a :obj:`signals.BeamProfile` object. If None,
            no beam convolution is applied.
        RADeg: If None, the signal will be inserted at the center of the generated map.
        decDeg: If None, the signal will be inserted at the center of the generated map.
        GNFWParams (dict, optional): Used to specify a different profile shape to the default (which follows
            Arnaud et al. 2010). If GNFWParams = 'default', then the default parameters as listed in
            gnfw.py are used, i.e., GNFWParams = {'gamma': 0.3081, 'alpha': 1.0510, 'beta': 5.4905,
            'tol': 1e-7, 'npts': 100}. Otherwise, give a dictionary that specifies the wanted values. This
            would usually be specified using the GNFWParams key in the .yml config used when running nemo
            (see the examples/ directory).
        amplitude (float, optional): Amplitude of the cluster, i.e., the central decrement (in map units,
            e.g., uK), or the central Comptonization parameter (dimensionless), before beam convolution.
            Not needed for generating filter kernels.
        maxSizeDeg (float, optional): Use to limit the region over which the beam convolution is done,
            for optimization purposes.
        convolveWithBeam (bool, optional): If False, no beam convolution is done.

    Returns:
        signalMap (:obj:`np.ndarray`).

    Note:
        The pixel window function is not applied here; use pixell.enmap.apply_window to do that (see
        nemo.filters.filterMaps).

    """

    if GNFWParams == 'default':
        # NOTE: These are Table 1 values from Battaglia+2012 for M500c
        GNFWParams={'P0': 7.49, 'gamma': 0.3, 'alpha': 1.0, 'beta': 4.49, 'c500': 1.408, 'tol': 1e-7, 'npts': 100}

    #----
    # Old
    ## Making the 1d profile itself is the slowest part (~1 sec)
    #signalDict=makeBattagliaModelProfile(z, M500, GNFWParams = GNFWParams)
    #tckP=signalDict['tckP']
    
    ## Make cluster map (unit-normalised profile)
    #rDeg=np.linspace(0.0, maxSizeDeg, 5000)
    #profile1d=interpolate.splev(rDeg, tckP, ext = 1)
    #if amplitude is not None:
        #profile1d=profile1d*amplitude
    #r2p=interpolate.interp1d(rDeg, profile1d, bounds_error=False, fill_value=0.0)
    #signalMap=r2p(degreesMap)
    
    #if convolveWithBeam == True:
        #signalMap=maps.convolveMapWithBeam(signalMap, wcs, beam, maxDistDegrees = maxSizeDeg)

    # New - using Sigurd object painter
    signalDict=makeBattagliaModelProfile(z, M500, GNFWParams = GNFWParams, cosmoModel = cosmoModel)
    tckP=signalDict['tckP']
    return _paintSignalMap(shape, wcs, tckP, beam = beam, RADeg = RADeg, decDeg = decDeg,
                           amplitude = amplitude, maxSizeDeg = maxSizeDeg,
                           convolveWithBeam = convolveWithBeam, omap = omap,
                           obsFrequencyGHz = obsFrequencyGHz, TCMBAlpha = TCMBAlpha, z = z)
    
    return signalMap

#------------------------------------------------------------------------------------------------------------
def getFRelWeights(config):
    """Returns a dictionary of frequency weights used in relativistic corrections for each tile. This is 
    cached in the selFn/ dir after the first time this routine is called.
    
    """
    
    if 'photFilter' not in config.parDict.keys() or config.parDict['photFilter'] is None:
        return {}
    
    fRelWeightsFileName=config.selFnDir+os.path.sep+"fRelWeights.fits"
    if os.path.exists(fRelWeightsFileName) == False:
        fRelTab=atpy.Table()
        fRelTab.add_column(atpy.Column(config.allTileNames, 'tileName'))
        for tileCount in range(len(config.allTileNames)):
            tileName=config.allTileNames[tileCount]
            filterFileName=config.diagnosticsDir+os.path.sep+tileName+os.path.sep+"filter_%s#%s.fits" % (config.parDict['photFilter'], tileName)
            with pyfits.open(filterFileName) as img:
                for i in range(1, 10):
                    if 'RW%d_GHZ' % (i) in img[0].header.keys():
                        freqGHz=str(img[0].header['RW%d_GHZ' % (i)])
                        if freqGHz == '':
                            freqGHz='148.0'
                            print(">>> WARNING: setting freqGHz = '%s' in getFRelWeights - this is okay if you're running on a TILe-C y-map" % (freqGHz))
                        if freqGHz not in fRelTab.keys():
                            fRelTab.add_column(atpy.Column(np.zeros(len(config.allTileNames)), freqGHz))
                        fRelTab[freqGHz][tileCount]=img[0].header['RW%d' % (i)]
        fRelTab.meta['NEMOVER']=nemo.__version__
        fRelTab.write(fRelWeightsFileName, overwrite = True)
    
    return loadFRelWeights(fRelWeightsFileName)

#------------------------------------------------------------------------------------------------------------
def loadFRelWeights(fRelWeightsFileName):
    """Returns a dictionary of frequency weights used in relativistic corrections for each tile (stored in
    a .fits table, made by getFRelWeights).
    
    """

    fRelTab=atpy.Table().read(fRelWeightsFileName)
    fRelWeightsDict={}
    for row in fRelTab:
        fRelWeightsDict[row['tileName']]={}
        for key in fRelTab.keys():
            if key != 'tileName':
                fRelWeightsDict[row['tileName']][float(key)]=row[key]
    
    return fRelWeightsDict

#------------------------------------------------------------------------------------------------------------
def fitQ(config):
    """Calculates the filter mismatch function *Q* on a grid of scale sizes for each tile in the map. The
    results are combined into a single file written under the `selFn` directory.
    
    The `GNFWParams` key in the `config` dictionary can be used to specify a different cluster profile shape.
    
    Args:
        config (:obj:`startUp.NemoConfig`): A NemoConfig object.
    
    Note:
        See :class:`QFit` for how to read in and use the file produced by this function.
        
    """
    
    cosmoModel=fiducialCosmoModel
        
    # Spin through the filter kernels
    photFilterLabel=config.parDict['photFilter']
    filterList=config.parDict['mapFilters']
    for f in filterList:
        if f['label'] == photFilterLabel:
            ref=f
    
    # This could be more general... but A10 model has no z-dependence, B12 model does
    # So Q is a function of (theta500, z) for the latter
    # We add a header keyword to the QFit.fits table to indicate if z-dependence important or not
    # Everything is then handled internally by QFit class
    if ref['class'].find("Arnaud") != -1:
        makeSignalModelMap=makeArnaudModelSignalMap
        zDepQ=0
    elif ref['class'].find("Battaglia") != -1:
        makeSignalModelMap=makeBattagliaModelSignalMap
        zDepQ=1
    else:
        raise Exception("Signal model for Q calculation should either be 'Arnaud' or 'Battaglia'")
    
    # M, z -> theta ranges for Q calc - what's most efficient depends on whether there is z-dependence, or not
    # NOTE: ref filter that sets scale we compare to must ALWAYS come first
    if zDepQ == 0:
        # To safely (numerically, at least) apply Q at z ~ 0.01, we need to go to theta500 ~ 500 arcmin (< 10 deg)
        MRange=[ref['params']['M500MSun']]
        zRange=[ref['params']['z']]
        # Old - to 30 arcmin
        # theta500Arcmin_wanted=np.logspace(np.log10(0.1), np.log10(30), 50)
        # zRange_wanted=[2.0]*10 + [1.0]*10 + [0.6]*10 + [0.3]*10 + [0.1]*10
        # New - same spacing over same range, with some extra scales
        theta500Arcmin_wanted=np.power(10, np.arange(np.log10(0.1), np.log10(50), 0.05055349))
        zRange_wanted=[2.0]*10 + [1.0]*10 + [0.6]*10 + [0.3]*10 + [0.1]*10 + [0.07]*4
        MRange_wanted=[]
        for theta500Arcmin, z in zip(theta500Arcmin_wanted, zRange_wanted):
            Ez=ccl.h_over_h0(cosmoModel, 1/(1+z))
            criticalDensity=ccl.physical_constants.RHO_CRITICAL*(Ez*cosmoModel['h'])**2
            R500Mpc=np.tan(np.radians(theta500Arcmin/60.0))*ccl.angular_diameter_distance(cosmoModel, 1/(1+z))
            M500=(4/3.0)*np.pi*np.power(R500Mpc, 3)*500*criticalDensity
            MRange_wanted.append(M500)
        MRange=MRange+MRange_wanted
        zRange=zRange+zRange_wanted
        signalMapSizeDeg=15.0
        # Old
        # minTheta500Arcmin=0.1
        # maxTheta500Arcmin=500
        # numPoints=50
        # theta500Arcmin_wanted=np.logspace(np.log10(minTheta500Arcmin), np.log10(maxTheta500Arcmin), numPoints)
        # zRange_wanted=np.zeros(numPoints)
        # zRange_wanted[np.less(theta500Arcmin_wanted, 3.0)]=2.0
        # zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 3.0), np.less(theta500Arcmin_wanted, 6.0))]=1.0
        # zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 6.0), np.less(theta500Arcmin_wanted, 10.0))]=0.5
        # zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 10.0), np.less(theta500Arcmin_wanted, 20.0))]=0.1
        # zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 20.0), np.less(theta500Arcmin_wanted, 30.0))]=0.05
        # zRange_wanted[np.greater(theta500Arcmin_wanted, 30.0)]=0.01
        # MRange_wanted=[]
        # for theta500Arcmin, z in zip(theta500Arcmin_wanted, zRange_wanted):
        #     Ez=ccl.h_over_h0(cosmoModel, 1/(1+z))
        #     criticalDensity=ccl.physical_constants.RHO_CRITICAL*(Ez*cosmoModel['h'])**2
        #     R500Mpc=np.tan(np.radians(theta500Arcmin/60.0))*ccl.angular_diameter_distance(cosmoModel, 1/(1+z))
        #     M500=(4/3.0)*np.pi*np.power(R500Mpc, 3)*500*criticalDensity
        #     MRange_wanted.append(M500)
        # MRange=MRange+MRange_wanted
        # zRange=zRange+zRange_wanted.tolist()
        # signalMapSizeDeg=15.0
    elif zDepQ == 1:
        # On a z grid for evolving profile models (e.g., Battaglia et al. 2012)
        MRange=[ref['params']['M500MSun']]
        zRange=[ref['params']['z']]
        zGrid=[0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0]
        minTheta500Arcmin=0.1
        maxTheta500Arcmin=100.0  # 55' corresponds to M500c = 1e16 MSun at z = 0.05
        numPoints=24
        theta500Arcmin_wanted=np.logspace(np.log10(minTheta500Arcmin), np.log10(maxTheta500Arcmin), numPoints)
        for z in zGrid:
            MRange_wanted=[]
            for theta500Arcmin in theta500Arcmin_wanted:
                Ez=ccl.h_over_h0(cosmoModel, 1/(1+z))
                criticalDensity=ccl.physical_constants.RHO_CRITICAL*(Ez*cosmoModel['h'])**2
                R500Mpc=np.tan(np.radians(theta500Arcmin/60.0))*ccl.angular_diameter_distance(cosmoModel, 1/(1+z))
                M500=(4/3.0)*np.pi*np.power(R500Mpc, 3)*500*criticalDensity
                MRange_wanted.append(M500)
            MRange=MRange+MRange_wanted
            zRange=zRange+([z]*len(MRange_wanted))
        signalMapSizeDeg=15.0
    else:
        raise Exception("valid values for zDepQ are 0 or 1")
            
    # Here we save the fit for each tile separately... 
    QTabDict={}
    for tileName in config.tileNames:
        t0=time.time()
        print("... fitting Q in tile %s" % (tileName))

        # Load reference scale filter (it may be in memory already)
        if tileName in config.cachedFilters.keys():
            filterObj=config.cachedFilters[tileName]
        else:
            # Otherwise, can still get from disk
            foundFilt=False
            for filt in config.parDict['mapFilters']:
                if filt['label'] == config.parDict['photFilter']:
                    foundFilt=True
                    break
            if foundFilt == False:
                raise Exception("couldn't find filter that matches photFilter")
            filterClass=eval('filters.%s' % (filt['class']))
            filterObj=filterClass(filt['label'], config.unfilteredMapsDictList, filt['params'], \
                                tileName = tileName,
                                diagnosticsDir = config.diagnosticsDir)
            filterObj.loadFilter()
        
        # Real space kernel or Fourier space filter?
        if issubclass(filterObj.__class__, filters.RealSpaceMatchedFilter) == True:
            realSpace=True
        else:
            realSpace=False
        
        # Set-up the beams
        beamsDict={}
        for mapDict in config.parDict['unfilteredMaps']:
            obsFreqGHz=mapDict['obsFreqGHz']
            beamsDict[obsFreqGHz]=mapDict['beamFileName']

        # Uncomment if debugging large scale Q business
        #MRange=[MRange[-3]]
        #zRange=[zRange[-3]]

        # Actually measuring Q...
        extMap=np.zeros(filterObj.shape, dtype = np.float32)
        wcs=filterObj.wcs
        nativeNPix=extMap.shape[0]*extMap.shape[1]
        # Uncomment to pad signal maps [but then need to adjust normalization]
        #extMap=enmap.enmap(extMap, wcs = wcs.AWCS)
        #yZoom=signalMapSizeDeg/wcs.getFullSizeSkyDeg()[1]
        #xZoom=signalMapSizeDeg/wcs.getFullSizeSkyDeg()[0]
        #yPad=int(extMap.shape[0]*yZoom-extMap.shape[0])
        #xPad=int(extMap.shape[1]*xZoom-extMap.shape[1])
        #extMap=enmap.pad(extMap, (yPad, xPad))
        #h=extMap.wcs.to_header()
        #h.insert(0, ('NAXIS2', extMap.shape[0]))
        #h.insert('NAXIS2', ('NAXIS1', extMap.shape[1]))
        #wcs=astWCS.WCS(h, mode = 'pyfits')
        #paddedNPix=extMap.shape[0]*extMap.shape[1]
        # Set centre coords
        shape=extMap.shape
        RADeg, decDeg=wcs.getCentreWCSCoords()
        x, y=wcs.wcs2pix(RADeg, decDeg)

        # For testing effects of adding CMB + noise
        #simCMBDict={}
        #seed=1234
        #noiseLevel=120.0
        #for obsFreqGHz in list(beamsDict.keys()):
            #simCMBDict[obsFreqGHz]=maps.simCMBMap(shape, wcs, noiseLevel = noiseLevel, beam = beamsDict[obsFreqGHz], seed = seed)
        #print("... adding CMB and noise = %.3f in Q models" % (noiseLevel))

        # Input signal maps to which we will apply filter(s)
        # We do this once and store in a dictionary for speed
        theta500ArcminDict={}
        signalMap=np.zeros(shape, dtype = np.float32)
        Q=[]
        QTheta500Arcmin=[]
        Qz=[]
        cubeStore={} # For debugging object painting
        for obsFreqGHz in list(beamsDict.keys()):
            cubeStore[obsFreqGHz]=[]
        for z, M500MSun in zip(zRange, MRange):
            key='%.2f_%.2f' % (z, np.log10(M500MSun))
            signalMaps=[]
            fSignalMaps=[]
            y0=2e-04
            for obsFreqGHz in list(beamsDict.keys()):
                if mapDict['obsFreqGHz'] is not None:   # Normal case
                    amplitude=maps.convertToDeltaT(y0, obsFreqGHz)
                else:                                   # TILe-C case
                    amplitude=y0
                # NOTE: Q is to adjust for mismatched filter shape
                # Yes, this should have the beam in it (certainly for TILe-C)
                # NOTE: CCL can blow up for some of the extreme masses we try to feed in here
                # (so we just skip those if it happens)
                # try:
                signalMap=makeSignalModelMap(z, M500MSun, shape, wcs, beam = beamsDict[obsFreqGHz],
                                             amplitude = amplitude, convolveWithBeam = True,
                                             GNFWParams = config.parDict['GNFWParams'])
                signalMap=enmap.apply_window(signalMap, pow = 1.0)
                # except:
                #     continue
                #signalMap=signalMap+simCMBDict[obsFreqGHz]
                #signalMap=signalMap+np.random.normal(0, noiseLevel, signalMap.shape)
                #cubeStore[obsFreqGHz].append(signalMap)
                if realSpace == True:
                    signalMaps.append(signalMap)
                else:
                    signalMaps.append(enmap.fft(signalMap))
            signalMaps=np.array(signalMaps)
            # Filter maps with ref kernel
            if len(signalMaps) == len(list(beamsDict.keys())):
                filteredSignal=filterObj.applyFilter(signalMaps)
                mapInterpolator=interpolate.RectBivariateSpline(np.arange(filteredSignal.shape[0]),
                                                                np.arange(filteredSignal.shape[1]),
                                                                filteredSignal, kx = 3, ky = 3)
                peakFilteredSignal=mapInterpolator(y, x)[0][0]
                # Below is if we wanted to simulate object-finding here
                #peakFilteredSignal=filteredSignal[int(y)-10:int(y)+10, int(x)-10:int(x)+10].max() # Avoids mess at edges
                if peakFilteredSignal not in Q:
                    Q.append(peakFilteredSignal)
                    QTheta500Arcmin.append(calcTheta500Arcmin(z, M500MSun, fiducialCosmoModel))
                    Qz.append(z)
        Q=np.array(Q)
        if abs(1-Q[0]/y0) > 1e-6:
            raise Exception("Q[0]/y0 outside tolerance")
        Q=Q/y0

        # Sort and make FITS table
        QTab=atpy.Table()
        QTab.add_column(atpy.Column(Q, 'Q'))
        QTab.add_column(atpy.Column(QTheta500Arcmin, 'theta500Arcmin'))
        QTab.add_column(atpy.Column(Qz, 'z'))
        QTab.sort('theta500Arcmin')
        QTab.meta['NEMOVER']=nemo.__version__
        QTab.meta['ZDEPQ']=zDepQ
        QTab.meta['TILENAME']=tileName
        QTabDict[tileName]=QTab
        del Q, filterObj #, clipDict

    if config.MPIEnabled == True:
        QTabDictList=config.comm.gather(QTabDict, root = 0)
        if config.rank == 0:
            print("... gathered Q fits")
            combQTabDict={}
            for QTabDict in QTabDictList:
                for key in QTabDict:
                    if key not in combQTabDict:
                        combQTabDict[key]=QTabDict[key]
            QTabDict=combQTabDict

    # Write output as MEF
    if config.rank == 0:
        outFileName=config.selFnDir+os.path.sep+"QFit.fits"
        QTabMEF=pyfits.HDUList()
        for tileName in config.allTileNames:
            if tileName in QTabDict.keys():
                QTabHDU=pyfits.table_to_hdu(QTabDict[tileName])
                QTabHDU.name=tileName
                QTabMEF.append(QTabHDU)
        QTabMEF.writeto(outFileName, overwrite = True)

    if config.rank == 0:
        print("... after Q fits completed: time since start = %.3f sec" % (time.time()-config._timeStarted))

#------------------------------------------------------------------------------------------------------------
def calcWeightedFRel(z, M500, Ez, fRelWeightsDict):
    """Return fRel for the given (z, M500), weighted by frequency according to fRelWeightsDict
    
    """
    
    fRels=[]
    freqWeights=[]
    for obsFreqGHz in fRelWeightsDict.keys():
        if fRelWeightsDict[obsFreqGHz] > 0:
            fRels.append(calcFRel(z, M500, Ez, obsFreqGHz = obsFreqGHz))
            freqWeights.append(fRelWeightsDict[obsFreqGHz])
    fRel=np.average(fRels, weights = freqWeights)
    
    return fRel
    
#------------------------------------------------------------------------------------------------------------
def calcFRel(z, M500, Ez, obsFreqGHz = 148.0):
    """Calculates relativistic correction to SZ effect at specified frequency, given z, M500 in MSun.
       
    This assumes the Arnaud et al. (2005) M-T relation, and applies formulae of Itoh et al. (1998)
    
    As for H13, we return fRel = 1 + delta_SZE (see also Marriage et al. 2011)

    """
    
    # NOTE: we should define constants somewhere else...
    h=6.63e-34
    kB=1.38e-23
    sigmaT=6.6524586e-29
    me=9.11e-31
    e=1.6e-19
    c=3e8
    
    # Using Arnaud et al. (2005) M-T to get temperature
    A=3.84e14
    B=1.71
    #TkeV=5.*np.power(((cosmoModel.efunc(z)*M500)/A), 1/B)   # HMF/Astropy
    #TkeV=5.*np.power(((cosmoModel.Ez(z)*M500)/A), 1/B)   # Colossus
    TkeV=5.*np.power(((Ez*M500)/A), 1/B)
    TKelvin=TkeV*((1000*e)/kB)

    # Itoh et al. (1998) eqns. 2.25 - 2.30
    thetae=(kB*TKelvin)/(me*c**2)
    X=(h*obsFreqGHz*1e9)/(kB*TCMB)
    Xtw=X*(np.cosh(X/2.)/np.sinh(X/2.))
    Stw=X/np.sinh(X/2.)

    Y0=-4+Xtw

    Y1=-10. + (47/2.)*Xtw - (42/5.)*Xtw**2 + (7/10.)*Xtw**3 + np.power(Stw, 2)*(-(21/5.) + (7/5.)*Xtw)

    Y2=-(15/2.) +  (1023/8.)*Xtw - (868/5.)*Xtw**2 + (329/5.)*Xtw**3 - (44/5.)*Xtw**4 + (11/30.)*Xtw**5 \
        + np.power(Stw, 2)*(-(434/5.) + (658/5.)*Xtw - (242/5.)*Xtw**2 + (143/30.)*Xtw**3) \
        + np.power(Stw, 4)*(-(44/5.) + (187/60.)*Xtw)

    Y3=(15/2.) + (2505/8.)*Xtw - (7098/5.)*Xtw**2 + (14253/10.)*Xtw**3 - (18594/35.)*Xtw**4 + (12059/140.)*Xtw**5 - (128/21.)*Xtw**6 + (16/105.)*Xtw**7 \
        + np.power(Stw, 2)*(-(7098/10.) + (14253/5.)*Xtw - (102267/35.)*Xtw**2 + (156767/140.)*Xtw**3 - (1216/7.)*Xtw**4 + (64/7.)*Xtw**5) \
        + np.power(Stw, 4)*(-(18594/35.) + (205003/280.)*Xtw - (1920/7.)*Xtw**2 + (1024/35.)*Xtw**3) \
        + np.power(Stw, 6)*(-(544/21.) + (992/105.)*Xtw)

    Y4=-(135/32.) + (30375/128.)*Xtw - (62391/10.)*Xtw**2 + (614727/40.)*Xtw**3 - (124389/10.)*Xtw**4 \
        + (355703/80.)*Xtw**5 - (16568/21.)*Xtw**6 + (7516/105.)*Xtw**7 - (22/7.)*Xtw**8 + (11/210.)*Xtw**9 \
        + np.power(Stw, 2)*(-(62391/20.) + (614727/20.)*Xtw - (1368279/20.)*Xtw**2 + (4624139/80.)*Xtw**3 - (157396/7.)*Xtw**4 \
        + (30064/7.)*Xtw**5 - (2717/7.)*Xtw**6 + (2761/210.)*Xtw**7) \
        + np.power(Stw, 4)*(-(124389/10.) + (6046951/160.)*Xtw - (248520/7.)*Xtw**2 + (481024/35.)*Xtw**3 - (15972/7.)*Xtw**4 + (18689/140.)*Xtw**5) \
        + np.power(Stw, 6)*(-(70414/21.) + (465992/105.)*Xtw - (11792/7.)*Xtw**2 + (19778/105.)*Xtw**3) \
        + np.power(Stw, 8)*(-(682/7.) + (7601/210.)*Xtw)

    deltaSZE=((X**3)/(np.exp(X)-1)) * ((thetae*X*np.exp(X))/(np.exp(X)-1)) * (Y0 + Y1*thetae + Y2*thetae**2 + Y3*thetae**3 + Y4*thetae**4)

    fRel=1+deltaSZE
    
    return fRel

#------------------------------------------------------------------------------------------------------------
def getM500FromP(P, log10M, calcErrors = True):
    """Returns M500 as the maximum likelihood value from given P(log10M) distribution, together with 
    1-sigma error bars (M500, -M500Err, +M500 err).
    
    """

    # Find max likelihood and integrate to get error bars
    tckP=interpolate.splrep(log10M, P)
    fineLog10M=np.linspace(log10M.min(), log10M.max(), 10000)
    fineP=interpolate.splev(fineLog10M, tckP)
    fineP=fineP/np.trapz(fineP, fineLog10M)
    index=np.argmax(fineP)
    
    clusterLogM500=fineLog10M[index]
    clusterM500=np.power(10, clusterLogM500)/1e14

    if calcErrors == True:
        for n in range(fineP.shape[0]):
            minIndex=index-n
            maxIndex=index+n
            if minIndex < 0 or maxIndex > fineP.shape[0]:
                # This shouldn't happen; if it does, probably y0 is in the wrong units
                # Previously we threw an exception here, but we can't if using this for forced photometry
                #print("WARNING: outside M500 range - check y0 units or for problem at cluster location in map (if not in forced photometry mode)")
                clusterM500MinusErr=0.
                clusterM500PlusErr=0.
                break
            p=np.trapz(fineP[minIndex:maxIndex], fineLog10M[minIndex:maxIndex])
            if p >= 0.6827:
                clusterLogM500Min=fineLog10M[minIndex]
                clusterLogM500Max=fineLog10M[maxIndex]
                clusterM500MinusErr=(np.power(10, clusterLogM500)-np.power(10, clusterLogM500Min))/1e14
                clusterM500PlusErr=(np.power(10, clusterLogM500Max)-np.power(10, clusterLogM500))/1e14
                break        
    else:
        clusterM500MinusErr=0.
        clusterM500PlusErr=0.
    
    return clusterM500, clusterM500MinusErr, clusterM500PlusErr

#------------------------------------------------------------------------------------------------------------
def y0FromLogM500(log10M500, z, tckQFit, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, sigma_int = 0.2,
                  cosmoModel = None, applyRelativisticCorrection = True, fRelWeightsDict = {148.0: 1.0}):
    """Predict y0~ given logM500 (in MSun) and redshift. Default scaling relation parameters are A10 (as in
    H13).
    
    Use cosmoModel (:obj:`pyccl.Cosmology`) to change/specify cosmological parameters.
    
    fRelWeightsDict is used to account for the relativistic correction when y0~ has been constructed
    from multi-frequency maps. Weights should sum to 1.0; keys are observed frequency in GHz.
    
    Returns y0~, theta500Arcmin, Q
    
    NOTE: Depreciated? Nothing we have uses this.
    
    """

    if type(Mpivot) == str:
        raise Exception("Mpivot is a string - check Mpivot in your .yml config file: use, e.g., 3.0e+14 (not 3e14 or 3e+14)")
        
    # Filtering/detection was performed with a fixed fiducial cosmology... so we don't need to recalculate Q
    # We just need to recalculate theta500Arcmin and E(z) only
    M500=np.power(10, log10M500)
    theta500Arcmin=calcTheta500Arcmin(z, M500, cosmoModel)
    Q=calcQ(theta500Arcmin, tckQFit)
    
    # Relativistic correction: now a little more complicated, to account for fact y0~ maps are weighted sum
    # of individual frequency maps, and relativistic correction size varies with frequency
    if applyRelativisticCorrection == True:
        fRels=[]
        freqWeights=[]
        for obsFreqGHz in fRelWeightsDict.keys():
            fRels.append(calcFRel(z, M500, cosmoModel, obsFreqGHz = obsFreqGHz))
            freqWeights.append(fRelWeightsDict[obsFreqGHz])
        fRel=np.average(np.array(fRels), axis = 0, weights = freqWeights)
    else:
        fRel=1.0
    
    # UPP relation according to H13
    # NOTE: m in H13 is M/Mpivot
    # NOTE: this goes negative for crazy masses where the Q polynomial fit goes -ve, so ignore those
    y0pred=tenToA0*np.power(cosmoModel.efunc(z), 2)*np.power(M500/Mpivot, 1+B0)*Q*fRel
    
    return y0pred, theta500Arcmin, Q
            
#------------------------------------------------------------------------------------------------------------
def calcMass(y0, y0Err, z, zErr, QFit, mockSurvey, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, 
             sigma_int = 0.2, Ez_gamma = 2, onePlusRedshift_power = 0.0, applyMFDebiasCorrection = True, applyRelativisticCorrection = True,
             calcErrors = True, fRelWeightsDict = {148.0: 1.0}, tileName = None):
    """Returns M500 +/- errors in units of 10^14 MSun, calculated assuming a y0 - M relation (default values
    assume UPP scaling relation from Arnaud et al. 2010), taking into account the steepness of the mass
    function. The approach followed is described in H13, Section 3.2.
    
    Here, mockSurvey is a MockSurvey object. We're using this to handle the halo mass function calculations
    (in turn using the Colossus module). Supplying mockSurvey is no longer optional (and handles setting the 
    cosmology anyway when initialised or updated).
    
    tckQFit is a set of spline knots, as returned by fitQ.
    
    If applyMFDebiasCorrection == True, apply correction that accounts for steepness of mass function.
    
    If applyRelativisticCorrection == True, apply relativistic correction (weighted by frequency using the
    contents of fRelWeightsDict).
    
    If calcErrors == False, error bars are not calculated, they are just set to zero.
    
    fRelWeightsDict is used to account for the relativistic correction when y0~ has been constructed
    from multi-frequency maps. Weights should sum to 1.0; keys are observed frequency in GHz.
        
    Returns dictionary with keys M500, M500_errPlus, M500_errMinus
    
    """
    
    if y0 < 0:
        raise Exception('y0 cannot be negative')
    if y0 > 1e-2:
        raise Exception('y0 is suspiciously large - probably you need to multiply by 1e-4')
            
    P, bestQ=calcPMass(y0, y0Err, z, zErr, QFit, mockSurvey, tenToA0 = tenToA0, B0 = B0, Mpivot = Mpivot,
                       sigma_int = sigma_int, Ez_gamma = Ez_gamma, onePlusRedshift_power = onePlusRedshift_power,
                       applyMFDebiasCorrection = applyMFDebiasCorrection,
                       applyRelativisticCorrection = applyRelativisticCorrection, fRelWeightsDict = fRelWeightsDict,
                       tileName = tileName, returnQ = True)
    
    M500, errM500Minus, errM500Plus=getM500FromP(P, mockSurvey.log10M, calcErrors = calcErrors)
    
    label=mockSurvey.mdefLabel
    
    return {'%s' % (label): M500, '%s_errPlus' % (label): errM500Plus, '%s_errMinus' % (label): errM500Minus,
            'Q': bestQ}

#------------------------------------------------------------------------------------------------------------
def calcPMass(y0, y0Err, z, zErr, QFit, mockSurvey, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, 
              sigma_int = 0.2, Ez_gamma = 2, onePlusRedshift_power = 0.0, applyMFDebiasCorrection = True, 
              applyRelativisticCorrection = True, fRelWeightsDict = {148.0: 1.0}, return2D = False,
              returnQ = False, tileName = None):
    """Calculates P(M500) assuming a y0 - M relation (default values assume UPP scaling relation from Arnaud 
    et al. 2010), taking into account the steepness of the mass function. The approach followed is described 
    in H13, Section 3.2. The binning for P(M500) is set according to the given mockSurvey, as are the assumed
    cosmological parameters.
    
    This routine is used by calcMass.
    
    Ez_gamma is E(z)^gamma factor (we assumed E(z)^2 previously)
    
    onePlusRedshift_power: added multiplication by (1+z)**onePlusRedshift_power (for checking evolution)
    
    If return2D == True, returns a grid of same dimensions / binning as mockSurvey.z, mockSurvey.log10M,
    normalised such that the sum of the values is 1.
    
    """
    
    # For marginalising over photo-z errors (we assume +/-5 sigma is accurate enough)
    if zErr > 0:
        zMin=z-zErr*5
        zMax=z+zErr*5
        zMask=np.logical_and(np.greater_equal(mockSurvey.z, zMin), np.less(mockSurvey.z, zMax))
        zRange=mockSurvey.z[zMask]
        #if zMin <= 0:
            #zMin=1e-3
        #zRange=np.arange(zMin, zMax, 0.005)
        Pz=np.exp(-np.power(z-zRange, 2)/(2*(np.power(zErr, 2))))
        Pz=Pz/np.trapz(Pz, zRange)
    else:
        zRange=[z]
        Pz=np.ones(len(zRange))

    log_y0=np.log(y0)
    log_y0Err=y0Err/y0
        
    # NOTE: Swap below if want to use bigger log10M range...
    log10Ms=mockSurvey.log10M
    #log10MStep=mockSurvey.log10M[1]-mockSurvey.log10M[0]
    #log10Ms=np.arange(-100.0, 100.0, log10MStep)
            
    PArr=[]
    for k in range(len(zRange)):
        
        zk=zRange[k]
        
        # We've generalised mockSurvey to be able to use e.g. M200m, but Q defined for theta500c
        # So, need a mapping between M500c and whatever mass definition used in mockSurvey
        # This only needed for extracting Q, fRel values
        if mockSurvey.delta != 500 or mockSurvey.rhoType != "critical":
            log10M500c_zk=np.log10(mockSurvey._transToM500c(mockSurvey.cosmoModel,
                                                            np.power(10, log10Ms),
                                                            1/(1+zk)))
        else:
            log10M500c_zk=log10Ms
                
        mockSurvey_zIndex=np.argmin(abs(mockSurvey.z-zk))
        theta500s=interpolate.splev(log10M500c_zk, mockSurvey.theta500Splines[mockSurvey_zIndex], ext = 3)
        Qs=QFit.getQ(theta500s, zk, tileName = tileName)
        fRels=interpolate.splev(log10M500c_zk, mockSurvey.fRelSplines[mockSurvey_zIndex], ext = 3)   
        fRels[np.less_equal(fRels, 0)]=1e-4   # For extreme masses (> 10^16 MSun) at high-z, this can dip -ve
        y0pred=tenToA0*np.power(mockSurvey.Ez[mockSurvey_zIndex], Ez_gamma)*np.power(np.power(10, log10Ms)/Mpivot, 1+B0)*Qs
        y0pred=y0pred*np.power(1+zk, onePlusRedshift_power)
        if applyRelativisticCorrection == True:
            y0pred=y0pred*fRels

        ###
        # Adjusted below to cope with objects with big z errors
        valid=y0pred > 0
        if valid.sum() == 0:
            continue
        log_y0pred=np.log(y0pred[valid])
        Py0GivenM=np.exp(-np.power(log_y0-log_y0pred, 2)/(2*(np.power(log_y0Err, 2)+np.power(sigma_int, 2))))
        Py0GivenM=Py0GivenM/np.trapz(Py0GivenM, log10Ms[valid])

        # Mass function de-bias
        if applyMFDebiasCorrection == True:
            PLog10M=mockSurvey.getPLog10M(zk)
            PLog10M=PLog10M/np.trapz(PLog10M, log10Ms)
        else:
            PLog10M=1.0

        P=np.zeros(log10Ms.shape[0])
        if type(PLog10M) == float:
            P[valid]=Py0GivenM*PLog10M*Pz[k]
        else:
            P[valid]=Py0GivenM*PLog10M[valid]*Pz[k]
        PArr.append(P)
        ###

        # Previous version
        # if np.less(y0pred, 0).sum() > 0:
        #     pass
        #     # This generally means we wandered out of where Q is defined (e.g., beyond mockSurvey log10M limits)
        #     # Or fRel can dip -ve for extreme mass at high-z (can happen with large Om0)
        #     # Or it can be due to a photo-z with a huge error bar, and only affect some part of the z grid
        #     # So... we could write a warning here. Previously we used to trigger an exception
        #     # raise Exception("Some predicted y0 values are zero.")
        # else:
        #     log_y0pred=np.log(y0pred)
        #
        #     Py0GivenM=np.exp(-np.power(log_y0-log_y0pred, 2)/(2*(np.power(log_y0Err, 2)+np.power(sigma_int, 2))))
        #     Py0GivenM=Py0GivenM/np.trapz(Py0GivenM, log10Ms)
        #
        #     # Mass function de-bias
        #     if applyMFDebiasCorrection == True:
        #         PLog10M=mockSurvey.getPLog10M(zk)
        #         PLog10M=PLog10M/np.trapz(PLog10M, log10Ms)
        #     else:
        #         PLog10M=1.0
        #
        #     P=Py0GivenM*PLog10M*Pz[k]
        #     PArr.append(P)

    # 2D PArr is what we would want to project onto (M, z) grid
    PArr=np.array(PArr)
        
    # Marginalised over z uncertainty
    P=np.sum(PArr, axis = 0)
    P=P/np.trapz(P, log10Ms)
    # If we want Q corresponding to mass (need more work to add errors if we really want them)
    PQ=P/np.trapz(P, Qs)
    fittedQ=Qs[np.argmax(PQ)]

    # Reshape to (M, z) grid - use this if use different log10M range to mockSurvey
    #tck=interpolate.splrep(log10Ms, P)
    #P=interpolate.splev(mockSurvey.log10M, tck, ext = 1)
        
    if return2D == True:
        P2D=np.zeros(mockSurvey.clusterCount.shape)
        if zErr == 0:
            P2D[np.argmin(abs(mockSurvey.z-z))]=PArr
        else:
            P2D[zMask]=PArr
        P=P2D/P2D.sum()
        #astImages.saveFITS("test.fits", P.transpose(), None)

    if returnQ == True:
        return P, fittedQ
    else:
        return P

#------------------------------------------------------------------------------------------------------------
def MDef1ToMDef2(mass, z, MDef1, MDef2, cosmoModel, c_m_relation = 'Bhattacharya13'):
    """Convert some mass at some z defined using MDef1 into a mass defined according to MDef2.

    Args:
        mass (float): Halo mass in MSun.
        z (float): Redshift of the halo.
        MDef1 (`obj`:ccl.halos.MassDef): CCL halo mass definition you want to convert from.
        MDef2 (`obj`:ccl.halos.MassDef): CCL halo mass definition you want to convert to.
    ,   c_m_relation ('obj':`str`): Name of the concentration -- mass relation to assume, as understood by CCL.

    """

    tolerance=1e-5
    scaleFactor=3.0
    ratio=1e6
    count=0
    try:
        trans1To2=ccl.halos.mass_translator(mass_in = MDef1, mass_out = MDef2, concentration = c_m_relation)
        massX=trans1To2(cosmoModel, mass, 1/(1+z))
    except:
        trans2To1=ccl.halos.mass_translator(mass_in = MDef2, mass_out = MDef1, concentration = c_m_relation)
        while abs(1.0-ratio) > tolerance:
            testMass=trans2To1(cosmoModel, scaleFactor*mass, 1/(1+z))
            ratio=mass/testMass
            scaleFactor=scaleFactor*ratio
            count=count+1
            if count > 10:
                raise Exception("MDef1 -> MDef2 mass conversion didn't converge quickly enough")
        massX=scaleFactor*mass

    return massX

#------------------------------------------------------------------------------------------------------------
def M500cToMdef(M500c, z, massDef, cosmoModel, c_m_relation = 'Bhattacharya13'):
    """Convert M500c to some other mass definition.
    
    massDef (`obj`:ccl.halos.MassDef): CCL halo mass definition
    
    """

    return MDef1ToMDef2(M500c, z, ccl.halos.MassDef(500, "critical"), massDef, cosmoModel,
                        c_m_relation = c_m_relation)
