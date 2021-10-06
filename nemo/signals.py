"""

This module contains routines for modeling cluster and source signals.

"""

import os
import sys
from pixell import enmap
import astropy
import astropy.wcs as enwcs
import astropy.io.fits as pyfits
import astropy.constants as constants
#from astropy.cosmology import FlatLambdaCDM
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
from . import plotSettings
import numpy as np
import numpy.fft as fft
import math
import pylab as plt
import pickle
import operator
import pyximport; pyximport.install()
import nemoCython
import nemo
import glob
import shutil
import yaml
import warnings
#import IPython
np.random.seed()

#------------------------------------------------------------------------------------------------------------
# Global constants (we could move others here but then need to give chunky obvious names, not just e.g. h)
TCMB=2.72548
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
    M200mDef=ccl.halos.MassDef(200, "matter", c_m_relation = 'Bhattacharya13')
    M200cDef=ccl.halos.MassDef(200, "critical", c_m_relation = 'Bhattacharya13')
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
    and the response in the second column), or can be set directly using arrays.
    
    Args:
        beamFileName(:obj:`str`, optional): Path to text file containing a beam profile in the ACT format.
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
            self.profile1d=beamData[1]
            self.rDeg=beamData[0]
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
        QFitFileName (:obj:`str`): Path to a FITS-table format file as made by :meth:`fitQ`.
        tileNames (:obj:`list`): If given, the Q-function will be defined only for these tiles (their names
            must appear in the file specified by `QFitFileName`).
    
    Attributes:
        fitDict (:obj:`dict`): Dictionary of interpolation objects, indexed by `tileName`. You should not
            need to access this directly - use :meth:`getQ` instead.
    
    """
        
    def __init__(self, QFitFileName = None, tileNames = None):
        self._zGrid=np.array([0.05, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1.0, 1.2, 1.6, 2.0])
        self._theta500ArcminGrid=np.logspace(np.log10(0.1), np.log10(55), 10)
        self.zMin=(self._zGrid).min()
        self.zMax=(self._zGrid).max()
        self.zDependent=None
        self.zDepThetaMax=None

        self.fitDict={}                
        if QFitFileName is not None:
            self.loadQ(QFitFileName, tileNames = tileNames)
        
        
    def loadQ(self, source, tileNames = None):
        """Load the filter mismatch function Q (see `Hasselfield et al. 2013 
        <https://ui.adsabs.harvard.edu/abs/2013JCAP...07..008H/abstract>`_) as a dictionary of spline fits.
        
        Args:
            source (:obj:`nemo.startUp.NemoConfig` or str): Either the path to a .fits table (containing Q fits
                for all tiles - this is normally ``selFn/QFit.fits``), or a :obj:`nemo.startUp.NemoConfig` object 
                (from which the path and tiles to use will be inferred).
            tileNames (optional, list): A list of tiles for which the Q function spline fit coefficients 
                will be extracted. If source is a :obj:`nemo.startUp.NemoConfig` object, this should be set to 
                ``None``.
        
        Returns:
            A dictionary (with tilNames as keys), containing spline knots for the Q function for each tile.
            Q values can then be obtained by using these with :func:`scipy.interpolate.splev`.
            
        """

        # Bit messy, but two modes here: 
        # - combined Q fit file for all tiles
        # - single Q fit for a single tile (interim stage, when under nemo MPI run)
        if type(source) == nemo.startUp.NemoConfig:
            tileNames=source.tileNames
            combinedQTabFileName=source.selFnDir+os.path.sep+"QFit.fits"
            loadMode=None
            if os.path.exists(combinedQTabFileName) == True:
                tileNamesInFile=[]
                with pyfits.open(combinedQTabFileName) as QTabFile:
                    for ext in QTabFile:
                        if type(ext) == astropy.io.fits.hdu.table.BinTableHDU:
                            tileNamesInFile.append(ext.name)
                tileNamesInFile.sort()
                if tileNames is None:
                    tileNames=tileNamesInFile
                loadMode="combined"
            else:
                globStr=source.selFnDir+os.path.sep+"QFit#*.fits"
                QTabFileNames=glob.glob(globStr)
                loadMode="single"
                if len(QTabFileNames) == 0:
                    raise Exception("could not find either '%s' or '%s' - needed to make QFit object" % (combinedQTabFileName, globStr))
            zMin=self._zGrid.max()
            zMax=self._zGrid.min()
            for tileName in tileNames:
                if loadMode == "combined":
                    QTab=atpy.Table().read(combinedQTabFileName, hdu = tileName)
                elif loadMode == "single":
                    QTab=atpy.Table().read(source.selFnDir+os.path.sep+"QFit#%s.fits" % (tileName))
                else:
                    raise Exception("loadMode is not defined")
                if QTab['z'].min() < zMin:
                    self.zMin=QTab['z'].min()
                if QTab['z'].max() > zMax:
                    self.zMax=QTab['z'].max()
                self.fitDict[tileName]=self._makeInterpolator(QTab)

        elif os.path.exists(source) == True:
            # Inspect file and get tile names if MEF
            if tileNames is None:
                tileNames=[]
                with pyfits.open(source) as QTab:
                    for ext in QTab:
                        if type(ext) == astropy.io.fits.hdu.table.BinTableHDU:
                            tileNames.append(ext.name)
            zMin=self._zGrid.max()
            zMax=self._zGrid.min()
            for tileName in tileNames:
                if tileName == '': # Individual, interim file name
                    assert(source.find("QFit#") > 0)
                    tileName=os.path.split(source)[-1].split("QFit#")[-1].split(".fits")[0]
                    QTab=atpy.Table().read(source)
                else:
                    QTab=atpy.Table().read(source, hdu = tileName)
                if QTab['z'].min() < zMin:
                    self.zMin=QTab['z'].min()
                if QTab['z'].max() > zMax:
                    self.zMax=QTab['z'].max()
                self.fitDict[tileName]=self._makeInterpolator(QTab)
    

    def _makeInterpolator(self, QTab):
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
            
        Returns:
            The value of *Q* (an array or a single float, depending on the input).
            
        Note:
            In the case where *Q* is a function of redshift, values outside of the range for which
            *Q* has been calculated will be filled with zeros (i.e., there is no extrapolation in
            redshift).
                    
        """
        
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
        wrtDensity=ccl.omega_x(cosmoModel, 1/(1+z), 'matter')
        #wrtDensity=cosmoModel.Om(z)*cosmoModel.critical_density(z).value
    else:
        raise Exception("wrt should be either 'critical' or 'mean'")
    #wrtDensity=(wrtDensity*np.power(Mpc_in_cm, 3))/MSun_in_g # NOTE: not needed for CCL units (MSun, Mpc etc.)
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
def makeArnaudModelProfile(z, M500, GNFWParams = 'default', cosmoModel = None):
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
    bRange=np.linspace(0, 30, 1000)
    cylPProfile=[]
    tol=1e-6
    for i in range(len(bRange)):
        b=bRange[i]
        cylPProfile.append(gnfw.integrated(b, params = GNFWParams))
        if i > 0 and abs(cylPProfile[i] - cylPProfile[i-1]) < tol:
            break
    cylPProfile=np.array(cylPProfile)
    bRange=bRange[:i+1]
    
    # Normalise to 1 at centre
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
    
    Used by ArnaudModelFilter
    
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
    
    # Adjust tol for speed vs. range of b covered
    bRange=np.linspace(0, 30, 1000)
    cylPProfile=[]
    tol=1e-6
    for i in range(len(bRange)):
        b=bRange[i]
        cylPProfile.append(gnfw.integrated(b, params = GNFWParams))
        if i > 0 and abs(cylPProfile[i] - cylPProfile[i-1]) < tol:
            break
    cylPProfile=np.array(cylPProfile)
    bRange=bRange[:i+1]
    
    # Normalise to 1 at centre
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
def makeArnaudModelSignalMap(z, M500, degreesMap, wcs, beam, GNFWParams = 'default', amplitude = None, 
                             maxSizeDeg = 15.0, convolveWithBeam = True):
    """Makes a 2d signal only map containing an Arnaud model cluster. 
    
    Args:
        z (float): Redshift; used for setting angular size.
        M500 (float): Mass within R500, defined with respect to critical density; units are solar masses.
        degreesMap (:obj:`numpy.ndarray`): A 2d array containing radial distance measured in degrees from 
            the centre of the model to be inserted. The output map will have the same dimensions and pixel
            scale (see nemoCython.makeDegreesDistanceMap).
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
        convolveWithBeam (bool, optional): If False, no beam convolution is done (it can be quicker to apply
            beam convolution over a whole source-injected map rather than per object).
    
    Returns:
        signalMap (:obj:`np.ndarray`).

    Note:
        The pixel window function is not applied here; use pixell.enmap.apply_window to do that (see 
        nemo.filters.filterMaps).    
        
    """

    # Making the 1d profile itself is the slowest part (~1 sec)
    signalDict=makeArnaudModelProfile(z, M500, GNFWParams = GNFWParams)
    tckP=signalDict['tckP']
    
    # Make cluster map (unit-normalised profile)
    rDeg=np.linspace(0.0, maxSizeDeg, 5000)
    profile1d=interpolate.splev(rDeg, tckP, ext = 1)
    if amplitude is not None:
        profile1d=profile1d*amplitude
    r2p=interpolate.interp1d(rDeg, profile1d, bounds_error=False, fill_value=0.0)
    signalMap=r2p(degreesMap)
    
    if convolveWithBeam == True:        
        signalMap=maps.convolveMapWithBeam(signalMap, wcs, beam, maxDistDegrees = maxSizeDeg)
    
    return signalMap

#------------------------------------------------------------------------------------------------------------
def makeBattagliaModelSignalMap(z, M500, degreesMap, wcs, beam, GNFWParams = 'default', amplitude = None, 
                                 maxSizeDeg = 15.0, convolveWithBeam = True):
    """Makes a 2d signal only map containing a Battaglia+2012 model cluster (taking into account the redshift
    evolution described in Table 1 and equation 11 there).
    
    Args:
        z (float): Redshift; used for setting angular size.
        M500 (float): Mass within R500, defined with respect to critical density; units are solar masses.
        degreesMap (:obj:`numpy.ndarray`): A 2d array containing radial distance measured in degrees from 
            the centre of the model to be inserted. The output map will have the same dimensions and pixel
            scale (see nemoCython.makeDegreesDistanceMap).
        GNFWParams (dict, optional): Used to specify a different profile shape to the default (which follows 
            Battaglia et al. 2012). If GNFWParams = 'default', then the default parameters as listed in 
            Battaglia et al. 2012 are used, i.e., GNFWParams = {'gamma': 0.3, 'alpha': 1.0, 'beta': 4.49,
            'c500': 1.408, 'tol': 1e-7, 'npts': 100}. Note that the definitions/sign convention is slightly
            different in Battaglia+2012 compared to Arnaud+2010 (we follow the latter). 
            Otherwise, give a dictionary that specifies the wanted values. This 
            would usually be specified using the GNFWParams key in the .yml config used when running nemo
            (see the examples/ directory).
        amplitude (float, optional): Amplitude of the cluster, i.e., the central decrement (in map units, 
            e.g., uK), or the central Comptonization parameter (dimensionless), before beam convolution. 
            Not needed for generating filter kernels.
        maxSizeDeg (float, optional): Use to limit the region over which the beam convolution is done, 
            for optimization purposes.
        convolveWithBeam (bool, optional): If False, no beam convolution is done (it can be quicker to apply
            beam convolution over a whole source-injected map rather than per object).
    
    Returns:
        signalMap (:obj:`np.ndarray`).

    Note:
        The pixel window function is not applied here; use pixell.enmap.apply_window to do that (see 
        nemo.filters.filterMaps).    
        
    """
    
    if GNFWParams == 'default':
        # NOTE: These are Table 1 values from Battaglia+2012 for M500c
        GNFWParams={'P0': 7.49, 'gamma': 0.3, 'alpha': 1.0, 'beta': 4.49, 'c500': 1.408, 'tol': 1e-7, 'npts': 100}

    # Making the 1d profile itself is the slowest part (~1 sec)
    signalDict=makeBattagliaModelProfile(z, M500, GNFWParams = GNFWParams)
    tckP=signalDict['tckP']
    
    # Make cluster map (unit-normalised profile)
    rDeg=np.linspace(0.0, maxSizeDeg, 5000)
    profile1d=interpolate.splev(rDeg, tckP, ext = 1)
    if amplitude is not None:
        profile1d=profile1d*amplitude
    r2p=interpolate.interp1d(rDeg, profile1d, bounds_error=False, fill_value=0.0)
    signalMap=r2p(degreesMap)
    
    if convolveWithBeam == True:        
        signalMap=maps.convolveMapWithBeam(signalMap, wcs, beam, maxDistDegrees = maxSizeDeg)
    
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
    results are initially cached (with a separate .fits table for each tile) under the `selFn` directory,
    before being combined into a single file at the end of a :ref:`nemoCommand` run. 
    
    The `GNFWParams` key in the `config` dictionary can be used to specify a different cluster profile shape.
    
    Args:
        config (:obj:`startUp.NemoConfig`): A NemoConfig object.
    
    Note:
        See :class:`QFit` for how to read in and use the output of this function.
        
    """
    
    t0=time.time()
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
        minTheta500Arcmin=0.1
        maxTheta500Arcmin=500.0
        numPoints=50
        theta500Arcmin_wanted=np.logspace(np.log10(minTheta500Arcmin), np.log10(maxTheta500Arcmin), numPoints)
        zRange_wanted=np.zeros(numPoints)
        zRange_wanted[np.less(theta500Arcmin_wanted, 3.0)]=2.0
        zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 3.0), np.less(theta500Arcmin_wanted, 6.0))]=1.0
        zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 6.0), np.less(theta500Arcmin_wanted, 10.0))]=0.5
        zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 10.0), np.less(theta500Arcmin_wanted, 20.0))]=0.1
        zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 20.0), np.less(theta500Arcmin_wanted, 30.0))]=0.05
        zRange_wanted[np.greater(theta500Arcmin_wanted, 30.0)]=0.01
        MRange_wanted=[]
        for theta500Arcmin, z in zip(theta500Arcmin_wanted, zRange_wanted):
            Ez=ccl.h_over_h0(cosmoModel, 1/(1+z))
            criticalDensity=ccl.physical_constants.RHO_CRITICAL*(Ez*cosmoModel['h'])**2
            R500Mpc=np.tan(np.radians(theta500Arcmin/60.0))*ccl.angular_diameter_distance(cosmoModel, 1/(1+z))
            M500=(4/3.0)*np.pi*np.power(R500Mpc, 3)*500*criticalDensity
            MRange_wanted.append(M500)         
        MRange=MRange+MRange_wanted
        zRange=zRange+zRange_wanted.tolist()
        signalMapSizeDeg=10.0
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
        signalMapSizeDeg=5.0
    else:
        raise Exception("valid values for zDepQ are 0 or 1")
            
    # Here we save the fit for each tile separately... 
    # completeness.tidyUp will put them into one file at the end of a nemo run
    for tileName in config.tileNames:
        tileQTabFileName=config.selFnDir+os.path.sep+"QFit#%s.fits" % (tileName)
        if os.path.exists(tileQTabFileName) == True:
            print("... already done Q fit for tile %s ..." % (tileName))
            continue
        print("... fitting Q in tile %s ..." % (tileName))

        # Load reference scale filter
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
                              diagnosticsDir = config.diagnosticsDir+os.path.sep+tileName)
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
        
        # A bit clunky but gets map pixel scale and shrinks map size we'll use for inserting signals
        # signalMapSizeDeg set according to lowest z model (see above), using smaller for z dependent to save RAM
        # (but then have a higher low-z cut where Q will be valid)
        extMap=np.zeros(filterObj.shape)
        wcs=filterObj.wcs
        RADeg, decDeg=wcs.getCentreWCSCoords()                
        clipDict=astImages.clipImageSectionWCS(extMap, wcs, RADeg, decDeg, signalMapSizeDeg)
        wcs=clipDict['wcs']
        extMap=clipDict['data']
        
        # Input signal maps to which we will apply filter(s)
        # We do this once and store in a dictionary for speed
        theta500ArcminDict={}
        signalMapDict={}
        signalMap=np.zeros(extMap.shape)
        degreesMap=np.ones(signalMap.shape, dtype = float)*1e6
        degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, signalMapSizeDeg)
        for z, M500MSun in zip(zRange, MRange):
            key='%.2f_%.2f' % (z, np.log10(M500MSun))
            signalMaps=[]
            fSignalMaps=[]
            y0=2e-4
            for obsFreqGHz in list(beamsDict.keys()):
                if mapDict['obsFreqGHz'] is not None:   # Normal case
                    amplitude=maps.convertToDeltaT(y0, obsFreqGHz)
                else:                                   # TILe-C case
                    amplitude=y0
                # NOTE: Q is to adjust for mismatched filter shape
                # Yes, this should have the beam in it (certainly for TILe-C)
                # NOTE: CCL can blow up for some of the extreme masses we try to feed in here
                # (so we just skip those if it happens)
                try:
                    signalMap=makeSignalModelMap(z, M500MSun, degreesMap, wcs, beamsDict[obsFreqGHz], 
                                                 amplitude = amplitude, convolveWithBeam = True, 
                                                 GNFWParams = config.parDict['GNFWParams'])
                except:
                    continue
                if realSpace == True:
                    signalMaps.append(signalMap)
                else:
                    signalMaps.append(enmap.fft(signalMap))
            signalMaps=np.array(signalMaps)
            # Skip any failed ones (see above - CCL blowing up for extreme masses)
            if len(signalMaps) == len(list(beamsDict.keys())):
                signalMapDict[key]=signalMaps
                theta500ArcminDict[key]=calcTheta500Arcmin(z, M500MSun, fiducialCosmoModel)
        
        # Filter maps with the ref kernel
        # NOTE: keep only unique values of Q, theta500Arcmin (or interpolation routines will fail)
        Q=[]
        QTheta500Arcmin=[]
        Qz=[]
        for z, M500MSun in zip(zRange, MRange):
            key='%.2f_%.2f' % (z, np.log10(M500MSun))
            if key in signalMapDict.keys():
                filteredSignal=filterObj.applyFilter(signalMapDict[key]) 
                peakFilteredSignal=filteredSignal.max()
                if peakFilteredSignal not in Q:
                    Q.append(peakFilteredSignal)      
                    QTheta500Arcmin.append(theta500ArcminDict[key])
                    Qz.append(z)
        Q=np.array(Q)
        Q=Q/Q[0]
            
        # Sort and save as FITS table (interim - all tile files gets combined at end of nemo run)
        QTab=atpy.Table()
        QTab.add_column(atpy.Column(Q, 'Q'))
        QTab.add_column(atpy.Column(QTheta500Arcmin, 'theta500Arcmin'))
        QTab.add_column(atpy.Column(Qz, 'z'))
        QTab.sort('theta500Arcmin')
        QTab.meta['NEMOVER']=nemo.__version__
        QTab.meta['ZDEPQ']=zDepQ
        QTab.write(tileQTabFileName, overwrite = True)
        
        # Test plot
        Q=QFit(tileQTabFileName)
        plotSettings.update_rcParams()
        plt.figure(figsize=(9,6.5))
        ax=plt.axes([0.12, 0.11, 0.86, 0.88])
        for z in [0.05, 0.1, 0.4, 1.0, 2.0]:
            mask=(QTab['z'] == z)
            if mask.sum() > 0:
                plt.plot(QTab['theta500Arcmin'][mask], QTab['Q'][mask], '.', label = "z = %.2f" % (z))
                thetaArr=np.logspace(np.log10(QTab['theta500Arcmin'][mask].min()), 
                                 np.log10(QTab['theta500Arcmin'][mask].max()), numPoints)
                plt.plot(thetaArr, Q.getQ(thetaArr, z, tileName = tileName), 'k-')
        plt.legend()
        plt.semilogx()
        plt.xlabel("$\\theta_{\\rm 500c}$ (arcmin)")
        plt.ylabel("$Q$ ($\\theta_{\\rm 500c}$, $z$)")
        plt.savefig(config.diagnosticsDir+os.path.sep+tileName+os.path.sep+"QFit_%s.pdf" % (tileName))
        plt.savefig(config.diagnosticsDir+os.path.sep+tileName+os.path.sep+"QFit_%s.png" % (tileName))
        plt.close()
        
        t1=time.time()
        print("... Q fit finished [tileName = %s, rank = %d, time taken = %.3f] ..." % (tileName, config.rank, t1-t0))
        

#------------------------------------------------------------------------------------------------------------
def makeCombinedQTable(config):
    """Writes dictionary of tables (containing individual tile Q fits) as a single .fits table.
    
    Returns combined Q astropy table object
    
    """

    outFileName=config.selFnDir+os.path.sep+"QFit.fits"    
    if os.path.exists(outFileName) == True:
        return atpy.Table().read(outFileName)
        
    QTabDict={}
    for tileName in config.allTileNames:
        QTabDict[tileName]=atpy.Table().read(config.selFnDir+os.path.sep+"QFit#%s.fits" % (tileName))
    
    #----
    # New - MEF
    QTabMEF=pyfits.HDUList()
    for tileName in config.allTileNames:
        with pyfits.open(config.selFnDir+os.path.sep+"QFit#%s.fits" % (tileName)) as QTab:
            QTab[1].name=tileName
            QTabMEF.append(QTab[1].copy())
    QTabMEF.writeto(outFileName, overwrite = True)
    combinedQTab=QTabMEF
    
    #----
    # Old
    #combinedQTab=atpy.Table()
    #for tabKey in list(QTabDict.keys()):
        #for colKey in QTabDict[tabKey].keys():
            #if colKey == 'theta500Arcmin':
                #if colKey not in combinedQTab.keys():
                    #combinedQTab.add_column(QTabDict[tabKey]['theta500Arcmin'], index = 0)
            #else:
                #combinedQTab.add_column(atpy.Column(QTabDict[tabKey][colKey].data, tabKey))
    #combinedQTab.meta['NEMOVER']=nemo.__version__
    #combinedQTab.write(outFileName, overwrite = True)
    
    return combinedQTab

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
            
    P=calcPMass(y0, y0Err, z, zErr, QFit, mockSurvey, tenToA0 = tenToA0, B0 = B0, Mpivot = Mpivot, 
                sigma_int = sigma_int, Ez_gamma = Ez_gamma, onePlusRedshift_power = onePlusRedshift_power, 
                applyMFDebiasCorrection = applyMFDebiasCorrection,
                applyRelativisticCorrection = applyRelativisticCorrection, fRelWeightsDict = fRelWeightsDict,
                tileName = tileName)
    
    M500, errM500Minus, errM500Plus=getM500FromP(P, mockSurvey.log10M, calcErrors = calcErrors)
    
    label=mockSurvey.mdefLabel
    
    return {'%s' % (label): M500, '%s_errPlus' % (label): errM500Plus, '%s_errMinus' % (label): errM500Minus}

#------------------------------------------------------------------------------------------------------------
def calcPMass(y0, y0Err, z, zErr, QFit, mockSurvey, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, 
              sigma_int = 0.2, Ez_gamma = 2, onePlusRedshift_power = 0.0, applyMFDebiasCorrection = True, 
              applyRelativisticCorrection = True, fRelWeightsDict = {148.0: 1.0}, return2D = False,
              tileName = None):
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
            log10M500c_zk=np.log10(mockSurvey.mdef.translate_mass(mockSurvey.cosmoModel, 
                                                                  np.power(10, log10Ms),
                                                                  1/(1+zk), mockSurvey._M500cDef))
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
        if np.less(y0pred, 0).sum() > 0:
            # This generally means we wandered out of where Q is defined (e.g., beyond mockSurvey log10M limits)
            # Or fRel can dip -ve for extreme mass at high-z (can happen with large Om0)
            raise Exception("Some predicted y0 values are negative.")
        log_y0pred=np.log(y0pred)

        Py0GivenM=np.exp(-np.power(log_y0-log_y0pred, 2)/(2*(np.power(log_y0Err, 2)+np.power(sigma_int, 2))))
        Py0GivenM=Py0GivenM/np.trapz(Py0GivenM, log10Ms)

        # Mass function de-bias
        if applyMFDebiasCorrection == True:
            PLog10M=mockSurvey.getPLog10M(zk)
            PLog10M=PLog10M/np.trapz(PLog10M, log10Ms)
        else:
            PLog10M=1.0
        
        P=Py0GivenM*PLog10M*Pz[k]
        PArr.append(P)
        
    # 2D PArr is what we would want to project onto (M, z) grid
    PArr=np.array(PArr)
        
    # Marginalised over z uncertainty
    P=np.sum(PArr, axis = 0)
    P=P/np.trapz(P, log10Ms)
    
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
        
    return P

#------------------------------------------------------------------------------------------------------------
# Mass conversion routines

# For getting x(f) - see Hu & Kravtsov
x=np.linspace(1e-3, 10, 1000)
fx=(x**3)*(np.log(1+1./x)-np.power(1+x, -1))
XF_TCK=interpolate.splrep(fx, x)
FX_TCK=interpolate.splrep(x, fx)

#------------------------------------------------------------------------------------------------------------
def gz(zIn, zMax = 1000, dz = 0.1):
    """Calculates linear growth factor at redshift z. Use Dz if you want normalised to D(z) = 1.0 at z = 0.
    
    See http://www.astronomy.ohio-state.edu/~dhw/A873/notes8.pdf for some notes on this.
    
    """
    
    zRange=np.arange(zIn, zMax, dz)
    HzPrime=[]
    for zPrime in zRange:
        HzPrime.append(astCalc.Ez(zPrime)*astCalc.H0)
    HzPrime=np.array(HzPrime)
    gz=astCalc.Ez(zIn)*np.trapz((np.gradient(zRange)*(1+zRange)) / np.power(HzPrime, 3), zRange)
    
    return gz

#------------------------------------------------------------------------------------------------------------
def calcDz(zIn):
    """Calculate linear growth factor, normalised to D(z) = 1.0 at z = 0.
    
    """
    return gz(zIn)/gz(0.0)

#------------------------------------------------------------------------------------------------------------
def criticalDensity(z):
    """Returns the critical density at the given z.
    
    """
    
    G=4.301e-9  # in MSun-1 km2 s-2 Mpc, see Robotham GAMA groups paper
    Hz=astCalc.H0*astCalc.Ez(z)
    rho_crit=((3*np.power(Hz, 2))/(8*np.pi*G))
    
    return rho_crit

#------------------------------------------------------------------------------------------------------------
def meanDensity(z):
    """Returns the mean density at the given z.
    
    """
    
    rho_mean=astCalc.OmegaMz(z)*criticalDensity(z)
    
    return rho_mean  

#------------------------------------------------------------------------------------------------------------
def M500cToMdef(M500c, z, massDef, cosmoModel):
    """Convert M500c to some other mass definition.
    
    massDef (`obj`:ccl.halos.MassDef): CCL halo mass definition
    
    """

    M500cDef=ccl.halos.MassDef(500, "critical")

    tolerance=1e-5
    scaleFactor=3.0
    ratio=1e6
    count=0
    while abs(1.0-ratio) > tolerance:
        testM500c=massDef.translate_mass(cosmoModel, scaleFactor*M500c, 1/(1+z), M500cDef)
        ratio=M500c/testM500c
        scaleFactor=scaleFactor*ratio
        count=count+1
        if count > 10:
            raise Exception("M500c -> massDef conversion didn't converge quickly enough")
        
    massX=scaleFactor*M500c
    
    return massX

#------------------------------------------------------------------------------------------------------------
def convertM200mToM500c(M200m, z):
    """Returns M500c (MSun), R500c (Mpc) for given M200m and redshift. Uses the Bhattacharya et al. c-M
    relation: http://adsabs.harvard.edu/abs/2013ApJ...766...32B

    See also Hu & Kravtsov: http://iopscience.iop.org/article/10.1086/345846/pdf
    
    """
        
    # c-M relation for full cluster sample
    Dz=calcDz(z)    # <--- this is the slow part. 3 seconds!
    nu200m=(1./Dz)*(1.12*np.power(M200m / (5e13 * np.power(astCalc.H0/100., -1)), 0.3)+0.53)
    c200m=np.power(Dz, 1.15)*9.0*np.power(nu200m, -0.29)
    
    rho_crit=criticalDensity(z)
    rho_mean=meanDensity(z)
    R200m=np.power((3*M200m)/(4*np.pi*200*rho_mean), 1/3.)

    rs=R200m/c200m
    
    f_rsOverR500c=((500*rho_crit) / (200*rho_mean)) * interpolate.splev(1./c200m, FX_TCK)
    x_rsOverR500c=interpolate.splev(f_rsOverR500c, XF_TCK)
    R500c=rs/x_rsOverR500c

    M500c=(4/3.0)*np.pi*R500c**3*(500*rho_crit)
        
    return M500c, R500c

#------------------------------------------------------------------------------------------------------------
def convertM500cToM200m(M500c, z):
    """Returns M200m given M500c
    
    """
    
    tolerance=1e-5
    scaleFactor=3.0
    ratio=1e6
    count=0
    while abs(1.0-ratio) > tolerance:
        testM500c, testR500c=convertM200mToM500c(scaleFactor*M500c, z)
        ratio=M500c/testM500c
        scaleFactor=scaleFactor*ratio
        count=count+1
        if count > 10:
            raise Exception("M500c -> M200m conversion didn't converge quickly enough")
        
    M200m=scaleFactor*M500c
    
    return M200m
