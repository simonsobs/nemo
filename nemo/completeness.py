"""

This module contains tools for calculating the completeness of cluster samples and handling area masks.

"""

import os
import sys
import resource
import glob
import numpy as np
import pylab as plt
import astropy.table as atpy
from astLib import *
from scipy import stats
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy import ndimage
from scipy import optimize
from scipy import integrate
from scipy import special
on_rtd=os.environ.get('READTHEDOCS', None)
if on_rtd is None:
    import pyccl as ccl
import nemo
from . import signals
from . import maps
from . import MockSurvey
from . import plotSettings
from . import startUp
from . import catalogs
from collections import OrderedDict
import colorcet
import types
import pickle
import astropy.io.fits as pyfits
import time
import shutil
import yaml
from decimal import Decimal

# If want to catch warnings as errors...
#import warnings
#warnings.filterwarnings('error')

#------------------------------------------------------------------------------------------------------------
class FootprintError(Exception):
    pass

#------------------------------------------------------------------------------------------------------------
class SelFn(object):
    """An object that describes the survey selection function. It uses the output in the ``selFn/`` directory
    (produced by the :ref:`nemoCommand` command) to calculate the survey completeness for a given
    signal-to-noise cut on a (log\ :sub:`10` mass, z) grid for a given set of cosmological and scaling
    relation parameters.
       
    Args:
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
        SNRCut (:obj:`float`): Completeness will be computed relative to this signal-to-noise selection
            cut (labelled as `fixed_SNR` in the catalogs produced by :ref:`nemoCommand`).
        configFileName (:obj:`str`, optional): Path to a Nemo configuration file. If not given, this
            will be read from ``selFnDir/config.yml``, which is the config file for the
            :ref:`nemoCommand` run that produced the ``selFn`` directory.
        footprint (:obj:`str`, optional): Use this to specify a footprint, if any are defined in the
            Nemo config used to produce the ``selFn`` dir (e.g., 'DES', 'HSC', 'KiDS' etc.). The default
            value of ``None`` uses the whole survey footprint.
        zStep (:obj:`float`, optional): Use this to set the binning in redshift for completeness
            calculations.
        tileNames (:obj:`list`, optional): If given, restrict the :class:`SelFn` object to use only these
            tiles.
        mockOversampleFactor (:obj:`float`, optional): Used only by :func:`generateMockSample`. Sets
            the oversampling level for the generated mock sample.
        downsampleRMS (:obj:`float`, optional): Downsample the resolution of the RMS (noise) tables by
            this factor. The RMS tables are generated from the noise maps, and are just a listing of noise
            level versus survey area. Downsampling speeds up completeness calculations considerably.
        applyMFDebiasCorrection (:obj:`bool`, optional): Set to `False` to disable the Eddington bias
            correction of mass estimates. Probably only useful for debugging.
        applyRelativisticCorrection (:obj:`bool`, optional): Set to `False` to disable inclusion of
            the relativistic correction in completeness calculations.
        setupAreaMask (:obj:`bool`, optional): If `True`, read in the area masks so that quick position
            checks can be done (e.g., by :meth:`SelFn.checkCoordsAreInMask`).
        enableCompletenessCalc (:obj:`bool`, optional): If `True`, set up the machinery needed to do
            completeness calculations.
        massFunction (:obj:`str`, optional): Name of the mass function to use, currently either
            'Tinker08' or 'Tinker10'. Mass function calculations are done by CCL.
        maxTheta500Arcmin (:obj:`float`, optional): If given, exclude clusters with expected angular size
            greater than this from the cluster counts (set their completeness to 0).
        method (:obj:`str`, optional): The method for calculating completeness. Options are: 'fast'
            (using a simple model based on the noise and the expected cluster signals), or 'injection'
            (directly using the results of end-to-end cluster injection and recovery sims).
        QSource (:obj:`str`, optional): The source to use for Q (the filter mismatch function) - either
            'fit' (to use results from the original Q-fitting routine) or 'injection' (to use Q derived
            from source injection simulations).
        overrideNoise (:obj:`float`, optional): If given, overwrite the contents of all noise maps with
            the given value (used for testing).
        maxFlags (:obj:`int`, optional): If given, apply the flags mask, cutting areas where the flags
            value is greater than this value.

    Attributes:
        SNRCut (:obj:`float`): Completeness will be computed relative to this signal-to-noise selection cut
            (labelled as `fixed_SNR` in the catalogs produced by :ref:`nemoCommand`).
        footprint (:obj:`str`): Use this to specify a footprint, if any are defined in the Nemo config
            used to produce the ``selFn`` dir (e.g., 'DES', 'HSC', 'KiDS' etc.). The default of ``None``
            uses the whole survey footprint.
        applyMFDebiasCorrection (:obj:`bool`): Set to `False` to disable the Eddington bias correction of
            mass estimates. Probably only useful for debugging.
        zStep (:obj:`float`): Use this to set the binning in redshift for completeness calculations.
        zMax (:obj:`float`): Sets the maximum z used for completeness calculations.
        tileNames (:obj:`list`): The list of tiles used by the SelFn object (default of None uses all tiles).
        WCSDict (:obj:`dict`): A dictionary indexed by `tileName`, containing :obj:`astWCS.WCS` objects
            that describe the mapping between pixel coords and (RA, dec) coords in each tile.
        areaMaskDict (:obj:`dict`): A dictionary containing the survey area masks, indexed by tileName.
            Values > 0 in these masks define the cluster or source search area.
        scalingRelationDict (:obj:`dict`): A dictionary of scaling relation parameters (see example Nemo
            config files for the format).
        Q (:class:`nemo.signals.QFit`): An object for calculating the filter mismatch function, referred
            to as `Q` in the ACT papers from `Hasselfield et al. (2013) <http://adsabs.harvard.edu/abs/2013JCAP...07..008H>`_
            onwards.
        RMSDict (:obj:`dict`): A dictionary of RMS tables, indexed by tileName. Each RMSTable contains 
            the noise level by area, as returned by :meth:`getRMSTab`.
        totalAreaDeg2 (:obj:`float`): The total area in square degrees, as measured from the survey mask,
            for the given set of tiles and footprint.
        fRelDict (:obj:`dict`): A dictionary of weights used for relativistic corrections, indexed by
            `tileName`.
        mockSurvey (:class:`nemo.MockSurvey.MockSurvey`): A :class:`MockSurvey` object, used for halo mass function 
            calculations and generating mock catalogs.

    Note: 
        Some of the methods of this class are experimental and not necessarily well tested.

    Note:
        Once SNRCut is set for this object, it cannot be changed later (well, it can... but changing
        self.SNRCut will not update anything, for the moment).
            
    """
        
    def __init__(self, selFnDir, SNRCut, configFileName = None, footprint = None, zStep = 0.01,
                 zMin = 0.0, zMax = 3.0, minMass = 5e13, numMassBins = 200, tileNames = None, mockOversampleFactor = 1.0,
                 downsampleRMS = 2, applyMFDebiasCorrection = True, applyRelativisticCorrection = True,
                 setUpAreaMask = False, enableCompletenessCalc = True, delta = 500, rhoType = 'critical',
                 massFunction = 'Tinker08', maxTheta500Arcmin = None, method = 'fast',
                 QSource = 'fit', useAverageQ = False, theoryCode = 'CCL', noiseCut = None, biasModel = None,
                 truncateDeltaSNR = 3.0, overrideNoise = None, massBinsTheory = 2000, zStepTheory = 0.001,
                 maxFlags = None):
        
        self.SNRCut=SNRCut
        self.biasModel=biasModel
        if footprint == 'full':
            footprint=None
        self.footprint=footprint
        self.downsampleRMS=downsampleRMS
        self.applyMFDebiasCorrection=applyMFDebiasCorrection
        self.applyRelativisticCorrection=applyRelativisticCorrection
        self.selFnDir=selFnDir
        self.zStep=zStep
        self.numMassBins=numMassBins
        self.maxTheta500Arcmin=maxTheta500Arcmin
        self.useAverageQ=useAverageQ
        self.massBinsTheory=massBinsTheory
        self.zStepTheory=zStepTheory
        self.truncateDeltaSNR=truncateDeltaSNR
        self.maxFlags=maxFlags

        if configFileName is None:
            configFileName=self.selFnDir+os.path.sep+"config.yml"
            if os.path.exists(configFileName) == False:
                raise Exception("No config .yml file found in selFnDir and no other location given.")
        self._config=startUp.NemoConfig(configFileName, makeOutputDirs = False, setUpMaps = False, verbose = False,
                                        selFnDir = self.selFnDir)
        parDict=self._config.parDict

        if tileNames is not None:
            self.tileNames=tileNames
        else:
            self.tileNames=self._config.tileNames
        
        deprecatedMethods=['montecarlo']
        if method in deprecatedMethods:
            raise Exception("Selection function completeness calculation method '%s' is deprecated - use a valid method (e.g., 'fast')" % (method))
        self.method=method

        # Needed for generating mock samples directly
        self.photFilterLabel=self._config.parDict['photFilter']
            
        # Check that any given footprint is defined - if not, give a useful error message
        if footprint is not None:
            if 'selFnFootprints' not in parDict.keys():
                raise Exception("No footprints defined in .yml config file")
            else:
                labelsList=[]
                for footprintDict in parDict['selFnFootprints']:
                    labelsList.append(footprintDict['label'])
                if footprint not in labelsList:
                    raise Exception("Footprint '%s' not found in selFnFootprints - check .yml config file" % (footprintLabel))
        
        # Load area masks
        if setUpAreaMask == True:
            self._setUpAreaMask()
        else:
            self.tileTab=None
            self.WCSDict=None
            self.areaMaskDict=None
                
        if enableCompletenessCalc == True:
            
            # We can now specify multiple scaling relations in the config - but we only use the first one here
            self.scalingRelationDict=parDict['massOptions']['scalingRelations'][0]

            # If the config didn't give cosmological parameters, put in defaults
            defaults={'H0': 70.0, 'Om0': 0.30, 'Ob0': 0.05, 'sigma8': 0.8, 'ns': 0.95}
            for key in defaults:
                if key not in parDict['massOptions'].keys():
                    parDict['massOptions'][key]=defaults[key]

            # For caching / faster updating
            self._last_H0=None
            self._last_Om0=None
            self._last_Ob0=None
            self._last_sigma8=None
            self._last_ns=None
            self._last_scalingRelationDict=None
            self._y0GridCache={}
            self._theta500GridCache={}

            # We should be able to do everything (except clustering) with this
            # NOTE: Some tiles may be empty, so we'll exclude them from tileNames list here
            RMSTabFileName=self.selFnDir+os.path.sep+"RMSTab.fits"
            if footprint is not None:
                RMSTabFileName=RMSTabFileName.replace(".fits", "_%s.fits" % (footprint))
            if maxFlags is not None:
                RMSTabFileName=RMSTabFileName.replace(".fits", "_maxFlags%d.fits" % (maxFlags))
            if os.path.exists(RMSTabFileName) == False:
                raise Exception("File %s not found - perhaps check maxFlags value given?" % (RMSTabFileName))
            self.RMSTab=atpy.Table().read(RMSTabFileName)
            if overrideNoise is not None:
                self.RMSTab['y0RMS']=overrideNoise
            # Sanitise just in case [this is an edge case that has happened on one sim]
            self.RMSTab=self.RMSTab[self.RMSTab['areaDeg2'] > 0]
            if noiseCut is not None:
                self.RMSTab=self.RMSTab[self.RMSTab['y0RMS'] < noiseCut]
            self.RMSDict={}
            tileNames=[]
            totalAreaDeg2=0.0 # Doing it this way so that tileNames can be chosen and fed into selFn
            for tileName in self.tileNames:
                tileTab=self.RMSTab[self.RMSTab['tileName'] == tileName]
                if downsampleRMS > 1 and len(tileTab) > 0:
                    tileTab=downsampleRMSTab(tileTab, downsampleRMS)
                if len(tileTab) > 0:    # We may have some blank tiles...
                    self.RMSDict[tileName]=tileTab
                    tileNames.append(tileName)
                    totalAreaDeg2=totalAreaDeg2+tileTab['areaDeg2'].sum()
            self.tileNames=tileNames
            self.totalAreaDeg2=totalAreaDeg2
            if totalAreaDeg2 == 0:
                raise FootprintError
            # If want a plot of noise distribution
            # plt.hist(self.RMSTab['y0RMS'], weights = self.RMSTab['areaDeg2'], bins  = 100, density=True)

            # For weighting - arrays where entries correspond with tileNames list
            tileAreas=[]    
            for tileName in self.tileNames:
                areaDeg2=self.RMSTab[self.RMSTab['tileName'] == tileName]['areaDeg2'].sum()
                tileAreas.append(areaDeg2)
            self.tileAreas=np.array(tileAreas)
            self.fracArea=self.tileAreas/self.totalAreaDeg2

            # Check of area consistency
            #checkAreaDeg2=0
            #for tileName in self.tileNames:
                #checkAreaDeg2=checkAreaDeg2+(self.areaMaskDict[tileName]*maps.getPixelAreaArcmin2Map(self.areaMaskDict[tileName].shape, self.WCSDict[tileName])/(60**2)).sum()

            # For quick sample generation
            self.mockOversampleFactor=mockOversampleFactor
            self.y0NoiseAverageDict={}
            for tileName in self.tileNames:
                RMSTab=self.RMSDict[tileName]
                areaWeights=RMSTab['areaDeg2'].data/RMSTab['areaDeg2'].data.sum()
                if areaWeights.sum() > 0:
                    self.y0NoiseAverageDict[tileName]=np.average(RMSTab['y0RMS'].data, 
                                                                 weights = areaWeights)            
            assert(len(self.tileNames) == len(self.tileAreas))
            assert(len(self.tileNames) == len(self.y0NoiseAverageDict.keys()))

            # For relativistic corrections
            self.fRelDict=signals.loadFRelWeights(self.selFnDir+os.path.sep+"fRelWeights.fits")

            # Stuff from the source injection sims (now required for completeness calculation)
            if self.method == 'injection':
                injDataPath=self.selFnDir+os.path.sep+"sourceInjectionData.fits"
                inputDataPath=self.selFnDir+os.path.sep+"sourceInjectionInputCatalog.fits"
                if os.path.exists(injDataPath) == False or os.path.exists(inputDataPath) == False:
                    raise Exception("%s not found - run a source injection test to generate (now required for completeness calculations)." % (injDataPath))
                injTab=atpy.Table().read(injDataPath)
                inputTab=atpy.Table().read(inputDataPath)
                if footprint is not None:
                    injTab=self.cutCatalogToSurveyArea(injTab)
                    inputTab=self.cutCatalogToSurveyArea(inputTab)
                theta500s, binCentres, compThetaGrid, thetaQ=_parseSourceInjectionData(injTab, inputTab, self.SNRCut)
                self.compThetaInterpolator=interpolate.RectBivariateSpline(theta500s, binCentres,
                                                                        compThetaGrid, kx = 3, ky = 3)

            # Q can now either be from the 'classic' fit or the source injection sims
            self.Q=signals.QFit(QSource = QSource, selFnDir = self.selFnDir, tileNames = tileNames)

            # Initial cosmology set-up
            H0=parDict['massOptions']['H0']
            Om0=parDict['massOptions']['Om0']
            Ob0=parDict['massOptions']['Ob0']
            sigma8=parDict['massOptions']['sigma8']
            ns=parDict['massOptions']['ns']
            self.mockSurvey=MockSurvey.MockSurvey(minMass, self.totalAreaDeg2, zMin, zMax, H0, Om0, Ob0, sigma8, ns,
                                                  zStep = self.zStepTheory, numMassBins = self.massBinsTheory,
                                                  delta = delta, rhoType = rhoType, massFunction = massFunction,
                                                  theoryCode = theoryCode)
            # This will usually be a coarser gridding than used in the MockSurvey object
            # NOTE: binning mass in log(M) rather than log10(M) here
            numRedshiftBins=int(round((zMax-zMin)/zStep))
            self.zBinEdges=np.linspace(zMin, zMax, numRedshiftBins+1)
            # self.zBinEdges=np.arange(zMin, zMax+zStep, zStep)
            self.logMBinEdges=np.linspace(np.log(self.mockSurvey.M).min(), np.log(self.mockSurvey.M).max(), self.numMassBins)
            self.clusterCount=np.zeros([self.zBinEdges.shape[0]-1, self.logMBinEdges.shape[0]-1])
            self.z=(self.zBinEdges[1:]+self.zBinEdges[:-1])/2
            self.a=1./(1+self.z)
            log10MBinEdges=np.log10(np.exp(self.logMBinEdges))
            self.log10M=(log10MBinEdges[1:]+log10MBinEdges[:-1])/2
            self.M=np.power(10, self.log10M)
            # Below used for Q calc only, changes only if cosmological parameters do
            self._theta500Grid=np.zeros([self.z.shape[0], self.log10M.shape[0]])
            self.update(H0, Om0, Ob0, sigma8, ns)


    def _setUpAreaMask(self):
        """Sets-up WCS info and loads area masks (taking into account the footprint, if specified)
        - needed for quick position checks etc.

        Note:
            This takes into account the maxFlags property.
        
        """
        
        # This takes ~20 sec to set-up - we could cache, or it's overhead when initialising SelFn
        # We could do a lot of this by just making the area mask, mass limit maps etc. MEFs
        # But then we wouldn't want to lazily load them (probably)
        self.tileTab=atpy.Table()
        self.tileTab.add_column(atpy.Column(list(self.tileNames), 'tileName'))
        self.tileTab.add_column(atpy.Column(np.zeros(len(self.tileNames)), 'RAMin'))
        self.tileTab.add_column(atpy.Column(np.zeros(len(self.tileNames)), 'RAMax'))
        self.tileTab.add_column(atpy.Column(np.zeros(len(self.tileNames)), 'decMin'))
        self.tileTab.add_column(atpy.Column(np.zeros(len(self.tileNames)), 'decMax'))
        self.WCSDict={}
        self.areaMaskDict={}
        for row in self.tileTab:
            if self.footprint is None:
                areaMap, wcs=loadAreaMask(row['tileName'], self.selFnDir)
            else:
                areaMap, wcs=loadIntersectionMask(row['tileName'], self.selFnDir, self.footprint)
            if self.maxFlags is not None:
                flagMap, wcs=loadFlagMask(row['tileName'], self.selFnDir)
                areaMap[flagMap > self.maxFlags]=0
            self.WCSDict[row['tileName']]=wcs.copy()
            self.areaMaskDict[row['tileName']]=areaMap
            ra0, dec0=self.WCSDict[row['tileName']].pix2wcs(0, 0)
            ra1, dec1=self.WCSDict[row['tileName']].pix2wcs(wcs.header['NAXIS1'], wcs.header['NAXIS2'])
            if ra1 > ra0:
                ra1=-(360-ra1)
            row['RAMin']=min([ra0, ra1])
            row['RAMax']=max([ra0, ra1])
            row['decMin']=min([dec0, dec1])
            row['decMax']=max([dec0, dec1])
            

    def cutCatalogToSurveyArea(self, catalog):
        """Checks that the coordinates of objects in the given catalog are within the survey area
        (taking into account a footprint, if specified) and cuts the catalog accordingly.

        Args:
            catalog (:obj:`astropy.table.Table`): The input catalog, as an astropy Table.

        Returns:
            A catalog (:obj:`astropy.table.Table`), cut to the survey area.

        """

        RAKey, decKey=catalogs.getTableRADecKeys(catalog)
        mask=self.checkCoordsInAreaMask(catalog[RAKey], catalog[decKey])
        return catalog[mask]


    def checkCoordsInAreaMask(self, RADeg, decDeg):
        """Checks if the given RA, dec coords are in valid regions of the map.
        
        Args:
            RADeg (:obj:`float` or :obj:`np.ndarray`): RA in decimal degrees.
            decDeg (:obj:`float` or :obj:`np.ndarray`): Dec in decimal degrees.
        
        Returns:
            `True` if the coordinates are in the area mask mask, `False` if not.
            
        """

        if self.tileTab is None:
            self._setUpAreaMask()

        RADeg=np.array(RADeg)
        decDeg=np.array(decDeg)
        if RADeg.shape == ():
            RADeg=np.array([RADeg])
        if decDeg.shape == ():
            decDeg=np.array([decDeg])
        inMask=np.zeros(len(RADeg), dtype = bool)
        tabList=[]
        for tileName in self.tileNames:
            wcs=self.WCSDict[tileName]
            areaMask=self.areaMaskDict[tileName]
            if areaMask.sum() > 0:
                coords=wcs.wcs2pix(RADeg, decDeg)
                coords=np.array(np.round(coords), dtype = int)
                mask1=np.logical_and(coords[:, 0] >= 0, coords[:, 1] >= 0)
                mask2=np.logical_and(coords[:, 0] < areaMask.shape[1], coords[:, 1] < areaMask.shape[0])
                mask=np.logical_and(mask1, mask2)
                inMask[mask]=inMask[mask]+areaMask[coords[:, 1][mask], coords[:, 0][mask]]

        return inMask
        

    def _doClusterCount(self, numRedshiftPoints = 25, numMassPoints = 25):
        """Update cluster counts (prior to applying selection function) on the gridding used by this object.
        This routine used to be in the MockSurvey class.

        Args:
            numRedshiftPoints (int): Number of points to subdivide each SelFn redshift bin into.
            numMassPoints (int): Number of points to subdivide each SelFn mass bin into.

        """

        # t0=time.time()
        summedOverBins=0
        for zi in range(self.zBinEdges.shape[0]-1):
            zMin=self.zBinEdges[zi]
            zMax=self.zBinEdges[zi+1]
            zPoints=np.linspace(zMin, zMax, numRedshiftPoints)
            for mi in range(self.logMBinEdges.shape[0]-1):
                mMin=self.logMBinEdges[mi]
                mMax=self.logMBinEdges[mi+1]
                mPoints=np.linspace(mMin, mMax, numMassPoints)
                dndz=integrate.simpson(self.mockSurvey.dndmdzInterpolator(mPoints, zPoints), x = mPoints, axis = 0)
                # dndm=integrate.simpson(dndmdz_interpolator(mPoints, zPoints), x = zPoints, axis = 1)
                numClusters_from_dndz=integrate.simpson(dndz, x = zPoints)
                # numClusters_from_dndm=integrate.simpson(dndm, x = mPoints)
                summedOverBins=summedOverBins+numClusters_from_dndz
                self.clusterCount[zi, mi]=numClusters_from_dndz
        res=abs(1-(summedOverBins/self.mockSurvey.numClusters))
        if res > 1e-2:
            raise Exception("Coarse (mass, z) binning used in SelFn does not accurately reproduce total theory cluster count (res = %.3e)." % (res))
        # t1=time.time()


    def update(self, H0, Om0, Ob0, sigma8, ns, scalingRelationDict = None):
        """Re-calculates the survey-average selection function for a given set of cosmological and scaling
        relation parameters.
        
        Returns:
            None - attributes such as :attr:`compMz`, the (log\ :sub:`10` mass, z) completeness grid, are
            updated in-place.

        """

        if scalingRelationDict is not None:
            self.scalingRelationDict=scalingRelationDict
        self.mockSurvey.update(H0, Om0, Ob0, sigma8, ns)

        self._doClusterCount()

        # Interpolator for erf, can be used by 'fast' completeness method
        # Only faster than direct calls of self._get_erf_diff for < 10000 points
        # Doesn't get any faster for < 1000 points
        self._compSNRLimInterp=None
        if self.method == 'faster':
            raise Exception("This has been disabled as it's not helpful unless you use very fine (M,z) gridding")
            qmin=0 # Should always be zero
            qmax=100
            qin=np.linspace(qmin, qmax, 2000)
            erf_compl=self._get_erf_diff(qin, self.SNRCut, qmax, self.SNRCut)
            try:
                truncIndex=np.where(erf_compl == 1)[0][0]
                self._compSNRLimInterp=interpolate.interp1d(qin[:truncIndex], erf_compl[:truncIndex], fill_value = 1, bounds_error = False)
            except:
                pass

        # theta500s for Q calc
        if self.mockSurvey._cosmoUpdated is True or self._theta500Grid.sum() == 0:
            zRange=self.z
            Ez=ccl.h_over_h0(self.mockSurvey.cosmoModel, self.a)
            Ez2=np.power(Ez, 2)
            DAz=ccl.angular_diameter_distance(self.mockSurvey.cosmoModel, self.a)
            criticalDensity=ccl.physical_constants.RHO_CRITICAL*(Ez*self.mockSurvey.cosmoModel['h'])**2
            for i in range(len(zRange)):
                if self.mockSurvey.delta != 500 or self.mockSurvey.rhoType != "critical":
                    M500s=self.mockSurvey._transToM500c(self.mockSurvey.cosmoModel, self.M, self.a[i])
                    self._log10M500s=np.log10(M500s)
                else:
                    M500s=self.M
                    self._log10M500s=self.log10M
                R500Mpc=np.power((3*M500s)/(4*np.pi*500*criticalDensity[i]), 1.0/3.0)
                self._theta500Grid[i]=np.degrees(np.arctan(R500Mpc/DAz[i]))*60.0

        if self.method == 'injection':
            # WARNING: Here there is no Q, and y0 is true y0, NOT y0~
            # This grid is purely mapping 'input' signals to M, z via the scaling relation
            y0Grid=self._makeSignalGrid(applyQ = False)
            compMz=np.zeros(y0Grid.shape)
            for i in range(len(self.z)):
                try:
                    compMz[i]=np.diag(self.compThetaInterpolator(self._theta500Grid[i], y0Grid[i]/1e-04))
                except:
                    # If above fails, it's because relativistic correction is on, and y0 doesn't always increase
                    for j in range(y0Grid[i].shape[0]):
                        compMz[i][j]=self.compThetaInterpolator(self._theta500Grid[i][j], y0Grid[i][j]/1e-04)
            self.compMz=compMz
            self.y0TildeGrid=self.Q.getQ(self._theta500Grid)*y0Grid
            if self.scalingRelationDict['sigma_int'] > 0:
                # Fiddling with Gaussian filter params has no effect
                mode='nearest'
                truncate=4.0
                logy0Grid=np.log(y0Grid)
                for i in range(logy0Grid.shape[0]):
                    dy=np.mean(np.gradient(logy0Grid[i]))
                    if dy > 0:
                        npix=self.scalingRelationDict['sigma_int']/dy
                        npix=npix*0.8   # Making this correction seems to work but not sure why yet
                        self.clusterCount[i]=ndimage.gaussian_filter1d(self.clusterCount[i], npix,
                                                                                  mode = mode, truncate = truncate)

        elif self.method == 'fast' or self.method == 'faster':
            compMzCube=np.zeros([len(self.tileNames), self.clusterCount.shape[0], self.clusterCount.shape[1]])
            t0=time.time()
            for tileIndex in range(len(self.tileNames)):
                tileName=self.tileNames[tileIndex]
                compMzCube[tileIndex]=self.calcFastCompletenessInTile(tileName, return_y0Grid = False)
            self.compMz=np.average(compMzCube, axis = 0, weights = self.fracArea)

        # Deals with corner at high S/N, high-z sometimes weirdly having lower than 1 completeness
        for i in range(self.compMz.shape[0]):
            minIndex=np.where(self.compMz[i] >= 1)[0]
            if len(minIndex) > 0:
                minIndex=minIndex[0]
                self.compMz[i][minIndex:]=1
        self.compMz[self.compMz < 0]=0

        self.predObsCount=self.compMz*self.clusterCount

        # For caching/update checks
        self._last_scalingRelationDict=self.scalingRelationDict
        self._last_H0=self.mockSurvey.H0
        self._last_Om0=self.mockSurvey.Om0
        self._last_Ob0=self.mockSurvey.Ob0
        self._last_sigma8=self.mockSurvey.sigma8
        self._last_ns=self.mockSurvey.ns


    def calcFastCompletenessInTile(self, tileName, return_y0Grid = False, RMSTab = None):
        """Calculate completeness on the (M, z) grid for the given tile using the fast method.

        Args:
            tileName (str): Name of the tile
            return_y0Grid (bool, optional): If given, return the ``fixed_y_c`` signal calculated
                on the (M, z) grid
            RMSTab (:obj:`astropy.table.Table`, optional): If given, this is used instead of the
                RMS noise information held by the ``SelFn`` object.

        Returns:
            2d array (completeness on the M,z grid for the tile with same dimensions as
            ``self.clusterCount``)

        """

        t00=time.time()
        # NOTE: y0Grid, theta500 grid are cached and only re-calculated if parameters change
        if self.useAverageQ == False:
            y0Grid=self._makeSignalGrid(tileName = tileName)
        else:
            y0Grid=self._makeSignalGrid(tileName = None)
        compMzTile=np.zeros(y0Grid.shape)
        # Calculate completeness using area-weighted average
        # NOTE: RMSTab that is fed in here can be downsampled in noise resolution for speed
        if RMSTab is None:
            RMSTab=self.RMSDict[tileName]
        areaWeights=RMSTab['areaDeg2']/RMSTab['areaDeg2'].sum()
        t22=time.time()
        # SOLikeT style - tiny loops are faster than doing this with array functions
        for i in range(len(RMSTab)):
            if self.biasModel is not None:
                trueSNR=y0Grid/RMSTab['y0RMS'][i]
                # Old model
                # corrFactors=self.biasModel['func'](trueSNR, self.biasModel['params'][0], self.biasModel['params'][1], self.biasModel['params'][2])
                # New model [if we do want the old model, should call it like this anyway]
                corrFactors=self.biasModel['func'](trueSNR, self.biasModel['params'])
                # Some models may give unphysically large correction factors when extrapolated S/N -> 0
                # So, we truncate this at some level (e.g. 3-sigma) below the S/N cut
                if self.truncateDeltaSNR is not None:
                    corrFactors[trueSNR < self.SNRCut-self.truncateDeltaSNR]=1.0
            else:
                corrFactors=np.ones(y0Grid.shape)
            # If we decide to bring back the 'faster' method, see below
            # if self._compSNRLimInterp is not None:
            #     compMzCube[tileIndex]=compMzCube[tileIndex]+self._compSNRLimInterp((y0Grid*corrFactors)/RMSTab['y0RMS'][i])*areaWeights[i]
            # else:
            #     compMzCube[tileIndex]=compMzCube[tileIndex]+self._get_erf_diff((y0Grid*corrFactors)/RMSTab['y0RMS'][i], self.SNRCut, 1e5, self.SNRCut)*areaWeights[i]
            if self.scalingRelationDict['sigma_int'] == 0:
                compMzTile=compMzTile+self._get_erf_diff((y0Grid*corrFactors)/RMSTab['y0RMS'][i], self.SNRCut, 1e5, self.SNRCut)*areaWeights[i]
            else:
                # SOLikeT style
                scatter=self.scalingRelationDict['sigma_int']
                lnyy=np.linspace(np.min(np.log(y0Grid)), np.max(np.log(y0Grid)), 44)
                yy0=np.exp(lnyy)
                mu=np.float32(np.log(y0Grid*corrFactors))
                fac=np.float32(1./np.sqrt(2.*np.pi*scatter**2))
                arg=self._get_erf_diff(yy0/RMSTab['y0RMS'][i], self.SNRCut, 1e5, self.SNRCut)
                cc=np.float32(arg*areaWeights[i])
                arg0=np.float32((lnyy[:, None,None]-mu)/(np.sqrt(2.)*scatter))
                args=fac*np.exp(np.float32(-arg0**2.)) * cc[:, None,None]
                compMzTile+=np.trapz(np.float32(args), x=lnyy, axis=0)

        if self.maxTheta500Arcmin is not None:
            compMzTile=compMzTile*np.array(self._theta500Grid < self.maxTheta500Arcmin, dtype = float)

        if return_y0Grid is True:
            return compMzTile, y0Grid
        else:
            return compMzTile


    def _get_erf_diff(self, qin, qmin, qmax, qcut):
        """Based on SOLikeT - qin == input SNR (y0/RMS on the grid)"""
        arg1 = (qin - qmax)/np.sqrt(2.)
        if qmin > qcut:
            qlim = qmin
        else:
            qlim = qcut
        arg2 = (qin - qlim)/np.sqrt(2.)
        erf_compl = (special.erf(arg2) - special.erf(arg1)) / 2.

        return erf_compl


    def _checkIfParametersUpdated(self):
        paramsUpdated=False
        if self._last_H0 != self.mockSurvey.H0:
            paramsUpdated=True
        if self._last_Om0 != self.mockSurvey.Om0:
            paramsUpdated=True
        if self._last_Ob0 != self.mockSurvey.Ob0:
            paramsUpdated=True
        if self._last_sigma8 != self.mockSurvey.sigma8:
            paramsUpdated=True
        if self._last_ns != self.mockSurvey.ns:
            paramsUpdated=True
        if self._last_scalingRelationDict != self.scalingRelationDict:
            paramsUpdated=True

        return paramsUpdated


    def _makeSignalGrid(self, applyQ = True, tileName = None):
        """Returns y0~ grid. tileName here is optional and only used for Q.

        Previous results may be cached and re-used

        """

        y0Grid=None

        if self._checkIfParametersUpdated() == False:
            if tileName in self._y0GridCache.keys():
                y0Grid=self._y0GridCache[tileName]

        if y0Grid is None:
            zRange=self.z
            tenToA0, B0, Mpivot, sigma_int=[self.scalingRelationDict['tenToA0'], self.scalingRelationDict['B0'],
                                            self.scalingRelationDict['Mpivot'], self.scalingRelationDict['sigma_int']]
            # Optional extras for z evolution
            if 'onePlusRedshift_power' not in self.scalingRelationDict.keys():
                onePlusRedshift_power=0.0
            else:
                onePlusRedshift_power=self.scalingRelationDict['onePlusRedshift_power']
            if 'Ez_gamma' not in self.scalingRelationDict.keys():
                Ez_gamma=2.0 # Default self-similar
            else:
                Ez_gamma=self.scalingRelationDict['Ez_gamma']

            y0Grid=np.zeros([zRange.shape[0], self.clusterCount.shape[1]])
            # NOTE: Still called Ez2, but now has gamma option enabled
            Ez2=np.power(ccl.h_over_h0(self.mockSurvey.cosmoModel, 1/(1+zRange)), Ez_gamma)
            for i in range(len(zRange)):
                zk=zRange[i]
                # NOTE: Now we have two z bin schemes (one in MockSurvey, one in SelFn) need to take care here with indices
                k=np.argmin(abs(self.mockSurvey.z-zk))
                Qs_zk=self.Q.getQ(self._theta500Grid[i], zk, tileName = tileName)
                #Qs_zk=self.compQInterpolator(theta500s_zk) # Survey-averaged Q from injection sims
                true_y0s_zk=tenToA0*Ez2[i]*np.power(np.power(10, self.log10M)/Mpivot, 1+B0)*np.power(1+zk, onePlusRedshift_power)
                if applyQ == True:
                    true_y0s_zk=true_y0s_zk*Qs_zk
                if self.applyRelativisticCorrection == True:
                    fRels_zk=interpolate.splev(self._log10M500s, self.mockSurvey.fRelSplines[k])
                    true_y0s_zk=true_y0s_zk*fRels_zk
                y0Grid[i]=true_y0s_zk
            # For some cosmological parameters, we can still get the odd -ve y0
            y0Grid[y0Grid <= 0] = 1e-9
            self._y0GridCache[tileName]=y0Grid

        return y0Grid


    def projectCatalogToMz(self, tab):
        """Project a Nemo cluster catalog (an astropy Table) into the (log\ :sub:`10` mass, z) grid, taking
        into account the uncertainties on y0, and redshift. Note that if the redshift error is non-zero, this
        is a lot slower.
        
        Args:
            tab (:obj:`astropy.table.Table`): A Nemo cluster catalog, containing `redshift` and `redshiftErr`
                columns.
                
        Returns:
            A 2d array containing the projection of the catalog on the (log\ :sub:`10` mass, z) grid.
        
        """
        
        catProjectedMz=np.zeros(self.clusterCount.shape)
        tenToA0, B0, Mpivot, sigma_int=self.scalingRelationDict['tenToA0'], self.scalingRelationDict['B0'], \
                                       self.scalingRelationDict['Mpivot'], self.scalingRelationDict['sigma_int']

        for row in tab:
            tileName=row['tileName']
            z=row['redshift']
            zErr=row['redshiftErr']
            y0=row['fixed_y_c']*1e-4
            y0Err=row['fixed_err_y_c']*1e-4
            P=signals.calcPMass(y0, y0Err, z, zErr, self.Q, self.mockSurvey, 
                                  tenToA0 = tenToA0, B0 = B0, Mpivot = Mpivot, sigma_int = sigma_int, 
                                  applyMFDebiasCorrection = self.applyMFDebiasCorrection, 
                                  fRelWeightsDict = self.fRelDict[tileName],
                                  return2D = True, tileName = tileName)
            # Paste into (M, z) grid
            catProjectedMz=catProjectedMz+P # For return2D = True, P is normalised such that 2D array sum is 1
        
        return catProjectedMz


    def projectCatalogToMz_simple(self, tab):
        """Project a Nemo cluster catalog (an astropy Table) into the (log\ :sub:`10` mass, z) grid. This version
        doesn't take into account any uncertainties (which may be okay if your binning is coarse enough).
                
        Args:
            tab (:obj:`astropy.table.Table`): A Nemo cluster catalog, containing `redshift` and `redshiftErr`
                columns.
                
        Returns:
            A 2d array containing the projection of the catalog on the (log\ :sub:`10` mass, z) grid.
        
        """
        
        tenToA0, B0, Mpivot, sigma_int=self.scalingRelationDict['tenToA0'], self.scalingRelationDict['B0'], \
                                       self.scalingRelationDict['Mpivot'], self.scalingRelationDict['sigma_int']
                                           
        obs_log10Ms=[]
        for row in tab:
            tileName=row['tileName']
            z=row['redshift']
            zErr=row['redshiftErr']
            y0=row['fixed_y_c']*1e-4
            y0Err=row['fixed_err_y_c']*1e-4
            massDict=signals.calcMass(y0, y0Err, z, zErr, self.Q, tenToA0 = tenToA0, B0 = B0, Mpivot = Mpivot, 
                                            sigma_int = sigma_int, mockSurvey = self.mockSurvey, 
                                            applyMFDebiasCorrection = self.applyMFDebiasCorrection,
                                            fRelWeightsDict = self.fRelDict[tileName],
                                            calcErrors = False, tileName = tileName)
            obs_log10Ms.append(14+np.log10(massDict['M500']))
        obsGrid, obs_log10MBinEdges, obs_zBinEdges=np.histogram2d(obs_log10Ms, tab['redshift'], 
                                                                  bins = [self.mockSurvey.log10MBinEdges, 
                                                                          self.mockSurvey.zBinEdges])
        obsGrid=obsGrid.transpose()
            
        return obsGrid


    def addPDetToCatalog(self, tab):
        """Given a catalog, add a column named `Pdet`, containing the detection probability.
        
        Args:
            tab (:obj:`astropy.table.Table`): A Nemo cluster catalog.
                
        Returns:
            Catalog with `Pdet` column added (:obj:`astropy.table.Table`)
        
        """

        log_y0Lim=np.log(self.SNRCut*tab['fixed_err_y_c']*1e-4)
        log_y0=np.log(tab['fixed_y_c']*1e-4)
        log_y0Err=1/tab['fixed_SNR']
        sigma_int=self.scalingRelationDict['sigma_int']
        log_totalErr=np.sqrt(log_y0Err**2 + sigma_int**2)
        Pdet=np.zeros(len(tab))
        for i in range(len(Pdet)):
            Pdet[i]=stats.norm.sf(log_y0Lim[i], loc = log_y0[i], scale = log_totalErr[i])
        tab['Pdet']=Pdet
        
        return tab

    
    def generateMockSample(self, mockOversampleFactor = None, applyPoissonScatter = True):
        """Returns a mock catalog (but with no object coordinate information).
        
        Args:
            mockOversampleFactor (:obj:`float`, optional): Scale the number of objects in the mock by this
                factor (e.g., this can be used to apply an overall area re-scaling). If given, this
                overrides self.mockOversampleFactor.
            applyPoissonScatter (:obj:`bool`, optional): If True, apply Poisson scatter to the cluster
                number counts when generating the mock catalog.
        
        """
        
        if mockOversampleFactor is None:
            mockOversampleFactor=self.mockOversampleFactor

        mockTabsList=[]
        for tileName, areaDeg2 in zip(self.tileNames, self.tileAreas):
            mockTab=self.mockSurvey.drawSample(self.RMSDict[tileName], self.scalingRelationDict,
                                               self.Q, wcs = None, 
                                               photFilterLabel = self.photFilterLabel, tileName = tileName, 
                                               makeNames = False,
                                               SNRLimit = self.SNRCut,
                                               areaDeg2 = areaDeg2*mockOversampleFactor,
                                               applyPoissonScatter = applyPoissonScatter,
                                               applyIntrinsicScatter = True,
                                               applyNoiseScatter = True,
                                               applyRelativisticCorrection = self.applyRelativisticCorrection,
                                               biasModel = self.biasModel)
            if mockTab is not None and len(mockTab) > 0:
                mockTabsList.append(mockTab)
        tab=atpy.vstack(mockTabsList)

        return tab
    
    
    def getMassLimit(self, completenessFraction, zBinEdges = None):
        """Return the mass limit (units of 10\ :sup:`14` MSun) as a function of redshift, for the given
        completeness level. 
        
        Args:
            completenessFraction (:obj:`float`): The completeness fraction (a number between 0 and 1) at
                which to return the mass limit.
            zBinEdges (:obj:`np.ndarray`, optional): The redshifts at which the completeness will be
                evaluated. If not given, :attr:`self.mockSurvey.z` will be used.
                
        Returns:
            Mass completeness by redshift (:obj:`np.ndarray`).
        
        """
        
        if zBinEdges is None:
            zBinEdges=self.zBinEdges

        massLimit=np.power(10, self.log10M[np.argmin(abs(self.compMz-completenessFraction), axis = 1)])/1e14
        if len(zBinEdges) > 0:
            binnedMassLimit=np.zeros(len(zBinEdges)-1)
            for i in range(len(zBinEdges)-1):
                binnedMassLimit[i]=np.average(massLimit[np.logical_and(self.z > zBinEdges[i], self.z <= zBinEdges[i+1])])
            massLimit=binnedMassLimit

        return massLimit

#------------------------------------------------------------------------------------------------------------
def optBiasModelFunc(snr, params):
    """Optimization bias model function, of the form ``corrFactor = 1 + p1/x + p2/x**2 + ... + pn/x**n``
    where p1...pn are fit coefficents given as the params array. This is for use with the `fast` completeness
    method of the ``SelFn`` class.

    Args:
        snr (:obj:`np.ndarray`): Array of true signal-to-noise ratio values
        params (:obj:`np.ndarray`): Fit coefficients.

    Returns:
        Array of correction factors.

    """

    model=np.ones(snr.shape)
    index=1
    for p in params:
        model=model+p/(snr**index)
        index=index+1

    return model

#------------------------------------------------------------------------------------------------------------
def optBiasPowerModelFunc(snr, param):
    """Optimization bias model function, of the form ``corrFactor = 1 + 1/x**p``
    where p is a fit paramer. This is for use with the `fast` completeness
    method of the ``SelFn`` class.

    Args:
        snr (:obj:`np.ndarray`): Array of true signal-to-noise ratio values
        params (:obj:`float`): Fit parameter.

    Returns:
        Array of correction factors.

    """

    model=1+1/np.power(snr, param)
    return model

#------------------------------------------------------------------------------------------------------------
def optBiasSeriesOffsetModelFunc(snr, params):
    """Optimization bias model function, of the form ``corrFactor = p0 + p1/x + p2/x**2 + ... + pn/x**n``
    where p0...pn are fit coefficents given as the params array. This is for use with the `fast` completeness
    method of the ``SelFn`` class.

    Args:
        snr (:obj:`np.ndarray`): Array of true signal-to-noise ratio values
        params (:obj:`np.ndarray`): Fit coefficients.

    Returns:
        Array of correction factors.

    """

    model=np.ones(snr.shape)*params[0]
    index=1
    for p in params[1:]:
        model=model+p/(snr**index)
        index=index+1
    return model

#------------------------------------------------------------------------------------------------------------
def _parseSourceInjectionData(injTab, inputTab, SNRCut):
    """Produce arrays for constructing interpolator objects from source injection test data.

    Args:
        injTab (:obj:`astropy.table.Table`): Output catalog produced by the source injection test.
        inputTab (:obj:`astropy.table.Table`): Input catalog produced by the source injection test.
        SNRCut (:obj:`float`): Selection threshold in S/N to apply.

    Returns:
        theta500s, ycBinCentres, compThetaGrid, thetaQ

    """

    # Completeness given y0 (NOT y0~) and theta500 and the S/N cut as 2D spline
    # We also derive survey-averaged Q here from the injection sim results [for y0 -> y0~ mapping]
    # NOTE: This is a survey-wide average, doesn't respect footprints at the moment
    # NOTE: This will need re-thinking for evolving, non-self-similar models?
    theta500s=np.unique(inputTab['theta500Arcmin'])
    binEdges=np.linspace(inputTab['inFlux'].min(), inputTab['inFlux'].max(), 101)
    binCentres=(binEdges[1:]+binEdges[:-1])/2
    compThetaGrid=np.zeros((theta500s.shape[0], binCentres.shape[0]))
    thetaQ=np.zeros(len(theta500s))
    #thetaQ_05=np.zeros(len(theta500s))
    #thetaQ_95=np.zeros(len(theta500s))
    for i in range(len(theta500s)):
        t=theta500s[i]
        injMask=np.logical_and(injTab['theta500Arcmin'] == t, injTab['SNR'] > SNRCut)
        inputMask=inputTab['theta500Arcmin'] == t
        injFlux=injTab['inFlux'][injMask]
        outFlux=injTab['outFlux'][injMask]
        inputFlux=inputTab['inFlux'][inputMask]
        recN, binEdges=np.histogram(injFlux, bins = binEdges)
        inpN, binEdges=np.histogram(inputFlux, bins = binEdges)
        valid=inpN > 0
        compThetaGrid[i][valid]=recN[valid]/inpN[valid]
        if len(outFlux) > 0:
            thetaQ[i]=np.median(outFlux/injFlux)
            #thetaQ_05[i]=np.percentile(outFlux/injFlux, 5)
            #thetaQ_95[i]=np.percentile(outFlux/injFlux, 95)

    return theta500s, binCentres, compThetaGrid, thetaQ

#------------------------------------------------------------------------------------------------------------
def loadAreaMask(tileName, selFnDir):
    """Loads the survey area mask, i.e., the area searched for sources and clusters, for the given tile.
    
    Args:
        tileName (:obj:`str`): The name of the tile for which the area mask will be loaded.
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
    
    Returns:
        Map array (2d :obj:`np.ndarray`), WCS object (:obj:`astWCS.WCS`)
    
    """
    
    areaMap, wcs=_loadTile(tileName, selFnDir, "areaMask", extension = 'fits')   
    return areaMap, wcs

#------------------------------------------------------------------------------------------------------------
def loadFlagMask(tileName, selFnDir):
    """Loads the flag mask, i.e., areas flagged for reasons such as point source subtraction, for the given
    tile.

    Args:
        tileName (:obj:`str`): The name of the tile for which the flag mask will be loaded.
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.

    Returns:
        Map array (2d :obj:`np.ndarray`), WCS object (:obj:`astWCS.WCS`)

    """

    areaMap, wcs=_loadTile(tileName, selFnDir, "flagMask", extension = 'fits')
    return areaMap, wcs

#------------------------------------------------------------------------------------------------------------
def loadRMSMap(tileName, selFnDir, photFilter):
    """Loads the RMS (noise) map for the given tile.
    
    Args:
        tileName (:obj:`str`): The name of the tile for which the RMS (noise) map will be loaded.
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
        photFilter (:obj:`str`): Name of the reference filter, as specified in, e.g., a :ref:`nemoCommand`
            config file (see :ref:`ConfigReference`).

    Returns:
        Map array (2d :obj:`np.ndarray`), WCS object (:obj:`astWCS.WCS`)
    
    """
    
    RMSMap, wcs=_loadTile(tileName, selFnDir, "RMSMap_%s" % (photFilter), extension = 'fits')   
    return RMSMap, wcs

#------------------------------------------------------------------------------------------------------------
def loadMassLimitMap(tileName, diagnosticsDir, z):
    """Loads the mass limit map for the given tile at the given redshift.

    Args:
        tileName (:obj:`str`): The name of the tile for which the mass limit map will be loaded.
        diagnosticsDir (:obj:`str`): Path to the ``diagnostics/`` directory, as produced by the 
            :ref:`nemoCommand` command.
        z (:obj:`float`): Redshift at which the mass limit map was made (should match an entry in the
            :ref:`nemoCommand` config file).

    Returns:
        Map array (2d :obj:`np.ndarray`), WCS object (:obj:`astWCS.WCS`)
    
    """
    
    massLimMap, wcs=_loadTile(tileName, diagnosticsDir, "massLimitMap_z%s" % (str(z).replace(".", "p")), 
                              extension = 'fits')   
    return massLimMap, wcs

#------------------------------------------------------------------------------------------------------------
def loadIntersectionMask(tileName, selFnDir, footprint):
    """Loads the intersection mask for the given tile and footprint.
    
    Args:
        tileName (:obj:`str`): The name of the tile for which the intersection mask will be loaded.
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
        footprint (:obj:`str`): The name of the footprint for which the intersection will be calculated,
            as defined in the :ref:`nemoCommand` config file (see :ref:`selFnFootprints`).
    
    Returns:
        Map array (2d :obj:`np.ndarray`), WCS object (:obj:`astWCS.WCS`)
    
    """           
            
    intersectMask, wcs=_loadTile(tileName, selFnDir, "intersect_%s" % (footprint), extension = 'fits')   
    return intersectMask, wcs

#------------------------------------------------------------------------------------------------------------
def _loadTile(tileName, baseDir, baseFileName, extension = 'fits'):
    """Generic function to load a tile image from either a multi-extension FITS file, or a file with 
    #tileName.fits type extension, whichever is found.
        
    Returns map array, wcs
    
    """
    
    # After tidyUp is run, this will be a MEF file... during first run, it won't be (written under MPI)
    if os.path.exists(baseDir+os.path.sep+"%s#%s.%s" % (baseFileName, tileName, extension)):
        fileName=baseDir+os.path.sep+"%s#%s.%s" % (baseFileName, tileName, extension)
    elif os.path.exists(baseDir+os.path.sep+tileName+os.path.sep+"%s#%s.%s" % (baseFileName, tileName, extension)):
        # Lustre-friendlier: each tile has its own directory
        fileName=baseDir+os.path.sep+tileName+os.path.sep+"%s#%s.%s" % (baseFileName, tileName, extension)
    else:
        fileName=baseDir+os.path.sep+"%s.%s" % (baseFileName, extension)
    with pyfits.open(fileName) as img:
        # If we find the tile - great. If not, we use first extension with data as it'll be compressed
        if tileName in img:
            extName=tileName
            data=img[extName].data
        else:
            data=None
        if data is None:
            for extName in img:
                data=img[extName].data
                if data is not None:
                    break
        data=img[extName].data
        wcs=astWCS.WCS(img[extName].header, mode = 'pyfits')
    
    return data, wcs

#------------------------------------------------------------------------------------------------------------
def getTileTotalAreaDeg2(tileName, selFnDir, masksList = [], footprintLabel = None):
    """Returns the total area of the tile given by `tileName` (taking into account masked regions).
    
    Args:
        tileName (:obj:`str`): The name of the tile for which the area will be calculated.
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
        masksList (:obj:`list`, optional): A list of paths to FITS-format mask images, which contain pixels
            with values of 1 to indicate valid survey area, 0 otherwise. If given, the area is calculated
            for the intersection of these masks with the survey area mask.
        footprintLabel (:obj:`str`, optional): The name of the footprint for which the intersection will
            be calculated, as defined in the :ref:`nemoCommand` config file (see :ref:`selFnFootprints`).
    
    Returns:
        Tile area, after masking, in square degrees.
    
    """
    
    areaMap, wcs=loadAreaMask(tileName, selFnDir)
    areaMapSqDeg=(maps.getPixelAreaArcmin2Map(areaMap.shape, wcs)*areaMap)/(60**2)
    totalAreaDeg2=areaMapSqDeg.sum()
    
    if footprintLabel is not None:  
        intersectMask=makeIntersectionMask(tileName, selFnDir, footprintLabel, masksList = masksList)
        totalAreaDeg2=(areaMapSqDeg*intersectMask).sum()        
        
    return totalAreaDeg2

#------------------------------------------------------------------------------------------------------------
def makeIntersectionMask(tileName, selFnDir, label, masksList = []):
    """Creates an intersection mask between the survey mask, and the mask files given in `masksList`.
    
    Args:
        tileName (:obj:`str`): The name of the tile for which the intersection mask will be made (or loaded
            from disk, if cached).
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
        label (:obj:`str`): The name of the footprint for which the intersection will be calculated,
            as defined in the :ref:`nemoCommand` config file (see :ref:`selFnFootprints`).
        masksList (:obj:`list`, optional): A list of paths to FITS-format mask images, which contain pixels
            with values of 1 to indicate valid survey area, 0 otherwise. If given, the area is calculated
            for the intersection of these masks with the survey area mask.
    
    Returns:
        Intersection mask (1 = valid area, 0 = outside of intersection area; 2d :obj:`np.ndarray`)
        
    Note:
        For speed, it is assumed that the declination axis is aligned with the vertical axis in each
        mask image. This routine caches the intersection masks in `selFnDir`.
            
    """
    
    # After tidyUp has run, there will be intersection mask MEF files
    intersectFileName=selFnDir+os.path.sep+"intersect_%s.fits" % (label)
    if os.path.exists(intersectFileName):
        intersectMask, wcs=loadIntersectionMask(tileName, selFnDir, label)
        return intersectMask
    
    # Otherwise, we may have a per-tile intersection mask
    intersectFileName=selFnDir+os.path.sep+tileName+os.path.sep+"intersect_%s#%s.fits" % (label, tileName)
    if os.path.exists(intersectFileName):
        intersectMask, wcs=loadIntersectionMask(tileName, selFnDir, label)
        return intersectMask

    # Otherwise... make it
    areaMap, wcs=loadAreaMask(tileName, selFnDir)
    RAMin, RAMax, decMin, decMax=wcs.getImageMinMaxWCSCoords()
    if masksList == []:
        raise Exception("didn't find previously cached intersection mask but makeIntersectionMask called with empty masksList")
    print("... creating %s intersection mask (%s) ..." % (label, tileName))         
    intersectMask=np.zeros(areaMap.shape)
    outRACoords=np.array(wcs.pix2wcs(np.arange(intersectMask.shape[1]), [0]*intersectMask.shape[1]))
    outDecCoords=np.array(wcs.pix2wcs([0]*np.arange(intersectMask.shape[0]), np.arange(intersectMask.shape[0])))
    outRA=outRACoords[:, 0]
    outDec=outDecCoords[:, 1]
    RAToX=interpolate.interp1d(outRA, np.arange(intersectMask.shape[1]), fill_value = 'extrapolate')
    DecToY=interpolate.interp1d(outDec, np.arange(intersectMask.shape[0]), fill_value = 'extrapolate')
    for fileName in masksList:
        with pyfits.open(fileName) as maskImg:
            for hdu in maskImg:
                if type(hdu) == pyfits.ImageHDU:
                    break
            maskWCS=astWCS.WCS(hdu.header, mode = 'pyfits')
            maskData=hdu.data
        # From sourcery tileDir stuff
        RAc, decc=wcs.getCentreWCSCoords()
        xc, yc=maskWCS.wcs2pix(RAc, decc)
        xc, yc=int(xc), int(yc)
        xIn=np.arange(maskData.shape[1])
        yIn=np.arange(maskData.shape[0])
        inRACoords=np.array(maskWCS.pix2wcs(xIn, [yc]*len(xIn)))
        inDecCoords=np.array(maskWCS.pix2wcs([xc]*len(yIn), yIn))
        inRA=inRACoords[:, 0]
        inDec=inDecCoords[:, 1]
        RAToX=interpolate.interp1d(inRA, xIn, fill_value = 'extrapolate')
        DecToY=interpolate.interp1d(inDec, yIn, fill_value = 'extrapolate')
        outRACoords=np.array(wcs.pix2wcs(np.arange(intersectMask.shape[1]), [0]*intersectMask.shape[1]))
        outDecCoords=np.array(wcs.pix2wcs([0]*np.arange(intersectMask.shape[0]), np.arange(intersectMask.shape[0])))
        outRA=outRACoords[:, 0]
        outDec=outDecCoords[:, 1]
        xIn=np.array(RAToX(outRA), dtype = int)
        yIn=np.array(DecToY(outDec), dtype = int)
        xMask=np.logical_and(xIn >= 0, xIn < maskData.shape[1])
        yMask=np.logical_and(yIn >= 0, yIn < maskData.shape[0])
        xOut=np.arange(intersectMask.shape[1])
        yOut=np.arange(intersectMask.shape[0])
        for i in yOut[yMask]:
            intersectMask[i][xMask]=intersectMask[i][xMask]+maskData[yIn[i], xIn[xMask]]
    intersectMask=np.array(np.greater(intersectMask, 0.5), dtype = int)
    maps.saveFITS(intersectFileName, intersectMask*areaMap, wcs, compressionType = 'PLIO_1')

    return intersectMask

#------------------------------------------------------------------------------------------------------------
def getRMSTab(tileName, photFilterLabel, selFnDir, footprintLabel = None, maxFlags = None):
    """Makes a table containing map area in the tile refered to by `tileName` against RMS (noise level)
    values, compressing the information in the RMS maps. Results are cached under `selFnDir`, and read from
    disk if found.
    
    Args:
        tileName (:obj:`str`): The name of the tile.
        photFilterLabel (:obj:`str`): Name of the reference filter, as specified in, e.g., a 
            :ref:`nemoCommand` config file (see :ref:`ConfigReference`).
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
        footprintLabel (:obj:`str`, optional): The name of the footprint in which the calculation will be
            done, as defined in the :ref:`nemoCommand` config file (see :ref:`selFnFootprints`).
    
    Returns:
        A table of RMS (noise level) values versus area in square degrees (:obj:`astropy.table.Table`).
    
    """
    
    # After tidyUp has run, we may have one global RMS table with an extra tileName column we can use
    RMSTabFileName=selFnDir+os.path.sep+"RMSTab.fits"
    if footprintLabel is not None:
        RMSTabFileName=RMSTabFileName.replace(".fits", "_%s.fits" % (footprintLabel))
    if maxFlags is not None:
        RMSTabFileName=RMSTabFileName.replace(".fits", "_maxFlags%d.fits" % (maxFlags))
    if os.path.exists(RMSTabFileName):
        tab=atpy.Table().read(RMSTabFileName)
        return tab[np.where(tab['tileName'] == tileName)]

    # Table doesn't exist, so make it...
    print(("... making RMS table for tile = %s, footprint = %s, maxFlags = %s" % (tileName, footprintLabel, str(maxFlags))))
    RMSMap, wcs=loadRMSMap(tileName, selFnDir, photFilterLabel)
    areaMap, wcs=loadAreaMask(tileName, selFnDir)
    if maxFlags is not None:
        flagMask, wcs=loadFlagMask(tileName, selFnDir)
        areaMap[flagMask > maxFlags]=0
    areaMapSqDeg=(maps.getPixelAreaArcmin2Map(areaMap.shape, wcs)*areaMap)/(60**2)

    if footprintLabel != None:  
        intersectMask=makeIntersectionMask(tileName, selFnDir, footprintLabel)
        areaMapSqDeg=areaMapSqDeg*intersectMask        
        RMSMap=RMSMap*intersectMask

    RMSValues=np.unique(RMSMap[np.nonzero(RMSMap)])
    totalAreaDeg2=areaMapSqDeg.sum()
    tileArea=np.zeros(len(RMSValues))
    for i in range(len(RMSValues)):
        tileArea[i]=areaMapSqDeg[np.equal(RMSMap, RMSValues[i])].sum()
    RMSTab=atpy.Table()
    RMSTab.add_column(atpy.Column(tileArea, 'areaDeg2'))
    RMSTab.add_column(atpy.Column(RMSValues, 'y0RMS'))
    # Checks - these should be impossible but we have seen (e.g., when messed up masks)
    tol=0.003
    if abs(RMSTab['areaDeg2'].sum()-areaMapSqDeg.sum()) > tol:
        raise Exception("Mismatch between area map and area in RMSTab for tile '%s'" % (tileName))
    if np.less(RMSTab['areaDeg2'], 0).sum() > 0:
        raise Exception("Negative area in tile '%s' - check your survey mask (and delete/remake tileDir files if necessary)." % (tileName))
    RMSTab.meta['NEMOVER']=nemo.__version__

    return RMSTab

#------------------------------------------------------------------------------------------------------------
def downsampleRMSTab(RMSTab, downsampleFactor = 2):
    """Downsamples `RMSTab` (see :meth:`getRMSTab`) in terms of noise resolution, binning by `stepSize`.
    
    Args:
        RMSTab (:obj:`astropy.table.Table`): An RMS table, as produced by :meth:`getRMSTab`.
        downsampleFactor (:obj:`float`, optional): Downsample the number of y0 bins by this factor.
            (e.g., downsample = 2 gives half the number of bins compared to the original table).
        
    Returns:
        A table of RMS (noise level) values versus area in square degrees (:obj:`astropy.table.Table`).        
    
    """
    
    if downsampleFactor < 1:
        raise Exception("downsampleFactor must be >= 1 - given %.1f" % (downsampleFactor))
    numBinEdges=int(len(RMSTab)/downsampleFactor)
    # We go slightly out of range here to allow overrideNoise option to work with changing anything else
    binEdges=np.linspace(RMSTab['y0RMS'].min()*0.999, RMSTab['y0RMS'].max()*1.001, numBinEdges+1)
    y0Binned=[]
    tileAreaBinned=[]
    binMins=[]
    binMaxs=[]
    for i in range(len(binEdges)-1):
        mask=np.logical_and(RMSTab['y0RMS'] >= binEdges[i], RMSTab['y0RMS'] < binEdges[i+1])
        if mask.sum() > 0:
            y0Binned.append(np.average(RMSTab['y0RMS'][mask], weights = RMSTab['areaDeg2'][mask]))
            tileAreaBinned.append(np.sum(RMSTab['areaDeg2'][mask]))
            binMins.append(binEdges[i])
            binMaxs.append(binEdges[i+1])
    newRMSTab=atpy.Table()
    newRMSTab.add_column(atpy.Column(y0Binned, 'y0RMS'))
    newRMSTab.add_column(atpy.Column(tileAreaBinned, 'areaDeg2'))
    RMSTab=newRMSTab
        
    return RMSTab

#------------------------------------------------------------------------------------------------------------
def calcTileWeightedAverageNoise(tileName, photFilterLabel, selFnDir, footprintLabel = None,
                                 maxFlags = None):
    """Returns the area weighted average ỹ\ :sub:`0` noise value in the tile.
    
    Args:
        tileName (:obj:`str`): The name of the tile.
        photFilterLabel (:obj:`str`): Name of the reference filter, as specified in, e.g., a 
            :ref:`nemoCommand` config file (see :ref:`ConfigReference`).
        selFnDir (:obj:`str`): Path to a ``selFn/`` directory, as produced by the :ref:`nemoCommand`
            command. This directory contains information such as the survey noise maps, area masks,
            and information needed to construct the filter mismatch function, `Q`, used in mass
            modeling.
        footprintLabel (:obj:`str`, optional): The name of the footprint in which the calculation will be
            done, as defined in the :ref:`nemoCommand` config file (see :ref:`selFnFootprints`).
    
    Returns:
        Area-weighted average noise in ỹ\ :sub:`0` (``fixed_y_c`` in Nemo cluster catalogs).
        
    """

    RMSTab=getRMSTab(tileName, photFilterLabel, selFnDir, footprintLabel = footprintLabel,
                     maxFlags = maxFlags)
    RMSValues=np.array(RMSTab['y0RMS'])
    tileArea=np.array(RMSTab['areaDeg2'])
    tileRMSValue=np.average(RMSValues, weights = tileArea)

    return tileRMSValue

#------------------------------------------------------------------------------------------------------------
def completenessByFootprint(config):
    """Write out the average (log\ :sub:`10` mass, z) grid over all survey footprints defined in the config,
    weighted by fraction of total survey area within the footprint. Also prints some useful statistics and
    produces some plots that are written to the `diagnosticsDir` directory.
    
    Args:
        config (:class:`nemo.startUp.NemoConfig`): A NemoConfig object.
    
    Returns:
        None
        
    Note:
        Output is written to files named e.g. ``diagnosticsDir/MzCompleteness_label.npz``, where `label`
        is the footprint name (i.e., a key in `selFnCollection`); 'full' is the default (survey-wide
        average).
    
    """

    zBinEdges=np.arange(0.05, 2.1, 0.1)
    zBinCentres=(zBinEdges[:-1]+zBinEdges[1:])/2.

    footprintLabels=[None]
    if 'selFnFootprints' in config.parDict.keys():
        for footprintDict in config.parDict['selFnFootprints']:
            footprintLabels.append(footprintDict['label'])

    for footprintLabel in footprintLabels:

        # We now support multiple relations, but here we take only the first entry
        scalingRelationDict=config.parDict['massOptions']['scalingRelations'][0]

        # Optimization bias can now be applied here
        biasModelDict=None
        if 'biasModel' in config.parDict['selFnOptions'].keys():
            if config.parDict['selFnOptions']['biasModel'] == 'power':
                biasModel=optBiasPowerModelFunc
            elif config.parDict['selFnOptions']['biasModel'] == 'series':
                biasModel=optBiasModelFunc
            elif config.parDict['selFnOptions']['biasModel'] == 'exp':
                biasModel=optBiasFuncExpModel
            else:
                raise Exception("biasModel must be 'series' or 'exp' - check selFnOptions in your config")
            try:
                biasModelParams=config.parDict['selFnOptions']['biasModelParams']
            except:
                raise Exception("If you specify biasModel, you must also give biasModelParams - check selFnOptions in your config")
            biasModelDict={'func': biasModel, 'params': biasModelParams}
            if footprintLabel == footprintLabels[0]:
                print(">>> Optimization bias model will be applied")

        try:
            selFn=SelFn(config.selFnDir, config.parDict['selFnOptions']['fixedSNRCut'],
                        footprint = footprintLabel, zStep = 0.02,
                        massBinsTheory = 200, zStepTheory = 0.01,
                        downsampleRMS = False,
                        applyRelativisticCorrection = scalingRelationDict['relativisticCorrection'],
                        delta = scalingRelationDict['delta'],
                        rhoType = scalingRelationDict['rhoType'],
                        method = config.parDict['selFnOptions']['method'],
                        QSource = config.parDict['selFnOptions']['QSource'],
                        maxFlags = config.parDict['selFnOptions']['maxFlags'],
                        biasModel = biasModelDict)
        except FootprintError:
            continue
            #print("... no overlapping area with footprint %s" % (footprintLabel))

        massLabel=selFn.mockSurvey.mdefLabel

        if footprintLabel is None:
            footprintLabel="full"

        outFileName=config.diagnosticsDir+os.path.sep+"MzCompleteness_%s_%s.npz" % (selFn.method, footprintLabel)
        np.savez(outFileName, z = selFn.z, log10M = selFn.log10M, completeness = selFn.compMz)
        makeMzCompletenessPlot(selFn.compMz, selFn.log10M, selFn.z, footprintLabel, massLabel,
                               config.diagnosticsDir+os.path.sep+"MzCompleteness_%s_%s.pdf" % (selFn.method, footprintLabel))

        # 90% mass completeness limit [we don't really need the plots any more - can use the Mz ones above]
        massLimit_90Complete=selFn.getMassLimit(0.9)
        # makeMassLimitVRedshiftPlot(massLimit_90Complete, selFn.z,
        #                            config.diagnosticsDir+os.path.sep+"completeness90Percent_%s_%s.pdf" % (selFn.method, footprintLabel),
        #                            title = "footprint: %s" % (footprintLabel))
        zMask=np.logical_and(selFn.z >= 0.2, selFn.z < 1.0)
        averageMassLimit_90Complete=np.average(massLimit_90Complete[zMask])
        print(">>> Survey-averaged results inside footprint %s [maxFlags = %s]:" % (footprintLabel, config.parDict['selFnOptions']['maxFlags']))
        print("... total survey area (after masking) = %.1f sq deg" % (selFn.totalAreaDeg2))
        print("... survey-averaged 90%% mass (%s) completeness limit (z = 0.5) = %.1f x 10^14 MSun" % (massLabel, massLimit_90Complete[np.argmin(abs(zBinCentres-0.5))]))

#------------------------------------------------------------------------------------------------------------
def calcCompletenessContour(compMz, log10M, z, level = 0.90):
    """Calculates a completeness contour on the (log\ :sub:`10` mass, z) plane.
    
    Args:
        compMz (:obj:`np.ndarray`): Map (2d array) of completeness on the (log\ :sub:`10` mass, z) plane.
        log10M (:obj:`np.ndarray`): One dimensional array of log\ :sub:`10` mass values corresponding to
            `compMz`.
        z (:obj:`np.ndarray`): One dimensional arra of redshifts corresponding to `compMz`.
        level (:obj:`float`, optional): Fractional completeness level (e.g., 0.90 is 90% completeness).
    
    Returns:
        Contour values for the given completeness level (a pair of arrays - redshifts, and log\ :sub:`10`
        mass values).
    
    """
    # Easiest way to get at contour for plotting later
    # The smoothing may only be necessary if compMz is made by montecarlo method
    contours=plt.contour(z, log10M, compMz.transpose(), levels = [level])
    cont_z=[]
    cont_log10M=[]
    for p in contours.collections[0].get_paths():
        v=p.vertices
        cont_z=cont_z+v[:, 0].tolist()
        cont_log10M=cont_log10M+v[:, 1].tolist()
    plt.close()
    contTab=atpy.Table()
    contTab.add_column(atpy.Column(cont_z, 'z'))
    contTab.add_column(atpy.Column(cont_log10M, 'log10M'))
    contTab.sort('z')
    cont_z=[]
    cont_log10M=[]
    for zi in z:
        mask=np.equal(contTab['z'], zi)
        if mask.sum() > 0:
            cont_z.append(zi)
            cont_log10M.append(np.median(contTab['log10M'][mask]))
    cont_z=np.array(cont_z)
    cont_log10M=np.array(cont_log10M)
    #cont_log10M=ndimage.uniform_filter(cont_log10M, int(np.ceil(len(cont_log10M)/20)))
    
    return cont_z, cont_log10M
    
#------------------------------------------------------------------------------------------------------------
def makeMzCompletenessPlot(compMz, log10M, z, title, massLabel, outFileName):
    """Makes a (log\ :sub:`10` mass, z) completeness plot.
    
    Args:
        compMz (:obj:`np.ndarray`): Map (2d array) of completeness on the (log\ :sub:`10` mass, z) plane.
        log10M (:obj:`np.ndarray`): One dimensional array of log\ :sub:`10` mass values corresponding to
            `compMz`.
        z (:obj:`np.ndarray`): One dimensional arra of redshifts corresponding to `compMz`.
        title (:obj:`str`): Title that will be written at the top of the plot.
        massLabel (:obj:`str`): Label for mass quantity (e.g., "M200c").
        outFileName (:obj:`str`): Path where the plot will be written to a file, with the format being
            determined by the file extension.
    
    Returns:
        None
        
    """
    
    cont_z, cont_log10M=calcCompletenessContour(compMz, log10M, z, level = 0.90)
    
    # Actual plot
    plotSettings.update_rcParams()
    fig, ax = plt.subplots(figsize=(9.5,6.5))

    plt.imshow((compMz*100).transpose(), cmap = colorcet.m_rainbow, origin = 'lower', aspect = 'auto')
    
    y_tck=interpolate.splrep(log10M, np.arange(log10M.shape[0]))
    plot_log10M=np.linspace(13.8, 15.4, 9)
    coords_log10M=interpolate.splev(plot_log10M, y_tck)
    labels_log10M=[]
    for lm in plot_log10M:
        labels_log10M.append("%.2f" % (lm))
    plt.yticks(interpolate.splev(plot_log10M, y_tck), labels_log10M)
    plt.ylim(coords_log10M.min(), coords_log10M.max())
    if massLabel[0] == "M":
        massLabel=massLabel[1:]
    plt.ylabel("log$_{10}$ ($M_{\\rm %s} / M_{\odot}$)" % (massLabel))
    
    x_tck=interpolate.splrep(z, np.arange(z.shape[0]))
    plot_z=np.linspace(0.0, 2.0, 11)
    coords_z=interpolate.splev(plot_z, x_tck)
    labels_z=[]
    for lz in plot_z:
        labels_z.append("%.1f" % (lz))
    plt.xticks(interpolate.splev(plot_z, x_tck), labels_z)
    plt.xlim(coords_z.min(), coords_z.max())
    plt.xlabel("$z$")
    
    coords_cont_z=interpolate.splev(cont_z, x_tck)
    coords_cont_log10M=interpolate.splev(cont_log10M, y_tck)
    plt.plot(coords_cont_z, coords_cont_log10M, 'k:', lw = 3)
    
    plt.colorbar(pad = 0.03)
    cbLabel="Completeness (%)" 
    plt.figtext(0.96, 0.52, cbLabel, ha="center", va="center", family = "sans-serif", rotation = "vertical")

    plt.tight_layout()

    if title != 'full':
        plt.title(title)    
    plt.savefig(outFileName)
    plt.close()
      
#------------------------------------------------------------------------------------------------------------
def makeMassLimitMapsAndPlots(config):
    """Makes maps of approximate mass completeness and associated plots. Output is written to the
    `diagnosticsDir` directory. The maps to make are controlled by the ``massLimitMaps`` parameter in the
    selFnOptions`` dictionary in the Nemo config file.
    
    Args:
        config (:class:`nemo.startUp.NemoConfig`): A NemoConfig object.

    Returns:
        None
    
    """


    # We now support multiple relations, but here we take only the first entry
    scalingRelationDict=config.parDict['massOptions']['scalingRelations'][0]

    # Optimization bias can now be applied here
    biasModelDict=None
    if 'biasModel' in config.parDict['selFnOptions'].keys():
        if config.parDict['selFnOptions']['biasModel'] == 'power':
            biasModel=optBiasPowerModelFunc
        elif config.parDict['selFnOptions']['biasModel'] == 'series':
            biasModel=optBiasModelFunc
        elif config.parDict['selFnOptions']['biasModel'] == 'exp':
            biasModel=optBiasFuncExpModel
        else:
            raise Exception("biasModel must be 'series' or 'exp' - check selFnOptions in your config")
        try:
            biasModelParams=config.parDict['selFnOptions']['biasModelParams']
        except:
            raise Exception("If you specify biasModel, you must also give biasModelParams - check selFnOptions in your config")
        biasModelDict={'func': biasModel, 'params': biasModelParams}
        # print(">>> Optimization bias model will be applied")

    selFn=SelFn(config.selFnDir, config.parDict['selFnOptions']['fixedSNRCut'],
                footprint = None, zStep = 0.1, setUpAreaMask = True,
                downsampleRMS = False,
                applyRelativisticCorrection = scalingRelationDict['relativisticCorrection'],
                delta = scalingRelationDict['delta'],
                rhoType = scalingRelationDict['rhoType'],
                method = config.parDict['selFnOptions']['method'],
                QSource = config.parDict['selFnOptions']['QSource'],
                maxFlags = config.parDict['selFnOptions']['maxFlags'],
                biasModel = biasModelDict)

    # Load in tiles RMS map and fill a mass limit tiled map, then stitch at the end
    for massLimDict in config.parDict['selFnOptions']['massLimitMaps']:
        t0=time.time()
        if 'completenessFraction' not in massLimDict.keys():
            completenessFraction=0.9
        else:
            completenessFraction=massLimDict['completenessFraction']
        zIndex=np.argmin(abs(selFn.z-massLimDict['z']))
        z=selFn.z[zIndex]
        outFileName=config.diagnosticsDir+os.path.sep+"massLimitMap_z%.2f_comp%.2f.fits" % (z, completenessFraction)
        inFileName=config.selFnDir+os.path.sep+"RMSMap_%s.fits" % (config.parDict['photFilter'])
        stitchedMapLimDict=maps.TileDict({}, tileCoordsDict = config.tileCoordsDict)
        with pyfits.open(inFileName) as img:
            count=0
            for tileName in config.allTileNames:
                count=count+1
                d=img[tileName].data
                noiseLevels=np.unique(d)
                noiseLevels=noiseLevels[noiseLevels > 0]
                for n in noiseLevels:
                    RMSTab=atpy.Table()
                    RMSTab['y0RMS']=[n]
                    RMSTab['areaDeg2']=[100] # doesn't matter, just for weighting
                    compMzTile, y0Grid=selFn.calcFastCompletenessInTile(tileName = tileName, return_y0Grid = True,
                                                                        RMSTab = RMSTab)
                    massLim=selFn.log10M[np.where(compMzTile[zIndex] > completenessFraction)[0][0]]
                    massLim=np.power(10, massLim)/1e14
                    d[d == n]=massLim
                stitchedMapLimDict[tileName]=d
        stitchedMapLimDict.saveStitchedFITS(outFileName, config.origWCS, compressionType = "RICE_1")
        t1=time.time()
        print("... made mass limit map '%s' (time taken = %.3f sec)" % (outFileName, t1-t0))

        # Plots
        plotSettings.update_rcParams()
        for key in list(stitchedMapLimDict.keys()):
            del stitchedMapLimDict[key]
        del stitchedMapLimDict
        with pyfits.open(outFileName) as img:
            for ext in img:
                if ext.data is not None:
                    massLimMap=np.float16(ext.data)
                    wcs=astWCS.WCS(ext.header, mode = 'pyfits')

        # Map plot [this is only set-up really for large area maps]
        colorBarPos='vertical'
        axesLabels="sexagesimal"
        fontSize=20.0
        figSize=(15.5, 4)
        axes=[0.08,0.07,0.88,0.96]
        cbShrink=0.8
        cbAspect=10
        fraction=0.05
        pad=0.01
        rotation=90
        labelpad=5
        # massLimMap=np.ma.masked_where(massLimMap <1e-6, massLimMap) # We run out of memory with this
        massLimMap[massLimMap < 1e-6]=np.nan
        cutLevels=[massLimMap[massLimMap > 0].min(), massLimMap[massLimMap > 0].max()]
        colorMapName=colorcet.m_rainbow
        fig=plt.figure(figsize = figSize)
        p=astPlots.ImagePlot(massLimMap, wcs, cutLevels = cutLevels, title = None, axes = axes,
                             axesLabels = axesLabels, colorMapName = colorMapName, axesFontFamily = 'sans-serif',
                             RATickSteps = {'deg': 30.0, 'unit': 'h'}, decTickSteps = {'deg': 20.0, 'unit': 'd'},
                             axesFontSize = fontSize)
        cbLabel="$M_{\\rm %d%s}$ (10$^{14}$ M$_{\odot}$)\n[%d%% complete]" % (selFn.mockSurvey.delta, selFn.mockSurvey.rhoType[0], int(100*completenessFraction))
        cb=plt.colorbar(p.axes.images[0], ax = p.axes, orientation = colorBarPos, fraction = fraction, pad = pad,
                shrink = cbShrink, aspect = cbAspect)
        cb.ax.get_yaxis().labelpad=labelpad
        cb.ax.set_ylabel(cbLabel, rotation = rotation)
        plt.savefig(outFileName.replace(".fits", ".pdf"), dpi = 300)
        plt.savefig(outFileName.replace(".fits", ".png"), dpi = 300)
        plt.close()
        del p
        del selFn

        # Cumulative area info [note, compression drives up number of 'mass limits', so we rebin]
        pixAreaMap=maps.getPixelAreaArcmin2Map(massLimMap.shape, wcs)
        pixAreaMap=pixAreaMap/(60**2) # Have to do this here to avoid overflow if 16 bit
        limits=np.unique(massLimMap[massLimMap > 0])
        limits=limits[limits > 0]
        binEdges=np.linspace(limits.min(), limits.max(), 30)
        binCentres=(binEdges[1:]+binEdges[:-1])/2
        areas=np.zeros(len(binCentres))
        for i in range(len(binEdges)-1):
            binMin=binEdges[i]
            binMax=binEdges[i+1]
            mask=np.logical_and(massLimMap > binMin, massLimMap <= binMax)
            if mask.sum() > 0:
                areas[i]=pixAreaMap[mask].sum()
        tab=atpy.Table()
        tab['MLim']=binCentres
        tab['areaDeg2']=areas
        tab.sort('MLim')
        del pixAreaMap, massLimMap, limits

        # Full survey cumulative area plot
        plt.figure(figsize=(9,6.5))
        ax=plt.axes([0.155, 0.12, 0.82, 0.86])
        plt.minorticks_on()
        plt.plot(tab['MLim'], np.cumsum(tab['areaDeg2']), 'k-')
        #plt.plot(plotMRange, plotCumArea, 'k-')
        # As in DR5 paper
        # plt.ylabel("survey area < $M_{\\rm %d%s}$ limit (deg$^2$)" % (scalingRelationDict['delta'],
        #                                                               scalingRelationDict['rhoType'][0]))
        # plt.xlabel("$M_{\\rm %d%s}$ (10$^{14}$ M$_{\odot}$) [%d%% complete]" % (scalingRelationDict['delta'],
        #                                                                         scalingRelationDict['rhoType'][0],
        #                                                                         int(100*completenessFraction)))
        # Suggested by Arthur in DR6 review process
        plt.ylabel("sky area (deg$^2$)")
        plt.xlabel("%d%% completeness limit $M_{\\rm %d%s}$ (10$^{14}$ M$_{\odot}$)" % (int(100*completenessFraction),
                                                                                        scalingRelationDict['delta'],
                                                                                        scalingRelationDict['rhoType'][0]))
        labelStr="total survey area = %.0f deg$^2$" % (np.cumsum(tab['areaDeg2']).max())
        plt.ylim(0.0, 1.2*np.cumsum(tab['areaDeg2']).max())
        plt.xlim(tab['MLim'].min(), tab['MLim'].max())
        plt.figtext(0.2, 0.9, labelStr, ha="left", va="center")
        plt.savefig(config.diagnosticsDir+os.path.sep+"cumulativeArea_massLimit_z%.2f_comp%.2f.pdf" % (z, completenessFraction))
        plt.savefig(config.diagnosticsDir+os.path.sep+"cumulativeArea_massLimit_z%.2f_comp%.2f.png" % (z, completenessFraction))
        plt.close()

        # Deepest 20% cumulative area plot - we show a bit beyond this
        totalAreaDeg2=tab['areaDeg2'].sum()
        deepTab=tab[np.where(np.cumsum(tab['areaDeg2']) < 0.25 * totalAreaDeg2)]
        plt.figure(figsize=(9,6.5))
        ax=plt.axes([0.155, 0.12, 0.82, 0.86])
        plt.minorticks_on()
        plt.plot(tab['MLim'], np.cumsum(tab['areaDeg2']), 'k-')
        # plt.ylabel("survey area < $M_{\\rm %d%s}$ limit (deg$^2$)" % (scalingRelationDict['delta'],
        #                                                               scalingRelationDict['rhoType'][0]))
        # plt.xlabel("$M_{\\rm %d%s}$ (10$^{14}$ M$_{\odot}$) [%d%% complete]" % (scalingRelationDict['delta'],
        #                                                                         scalingRelationDict['rhoType'][0],
        #                                                                         int(100*completenessFraction)))
        # Suggested by Arthur in DR6 review process
        plt.ylabel("sky area (deg$^2$)")
        plt.xlabel("%d%% completeness limit $M_{\\rm %d%s}$ (10$^{14}$ M$_{\odot}$)" % (int(100*completenessFraction),
                                                                                        scalingRelationDict['delta'],
                                                                                        scalingRelationDict['rhoType'][0]))
        labelStr="area of deepest 20%% = %.0f deg$^2$" % (0.2 * totalAreaDeg2)
        plt.ylim(0.0, 1.2*np.cumsum(deepTab['areaDeg2']).max())
        plt.xlim(deepTab['MLim'].min(), deepTab['MLim'].max())
        plt.figtext(0.2, 0.9, labelStr, ha="left", va="center")
        plt.savefig(config.diagnosticsDir+os.path.sep+"cumulativeArea_massLimit_z%.2f_comp%.2f_deepest20percent.pdf" % (z, completenessFraction))
        plt.savefig(config.diagnosticsDir+os.path.sep+"cumulativeArea_massLimit_z%.2f_comp%.2f_deepest20percent.png" % (z, completenessFraction))
        plt.close()

#------------------------------------------------------------------------------------------------------------
def makeMassLimitVRedshiftPlot(massLimit_90Complete, zRange, outFileName, title = None):
    """Makes a plot of 90%-completeness mass limit versus redshift. Uses spline interpolation.
    
    Args:
        massLimit_90Complete (:obj:`np.ndarray`): Mass limit at each redshift, corresponding to 90%
            completeness.
        zRange (:obj:`np.ndarray`): Redshifts at which mass completeness was evaluted.
        outFileName (:obj:`str`): Path to which the plot file will be written. The format is determined by
            the file extension.
        title (:obj:`str`): The title that will be written at the top of the plot.
    
    Returns:
        None
    
    """
    
    plotSettings.update_rcParams()
    plt.figure(figsize=(9,6.5))
    ax=plt.axes([0.10, 0.11, 0.87, 0.86])
    if title is not None:
        plt.figtext(0.15, 0.2, title, ha="left", va="center")
    tck=interpolate.splrep(zRange, massLimit_90Complete)
    plotRange=np.linspace(0, 2, 100)
    plt.plot(plotRange, interpolate.splev(plotRange, tck), 'k-')
    plt.plot(zRange, massLimit_90Complete, 'D', ms = 8)
    plt.xlabel("$z$")
    plt.ylim(0.5, 8)
    plt.xticks(np.arange(0, 2.2, 0.2))
    plt.xlim(0, 2)
    labelStr="$M_{\\rm 500c}$ (10$^{14}$ M$_{\odot}$) [90% complete]"
    plt.ylabel(labelStr)
    plt.savefig(outFileName)
    if outFileName.find(".pdf") != -1:
        plt.savefig(outFileName.replace(".pdf", ".png"))
    plt.close()   
    
        
#------------------------------------------------------------------------------------------------------------
def makeFullSurveyMassLimitMapPlot(z, config):
    """Makes full area mass limit map by reprojecting the tile mass limit maps onto the full map pixelisation.
    Creates both a plot and a FITS image.
    
    Args:
        z (:obj:`float`): The redshift at which the mass completeness is evaluated.
        config (:class:`nemo.startUp.NemoConfig`): A NemoConfig object.
    
    Note:
        Output is written to the `diagnostics` directory, as produced by the :ref:`nemoCommand` command.
    
    Returns:
        None
    
    """
    
    if 'makeQuickLookMaps' not in config.parDict.keys():
        config.quicklookScale=0.25
        config.quicklookShape, config.quicklookWCS=maps.shrinkWCS(config.origShape, config.origWCS, config.quicklookScale)

    outFileName=config.diagnosticsDir+os.path.sep+"reproj_massLimitMap_z%s.fits" % (str(z).replace(".", "p"))
    maps.stitchTilesQuickLook(config.diagnosticsDir+os.path.sep+"*"+os.path.sep+"massLimitMap_z%s#*.fits" % (str(z).replace(".", "p")),
                              outFileName, config.quicklookWCS, config.quicklookShape,
                              fluxRescale = config.quicklookScale)

    # Make plot
    if os.path.exists(outFileName) == True:
        with pyfits.open(outFileName) as img:
            for hdu in img:
                if hdu.shape != ():
                    reproj=np.nan_to_num(hdu.data)
                    reproj=np.ma.masked_where(reproj <1e-6, reproj)
                    wcs=astWCS.WCS(hdu.header, mode = 'pyfits')
            plotSettings.update_rcParams()
            fontSize=20.0
            figSize=(16, 5.7)
            axesLabels="sexagesimal"
            axes=[0.08,0.15,0.91,0.88]
            cutLevels=[2, int(np.median(reproj[np.nonzero(reproj)]))+2]
            colorMapName=colorcet.m_rainbow
            fig=plt.figure(figsize = figSize)
            p=astPlots.ImagePlot(reproj, wcs, cutLevels = cutLevels, title = None, axes = axes, 
                                axesLabels = axesLabels, colorMapName = colorMapName, axesFontFamily = 'sans-serif', 
                                RATickSteps = {'deg': 30.0, 'unit': 'h'}, decTickSteps = {'deg': 20.0, 'unit': 'd'},
                                axesFontSize = fontSize)
            cbLabel="$M_{\\rm 500c}$ (10$^{14}$ M$_{\odot}$) [90% complete]"
            cbShrink=0.7
            cbAspect=40
            cb=plt.colorbar(p.axes.images[0], ax = p.axes, orientation="horizontal", fraction = 0.05, pad = 0.18, 
                            shrink = cbShrink, aspect = cbAspect)
            plt.figtext(0.53, 0.04, cbLabel, ha="center", va="center", fontsize = fontSize, family = "sans-serif")
            plt.savefig(outFileName.replace(".fits", ".pdf"), dpi = 300)
            plt.savefig(outFileName.replace(".fits", ".png"), dpi = 300)
            plt.close()

#------------------------------------------------------------------------------------------------------------
def tidyUp(config):
    """Tidies up the `selFn` directory, constructing multi-extension FITS files from individual tile images
    and tables, deleting the individual tile files when complete. This routine also copies the given Nemo
    configuration into the `selFn` directory, and writes a plain text file that lists the tile names, and the
    areas covered by each tile.
    
    Args:
        config (:class:`nemo.startUp.NemoConfig`): A NemoConfig object.
        
    Returns:
        None
    
    """
    
    shutil.copy(config.configFileName, config.selFnDir+os.path.sep+"config.yml")

    # Make MEFs
    MEFsToBuild=["RMSMap_%s" % (config.parDict['photFilter'])]
    compressionTypes=["RICE_1"]
    dtypes=[np.float32]
    if 'selFnFootprints' in config.parDict.keys():
        for footprintDict in config.parDict['selFnFootprints']:
            MEFsToBuild.append("intersect_%s" % footprintDict['label'])
            compressionTypes.append("PLIO_1")
            dtypes.append(np.uint8)
    assert(len(MEFsToBuild) == len(compressionTypes))
    assert(len(MEFsToBuild) == len(dtypes))
    for MEFBaseName, compressionType, dtype in zip(MEFsToBuild, compressionTypes, dtypes):
        outFileName=config.selFnDir+os.path.sep+MEFBaseName+".fits"
        newImg=pyfits.HDUList()
        filesToRemove=[]
        for tileName in config.allTileNames:
            fileName=config.selFnDir+os.path.sep+tileName+os.path.sep+MEFBaseName+"#"+tileName+".fits"
            if os.path.exists(fileName):
                with pyfits.open(fileName) as img:
                    for extName in img:
                        if img[extName].data is not None:
                            break
                    hdu=pyfits.CompImageHDU(np.array(img[extName].data, dtype = dtype), img[extName].header, 
                                            name = tileName, compression_type = compressionType)
                filesToRemove.append(fileName)
                newImg.append(hdu)
        if len(newImg) > 0:
            newImg.writeto(outFileName, overwrite = True)
            for f in filesToRemove:
                os.remove(f)

    # Write a table of tile areas for those that want it
    with open(config.selFnDir+os.path.sep+"tileAreas.txt", "w") as outFile:
        outFile.write("#tileName areaDeg2\n")
        for tileName in config.allTileNames:
            tileAreaDeg2=getTileTotalAreaDeg2(tileName, config.selFnDir)
            outFile.write("%s %.6f\n" % (tileName, tileAreaDeg2))

    # Clean out any per-tile directories
    for tileName in config.allTileNames:
        if os.path.exists(config.selFnDir+os.path.sep+tileName) == True:
            shutil.rmtree(config.selFnDir+os.path.sep+tileName)


    
