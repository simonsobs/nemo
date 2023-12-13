"""

This module contains tools for manipulating maps.

"""

from astLib import *
from scipy import ndimage
from scipy import interpolate
from scipy.signal import convolve as scipy_convolve
from scipy import optimize
import astropy.io.fits as pyfits
import astropy.table as atpy
import astropy.stats as apyStats
import mahotas
import colorcet
import numpy as np
import pylab as plt
import glob
import os
import sys
import math
import time
import shutil
import copy
import yaml
import pickle
from pixell import enmap, curvedsky, utils, powspec
import nemo
try:
    import reproject
except:
    pass
from . import catalogs
from . import signals
from . import photometry
from . import plotSettings
from . import pipelines
from . import completeness
np.random.seed()

#------------------------------------------------------------------------------------------------------------
class MapDict(dict):
    """A dictionary for managing a sky map (a 2d array with an associated WCS) within Nemo. Keys within the
    dictionary can be set to values that control preprocessing of the map (usually done before filtering).
    Many of the keys in the dictionary map to entries in the .yml config file used by Nemo.

    Args:
        inputDict (:obj:`dict`): Input dictionary (usually this mirrors the contents of `unfilteredMaps` in
            Nemo .yml config files).
        tileCoordsDict (:obj:`dict`, optional): A dictionary that describes the tiling of a large map, as
            produced by :meth:`startUp.NemoConfig.getTileCoordsDict`.

    Attributes:
        tileCoordsDict (:obj:`dict`): A dictionary that describes the tiling of a large map, as
            produced by :meth:`startUp.NemoConfig.getTileCoordsDict`.
        validMapKeys (:obj:`list`): A list of keys that may contain a path to a map in FITS image format.
            These are: ['mapFileName', 'weightsFileName', 'pointSourceMask', 'surveyMask', 'flagMask'].

    """

    def __init__(self, inputDict, tileCoordsDict = None):
        super(MapDict, self).__init__(inputDict)
        self.tileCoordsDict=tileCoordsDict
        self._maskKeys=['pointSourceMask', 'surveyMask', 'flagMask', 'extendedMask']
        self.validMapKeys=['mapFileName', 'weightsFileName']+self._maskKeys


    def copy(self):
        """Make a copy of this :class:`MapDict` object.

        Returns:
            A deep copy of the :class:`MapDict` object.

        """
        return MapDict(self, tileCoordsDict = self.tileCoordsDict)


    def loadTile(self, mapKey, tileName, returnWCS = False):
        """Given a key in the map dictionary that contains the path to a FITS image, return the map
        as a 2d array and (optionally) the WCS.

        Args:
            mapKey (:obj:`str`): Name of a key in a map dictionary that contains the path to map FITS image.
                See self.`validMapKeys` for a list.
            tileName (:obj:`str`): The name of the tile to load.

        Returns:
            Map data as a 2d array (and optionally a WCS)

        Note:
            Tiles can be re-projected from CAR to TAN on the fly if the 'reprojectToTan' is set in the
            Nemo config.

        """

        if mapKey not in self.validMapKeys:
            raise Exception("mapKey must be one of %s - given mapKey = '%s'." % (self.validMapKeys, mapKey))

        pathToTileImages=self.get(mapKey)
        if os.path.isdir(pathToTileImages) == True:
            # Directory full of tile images (used by, e.g., on-the-fly extended source masking)
            with pyfits.open(pathToTileImages+os.path.sep+tileName+".fits") as img:
                extName=0
                tileData=img[extName].data
                if tileData is None:
                    for extName in img:
                        tileData=img[extName].data
                        if tileData is not None:
                            break
                assert tileData is not None
                if returnWCS == True or self['reprojectToTan'] == True:
                    # Zapping keywords in old ACT maps that confuse astropy.wcs
                    wcs=astWCS.WCS(img[extName].header, mode = 'pyfits', zapKeywords = ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'])
                data=tileData
        elif type(pathToTileImages) == np.ndarray:
            # We no longer want to support this kind of thing... clean this up later
            raise Exception("Expected a path but got an array instead (image already loaded).")
        else:
            # On-the-fly tile clipping
            with pyfits.open(pathToTileImages) as img:
                for ext in img:
                    if img[ext].data is not None:
                        break
                if returnWCS == True or self['reprojectToTan'] == True:
                    wcs=astWCS.WCS(self.tileCoordsDict[tileName]['header'], mode = 'pyfits')
                minX, maxX, minY, maxY=self.tileCoordsDict[tileName]['clippedSection']
                if img[ext].data.ndim == 3:
                    data=img[ext].data[0, minY:maxY, minX:maxX]
                elif img[ext].data.ndim == 2:
                    data=img[ext].data[minY:maxY, minX:maxX]
                else:
                    raise Exception("Map data has %d dimensions - only ndim = 2 or ndim = 3 are currently handled." % (img[ext].data.ndim))
                # Avoiding potential for read-only weirdness
                data=data[:]
                data=data.copy()

        # Convert any mask to 8-bit unsigned ints to save memory
        if mapKey in self._maskKeys:
            if data.dtype != np.uint8:
                data=np.array(data, dtype = np.uint8)

        # Survey masks are special: we need to zap the border overlap area or area calculations will be wrong
        if mapKey == 'surveyMask':
            minX, maxX, minY, maxY=self.tileCoordsDict[tileName]['areaMaskInClipSection']
            data[:minY, :]=0
            data[maxY:, :]=0
            data[:, :minX]=0
            data[:, maxX:]=0

        # Optional TAN reprojection - may help avoid biases due to distortion at high dec in CAR
        # WARNING: Probably introduces a new pixel window if we're not careful
        if self['reprojectToTan'] == True:
            if mapKey in self._maskKeys:
                order=0
            else:
                order='bicubic'
            tanWCS=_makeTanWCS(wcs)
            ySizePix, xSizePix=tanWCS.header['NAXIS2'], tanWCS.header['NAXIS1']
            tanData, footprint=reproject.reproject_interp((data, wcs.AWCS), tanWCS.AWCS, shape_out = [ySizePix, xSizePix],
                                                          order = order, return_footprint = True)
            tanData[footprint == 0]=0 # get rid of nans which will be in borders anyway
            # checkData=reproject.reproject_interp((tanData, tanWCS.AWCS), wcs.AWCS, shape_out = data.shape, order = 'bicubic',
                                                 # return_footprint = False)
            wcs=tanWCS
            data=tanData

        if returnWCS == True:
            return data, wcs
        else:
            return data


    def preprocess(self, tileName = 'PRIMARY', diagnosticsDir = None):
        """Applies a number of pre-processing steps to the map described by this :class:`MapDict` object,
        typically used before filtering.

        The first step is to load the map itself and the associated weights. Some other operations that may be
        applied are controlled by keys added to the MapDict. Some of these may be specified in the .yml
        configuration file, while others are applied by particular filter objects or by routines that generate
        simulated data. The following keys are understood:

        surveyMask (:obj:`str`)
            Path to a mask (.fits image; 1 = valid, 0 = masked) that defines the valid object search area.

        pointSourceMask (:obj:`str`)
            Path to a mask (.fits image; 1 = valid, 0 = masked) that contains holes at the locations of point
            sources, defining regions that are excluded from the object search area.

        RADecSection (:obj:`list`)
            Defines a region to extract from the map. Use the format [RAMin, RAMax, decMin, decMax] (units:
            decimal degrees).

        CMBSimSeed (:obj:`int`)
            If present, replace the map with a source-free simulated CMB realisation, generated using the given
            seed number. Used by :meth:`estimateContaminationFromSkySim`.

        applyBeamConvolution (:obj:`bool`)
            If True, the map is convolved with the beam given in the beamFileName key. This should only be
            needed when using preliminary y-maps made by tILe-C.

        Args:
            tileName (:obj:`str`): Name of the map tile (extension name) to operate on.
            diagnosticsDir (:obj:`str`): Path to a directory where miscellaneous diagnostic data may be written.

        Returns:
            None - values in the map dictionary are updated in-place, and additional keys may be added.

        """

        data, wcs=self.loadTile('mapFileName', tileName, returnWCS = True)

        # Optional calibration factor
        if 'calibFactor' in self.keys():
            data=data*self['calibFactor']

        if self['units'] == 'Jy/sr':
            if self['obsFreqGHz'] == 148:
                data=(data/1.072480e+09)*2.726*1e6
            elif self['obsFreqGHz'] == 219:
                data=(data/1.318837e+09)*2.726*1e6
            else:
                raise Exception("no code added to support conversion to uK from Jy/sr for freq = %.0f GHz" \
                        % (self['obsFreqGHz']))

        # Load weight map if given
        if 'weightsFileName' in list(self.keys()) and self['weightsFileName'] is not None:
            weights=self.loadTile('weightsFileName', tileName)
            # For Enki maps... take only I (temperature) for now, add options for this later
            if weights.ndim == 3:       # I, Q, U
                weights=weights[0, :]
            elif weights.ndim == 4:     # I, Q, U and also a covariance matrix
                weights=weights[0, 0, :]
        else:
            weights=np.ones(data.shape)

        # We rely on pixels with zero weight having zero value in actual maps later (automated edge trimming)
        # This might not be the case if the map has been filtered slightly before being fed into nemo
        data[weights == 0]=0

        # Load survey and point source masks, if given
        if 'surveyMask' in list(self.keys()) and self['surveyMask'] is not None:
            surveyMask=self.loadTile('surveyMask', tileName)
        else:
            surveyMask=np.ones(data.shape, dtype = np.uint8)
            surveyMask[weights == 0]=0

        # Some apodisation of the data outside the survey mask
        # NOTE: should add adjustable parameter for this somewhere later
        if 'apodizeUsingSurveyMask' in list(self.keys()) and self['apodizeUsingSurveyMask'] == True:
            # We need to remain unapodized to at least noiseGridArcmin beyond the edge of the survey mask
            # We'll need to make these adjustable parameters
            apodMask=np.array(surveyMask, dtype = bool)
            for i in range(120):
                apodMask=mahotas.dilate(apodMask)
            apodMask=ndimage.gaussian_filter(np.array(apodMask, dtype = float), 20)
            data=data*apodMask
            del apodMask

        if 'pointSourceMask' in list(self.keys()) and self['pointSourceMask'] is not None:
            psMask=self.loadTile('pointSourceMask', tileName)
        else:
            psMask=np.ones(data.shape, dtype = np.uint8)

        # Use for tracking regions where subtraction/in-painting took place to make flags in catalog
        # We can also supply a flag mask at the start, e.g., for marking dusty regions without zapping them
        # NOTE: flag masks for each frequency map get combined within filter objects
        if 'flagMask' in list(self.keys()) and self['flagMask'] is not None:
            flagMask=self.loadTile('flagMask', tileName)*surveyMask
        else:
            flagMask=np.zeros(data.shape, dtype = np.uint8)

        # Optional map clipping
        if 'RADecSection' in list(self.keys()) and self['RADecSection'] is not None:
            RAMin, RAMax, decMin, decMax=self['RADecSection']
            clip=astImages.clipUsingRADecCoords(data, wcs, RAMin, RAMax, decMin, decMax)
            data=clip['data']
            whtClip=astImages.clipUsingRADecCoords(weights, wcs, RAMin, RAMax, decMin, decMax)
            weights=whtClip['data']
            psClip=astImages.clipUsingRADecCoords(psMask, wcs, RAMin, RAMax, decMin, decMax)
            psMask=psClip['data']
            surveyClip=astImages.clipUsingRADecCoords(surveyMask, wcs, RAMin, RAMax, decMin, decMax)
            surveyMask=surveyClip['data']
            flagClip=astImages.clipUsingRADecCoords(flagMask, wcs, RAMin, RAMax, decMin, decMax)
            flagMask=flagClip['data']
            wcs=clip['wcs']
            if len(clip['data']) == 0:
                raise Exception("Clipping using RADecSection returned empty array - check RADecSection in config .yml file is in map")

        # For source-free simulations (contamination tests)
        if 'CMBSimSeed' in list(self.keys()):
            randMap=simCMBMap(data.shape, wcs, noiseLevel = 0, beam = self['beamFileName'],
                              seed = self['CMBSimSeed'])
            randMap[np.equal(weights, 0)]=0
            # Add white noise that varies according to inv var map...
            # Noise needed is the extra noise we need to add to match the real data, scaled by inv var map
            # This initial estimate is too high, so we use a grid search to get a better estimate
            mask=np.nonzero(data)
            dataSigma=data[mask].std()
            whiteNoiseLevel=np.zeros(weights.shape)
            whiteNoiseLevel[mask]=1/np.sqrt(weights[mask])
            noiseNeeded=np.sqrt(data[mask].var()-randMap[mask].var()-np.median(whiteNoiseLevel[mask])**2)
            noiseBoostFactor=noiseNeeded/np.median(whiteNoiseLevel[mask])
            # NOTE: disabled finding boost factor below for now...
            bestBoostFactor=1.
            # --- disabled
            #bestDiff=1e6
            #bestBoostFactor=noiseBoostFactor
            #simNoiseValues=simNoise[mask]
            #for boostFactor in np.linspace(noiseBoostFactor*0.5, noiseBoostFactor, 10):
                #diff=abs(dataSigma-(simNoiseValues+generatedNoise*boostFactor).std())
                #if diff < bestDiff:
                    #bestBoostFactor=boostFactor
                    #bestDiff=diff
            # ---
            data[mask]=np.random.normal(randMap[mask], bestBoostFactor*whiteNoiseLevel[mask],
                                        whiteNoiseLevel[mask].shape)
            outFileName=diagnosticsDir+os.path.sep+"CMBSim_%d#%s.fits" % (self['obsFreqGHz'], tileName)
            saveFITS(outFileName, data, wcs)

        # For position recovery tests, completeness calculations
        if 'injectSources' in list(self.keys()):
            # NOTE: Need to add varying GNFWParams here
            if 'GNFWParams' in self['injectSources'].keys():
                GNFWParams=self['injectSources']['GNFWParams']
                obsFreqGHz=self['obsFreqGHz']
            else:
                GNFWParams=None
                obsFreqGHz=None
            # Unsure if we actually want/need the below...
            # source injection sim tiles are processed independently (so shouldn't be double counting in overlaps anyway)
            validAreaSection=self.tileCoordsDict[tileName]['areaMaskInClipSection']
            modelMap=makeModelImage(data.shape, wcs, self['injectSources']['catalog'],
                                    self['beamFileName'], obsFreqGHz = obsFreqGHz,
                                    GNFWParams = GNFWParams, profile = self['injectSources']['profile'],
                                    validAreaSection = validAreaSection,
                                    override = self['injectSources']['override'])
            if modelMap is not None:
                modelMap[weights == 0]=0
                data=data+modelMap

        # Should only be needed for handling preliminary tILe-C maps
        if 'applyBeamConvolution' in self.keys() and self['applyBeamConvolution'] == True:
            data=convolveMapWithBeam(data, wcs, self['beamFileName'], maxDistDegrees = 1.0)
            if diagnosticsDir is not None:
                saveFITS(diagnosticsDir+os.path.sep+"beamConvolved#%s.fits" % (tileName), data, wcs)

        # Smoothing with some kernel (used, e.g., in PSF-matching between maps in nemoSpec)
        if 'smoothKernel' in self.keys():
            if 'smoothAttenuationFactor' in self.keys():
                data=data*self['smoothAttenuationFactor']
            data=convolveMapWithBeam(data, wcs, self['smoothKernel'], maxDistDegrees = 1.0)

        # Check if we're going to need to fill holes - if so, set-up smooth background only once
        # NOTE: If this needs changing, needs parametrizing as used in e.g. ACT DR5 results
        holeFillingKeys=['maskPointSourcesFromCatalog', 'maskAndFillFromCatalog', 'extendedMask']
        holeFilling=False
        for h in holeFillingKeys:
            if h in list(self.keys()):
                holeFilling=True
                break
        if holeFilling == True:
            pixRad=(10.0/60.0)/wcs.getPixelSizeDeg()
            bckData=ndimage.median_filter(data, int(pixRad))

        if 'extendedMask' in list(self.keys()):
            # Filling with white noise + smooth large scale image
            # WARNING: Assumes weights are ivar maps [true for Sigurd's maps]
            extendedMask=self.loadTile('extendedMask', tileName = tileName)
            mask=np.nonzero(weights)
            whiteNoiseLevel=np.zeros(weights.shape)
            whiteNoiseLevel[mask]=1/np.sqrt(weights[mask])
            data[extendedMask == 1]=bckData[extendedMask == 1]+np.random.normal(0, whiteNoiseLevel[extendedMask == 1])
            surveyMask=surveyMask*(1-extendedMask)
            #flagMask=flagMask+extendedMask

        # Optional masking of point sources from external catalog
        # Especially needed if using Fourier-space matched filter (and maps not already point source subtracted)
        if 'maskPointSourcesFromCatalog' in list(self.keys()) and self['maskPointSourcesFromCatalog'] is not None:
            # This is fast enough if using small tiles and running in parallel...
            # If our masking/filling is effective enough, we may not need to mask so much here...
            if type(self['maskPointSourcesFromCatalog']) is not list:
                self['maskPointSourcesFromCatalog']=[self['maskPointSourcesFromCatalog']]
            psMask=np.ones(data.shape, dtype = np.uint8)
            #pixRad=(10.0/60.0)/wcs.getPixelSizeDeg()
            #bckData=ndimage.median_filter(data, int(pixRad))
            rDegMap=np.ones(data.shape, dtype = float)*1e6
            for catalogInfo in self['maskPointSourcesFromCatalog']:
                if type(catalogInfo) == str:
                    catalogPath=catalogInfo
                    fluxCutJy=0.0
                elif type(catalogInfo) == dict:
                    catalogPath=catalogInfo['path']
                    fluxCutJy=catalogInfo['fluxCutJy']
                else:
                    raise Exception("Didn't understand contents of 'maskPointSourcesFromCatalog' - should be a path, or a dict with 'path' key.")
                tab=atpy.Table().read(catalogPath)
                if 'fluxJy' in tab.keys():
                    tab=tab[tab['fluxJy'] > fluxCutJy]
                tab=catalogs.getCatalogWithinImage(tab, data.shape, wcs)
                # If we're given a catalog that already has rArcmin in it, we use that to set hole size
                # Otherwise, if we have shape measurements (ellipse_A at least), we can use that
                for row in tab:
                    # Extended sources - identify by measured size > masking radius
                    # These will mess up noise term in filter, so add to psMask also and fill + smooth
                    # We won't fiddle with PA here, we'll just maximise based on x-pixel scale (because CAR)
                    if 'rArcmin' in tab.keys():
                        maskRadiusArcmin=row['rArcmin']
                    elif 'ellipse_A' in tab.keys():
                        xPixSizeArcmin=(wcs.getXPixelSizeDeg()/np.cos(np.radians(row['decDeg'])))*60
                        ASizeArcmin=row['ellipse_A']/xPixSizeArcmin
                        maskRadiusArcmin=ASizeArcmin/2
                    else:
                        raise Exception("To mask sources in a catalog, need either 'rArcmin' or 'ellipse_A' column to be present.")
                    rDegMap, xBounds, yBounds=makeDegreesDistanceMap(rDegMap, wcs,
                                                                     row['RADeg'], row['decDeg'],
                                                                     maskRadiusArcmin/60)
                    surveyMask[rDegMap < maskRadiusArcmin/60.0]=0
                    psMask[rDegMap < maskRadiusArcmin/60.0]=0
                    data[rDegMap < maskRadiusArcmin/60.0]=bckData[rDegMap < maskRadiusArcmin/60.0]

        if 'subtractModelFromCatalog' in list(self.keys()) and self['subtractModelFromCatalog'] is not None:
            if type(self['subtractModelFromCatalog']) is not list:
                self['subtractModelFromCatalog']=[self['subtractModelFromCatalog']]
            for tab in self['subtractModelFromCatalog']:
                if type(tab) != atpy.Table:
                    tab=atpy.Table().read(catalogPath)
                tab=catalogs.getCatalogWithinImage(tab, data.shape, wcs)
                model=makeModelImage(data.shape, wcs, tab, self['beamFileName'], obsFreqGHz = self['obsFreqGHz'])
                if model is not None:
                    data=data-model
                    # Threshold of > 1 uK here should be made adjustable in config
                    flagMask=flagMask+np.greater(model, 1)

        if 'maskAndFillFromCatalog' in list(self.keys()) and self['maskAndFillFromCatalog'] is not None:
            if type(self['maskAndFillFromCatalog']) is not list:
                self['maskAndFillFromCatalog']=[self['maskAndFillFromCatalog']]
            for tab in self['maskAndFillFromCatalog']:
                if type(tab) != atpy.Table:
                    tab=atpy.Table().read(catalogPath)
                tab=catalogs.getCatalogWithinImage(tab, data.shape, wcs)
                if len(tab) > 0 and 'ellipse_A' not in tab.keys():
                    raise Exception("Need to set measureShapes: True to use maskAndFillFromCatalog")
                for row in tab:
                    x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
                    rArcminMap=np.ones(data.shape, dtype = float)*1e6
                    if 'ellipse_A' and 'ellipse_B' in tab.keys():
                        xPixSizeArcmin=(wcs.getXPixelSizeDeg()/np.cos(np.radians(row['decDeg'])))*60
                        maskRadiusArcmin=(row['ellipse_A']/xPixSizeArcmin)/2
                    if 'maskHoleDilationFactor' in self.keys() and self['maskHoleDilationFactor'] is not None:
                        maskRadiusArcmin=maskRadiusArcmin*self['maskHoleDilationFactor']
                    rArcminMap, xBounds, yBounds=makeDegreesDistanceMap(rArcminMap, wcs,
                                                                        row['RADeg'], row['decDeg'],
                                                                        maskRadiusArcmin/60)
                    rArcminMap=rArcminMap*60
                    surveyMask[rArcminMap < maskRadiusArcmin]=0
                    psMask[rArcminMap < maskRadiusArcmin]=0
                    data[rArcminMap < maskRadiusArcmin]=bckData[rArcminMap < maskRadiusArcmin]

        # Add the map data to the dict
        self['data']=data
        self['weights']=weights
        self['wcs']=wcs
        self['surveyMask']=surveyMask
        self['pointSourceMask']=psMask
        self['flagMask']=flagMask
        self['tileName']=tileName

        # No point continuing if masks are different shape to map (easier to tell user here)
        if self['data'].shape != self['pointSourceMask'].shape:
            raise Exception("Map and point source mask dimensions are not the same (they should also have same WCS)")
        if self['data'].shape != self['surveyMask'].shape:
            raise Exception("Map and survey mask dimensions are not the same (they should also have same WCS)")


#------------------------------------------------------------------------------------------------------------
class MapDictList(object):
    """Blah. We want this to iterate over the mapDictList and be indexable.

    """

    def __init__(self, mapDictList, tileCoordsDict = None):
        """Blah.

        """

        self.mapDictList=[]
        for mapDict in mapDictList:
            self.mapDictList.append(MapDict(mapDict, tileCoordsDict))


    def __iter__(self):
        yield from self.mapDictList


    def __getitem__(self, item):
        return self.mapDictList[item]


#------------------------------------------------------------------------------------------------------------
class TileDict(dict):
    """A dictionary for collecting tile images, for later saving as multi-extension FITS or outputting as a
    single monolithic FITS image. Keys within the dictionary map to tile names. Handles on-the-fly
    reprojection between TAN and CAR if specified in the Nemo config.

    Args:
        inputDict (:obj:`dict`): Input dictionary (keys map to tile names).
        tileCoordsDict (:obj:`dict`, optional): A dictionary that describes the tiling of a large map, as
            produced by :meth:`startUp.NemoConfig.getTileCoordsDict`.

    Attributes:
        tileCoordsDict (:obj:`dict`): A dictionary that describes the tiling of a large map, as
            produced by :meth:`startUp.NemoConfig.getTileCoordsDict`.

    """

    def __init__(self, inputDict, tileCoordsDict = None):
        super(TileDict, self).__init__(inputDict)
        self.tileCoordsDict=tileCoordsDict


    def copy(self):
        """Make a copy of this :class:`TileDict` object.

        Returns:
            A deep copy of the :class:`TileDict` object.

        """
        return TileDict(self, tileCoordsDict = self.tileCoordsDict)


    def saveMEF(self, outFileName, compressionType = None):
        """Save the tile images as a multi-extension FITS file.

        Args:
            outFileName (:obj:`str`): Path where the MEF file will be written.
            compressionType (:obj:`str`): If given, the data will be compressed using the given
                method (as understood by :mod:`astropy.io.fits`). Use `PLIO_1` for masks,
                and `RICE_1` for other image data that can stand lossy compression. If None,
                the image data is not compressed.

        Returns:
            None

        """
        newImg=pyfits.HDUList()
        for tileName in self.keys():
            if self.tileCoordsDict[tileName]['reprojectToTan'] == True:
                wcs=astWCS.WCS(self.tileCoordsDict[tileName]['header'], mode = 'pyfits')
                tanWCS=_makeTanWCS(wcs)
                header=tanWCS.header
            else:
                header=self.tileCoordsDict[tileName]['header']
            if compressionType is not None:
                if compressionType == 'PLIO_1':
                    dtype=np.uint8
                else:
                    dtype=np.float32
                hdu=pyfits.CompImageHDU(np.array(self[tileName], dtype = dtype),
                                        header, name = tileName,
                                        compression_type = compressionType)
            else:
                hdu=pyfits.ImageHDU(self[tileName], header, name = tileName)
            newImg.append(hdu)
        newImg.writeto(outFileName, overwrite = True)


    def saveStitchedFITS(self, outFileName, stitchedWCS, compressionType = None):
        """Stitch together the tiles into a monolithic image and save in a FITS file.

        Args:
            outFileName (:obj:`str`): Path where the stitched image FITS file will be written.
            stitchedWCS (:obj:`astWCS.WCS`): WCS object corresponding to the stitched map
                that will be produced.
            compressionType (:obj:`str`): If given, the data will be compressed using the given
                method (as understood by :mod:`astropy.io.fits`). Use `PLIO_1` for masks,
                and `RICE_1` for other image data that can stand lossy compression. If None,
                the image data is not compressed.

        Returns:
            None

        """

        wcs=stitchedWCS
        d=np.zeros([stitchedWCS.header['NAXIS2'], stitchedWCS.header['NAXIS1']], dtype = np.float32)
        for tileName in self.keys():
            if self.tileCoordsDict[tileName]['reprojectToTan'] == True:
                carWCS=astWCS.WCS(self.tileCoordsDict[tileName]['header'], mode = 'pyfits')
                tanWCS=_makeTanWCS(carWCS)
                shape=[self.tileCoordsDict[tileName]['header']['NAXIS2'],
                       self.tileCoordsDict[tileName]['header']['NAXIS1']]
                if compressionType == 'PLIO_1':
                    order=0
                else:
                    order='bicubic'
                carData, footprint=reproject.reproject_interp((self[tileName], tanWCS.AWCS), carWCS.AWCS, shape_out = shape, order = order,
                                                              return_footprint = True)
                carData[footprint == 0]=0 # get rid of nans which will be in borders anyway
            else:
                carData=self[tileName]
            minX, maxX, minY, maxY=self.tileCoordsDict[tileName]['clippedSection']
            d[minY:maxY, minX:maxX]=d[minY:maxY, minX:maxX]+carData
        saveFITS(outFileName, d, wcs, compressionType = compressionType)

#-------------------------------------------------------------------------------------------------------------
def _makeTanWCS(wcs, pixScale = 0.5/60.0):
    """Generate a TAN WCS.

    Returns:
        TAN WCS

    """

    RADeg, decDeg=wcs.getCentreWCSCoords()
    CRVAL1, CRVAL2=RADeg, decDeg
    xSizeDeg, ySizeDeg=wcs.getFullSizeSkyDeg()
    xSizePix, ySizePix=int(xSizeDeg/pixScale), int(ySizeDeg/pixScale)
    xRefPix=xSizePix/2.0
    yRefPix=ySizePix/2.0
    xOutPixScale=xSizeDeg/xSizePix
    yOutPixScale=ySizeDeg/ySizePix
    newHead=pyfits.Header()
    newHead['NAXIS']=2
    newHead['NAXIS1']=xSizePix
    newHead['NAXIS2']=ySizePix
    newHead['CTYPE1']='RA---TAN'
    newHead['CTYPE2']='DEC--TAN'
    newHead['CRVAL1']=CRVAL1
    newHead['CRVAL2']=CRVAL2
    newHead['CRPIX1']=xRefPix+1
    newHead['CRPIX2']=yRefPix+1
    newHead['CDELT1']=-xOutPixScale
    newHead['CDELT2']=xOutPixScale    # Makes more sense to use same pix scale
    newHead['CUNIT1']='DEG'
    newHead['CUNIT2']='DEG'
    tanWCS=astWCS.WCS(newHead, mode='pyfits')

    return tanWCS

    import reproject
    tanData, footprint=reproject.reproject_interp((data, wcs.AWCS), tanWCS.AWCS, shape_out = [ySizePix, xSizePix],
                                                    order = 'bicubic', return_footprint = True)
    tanData[footprint == 0]=0 # get rid of nans which will be in borders anyway
    # checkData=reproject.reproject_interp((tanData, tanWCS.AWCS), wcs.AWCS, shape_out = data.shape, order = 'bicubic',
                                            # return_footprint = False)
    wcs=tanWCS
    data=tanData

#-------------------------------------------------------------------------------------------------------------
def convertToY(mapData, obsFrequencyGHz = 148):
    """Converts an array (e.g., a map) in ΔTemperature (μK) with respect to the CMB to Compton y parameter
    values at the given frequency.
    
    Args:
        mapData (:obj:`np.ndarray`): An array containing delta T (micro Kelvin, with respect to CMB) values.
        obsFrequencyGHz (:obj:`float`): Frequency in GHz at which to do the conversion.
    
    Returns:
        An array of Compton y parameter values.
    
    """
    fx=signals.fSZ(obsFrequencyGHz)    
    mapData=(mapData/(signals.TCMB*1e6))/fx # remember, map is in deltaT uK 
    
    return mapData

#-------------------------------------------------------------------------------------------------------------
def convertToDeltaT(mapData, obsFrequencyGHz = 148, TCMBAlpha = 0.0, z = None):
    """Converts an array (e.g., a map) of Compton y parameter values to ΔTemperature (μK) with respect to the
    CMB at the given frequency.
    
    Args:
        mapData (:obj:`np.ndarray`): An array containing Compton y parameter values.
        obsFrequencyGHz (:obj:`float`): Frequency in GHz at which to do the conversion.
        TCMBAlpha (:obj:`float`, optional): This should always be zero unless you really do want to make a model
            where CMB temperature evolves as T\ :sub:`0` * (1+z)\ :sup:`1-TCMBAlpha`.
        z (:obj:`float`, optional): Redshift - needed only if TCMBAlpha is non-zero.
    
    Returns:
        An array of ΔT (μK) values.
    
    """
    fx=signals.fSZ(obsFrequencyGHz, TCMBAlpha = TCMBAlpha, z = z)
    mapData=mapData*fx*(signals.TCMB*1e6)   # into uK
    
    return mapData

#-------------------------------------------------------------------------------------------------------------
def autotiler(surveyMask, wcs, targetTileWidth, targetTileHeight):
    """Given a survey mask (where values > 0 indicate valid area, and 0 indicates area to be ignored), 
    figure out an optimal tiling strategy to accommodate tiles of the given dimensions. The survey mask need
    not be contiguous (e.g., AdvACT and SO maps, using the default pixelization, can be segmented into three
    or more different regions).
    
    Args:
        surveyMask (numpy.ndarray): Survey mask image (2d array). Values > 0 will be taken to define valid 
            area.
        wcs (astWCS.WCS): WCS associated with survey mask image.
        targetTileWidth (float): Desired tile width, in degrees (RA direction for CAR).
        targetTileHeight (float): Desired tile height, in degrees (dec direction for CAR).
    
    Returns:
        Dictionary list defining tiles in same format as config file.
    
    Note:
        While this routine will try to match the target file sizes, it may not match exactly. Also,
        :meth:`startUp.NemoConfig.getTileCoordsDict` will expand tiles by a user-specified amount such that
        they overlap.
    
    """
    
    # This deals with identifying boss vs. full AdvACT footprint maps 
    mapCentreRA, mapCentreDec=wcs.getCentreWCSCoords()    
    skyWidth, skyHeight=wcs.getFullSizeSkyDeg()
    if mapCentreRA < 0.1 and skyWidth < 0.1 or skyWidth > 359.9:
        handle180Wrap=True
    else:
        handle180Wrap=False
    
    segMap=surveyMask
    try:
        numObjects=ndimage.label(segMap, output = segMap)
    except:
        raise Exception("surveyMask given for autotiler is probably too complicated (breaks into > 256 regions) - check your mask and/or config file.")

    # More memory efficient than previous version
    fieldIDs=np.arange(1, numObjects+1, dtype = segMap.dtype)
    maskSections=ndimage.find_objects(segMap)
    tileList=[]
    for maskSection, f in zip(maskSections, fieldIDs):
        yMin=maskSection[0].start
        yMax=maskSection[0].stop-1
        if yMax-yMin < 1000:  # In case of stray individual pixels (e.g., combined with extended sources mask)
            continue
        xc=int((maskSection[1].start+(maskSection[1].stop-1))/2)

        # Some people want to run on full sky CAR ... so we have to avoid that blowing up at the poles
        decMin, decMax=np.nan, np.nan
        deltaY=0
        while np.isnan(decMin) and np.isnan(decMax):
            RAc, decMin=wcs.pix2wcs(xc, yMin+deltaY)
            RAc, decMax=wcs.pix2wcs(xc, yMax-deltaY)
            deltaY=deltaY+0.01

        numRows=int((decMax-decMin)/targetTileHeight)
        if numRows == 0:
            raise Exception("targetTileHeight is larger than the height of the map - edit your config file accordingly.")
        tileHeight=np.ceil(((decMax-decMin)/numRows)*100)/100

        for i in range(numRows):
            decBottom=decMin+i*tileHeight
            decTop=decMin+(i+1)*tileHeight
            xc, yBottom=wcs.wcs2pix(RAc, decBottom)
            xc, yTop=wcs.wcs2pix(RAc, decTop)
            yBottom=int(yBottom)
            yTop=int(yTop)
            yc=int((yTop+yBottom)/2)

            strip=segMap[yBottom:yTop]
            ys, xs=np.where(strip == f)
            xMin=xs.min()
            xMax=xs.max()
            del ys, xs, strip
            stripWidthDeg=(xMax-xMin)*wcs.getXPixelSizeDeg()
            RAMax, decc=wcs.pix2wcs(xMin, yc)
            RAMin, decc=wcs.pix2wcs(xMax, yc)
            numCols=int(stripWidthDeg/targetTileWidth)
            tileWidth=np.ceil((stripWidthDeg/numCols)*100)/100
            #assert(tileWidth < targetTileWidth*1.1)

            stretchFactor=1/np.cos(np.radians(decTop))
            numCols=int(stripWidthDeg/(targetTileWidth*stretchFactor))
            for j in range(numCols):
                tileWidth=np.ceil((stripWidthDeg/numCols)*100)/100
                RALeft=RAMax-j*tileWidth
                RARight=RAMax-(j+1)*tileWidth
                if RALeft < 0:
                    RALeft=RALeft+360
                if RARight < 0:
                    RARight=RARight+360
                # HACK: Edge-of-map handling
                if handle180Wrap == True:
                    if RARight < 180.01 and RALeft < 180+tileWidth and RALeft > 180.01:
                        RARight=180.01
                # NOTE: floats here to make tileDefinitions.yml readable
                tileList.append({'tileName': '%d_%d_%d' % (f, i, j),
                                 'RADecSection': [float(RARight), float(RALeft), float(decBottom), float(decTop)]})

    return tileList

#------------------------------------------------------------------------------------------------------------
def saveTilesDS9RegionsFile(parDict, DS9RegionFileName):
    """Writes a DS9 .reg file containing the locations of tiles defined in parDict.
    
    Args:
        parDict (:obj:`dict`): Dictionary containing the contents of the Nemo config file.
        DS9RegionFileName (str): Path to DS9 regions file to be written.
    
    """

    if type(parDict['tileDefinitions']) is not list:
        raise Exception("parDict did not contain a list of tile definitions.")
    
    tileNames=[]
    coordsList=[]
    for tileDict in parDict['tileDefinitions']:
        ra0, ra1, dec0, dec1=tileDict['RADecSection']
        coordsList.append([ra0, ra1, dec0, dec1])
        tileNames.append(tileDict['tileName'])   
    with open(DS9RegionFileName, "w") as outFile:
        outFile.write("# Region file format: DS9 version 4.1\n")
        outFile.write('global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
        outFile.write("fk5\n")
        for c, name in zip(coordsList, tileNames):
            outFile.write('polygon(%.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f, %.6f) # text="%s"\n' % (c[0], c[2], c[0], c[3], c[1], c[3], c[1], c[2], name))  

#-------------------------------------------------------------------------------------------------------------
def shrinkWCS(origShape, origWCS, scaleFactor):
    """Given an astWCS object and corresponding image shape, scale the WCS by scaleFactor. Used for making 
    downsampled quicklook images (using stitchMaps).
    
    Args:
        origShape (tuple): Shape of the original image.
        origWCS (astWCS.WCS object): WCS for the original image.
        scaleFactor (float): The factor by which to scale the image WCS.
    Returns:
        shape (tuple), WCS (astWCS.WCS object)
    
    """
    
    scaledShape=[int(origShape[0]*scaleFactor), int(origShape[1]*scaleFactor)]
    scaledData=np.zeros(scaledShape)
    
    trueScaleFactor=np.array(scaledData.shape, dtype = float) / np.array(origShape, dtype = float)
    offset=0.
    imageWCS=origWCS.copy()
    try:
        oldCRPIX1=imageWCS.header['CRPIX1']
        oldCRPIX2=imageWCS.header['CRPIX2']
        CD11=imageWCS.header['CD1_1']
        CD21=imageWCS.header['CD2_1']
        CD12=imageWCS.header['CD1_2']
        CD22=imageWCS.header['CD2_2'] 
    except KeyError:
        oldCRPIX1=imageWCS.header['CRPIX1']
        oldCRPIX2=imageWCS.header['CRPIX2']
        CD11=imageWCS.header['CDELT1']
        CD21=0
        CD12=0
        CD22=imageWCS.header['CDELT2']

    CDMatrix=np.array([[CD11, CD12], [CD21, CD22]], dtype=np.float64)
    scaleFactorMatrix=np.array([[1.0/trueScaleFactor[1], 0], [0, 1.0/trueScaleFactor[0]]])
    scaleFactorMatrix=np.array([[1.0/trueScaleFactor[1], 0], [0, 1.0/trueScaleFactor[0]]])
    scaledCDMatrix=np.dot(scaleFactorMatrix, CDMatrix)

    scaledWCS=imageWCS.copy()
    scaledWCS.header['NAXIS1']=scaledData.shape[1]
    scaledWCS.header['NAXIS2']=scaledData.shape[0]
    scaledWCS.header['CRPIX1']=oldCRPIX1*trueScaleFactor[1]
    scaledWCS.header['CRPIX2']=oldCRPIX2*trueScaleFactor[0]
    scaledWCS.header['CD1_1']=scaledCDMatrix[0][0]
    scaledWCS.header['CD2_1']=scaledCDMatrix[1][0]
    scaledWCS.header['CD1_2']=scaledCDMatrix[0][1]
    scaledWCS.header['CD2_2']=scaledCDMatrix[1][1]
    scaledWCS.updateFromHeader()
    
    return scaledShape, scaledWCS

#-------------------------------------------------------------------------------------------------------------
def chunkLoadMask(fileName, numChunks = 8, dtype = np.uint8, returnWCS = True):
    """Load a FITS-format mask file (with default 8-bit integer values) in chunks, for memory efficiency,
    at the expense of some speed. Masks in compressed format (see :meth:`saveFITS`) are supported.

    Args:
        fileName (:obj:`str`): Path to the FITS-format mask file.
        numChunks (:obj:`int`): Number of chunks in which to load the file. Largers numbers use less memory,
            but it takes a little longer for the mask to load.
        returnWCS (:obj:`bool`, optional): If given, return the WCS of the mask.

    Returns:
        Mask image (2d array of 8-bit unsigned integers), and optionally a WCS object.

    Note:
        This can also be used to load large compressed maps in a memory-efficient way by setting
        ``dtype = np.float32``.

    """

    shape=None
    with pyfits.open(fileName) as img:
        for hdu in img:
            if hdu.data is not None:
                shape=hdu.data.shape
                if returnWCS == True:
                    wcs=astWCS.WCS(hdu.header, mode = 'pyfits')
                break
    del img

    height=shape[0]
    chunkSize=int(height/numChunks)
    maskArr=np.zeros(shape, dtype = dtype)
    for i in range(numChunks):
        with pyfits.open(fileName) as img:
            for hdu in img:
                if hdu.data is not None:
                    start=i*chunkSize
                    end=(i+1)*chunkSize
                    if end >= height:
                        end=height-1
                    chunk=hdu.data[start:end]
                    maskArr[start:end]=chunk
                    del chunk
            del hdu.data
        del img

    if returnWCS == True:
        return maskArr, wcs
    else:
        return maskArr

#-------------------------------------------------------------------------------------------------------------
def checkMask(fileName, numChunks = 8):
    """Checks whether a mask contains negative values (invalid) and throws an exception if this is the case.
    
    Args:
        fileName (str): Name of the FITS format mask file to check.
        
    """

    # This is now more horrid looking to save memory, at the expense of speed
    height=None
    with pyfits.open(fileName) as img:
        for hdu in img:
            if hdu.data is not None:
                height=hdu.data.shape[0]
    del img

    chunkSize=int(height/numChunks)
    for i in range(numChunks):
        with pyfits.open(fileName) as img:
            for hdu in img:
                if hdu.data is not None:
                    start=i*chunkSize
                    end=(i+1)*chunkSize
                    if end >= height:
                        end=height-1
                    chunk=hdu.data[start:end].flatten()
                    if np.any(chunk < 0) == True:
                        raise Exception("Mask file '%s' contains negative values - please fix your mask." % (fileName))
                    del chunk
            del hdu.data
        del img

#-------------------------------------------------------------------------------------------------------------
def maskOutSources(mapData, wcs, catalog, radiusArcmin = 7.0, mask = 0.0, growMaskedArea = 1.0):
    """Given a mapData array and a catalog of source positions, replace the values at the object positions 
    in the map within radiusArcmin with replacement values. If mask == 'whiteNoise', this will be white
    noise with mean and sigma set by the pixel values in an annulus of 1 < r < 2 * radiusArcmin.
    
    growMaskedArea sets factor larger than radiusArcmin to set masked area to in returned mask. This can
    avoid any weird artefacts making it into source lists.
    
    Returns a dictionary with keys 'data' (mapData with mask applied), 'mask' (0-1 mask of areas masked).
    
    """
        
    maskMap=np.zeros(mapData.shape)
    maskedMapData=np.zeros(mapData.shape, dtype=np.float64)+mapData    # otherwise, gets modified in place.
    
    bckSubbed=subtractBackground(mapData, wcs, smoothScaleDeg = 1.4/60.0) # for source subtracting
    
    mapInterpolator=interpolate.RectBivariateSpline(np.arange(mapData.shape[0]), 
                                                np.arange(mapData.shape[1]), 
                                                bckSubbed, kx = 1, ky = 1)

    for obj in catalog:
        if wcs.coordsAreInImage(obj['RADeg'], obj['decDeg']) == True:
            degreesMap=np.ones(mapData.shape, dtype = float)*1e6
            rRange, xBounds, yBounds=makeDegreesDistanceMap(degreesMap, wcs,
                                                            obj['RADeg'], obj['decDeg'],
                                                            20.0/60.0)
            circleMask=np.less(rRange, radiusArcmin/60.0)
            grownCircleMask=np.less(rRange, (radiusArcmin*growMaskedArea)/60.0)
            maskMap[grownCircleMask]=1.0
            if type(mask) == float or type(mask) == int:
                maskedMapData[circleMask]=mask

            elif mask == 'shuffle':
                # How about copying random pixels from the vicinity into the area to be masked?
                annulusMask=np.logical_and(np.greater(rRange, 5.0/60.0), \
                                            np.less(rRange, 10.0/60.0))
                annulusValues=mapData[annulusMask].flatten()
                indices=np.random.randint(0, annulusValues.shape[0], circleMask.flatten().nonzero()[0].shape[0])
                maskedMapData[circleMask]=annulusValues[indices]
                
            elif mask == 'subtract':         
                peakValue=mapData[int(round(obj['y'])), int(round(obj['x']))]
                sigmaDeg=(1.4/60.0)/np.sqrt(8.0*np.log(2.0))            
                profRDeg=np.linspace(0.0, 30.0/60.0, 5000)
                profile1d=peakValue*np.exp(-((profRDeg**2)/(2*sigmaDeg**2)))                
                r2p=interpolate.interp1d(profRDeg, profile1d, bounds_error=False, fill_value=0.0)
                profile2d=np.zeros(rRange.shape)
                profMask=np.less(rRange, 1.0)
                profile2d[profMask]=r2p(rRange[profMask])
                maskedMapData[profMask]=maskedMapData[profMask]-profile2d[profMask]
                
                # NOTE: below old, replaced Jul 2015 but not deleted as yet...
                # 1.3197 is a correction factor for effect of filtering on bckSubbed
                # Worked out by comparing peak value of bckSubbed profile2d only map
                #peakValue=mapInterpolator(obj['y'], obj['x'])[0][0]*1.3197   
                #sigmaDeg=(1.4/60.0)/np.sqrt(8.0*np.log(2.0))            
                #profRDeg=np.linspace(0.0, 30.0/60.0, 5000)
                #profile1d=peakValue*np.exp(-((profRDeg**2)/(2*sigmaDeg**2)))                
                #r2p=interpolate.interp1d(profRDeg, profile1d, bounds_error=False, fill_value=0.0)
                #profile2d=np.zeros(rRange.shape)
                #profMask=np.less(rRange, 1.0)
                #profile2d[profMask]=r2p(rRange[profMask])
                #maskedMapData[profMask]=maskedMapData[profMask]-profile2d[profMask]
            
                
            elif mask == "whiteNoise":
                # Get pedestal level and white noise level from average between radiusArcmin and  2*radiusArcmin
                annulusMask=np.logical_and(np.greater(rRange, 2*radiusArcmin/60.0), \
                                            np.less(rRange, 4*radiusArcmin/60.0))
                maskedMapData[circleMask]=np.random.normal(mapData[annulusMask].mean(), \
                                                            mapData[annulusMask].std(),  \
                                                            mapData[circleMask].shape)
    
    return {'data': maskedMapData, 'mask': maskMap}

#-------------------------------------------------------------------------------------------------------------
def applyPointSourceMask(maskFileName, mapData, mapWCS, mask = 0.0, radiusArcmin = 2.8):
    """Given file name pointing to a point source mask (as made by maskOutSources), apply it to given mapData.
    
    """
    
    img=pyfits.open(maskFileName)
    maskData=img[0].data

    maskedMapData=np.zeros(mapData.shape)+mapData    # otherwise, gets modified in place.
    
    # Thresholding to identify significant pixels
    threshold=0
    sigPix=np.array(np.greater(maskData, threshold), dtype=int)
    sigPixMask=np.equal(sigPix, 1)
    
    # Fast, simple segmentation - don't know about deblending, but doubt that's a problem for us
    segmentationMap, numObjects=ndimage.label(sigPix)
    
    # Get object positions, number of pixels etc.
    objIDs=np.unique(segmentationMap)
    objPositions=ndimage.center_of_mass(maskData, labels = segmentationMap, index = objIDs)
    objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
    
    for objID, pos, numPix in zip(objIDs, objPositions, objNumPix):
        circleMask=np.equal(segmentationMap, objID)
        if type(mask) == float or type(mask) == int:
            maskedMapData[circleMask]=mask
        elif mask == "subtract":
            print("Add code to subtract point sources")
            ipshell()
            sys.exit()
        elif mask == "whiteNoise":
            RADeg, decDeg=mapWCS.pix2wcs(pos[1], pos[0])
            if np.isnan(RADeg) == False and np.isnan(decDeg) == False:
                degreesMap=np.ones(mapData.shape, dtype = float)*1e6
                rRange, xBounds, yBounds=makeDegreesDistanceMap(degreesMap, mapWCS,
                                                                RADeg, decDeg,
                                                                (radiusArcmin*4)/60.0)
                # Get pedestal level and white noise level from average between radiusArcmin and  2*radiusArcmin
                annulusMask=np.logical_and(np.greater(rRange, radiusArcmin/60.0), \
                                              np.less(rRange, 2*radiusArcmin/60.0))
                # Below just does a quick sanity check - we don't bother masking if std == 0, because we're
                # most likely applying this in the middle of a fake source sim with map set to zero for testing
                sigma=mapData[annulusMask].std()
                if sigma > 0:
                    maskedMapData[circleMask]=np.random.normal(mapData[annulusMask].mean(), \
                                                                  sigma,  \
                                                                  mapData[circleMask].shape)
    
    return maskedMapData
                                                             
#-------------------------------------------------------------------------------------------------------------
def addWhiteNoise(mapData, noisePerPix):
    """Adds Gaussian distributed white noise to mapData.
    
    """
    
    noise=np.random.normal(0, noisePerPix, mapData.shape)
    mapData=mapData+noise
    
    return mapData

#------------------------------------------------------------------------------------------------------------
def simCMBMap(shape, wcs, noiseLevel = None, beam = None, seed = None):
    """Generate a simulated CMB map, optionally convolved with the beam and with (white) noise added.
    
    Args:
        shape (:obj:`tuple`): A tuple describing the map (numpy array) shape in pixels (height, width).
        wcs (:obj:`astWCS.WCS`): An astWCS object.
        noiseLevel (:obj:`numpy.ndarray` or float): If a single number, this is taken as sigma (in map units,
            usually uK) for generating white noise that is added across the whole map. Alternatively, an array
            with the same dimensions as shape may be used, specifying sigma (in map units) per corresponding 
            pixel. Noise will only be added where non-zero values appear in noiseLevel.
        beam (:obj:`str` or :obj:`signals.BeamProfile`): Either the file name of the text file that describes
            the beam with which the map will be convolved, or a :obj:`signals.BeamProfile` object. If None,
            no beam convolution is applied.
        seed (:obj:`int`): The seed used for the random CMB realisation.
            
    Returns:
        A map (:obj:`numpy.ndarray`)
    
    """

    # Power spectrum array ps here is indexed by ell, starting from 0
    # i.e., each element corresponds to the power at ell = 0, 1, 2 ... etc.
    ps=powspec.read_spectrum(nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"planck_lensedCls.dat",
                             scale = True, expand = None)
    ps=ps[0]
    lps=np.arange(0, len(ps))

    if beam is not None:
        if type(beam) == str:
            beam=signals.BeamProfile(beamFileName = beam)
        assert(type(beam) == signals.BeamProfile)
        lbeam=np.interp(lps, beam.ell, beam.Bell)
        ps*=np.power(lbeam, 2)

    randMap=curvedsky.rand_map(shape, wcs.AWCS, ps = ps, spin = [0,2], seed = seed)

    if noiseLevel is not None:
        randMap=randMap+simNoiseMap(shape, noiseLevel)

    np.random.seed()
    
    return randMap

#-------------------------------------------------------------------------------------------------------------
def simNoiseMap(shape, noiseLevel, wcs = None, lKnee = None, alpha = -3, noiseMode = 'perPixel'):
    """Generate a simulated noise map. This may contain just white noise, or optionally a 1/f noise component
    can be generated.

    Args:
        shape (:obj:`tuple`): A tuple describing the map (numpy array) shape in pixels (height, width).
        noiseLevel (:obj:`numpy.ndarray` or float): If a single number, this is taken as sigma (in map units,
            usually uK) for generating white noise that is added across the whole map. Alternatively, an array
            with the same dimensions as shape may be used, specifying sigma (in map units) per corresponding
            pixel. Noise will only be added where non-zero values appear in noiseLevel.
        wcs (:obj:`astWCS.WCS`, optional): WCS corresponding to the map shape.
        lKnee (:obj:`float`, optional): If given, 1/f noise will be generated using the power spectrum
            N_l = (1 + l/lknee)^-alpha) - see Appendix A of MacCrann et al. 2023.
        alpha (:obj:`float`, optional): Power-law exponent in the power spectrum used for generating 1/f noise. Has
            no effect unless lKnee is also given.
        noiseMode(:obj:`str`, optional): Either 'perPixel', or 'perSquareArcmin' - if the latter, constant noise in terms
            of surface brightness will be added (accounts for varying pixel scale, if present - which requires
            `wcs` to be supplied).

    Returns:
        A map (:obj:`numpy.ndarray`)

    """

    np.random.seed()

    assert(noiseMode in ['perPixel', 'perSquareArcmin'])
    if noiseMode == 'perSquareArcmin' and lKnee is not None:
        raise Exception("Adding 1/f noise when noiseMode != 'perPixel' is not supported yet")
    if noiseMode == 'perSquareArcmin' and type(noiseLevel) == np.ndarray:
        raise Exception("noiseLevel is a map - this is only currently supported if noiseMode = 'perPixel' (noiseMode = 'perSquareArcmin' given)")

    if lKnee is None:
        # White noise only
        randMap=np.zeros(shape)
        generatedNoise=np.zeros(randMap.shape)
        if type(noiseLevel) == np.ndarray:
            mask=np.nonzero(noiseLevel)
            generatedNoise[mask]=np.random.normal(0, noiseLevel[mask], noiseLevel[mask].shape)
        else:
            if noiseLevel > 0:
                if noiseMode == 'perPixel':
                    generatedNoise=np.random.normal(0, noiseLevel, randMap.shape)
                else:
                    arcmin2Map=getPixelAreaArcmin2Map(shape, wcs)
                    generatedNoise=np.random.normal(0, noiseLevel/arcmin2Map, randMap.shape)
        randMap=randMap+generatedNoise

    else:
        # 1/f noise + white noise, using Niall's routines
        mlmax=6000 # following config in Niall's code, could be made a parameter
        if wcs is None:
            raise Exception("wcs is None - need to supply a wcs to generate a noise map with 1/f noise included.")
        if type(noiseLevel) == np.ndarray:
            mask=noiseLevel > 1e-07
            ivarMap=np.zeros(shape)
            ivarMap[mask]=1/noiseLevel[mask]**2
            ivarMap=enmap.enmap(ivarMap, wcs.AWCS)
        else:
            ivarMap=enmap.enmap(np.ones(shape)*(1/noiseLevel**2), wcs.AWCS)

        def _mod_noise_map(ivar, Nl):
            map1 = enmap.rand_gauss(ivar.shape, ivar.wcs)
            lmax = len(Nl)-1
            ainfo = curvedsky.alm_info(lmax)
            alm = curvedsky.map2alm(map1, ainfo=ainfo)
            map2 = curvedsky.alm2map(alm, np.zeros_like(map1))
            map1 -= map2
            ainfo.lmul(alm, Nl**0.5, alm)
            map2 = curvedsky.alm2map(alm, np.zeros_like(map1))
            map1 += map2
            map1 *= ivar**-0.5
            ivar_nonzero = ivar>0.
            ivar_median = np.median(ivar[ivar_nonzero])
            valid_ivar = ivar_nonzero*(ivar>ivar_median/1.e6)
            map1[~(valid_ivar)] = 0.
            return map1

        assert np.all(np.isfinite(ivarMap))
        shape, wcs=ivarMap.shape, ivarMap.wcs
        ells = np.arange(mlmax+1)
        Nl_atm = (lKnee/ells)**-alpha + 1
        Nl_atm[~np.isfinite(Nl_atm)] = 0.
        assert np.all(np.isfinite(Nl_atm))
        randMap = _mod_noise_map(ivarMap, Nl_atm)
        assert np.all(np.isfinite(randMap))

    return randMap

#-------------------------------------------------------------------------------------------------------------
def subtractBackground(data, wcs, RADeg = 'centre', decDeg = 'centre', smoothScaleDeg = 30.0/60.0):
    """Smoothes map with Gaussian of given scale and subtracts it, to get rid of large scale power.
    
    If RADeg, decDeg = 'centre', then the pixel scales used to set the kernel shape will be set from that at the
    centre of the WCS. Otherwise, they will be taken at the given coords.
    
    Note that wcs is only used to figure out the pixel scales here.
    
    """
            
    smoothedData=smoothMap(data, wcs, RADeg = RADeg, decDeg = decDeg, smoothScaleDeg = smoothScaleDeg)
    data=data-smoothedData
    
    return data

#------------------------------------------------------------------------------------------------------------
def convolveMapWithBeam(data, wcs, beam, maxDistDegrees = 1.0):
    """Convolves map defined by data, wcs with the beam.
    
    Args:
        data (:obj:`numpy.ndarray`): Map to convolve, as 2d array.
        wcs (:obj:`astWCS.WCS`): WCS corresponding to data (i.e., the map).
        beam (:obj:`BeamProfile` or str): Either a BeamProfile object, or a string that gives the path to a 
            text file that describes the beam profile.
        maxDistDegrees (float): Sets the size of the convolution kernel, for optimization purposes.
    
    Returns:
        Beam-convolved map (numpy array).
    
    Note:
        The pixel scale used to define the convolution kernel is evaluated at the central map pixel. So, 
        this routine should only be used with either pixelisations where the scale is constant or on 
        relatively small tiles.
        
    """
    
    if type(beam) == str:
        beam=signals.BeamProfile(beamFileName = beam)

    # Pad the beam kernel to odd number of pixels (so we know shift to apply)
    # We're only really using WCS info here for the pixel scale at the centre of the map
    if data.shape[0] % 2 == 0:
        yPad=1
    else:
        yPad=0
    if data.shape[1] % 2 == 0:
        xPad=1
    else:
        xPad=0
    degreesMap=np.ones([data.shape[0]+yPad, data.shape[1]+xPad], dtype = float)*1e6
    RADeg, decDeg=wcs.pix2wcs(int(degreesMap.shape[1]/2)+1, int(degreesMap.shape[0]/2)+1)
    degreesMap, xBounds, yBounds=makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg,
                                                        maxDistDegrees)
    beamMap=signals.makeBeamModelSignalMap(degreesMap, wcs, beam)
    if (yBounds[1]-yBounds[0]) > beamMap.shape[1] and (yBounds[1]-yBounds[0]) % 2 == 0:
        yBounds[0]=yBounds[0]-1
    if (xBounds[1]-xBounds[0]) > beamMap.shape[0] and (xBounds[1]-xBounds[0]) % 2 == 0:
        xBounds[0]=xBounds[0]-1    
    beamMap=beamMap[yBounds[0]:yBounds[1], xBounds[0]:xBounds[1]]
    beamMap=beamMap/np.sum(beamMap)

    # For testing for shift
    # This shows we get (-1, -1) shift with scipy_convolve and odd-shaped kernel
    #testMap=np.zeros([301, 301])
    #yc1=151
    #xc1=151
    #testMap[yc1, xc1]=1.
    #outMap=scipy_convolve(testMap, beamMap, mode = 'same')
    #yc2, xc2=np.where(outMap == outMap.max())
    #yc2=int(yc2)
    #xc2=int(xc2)
    #outMap=ndimage.shift(outMap, [yc1-yc2, xc1-xc2])
    
    outMap=ndimage.shift(scipy_convolve(data, beamMap, mode = 'same'), [-1, -1])
        
    return outMap

#-------------------------------------------------------------------------------------------------------------
def smoothMap(data, wcs, RADeg = 'centre', decDeg = 'centre', smoothScaleDeg = 5.0/60.0):
    """Smoothes map with Gaussian of given scale.
    
    If RADeg, decDeg = 'centre', then the pixel scales used to set the kernel shape will be set from that at the
    centre of the WCS. Otherwise, they will be taken at the given coords.
    
    Note that wcs is only used to figure out the pixel scales here.
    
    """
    
    ra0, dec0=wcs.getCentreWCSCoords()
    if RADeg != 'centre':
        ra0=float(RADeg)
    if decDeg != 'centre':
        dec0=float(decDeg)
    x0, y0=wcs.wcs2pix(ra0, dec0)
    x1=x0+1
    y1=y0+1
    ra1, dec1=wcs.pix2wcs(x1, y1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    xSmoothScalePix=smoothScaleDeg/xPixScale
    ySmoothScalePix=smoothScaleDeg/yPixScale
    smoothedData=ndimage.gaussian_filter(data, (ySmoothScalePix, xSmoothScalePix))
    
    return smoothedData
    
#-------------------------------------------------------------------------------------------------------------
def getPixelAreaArcmin2Map(shape, wcs):
    """Returns a map of pixel area in arcmin2
    
    """
    
    # Get pixel size as function of position
    pixAreasDeg2=[]
    RACentre, decCentre=wcs.getCentreWCSCoords()
    x0, y0=wcs.wcs2pix(RACentre, decCentre)
    x1=x0+1
    for y0 in range(shape[0]):
        y1=y0+1
        ra0, dec0=wcs.pix2wcs(x0, y0)
        ra1, dec1=wcs.pix2wcs(x1, y1)
        xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
        yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
        pixAreasDeg2.append(xPixScale*yPixScale)
    pixAreasDeg2=np.array(pixAreasDeg2)
    pixAreasArcmin2=pixAreasDeg2*(60**2)
    pixAreasArcmin2Map=np.array([pixAreasArcmin2]*shape[1]).transpose()
    
    return pixAreasArcmin2Map    
    
#-------------------------------------------------------------------------------------------------------------
def estimateContaminationFromSkySim(config, imageDict):
    """Estimate contamination by running on source-free sky simulations (CMB plus noise that we generate here
    on the fly).
    
    This uses the same kernels that were constructed and used on the real maps. The whole filtering and object
    detection pipeline is run on the simulated maps repeatedly. The number of sky sims used (set by numSkySims
    in the .yml config file) should be fairly large (~100) for the results to be robust (results on individual
    sims can vary by a lot).
    
    Args:
        config (:obj:`startUp.NemoConfig`): Nemo configuration object.
        imageDict (:obj:`dict`): A dictionary containing the output filtered maps and catalogs from running on 
            the real data (i.e., the output of pipelines.filterMapsAndMakeCatalogs). This will not be modified,
            but is used for estimating the contamination rate by comparison to the source-free sims.
    
    Returns:
        A dictionary where each key points to an astropy Table object containing the average contamination 
        estimate corresponding to SNR (maximal estimate) and fixed_SNR (for the chosen reference filter 
        scale).
    
    """

    simRootOutDir=config.diagnosticsDir+os.path.sep+"skySim_rank%d" % (config.rank)
    SNRKeys=['fixed_SNR']        
    numSkySims=config.parDict['numSkySims']
    resultsList=[]
    for i in range(numSkySims):
        
        # NOTE: we throw the first sim away on figuring out noiseBoostFactors
        print(">>> Sky sim %d/%d [rank = %d] ..." % (i+1, numSkySims, config.rank))
        t0=time.time()

        # We don't copy this, because it's complicated due to containing MPI-related things (comm)
        # So... we modify the config parameters in-place, and restore them before exiting this method
        simConfig=config
        
        # We use the seed here to keep the CMB sky the same across frequencies...
        CMBSimSeed=np.random.randint(16777216)
        
        # NOTE: This block below should be handled when parsing the config file - fix/remove
        # Optional override of default GNFW parameters (used by Arnaud model), if used in filters given
        if 'GNFWParams' not in list(simConfig.parDict.keys()):
            simConfig.parDict['GNFWParams']='default'
        for filtDict in simConfig.parDict['mapFilters']:
            filtDict['params']['GNFWParams']=simConfig.parDict['GNFWParams']
        
        # Delete all non-reference scale filters (otherwise we'd want to cache all filters for speed)
        for filtDict in simConfig.parDict['mapFilters']:
            if filtDict['label'] == simConfig.parDict['photFilter']:
                break
        simConfig.parDict['mapFilters']=[filtDict] 
        
        # Filling in with sim will be done when maps.preprocessMapDict is called by the filter object
        for mapDict in simConfig.unfilteredMapsDictList:
            mapDict['CMBSimSeed']=CMBSimSeed
                    
        # NOTE: we need to zap ONLY specific maps for when we are running in parallel
        for tileName in simConfig.tileNames:
            mapFileNames=glob.glob(simRootOutDir+os.path.sep+"filteredMaps"+os.path.sep+"*#%s_*.fits" % (tileName))
            for m in mapFileNames:
                os.remove(m)
                
        simImageDict=pipelines.filterMapsAndMakeCatalogs(simConfig, 
                                                         rootOutDir = simRootOutDir,
                                                         useCachedFilters = True)
        
        # Write out the last sim map catalog for debugging
        # NOTE: tileName here makes no sense - this should be happening in the pipeline call above
        #optimalCatalogFileName=simRootOutDir+os.path.sep+"CMBSim_optimalCatalog#%s.csv" % (tileName)    
        #optimalCatalog=simImageDict['optimalCatalog']
        #if len(optimalCatalog) > 0:
            #catalogs.writeCatalog(optimalCatalog, optimalCatalogFileName.replace(".csv", ".fits"), constraintsList = ["SNR > 0.0"])
        
        # Contamination estimate...
        contaminTabDict=estimateContamination(simImageDict, imageDict, SNRKeys, 'skySim', config.diagnosticsDir)
        resultsList.append(contaminTabDict)
        t1=time.time()
        print("... time taken for sky sim run = %.3f sec" % (t1-t0))

    # Average results
    avContaminTabDict={}
    for k in list(resultsList[0].keys()):
        avContaminTabDict[k]=atpy.Table()
        for kk in list(resultsList[0][k].keys()):
            avContaminTabDict[k].add_column(atpy.Column(np.zeros(len(resultsList[0][k])), kk))
            for i in range(len(resultsList)):
                avContaminTabDict[k][kk]=avContaminTabDict[k][kk]+resultsList[i][k][kk]
            avContaminTabDict[k][kk]=avContaminTabDict[k][kk]/float(len(resultsList))
    
    # For writing separate contamination .fits tables if running in parallel
    # (if we're running in serial, then we'll get a giant file name with full tileNames list... fix later)
    tileNamesLabel="#"+str(config.tileNames).replace("[", "").replace("]", "").replace("'", "").replace(", ", "#")
    for k in list(avContaminTabDict.keys()):
        fitsOutFileName=config.diagnosticsDir+os.path.sep+"%s_contaminationEstimate_%s.fits" % (k, tileNamesLabel)
        contaminTab=avContaminTabDict[k]
        contaminTab.meta['NEMOVER']=nemo.__version__
        contaminTab.write(fitsOutFileName, overwrite = True)
    
    # Restore the original config parameters (which we overrode to make the sims here)
    config.restoreConfig()
    
    return avContaminTabDict

#-------------------------------------------------------------------------------------------------------------
def estimateContaminationFromInvertedMaps(config, imageDict):
    """Run the whole filtering set up again, on inverted maps.
    
    Writes a DS9. reg file, which contains only the highest SNR contaminants (since these
    are most likely to be associated with artefacts in the map - e.g., point source masking).
    
    Writes a plot and a .fits table to the diagnostics dir.
    
    Runs over both SNR and fixed_SNR values.
    
    Returns a dictionary containing the results
    
    """
    
    invertedDict={}
    ignoreKeys=['optimalCatalog', 'mergedCatalog']
    for key in imageDict:
        if key not in ignoreKeys:
            invertedDict[key]=imageDict[key]
            
    invertedDict=pipelines.filterMapsAndMakeCatalogs(config, measureFluxes = False, invertMap = True)
    
    SNRKeys=['SNR', 'fixed_SNR']
    contaminTabDict=estimateContamination(invertedDict, imageDict, SNRKeys, 'invertedMap', config.diagnosticsDir)

    for k in list(contaminTabDict.keys()):
        fitsOutFileName=config.diagnosticsDir+os.path.sep+"%s_contaminationEstimate.fits" % (k)
        contaminTab=contaminTabDict[k]
        contaminTab.write(fitsOutFileName, overwrite = True)
        
    return contaminTabDict

#------------------------------------------------------------------------------------------------------------
def plotContamination(contaminTabDict, diagnosticsDir):
    """Makes contamination rate plots, output stored under diagnosticsDir
    
    While we're at it, we write out a text file containing interpolated values for e.g., 5%, 10% 
    contamination levels
    
    """

    plotSettings.update_rcParams()

    for k in list(contaminTabDict.keys()):
        if k.find('fixed') != -1:
            SNRKey="fixed_SNR"
            SNRLabel="SNR$_{\\rm 2.4}$"
        else:
            SNRKey="SNR"
            SNRLabel="SNR"
        binEdges=contaminTabDict[k][SNRKey]
        cumContamination=contaminTabDict[k]['cumContamination']
        plt.figure(figsize=(9,6.5))
        ax=plt.axes([0.10, 0.11, 0.87, 0.87])  
        plt.plot(binEdges, cumContamination, 'k-')# % (l))#, label = legl)
        plt.xlabel("%s" % (SNRLabel))#, fontdict = fontDict)
        plt.ylabel("Contamination fraction > %s" % (SNRLabel))#, fontdict = fontDict)
        allLabels=['4.0', '', '', '', '', '5.0', '', '', '', '', '6.0', '', '', '', '', '7.0', '', '', '', '', '8.0']
        allTicks=np.arange(4.0, 8.2, 0.2)
        plt.xticks(allTicks, allLabels)
        plt.xlim(4, 8)
        #plt.xlim(binMin, 10.01)#binMax)
        plt.ylim(-0.05, 0.6)
        #plt.legend()
        plt.savefig(diagnosticsDir+os.path.sep+"%s_contaminationEstimate.pdf" % (k))
        plt.close()  
        
        tck=interpolate.splrep(binEdges, contaminTabDict[k]['cumContamination'])
        fineSNRs=np.linspace(binEdges.min(), binEdges.max(), 1000)
        fineContamination=interpolate.splev(fineSNRs, tck, ext = 1)
        with open(diagnosticsDir+os.path.sep+"%s_contaminationEstimate_usefulFractions.txt" % (k), "w") as outFile:
            fracs=[0.4, 0.3, 0.2, 0.1, 0.05, 0.01]
            for f in fracs:
                SNRf=fineSNRs[np.argmin(abs(fineContamination-f))]
                logStr="... contamination fraction = %.2f for %s > %.3f ..." % (f, SNRKey, SNRf)
                print(logStr)
                outFile.write(logStr+"\n")
        
#------------------------------------------------------------------------------------------------------------
def estimateContamination(contamSimDict, imageDict, SNRKeys, label, diagnosticsDir):
    """Performs the actual contamination estimate, makes output under diagnosticsDir.
        
    Use label to set a prefix for output (plots / .fits tables), e.g., label = "skySim"
    
    """
    
    invertedDict=contamSimDict
    contaminTabDict={}
    for SNRKey in SNRKeys:
        #catalogs.catalog2DS9(invertedDict['optimalCatalog'], rootOutDir+os.path.sep+"skySimCatalog_%s_gtr_5.reg" % (SNRKey), 
                                 #constraintsList = ['%s > 5' % (SNRKey)])
        
        invertedSNRs=[]
        for obj in invertedDict['optimalCatalog']:
            invertedSNRs.append(obj[SNRKey])
        invertedSNRs=np.array(invertedSNRs)
        invertedSNRs.sort()
        numInverted=np.arange(len(invertedSNRs))+1
        
        candidateSNRs=[]
        for obj in imageDict['optimalCatalog']:
            candidateSNRs.append(obj[SNRKey])
        candidateSNRs=np.array(candidateSNRs)
        candidateSNRs.sort()
        numCandidates=np.arange(len(candidateSNRs))+1
        
        binMin=4.0
        binMax=20.0
        binStep=0.2
        binEdges=np.linspace(binMin, binMax, int((binMax-binMin)/binStep+1))
        binCentres=(binEdges+binStep/2.0)[:-1]
        candidateSNRHist=np.histogram(candidateSNRs, bins = binEdges)
        invertedSNRHist=np.histogram(invertedSNRs, bins = binEdges)    
        
        cumSumCandidates=[]
        cumSumInverted=[]
        for i in range(binCentres.shape[0]):
            cumSumCandidates.append(candidateSNRHist[0][i:].sum())
            cumSumInverted.append(invertedSNRHist[0][i:].sum())
        cumSumCandidates=np.array(cumSumCandidates, dtype = float)
        cumSumInverted=np.array(cumSumInverted, dtype = float)
        
        # Plot cumulative contamination estimate (this makes more sense than plotting purity, since we don't know
        # that from what we're doing here, strictly speaking)
        cumContamination=np.zeros(cumSumCandidates.shape)
        mask=np.greater(cumSumCandidates, 0)
        cumContamination[mask]=cumSumInverted[mask]/cumSumCandidates[mask]
        
        # Remember, this is all cumulative (> SNR, so lower bin edges)
        contaminDict={}
        contaminDict['%s' % (SNRKey)]=binEdges[:-1]
        contaminDict['cumSumRealCandidates']=cumSumCandidates
        contaminDict['cumSumSimCandidates']=cumSumInverted
        contaminDict['cumContamination']=cumContamination       
        
        # Convert to .fits table
        contaminTab=atpy.Table()
        for key in list(contaminDict.keys()):
            contaminTab.add_column(atpy.Column(contaminDict[key], key))
        
        contaminTabDict['%s_%s' % (label, SNRKey)]=contaminTab
        
    return contaminTabDict

#------------------------------------------------------------------------------------------------------------
def makeModelImage(shape, wcs, catalog, beamFileName, obsFreqGHz = None, GNFWParams = 'default',\
                   profile = 'A10', cosmoModel = None, applyPixelWindow = True, override = None,\
                   validAreaSection = None, minSNR = -99, TCMBAlpha = 0, reportTimingInfo = False):
    """Make a map with the given dimensions (shape) and WCS, containing model clusters or point sources, 
    with properties as listed in the catalog. This can be used to either inject or subtract sources
    from real maps.
    
    Args:
        shape (tuple): The dimensions of the output map (height, width) that will contain the model sources.
        wcs (:obj:`astWCS.WCS`): A WCS object that defines the coordinate system of the map. 
        catalog (:obj:`astropy.table.Table` or str): An astropy Table object containing the catalog,
            or a string containing the path to a catalog that astropy Table understands. The catalog must
            include columns named 'RADeg', 'decDeg' that give object coordinates. For point sources, the 
            amplitude in uK must be given in a column named 'deltaT_c'. For clusters, either 'M500' (in 
            units of 10^14 MSun), 'z', and 'fixed_y_c' must be given (as in a mock catalog), OR the 
            catalog must contain a 'template' column, with templates named like, e.g., Arnaud_M1e14_z0p2
            (for a z = 0.2, M500 = 1e14 MSun cluster; see the example .yml config files included with nemo).
        beamFileName: Path to a text file that describes the beam.
        obsFreqGHz (float, optional): Used only by cluster catalogs - if given, the returned map will be 
            converted into delta T uK, assuming the given frequency. Otherwise, a y0 map is returned.
        GNFWParams (str or dict, optional): Used only by cluster catalogs. If 'default', the Arnaud et al. 
            (2010) Universal Pressure Profile is assumed. Otherwise, a dictionary that specifies the profile
            parameters can be given here (see gnfw.py).
        profile (str, optional): Used by cluster models only - sets the profile shape to use: 'A10'
            for Arnaud et al. (2010) UPP models, or 'B12' for Battaglia et al. (2012) models.
        override (dict, optional): Used only by cluster catalogs. If a dictionary containing keys
            {'M500', 'redshift'} is given, all objects in the model image are forced to have the 
            corresponding angular size. Used by :meth:`sourceInjectionTest`.
        applyPixelWindow (bool, optional): If True, apply the pixel window function to the map.
        validAreaSection (list, optional): Pixel coordinates within the wcs in the format
            [xMin, xMax, yMin, yMax] that define valid area within the model map. Pixels outside this 
            region will be set to zero. Use this to remove overlaps between tile boundaries.
        minSNR (float, optional): Only include objects with SNR (or fixed_SNR) > this value in the model.
            If found, the 'SNR' column will be used, otherwise the 'fixed_SNR' column will be used. If
            neither is present, no cuts on the catalog will be performed.
        TCMBAlpha (float, optional): This should always be zero unless you really do want to make a
            cluster model image where CMB temperature evolves as T0*(1+z)^{1-TCMBAlpha}.
        reportTimingInfo (bool, optional): If True, report how long each step takes.
        
    Returns:
        Map containing injected sources, or None if there are no objects within the map dimensions.
    
    """

    modelMap=enmap.zeros(shape, dtype = np.float32) #np.zeros(shape, dtype = np.float32)

    if type(catalog) == str:
        catalog=atpy.Table().read(catalog)
    
    # This works per-tile, so throw out objects that aren't in it
    t0=time.time()
    catalog=catalogs.getCatalogWithinImage(catalog, shape, wcs)
    t1=time.time()
    if reportTimingInfo: print("makeModelImage - getting catalog within image - took %.3f sec" % (t1-t0))

    # Optional SNR cuts
    if 'SNR' in catalog.keys():
        SNRKey='SNR'
    elif 'fixed_SNR' in catalog.keys():
        SNRKey='fixed_SNR'
    else:
        SNRKey=None
    if SNRKey is not None:
        catalog=catalog[catalog[SNRKey] > minSNR]

    # If we want to restrict painting to just area mask within in a tile
    # (avoids double painting of objects in overlap areas)
    if validAreaSection is not None and len(catalog) > 0:
        t0=time.time()
        x0, x1, y0, y1=validAreaSection
        coords=wcs.wcs2pix(catalog['RADeg'], catalog['decDeg'])
        x=np.array(coords)[:, 0]
        y=np.array(coords)[:, 1]
        xMask=np.logical_and(x >= x0, x < x1)
        yMask=np.logical_and(y >= y0, y < y1)
        cMask=np.logical_and(xMask, yMask)
        catalog=catalog[cMask]
        t1=time.time()
        if reportTimingInfo: print("makeModelImage - cutting catalog to area mask - took %.3f sec" % (t1-t0))

    if len(catalog) == 0:
        return None

    if cosmoModel is None:
        cosmoModel=signals.fiducialCosmoModel

    # Set initial max size in degrees from beam file (used for sources; clusters adjusted for each object)
    t0=time.time()
    numFWHM=5.0
    beam=signals.BeamProfile(beamFileName = beamFileName)
    maxSizeDeg=(beam.FWHMArcmin*numFWHM)/60
    t1=time.time()
    if reportTimingInfo: print("makeModelImage - set up beam - took %.3f sec" % (t1-t0))
    
    t0=time.time()
    if 'y_c' in catalog.keys() or 'true_y_c' in catalog.keys():
        # Clusters - insert one at a time (with different scales etc.)
        # We could use this to replace how GNFWParams are fed in also (easier for nemoModel script)
        if profile == 'A10':
            makeClusterSignalMap=signals.makeArnaudModelSignalMap
        elif profile == 'B12':
            makeClusterSignalMap=signals.makeBattagliaModelSignalMap
        else:
            raise Exception("Didn't understand profile - should be A10 or B12. This would be an excellent place\
                            to accept a string of GNFW parameters, but that is not implemented yet.")
        count=0
        # First bit here (override) is for doing injection sims faster
        if override is not None:
            z=override['redshift']
            M500=override['M500']
            y0ToInsert=catalog['y_c'].data*1e-4
            RAs=catalog['RADeg'].data
            decs=catalog['decDeg'].data
            theta500Arcmin=signals.calcTheta500Arcmin(z, M500, cosmoModel)
            maxSizeDeg=5*(theta500Arcmin/60)
            modelMap=makeClusterSignalMap(z, M500, modelMap.shape, wcs, RADeg = RAs,
                                          decDeg = decs, beam = beam,
                                          GNFWParams = GNFWParams, amplitude = y0ToInsert,
                                          maxSizeDeg = maxSizeDeg, convolveWithBeam = True,
                                          cosmoModel = cosmoModel)
            if obsFreqGHz is not None:
                modelMap=convertToDeltaT(modelMap, obsFrequencyGHz = obsFreqGHz,
                                         TCMBAlpha = TCMBAlpha, z = z)
        else:
            for row in catalog:
                count=count+1
                if 'true_M500c' in catalog.keys():
                    # This case is for when we're running from nemoMock output
                    # Since the idea of this is to create noise-free model images, we must use true values here
                    # (to avoid any extra scatter/selection effects after adding model clusters to noise maps).
                    M500=row['true_M500c']*1e14
                    z=row['redshift']
                    y0ToInsert=row['true_y_c']*1e-4
                else:
                    # NOTE: This case is for running from nemo output
                    # We need to adapt this for when the template names are not in this format
                    if 'template' not in catalog.keys():
                        raise Exception("No M500, z, or template column found in catalog.")
                    bits=row['template'].split("#")[0].split("_")
                    M500=float(bits[1][1:].replace("p", "."))
                    z=float(bits[2][1:].replace("p", "."))
                    y0ToInsert=row['y_c']*1e-4  # or fixed_y_c...
                theta500Arcmin=signals.calcTheta500Arcmin(z, M500, cosmoModel)
                maxSizeDeg=5*(theta500Arcmin/60)
                # Updated in place
                makeClusterSignalMap(z, M500, modelMap.shape, wcs, RADeg = row['RADeg'],
                                     decDeg = row['decDeg'], beam = beam,
                                     GNFWParams = GNFWParams, amplitude = y0ToInsert,
                                     maxSizeDeg = maxSizeDeg, convolveWithBeam = True,
                                     cosmoModel = cosmoModel, omap = modelMap,
                                     obsFrequencyGHz = obsFreqGHz, TCMBAlpha = TCMBAlpha)
    else:
        # Sources - slower but more accurate way
        for row in catalog:
            if validAreaSection is not None:
                x0, x1, y0, y1=validAreaSection
                x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
                if (x >= x0 and x < x1 and y >= y0 and y < y1) == False:
                    continue
            degreesMap=np.ones(modelMap.shape, dtype = float)*1e6 # NOTE: never move this
            degreesMap, xBounds, yBounds=makeDegreesDistanceMap(degreesMap, wcs,
                                                                row['RADeg'], row['decDeg'],
                                                                maxSizeDeg)
            signalMap=signals.makeBeamModelSignalMap(degreesMap, wcs, beam)*row['deltaT_c']
            modelMap=modelMap+signalMap
    t1=time.time()
    if reportTimingInfo: print("makeModelImage - painting objects - took %.3f sec" % (t1-t0))

    # Optional: apply pixel window function - generally this should be True
    # (because the source-insertion routines in signals.py interpolate onto the grid rather than average)
    if applyPixelWindow == True:
        t0=time.time()
        modelMap=enmap.apply_window(modelMap, pow = 1.0)
        t1=time.time()
        if reportTimingInfo: print("makeModelImage - pix win application - took %.3f sec" % (t1-t0))

    return modelMap
        
#------------------------------------------------------------------------------------------------------------
def sourceInjectionTest(config):
    """Insert sources with known positions and properties into the map, apply the filter, and record their
    offset with respect to the true location as a function of S/N (for the fixed reference scale only).
    If the inserted sources are clusters, the Q function will be applied to the output fluxes, to account 
    for any mismatch between the reference filter scale and the inserted clusters.
    
    Writes output to the diagnostics/ directory.
    
    Args:
        config (:obj:`nemo.startUp.NemoConfig`): Nemo configuration object.

    Returns:
        An astropy Table containing recovered position offsets and fluxes versus fixed_SNR for inserted
        sources.

    Note:
        Injection tests for clusters use the reference filter only (set with the `photFilter` keyword
        in the config). Input amplitudes for clusters are in `y_c`, while output is in `fixed_y_c`
        (because the reference filter is used). Similarly, output SNR is `fixed_SNR`, although the
        output column is labelled as `SNR`.
    
    """

    # WARNING: For multi-pass mode, this has the desired behaviour IF this is called after a nemo run
    # (i.e., the config is already set to the last filterSet)
    # But, can we guarantee that? Probably not. But we could put a warning in the docs for now?

    # This should perhaps be a config parameter
    realExclusionRadiusArcmin=5.0

    # This should make it quicker to generate test catalogs (especially when using tiles)
    selFn=completeness.SelFn(config.selFnDir, 4.0, configFileName = config.configFileName,
                             enableCompletenessCalc = False, setUpAreaMask = True,
                             tileNames = config.allTileNames)
    
    print(">>> Position recovery test [rank = %d] ..." % (config.rank))

    if 'sourceInjectionIterations' not in config.parDict.keys():
        numIterations=1
    else:
        numIterations=config.parDict['sourceInjectionIterations']

    # Change to previous behavior - if config doesn't specify models to use, assume it's point sources
    if 'sourceInjectionModels' in config.parDict.keys():
        clusterMode=True
        sourceInjectionModelList=config.parDict['sourceInjectionModels']
        SNRCol='SNR'
        fluxCol='y_c'
        noiseLevelCol='err_y_c'
        for sourceInjectionModel in sourceInjectionModelList:
            theta500Arcmin=signals.calcTheta500Arcmin(sourceInjectionModel['redshift'],
                                                      sourceInjectionModel['M500'],
                                                      signals.fiducialCosmoModel)
            label='%.2f' % (theta500Arcmin)
            sourceInjectionModel['label']=label
            sourceInjectionModel['theta500Arcmin']=theta500Arcmin
    else:
        # Sources
        clusterMode=False
        sourceInjectionModelList=[{'label': 'pointSource'}]
        SNRCol='SNR'
        fluxCol='deltaT_c'
        noiseLevelCol='err_deltaT_c'

    # This isn't really important as avoidance radius will stop us putting in too many sources
    if 'sourcesPerTile' not in config.parDict.keys():
        numSourcesPerTile=300
    else:
        numSourcesPerTile=config.parDict['sourcesPerTile']
    
    # We need the actual catalog to throw out spurious 'recoveries'
    # i.e., we only want to cross-match with objects we injected
    catFileName=config.rootOutDir+os.path.sep+"%s_optimalCatalog.fits" % (os.path.split(config.rootOutDir)[-1])
    if os.path.exists(catFileName) == False:
        raise Exception("Catalog file '%s' not found - needed to do source injection test." % (catFileName))
    realCatalog=atpy.Table().read(catFileName)
    
    # Run each scale / model and then collect everything into one table afterwards
    # NOTE: These dictionaries contain recovered measurements from running the finder
    RADegDict={}
    decDegDict={}
    SNRDict={}
    rArcminDict={}
    inFluxDict={}
    outFluxDict={}
    noiseLevelDict={}
    tileNamesDict={}
    # NOTE: This list collects all the input catalogs
    allInputCatalogs=[]
    modelCount=0
    for sourceInjectionModel in sourceInjectionModelList:
        modelCount=modelCount+1
        print(">>> Source injection model: %d/%d" % (modelCount, len(sourceInjectionModelList)))
        RADegDict[sourceInjectionModel['label']]=[]
        decDegDict[sourceInjectionModel['label']]=[]
        SNRDict[sourceInjectionModel['label']]=[]
        rArcminDict[sourceInjectionModel['label']]=[]
        inFluxDict[sourceInjectionModel['label']]=[]
        outFluxDict[sourceInjectionModel['label']]=[]
        noiseLevelDict[sourceInjectionModel['label']]=[]
        tileNamesDict[sourceInjectionModel['label']]=[]
        for i in range(numIterations):        
            print(">>> Source injection and recovery test %d/%d [rank = %d]" % (i+1, numIterations, config.rank))

            # NOTE: This block below should be handled when parsing the config file - fix/remove
            # Optional override of default GNFW parameters (used by Arnaud model), if used in filters given
            if 'GNFWParams' not in list(config.parDict.keys()):
                config.parDict['GNFWParams']='default'
            for filtDict in config.parDict['mapFilters']:
                filtDict['params']['GNFWParams']=config.parDict['GNFWParams']
            
            # We don't want to save/cache position recovery test maps
            for filtDict in config.parDict['mapFilters']:
                keysToFalsify=['saveFilteredMaps', 'savePlots']
                for key in keysToFalsify:
                    filtDict['params'][key]=False
            
            # Delete all non-reference scale filters (otherwise we'd want to cache all filters for speed)
            # NOTE: As it stands, point-source only runs may not define photFilter - we need to handle that
            # That should be obvious, as mapFilters will only have one entry
            for filtDict in config.parDict['mapFilters']:
                if filtDict['label'] == config.parDict['photFilter']:
                    break
            config.parDict['mapFilters']=[filtDict]
                
            # Filling maps with injected sources will be done when maps.preprocessMapDict is called by the filter object
            # So, we only generate the catalog here
            print("... generating mock catalog")
            if config.rank == 0:
                if filtDict['class'].find("ArnaudModel") != -1:
                    if 'sourceInjectionAmplitudeRange' not in config.parDict.keys():
                        amplitudeRange=[0.001, 10]
                    else:
                        amplitudeRange=config.parDict['sourceInjectionAmplitudeRange']
                        if amplitudeRange == 'auto':
                            amplitudeRange=[realCatalog['fixed_y_c'].min()*0.5, realCatalog['fixed_y_c'].max()]
                    if 'sourceInjectionDistribution' not in config.parDict.keys():
                        distribution='linear'
                    else:
                        distribution=config.parDict['sourceInjectionDistribution']
                    # Quick test catalog - takes < 1 sec to generate
                    mockCatalog=catalogs.generateTestCatalog(config, numSourcesPerTile,
                                                            amplitudeColumnName = fluxCol,
                                                            amplitudeRange = amplitudeRange,
                                                            amplitudeDistribution = distribution,
                                                            selFn = selFn, maskDilationPix = 20)
                    # Or... proper mock, but this takes ~24 sec for E-D56
                    #mockCatalog=pipelines.makeMockClusterCatalog(config, writeCatalogs = False, verbose = False)[0]
                    injectSources={'catalog': mockCatalog, 'GNFWParams': config.parDict['GNFWParams'],
                                   'override': sourceInjectionModel, 'profile': 'A10'}
                elif filtDict['class'].find("Beam") != -1:
                    if 'sourceInjectionAmplitudeRange' not in config.parDict.keys():
                        amplitudeRange=[1, 1000]
                    else:
                        amplitudeRange=config.parDict['sourceInjectionAmplitudeRange']
                    if 'sourceInjectionDistribution' not in config.parDict.keys():
                        distribution='log'
                    else:
                        distribution=config.parDict['sourceInjectionDistribution']
                    mockCatalog=catalogs.generateTestCatalog(config, numSourcesPerTile,
                                                            amplitudeColumnName = fluxCol,
                                                            amplitudeRange = amplitudeRange,
                                                            amplitudeDistribution = distribution,
                                                            selFn = selFn, maskDilationPix = 20)
                    injectSources={'catalog': mockCatalog, 'override': sourceInjectionModel, 'profile': None}
                else:
                    raise Exception("Don't know how to generate injected source catalogs for filterClass '%s'" % (filtDict['class']))
                if 'theta500Arcmin' in sourceInjectionModel.keys():
                    mockCatalog['theta500Arcmin']=sourceInjectionModel['theta500Arcmin']
                allInputCatalogs.append(mockCatalog)
            else:
                injectSources=None
                mockCatalog=None

            if config.MPIEnabled == True:
                bcastInjectSources=config.comm.bcast(injectSources, root = 0)
                config.comm.barrier()
                if config.rank > 0:
                    injectSources=bcastInjectSources
                    mockCatalog=bcastInjectSources['catalog']

            for mapDict in config.unfilteredMapsDictList:
                mapDict['injectSources']=injectSources
            
            # Ideally we shouldn't have blank tiles... but if we do, skip
            if len(mockCatalog) > 0:

                # Uncomment line below if want to save filtered maps for quick and dirty debugging
                # Overwrites the original filtered map, but can compare to 'stitched' map
                # config.parDict['mapFilters'][0]['params']['saveFilteredMaps']=True

                recCatalog=pipelines.filterMapsAndMakeCatalogs(config, useCachedFilters = True,
                                                               useCachedRMSMap = True, writeAreaMask = False,
                                                               writeFlagMask = False)

                # NOTE: Below here only rank 0 really needed (could then broadcast result)

                # We should be conservative in removing potential matches with real objects
                # Because we have a huge sky area and there's no reason to risk contamination of this kind
                # Effectively this is the same as using 5' circular holes in the survey mask on real objects
                # (but actually adding the avoidance radius parameter to the test catalogs really solved this)
                if len(recCatalog) > 0:
                    recCatalog=catalogs.removeCrossMatched(recCatalog, realCatalog,
                                                           radiusArcmin = realExclusionRadiusArcmin)
                if len(recCatalog) > 0:
                    try:
                        x_mockCatalog, x_recCatalog, rDeg=catalogs.crossMatch(mockCatalog, recCatalog,
                                                                              radiusArcmin = realExclusionRadiusArcmin)
                    except:
                        raise Exception("Source injection test: cross match failed on tileNames = %s; mockCatalog length = %d; recCatalog length = %d" % (str(config.tileNames), len(mockCatalog), len(recCatalog)))

                    # Catching any crazy mismatches, writing output for debugging
                    if clusterMode == False and np.logical_and(rDeg > 1.5/60, x_recCatalog['SNR'] > 10).sum() > 0:
                        mask=np.logical_and(rDeg > 1.5/60, x_recCatalog['SNR'] > 10)
                        config.parDict['mapFilters'][0]['params']['saveFilteredMaps']=True
                        recCatalog2=pipelines.filterMapsAndMakeCatalogs(config, useCachedFilters = True,
                                                                        writeAreaMask = False, writeFlagMask = False)
                        recCatalog2=catalogs.removeCrossMatched(recCatalog2, realCatalog,
                                                                radiusArcmin = realExclusionRadiusArcmin)
                        catalogs.catalog2DS9(x_recCatalog[mask],
                                             simRootOutDir+os.path.sep+"filteredMaps"+os.path.sep+tileName+os.path.sep+"mismatch-rec.reg")
                        catalogs.catalog2DS9(x_mockCatalog[mask],
                                             simRootOutDir+os.path.sep+"filteredMaps"+os.path.sep+tileName+os.path.sep+"mismatch-input.reg",
                                             color = 'red')
                        msg="Caught recovered source at large offset - check output under %s" % (simRootOutDir+os.path.sep+"filteredMaps"+os.path.sep+tileName)
                        if config.parDict['haltOnPositionRecoveryProblem'] == True:
                            raise Exception(msg)
                        else:
                            print("... Warning: %s ..." % (msg))

                    # Store everything - analyse later
                    RADegDict[sourceInjectionModel['label']]=RADegDict[sourceInjectionModel['label']]+x_recCatalog['RADeg'].tolist()
                    decDegDict[sourceInjectionModel['label']]=decDegDict[sourceInjectionModel['label']]+x_recCatalog['decDeg'].tolist()
                    SNRDict[sourceInjectionModel['label']]=SNRDict[sourceInjectionModel['label']]+x_recCatalog[SNRCol].tolist()
                    rArcminDict[sourceInjectionModel['label']]=rArcminDict[sourceInjectionModel['label']]+(rDeg*60).tolist()
                    inFluxDict[sourceInjectionModel['label']]=inFluxDict[sourceInjectionModel['label']]+x_mockCatalog[fluxCol].tolist()
                    outFluxDict[sourceInjectionModel['label']]=outFluxDict[sourceInjectionModel['label']]+x_recCatalog[fluxCol].tolist()
                    noiseLevelDict[sourceInjectionModel['label']]=noiseLevelDict[sourceInjectionModel['label']]+x_recCatalog[noiseLevelCol].tolist()
                    tileNamesDict[sourceInjectionModel['label']]=tileNamesDict[sourceInjectionModel['label']]+x_recCatalog['tileName'].tolist()

        RADegDict[sourceInjectionModel['label']]=np.array(RADegDict[sourceInjectionModel['label']])
        decDegDict[sourceInjectionModel['label']]=np.array(decDegDict[sourceInjectionModel['label']])
        SNRDict[sourceInjectionModel['label']]=np.array(SNRDict[sourceInjectionModel['label']])
        rArcminDict[sourceInjectionModel['label']]=np.array(rArcminDict[sourceInjectionModel['label']])
        inFluxDict[sourceInjectionModel['label']]=np.array(inFluxDict[sourceInjectionModel['label']])
        outFluxDict[sourceInjectionModel['label']]=np.array(outFluxDict[sourceInjectionModel['label']])
        noiseLevelDict[sourceInjectionModel['label']]=np.array(noiseLevelDict[sourceInjectionModel['label']])
        tileNamesDict[sourceInjectionModel['label']]=np.array(tileNamesDict[sourceInjectionModel['label']])
        
    # Collecting all results into one giant table
    models=[]
    theta500s=[]
    RAs=[]
    decs=[]
    SNRs=[]
    rArcmin=[]
    inFlux=[]
    outFlux=[]
    noiseLevel=[]
    tileNames=[]
    for sourceInjectionModel in sourceInjectionModelList:
        label=sourceInjectionModel['label']
        if 'theta500Arcmin' in sourceInjectionModel.keys():
            theta500s=theta500s+[sourceInjectionModel['theta500Arcmin']]*len(SNRDict[label])
        models=models+[label]*len(SNRDict[label])
        RAs=RAs+RADegDict[label].tolist()
        decs=decs+decDegDict[label].tolist()
        SNRs=SNRs+SNRDict[label].tolist()
        rArcmin=rArcmin+rArcminDict[label].tolist()
        inFlux=inFlux+inFluxDict[label].tolist()
        outFlux=outFlux+outFluxDict[label].tolist()
        noiseLevel=noiseLevel+noiseLevelDict[label].tolist()
        tileNames=tileNames+tileNamesDict[label].tolist()
    resultsTable=atpy.Table()
    resultsTable.add_column(atpy.Column(RAs, 'RADeg'))
    resultsTable.add_column(atpy.Column(decs, 'decDeg'))
    resultsTable.add_column(atpy.Column(models, 'sourceInjectionModel'))
    if len(theta500s) == len(resultsTable):
        resultsTable.add_column(atpy.Column(theta500s, 'theta500Arcmin'))
    resultsTable.add_column(atpy.Column(SNRs, SNRCol))
    resultsTable.add_column(atpy.Column(rArcmin, 'rArcmin'))
    resultsTable.add_column(atpy.Column(inFlux, 'inFlux'))
    resultsTable.add_column(atpy.Column(outFlux, 'outFlux'))
    resultsTable.add_column(atpy.Column(noiseLevel, 'noiseLevel'))
    resultsTable.add_column(atpy.Column(tileNames, 'tileName'))

    # Store the giant combined input catalog as well, for completeness calculations
    # NOTE: Not all objects in this may have been injected (masking, avoiding overlap etc.)
    if config.rank == 0:
        allInputTab=atpy.vstack(allInputCatalogs)
        allInputTab.rename_column(fluxCol, "inFlux")
        allInputTab=catalogs.removeCrossMatched(allInputTab, realCatalog, radiusArcmin = realExclusionRadiusArcmin)
        allInputTab.write(config.selFnDir+os.path.sep+"sourceInjectionInputCatalog.fits", overwrite = True)

    # Restore the original config parameters (which we overrode here)
    config.restoreConfig()

    return resultsTable

#------------------------------------------------------------------------------------------------------------
def positionRecoveryAnalysis(posRecTable, plotFileName, percentiles = [50, 95, 99.7],
                             sourceInjectionModel = None, plotRawData = True, rawDataAlpha = 1,
                             pickleFileName = None, selFnDir = None):
    """Estimate and plot position recovery accuracy as function of fixed filter scale S/N (fixed_SNR), using
    the contents of posRecTable (see positionRecoveryTest).
    
    Args:
        posRecTable (:obj:`astropy.table.Table`): Table containing recovered position offsets versus SNR
            or fixed_SNR for various cluster/source models (produced by sourceInjectionTest).
        plotFileName (str): Path where the plot file will be written.
        percentiles (list, optional): List of percentiles to plot (some interpolation will be done) and
            for which corresponding model fit parameters will be saved (if selFnDir is not None).
        sourceInjectionModel (str, optional): If given, select only objects matching the given source
            injection model name from the input table. This can be used to get results for individual
            cluster scales, for example.
        plotRawData (bool, optional): Plot the raw (fixed_SNR, positional offset) data in the background.
        pickleFileName (string, optional): Saves the percentile contours data as a pickle file if not None.
            This is saved as a dictionary with top-level keys named according to percentilesToPlot.
        selFnDir (string, optional): If given, model fit parameters will be written to a file named
            posRecModelParameters.txt under the given selFn directory path.
            
    """

    # Sources or clusters table?
    tab=posRecTable
    if len(tab) == 0:
        return None
    if np.unique(tab['sourceInjectionModel'])[0] == 'pointSource':
        SNRCol='SNR'
        plotSNRLabel="SNR"
        rArcminThreshold=np.linspace(0, 5, 201)
        plotUnits="arcsec"
        plotUnitsMultiplier=60
        plotUnitsLabel="$^{\prime\prime}$"
    else:
        # Clusters - SNR is really fixed_SNR here, because injection sims use only ref filter
        SNRCol='SNR'
        plotSNRLabel="fixed_SNR"
        rArcminThreshold=np.linspace(0, 10, 101)
        plotUnits="arcmin"
        plotUnitsMultiplier=1
        plotUnitsLabel="$^\prime$"
    
    # Optional cut on injected signal model
    if sourceInjectionModel is not None:
        tab=tab[tab['sourceInjectionModel'] == str(sourceInjectionModel)]
    
    # Evaluate %-age of sample in bins of SNR within some rArcmin threshold
    # No longer separating by input model (clusters are all shapes anyway)
    SNREdges=np.linspace(3.0, 10.0, 36)#np.linspace(0, 10, 101)
    SNRCentres=(SNREdges[1:]+SNREdges[:-1])/2.
    grid=np.zeros([rArcminThreshold.shape[0], SNREdges.shape[0]-1])
    totalGrid=np.zeros(grid.shape)
    withinRGrid=np.zeros(grid.shape)
    for i in range(SNREdges.shape[0]-1):
        SNRMask=np.logical_and(tab[SNRCol] >= SNREdges[i], tab[SNRCol] < SNREdges[i+1])
        for j in range(rArcminThreshold.shape[0]):
            total=SNRMask.sum()
            withinR=(tab['rArcmin'][SNRMask] < rArcminThreshold[j]).sum()
            totalGrid[j, i]=total
            withinRGrid[j, i]=withinR
            if total > 0:
                grid[j, i]=withinR/total
    
    # What we want are contours of constant prob - easiest to get this via matplotlib
    levelsList=np.array(percentiles)/100.
    contours=plt.contour(SNRCentres, rArcminThreshold, grid, levels = levelsList)
    minSNR=SNRCentres[np.sum(grid, axis = 0) > 0].min()
    maxSNR=SNRCentres[np.sum(grid, axis = 0) > 0].max()
    plt.close()
    
    # We make our own plot so we use consistent colours, style (haven't fiddled with contour rc settings)
    plotSettings.update_rcParams()
    plt.figure(figsize=(9,6.5))
    ax=plt.axes([0.11, 0.11, 0.88, 0.87])
    if plotRawData == True:
        plt.plot(posRecTable[SNRCol], posRecTable['rArcmin']*plotUnitsMultiplier, 
                 '.', color = '#A0A0A0', alpha = rawDataAlpha)
    contoursDict={}
    for i in range(len(levelsList)):
        vertices=contours.collections[i].get_paths()[0].vertices
        SNRs=vertices[:, 0]
        rArcminAtProb=vertices[:, 1]
        labelStr="%.1f" % (percentiles[i]) + "%"
        contoursDict[labelStr]={SNRCol: SNRs, 'rArcmin': rArcminAtProb}
        plt.plot(SNRs, rArcminAtProb*plotUnitsMultiplier, label = labelStr, lw = 3)
    plt.xlim(minSNR, maxSNR)
    #plt.ylim(0, 5)
    #plt.ylim(0,3)
    plt.legend(loc = 'upper right')
    plt.xlabel(plotSNRLabel)
    plt.ylabel("Recovered Position Offset (%s)" % (plotUnitsLabel))
    plt.savefig(plotFileName)
    plt.close()
    
    # Save %-ile contours in case we want to use them in some modelling later
    if pickleFileName is not None:
        with open(pickleFileName, "wb") as pickleFile:
            pickler=pickle.Pickler(pickleFile)
            pickler.dump(contoursDict)
    
    # Fit and save a position recovery model under selFn directory
    if selFnDir is not None:
        # This extra plot isn't really necessary
        outDir, fileName=os.path.split(os.path.abspath(plotFileName))
        fitPlotFileName=outDir+os.path.sep+"modelFits_"+fileName
        keys=contoursDict.keys()
        fitParamsDict={}
        plotSettings.update_rcParams()
        plt.figure(figsize=(9,6.5), dpi = 300)
        ax=plt.axes([0.11, 0.11, 0.88, 0.87])
        if plotRawData == True:
            posRecTable=tab
            plt.plot(posRecTable[SNRCol], posRecTable['rArcmin']*plotUnitsMultiplier, 
                     '.', color = '#A0A0A0', alpha = rawDataAlpha)
        for key in keys:
            a=contoursDict[key]
            valid=np.where(a[SNRCol] >= 4.1)
            snr=a[SNRCol][valid]
            rArcmin=a['rArcmin'][valid]
            try:
                results=optimize.curve_fit(catalogs._posRecFitFunc, snr, rArcmin)
            except:
                print("... WARNING: curve_fit failed for key = %s ..." % (key))
                continue
            bestFitSNRFold, bestFitPedestal, bestFitNorm=results[0]
            fitParamsDict[key]=np.array([bestFitSNRFold, bestFitPedestal, bestFitNorm])
            fitSNRs=np.linspace(4, 10, 100)
            plt.plot(fitSNRs, 
                     catalogs._posRecFitFunc(fitSNRs, bestFitSNRFold, bestFitPedestal, bestFitNorm)*plotUnitsMultiplier, 
                     '-', label = key)
        #plt.ylim(0, 3)
        plt.legend(loc = 'upper right')
        plt.xlim(snr.min(), snr.max())
        plt.xlabel(plotSNRLabel)
        plt.ylabel("Recovered Position Offset (%s)" % (plotUnitsLabel))
        plt.savefig(fitPlotFileName)
        plt.close()
        # Save the fits
        outFileName=selFnDir+os.path.sep+"posRecModelFits.pkl"
        with open(outFileName, "wb") as pickleFile:
            pickler=pickle.Pickler(pickleFile)
            pickler.dump(fitParamsDict)

#------------------------------------------------------------------------------------------------------------
def noiseBiasAnalysis(sourceInjTable, plotFileName, sourceInjectionModel = None):
    """Estimate the noise bias from the ratio of input to recovered flux as a function of signal-to-noise.
    
    Args:
        posRecTable (:obj:`astropy.table.Table`): Table containing recovered position offsets versus fixed_SNR 
            for various cluster/source models (produced by sourceInjectionTest).
        plotFileName (str): Path where the plot file will be written.
        clipPercentile (float, optional): Clips offset values outside of this percentile of the whole 
            position offsets distribution, to remove a small number of outliers (spurious next-neighbour 
            cross matches) that otherwise bias the contours high for large (99%+) percentile cuts in 
            individual fixed_SNR bins.
        sourceInjectionModel (str, optional): If given, restrict analysis to only objects matching this.
    
    Notes:
        For clusters, bear in mind this only makes sense if any mismatch between the inserted cluster's 
        shape and the signal assumed by the filter is taken into account. This is done using the Q-function
        in sourceInjectionTest.
            
    """

    print("Work in progress - skipped")
    return None
        
#---------------------------------------------------------------------------------------------------
def saveFITS(outputFileName, mapData, wcs, compressionType = None):
    """Writes a map (2d image array) to a new FITS file.
    
    Args:
        outputFileName (str): Filename of output FITS image.
        mapData (:obj:`np.ndarray`): Map data array.
        wcs (:obj:`astWCS.WCS`): Map WCS object.
        compressionType (str, optional): If given, the data will be compressed using the given
            method (as understood by :mod:`astropy.io.fits`). Use `PLIO_1` for masks,
            and `RICE_1` for other image data that can stand lossy compression. If None,
            the image data is not compressed.
    
    """
    
    wcs.header['NEMOVER']=nemo.__version__
    
    if os.path.exists(outputFileName):
        os.remove(outputFileName)
    
    if compressionType is None:
        if wcs is not None:
            hdu=pyfits.PrimaryHDU(mapData, wcs.header)
        else:
            hdu=pyfits.PrimaryHDU(mapData, None)
    
    if compressionType is not None:
        if wcs is not None:
            if compressionType == 'PLIO_1':
                dtype=np.uint8
            else:
                dtype=np.float32
            hdu=pyfits.CompImageHDU(np.array(mapData, dtype = dtype), wcs.header, 
                                    compression_type = compressionType)
        else:
            hdu=pyfits.CompImageHDU(np.array(mapData, dtype = dtype), None,
                                    compression_type = compressionType)
            
    newImg=pyfits.HDUList()
    newImg.append(hdu)
    newImg.writeto(outputFileName)
    newImg.close()

#---------------------------------------------------------------------------------------------------
def makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, maxDistDegrees):
    """Fills (in place) the 2d array degreesMap with distance in degrees from the given position,
    out to some user-specified maximum distance.

    Args:
        degreesMap (:obj:`np.ndarray`): Map (2d array) that will be filled with angular distance
            from the given coordinates. Probably you should feed in an array set to some extreme
            initial value (e.g., 1e6 everywhere) to make it easy to filter for pixels near the
            object coords afterwards.
        wcs (:obj:`astWCS.WCS`): WCS corresponding to degreesMap.
        RADeg (float): RA in decimal degrees of position of interest (e.g., object location).
        decDeg (float): Declination in decimal degrees of position of interest (e.g., object
            location).
        maxDistDegrees: The maximum radius out to which distance will be calculated.

    Returns:
        A map (2d array) of distance in degrees from the given position,
        (min x, max x) pixel coords corresponding to maxDistDegrees box,
        (min y, max y) pixel coords corresponding to maxDistDegrees box

    Note:
        This routine measures the pixel scale local to the given position, then assumes that it
        does not change. So, this routine may only be accurate close to the given position,
        depending upon the WCS projection used.

    """

    x0, y0=wcs.wcs2pix(RADeg, decDeg)
    ra0, dec0=RADeg, decDeg
    ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
    xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)

    xDistPix=int(round((maxDistDegrees)/xPixScale))
    yDistPix=int(round((maxDistDegrees)/yPixScale))

    Y=degreesMap.shape[0]
    X=degreesMap.shape[1]

    minX=int(round(x0))-xDistPix
    maxX=int(round(x0))+xDistPix
    minY=int(round(y0))-yDistPix
    maxY=int(round(y0))+yDistPix
    if minX < 0:
        minX=0
    if maxX > X:
        maxX=X
    if minY < 0:
        minY=0
    if maxY > Y:
        maxY=Y

    xDeg=(np.arange(degreesMap.shape[1])-x0)*xPixScale
    yDeg=(np.arange(degreesMap.shape[0])-y0)*yPixScale
    for i in range(minY, maxY):
        degreesMap[i][minX:maxX]=np.sqrt(yDeg[i]**2+xDeg[minX:maxX]**2)

    return degreesMap, [minX, maxX], [minY, maxY]

#---------------------------------------------------------------------------------------------------
def makeExtendedSourceMask(config, tileName):
    """Find extended sources in all maps, adding an extended mask to the Nemo config. Each frequency
    map will then have extended mask holes filled when preprocess is called.

    """

    settings=config.parDict['findAndMaskExtended']

    maskCube=[]
    for mapDict in config.unfilteredMapsDictList:
        data, wcs=mapDict.loadTile('mapFileName', tileName, returnWCS = True)
        weights=mapDict.loadTile('weightsFileName', tileName)
        validMask=np.nonzero(weights)
        whiteNoiseLevel=np.zeros(weights.shape)
        whiteNoiseLevel[validMask]=1/np.sqrt(weights[validMask]) # Assumed inverse variance
        # Isolate a scale that's extended
        s1=subtractBackground(data, wcs, smoothScaleDeg = settings['bigScaleDeg'])
        s2=subtractBackground(data, wcs, smoothScaleDeg = settings['smallScaleDeg'])
        s=s1-s2
        del s1, s2
        # Make a simple global 3-sigma clipped noise estimate from the filtered map
        # Then scale that according to the white noise level map from the map maker
        # Assume that median white noise level there should correspond with our global clipped noise estimate
        # (we were using mean but that blows up in edge tiles)
        mean=0
        sigma=1e6
        vals=s.flatten()
        for i in range(10):
            mask=np.less(abs(vals-mean), 3*sigma)
            mean=np.mean(vals[mask])
            sigma=np.std(vals[mask])
        scaleFactor=sigma/np.median(whiteNoiseLevel[validMask])
        whiteNoiseLevel[validMask]=whiteNoiseLevel[validMask]*scaleFactor
        snr=np.zeros(s.shape)
        snr[validMask]=s[validMask]/whiteNoiseLevel[validMask]
        # Mask set such that 1 = masked, 0 = not masked
        extendedMask=np.array(np.greater(snr, settings['thresholdSigma']), dtype = np.uint8)
        if 'dilationPix' in settings.keys() and settings['dilationPix'] > 0:
            for i in range(settings['dilationPix']):
                extendedMask=mahotas.dilate(extendedMask)
        extendedMask[extendedMask > 0]=1
        maskCube.append(extendedMask)
    maskCube=np.array(maskCube, dtype = np.uint8)
    extendedMask=maskCube.sum(axis = 0)
    extendedMask[extendedMask > 0]=1

    # Optionally cut any small objects
    if 'minSizeArcmin2' in settings.keys() and settings['minSizeArcmin2'] > 0:
        arcmin2Map=getPixelAreaArcmin2Map(extendedMask.shape, wcs)
        segMap, numObjects=ndimage.label(extendedMask)
        for i in range(1, numObjects+1):
            if arcmin2Map[segMap == i].sum() < settings['minSizeArcmin2']:
                extendedMask[segMap == i]=0

    os.makedirs(config.diagnosticsDir+os.path.sep+"extendedMask", exist_ok = True)
    outFileName=config.diagnosticsDir+os.path.sep+"extendedMask"+os.path.sep+tileName+".fits"
    saveFITS(outFileName, extendedMask, wcs, compressionType = 'PLIO_1')

    for mapDict in config.unfilteredMapsDictList:
        mapDict['extendedMask']=config.diagnosticsDir+os.path.sep+"extendedMask"

#------------------------------------------------------------------------------------------------------------
def makeMaskFromDS9PolyRegionFile(regionFileName, shape, wcs):
    """Make a mask from a DS9 region file. The region file must have been created with RA, dec coordinates
    given in decimal degrees, and the shapes defining the mask must consist of polygon regions only.

    Args:
        regionFileName (:obj:`str`): Path to SAOImage DS9 region file.
        origShape (:obj:`tuple`): Shape of the output mask.
        origWCS (:obj:`astWCS.WCS object`): WCS for the output mask.

    Returns:
        Mask (2d array)

    """

    with open(regionFileName, "r") as inFile:
        lines=inFile.readlines()
    polyList=[]
    for line in lines:
        if line.find("polygon") != -1:
            polyPoints=[]
            coords=line.split("polygon(")[-1].split(") ")[0].split(",")
            for i in range(0, len(coords), 2):
                try:
                    RADeg, decDeg=[float(coords[i]), float(coords[i+1])]
                except:
                    raise Exception("failed to parse coords in region file %s - problem at: %s" % (regionFileName, coords))
                x, y=wcs.wcs2pix(RADeg, decDeg)
                polyPoints.append((int(round(y)), int(round(x))))
            polyList.append(polyPoints)
    surveyMask=np.zeros(shape, dtype = int)
    for polyPoints in polyList:
        mahotas.polygon.fill_polygon(polyPoints, surveyMask)

    return surveyMask
