"""

This module contains tools for manipulating maps (e.g., conversion of units etc.).

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
import pyximport; pyximport.install()
import nemoCython
import time
import shutil
import copy
import yaml
import pickle
#import IPython
from pixell import enmap
import nemo
from . import catalogs
from . import signals
from . import photometry
from . import plotSettings
from . import pipelines
from . import completeness
np.random.seed()
              
#-------------------------------------------------------------------------------------------------------------
def convertToY(mapData, obsFrequencyGHz = 148):
    """Converts mapData (in delta T) at given frequency to y.
    
    """
    fx=signals.fSZ(obsFrequencyGHz)    
    mapData=(mapData/(signals.TCMB*1e6))/fx # remember, map is in deltaT uK 
    
    return mapData

#-------------------------------------------------------------------------------------------------------------
def convertToDeltaT(mapData, obsFrequencyGHz = 148):
    """Converts mapData (in yc) to deltaT (micro Kelvin) at given frequency.
    
    """
    fx=signals.fSZ(obsFrequencyGHz)   
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
        makeTileDir will expand tiles by a user-specified amount such that they overlap.
    
    """
    
    # This deals with identifying boss vs. full AdvACT footprint maps 
    mapCentreRA, mapCentreDec=wcs.getCentreWCSCoords()    
    skyWidth, skyHeight=wcs.getFullSizeSkyDeg()
    if mapCentreRA < 0.1 and skyWidth < 0.1 or skyWidth > 359.9:
        handle180Wrap=True
    else:
        handle180Wrap=False
    
    segMap=surveyMask
    segMap, numObjects=ndimage.label(np.greater(segMap, 0))
    fieldIDs=np.arange(1, numObjects+1)
        
    tileList=[]
    for f in fieldIDs:
        ys, xs=np.where(segMap == f)
        if len(ys) < 1000:  # In case of stray individual pixels (e.g., combined with extended sources mask)
            continue
        yMin=ys.min()
        yMax=ys.max()
        xc=int((xs.min()+xs.max())/2)
        RAc, decMin=wcs.pix2wcs(xc, yMin)
        RAc, decMax=wcs.pix2wcs(xc, yMax)
        
        numRows=int((decMax-decMin)/targetTileHeight)
        tileHeight=np.ceil(((decMax-decMin)/numRows)*100)/100
        assert(tileHeight < 10)
        
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
            outFile.write('polygon(%d, %d, %d, %d, %d, %d, %d, %d) # text="%s"\n' % (c[0], c[2], c[0], c[3], c[1], c[3], c[1], c[2], name))  
                
#------------------------------------------------------------------------------------------------------------
def addAutoTileDefinitions(parDict, DS9RegionFileName = None, cacheFileName = None):
    """Runs autotiler to add automatic tile definitions into the parameters dictionary in-place.
    
    Args:
        parDict (:obj:`dict`): Dictionary containing the contents of the Nemo config file.
        DS9RegionFileName (str, optional): Path to DS9 regions file to be written.
        cacheFileName (str, optional): Path to output a cached .yml file which will can be read instead on 
            repeated runs (for speed).

    """
    
    if cacheFileName is not None and os.path.exists(cacheFileName):
        with open(cacheFileName, "r") as stream:
            parDict['tileDefinitions']=yaml.safe_load(stream)
        return None

    if 'tileDefinitions' in parDict.keys() and type(parDict['tileDefinitions']) == dict:
        # If we're not given a survey mask, we'll make one up from the map image itself
        if 'mask' in parDict['tileDefinitions'].keys() and parDict['tileDefinitions']['mask'] is not None:
            surveyMaskPath=parDict['tileDefinitions']['mask']
        else:
            surveyMaskPath=parDict['unfilteredMaps'][0]['mapFileName']
        with pyfits.open(surveyMaskPath) as img:    
            # Just in case RICE-compressed or similar
            if img[0].data is None:
                surveyMask=img['COMPRESSED_IMAGE'].data
                wcs=astWCS.WCS(img['COMPRESSED_IMAGE'].header, mode = 'pyfits')
            else:
                surveyMask=img[0].data
                wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
            # One day we will write a routine to deal with the multi-plane thing sensibly...
            # But today is not that day
            if surveyMask.ndim == 3:
                surveyMask=surveyMask[0, :]
            assert(surveyMask.ndim == 2)   
            surveyMask[surveyMask != 0]=1
        parDict['tileDefinitions']=autotiler(surveyMask, wcs, 
                                             parDict['tileDefinitions']['targetTileWidthDeg'],
                                             parDict['tileDefinitions']['targetTileHeightDeg'])
        print("... breaking map into %d tiles ..." % (len(parDict['tileDefinitions'])))
        if DS9RegionFileName is not None:
            saveTilesDS9RegionsFile(parDict, DS9RegionFileName)
    
        if cacheFileName is not None:
            stream=yaml.dump(parDict['tileDefinitions'])
            with open(cacheFileName, "w") as outFile: 
                outFile.write(stream) 

#-------------------------------------------------------------------------------------------------------------
def loadTile(pathToTileImages, tileName, returnWCS = False):
    """Given a path to a directory full of tiles, or a .fits file, return the map array and (optionally)
    the WCS.
    
    Args:
        pathToTileImages(str): Path to either a .fits file, or a directory full of .fits files named by tileName.
        tileName(str): The name of the tile to load.
        
    Returns:
        Map data (and optionally wcs)
    
    """

    if os.path.isdir(pathToTileImages) == True:
        with pyfits.open(pathToTileImages+os.path.sep+tileName+".fits") as img:
            if returnWCS == True:
                # Zapping keywords in old ACT maps that confuse astropy.wcs
                wcs=astWCS.WCS(img[0].header, mode = 'pyfits', zapKeywords = ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'])
            data=img[0].data
    else:
        with pyfits.open(pathToTileImages) as img:
            if returnWCS == True:
                wcs=astWCS.WCS(img[tileName].header, mode = 'pyfits', zapKeywords = ['PC1_1', 'PC1_2', 'PC2_1', 'PC2_2'])
            data=img[tileName].data

    if returnWCS == True:
        return data, wcs
    else:
        return data

#-------------------------------------------------------------------------------------------------------------
def makeTileDir(parDict, writeToDisk = True):
    """Update this later. Revised version. Instead of making a MEF, makes a directory for each map and puts
    individual tile images there. We'll only need to edit preprocessMapDict to handle the difference. Why
    are we doing this? Just in case reading from the same file is gumming up MPI runs on hippo when using lots
    of nodes.
    
    Makes a tileDir (directory containing .fits images, one per tile) file, if the needed parameters are 
    given in parDict. Adjusts unfilteredMapsDictList accordingly and returns it.
    
    If the options for making a tileDir image aren't given in parDict, then we pass through a standard
    single extension file (or rather the path to it, as originally given)
    
    NOTE: If the map given in unfilteredMaps is 3d (enki gives I, Q, U as a datacube), then this will extract
    only the I (temperature) part and save that in the tileDir file. This will need changing if hunting for
    polarized sources...
    
    Returns unfilteredMapsDictList [input for filterMaps], list of extension names, dictionary of clip coords
    
    NOTE: Under MPI, this should only be called by the rank = 0 process
    
    """
    
    if 'makeTileDir' not in list(parDict.keys()):
        parDict['makeTileDir']=False
    
    # Some of this is rather clunky...
    unfilteredMapsDictList=[]
    clipCoordsDict={}
    if parDict['makeTileDir'] == False:
        tileNames=[]        
        for mapDict in parDict['unfilteredMaps']:
            unfilteredMapsDictList.append(mapDict.copy())
            img=pyfits.open(mapDict['mapFileName'])
            if tileNames == []:
                for ext in img:
                    tileNames.append(ext.name)
                    clipCoordsDict[ext.name]={'clippedSection': [0, ext.header['NAXIS1'], 0, ext.header['NAXIS2']],
                                              'header': ext.header,
                                              'areaMaskInClipSection': [0, ext.header['NAXIS1'], 0, ext.header['NAXIS2']]}
            else:
                for ext in img:
                    if ext.name not in tileNames:
                        raise Exception("extension names do not match between all maps in unfilteredMapsDictList")
            img.close()
    else:
        tileNames=[]
        wcs=None
        for mapDict in parDict['unfilteredMaps']:
            if 'tileDefLabel' in list(parDict.keys()):
                tileDefLabel=parDict['tileDefLabel']
            else:
                tileDefLabel='userDefined'
            tileDirFileNameLabel="%s_%.1f" % (tileDefLabel, parDict['tileOverlapDeg'])
            
            # Figure out what the input / output files will be called
            # NOTE: we always need to make a survey mask if none exists, as used to zap over regions, so that gets special treatment
            fileNameKeys=['mapFileName', 'weightsFileName', 'pointSourceMask', 'surveyMask']
            inFileNames=[]
            outFileNames=[]
            mapTypeList=[]
            for f in fileNameKeys:
                if f in list(mapDict.keys()) and mapDict[f] is not None:
                    inFileNames.append(mapDict[f])
                    mapDir, mapFileName=os.path.split(mapDict[f])
                    if mapDir != '':
                        mapDirStr=mapDir+os.path.sep
                    else:
                        mapDirStr=''
                    outFileNames.append(mapDirStr+"tileDir_%s_" % (tileDirFileNameLabel)+mapFileName)
                    mapTypeList.append(f)

            if wcs is None:
                wcsPath=parDict['unfilteredMaps'][0]['mapFileName']
                with pyfits.open(wcsPath) as img:
                    wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
                
            # Extract tile definitions (may have been inserted by autotiler before calling here)
            tileNames=[]
            coordsList=[]
            for tileDict in parDict['tileDefinitions']:
                ra0, ra1, dec0, dec1=tileDict['RADecSection']
                x0, y0=wcs.wcs2pix(ra0, dec0)
                x1, y1=wcs.wcs2pix(ra1, dec1)
                xMin=min([x0, x1])
                xMax=max([x0, x1])
                yMin=min([y0, y1])
                yMax=max([y0, y1])
                coordsList.append([xMin, xMax, yMin, yMax])
                tileNames.append(tileDict['tileName'])   
            
            # Make tiles
            # NOTE: we accommodate having user-defined regions for calculating noise power in filters here
            # Since we would only use such an option with tileDir files, this should be okay
            # Although since we do this by modifying headers, would need to remake tileDir files each time adjusted in .par file
            # NOTE: now treating surveyMask as special, and zapping overlap regions there (simplify selection function stuff later)
            tileOverlapDeg=parDict['tileOverlapDeg']
            for mapType, inMapFileName, outMapFileName in zip(mapTypeList, inFileNames, outFileNames):
                mapData=None    # only load the map if we have to
                if writeToDisk == True:
                    os.makedirs(outMapFileName, exist_ok = True)
                for c, name in zip(coordsList, tileNames):
                    tileFileName=outMapFileName+os.path.sep+name+".fits"
                    if mapData is None:
                        #deckImg=pyfits.HDUList()
                        # Special handling for case where surveyMask = None in the .par file (tidy later...)
                        if mapType == 'surveyMask' and inMapFileName is None:
                            with pyfits.open(inFileNames[0]) as img:
                                mapData=np.ones(img[0].data.shape)
                        else:
                            with pyfits.open(inMapFileName) as img:
                                mapData=img[0].data
                        # Deal with Sigurd's maps which have T, Q, U as one 3d array
                        # If anyone wants to find polarized sources, this will need changing...
                        if mapData.ndim == 3:
                            mapData=mapData[0, :]
                    # Defining clip region
                    y0=c[2]
                    y1=c[3]
                    x0=c[0]
                    x1=c[1]
                    ra0, dec0=wcs.pix2wcs(x0, y0)
                    ra1, dec1=wcs.pix2wcs(x1, y1)
                    # Be careful with signs here... and we're assuming approx pixel size is ok
                    if x0-tileOverlapDeg/wcs.getPixelSizeDeg() > 0:
                        ra0=ra0+tileOverlapDeg
                    if x1+tileOverlapDeg/wcs.getPixelSizeDeg() < mapData.shape[1]:
                        ra1=ra1-tileOverlapDeg
                    if y0-tileOverlapDeg/wcs.getPixelSizeDeg() > 0:
                        dec0=dec0-tileOverlapDeg
                    if y1+tileOverlapDeg/wcs.getPixelSizeDeg() < mapData.shape[0]:
                        dec1=dec1+tileOverlapDeg
                    if ra1 > ra0:
                        ra1=-(360-ra1)
                    # This bit is necessary to avoid Q -> 0.2 ish problem with Fourier filter
                    # (which happens if image dimensions are both odd)
                    # I _think_ this is related to the interpolation done in signals.fitQ
                    ddec=0.5/60.
                    count=0
                    clip=astImages.clipUsingRADecCoords(mapData, wcs, ra1, ra0, dec0, dec1)
                    while clip['data'].shape[0] % 2 != 0:
                        clip=astImages.clipUsingRADecCoords(mapData, wcs, ra1, ra0, dec0-ddec*count, dec1)
                        count=count+1
                        if count > 100:
                            raise Exception("Triggered stupid bug in makeTileDir... this should be fixed properly")
                    # Storing clip coords etc. so can stitch together later
                    # areaMaskSection here is used to define the region that would be kept (takes out overlap)
                    ra0, dec0=wcs.pix2wcs(x0, y0)
                    ra1, dec1=wcs.pix2wcs(x1, y1)
                    clip_x0, clip_y0=clip['wcs'].wcs2pix(ra0, dec0)
                    clip_x1, clip_y1=clip['wcs'].wcs2pix(ra1, dec1)
                    clip_x0=int(round(clip_x0))
                    clip_x1=int(round(clip_x1))
                    clip_y0=int(round(clip_y0))
                    clip_y1=int(round(clip_y1))
                    if name not in clipCoordsDict:
                        clipCoordsDict[name]={'clippedSection': clip['clippedSection'], 'header': clip['wcs'].header,
                                              'areaMaskInClipSection': [clip_x0, clip_x1, clip_y0, clip_y1]}
                    else:
                        assert(clipCoordsDict[name]['clippedSection'] == clip['clippedSection'])
                    # Old
                    #clip=astImages.clipUsingRADecCoords(mapData, wcs, ra1, ra0, dec0, dec1)
                    print("... adding %s [%d, %d, %d, %d ; %d, %d] ..." % (name, ra1, ra0, dec0, dec1, ra0-ra1, dec1-dec0))
                    header=clip['wcs'].header#.copy()
                    if 'tileNoiseRegions' in list(parDict.keys()):
                        if name in list(parDict['tileNoiseRegions'].keys()):
                            noiseRAMin, noiseRAMax, noiseDecMin, noiseDecMax=parDict['tileNoiseRegions'][name]
                        else:
                            if 'autoBorderDeg' in parDict['tileNoiseRegions']:
                                autoBorderDeg=parDict['tileNoiseRegions']['autoBorderDeg']
                                for tileDef in parDict['tileDefinitions']:
                                    if tileDef['tileName'] == name:
                                        break
                                noiseRAMin, noiseRAMax, noiseDecMin, noiseDecMax=tileDef['RADecSection']
                                noiseRAMin=noiseRAMin+autoBorderDeg
                                noiseRAMax=noiseRAMax-autoBorderDeg
                                noiseDecMin=noiseDecMin+autoBorderDeg
                                noiseDecMax=noiseDecMax-autoBorderDeg
                            else:
                                raise Exception("No entry in tileNoiseRegions in config file for tileName '%s' - either add one, or add 'autoBorderDeg': 0.5 (or similar) to tileNoiseRegions" % (name))
                        print("... adding noise region [%.3f, %.3f, %.3f, %.3f] to header %s ..." % (noiseRAMin, noiseRAMax, noiseDecMin, noiseDecMax, name))
                        header['NRAMIN']=noiseRAMin
                        header['NRAMAX']=noiseRAMax
                        header['NDEMIN']=noiseDecMin
                        header['NDEMAX']=noiseDecMax
                    # Survey mask is special: zap overlap regions outside of tile definitions
                    if mapType == 'surveyMask':
                        clip_x0, clip_x1, clip_y0, clip_y1=clipCoordsDict[name]['areaMaskInClipSection']
                        zapMask=np.zeros(clip['data'].shape)
                        zapMask[clip_y0:clip_y1, clip_x0:clip_x1]=1.
                        clip['data']=clip['data']*zapMask
                    if os.path.exists(tileFileName) == False and writeToDisk == True:
                        astImages.saveFITS(tileFileName, clip['data'], clip['wcs'])
                                
            # Replace entries in unfilteredMapsDictList in place
            for key, outFileName in zip(mapTypeList, outFileNames):
                mapDict[key]=outFileName
            unfilteredMapsDictList.append(mapDict.copy())
        
    return unfilteredMapsDictList, tileNames, clipCoordsDict

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
def checkMask(fileName):
    """Checks whether a mask contains negative values (invalid) and throws an exception if this is the case.
    
    Args:
        fileName (str): Name of the .fits format mask file to check
        
    """
    
    with pyfits.open(fileName) as img:
        for hdu in img:
            if hdu.data is not None:
                if np.less(hdu.data, 0).sum() > 0:
                    raise Exception("Mask file '%s' contains negative values - please fix your mask." % (fileName))

#-------------------------------------------------------------------------------------------------------------
def stitchTiles(config):
    """Stitches together full size filtered maps, SN maps, area maps, and noise maps that have been previously
    been saved as tiles.
    
    """
    
    # Defining the maps to stitch together and where they will go
    stitchDictList=[{'pattern': config.filteredMapsDir+os.path.sep+"$TILENAME"+os.path.sep+"$FILTLABEL#$TILENAME_filteredMap.fits",
                     'outFileName': config.filteredMapsDir+os.path.sep+"stitched_$FILTLABEL_filteredMap.fits",
                     'compressed': False},
                    {'pattern': config.filteredMapsDir+os.path.sep+"$TILENAME"+os.path.sep+"$FILTLABEL#$TILENAME_SNMap.fits",
                     'outFileName': config.filteredMapsDir+os.path.sep+"stitched_$FILTLABEL_SNMap.fits",
                     'compressed': False},
                    {'pattern': config.selFnDir+os.path.sep+"areaMask#$TILENAME.fits",
                     'outFileName': config.selFnDir+os.path.sep+"stitched_areaMask.fits",
                     'compressed': True},
                    {'pattern': config.selFnDir+os.path.sep+"RMSMap_$FILTLABEL#$TILENAME.fits",
                     'outFileName': config.selFnDir+os.path.sep+"stitched_RMSMap_$FILTLABEL.fits",
                     'compressed': True}]

    tileCoordsDict=config.tileCoordsDict
    for filterDict in config.parDict['mapFilters']:
        if 'saveFilteredMaps' in filterDict['params'].keys() and filterDict['params']['saveFilteredMaps'] == True:
            for stitchDict in stitchDictList:
                pattern=stitchDict['pattern']
                outFileName=stitchDict['outFileName'].replace("$FILTLABEL", filterDict['label'])
                compressed=stitchDict['compressed']
                d=np.zeros([config.origWCS.header['NAXIS2'], config.origWCS.header['NAXIS1']])
                wcs=config.origWCS
                for tileName in tileCoordsDict.keys():
                    f=pattern.replace("$TILENAME", tileName).replace("$FILTLABEL", filterDict['label'])
                    if os.path.exists(f) == True:
                        # Handle compressed or otherwise
                        with pyfits.open(f) as img:
                            tileData=img[0].data
                            if tileData is None:
                                for extName in img:
                                    tileData=img[extName].data
                                    if tileData is not None:
                                        break
                            assert tileData is not None
                    else:
                        continue
                    areaMask, areaWCS=completeness.loadAreaMask(tileName, config.selFnDir)
                    minX, maxX, minY, maxY=config.tileCoordsDict[tileName]['clippedSection']
                    d[minY:maxY, minX:maxX]=d[minY:maxY, minX:maxX]+areaMask*tileData
                saveFITS(outFileName, d, wcs, compressed = compressed)

#-------------------------------------------------------------------------------------------------------------
def stitchTilesQuickLook(filePattern, outFileName, outWCS, outShape, fluxRescale = 1.0):
    """Fast routine for stitching map tiles back together. Since this uses interpolation, you probably don't 
    want to do analysis on the output - this is just for checking / making plots etc.. This routine sums 
    images as it pastes them into the larger map grid. So, if the zeroed (overlap) borders are not handled,
    correctly, this will be obvious in the output.

    NOTE: This assumes RA in x direction, dec in y direction (and CAR projection).
    
    NOTE: This routine only writes output if there are multiple files that match filePattern (to save needless
    duplicating maps if nemo was not run in tileDir mode).
    
    Output map will be multiplied by fluxRescale (this is necessary if downsampling in resolution).

    Takes 10 sec for AdvACT S16-sized downsampled by a factor of 4 in resolution.
    
    """
    
    # Set-up template blank map into which we'll dump tiles
    outData=np.zeros(outShape)
    outRACoords=np.array(outWCS.pix2wcs(np.arange(outData.shape[1]), [0]*outData.shape[1]))
    outDecCoords=np.array(outWCS.pix2wcs([0]*np.arange(outData.shape[0]), np.arange(outData.shape[0])))
    outRA=outRACoords[:, 0]
    outDec=outDecCoords[:, 1]
    RAToX=interpolate.interp1d(outRA, np.arange(outData.shape[1]), fill_value = 'extrapolate')
    DecToY=interpolate.interp1d(outDec, np.arange(outData.shape[0]), fill_value = 'extrapolate')

    # Splat tiles into output map
    inFiles=glob.glob(filePattern)
    if len(inFiles) < 1:
        return None # We could raise an Exception here instead
    count=0
    for f in inFiles:
        count=count+1
        #print("... %d/%d ..." % (count, len(inFiles)))
        with pyfits.open(f) as img:
            for hdu in img:
                if hdu.shape is not ():
                    d=hdu.data
                    inWCS=astWCS.WCS(hdu.header, mode = 'pyfits')
                    break
        xIn=np.arange(d.shape[1])
        yIn=np.arange(d.shape[0])
        inRACoords=np.array(inWCS.pix2wcs(xIn, [0]*len(xIn)))
        inDecCoords=np.array(inWCS.pix2wcs([0]*len(yIn), yIn))
        inRA=inRACoords[:, 0]
        inDec=inDecCoords[:, 1]
        xOut=np.array(RAToX(inRA), dtype = int)
        yOut=np.array(DecToY(inDec), dtype = int)
        for i in range(len(yOut)):
            try:
                outData[yOut[i]][xOut]=outData[yOut[i]][xOut]+d[yIn[i], xIn]
            except:
                raise Exception("Output pixel coords invalid - if you see this, probably outWCS.header has keywords that confuse astropy.wcs (PC1_1 etc. - in old ACT maps)")
    saveFITS(outFileName, outData*fluxRescale, outWCS, compressed = True)

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
            rRange, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, 
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
                rRange, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, mapWCS, 
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
    
#-------------------------------------------------------------------------------------------------------------
def preprocessMapDict(mapDict, tileName = 'PRIMARY', diagnosticsDir = None):
    """Applies a number of pre-processing steps to the map described by `mapDict`, prior to filtering.
    
    The first step is to load the map itself and the associated weights. Some other operations that may be 
    applied are controlled by keys added to `mapDict`. Some of these may be specified in the .yml configuration
    file, while others are applied by particular filter objects or by routines that generate simulated data. 
    The following keys are understood:
    
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
    
    applyBeamConvolution (:obj: `str`)
        If True, the map is convolved with the beam given in the beamFileName key. This should only be 
        needed when using preliminary y-maps made by tILe-C.
            
    Args:
        mapDict (:obj:`dict`): A dictionary with the same keys as given in the unfilteredMaps list in the 
            .yml configuration file (i.e., it must contain at least the keys 'mapFileName', 'units', and
            'weightsFileName', and may contain some of the optional keys listed above).
        tileName (:obj:`str`): Name of the map tile (extension name) to operate on.
        diagnosticsDir (:obj:`str`): Path to a directory where miscellaneous diagnostic data are written.
    
    Returns:
        A dictionary with keys that point to the map itself ('data'), weights ('weights'), masks 
        ('surveyMask', 'pointSourceMask'), and WCS object ('wcs').
    
    """
    
    data, wcs=loadTile(mapDict['mapFileName'], tileName, returnWCS = True)
    
    # Optional calibration factor
    if 'calibFactor' in mapDict.keys():
        data=data*mapDict['calibFactor']
    
    # For Enki maps... take only I (temperature) for now, add options for this later
    if data.ndim == 3:
        data=data[0, :]
    if mapDict['units'] == 'Jy/sr':
        if mapDict['obsFreqGHz'] == 148:
            data=(data/1.072480e+09)*2.726*1e6
        elif mapDict['obsFreqGHz'] == 219:
            data=(data/1.318837e+09)*2.726*1e6
        else:
            raise Exception("no code added to support conversion to uK from Jy/sr for freq = %.0f GHz" \
                    % (mapDict['obsFreqGHz']))

    # Load weight map if given
    if 'weightsFileName' in list(mapDict.keys()) and mapDict['weightsFileName'] is not None:
        weights=loadTile(mapDict['weightsFileName'], tileName)
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
    if 'surveyMask' in list(mapDict.keys()) and mapDict['surveyMask'] is not None:
        surveyMask=loadTile(mapDict['surveyMask'], tileName)
    else:
        surveyMask=np.ones(data.shape)
        surveyMask[weights == 0]=0
    
    # Some apodisation of the data outside the survey mask
    # NOTE: should add adjustable parameter for this somewhere later
    if 'apodizeUsingSurveyMask' in list(mapDict.keys()) and mapDict['apodizeUsingSurveyMask'] == True:
        # We need to remain unapodized to at least noiseGridArcmin beyond the edge of the survey mask
        # We'll need to make these adjustable parameters
        apodMask=np.array(surveyMask, dtype = bool)
        for i in range(120):
            apodMask=mahotas.dilate(apodMask)
        apodMask=ndimage.gaussian_filter(np.array(apodMask, dtype = float), 20)
        data=data*apodMask
        del apodMask

    if 'pointSourceMask' in list(mapDict.keys()) and mapDict['pointSourceMask'] is not None:
        psMask=loadTile(mapDict['pointSourceMask'], tileName)
    else:
        psMask=np.ones(data.shape)
                
    # Optional map clipping
    if 'RADecSection' in list(mapDict.keys()) and mapDict['RADecSection'] is not None:
        RAMin, RAMax, decMin, decMax=mapDict['RADecSection']
        clip=astImages.clipUsingRADecCoords(data, wcs, RAMin, RAMax, decMin, decMax)
        data=clip['data']
        whtClip=astImages.clipUsingRADecCoords(weights, wcs, RAMin, RAMax, decMin, decMax)
        weights=whtClip['data']
        psClip=astImages.clipUsingRADecCoords(psMask, wcs, RAMin, RAMax, decMin, decMax)
        psMask=psClip['data']
        surveyClip=astImages.clipUsingRADecCoords(surveyMask, wcs, RAMin, RAMax, decMin, decMax)
        surveyMask=surveyClip['data']
        wcs=clip['wcs']
        if len(clip['data']) == 0:
            raise Exception("Clipping using RADecSection returned empty array - check RADecSection in config .yml file is in map")
        #astImages.saveFITS(diagnosticsDir+os.path.sep+'%d' % (mapDict['obsFreqGHz'])+"_weights.fits", weights, wcs)
    
    # For source-free simulations (contamination tests)
    if 'CMBSimSeed' in list(mapDict.keys()):
        randMap=simCMBMap(data.shape, wcs, noiseLevel = 0, beamFileName = mapDict['beamFileName'], 
                          seed = mapDict['CMBSimSeed'])
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
        # Sanity check
        outFileName=diagnosticsDir+os.path.sep+"CMBSim_%d#%s.fits" % (mapDict['obsFreqGHz'], tileName) 
        astImages.saveFITS(outFileName, data, wcs)
    
    # For position recovery tests
    if 'injectSources' in list(mapDict.keys()):
        # NOTE: Need to add varying GNFWParams here
        modelMap=makeModelImage(data.shape, wcs, mapDict['injectSources']['catalog'], 
                                mapDict['beamFileName'], obsFreqGHz = mapDict['obsFreqGHz'], 
                                GNFWParams = mapDict['injectSources']['GNFWParams'],
                                override = mapDict['injectSources']['override'])
        modelMap[weights == 0]=0
        data=data+modelMap

    # Should only be needed for handling preliminary tILe-C maps
    if 'applyBeamConvolution' in mapDict.keys() and mapDict['applyBeamConvolution'] == True:
        data=convolveMapWithBeam(data, wcs, mapDict['beamFileName'], maxDistDegrees = 1.0)
        if diagnosticsDir is not None:
            astImages.saveFITS(diagnosticsDir+os.path.sep+"beamConvolved#%s.fits" % (tileName), data, wcs)
        
    # Optional masking of point sources from external catalog
    # Especially needed if using Fourier-space matched filter (and maps not already point source subtracted)
    if 'maskPointSourcesFromCatalog' in list(mapDict.keys()) and mapDict['maskPointSourcesFromCatalog'] is not None:  
        # This is fast enough if using small tiles and running in parallel...
        # If our masking/filling is effective enough, we may not need to mask so much here...
        if type(mapDict['maskPointSourcesFromCatalog']) is not list:
            mapDict['maskPointSourcesFromCatalog']=[mapDict['maskPointSourcesFromCatalog']]
        psMask=np.ones(data.shape)
        for catalogPath in mapDict['maskPointSourcesFromCatalog']:
            tab=atpy.Table().read(catalogPath)
            tab=catalogs.getCatalogWithinImage(tab, data.shape, wcs)
            # Variable sized holes: based on inspecting sources by deltaT in f150 maps
            # To avoid the problem of rings around sources, we sort the list and mask the brightest first
            # If an object lands in an already masked area, we remove it from the catalog
            # If we're given a catalog that already has rArcmin in it, we use that instead
            if 'rArcmin' not in tab.keys():
                tab.sort('deltaT_c', reverse = True)
                tab.add_column(atpy.Column(np.zeros(len(tab)), 'rArcmin'))
                tab['rArcmin'][tab['deltaT_c'] < 500]=3.0
                tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 500, tab['deltaT_c'] < 1000)]=4.0
                tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 1000, tab['deltaT_c'] < 2000)]=5.0
                tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 2000, tab['deltaT_c'] < 3000)]=5.5
                tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 3000, tab['deltaT_c'] < 10000)]=6.0
                tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 10000, tab['deltaT_c'] < 40000)]=8.0
                tab['rArcmin'][tab['deltaT_c'] >= 40000]=12.0
            selectedRows=[]
            for row in tab:
                x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
                if psMask[int(y), int(x)] != 0:
                    rArcminMap=np.ones(data.shape, dtype = float)*1e6
                    rArcminMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(rArcminMap, wcs, 
                                                                                   row['RADeg'], row['decDeg'], 
                                                                                   row['rArcmin']/60.)
                    rArcminMap=rArcminMap*60
                    psMask[rArcminMap < row['rArcmin']]=0
                    selectedRows.append(True)
                else:
                    selectedRows.append(False)
            tab=tab[np.where(selectedRows)]
            
            # Subtract sources (if there are small residuals, doesn't matter, as masked later anyway)
            # NOTE: Doesn't work very well at the moment
            #model=makeModelImage(data.shape, wcs, tab, mapDict['beamFileName'])
            #data=data-model
            # Or fill holes with smoothed map
            # NOTE: We were filling with white noise - this is bad for repeatability (since we use as noise term in filter)
            pixRad=(10.0/60.0)/wcs.getPixelSizeDeg()
            bckData=ndimage.median_filter(data, int(pixRad)) # Size chosen for max hole size... slow... but quite good
            if mapDict['weightsType'] =='invVar':
                rms=np.zeros(weights.shape)
                rms[np.nonzero(weights)]=1.0/np.sqrt(weights[np.nonzero(weights)])
            else:
                raise Exception("Not implemented white noise estimate for non-inverse variance weights for masking sources from catalog")
            data[np.where(psMask == 0)]=bckData[np.where(psMask == 0)]#+np.random.normal(0, rms[np.where(psMask == 0)]) 
            #astImages.saveFITS("test_%s.fits" % (tileName), data, wcs)
    
    # Optional subtraction of point sources based on a catalog
    # We'll also (optionally) add a small mask at these locations to the survey mask
    if 'subtractPointSourcesFromCatalog' in list(mapDict.keys()) and mapDict['subtractPointSourcesFromCatalog'] is not None:
        if type(mapDict['subtractPointSourcesFromCatalog']) is not list:
            mapDict['subtractPointSourcesFromCatalog']=[mapDict['subtractPointSourcesFromCatalog']]
        for tab in mapDict['subtractPointSourcesFromCatalog']:
            if type(tab) != atpy.Table:
                tab=atpy.Table().read(catalogPath)
            tab=catalogs.getCatalogWithinImage(tab, data.shape, wcs) 
            model=makeModelImage(data.shape, wcs, tab, mapDict['beamFileName'], obsFreqGHz = None)
            if model is not None:
                data=data-model
            # Optionally blank small exclusion zone around these sources in survey mask
            # (this is only needed for SZ searches, to avoid any issue with possible oversubtraction)
            if 'maskSubtractedPointSources' in list(mapDict.keys()) and mapDict['maskSubtractedPointSources'] == True: 
                for row in tab:
                    x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
                    if surveyMask[int(y), int(x)] != 0:
                        rArcminMap=np.ones(data.shape, dtype = float)*1e6
                        if row['SNR'] > 1000:
                            maskRadiusArcmin=10.0
                        else:
                            maskRadiusArcmin=4.0
                        rArcminMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(rArcminMap, wcs, 
                                                                                    row['RADeg'], row['decDeg'], 
                                                                                    maskRadiusArcmin/60)
                        rArcminMap=rArcminMap*60
                        surveyMask[rArcminMap < maskRadiusArcmin]=0
    
    # Add the map data to the dict
    mapDict['data']=data
    mapDict['weights']=weights
    mapDict['wcs']=wcs
    mapDict['surveyMask']=surveyMask
    mapDict['psMask']=psMask
    mapDict['tileName']=tileName
    
    # Sanity check - no point continuing if masks are different shape to map (easier to tell user here)
    if mapDict['data'].shape != mapDict['psMask'].shape:
        raise Exception("Map and point source mask dimensions are not the same (they should also have same WCS)")
    if mapDict['data'].shape != mapDict['surveyMask'].shape:
        raise Exception("Map and survey mask dimensions are not the same (they should also have same WCS)")
    
    ## Save trimmed weights - this isn't necessary
    #if os.path.exists(diagnosticsDir+os.path.sep+"weights#%s.fits" % (tileName)) == False:
        #astImages.saveFITS(diagnosticsDir+os.path.sep+"weights#%s.fits" % (tileName), weights, wcs)
    
    return mapDict

#------------------------------------------------------------------------------------------------------------
def simCMBMap(shape, wcs, noiseLevel = 0.0, beamFileName = None, seed = None):
    """Generate a simulated CMB map, optionally convolved with the beam and with (white) noise added.
    
    Args:
        shape (:obj:`tuple`): A tuple describing the map (numpy array) shape in pixels (height, width).
        wcs (:obj:`astWCS.WCS`): An astWCS object.
        noiseLevel (:obj:`numpy.ndarray` or float): If a single number, this is taken as sigma (in map units,
            usually uK) for generating white noise that is added across the whole map. Alternatively, an array
            with the same dimensions as shape may be used, specifying sigma (in map units) per corresponding 
            pixel. Noise will only be added where non-zero values appear in noiseLevel.
        beamFileName (:obj:`str`): The file name of the text file that describes the beam with which the map will be
            convolved. If None, no beam convolution is applied.
        seed (:obj:`int`): The seed used for the random CMB realisation.
            
    Returns:
        A map (:obj:`numpy.ndarray`)
    
    """
    
    from pixell import enmap, utils, powspec
    import astropy.wcs as apywcs
    enlibWCS=apywcs.WCS(wcs.header)
    ps=powspec.read_spectrum(nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"planck_lensedCls.dat", 
                             scale = True)
    randMap=enmap.rand_map(shape, enlibWCS, ps, seed = seed)
    np.random.seed()    # Otherwise, we will end up with identical white noise...
    
    if beamFileName != None:
        randMap=convolveMapWithBeam(randMap, wcs, beamFileName)

    if type(noiseLevel) == np.ndarray:
        mask=np.nonzero(noiseLevel)
        generatedNoise=np.zeros(randMap.shape)
        generatedNoise[mask]=np.random.normal(0, noiseLevel[mask], noiseLevel[mask].shape)
        randMap=randMap+generatedNoise
    else:
        if noiseLevel > 0:
            generatedNoise=np.random.normal(0, noiseLevel, randMap.shape)
            randMap=randMap+generatedNoise

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
    degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, 
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
def getPixelAreaArcmin2Map(mapData, wcs):
    """Returns a map of pixel area in arcmin2
    
    """
    
    # Get pixel size as function of position
    pixAreasDeg2=[]
    RACentre, decCentre=wcs.getCentreWCSCoords()
    x0, y0=wcs.wcs2pix(RACentre, decCentre)
    x1=x0+1
    for y0 in range(mapData.shape[0]):
        y1=y0+1
        ra0, dec0=wcs.pix2wcs(x0, y0)
        ra1, dec1=wcs.pix2wcs(x1, y1)
        xPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
        yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
        pixAreasDeg2.append(xPixScale*yPixScale)
    pixAreasDeg2=np.array(pixAreasDeg2)
    pixAreasArcmin2=pixAreasDeg2*(60**2)
    pixAreasArcmin2Map=np.array([pixAreasArcmin2]*mapData.shape[1]).transpose()
    
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
                                                         copyFilters = True)
        
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
def makeModelImage(shape, wcs, catalog, beamFileName, obsFreqGHz = None, GNFWParams = 'default', 
                   cosmoModel = None, applyPixelWindow = True, override = None,
                   validAreaSection = None):
    """Make a map with the given dimensions (shape) and WCS, containing model clusters or point sources, 
    with properties as listed in the catalog. This can be used to either inject or subtract sources
    from real maps.
    
    Args:
        shape (tuple): The dimensions of the output map (height, width) that will contain the model sources.
        wcs (:obj:`astWCS.WCS`): A WCS object that defines the coordinate system of the map. 
        catalog (:obj:`astropy.table.Table`): An astropy Table object containing the catalog. This must 
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
        override (dict, optional): Used only by cluster catalogs. If a dictionary containing keys
            {'M500', 'redshift'} is given, all objects in the model image are forced to have the 
            corresponding angular size. Used by :meth:`sourceInjectionTest`.
        applyPixelWindow (bool, optional): If True, apply the pixel window function to the map.
        validAreaSection (list, optional): Pixel coordinates within the wcs in the format
            [xMin, xMax, yMin, yMax] that define valid area within the model map. Pixels outside this 
            region will be set to zero. Use this to remove overlaps between tile boundaries.
            
    Returns:
        Map containing injected sources, or None if there are no objects within the map dimensions.
    
    """
    
    modelMap=np.zeros(shape, dtype = float)
    
    # This works per-tile, so throw out objects that aren't in it
    catalog=catalogs.getCatalogWithinImage(catalog, shape, wcs)
    if len(catalog) == 0:
        return None

    if cosmoModel is None:
        cosmoModel=signals.FlatLambdaCDM(H0 = 70.0, Om0 = 0.3, Ob0 = 0.05, Tcmb0 = signals.TCMB)
    
    # Set initial max size in degrees from beam file (used for sources; clusters adjusted for each object)
    numFWHM=5.0
    beam=signals.BeamProfile(beamFileName = beamFileName)
    maxSizeDeg=beam.rDeg[np.argmin(abs(beam.profile1d-0.5))]*2*numFWHM 
    
    # Map of distance(s) from objects - this will get updated in place (fast)
    degreesMap=np.ones(modelMap.shape, dtype = float)*1e6
    
    if 'fixed_y_c' in catalog.keys():
        # Clusters: for speed - assume all objects are the same shape
        if override is not None:
            fluxScaleMap=np.zeros(modelMap.shape)
            for row in catalog:
                degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, 
                                                                            row['RADeg'], row['decDeg'], 
                                                                            maxSizeDeg)
                fluxScaleMap[yBounds[0]:yBounds[1], xBounds[0]:xBounds[1]]=row['fixed_y_c']*1e-4
            theta500Arcmin=signals.calcTheta500Arcmin(override['redshift'], override['M500'], cosmoModel)
            maxSizeDeg=5*(theta500Arcmin/60)
            modelMap=signals.makeArnaudModelSignalMap(override['redshift'], override['M500'], degreesMap, 
                                                      wcs, beam, GNFWParams = GNFWParams,
                                                      maxSizeDeg = maxSizeDeg, convolveWithBeam = False)
            modelMap=modelMap*fluxScaleMap
            modelMap=convolveMapWithBeam(modelMap, wcs, beam, maxDistDegrees = 1.0)

        # Clusters - insert one at a time (with different scales etc.) - currently taking ~1.6 sec per object
        else:
            count=0
            for row in catalog:
                # This should avoid overlaps if tiled - we only add cluster if inside areaMask region
                # NOTE: Should move this out of this switch so applied to all catalog types
                if validAreaSection is not None:
                    x0, x1, y0, y1=validAreaSection
                    x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
                    if (x >= x0 and x < x1 and y >= y0 and y < y1) == False:
                        continue
                count=count+1
                # NOTE: We need to think about this a bit more, for when we're not working at fixed filter scale
                if 'true_M500' in catalog.keys():
                    M500=row['true_M500']*1e14
                    z=row['redshift']
                    y0ToInsert=row['fixed_y_c']*1e-4
                else:
                    if 'template' not in catalog.keys():
                        raise Exception("No M500, z, or template column found in catalog.")
                    bits=row['template'].split("#")[0].split("_")
                    M500=float(bits[1][1:].replace("p", "."))
                    z=float(bits[2][1:].replace("p", "."))
                    y0ToInsert=row['y_c']*1e-4  # or fixed_y_c...
                theta500Arcmin=signals.calcTheta500Arcmin(z, M500, cosmoModel)
                maxSizeDeg=5*(theta500Arcmin/60)
                degreesMap=np.ones(modelMap.shape, dtype = float)*1e6 # NOTE: never move this
                degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, 
                                                                            row['RADeg'], row['decDeg'], 
                                                                            maxSizeDeg)
                signalMap=signals.makeArnaudModelSignalMap(z, M500, degreesMap, wcs, beam, 
                                                           GNFWParams = GNFWParams, amplitude = y0ToInsert,
                                                           maxSizeDeg = maxSizeDeg, convolveWithBeam = False)
                modelMap=modelMap+signalMap
            modelMap=convolveMapWithBeam(modelMap, wcs, beam, maxDistDegrees = 1.0)

    else:
        # Sources - note this is extremely fast, but will be wrong for sources close enough to blend
        fluxScaleMap=np.zeros(modelMap.shape)
        for row in catalog:
            degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, 
                                                                        row['RADeg'], row['decDeg'], 
                                                                        maxSizeDeg)
            fluxScaleMap[yBounds[0]:yBounds[1], xBounds[0]:xBounds[1]]=row['deltaT_c']
        modelMap=signals.makeBeamModelSignalMap(degreesMap, wcs, beam)
        modelMap=modelMap*fluxScaleMap

    # Optional: convert map to deltaT uK
    # This should only be used if working with clusters - source amplitudes are fed in as delta T uK already
    if obsFreqGHz is not None:
        modelMap=convertToDeltaT(modelMap, obsFrequencyGHz = obsFreqGHz)
    
    # Optional: apply pixel window function - generally this should be True
    # (because the source-insertion routines in signals.py interpolate onto the grid rather than average)
    if applyPixelWindow == True:
        modelMap=enmap.apply_window(modelMap, pow = 1.0)
        
    return modelMap
        
#------------------------------------------------------------------------------------------------------------
def sourceInjectionTest(config, writeRankTable = True):
    """Insert sources with known positions and properties into the map, apply the filter, and record their
    offset with respect to the true location as a function of S/N (for the fixed reference scale only).
    If the inserted sources are clusters, the Q function will be applied to the output fluxes, to account 
    for any mismatch between the reference filter scale and the inserted clusters.
    
    Writes output to the diagnostics/ directory.
    
    Args:
        config (:obj:`nemo.startUp.NemoConfig`): Nemo configuration object.
        writeRankTable (bool, optional): If True, saves a table as output for this MPI rank under the 
            diagnostics/ directory. Useful for MPI debugging only.
    
    Returns:
        An astropy Table containing recovered position offsets and fluxes versus fixed_SNR for inserted
        sources.
    
    """
    
    simRootOutDir=config.diagnosticsDir+os.path.sep+"sourceInjection_rank%d" % (config.rank)
    SNRKeys=['fixed_SNR']
    
    # We don't copy this, because it's complicated due to containing MPI-related things (comm)
    # So... we modify the config parameters in-place, and restore them before exiting this method
    simConfig=config
    
    # This should make it quicker to generate test catalogs (especially when using tiles)
    selFn=completeness.SelFn(config.selFnDir, 4.0, configFileName = config.configFileName,
                             enableCompletenessCalc = False, setUpAreaMask = True,
                             tileNames = config.tileNames)
    
    print(">>> Position recovery test [rank = %d] ..." % (config.rank))

    if 'sourceInjectionIterations' not in config.parDict.keys():
        numIterations=1
    else:
        numIterations=config.parDict['sourceInjectionIterations']
    
    # For clusters, we may want to run multiple scales
    # We're using theta500Arcmin as the label here
    filtDict=simConfig.parDict['mapFilters'][0]
    if filtDict['class'].find("ArnaudModel") != -1:
        clusterMode=True
        if 'sourceInjectionModels' not in config.parDict.keys():
            sourceInjectionModelList=[{'redshift': 0.4, 'M500': 2e14}]
        else:
            sourceInjectionModelList=config.parDict['sourceInjectionModels']
        for sourceInjectionModel in sourceInjectionModelList:
            label='%.2f' % (signals.calcTheta500Arcmin(sourceInjectionModel['redshift'], 
                                                       sourceInjectionModel['M500'], signals.fiducialCosmoModel))
            sourceInjectionModel['label']=label
        # We need Q for flux recovery stuff...
        tckQFitDict=signals.loadQ(config.selFnDir+os.path.sep+"QFit.fits", tileNames = config.tileNames)
    else:
        # Sources
        clusterMode=False
        sourceInjectionModelList=[{'label': 'pointSource'}]
    #
    if 'sourceInjectionSourcesPerTile' not in config.parDict.keys():
        numSourcesPerTile=300
    else:
        numSourcesPerTile=config.parDict['sourceInjectionSourcesPerTile']
    
    # We need the actual catalog to throw out spurious 'recoveries'
    # i.e., we only want to cross-match with objects we injected
    catFileName=config.rootOutDir+os.path.sep+"%s_optimalCatalog.fits" % (os.path.split(config.rootOutDir)[-1])
    if os.path.exists(catFileName) == False:
        raise Exception("Catalog file '%s' not found - needed to do source injection test." % (catFileName))
    realCatalog=atpy.Table().read(catFileName)
    
    # Run each scale / model and then collect everything into one table afterwards
    # NOTE: raw flux error in catalogs is from RMS map, but e.g. outFlux will have Q applied here if cluster
    SNRDict={}
    rArcminDict={}
    inFluxDict={}
    outFluxDict={}
    noiseLevelDict={}
    tileNamesDict={}
    for sourceInjectionModel in sourceInjectionModelList:
        SNRDict[sourceInjectionModel['label']]=[]
        rArcminDict[sourceInjectionModel['label']]=[]
        inFluxDict[sourceInjectionModel['label']]=[]
        outFluxDict[sourceInjectionModel['label']]=[]
        noiseLevelDict[sourceInjectionModel['label']]=[]
        tileNamesDict[sourceInjectionModel['label']]=[]
        for i in range(numIterations):        
            print(">>> Source injection and recovery test %d/%d [rank = %d] ..." % (i+1, numIterations, config.rank))
                        
            # NOTE: This block below should be handled when parsing the config file - fix/remove
            # Optional override of default GNFW parameters (used by Arnaud model), if used in filters given
            if 'GNFWParams' not in list(simConfig.parDict.keys()):
                simConfig.parDict['GNFWParams']='default'
            for filtDict in simConfig.parDict['mapFilters']:
                filtDict['params']['GNFWParams']=simConfig.parDict['GNFWParams']
            
            # We don't want to save/cache position recovery test maps
            for filtDict in simConfig.parDict['mapFilters']:
                keysToFalsify=['saveFilteredMaps', 'savePlots']
                for key in keysToFalsify:
                    filtDict['params'][key]=False
            
            # Delete all non-reference scale filters (otherwise we'd want to cache all filters for speed)
            # NOTE: As it stands, point-source only runs may not define photFilter - we need to handle that
            # That should be obvious, as mapFilters will only have one entry
            for filtDict in simConfig.parDict['mapFilters']:
                if filtDict['label'] == simConfig.parDict['photFilter']:
                    break
            simConfig.parDict['mapFilters']=[filtDict] 
                
            # Filling maps with injected sources will be done when maps.preprocessMapDict is called by the filter object
            # So, we only generate the catalog here
            print("... generating mock catalog ...")
            if filtDict['class'].find("ArnaudModel") != -1:
                fluxCol='fixed_y_c'
                noiseLevelCol='fixed_err_y_c'
                # Quick test catalog - takes < 1 sec to generate
                mockCatalog=catalogs.generateTestCatalog(config, numSourcesPerTile, 
                                                         amplitudeColumnName = 'fixed_y_c', 
                                                         amplitudeRange = [0.001, 1], 
                                                         amplitudeDistribution = 'linear',
                                                         selFn = selFn)
                # Or... proper mock, but this takes ~24 sec for E-D56
                #mockCatalog=pipelines.makeMockClusterCatalog(config, writeCatalogs = False, verbose = False)[0]                
                injectSources={'catalog': mockCatalog, 'GNFWParams': config.parDict['GNFWParams'], 
                               'override': sourceInjectionModel}
            elif filtDict['class'].find("BeamModel") != -1:
                fluxCol='deltaT_c'
                noiseLevelCol='err_deltaT_c'
                raise Exception("Haven't implemented generating mock source catalogs here yet")
            else:
                raise Exception("Don't know how to generate injected source catalogs for filterClass '%s'" % (filtDict['class']))
            for mapDict in simConfig.unfilteredMapsDictList:
                mapDict['injectSources']=injectSources
            
            # NOTE: we need to zap ONLY specific maps for when we are running in parallel
            for tileName in simConfig.tileNames:
                mapFileNames=glob.glob(simRootOutDir+os.path.sep+"filteredMaps"+os.path.sep+"*#%s_*.fits" % (tileName))
                for m in mapFileNames:
                    os.remove(m)

            # Ideally we shouldn't have blank tiles... but if we do, skip
            if len(mockCatalog) > 0:
                recCatalog=pipelines.filterMapsAndMakeCatalogs(simConfig, rootOutDir = simRootOutDir,
                                                               copyFilters = True, useCachedMaps = False)
                # We should be conservative in removing potential matches with real objects
                # Because we have a huge sky area and there's no reason to risk contamination of this kind
                # Effectively this is the same as using 5' circular holes in the survey mask on real objects
                # (but actually adding the avoidance radius parameter to the test catalogs really solved this)
                if len(recCatalog) > 0:
                    recCatalog=catalogs.removeCrossMatched(recCatalog, realCatalog, radiusArcmin = 5.0)
                if len(recCatalog) > 0:
                    try:
                        x_mockCatalog, x_recCatalog, rDeg=catalogs.crossMatch(mockCatalog, recCatalog, radiusArcmin = 5.0)
                    except:
                        raise Exception("Position recovery test: cross match failed on tileNames = %s; mockCatalog length = %d; recCatalog length = %d" % (str(simConfig.tileNames), len(mockCatalog), len(recCatalog)))
                    # If we're using clusters, we need to put in the Q correction
                    # NOTE: This assumes the model name gives theta500c in arcmin!
                    if clusterMode == True:
                        for tileName in np.unique(x_recCatalog['tileName']):
                            theta500Arcmin=float(sourceInjectionModel['label'])
                            Q=interpolate.splev(theta500Arcmin, tckQFitDict[tileName])
                            mask=(x_recCatalog['tileName'] == tileName)
                            x_recCatalog[fluxCol][mask]=x_recCatalog[fluxCol][mask]/Q
                    # Store everything - analyse later
                    SNRDict[sourceInjectionModel['label']]=SNRDict[sourceInjectionModel['label']]+x_recCatalog['fixed_SNR'].tolist()
                    rArcminDict[sourceInjectionModel['label']]=rArcminDict[sourceInjectionModel['label']]+(rDeg*60).tolist()
                    inFluxDict[sourceInjectionModel['label']]=inFluxDict[sourceInjectionModel['label']]+x_mockCatalog[fluxCol].tolist()
                    outFluxDict[sourceInjectionModel['label']]=outFluxDict[sourceInjectionModel['label']]+x_recCatalog[fluxCol].tolist()
                    noiseLevelDict[sourceInjectionModel['label']]=noiseLevelDict[sourceInjectionModel['label']]+x_recCatalog[noiseLevelCol].tolist()
                    tileNamesDict[sourceInjectionModel['label']]=tileNamesDict[sourceInjectionModel['label']]+x_recCatalog['tileName'].tolist()

        SNRDict[sourceInjectionModel['label']]=np.array(SNRDict[sourceInjectionModel['label']])
        rArcminDict[sourceInjectionModel['label']]=np.array(rArcminDict[sourceInjectionModel['label']])
        inFluxDict[sourceInjectionModel['label']]=np.array(inFluxDict[sourceInjectionModel['label']])
        outFluxDict[sourceInjectionModel['label']]=np.array(outFluxDict[sourceInjectionModel['label']])
        noiseLevelDict[sourceInjectionModel['label']]=np.array(noiseLevelDict[sourceInjectionModel['label']])
        tileNamesDict[sourceInjectionModel['label']]=np.array(tileNamesDict[sourceInjectionModel['label']])
        
    # Just collect results as long tables (model, SNR, rArcmin, inFlux, outFlux) that we can later stack and average etc.
    # (see positionRecoveryAnalysis below)
    models=[]
    SNRs=[]
    rArcmin=[]
    inFlux=[]
    outFlux=[]
    noiseLevel=[]
    tileNames=[]
    for sourceInjectionModel in sourceInjectionModelList:
        label=sourceInjectionModel['label']
        models=models+[label]*len(SNRDict[label])
        SNRs=SNRs+SNRDict[label].tolist()
        rArcmin=rArcmin+rArcminDict[label].tolist()
        inFlux=inFlux+inFluxDict[label].tolist()
        outFlux=outFlux+outFluxDict[label].tolist()
        noiseLevel=noiseLevel+noiseLevelDict[label].tolist()
        tileNames=tileNames+tileNamesDict[label].tolist()
    resultsTable=atpy.Table()
    resultsTable.add_column(atpy.Column(models, 'sourceInjectionModel'))
    resultsTable.add_column(atpy.Column(SNRs, 'fixed_SNR'))
    resultsTable.add_column(atpy.Column(rArcmin, 'rArcmin'))
    resultsTable.add_column(atpy.Column(inFlux, 'inFlux'))
    resultsTable.add_column(atpy.Column(outFlux, 'outFlux'))
    resultsTable.add_column(atpy.Column(noiseLevel, 'noiseLevel'))
    resultsTable.add_column(atpy.Column(tileNames, 'tileName'))

    # Shouldn't be necessary BUT seems we have trouble gathering very large runs
    # So it's actually more reliable to write/read from disk
    if writeRankTable == True:
        fitsOutFileName=config.diagnosticsDir+os.path.sep+"sourceInjection_rank%d.fits" % (config.rank)
        resultsTable.write(fitsOutFileName, overwrite = True)
    
    # Restore the original config parameters (which we overrode here)
    config.restoreConfig()

    return resultsTable

#------------------------------------------------------------------------------------------------------------
def positionRecoveryAnalysis(posRecTable, plotFileName, percentiles = [50, 95, 99.7],
                             plotRawData = True, rawDataAlpha =  0.01, pickleFileName = None,
                             selFnDir = None):
    """Estimate and plot position recovery accuracy as function of fixed filter scale S/N (fixed_SNR), using
    the contents of posRecTable (see positionRecoveryTest).
    
    Args:
        posRecTable (:obj:`astropy.table.Table`): Table containing recovered position offsets versus fixed_SNR
            for various cluster/source models (produced by sourceInjectionTest).
        plotFileName (str): Path where the plot file will be written.
        percentiles (list, optional): List of percentiles to plot (some interpolation will be done) and
            for which corresponding model fit parameters will be saved (if selFnDir is not None).
        plotRawData (bool, optional): Plot the raw (fixed_SNR, positional offset) data in the background.
        pickleFileName (string, optional): Saves the percentile contours data as a pickle file if not None.
            This is saved a dictionary with top-level keys named according to percentilesToPlot.
        selFnDir (string, optional): If given, model fit parameters will be written to a file named
            posRecModelParameters.txt under the given selFn directory path.
            
    """

    # Evaluate %-age of sample in bins of SNR within some rArcmin threshold
    # No longer separating by input model (clusters are all shapes anyway)
    tab=posRecTable
    SNREdges=np.linspace(3.0, 10.0, 36)#np.linspace(0, 10, 101)
    SNRCentres=(SNREdges[1:]+SNREdges[:-1])/2.
    rArcminThreshold=np.linspace(0, 10, 101)
    grid=np.zeros([rArcminThreshold.shape[0], SNREdges.shape[0]-1])
    totalGrid=np.zeros(grid.shape)
    withinRGrid=np.zeros(grid.shape)
    for i in range(SNREdges.shape[0]-1):
        SNRMask=np.logical_and(tab['fixed_SNR'] >= SNREdges[i], tab['fixed_SNR'] < SNREdges[i+1])
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
        plt.plot(posRecTable['fixed_SNR'], posRecTable['rArcmin'], '.', color = '#A0A0A0', alpha = rawDataAlpha)
    contoursDict={}
    for i in range(len(levelsList)):
        vertices=contours.collections[i].get_paths()[0].vertices
        SNRs=vertices[:, 0]
        rArcminAtProb=vertices[:, 1]
        labelStr="%.1f" % (percentiles[i]) + "%"
        contoursDict[labelStr]={'fixed_SNR': SNRs, 'rArcmin': rArcminAtProb}
        plt.plot(SNRs, rArcminAtProb, label = labelStr, lw = 3)
    plt.xlim(minSNR, maxSNR)
    #plt.ylim(0, 5)
    plt.ylim(0,3)
    plt.legend(loc = 'upper right')
    plt.xlabel("SNR$_{2.4}$")
    plt.ylabel("Recovered Position Offset ($^\prime$)")
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
            plt.plot(posRecTable['fixed_SNR'], posRecTable['rArcmin'], '.', color = '#A0A0A0', alpha = rawDataAlpha)
        for key in keys:
            a=contoursDict[key]
            valid=np.where(a['fixed_SNR'] >= 4.1)
            snr=a['fixed_SNR'][valid]
            rArcmin=a['rArcmin'][valid]
            try:
                results=optimize.curve_fit(catalogs._posRecFitFunc, snr, rArcmin)
            except:
                print("... WARNING: curve_fit failed for key = %s ..." % (key))
                continue
            bestFitSNRFold, bestFitPedestal, bestFitNorm=results[0]
            fitParamsDict[key]=np.array([bestFitSNRFold, bestFitPedestal, bestFitNorm])
            fitSNRs=np.linspace(4, 10, 100)
            plt.plot(fitSNRs, catalogs._posRecFitFunc(fitSNRs, bestFitSNRFold, bestFitPedestal, bestFitNorm), '-', label = key)
        plt.ylim(0, 3)
        plt.legend(loc = 'upper right')
        plt.xlim(snr.min(), snr.max())
        plt.xlabel("SNR$_{2.4}$")
        plt.ylabel("Recovered Position Offset ($^\prime$)")
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
def saveFITS(outputFileName, mapData, wcs, compressed = False):
    """Writes a map (2d image array) to a new .fits file.
    
    Args:
        outputFileName (str): Filename of output FITS image.
        mapData (:obj:`np.ndarray`): Map data array.
        wcs (:obj:`astWCS.WCS`): Map WCS object.
        compressed (bool): If True, writes a compressed image.
    
    """
    
    if os.path.exists(outputFileName):
        os.remove(outputFileName)
    
    if compressed == False:
        if wcs is not None:
            hdu=pyfits.PrimaryHDU(mapData, wcs.header)
        else:
            hdu=pyfits.PrimaryHDU(mapData, None)
    
    if compressed == True:
        if wcs is not None:
            hdu=pyfits.CompImageHDU(np.array(mapData, dtype = float), wcs.header)
        else:
            hdu=pyfits.CompImageHDU(np.array(mapData, dtype = float), None)
            
    newImg=pyfits.HDUList()
    newImg.append(hdu)
    newImg.writeto(outputFileName)
    newImg.close()
