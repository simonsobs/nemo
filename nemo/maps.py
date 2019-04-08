"""

This module contains tools for manipulating maps (e.g., conversion of units etc.).

"""

from astLib import *
from scipy import ndimage
from scipy import interpolate
from scipy.signal import convolve as scipy_convolve
import astropy.io.fits as pyfits
import astropy.table as atpy
import astropy.stats as apyStats
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
import IPython
import nemo
from . import catalogs
from . import signals
from . import photometry
from . import plotSettings
from . import pipelines
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
def makeTileDeck(parDict):
    """Makes a tileDeck multi-extension .fits file, if the needed parameters are given in parDict, or
    will handle setting up such a file if given directly in unfilteredMapsDictList in parDict (and the .par
    file). Adjusts unfilteredMapsDictList accordingly and returns it.
    
    If the options for making a tileDeck image aren't given in parDict, then we pass through a standard
    single extension file (or rather the path to it, as originally given)
    
    NOTE: If the map given in unfilteredMaps is 3d (enki gives I, Q, U as a datacube), then this will extract
    only the I (temperature) part and save that in the tileDeck file. This will need changing if hunting for
    polarized sources...
    
    Returns unfilteredMapsDictList [input for filterMaps], list of extension names
    
    """
    
    if 'makeTileDeck' not in list(parDict.keys()):
        parDict['makeTileDeck']=False

    # Some of this is rather clunky...
    unfilteredMapsDictList=[]
    if parDict['makeTileDeck'] == False:
        extNames=[]        
        for mapDict in parDict['unfilteredMaps']:
            unfilteredMapsDictList.append(mapDict.copy())
            img=pyfits.open(mapDict['mapFileName'])
            if extNames == []:
                for ext in img:
                    extNames.append(ext.name)
            else:
                for ext in img:
                    if ext.name not in extNames:
                        raise Exception("extension names do not match between all maps in unfilteredMapsDictList")
            img.close()
    else:
        extNames=[]        
        for mapDict in parDict['unfilteredMaps']:
                        
            # Added an option to define tiles in the .par file... otherwise, we will do the automatic tiling
            if 'tileDefinitions' in list(parDict.keys()):
                if 'tileDefLabel' in list(parDict.keys()):
                    tileDefLabel=parDict['tileDefLabel']
                else:
                    tileDefLabel='userDefined'
                tileDeckFileNameLabel="%s_%.1f" % (tileDefLabel, parDict['tileOverlapDeg'])
                defineTilesAutomatically=False
            else:
                tileDeckFileNameLabel="%dx%d_%.1f" % (parDict['numHorizontalTiles'],
                                                      parDict['numVerticalTiles'], 
                                                      parDict['tileOverlapDeg'])
                defineTilesAutomatically=True
            
            # Figure out what the input / output files will be called
            # NOTE: we always need to make a survey mask if none exists, as used to zap over regions, so that gets special treatment
            fileNameKeys=['mapFileName', 'weightsFileName', 'pointSourceMask', 'surveyMask']
            inFileNames=[]
            outFileNames=[]
            mapTypeList=[]
            for f in fileNameKeys:
                if f in list(mapDict.keys()) and mapDict[f] != None:
                    inFileNames.append(mapDict[f])
                    mapDir, mapFileName=os.path.split(mapDict[f])
                    if mapDir != '':
                        mapDirStr=mapDir+os.path.sep
                    else:
                        mapDirStr=''
                    outFileNames.append(mapDirStr+"tileDeck_%s_" % (tileDeckFileNameLabel)+mapFileName)
                    mapTypeList.append(f)
                if f == 'surveyMask' and mapDict[f] == None:
                    inFileNames.append(None)
                    outFileNames.append(outFileNames[0].replace(".fits", "_surveyMask.fits"))
                    mapTypeList.append(f)

            allFilesMade=True
            for f in outFileNames:
                if os.path.exists(f) == False:
                    allFilesMade=False
            
            if allFilesMade == True:
                # We need the extension names only here...
                img=pyfits.open(outFileNames[0])
                if extNames == []:
                    for ext in img:
                        extNames.append(ext.name)
                else:
                    for ext in img:
                        if ext.name not in extNames:
                            raise Exception("extension names do not match between all maps in unfilteredMapsDictList")
            else:
                
                # Whether we make tiles automatically or not, we need the WCS from somewhere...
                if 'surveyMask' in list(mapDict.keys()) and mapDict['surveyMask'] != None:
                    wht=pyfits.open(mapDict['surveyMask'])
                    print(">>> Using survey mask to determine tiling ...")
                else:
                    wht=pyfits.open(mapDict['weightsFileName'])
                    print(">>> Using weight map to determine tiling ...")
                wcs=astWCS.WCS(wht[0].header, mode = 'pyfits')
                tileOverlapDeg=parDict['tileOverlapDeg']
   
                if defineTilesAutomatically == True:
                    
                    if 'surveyMask' in mapDict.keys() and mapDict['surveyMask'] == None:
                        print("... WARNING: same tiling not guaranteed across multiple frequencies ...")
                    
                    # NOTE: here we look at surveyMask first to determine where to put down tiles
                    # since this will ensure algorithm uses same tiles for multi-freq data
                    # Otherwise, we use the wht image (but then can't guarantee f090 and f150 have same tiles)
                    deckWht=pyfits.HDUList()
                    whtData=wht[0].data
                    mapWidth=whtData.shape[1]
                    mapHeight=whtData.shape[0]

                    # Figure out where edges are
                    edges={}
                    for y in range(mapHeight):
                        xIndices=np.where(whtData[y] != 0)[0]
                        if len(xIndices) > 0:
                            xMin=xIndices.min()
                            xMax=xIndices.max()
                            edges[y]=[xMin, xMax]
                    
                    # Starting from bottom left, work our way around the map adding tiles, ignoring blank regions
                    numHorizontalTiles=parDict['numHorizontalTiles']
                    numVerticalTiles=parDict['numVerticalTiles']
                    ys=list(edges.keys())
                    ys.sort()
                    ys=np.array(ys)
                    coordsList=[]
                    extNames=[]
                    tileHeightPix=int(np.ceil((ys.max()-ys.min())/float(numVerticalTiles)))
                    for i in range(numVerticalTiles):
                        yMin=ys.min()+i*tileHeightPix
                        yMax=ys.min()+(i+1)*tileHeightPix
                        keys=np.arange(yMin, yMax)
                        minXMin=1e6
                        maxXMax=0
                        for k in keys:
                            if k in list(edges.keys()):
                                xMin, xMax=edges[k]
                                if xMin < minXMin:
                                    minXMin=xMin
                                if xMax > maxXMax:
                                    maxXMax=xMax
                        tileWidthPix=int(np.ceil(maxXMax-minXMin)/float(numHorizontalTiles))
                        for j in range(numHorizontalTiles):
                            xMin=minXMin+j*tileWidthPix
                            xMax=minXMin+(j+1)*tileWidthPix
                            coordsList.append([xMin, xMax, yMin, yMax])
                            extNames.append("%d_%d" % (j, i))
                    
                    # Not sure if this will actually tidy up...
                    wht.close()
                    del whtData
                
                else:
                    # Use user-defined tiles - this is a bit of a faff, to avoid re-writing below bit where we make the tiles...
                    extNames=[]
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
                        extNames.append(tileDict['extName'])   
                
                # Output a .reg file for debugging (pixel coords)
                outFile=open(outFileNames[0].replace(".fits", "_tiles.reg"), "w")
                outFile.write("# Region file format: DS9 version 4.1\n")
                outFile.write('global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                outFile.write("image\n")
                for c, name in zip(coordsList, extNames):
                    outFile.write('polygon(%d, %d, %d, %d, %d, %d, %d, %d) # text="%s"\n' % (c[0], c[2], c[0], c[3], c[1], c[3], c[1], c[2], name))
                outFile.close()
                #print("check tiles")
                #sys.exit()
                
                # Make tiles
                # NOTE: we accommodate having user-defined regions for calculating noise power in filters here
                # Since we would only use such an option with tileDeck files, this should be okay
                # Although since we do this by modifying headers, would need to remake tileDeck files each time adjusted in .par file
                # NOTE: now treating surveyMask as special, and zapping overlap regions there (simplify selection function stuff later)
                for mapType, inMapFileName, outMapFileName in zip(mapTypeList, inFileNames, outFileNames):
                    if os.path.exists(outMapFileName) == False:
                        print(">>> Writing tileDeck file %s ..." % (outMapFileName))
                        deckImg=pyfits.HDUList()
                        # Special handling for case where surveyMask = None in the .par file (tidy later...)
                        if mapType == 'surveyMask' and inMapFileName == None:
                            img=pyfits.open(inFileNames[0])
                            mapData=np.ones(img[0].data.shape)
                        else:
                            img=pyfits.open(inMapFileName)
                            mapData=img[0].data

                        # Deal with Sigurd's maps which have T, Q, U as one 3d array
                        # If anyone wants to find polarized sources, this will need changing...
                        if mapData.ndim == 3:
                            mapData=mapData[0, :]
                        for c, name in zip(coordsList, extNames):
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
                                clip=astImages.clipUsingRADecCoords(mapData, wcs, ra1, ra0, dec0, dec1+ddec*count)
                                count=count+1
                            # Old
                            #clip=astImages.clipUsingRADecCoords(mapData, wcs, ra1, ra0, dec0, dec1)
                            print("... adding %s [%d, %d, %d, %d ; %d, %d] ..." % (name, ra1, ra0, dec0, dec1, ra0-ra1, dec1-dec0))
                            header=clip['wcs'].header.copy()
                            if 'tileNoiseRegions' in list(parDict.keys()):
                                if name in list(parDict['tileNoiseRegions'].keys()):
                                    noiseRAMin, noiseRAMax, noiseDecMin, noiseDecMax=parDict['tileNoiseRegions'][name]
                                else:
                                    if 'autoBorderDeg' in parDict['tileNoiseRegions']:
                                        autoBorderDeg=parDict['tileNoiseRegions']['autoBorderDeg']
                                        for tileDef in parDict['tileDefinitions']:
                                            if tileDef['extName'] == name:
                                                break
                                        noiseRAMin, noiseRAMax, noiseDecMin, noiseDecMax=tileDef['RADecSection']
                                        noiseRAMin=noiseRAMin+autoBorderDeg
                                        noiseRAMax=noiseRAMax-autoBorderDeg
                                        noiseDecMin=noiseDecMin+autoBorderDeg
                                        noiseDecMax=noiseDecMax-autoBorderDeg
                                    else:
                                        raise Exception("No entry in tileNoiseRegions in config file for extName '%s' - either add one, or add 'autoBorderDeg': 0.5 (or similar) to tileNoiseRegions" % (name))
                                print("... adding noise region [%.3f, %.3f, %.3f, %.3f] to header %s ..." % (noiseRAMin, noiseRAMax, noiseDecMin, noiseDecMax, name))
                                header['NRAMIN']=noiseRAMin
                                header['NRAMAX']=noiseRAMax
                                header['NDEMIN']=noiseDecMin
                                header['NDEMAX']=noiseDecMax
                            # Survey mask is special: zap overlap regions outside of tile definitions
                            if mapType == 'surveyMask':
                                ra0, dec0=wcs.pix2wcs(x0, y0)
                                ra1, dec1=wcs.pix2wcs(x1, y1)
                                clip_x0, clip_y0=clip['wcs'].wcs2pix(ra0, dec0)
                                clip_x1, clip_y1=clip['wcs'].wcs2pix(ra1, dec1)
                                clip_x0=int(round(clip_x0))
                                clip_x1=int(round(clip_x1))
                                clip_y0=int(round(clip_y0))
                                clip_y1=int(round(clip_y1))
                                zapMask=np.zeros(clip['data'].shape)
                                zapMask[clip_y0:clip_y1, clip_x0:clip_x1]=1.
                                clip['data']=clip['data']*zapMask
                                #astImages.saveFITS("test.fits", zapMask, clip['wcs'])
                            hdu=pyfits.ImageHDU(data = clip['data'].copy(), header = header, name = name)
                            deckImg.append(hdu)    
                        deckImg.writeto(outMapFileName)
                        deckImg.close()
                                
            # Replace entries in unfilteredMapsDictList in place
            for key, outFileName in zip(mapTypeList, outFileNames):
                mapDict[key]=outFileName
            unfilteredMapsDictList.append(mapDict.copy())
    
    return unfilteredMapsDictList, extNames

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
def stitchTiles(filePattern, outFileName, outWCS, outShape, fluxRescale = 1.0):
    """Fast routine for stitching map tiles back together. Since this uses interpolation, you probably don't 
    want to do analysis on the output - this is just for checking / making plots etc.. This routine sums 
    images as it pastes them into the larger map grid. So, if the zeroed (overlap) borders are not handled,
    correctly, this will be obvious in the output.

    NOTE: This assumes RA in x direction, dec in y direction (and CAR projection).
    
    NOTE: This routine only writes output if there are multiple files that match filePattern (to save needless
    duplicating maps if nemo was not run in tileDeck mode).
    
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
    if len(inFiles) < 2:
        return None # No point stitching together 1 tile (probably means didn't run in tiles mode)
    count=0
    for f in inFiles:
        count=count+1
        #print("... %d/%d ..." % (count, len(inFiles)))
        with pyfits.open(f) as img:
            d=img[0].data
            inWCS=astWCS.WCS(img[0].header, mode = 'pyfits')
        xIn=np.arange(d.shape[1])
        yIn=np.arange(d.shape[0])
        inRACoords=np.array(inWCS.pix2wcs(xIn, [0]*len(xIn)))
        inDecCoords=np.array(inWCS.pix2wcs([0]*len(yIn), yIn))
        inRA=inRACoords[:, 0]
        inDec=inDecCoords[:, 1]
        xOut=np.array(RAToX(inRA), dtype = int)
        yOut=np.array(DecToY(inDec), dtype = int)
        for i in range(len(yOut)):
            outData[yOut[i]][xOut]=outData[yOut[i]][xOut]+d[yIn[i], xIn]
    astImages.saveFITS(outFileName, outData*fluxRescale, outWCS)

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
            rRange=nemoCython.makeDegreesDistanceMap(maskedMapData, wcs, obj['RADeg'], obj['decDeg'], 
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
                # NOTE: This only makes sense to do on an unfiltered map...
                if obj['id'] == 1445:
                    print("Fix oversubtraction... peakValue is pointSource + CMB...")
                    IPython.embed()
                    sys.exit()
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
                rRange=nemoCython.makeDegreesDistanceMap(mapData, mapWCS, RADeg, decDeg, (radiusArcmin*4)/60.0)        
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
def preprocessMapDict(mapDict, extName = 'PRIMARY', diagnosticsDir = None):
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
            
    Args:
        mapDict (:obj:`dict`): A dictionary with the same keys as given in the unfilteredMaps list in the 
            .yml configuration file (i.e., it must contain at least the keys 'mapFileName', 'units', and
            'weightsFileName', and may contain some of the optional keys listed above).
        extName (:obj:`str`): Name of the map tile (extension name) to operate on.
        diagnosticsDir (:obj:`str`): Path to a directory where miscellaneous diagnostic data are written.
    
    Returns:
        A dictionary with keys that point to the map itself ('data'), weights ('weights'), masks 
        ('surveyMask', 'pointSourceMask'), and WCS object ('wcs').
    
    """
            
    with pyfits.open(mapDict['mapFileName'], memmap = True) as img:
        wcs=astWCS.WCS(img[extName].header, mode = 'pyfits')
        data=img[extName].data
    
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
        with pyfits.open(mapDict['weightsFileName'], memmap = True) as wht:
            weights=wht[extName].data
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
        with pyfits.open(mapDict['surveyMask'], memmap = True) as smImg:
            surveyMask=smImg[extName].data
    else:
        surveyMask=np.ones(data.shape)
        surveyMask[weights == 0]=0

    if 'pointSourceMask' in list(mapDict.keys()) and mapDict['pointSourceMask'] is not None:
        with pyfits.open(mapDict['pointSourceMask'], memmap = True) as psImg:
            psMask=psImg[extName].data
    else:
        psMask=np.ones(data.shape)
            
    #print("... opened map %s ..." % (mapDict['mapFileName']))
    
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
        outFileName=diagnosticsDir+os.path.sep+"CMBSim_%d#%s.fits" % (mapDict['obsFreqGHz'], extName) 
        astImages.saveFITS(outFileName, data, wcs)
    
    # For position recovery tests
    if 'injectSources' in list(mapDict.keys()):
        # NOTE: Need to add varying GNFWParams here
        # NOTE: We should also stop loading the beam from disk each time
        data=injectSources(data, wcs, mapDict['injectSources'], mapDict['beamFileName'], 
                           obsFreqGHz = mapDict['obsFreqGHz'], GNFWParams = 'default')
    
    # Optional masking of point sources from external catalog
    # Especially needed if using Fourier-space matched filter (and maps not already point source subtracted)
    if 'maskPointSourcesFromCatalog' in list(mapDict.keys()) and mapDict['maskPointSourcesFromCatalog'] is not None:  
        # This is fast enough if using small tiles and running in parallel...
        # If our masking/filling is effective enough, we may not need to mask so much here...     
        tab=atpy.Table().read(mapDict['maskPointSourcesFromCatalog'])
        xyCoords=wcs.wcs2pix(tab['RADeg'].tolist(), tab['decDeg'].tolist()) 
        xyCoords=np.array(xyCoords, dtype = int)
        mask=[]
        for i in range(len(tab)):
            x, y=xyCoords[i][0], xyCoords[i][1]
            if x >= 0 and x < data.shape[1]-1 and y >= 0 and y < data.shape[0]-1:
                mask.append(True)
            else:
                mask.append(False)
        tab=tab[mask]
        # Variable sized holes: based on inspecting sources by deltaT in f150 maps
        tab.add_column(atpy.Column(np.zeros(len(tab)), 'rArcmin'))
        tab['rArcmin'][tab['deltaT_c'] < 500]=3.0
        tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 500, tab['deltaT_c'] < 1000)]=4.0
        tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 1000, tab['deltaT_c'] < 2000)]=5.0
        tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 2000, tab['deltaT_c'] < 3000)]=5.5
        tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 3000, tab['deltaT_c'] < 10000)]=6.0
        tab['rArcmin'][np.logical_and(tab['deltaT_c'] >= 10000, tab['deltaT_c'] < 50000)]=8.0
        tab['rArcmin'][tab['deltaT_c'] >= 50000]=10.0
        psMask=np.ones(data.shape)
        for row in tab:
            rArcminMap=nemoCython.makeDegreesDistanceMap(psMask, wcs, row['RADeg'], row['decDeg'], row['rArcmin']/60.)*60
            psMask[rArcminMap < row['rArcmin']]=0
        # Fill holes with smoothed map + white noise
        pixRad=(10.0/60.0)/wcs.getPixelSizeDeg()
        bckData=ndimage.median_filter(data, int(pixRad)) # Size chosen for max hole size... slow... but quite good
        if mapDict['weightsType'] =='invVar':
            rms=np.zeros(weights.shape)
            rms[np.nonzero(weights)]=1.0/np.sqrt(weights[np.nonzero(weights)])
        else:
            raise Exception("Not implemented white noise estimate for non-inverse variance weights for masking sources from catalog")
        data[np.where(psMask == 0)]=bckData[np.where(psMask == 0)]+np.random.normal(0, rms[np.where(psMask == 0)]) 
        #astImages.saveFITS("test_%s.fits" % (extName), data, wcs)
    
    # Add the map data to the dict
    mapDict['data']=data
    mapDict['weights']=weights
    mapDict['wcs']=wcs
    mapDict['surveyMask']=surveyMask
    mapDict['psMask']=psMask
    mapDict['extName']=extName
    
    # Sanity check - no point continuing if masks are different shape to map (easier to tell user here)
    if mapDict['data'].shape != mapDict['psMask'].shape:
        raise Exception("Map and point source mask dimensions are not the same (they should also have same WCS)")
    if mapDict['data'].shape != mapDict['surveyMask'].shape:
        raise Exception("Map and survey mask dimensions are not the same (they should also have same WCS)")
    
    ## Save trimmed weights - this isn't necessary
    #if os.path.exists(diagnosticsDir+os.path.sep+"weights#%s.fits" % (extName)) == False:
        #astImages.saveFITS(diagnosticsDir+os.path.sep+"weights#%s.fits" % (extName), weights, wcs)
        
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
def convolveMapWithBeam(data, wcs, beamFileName, maxDistDegrees = 1.0):
    """Convolves map given by data, wcs with the beam.
    
    maxDistDegrees sets the size of the convolution kernel.
    
    Returns map (numpy array)
    
    """
    
    RADeg, decDeg=wcs.getCentreWCSCoords()
    degreesMap=nemoCython.makeDegreesDistanceMap(data, wcs, RADeg, decDeg, maxDistDegrees)

    beamMap, sigDict=signals.makeBeamModelSignalMap(degreesMap, wcs, beamFileName)
    ys, xs=np.where(degreesMap < maxDistDegrees)
    yMin=ys.min()
    yMax=ys.max()+1
    xMin=xs.min()
    xMax=xs.max()+1
    beamMap=beamMap[yMin:yMax, xMin:xMax]
    beamMap=beamMap/beamMap.sum()
    
    return scipy_convolve(data, beamMap, mode = 'same')

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
        for extName in simConfig.extNames:
            mapFileNames=glob.glob(simRootOutDir+os.path.sep+"filteredMaps"+os.path.sep+"*#%s_*.fits" % (extName))
            for m in mapFileNames:
                os.remove(m)
                
        simImageDict=pipelines.filterMapsAndMakeCatalogs(simConfig, 
                                                         rootOutDir = simRootOutDir,
                                                         copyFilters = True)
        
        # Write out the last sim map catalog for debugging
        # NOTE: extName here makes no sense - this should be happening in the pipeline call above
        #optimalCatalogFileName=simRootOutDir+os.path.sep+"CMBSim_optimalCatalog#%s.csv" % (extName)    
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
    # (if we're running in serial, then we'll get a giant file name with full extNames list... fix later)
    extNamesLabel="#"+str(config.extNames).replace("[", "").replace("]", "").replace("'", "").replace(", ", "#")
    for k in list(avContaminTabDict.keys()):
        fitsOutFileName=config.diagnosticsDir+os.path.sep+"%s_contaminationEstimate_%s.fits" % (k, extNamesLabel)
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
        binEdges=np.linspace(binMin, binMax, (binMax-binMin)/binStep+1)
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
def injectSources(data, wcs, catalog, beamFileName, obsFreqGHz = 148.0, GNFWParams = 'default'):
    """Inject sources (clusters or point sources) with properties listed in the catalog into the map defined
    by data, wcs.
    
    Args:
        data (:obj:`numpy.ndarray`): The (2d) pixel data of a map in units of uK (delta T CMB).
        wcs (:obj:`astWCS.WCS`): A WCS object that defines the coordinate system of the map. 
        catalog (:obj:`astropy.table.Table`): An astropy Table object containing the catalog. This must 
            include columns named RADeg, decDeg that give object coordinates. For point sources, the 
            amplitude in uK must be given in a column named deltaT_c. For clusters, M500 (in units of
            10^14 MSun) and z must be given (GNFW profile assumed).
        beamFileName: Path to a text file that describes the beam.
        obsFreqGHz (float, optional): Used only by cluster catalogs - for converting SZ y into delta T uK.
        GNFWParams (str or dict, optional): Used only by cluster catalogs. If 'default', the Arnaud et al. 
            (2010) Universal Pressure Profile is assumed. Otherwise, a dictionary that specifies the profile
            parameters can be given here (see gnfw.py).
            
    Returns:
        Map containing injected sources.
    
    """
    
    # Inspect the catalog - are we dealing with point sources or clusters?
    print("inject sources")
    IPython.embed()
    sys.exit()
    #for row in catalog:
        
    #if 'fixed_y_c' in catalog.keys():
    #signals.makeArnaudModelSignalMap(z, M500, obsFreqGHz, degreesMap, wcs, beamFileName, GNFWParams = 'default',
                             #deltaT0 = None, maxSizeDeg = 15.0, convolveWithBeam = True)
    #signals.makeBeamModelSignalMap(degreesMap, wcs, beamFileName, deltaT0 = None)
    
#------------------------------------------------------------------------------------------------------------
def positionRecoveryTest(config, imageDict):
    """Insert sources with known positions and properties into the map, apply the filter, and record their
    offset with respect to the true location as a function of S/N (for the fixed reference scale only).
    
    Writes output to the diagnostics/ directory.
    
    Args:
        config (:obj:`startUp.NemoConfig`): Nemo configuration object.
        imageDict (:obj:`dict`): A dictionary containing the output filtered maps and catalogs from running 
        on the real data (i.e., the output of pipelines.filterMapsAndMakeCatalogs). This will not be 
        modified.
    
    """
    
    simRootOutDir=config.diagnosticsDir+os.path.sep+"posRec_rank%d" % (config.rank)
    SNRKeys=['fixed_SNR']

    print(">>> Position recovery test [rank = %d] ..." % (config.rank))
    t0=time.time()

    # We don't copy this, because it's complicated due to containing MPI-related things (comm)
    # So... we modify the config parameters in-place, and restore them before exiting this method
    simConfig=config
    
    # NOTE: This block below should be handled when parsing the config file - fix/remove
    # Optional override of default GNFW parameters (used by Arnaud model), if used in filters given
    if 'GNFWParams' not in list(simConfig.parDict.keys()):
        simConfig.parDict['GNFWParams']='default'
    for filtDict in simConfig.parDict['mapFilters']:
        filtDict['params']['GNFWParams']=simConfig.parDict['GNFWParams']
    
    # Delete all non-reference scale filters (otherwise we'd want to cache all filters for speed)
    # NOTE: As it stands, point-source only runs may not define photFilter - we need to handle that
    # That should be obvious, as mapFilters will only have one entry
    for filtDict in simConfig.parDict['mapFilters']:
        if filtDict['label'] == simConfig.parDict['photFilter']:
            break
    simConfig.parDict['mapFilters']=[filtDict] 
        
    # Filling maps with injected sources will be done when maps.preprocessMapDict is called by the filter object
    # Generate catalog here
    if filtDict['class'].find("ArnaudModel") != -1:
        mockCatalog=pipelines.makeMockClusterCatalog(config, writeCatalogs = False, verbose = False)[0]                
    elif filtDict['class'].find("BeamModel") != -1:
        raise Exception("Haven't implemented generating mock source catalogs here yet")
    else:
        raise Exception("Don't know how to generate injected source catalogs for filterClass '%s'" % (filtDict['class']))
    for mapDict in simConfig.unfilteredMapsDictList:
        mapDict['injectSources']=mockCatalog
                
    # NOTE: we need to zap ONLY specific maps for when we are running in parallel
    for extName in simConfig.extNames:
        mapFileNames=glob.glob(simRootOutDir+os.path.sep+"filteredMaps"+os.path.sep+"*#%s_*.fits" % (extName))
        for m in mapFileNames:
            os.remove(m)
            
    simImageDict=pipelines.filterMapsAndMakeCatalogs(simConfig, 
                                                     rootOutDir = simRootOutDir,
                                                     copyFilters = True)
    
    # Cross match the output with the input catalog - how close were we?
    print("pos recovery results")
    IPython.embed()
    sys.exit()
    
        
