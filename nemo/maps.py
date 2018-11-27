"""

This module contains tools for manipulating maps (e.g., conversion of units etc.).

"""

from astLib import *
from scipy import ndimage
from scipy import interpolate
import astropy.io.fits as pyfits
import numpy as np
import glob
import os
import sys
import math
import pyximport; pyximport.install()
import nemoCython
import time
import IPython
import nemo
from . import catalogs
np.random.seed()

#-------------------------------------------------------------------------------------------------------------
def convertToY(mapData, obsFrequencyGHz = 148):
    """Converts mapData (in delta T) at given frequency to yc 
    
    """
    # Convert into y: (deltaT/TCMB) = yf(x), y=(kB*sigmaT)/(me*c^2) * integral of (dl neTe)
    # f(x) = [x*coth(x/2) -4][1+relativistic correction for high T] (coth = 1/tanh)
    # x=h*freq/kB*TCMB
    # Maps are in delta T
    # Let's assume we're non-relativistic and isothermal
    # all units in kelvin, metres, Hz, kg, etc
    h=6.63e-34
    kB=1.38e-23
    sigmaT=6.6524586e-29
    me=9.11e-31
    c=3e8
    TCMB=2.726
    x=(h*obsFrequencyGHz*1e9)/(kB*TCMB)
    fx=x*(1.0/math.tanh(x/2))-4.0
    mapData=(mapData/(TCMB*1e6))/fx # remember, map is in uK  
    
    return mapData

#-------------------------------------------------------------------------------------------------------------
def convertToDeltaT(mapData, obsFrequencyGHz = 148):
    """Converts mapData (in yc) to deltaT (micro Kelvin) at given frequency.
    
    """
    # Convert into y: (deltaT/TCMB) = yf(x), y=(kB*sigmaT)/(me*c^2) * integral of (dl neTe)
    # f(x) = [x*coth(x/2) -4][1+relativistic correction for high T] (coth = 1/tanh)
    # x=h*freq/kB*TCMB
    # Maps are in delta T
    # Let's assume we're non-relativistic and isothermal
    # all units in kelvin, metres, Hz, kg, etc
    h=6.63e-34
    kB=1.38e-23
    sigmaT=6.6524586e-29
    me=9.11e-31
    c=3e8
    TCMB=2.726
    x=(h*obsFrequencyGHz*1e9)/(kB*TCMB)
    fx=x*(1.0/math.tanh(x/2))-4.0
    mapData=mapData*fx*(TCMB*1e6)   # into uK
    
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
                tileDeckFileNameLabel="userDefined_%.1f" % (parDict['tileOverlapDeg'])
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
                outFile.write('global color=green dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
                outFile.write("image\n")
                for c, name in zip(coordsList, extNames):
                    outFile.write('polygon(%d, %d, %d, %d, %d, %d, %d, %d) # text="%s"\n' % (c[0], c[2], c[0], c[3], c[1], c[3], c[1], c[2], name))
                outFile.close()
                    
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
                            clip=astImages.clipUsingRADecCoords(mapData, wcs, ra1, ra0, dec0, dec1)
                            print("... adding %s [%d, %d, %d, %d ; %d, %d] ..." % (name, ra1, ra0, dec0, dec1, ra0-ra1, dec1-dec0))
                            header=clip['wcs'].header.copy()
                            if 'tileNoiseRegions' in list(parDict.keys()) and name in list(parDict['tileNoiseRegions'].keys()):
                                noiseRAMin, noiseRAMax, noiseDecMin, noiseDecMax=parDict['tileNoiseRegions'][name]
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
    """Does preprocessing steps according to parameters in mapDict, returns the mapDict with additional
    keys added ['data', 'weights', 'wcs'].
    
    """
            
    img=pyfits.open(mapDict['mapFileName'], memmap = True)
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
    if 'weightsFileName' in list(mapDict.keys()) and mapDict['weightsFileName'] != None:
        wht=pyfits.open(mapDict['weightsFileName'], memmap = True)
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
    if 'surveyMask' in list(mapDict.keys()) and mapDict['surveyMask'] !=  None:
        smImg=pyfits.open(mapDict['surveyMask'])
        surveyMask=smImg[extName].data
    else:
        surveyMask=np.ones(data.shape)
        surveyMask[weights == 0]=0

    if 'pointSourceMask' in list(mapDict.keys()) and mapDict['pointSourceMask'] != None:
        psImg=pyfits.open(mapDict['pointSourceMask'])
        psMask=psImg[extName].data
    else:
        psMask=np.ones(data.shape)
            
    print("... opened map %s ..." % (mapDict['mapFileName']))
    
    # Optional map clipping
    if 'RADecSection' in list(mapDict.keys()) and mapDict['RADecSection'] != None:
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
        #astImages.saveFITS(diagnosticsDir+os.path.sep+'%d' % (mapDict['obsFreqGHz'])+"_weights.fits", weights, wcs)
    
    if 'CMBSimSeed' in list(mapDict.keys()):
        print("... filling map with CMB + noise sim ...")
        # This is the only part of this module that depends on signals
        from . import signals
        
        # The old flipper-based routine that did this took 190 sec versus 0.7 sec for enlib
        # NOTE: enlib/pixell imports here for now, to save having to install if we're not using it
        from pixell import enmap, utils, powspec
        import astropy.wcs as apywcs
        enlibWCS=apywcs.WCS(wcs.header)
        ps=powspec.read_spectrum(nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"planck_lensedCls.dat", scale = True)
        randMap=enmap.rand_map(data.shape, enlibWCS, ps, seed = mapDict['CMBSimSeed'])
        np.random.seed()    # Otherwise, we will end up with identical white noise...
        # Convolve with beam - NOTE: beam map only 4 x HWHM on a side (otherwise convolution very slow)
        beamData=np.loadtxt(mapDict['beamFileName']).transpose()
        profile1d=beamData[1]
        rArcmin=beamData[0]*60.0
        HWHMDeg=rArcmin[np.where(profile1d <= 0.5)[0][0]]/60.
        radiusDeg=HWHMDeg*4
        RADeg, decDeg=wcs.getCentreWCSCoords()
        clip=astImages.clipUsingRADecCoords(np.array(randMap), wcs, RADeg+radiusDeg, RADeg-radiusDeg, decDeg-radiusDeg, decDeg+radiusDeg)
        RADeg, decDeg=clip['wcs'].getCentreWCSCoords()
        degreesMap=nemoCython.makeDegreesDistanceMap(clip['data'], clip['wcs'], RADeg, decDeg, 0.5)
        beamSignalMap, inputSignalProperties=signals.makeBeamModelSignalMap(mapDict['obsFreqGHz'], degreesMap, 
                                                                    clip['wcs'], mapDict['beamFileName'])
        randMap=ndimage.convolve(np.array(randMap), beamSignalMap/beamSignalMap.sum()) 
        # Add white noise that varies according to inv var map...
        # Noise needed is the extra noise we need to add to match the real data, scaled by inv var map
        # This initial estimate is too high, so we use a grid search to get a better estimate
        mask=np.nonzero(data)
        dataSigma=data[mask].std()
        whiteNoiseLevel=np.zeros(weights.shape)
        whiteNoiseLevel[mask]=1/np.sqrt(weights[mask])
        noiseNeeded=np.sqrt(data[mask].var()-randMap[mask].var()-np.median(whiteNoiseLevel[mask])**2)
        noiseBoostFactor=noiseNeeded/np.median(whiteNoiseLevel[mask])
        simNoise=np.zeros(data.shape)+np.array(randMap)
        simNoise[np.equal(data, 0)]=0.
        generatedNoise=np.random.normal(0, whiteNoiseLevel[mask], whiteNoiseLevel[mask].shape)
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
        # --- disabled
        simNoise[mask]=simNoise[mask]+generatedNoise*bestBoostFactor
        # Output
        data=np.array(simNoise)
        # Sanity check
        outFileName=diagnosticsDir+os.path.sep+"CMBSim_%d#%s.fits" % (mapDict['obsFreqGHz'], extName) 
        astImages.saveFITS(outFileName, data, wcs)
        
    # Optional adding of white noise
    if 'addNoise' in list(mapDict.keys()) and mapDict['addNoise'] != None:
        data=addWhiteNoise(data, mapDict['addNoise'])
        if diagnosticsDir != None:
            astImages.saveFITS(diagnosticsDir+os.path.sep+"simMapPlusNoise_%d.fits" \
                            % (mapDict['obsFreqGHz']), data, wcs)    
            
    # Optional background subtraction - subtract smoothed version of map, this is like high pass filtering
    # or a wavelet decomposition scale image
    if 'bckSubScaleArcmin' in list(mapDict.keys()) and mapDict['bckSubScaleArcmin'] != None:
        data=subtractBackground(data, wcs, smoothScaleDeg = mapDict['bckSubScaleArcmin']/60.)
    
    # Optional masking of point sources from external catalog - needed, e.g., for point source subtracted
    # maps from Jon's pipeline, because otherwise we get negative bits that are detected as spurious 
    # clusters
    if 'maskPointSourcesFromCatalog' in list(mapDict.keys()) and mapDict['maskPointSourcesFromCatalog'] != None:
        print("Add code for masking point sources in a catalog")
        IPython.embed()
        sys.exit()        
    
    # Apply a ready-made point source and in-paint - added Feb 2018 (should then remove other stuff above)
    # This only reduced bright edges around masked regions... which is nice, but not essential
    # May allow us to use smaller holes in PS mask though
    #smoothedMap=ndimage.median_filter(data, 13) # size chosen for typical hole size... slow... but quite good
    #data[np.where(psMask == 0)]=smoothedMap[np.where(psMask == 0)]
    
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
    
    # Save trimmed weights
    if os.path.exists(diagnosticsDir+os.path.sep+"weights.fits") == False:
        astImages.saveFITS(diagnosticsDir+os.path.sep+"weights.fits", weights, wcs)
        
    return mapDict
    
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
def estimateContaminationFromSkySim(imageDict, extNames, parDictFileName, numSkySims, diagnosticsDir):
    """Run the whole filtering set up again, on sky simulation, with noise, that we generate here.
    
    Turns out we need to run many realisations, as results can vary by a lot.
    
    We may want to combine the output from this with the inverted maps test (which is quicker and easier,
    but has different problems - e.g., how aggressive we are at masking point sources).
    
    Writes a plot and a .fits table to the diagnostics dir.
    
    Runs over both SNR and fixed_SNR values.
    
    Returns a dictionaries containing the results
    
    """
    
    # These may break in python3, so only here to limit damage...
    from . import filters
    from . import startUp
    
    # To ensure we use the same kernel for filtering the sim maps as was used on the real data, copy kernels to sims dir
    # The kernel will then be loaded automatically when filterMaps is called 
    # Yes, this is a bit clunky...
    rootOutDir=diagnosticsDir+os.path.sep+"skySim"
    kernelCopyDestDir=rootOutDir+os.path.sep+"diagnostics"
    dirList=[rootOutDir, kernelCopyDestDir]
    for d in dirList:
        if os.path.exists(d) == False:
            os.makedirs(d)
    for extName in extNames:
        fileNames=glob.glob(diagnosticsDir+os.path.sep+"kern2d*#%s*.fits" % (extName))
        for f in fileNames:
            shutil.copyfile(f, kernelCopyDestDir+os.path.sep+os.path.split(f)[-1]) 
                
    resultsList=[]
    for i in range(numSkySims):
        
        # NOTE: we throw the first sim away on figuring out noiseBoostFactors
        print(">>> sky sim %d/%d ..." % (i+1, numSkySims))
        t0=time.time()
        
        # We use the seed here to keep the CMB sky the same across frequencies...
        CMBSimSeed=np.random.randint(16777216)
        
        simParDict=startUp.parseConfigFile(parDictFileName)
            
        # Optional override of default GNFW parameters (used by Arnaud model), if used in filters given
        if 'GNFWParams' not in list(simParDict.keys()):
            simParDict['GNFWParams']='default'
        for filtDict in simParDict['mapFilters']:
            filtDict['params']['GNFWParams']=simParDict['GNFWParams']
        
        # We're feeding in extNames to be MPI-friendly (each process takes its own set of extNames)
        unfilteredMapsDictList, ignoreThis=makeTileDeck(simParDict)
        
        # Filling in with sim will be done when maps.preprocessMapDict is called by the filter object
        for mapDict in unfilteredMapsDictList:
            mapDict['CMBSimSeed']=CMBSimSeed
                    
        # NOTE: we need to zap ONLY specific maps for when we are running in parallel
        for extName in extNames:
            mapFileNames=glob.glob(rootOutDir+os.path.sep+"filteredMaps"+os.path.sep+"*#%s_*.fits" % (extName))
            for m in mapFileNames:
                os.remove(m)
        
        simImageDict=filters.filterMaps(unfilteredMapsDictList, simParDict['mapFilters'], extNames = extNames, rootOutDir = rootOutDir)
            
        # Below here is same as inverted maps right now....
        # If we have makeDS9Regions = True here, we overwrite the existing .reg files from when we ran on the non-inverted maps
        photometry.findObjects(simImageDict, threshold = simParDict['thresholdSigma'], minObjPix = simParDict['minObjPix'],
                               rejectBorder = simParDict['rejectBorder'], diagnosticsDir = diagnosticsDir,
                               invertMap = False, makeDS9Regions = True, useInterpolator = simParDict['useInterpolator'])    

        # For fixed filter scale
        # Adds fixed_SNR values to catalogs for all maps
        photometryOptions=simParDict['photometryOptions']
        if 'photFilter' in list(photometryOptions.keys()):
            photFilter=photometryOptions['photFilter']
        else:
            photFilter=None
        if photFilter != None:
            photometry.getSNValues(simImageDict, SNMap = 'file', prefix = 'fixed_', template = photFilter, invertMap = False)   
            SNRKeys=['SNR', 'fixed_SNR']
        else:
            SNRKeys=['SNR']
        
        minSNToIncludeInOptimalCatalog=simParDict['catalogCuts']
        catalogs.mergeCatalogs(simImageDict)
        catalogs.makeOptimalCatalog(simImageDict, minSNToIncludeInOptimalCatalog)

        # Write out sim map catalogs for debugging
        skySimDir=diagnosticsDir+os.path.sep+"skySim"
        if len(simImageDict['optimalCatalog']) > 0:
            tab=catalogs.catalogToTab(simImageDict['optimalCatalog'], catalogs.COLUMN_NAMES, ["SNR > 0.0"])    
            optimalCatalogFileName=skySimDir+os.path.sep+"skySim%d_optimalCatalog.csv" % (i)           
            catalogs.writeCatalogFromTab(tab, optimalCatalogFileName, \
                                            catalogs.COLUMN_NAMES, catalogs.COLUMN_FORMATS, constraintsList = ["SNR > 0.0"], 
                                            headings = True)
            addInfo=[{'key': 'SNR', 'fmt': '%.1f'}, {'key': 'fixed_SNR', 'fmt': '%.1f'}]
            catalogs.catalog2DS9(tab, optimalCatalogFileName.replace(".csv", ".reg"), constraintsList = ["SNR > 0.0"], \
                                    addInfo = addInfo, color = "cyan") 

        # Contamination estimate...
        contaminTabDict=estimateContamination(simImageDict, imageDict, SNRKeys, 'skySim', diagnosticsDir)
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
    extNamesLabel="#"+str(extNames).replace("[", "").replace("]", "").replace("'", "").replace(", ", "#")
    for k in list(avContaminTabDict.keys()):
        fitsOutFileName=diagnosticsDir+os.path.sep+"%s_contaminationEstimate_%s.fits" % (k, extNamesLabel)
        if os.path.exists(fitsOutFileName) == True:
            os.remove(fitsOutFileName)
        contaminTab=avContaminTabDict[k]
        contaminTab.write(fitsOutFileName)
        
    return avContaminTabDict

#-------------------------------------------------------------------------------------------------------------
def estimateContaminationFromInvertedMaps(imageDict, extNames, thresholdSigma, minObjPix, rejectBorder, 
                                          minSNToIncludeInOptimalCatalog, photometryOptions, diagnosticsDir, findCenterOfMass = True):
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
    
    # If we have makeDS9Regions = True here, we overwrite the existing .reg files from when we ran on the non-inverted maps
    photometry.findObjects(invertedDict, threshold = thresholdSigma, minObjPix = minObjPix,
                           rejectBorder = rejectBorder, diagnosticsDir = diagnosticsDir,
                           invertMap = True, makeDS9Regions = False, findCenterOfMass = findCenterOfMass)    

    # For fixed filter scale
    # Adds fixed_SNR values to catalogs for all maps
    if 'photFilter' in list(photometryOptions.keys()):
        photFilter=photometryOptions['photFilter']
    else:
        photFilter=None
    if photFilter != None:
        photometry.getSNValues(invertedDict, SNMap = 'file', prefix = 'fixed_', template = photFilter, invertMap = True)
        SNRKeys=['SNR', 'fixed_SNR']
    else:
        SNRKeys=['SNR']
        
    catalogs.mergeCatalogs(invertedDict)
    catalogs.makeOptimalCatalog(invertedDict, minSNToIncludeInOptimalCatalog)
    
    contaminTabDict=estimateContamination(invertedDict, imageDict, SNRKeys, 'invertedMap', diagnosticsDir)

    for k in list(contaminTabDict.keys()):
        fitsOutFileName=diagnosticsDir+os.path.sep+"%s_contaminationEstimate.fits" % (k)
        if os.path.exists(fitsOutFileName) == True:
            os.remove(fitsOutFileName)
        contaminTab=contaminTabDict[k]
        contaminTab.write(fitsOutFileName)
        
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