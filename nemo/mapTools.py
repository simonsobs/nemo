"""

This module contains tools for manipulating maps (e.g. filtering, conversion of units etc.).

"""

from nemo import mapFilters
from nemo import photometry
from nemo import catalogTools
from nemo import simsTools
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
                    print("... WARNING: same tiling not guaranteed across multiple frequencies ...")
                wcs=astWCS.WCS(wht[0].header, mode = 'pyfits')
                tileOverlapDeg=parDict['tileOverlapDeg']
   
                if defineTilesAutomatically == True:
                    
                    if 'surveyMask' in mapDict.keys() and mapDict['surveyMask'] == None:
                        print "... WARNING: same tiling not guaranteed across multiple frequencies ..."
                    
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
def filterMaps(unfilteredMapsDictList, filtersList, extNames = ['PRIMARY'], rootOutDir = ".", verbose = True):
    """Build and applies filters to the unfiltered maps(s). The output is a filtered map in yc. All filter
    operations are done in the filter objects, even if multifrequency (a change from previous behaviour).
   
    Filtered maps are written to rootOutDir/filteredMaps
    Filters, if stored, are written to rootOutDir/filters
    
    Returns a dictionary containing a map of filtered maps to keys in filterDict. We'll use this dictionary
    for keeping track of catalogs etc. subsequently.
    
    """
    
    # Storage, in case it doesn't already exist
    filteredMapsDir=rootOutDir+os.path.sep+"filteredMaps"
    diagnosticsDir=rootOutDir+os.path.sep+"diagnostics"
    dirList=[filteredMapsDir, diagnosticsDir]
    for d in dirList:
        if os.path.exists(d) == False:
            os.makedirs(d)
            
    # Dictionary to keep track of images we're going to make
    imageDict={}
    
    # For handling tileDeck style .fits files
    imageDict['extNames']=extNames
    
    # Since we're putting stuff like extNames in the top level, let's keep a separate list of mapDicts
    imageDict['mapKeys']=[]
    
    # Make filtered maps for each filter
    if verbose == True: print(">>> Making filtered maps and S/N maps ...")
    for f in filtersList:
        
        # Iterate over all extensions (for tileDeck files)...
        for extName in extNames:
            
            print("--> extName = %s ..." % (extName))
            label=f['label']+"#"+extName
            
            filteredMapFileName=filteredMapsDir+os.path.sep+"%s_filteredMap.fits"  % (label)
            SNMapFileName=filteredMapsDir+os.path.sep+"%s_SNMap.fits" % (label)
            signalMapFileName=diagnosticsDir+os.path.sep+"%s_signalMap.fits" % (label)
            #transferFnFileName=filteredMapsDir+os.path.sep+"%s_transferFunction.fits" % (f['label'])

            if os.path.exists(filteredMapFileName) == False:
                
                print("... making filtered map %s ..." % (label)) 
                filterClass=eval('mapFilters.%s' % (f['class']))
                filterObj=filterClass(label, unfilteredMapsDictList, f['params'], \
                                      extName = extName, 
                                      diagnosticsDir = diagnosticsDir)
                filteredMapDict=filterObj.buildAndApply()
                    
                # Keywords we need for photometry later
                #filteredMapDict['wcs'].header['BBIAS']=filteredMapDict['beamDecrementBias']
                #filteredMapDict['wcs'].header['ASCALING']=filteredMapDict['signalAreaScaling']
                filteredMapDict['wcs'].header['BUNIT']=filteredMapDict['mapUnits']
                filteredMapDict['wcs'].updateFromHeader()

                #if filteredMapDict['obsFreqGHz'] != 'yc':
                    #filteredMapDict['data']=convertToY(filteredMapDict['data'], \
                                                    #obsFrequencyGHz = filteredMapDict['obsFreqGHz'])
                                                
                astImages.saveFITS(filteredMapFileName, filteredMapDict['data'], filteredMapDict['wcs'])
                astImages.saveFITS(SNMapFileName, filteredMapDict['SNMap'], filteredMapDict['wcs'])            
                #astImages.saveFITS(signalMapFileName, filteredMapDict['signalMap'], filteredMapDict['wcs'])            

            else:
                print("... filtered map %s already made ..." % (label)) 
            
            # Add file names to imageDict
            if label not in imageDict:
                imageDict[label]={}
            imageDict[label]['filteredMap']=filteredMapFileName
            imageDict[label]['SNMap']=SNMapFileName
            imageDict[label]['signalMap']=signalMapFileName
            
            # Track e.g. reference filter scale with this key
            imageDict[label]['template']=f['label']
            
            # Track which keys have filtered maps that we might want to iterate over
            imageDict['mapKeys'].append(label)
            
            # May be handy to keep track of for plotting etc. later
            imageDict[label]['unfilteredMapsDictList']=unfilteredMapsDictList  
            
    return imageDict
    
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
def clipUsingRADecCoords(imageData, imageWCS, RAMin, RAMax, decMin, decMax, returnWCS = True):
    """Customised version of the astLib routine, because there is a wraparound problem for 10h field in the
    equator. Not sure if the fix used here will be generally stable for e.g. the south...
    
    This *might* add a shift of +/- 1 pixel - we know for sure that the equivalent astLib routine doesn't on 
    unfiltered equatorial maps, so use that instead (should be fine now the map WCS has been tweaked).
    
    """
    
    # No idea if this is stable, but it works for the ACT season23 equatorial map
    wcs=imageWCS.copy()
    rac, decc=wcs.getCentreWCSCoords()
    xc, yc=wcs.wcs2pix(rac, decc)
    wcs.header.update('CRPIX1', xc)
    wcs.header.update('CRVAL1', rac)
    wcs.updateFromHeader()
        
    imHeight=imageData.shape[0]
    imWidth=imageData.shape[1]
    
    xMin, yMin=wcs.wcs2pix(RAMin, decMin)
    xMax, yMax=wcs.wcs2pix(RAMax, decMax)
    
    xMin=int(round(xMin))
    xMax=int(round(xMax))
    yMin=int(round(yMin))
    yMax=int(round(yMax))
    X=[xMin, xMax]
    X.sort()
    Y=[yMin, yMax]
    Y.sort()
    
    if X[0] < 0:
        X[0]=0
    if X[1] > imWidth:
        X[1]=imWidth
    if Y[0] < 0:
        Y[0]=0
    if Y[1] > imHeight:
        Y[1]=imHeight   
    
    clippedData=imageData[Y[0]:Y[1],X[0]:X[1]]

    # Update WCS
    if returnWCS == True:
        try:
            oldCRPIX1=wcs.header['CRPIX1']
            oldCRPIX2=wcs.header['CRPIX2']
            clippedWCS=wcs.copy()
            clippedWCS.header.update('NAXIS1', clippedData.shape[1])
            clippedWCS.header.update('NAXIS2', clippedData.shape[0])
            clippedWCS.header.update('CRPIX1', oldCRPIX1-X[0])
            clippedWCS.header.update('CRPIX2', oldCRPIX2-Y[0])
            clippedWCS.updateFromHeader()
            
        except KeyError:
            
            if REPORT_ERRORS == True:
                
                print("WARNING: astImages.clipUsingRADecCoords() : no CRPIX1, CRPIX2 keywords found - not updating clipped image WCS.")
                
                clippedData=imageData[Y[0]:Y[1],X[0]:X[1]]
                clippedWCS=wcs.copy()
    else:
        clippedWCS=None
    
    return {'data': clippedData, 'wcs': clippedWCS, 'clippedSection': [X[0], X[1], Y[0], Y[1]]}
    
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
        if weights.ndim == 3:
            weights=weights[0, :]
    else:
        weights=np.ones(data.shape)

    # Load survey and point source masks, if given
    if 'surveyMask' in list(mapDict.keys()) and mapDict['surveyMask'] !=  None:
        smImg=pyfits.open(mapDict['surveyMask'])
        surveyMask=smImg[extName].data
    else:
        surveyMask=np.ones(data.shape)
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
        # The old flipper-based routine that did this took 190 sec versus 0.7 sec for enlib
        # NOTE: enlib imports here for now, to save having to install if we're not using it
        from enlib import enmap, utils, powspec
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
        beamSignalMap, inputSignalProperties=simsTools.makeBeamModelSignalMap(mapDict['obsFreqGHz'], degreesMap, 
                                                                    clip['wcs'], mapDict['beamFileName'])
        randMap=ndimage.convolve(np.array(randMap), beamSignalMap.data/beamSignalMap.data.sum()) 
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
    
    # Optional removal of point sources, using GaussianWienerFilter to find them
    # We only do this once, and cache the result in the diagnostics dir
    # NOTE: see above for new 'pointSourceMask' option - eventually we may remove the below...
    if 'pointSourceRemoval' in list(mapDict.keys()) and mapDict['pointSourceRemoval'] != None:
        outFileName=diagnosticsDir+os.path.sep+"psRemoved_%d.fits" % (mapDict['obsFreqGHz'])
        if os.path.exists(outFileName) == False:
            psRemovalMapDict={}
            psRemovalMapDict['data']=subtractBackground(data, wcs, smoothScaleDeg = 3.0/60.0)
            psRemovalMapDict['wcs']=wcs
            psRemovalMapDict['weights']=weights
            psRemovalMapDict['obsFreqGHz']=mapDict['obsFreqGHz']
            if 'beamFWHMArcmin' in list(mapDict.keys()):
                psRemovalMapDict['beamFWHMArcmin']=mapDict['beamFWHMArcmin']
                psRemovalClass=mapFilters.GaussianMatchedFilter
            if 'beamFileName' in list(mapDict.keys()):
                psRemovalMapDict['beamFileName']=mapDict['beamFileName']
                psRemovalClass=mapFilters.BeamMatchedFilter
            psRemovalParams=mapDict['pointSourceRemoval']
            psRemovalParams['FWHMArcmin']=1.4
            psRemovalParams['noiseParams']={'method': 'dataMap'}
            gaussFilter=psRemovalClass('psremoval-%d' % (psRemovalMapDict['obsFreqGHz']), \
                                                            [psRemovalMapDict], psRemovalParams, \
                                                            diagnosticsDir = diagnosticsDir)
            psRemoved=gaussFilter.buildAndApply()
            SNMap=psRemoved['SNMap']
            imageDict={'psRemoved': {}}
            imageDict['psRemoved']['SNMap']=SNMap
            imageDict['psRemoved']['wcs']=psRemoved['wcs']
            photometry.findObjects(imageDict, SNMap = 'array', 
                                threshold = mapDict['pointSourceRemoval']['threshold'],
                                minObjPix = mapDict['pointSourceRemoval']['minObjPix'], 
                                rejectBorder = mapDict['pointSourceRemoval']['rejectBorder'], 
                                makeDS9Regions = False, 
                                writeSegmentationMap = False)            
            if diagnosticsDir != None:
                outFileName=diagnosticsDir+os.path.sep+"pointSources-%d.reg" % \
                            (psRemovalMapDict['obsFreqGHz'])
                catalogTools.catalog2DS9(imageDict['psRemoved']['catalog'], outFileName)
                baseKeys=['name', 'RADeg', 'decDeg']
                baseFormats=["%s", "%.6f", "%.6f"]
                catalogTools.writeCatalog(imageDict['psRemoved']['catalog'], 
                                            outFileName.replace(".reg", ".csv"), baseKeys, baseFormats, 
                                            [], headings = True, writeNemoInfo = False)
            maskedDict=maskOutSources(data, psRemoved['wcs'], imageDict['psRemoved']['catalog'],
                                        radiusArcmin = psRemovalParams['radiusArcmin'], 
                                        mask = psRemovalParams['masking'])
            if diagnosticsDir != None:
                astImages.saveFITS(diagnosticsDir+os.path.sep+"psRemoved_%d.fits" \
                            % (mapDict['obsFreqGHz']), maskedDict['data'], wcs)
                astImages.saveFITS(diagnosticsDir+os.path.sep+"psMask_%d.fits" \
                            % (mapDict['obsFreqGHz']), maskedDict['mask'], wcs)
            data=maskedDict['data']
        else:
            # Do we need this? Probably not.
            img=pyfits.open(outFileName)
            data=img[0].data
    
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
    
