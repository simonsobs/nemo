# -*- coding: utf-8 -*-
"""This module contains tools for manipulating maps (e.g. filtering, conversion of units etc.).

"""

import mapFilters
import photometry
import catalogTools
from astLib import *
from scipy import ndimage
from scipy import interpolate
import pyfits
import numpy as np
import os
import sys
import math
import pyximport; pyximport.install()
import nemoCython
import time
import IPython
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
def filterMaps(unfilteredMapsDictList, filtersList, rootOutDir = ".", verbose = True):
    """Build and applies filters to the unfiltered maps(s). The output is a filtered map in yc. All filter
    operations are done in the filter objects, even if multifrequency (a change from previous behaviour).
   
    Filtered maps are written to rootOutDir/filteredMaps
    Filters, if stored, are written to rootOutDir/filters
    
    Returns a dictionary containing a map of filtered maps to keys in filterDict. We'll use this dictionary
    for keeping track of catalogs etc. subsequently.
    
    """
    
    # Storage, in case it doesn't already exist
    filteredMapsDir=rootOutDir+os.path.sep+"filteredMaps"
    filtersDir=rootOutDir+os.path.sep+"filters"
    diagnosticsDir=rootOutDir+os.path.sep+"diagnostics"
    dirList=[filteredMapsDir, filtersDir, diagnosticsDir]
    for d in dirList:
        if os.path.exists(d) == False:
            os.makedirs(d)
            
    # Dictionary to keep track of images we're going to make
    imageDict={}
        
    # Make filtered maps for each filter
    if verbose == True: print ">>> Making filtered maps and S/N maps ..."
    for f in filtersList:
        
        filteredMapFileName=filteredMapsDir+os.path.sep+"%s_filteredMap.fits"  % (f['label'])
        SNMapFileName=filteredMapsDir+os.path.sep+"%s_SNMap.fits" % (f['label'])
        signalMapFileName=diagnosticsDir+os.path.sep+"%s_signalMap.fits" % (f['label'])
        #transferFnFileName=filteredMapsDir+os.path.sep+"%s_transferFunction.fits" % (f['label'])

        if os.path.exists(filteredMapFileName) == False:
            
            print "... making filtered map %s ..." % (f['label']) 
            filterClass=eval('mapFilters.%s' % (f['class']))
            filterObj=filterClass(f['label'], unfilteredMapsDictList, f['params'], \
                                  outDir = filtersDir, diagnosticsDir = diagnosticsDir)
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
            print "... filtered map %s already made ..." % (f['label']) 
        
        # Add file names to imageDict
        if f['label'] not in imageDict:
            imageDict[f['label']]={}
        imageDict[f['label']]['filteredMap']=filteredMapFileName
        imageDict[f['label']]['SNMap']=SNMapFileName
        imageDict[f['label']]['signalMap']=signalMapFileName
        
        # May be handy to keep track of for plotting etc. later
        imageDict[f['label']]['unfilteredMapsDictList']=unfilteredMapsDictList  
        
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
                    print "Fix oversubtraction... peakValue is pointSource + CMB..."
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
            print "Add code to subtract point sources"
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
                
                print "WARNING: astImages.clipUsingRADecCoords() : no CRPIX1, CRPIX2 keywords found - not updating clipped image WCS."
                
                clippedData=imageData[Y[0]:Y[1],X[0]:X[1]]
                clippedWCS=wcs.copy()
    else:
        clippedWCS=None
    
    return {'data': clippedData, 'wcs': clippedWCS, 'clippedSection': [X[0], X[1], Y[0], Y[1]]}
    
#-------------------------------------------------------------------------------------------------------------
def preprocessMapDict(mapDict, diagnosticsDir = None):
    """Does preprocessing steps according to parameters in mapDict, returns the mapDict with additional
    keys added ['data', 'weights', 'wcs'].
    
    """
    
    # Load data map, if we don't already have one (e.g. if removing point sources)
    if 'data' not in mapDict.keys():   
        img=pyfits.open(mapDict['mapFileName'], memmap = True)
        wcs=astWCS.WCS(mapDict['mapFileName'])
        data=img[0].data
        if mapDict['units'] == 'Jy/sr':
            if mapDict['obsFreqGHz'] == 148:
                data=(data/1.072480e+09)*2.726*1e6
            elif mapDict['obsFreqGHz'] == 219:
                data=(data/1.318837e+09)*2.726*1e6
            else:
                raise Exception, "no code added to support conversion to uK from Jy/sr for freq = %.0f GHz" \
                        % (mapDict['obsFreqGHz'])

        # Load weight map if given
        if 'weightsFileName' in mapDict.keys() and mapDict['weightsFileName'] != None:
            wht=pyfits.open(mapDict['weightsFileName'], memmap = True)
            weights=wht[0].data
        else:
            weights=np.ones(data.shape)

        # Load survey and point source masks, if given
        if 'surveyMask' in mapDict.keys() and mapDict['surveyMask'] !=  None:
            smImg=pyfits.open(mapDict['surveyMask'])
            surveyMask=smImg[0].data
        else:
            surveyMask=np.ones(data.shape)
        if 'pointSourceMask' in mapDict.keys() and mapDict['pointSourceMask'] != None:
            psImg=pyfits.open(mapDict['pointSourceMask'])
            psMask=psImg[0].data
        else:
            psMask=np.ones(data.shape)
                
        print "... opened map %s ..." % (mapDict['mapFileName'])
        
        # Optional map clipping
        if 'RADecSection' in mapDict.keys() and mapDict['RADecSection'] != None:
            RAMin, RAMax, decMin, decMax=mapDict['RADecSection']
            clip=astImages.clipUsingRADecCoords(data, wcs, RAMin, RAMax, decMin, decMax)
            data=clip['data']
            whtClip=astImages.clipUsingRADecCoords(weights, wcs, RAMin, RAMax, decMin, decMax)
            weights=whtClip['data']
            wcs=clip['wcs']
            psClip=astImages.clipUsingRADecCoords(psMask, wcs, RAMin, RAMax, decMin, decMax)
            psMask=psClip['data']
            surveyClip=astImages.clipUsingRADecCoords(surveyMask, wcs, RAMin, RAMax, decMin, decMax)
            surveyMask=surveyClip['data']
            #astImages.saveFITS(diagnosticsDir+os.path.sep+'%d' % (mapDict['obsFreqGHz'])+"_weights.fits", weights, wcs)
        
        # Optional adding of white noise
        if 'addNoise' in mapDict.keys() and mapDict['addNoise'] != None:
            data=addWhiteNoise(data, mapDict['addNoise'])
            if diagnosticsDir != None:
                astImages.saveFITS(diagnosticsDir+os.path.sep+"simMapPlusNoise_%d.fits" \
                                % (mapDict['obsFreqGHz']), data, wcs)    
                
        # Optional background subtraction - subtract smoothed version of map, this is like high pass filtering
        # or a wavelet decomposition scale image
        if 'bckSubScaleArcmin' in mapDict.keys() and mapDict['bckSubScaleArcmin'] != None:
            data=subtractBackground(data, wcs, smoothScaleDeg = mapDict['bckSubScaleArcmin']/60.)
        
        # Added March 2017 as replacement for older 'pointSourceRemoval' code
        # Break out from here / tidy up later
        # NOTE: If we just mask big enough areas, we don't need to bother with this!
        #if 'pointSourceMask' in mapDict.keys() and mapDict['pointSourceMask'] != None:
            #psRemovedFileName=diagnosticsDir+os.path.sep+"psMasked_%d.fits" % (mapDict['obsFreqGHz'])
            #if os.path.exists(psRemovedFileName) == True:
                #print "... loading cached map %s which has point source masking applied ..." % (psRemovedFileName)
                #psRemovedImg=pyfits.open(psRemovedFileName)
                #data=psRemovedImg[0].data
            #else:
                #print "... filling in map at masked point source locations ..."
                #t0=time.time()
                ## The big smooth method
                ##smoothed=smoothMap(data, wcs, smoothScaleDeg = 10.0/60.0)
                ## Find hole locations - we could swap out the mask here for a catalog instead, or add as option
                ## Measure the average value in an annulus around each hole location
                #ptSrcData=1-psMask
                #sigPix=np.array(np.greater(ptSrcData, 0), dtype=int)
                #sigPixMask=np.equal(sigPix, 1)
                #segmentationMap, numObjects=ndimage.label(sigPix)
                #objIDs=np.unique(segmentationMap)
                #objPositions=ndimage.center_of_mass(ptSrcData, labels = segmentationMap, index = objIDs)
                #objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)
                #psCat=[]
                #idNumCount=0
                #minObjPix=5
                #for i in range(len(objIDs)):
                    #if type(objNumPix) != float and objNumPix[i] > minObjPix:
                        #idNumCount=idNumCount+1
                        #x=int(round(objPositions[i][1]))
                        #y=int(round(objPositions[i][0]))
                        #objDict={}
                        #objDict['id']=idNumCount
                        #objDict['x']=x
                        #objDict['y']=y
                        #objDict['RADeg'], objDict['decDeg']=wcs.pix2wcs(objDict['x'], objDict['y'])
                        #objDict['numPix']=objNumPix[i]
                        #rPix=int(round(np.sqrt(objDict['numPix']/np.pi)))
                        #annulus=photometry.makeAnnulus(rPix, rPix+4)
                        #yMin=objDict['y']-annulus.shape[0]/2
                        #yMax=objDict['y']+annulus.shape[0]/2
                        #xMin=objDict['x']-annulus.shape[1]/2
                        #xMax=objDict['x']+annulus.shape[1]/2
                        ## Just skip if too close to map edge - will be cropped out later anyway...
                        #if yMax > data.shape[0] or yMin < 0 or xMax > data.shape[1] or xMin < 0:
                            #continue        
                        ## for testing
                        ##if x > 6900 and x < 6940 and y >  1240 and y < 1270:
                        ##if x > 6070 and x < 6133 and y > 601  and y < 648:
                        #psClip=data[yMin:yMax, xMin:xMax]
                        #sigma=np.std(psClip[np.not_equal(psClip, 0)*annulus])
                        #med=np.median(psClip[np.not_equal(psClip, 0)*annulus])
                        #if sigma > 0:
                            #pixMask=np.where(segmentationMap == objIDs[i])
                            ##data[pixMask]=smoothed[pixMask]+np.random.normal(0, sigma, len(smoothed[pixMask]))
                            #data[pixMask]=np.random.normal(med, sigma, len(data[pixMask]))
                            #psCat.append(objDict)
                #t1=time.time()
                #print "... took %.3f sec ..." % (t1-t0)
                #astImages.saveFITS(psRemovedFileName, data, wcs)
        
        # Optional removal of point sources, using GaussianWienerFilter to find them
        # We only do this once, and cache the result in the diagnostics dir
        # NOTE: see above for new 'pointSourceMask' option - eventually we may remove the below...
        if 'pointSourceRemoval' in mapDict.keys() and mapDict['pointSourceRemoval'] != None:
            outFileName=diagnosticsDir+os.path.sep+"psRemoved_%d.fits" % (mapDict['obsFreqGHz'])
            if os.path.exists(outFileName) == False:
                psRemovalMapDict={}
                psRemovalMapDict['data']=subtractBackground(data, wcs, smoothScaleDeg = 3.0/60.0)
                psRemovalMapDict['wcs']=wcs
                psRemovalMapDict['weights']=weights
                psRemovalMapDict['obsFreqGHz']=mapDict['obsFreqGHz']
                if 'beamFWHMArcmin' in mapDict.keys():
                    psRemovalMapDict['beamFWHMArcmin']=mapDict['beamFWHMArcmin']
                    psRemovalClass=mapFilters.GaussianMatchedFilter
                if 'beamFileName' in mapDict.keys():
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
        if 'maskPointSourcesFromCatalog' in mapDict.keys() and mapDict['maskPointSourcesFromCatalog'] != None:
            print "Add code for masking point sources in a catalog"
            IPython.embed()
            sys.exit()
            
        # Add the map data to the dict
        mapDict['data']=data
        mapDict['weights']=weights
        mapDict['wcs']=wcs
        mapDict['surveyMask']=surveyMask
        mapDict['psMask']=psMask
               
        # Save trimmed weights
        if os.path.exists(diagnosticsDir+os.path.sep+"weights.fits") == False:
            astImages.saveFITS(diagnosticsDir+os.path.sep+"weights.fits", weights, wcs)
        
    return mapDict
    
#-------------------------------------------------------------------------------------------------------------
def subtractBackground(data, wcs, smoothScaleDeg = 30.0/60.0):
    """Smoothes map with Gaussian of given scale and subtracts it, to get rid of large scale power.
    
    """
            
    smoothedData=smoothMap(data, wcs, smoothScaleDeg=smoothScaleDeg)
    data=data-smoothedData
    
    return data

#-------------------------------------------------------------------------------------------------------------
def smoothMap(data, wcs, smoothScaleDeg = 5.0/60.0):
    """Smoothes map with Gaussian of given scale.
    
    """
            
    ra0, dec0=wcs.getCentreWCSCoords()
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
    
