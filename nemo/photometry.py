# -*- coding: utf-8 -*-
"""This module contains source finding and photometry routines.

"""

import pyfits
import numpy as np
import numpy.fft as fft
from astLib import *
import os
import pylab
import math
from scipy import ndimage
from scipy import interpolate
import catalogTools
import mapTools
import simsTools
import sys
import IPython
np.random.seed()

#------------------------------------------------------------------------------------------------------------
def findObjects(imageDict, SNMap = 'file', threshold = 3.0, minObjPix = 3, rejectBorder = 10, 
                makeDS9Regions = True, writeSegmentationMap = False, diagnosticsDir = None,
                invertMap = False):
    """Finds objects in the filtered maps pointed to by the imageDict. Threshold is in units of sigma 
    (as we're using S/N images to detect objects). Catalogs get added to the imageDict.
    
    SNMap == 'file' means that the SNMap will be loaded from disk, SNMap == 'array' means it is a np
    array already in memory. SNMap == 'local', means use the a ranked filter to get noise in annulus around
    given position. This doesn't seem to work as well...
    
    Throw out objects within rejectBorder pixels of edge of frame if rejectBorder != None
    
    Set invertMap == True to do a test of estimating fraction of spurious sources
    
    """
    
    print ">>> Finding objects ..."
    
    if rejectBorder == None:
        rejectBorder=0
  
    # Load point source mask - this is hacked for 148 only at the moment
    if diagnosticsDir != None and os.path.exists(diagnosticsDir+os.path.sep+"psMask_148.fits") == True:
        img=pyfits.open(diagnosticsDir+os.path.sep+"psMask_148.fits")
        psMaskMap=img[0].data
    else:
        psMaskMap=None
    
    # Load area mask
    if diagnosticsDir != None:
        img=pyfits.open(diagnosticsDir+os.path.sep+"areaMask.fits")
        areaMask=img[0].data
    else:
        areaMask=None
            
    # Do search on each filtered map separately
    for key in imageDict.keys():
        
        print "... searching %s ..." % (key)
        
        if SNMap == 'file':
            img=pyfits.open(imageDict[key]['SNMap'])
            wcs=astWCS.WCS(imageDict[key]['SNMap'])
            data=img[0].data
        elif SNMap == 'array':
            data=imageDict[key]['SNMap']
            wcs=imageDict[key]['wcs']
        elif SNMap == 'local':
            # This ignores the area mask, for now
            img=pyfits.open(imageDict[key]['ycFilteredMap'])
            wcs=astWCS.WCS(imageDict[key]['ycFilteredMap'])
            data=img[0].data
            bckInnerRadiusArcmin=10.0
            bckOuterRadiusArcmin=15.0
            ra0, dec0=wcs.getCentreWCSCoords()
            x0, y0=wcs.wcs2pix(ra0, dec0)
            x1=x0+1
            y1=y0+1
            ra1, dec1=wcs.pix2wcs(x1, y1)
            yPixScale=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
            bckInnerScalePix=int(round((bckInnerRadiusArcmin/60.0)/yPixScale))
            bckOuterScalePix=int(round((bckOuterRadiusArcmin/60.0)/yPixScale))
            annulus=np.zeros([bckOuterScalePix*2, bckOuterScalePix*2])
            xRange=np.array([np.arange(annulus.shape[1])-bckOuterScalePix]*annulus.shape[0])
            yRange=np.array([np.arange(annulus.shape[0])-bckOuterScalePix]*annulus.shape[1]).transpose()
            rRange=np.sqrt(xRange**2+yRange**2)
            annulus[np.logical_and(np.greater(rRange, bckInnerScalePix),
                                        np.less(rRange, bckOuterScalePix))]=1.0
            numValues=annulus.flatten().nonzero()[0].shape[0]
                    
            noise=ndimage.rank_filter(abs(data), int(round(0.683*numValues)), footprint = annulus)
            data=data/noise # now data is SNMap
            data[np.isnan(data)]=0.0      
            data[np.isinf(data)]=0.0            
        else:
            raise Exception, "Didn't understand SNMap value '%s'" % (str(SNMap))

        # This is for checking contamination
        if invertMap == True:
            data=data*-1
            
        # Thresholding to identify significant pixels
        sigPix=np.array(np.greater(data, threshold), dtype=int)
        sigPixMask=np.equal(sigPix, 1)
        
        # Fast, simple segmentation - don't know about deblending, but doubt that's a problem for us
        segmentationMap, numObjects=ndimage.label(sigPix)
        if writeSegmentationMap == True:
            segMapFileName=imageDict[key]['SNMap'].replace("SNMap", "segmentationMap")
            astImages.saveFITS(segMapFileName, segmentationMap, wcs)
            imageDict[key]['segmentationMap']=segMapFileName
        
        # Get object positions, number of pixels etc.
        objIDs=np.unique(segmentationMap)
        objPositions=ndimage.center_of_mass(data, labels = segmentationMap, index = objIDs)
        objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)

        mapInterpolator=interpolate.RectBivariateSpline(np.arange(data.shape[0]), 
                                                        np.arange(data.shape[1]), 
                                                        data, kx = 1, ky = 1)
                                                        
        # Border around edge where we might throw stuff out just to cut down contamination
        if areaMask != None:
            minX=np.where(areaMask > 0)[1].min()
            maxX=np.where(areaMask > 0)[1].max()
            minY=np.where(areaMask > 0)[0].min()
            maxY=np.where(areaMask > 0)[0].max()
        else:
            minX=0
            maxX=segmentationMap.shape[1]-1
            minY=0
            maxY=segmentationMap.shape[0]-1
        minX=minX+rejectBorder
        maxX=maxX-rejectBorder
        minY=minY+rejectBorder
        maxY=maxY-rejectBorder       
        
        # Build catalog
        catalog=[]
        idNumCount=1
        for i in range(len(objIDs)):
            if type(objNumPix) != float and objNumPix[i] > minObjPix:
                objDict={}
                objDict['id']=idNumCount
                objDict['x']=objPositions[i][1]
                objDict['y']=objPositions[i][0]
                objDict['RADeg'], objDict['decDeg']=wcs.pix2wcs(objDict['x'], objDict['y'])
                if objDict['RADeg'] < 0:
                    objDict['RADeg']=360.0+objDict['RADeg']
                galLong, galLat=astCoords.convertCoords("J2000", "GALACTIC", objDict['RADeg'], objDict['decDeg'], 2000)
                objDict['galacticLatDeg']=galLat
                objDict['name']=catalogTools.makeACTName(objDict['RADeg'], objDict['decDeg'])
                objDict['numSigPix']=objNumPix[i]
                objDict['template']=key
                objDict['SNR']=mapInterpolator(objDict['y'], objDict['x'])[0][0]              
                if objDict['x'] > minX and objDict['x'] < maxX and \
                    objDict['y'] > minY and objDict['y'] < maxY:
                    masked=False
                    if psMaskMap != None:
                        if psMaskMap[int(round(objDict['y'])), int(round(objDict['x']))] > 0:
                            masked=True
                    if masked == False:
                        catalog.append(objDict)
                        idNumCount=idNumCount+1     
        if makeDS9Regions == True:
            catalogTools.catalog2DS9(catalog, imageDict[key]['ycFilteredMap'].replace(".fits", ".reg"))
        imageDict[key]['catalog']=catalog

#------------------------------------------------------------------------------------------------------------
def getSNValues(imageDict, SNMap = 'file', invertMap = False):
    """Measures S/N values in maps at catalog positions. Like above, but used by e.g. nemoSPTCheck script
    
    Set invertMap == True, to do a test for estimating the spurious source fraction
    
    """
    
    print ">>> Getting S/N values ..."
          
    # Do search on each filtered map separately
    for key in imageDict.keys():
        
        print "... searching %s ..." % (key)

        if SNMap == 'file':
            img=pyfits.open(imageDict[key]['SNMap'])
            wcs=astWCS.WCS(imageDict[key]['SNMap'])
            data=img[0].data
        elif SNMap == 'array':
            data=imageDict[key]['SNMap']
            wcs=imageDict[key]['wcs']
        else:
            raise Exception, "Didn't understand SNMap value '%s'" % (str(SNMap))

        mapInterpolator=interpolate.RectBivariateSpline(np.arange(data.shape[0]), 
                                                        np.arange(data.shape[1]), 
                                                        data, kx = 1, ky = 1)
                                            
        for obj in imageDict[key]['catalog']:
            obj['x'], obj['y']=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
            if obj['x'] > 0 and obj['x'] < data.shape[1] and obj['y'] > 0 and obj['y'] < data.shape[0]:
                #obj['SNR']=data[obj['y'], obj['x']] # read directly off of S/N map
                obj['SNR']=mapInterpolator(obj['y'], obj['x'])[0][0]
            else:
                obj['SNR']=None
       
#------------------------------------------------------------------------------------------------------------
def measureFluxes(imageDict, photometryOptions, diagnosticsDir, unfilteredMapsDict = None):
    """Add flux measurements to each catalog pointed to in the imageDict. Measured in filtered maps in yc
    units. The interface for using the unfiltered maps option here is somewhat confused.
    
    NOTE: Jan 2015, flux is now integrated Y500 assuming R500 as given by the template for each object
    in the catalog if using the Arnaud profile (i.e., we normalised things there to give Y500).
    
    Use beamDecrementBias to undo beam smoothing for central yc, deltaTc values
    Use signalAreaScaling to multiply up map such that pixel at object location gives Y500 
    (if the cluster did have that mass and redshift)
    
    """

    # NOTE: I prefer sr to arcmin2 (but arcmin2 is how we define signalAreaScaling for Arnaud model)
    srToArcmin2=np.power(np.radians(1.0/60.0), 2)
    
    print ">>> Measuring SZ decrements and Y500s ..."
    for key in imageDict.keys():
        
        print "--> Map: %s ..." % (imageDict[key]['ycFilteredMap'])
        catalog=imageDict[key]['catalog']        
        
        img=pyfits.open(imageDict[key]['ycFilteredMap'])
        wcs=astWCS.WCS(imageDict[key]['ycFilteredMap'])
        mapData=img[0].data

        beamDecrementBias=wcs.header['BBIAS']       # multiply map values by this to get intrinsic yc, delta Tc
        signalAreaScaling=wcs.header['ASCALING']    # multiply map values by this to turn single pixel at object location to Y500
        
        mapInterpolator=interpolate.RectBivariateSpline(np.arange(mapData.shape[0]), 
                                                        np.arange(mapData.shape[1]), 
                                                        mapData, kx = 1, ky = 1) 
                
        # Maps must have been normalised such that peak pixel value for an object gives Y500 in arcmin2
        for obj in catalog:
            
            if key == obj['template']:
                
                # NOTE: We might want to avoid 2d interpolation here because that was found not to be robust elsewhere
                # i.e., avoid using interpolate.RectBivariateSpline
                mapValue=mapData[int(round(obj['y'])), int(round(obj['x']))]
                #mapValue=mapInterpolator(obj['y'], obj['x'])[0][0]
                #print "photometry"
                #IPython.embed()
                #sys.exit()
                Y500_sr=mapValue*signalAreaScaling*srToArcmin2
                yc=mapValue*beamDecrementBias
                deltaTc=mapTools.convertToDeltaT(yc, obsFrequencyGHz = 148.0)
                
                obj['Y500_sr']=Y500_sr
                obj['err_Y500_sr']=Y500_sr/obj['SNR']
                obj['y_c']=yc
                obj['err_y_c']=yc/obj['SNR']
                obj['deltaT_c']=deltaTc
                obj['err_deltaT_c']=deltaTc/obj['SNR']
        
#------------------------------------------------------------------------------------------------------------ 
def measureApertureFluxes(catalog, apertureRadiusArcmin, mapData, wcs, fluxCorrectionFactor, numBackgrounds = 10, 
                          maxBackgroundSepArcmin = 30.0, minBackgroundSepArcmin = 10.0,
                          allowOverlappingBackgrounds = True):
    """This is now doing the 'aperture-flux-of-signal-template-multiplied-by-norm-from-central-pixel thing.
    We assume here we already have SNR (we can use this and peak flux to estimate statistical uncertainty).
    
    Also measure delta T here, so should rename this
    """

    mapInterpolator=interpolate.RectBivariateSpline(np.arange(mapData.shape[0]), 
                                                    np.arange(mapData.shape[1]), 
                                                    mapData, kx = 1, ky = 1)                                                    
    for obj in catalog:

        # y_c, delta T c - assuming 148 GHz here
        yc=mapInterpolator(obj['y'], obj['x'])[0][0]
        deltaT=mapTools.convertToDeltaT(yc, obsFrequencyGHz = 148.0)
        obj['y_c']=yc
        obj['deltaT_c']=deltaT
        
        #if obj['name'] == 'ACT-CL J0059.1-0049':
            #print "WARNING: check in MatchedFilter code about normalisation of signalMap"
            #print "Here, check y_c stuff - make it y0 ish like Hass"
            #IPython.embed()
            #sys.exit()
            # y0 in Hasselfield_upp .fits table is corrected y0, i.e., NOT y0~
            # y0 = y0~ / (Q*fRel)
            # So, e.g., mapTools.convertToY(-336)/0.53 where -336 is deltaT0_fixed in Hass catalogue, reproduces y0
            #(obj['y_c']*Q)/1e-4
        
        # Integrated Y, summed flux in aperture method (works in white noise case)
        fluxArcmin2=objectFluxInAperture(obj, apertureRadiusArcmin, mapData, wcs, fluxCorrectionFactor)
        
        # Integrated Y, integrate profile method - insert it here, then comment out above
        
        # Error is a fudge, and wrong no doubt
        fluxErrArcmin2=fluxArcmin2/obj['SNR']
        
        # Add stuff to catalog
        obj['flux_arcmin2']=fluxArcmin2
        obj['fluxErr_arcmin2']=fluxErrArcmin2
        obj['fluxStatus']="Okay"
        obj['fluxRadius_arcmin']=apertureRadiusArcmin    
        
#------------------------------------------------------------------------------------------------------------
def getRadialDistanceMap(objDict, data, wcs):
    """Returns an array of same dimensions as data but containing radial distance to the object specified
    in the objDict in units degrees on the sky
    
    """
    
    x0=objDict['x']
    y0=objDict['y']
    x1=x0+1
    y1=y0+1
    ra1, dec1=wcs.pix2wcs(x1, y1)
    xPixScale=astCoords.calcAngSepDeg(objDict['RADeg'], objDict['decDeg'], ra1, objDict['decDeg'])
    yPixScale=astCoords.calcAngSepDeg(objDict['RADeg'], objDict['decDeg'], objDict['RADeg'], dec1)       
    xRange=np.array([np.arange(0, data.shape[1])-x0]*data.shape[0])*xPixScale
    yRange=(np.array([np.arange(0, data.shape[0])-y0]*data.shape[1])*yPixScale).transpose()
    rRange=np.sqrt(xRange**2+yRange**2)       

    return rRange

#------------------------------------------------------------------------------------------------------------
def getPixelsDistanceMap(objDict, data):
    """Returns an array of same dimensions as data but containing radial distance to the object specified
    in the objDict in units of pixels
    
    """
    
    x0=objDict['x']
    y0=objDict['y']
    x1=x0+1
    y1=y0+1     
    xRange=np.array([np.arange(0, data.shape[1])-x0]*data.shape[0])
    yRange=(np.array([np.arange(0, data.shape[0])-y0]*data.shape[1])).transpose()
    rRange=np.sqrt(xRange**2+yRange**2)       

    return rRange

#------------------------------------------------------------------------------------------------------------ 
def objectFluxInAperture(objDict, apertureRadiusArcmin, mapData, wcs, fluxCorrectionFactor = 1.0):
    """This routine simply sums the flux at the position of the objDict in the map, without doing
    any kind of background subtraction. Returned flux is in arcmin2.
                
    """
                
    # Clip out region big enough to contain object + background
    # This can handle objects near map edges
    ra0=objDict['RADeg']
    dec0=objDict['decDeg']
    x, y=wcs.wcs2pix(ra0, dec0)
    ra1, dec1=wcs.pix2wcs(x+1, y+1)    
    xLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    arcmin2PerPix=xLocalDegPerPix*yLocalDegPerPix*60.0**2
    clipSizeDeg=(apertureRadiusArcmin/60.0)*3
    xClipSizePix=int(round(clipSizeDeg/xLocalDegPerPix))
    yClipSizePix=int(round(clipSizeDeg/yLocalDegPerPix))
    left=int(round(x-xClipSizePix/2))
    right=int(round(x+xClipSizePix/2))
    top=int(round(y+yClipSizePix/2))
    bottom=int(round(y-yClipSizePix/2))
    yCorrection=0
    xCorrection=0
    if bottom < 0:
        yCorrection=bottom
        bottom=0
    if top > mapData.shape[0]:
        yCorrection=top-mapData.shape[0]
        top=mapData.shape[0]
    if left < 0:
        xCorrection=left
        left=0
    if right > mapData.shape[1]:
        xCorrection=right-mapData.shape[1]
        right=mapData.shape[1]
    data=mapData[bottom:top, left:right]
    if xCorrection < 0:
        x=xClipSizePix/2+xCorrection
    else:
        x=data.shape[1]-(xClipSizePix/2-xCorrection)
    if yCorrection < 0:
        y=yClipSizePix/2+yCorrection
    else:
        y=data.shape[0]-(yClipSizePix/2-yCorrection)
    xPix=np.array([np.arange(0, data.shape[1], dtype=float)]*data.shape[0])-x
    yPix=(np.array([np.arange(0, data.shape[0], dtype=float)]*data.shape[1])-y).transpose()
    xDeg=xPix*xLocalDegPerPix
    yDeg=yPix*yLocalDegPerPix
    rDeg=np.sqrt(xDeg**2+yDeg**2)  
    
    # Note flux correction factor in below
    fluxMask=np.less(rDeg, apertureRadiusArcmin/60.0)
    fluxPixels=data[fluxMask]
    flux_arcmin2=np.sum(fluxPixels)*arcmin2PerPix*fluxCorrectionFactor
        
    return flux_arcmin2
            
