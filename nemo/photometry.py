"""

This module contains source finding and photometry routines.

"""

import astropy.io.fits as pyfits
import astropy.table as atpy
import numpy as np
import numpy.fft as fft
from astLib import *
import os
import pylab
import math
from scipy import ndimage
from scipy import interpolate
from . import catalogs
from . import maps
import sys
import IPython
np.random.seed()

#------------------------------------------------------------------------------------------------------------
def findObjects(imageDict, SNMap = 'file', threshold = 3.0, minObjPix = 3, rejectBorder = 10, 
                findCenterOfMass = True, makeDS9Regions = True, writeSegmentationMap = False, 
                selFnDir = None, invertMap = False, objIdent = 'ACT-CL', longNames = False, 
                verbose = True, useInterpolator = True, measureShapes = False):
    """Finds objects in the filtered maps pointed to by the imageDict. Threshold is in units of sigma 
    (as we're using S/N images to detect objects). Catalogs get added to the imageDict.
    
    SNMap == 'file' means that the SNMap will be loaded from disk, SNMap == 'array' means it is a np
    array already in memory.
    
    Throw out objects within rejectBorder pixels of edge of frame if rejectBorder != None
    
    Set invertMap == True to do a test of estimating fraction of spurious sources
    
    Set measureShapes == True to fit ellipses in a similar style to how SExtractor does it (may be useful for
    automating finding of extended sources).
    
    """
    
    if verbose == True:
        print(">>> Finding objects ...")
    
    if rejectBorder == None:
        rejectBorder=0
                  
    # Do search on each filtered map separately
    for key in imageDict['mapKeys']:
        
        if verbose == True:
            print("... searching %s ..." % (key))
        
        # Load area mask
        tileName=key.split("#")[-1]
        if selFnDir != None:
            img=pyfits.open(selFnDir+os.path.sep+"areaMask#%s.fits.gz" % (tileName))
            areaMask=img[0].data
        else:
            areaMask=None
        
        if SNMap == 'file':
            img=pyfits.open(imageDict[key]['SNMap'])
            wcs=astWCS.WCS(imageDict[key]['SNMap'])
            data=img[0].data
        elif SNMap == 'array':
            data=imageDict[key]['SNMap']
            wcs=imageDict[key]['wcs']

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
        if findCenterOfMass == True:
           print("... working with center of mass method ...")
           objPositions=ndimage.center_of_mass(data, labels = segmentationMap, index = objIDs)
        else:
           print("... working with maximum position method ...")
           objPositions=ndimage.maximum_position(data, labels = segmentationMap, index = objIDs)

        objNumPix=ndimage.sum(sigPixMask, labels = segmentationMap, index = objIDs)

        # In the past this had problems - Simone using happily in his fork though, so now optional
        if useInterpolator == True:
            mapInterpolator=interpolate.RectBivariateSpline(np.arange(data.shape[0]), 
                                                            np.arange(data.shape[1]), 
                                                            data, kx = 3, ky = 3)
        else:
            print("... not using sub-pixel interpolation for SNRs ...")
                                                        
        # Border around edge where we might throw stuff out just to cut down contamination
        if type(areaMask) == np.ndarray and areaMask.sum() > 0:
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
        
        # Build catalog - initially as a list of dictionaries, then converted to astropy Table
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
                if longNames == False:
                    objDict['name']=catalogs.makeACTName(objDict['RADeg'], objDict['decDeg'], prefix = objIdent)
                else:
                    objDict['name']=catalogs.makeLongName(objDict['RADeg'], objDict['decDeg'], prefix = objIdent)                    
                objDict['numSigPix']=objNumPix[i]
                objDict['template']=key
                if useInterpolator == True:
                    objDict['SNR']=mapInterpolator(objDict['y'], objDict['x'])[0][0]
                else:
                    objDict['SNR']=data[int(round(objDict['y'])), int(round(objDict['x']))]
                # Optional SExtractor style shape measurements                
                if measureShapes == True:
                    doubleCheck=False
                    if objDict['numSigPix'] > 9:
                        mask=np.equal(segmentationMap, objIDs[i])
                        ys, xs=np.where(segmentationMap == objIDs[i])
                        yMin=ys.min()
                        xMin=xs.min()
                        yMax=ys.max()
                        xMax=xs.max()
                        # Centres (1st order moments) - cx2, cy2 here to avoid overwriting whatever catalog cx, cy is (see above)
                        xs=xs-xMin
                        ys=ys-yMin
                        cx2=(xs*data[mask]).sum()/data[mask].sum()
                        cy2=(ys*data[mask]).sum()/data[mask].sum()
                        # Spread (2nd order moments)
                        x2=(((xs**2)*data[mask]).sum()/data[mask].sum())-cx2**2
                        y2=(((ys**2)*data[mask]).sum()/data[mask].sum())-cy2**2
                        xy=(((xs*ys)*data[mask]).sum()/data[mask].sum())-cx2*cy2
                        # A, B, theta from above - correct theta has same sign as xy (see SExtractor manual) 
                        theta=np.degrees(np.arctan(2*(xy/(x2-y2)))/2.0)
                        if xy > 0 and theta < 0:
                            theta=theta+90
                        elif xy < 0 and theta > 0:
                            theta=theta-90
                        if theta > 0 and xy > 0:
                            doubleCheck=True
                        if theta < 0 and xy < 0:
                            doubleCheck=True
                        if doubleCheck == True:
                            A=np.sqrt((x2+y2)/2.0 + np.sqrt( ((x2-y2)/2)**2 + xy**2))
                            B=np.sqrt((x2+y2)/2.0 - np.sqrt( ((x2-y2)/2)**2 + xy**2))
                            # Moments work terribly for low surface brightness, diffuse things which aren't strongly peaked
                            # Shape measurement is okay though - so just scale A, B to match segMap area
                            segArea=float(np.count_nonzero(np.equal(segmentationMap, objIDs[i])))
                            curArea=A*B*np.pi
                            scaleFactor=np.sqrt(segArea/curArea)
                            A=A*scaleFactor
                            B=B*scaleFactor  
                            ecc=np.sqrt(1-B**2/A**2)
                            objDict['ellipse_PA']=theta
                            objDict['ellipse_A']=A
                            objDict['ellipse_B']=B
                            objDict['ellipse_x0']=cx2+xMin
                            objDict['ellipse_y0']=cy2+yMin
                            objDict['ellipse_e']=ecc
                    if objDict['numSigPix'] <= 9 or doubleCheck == False:
                        objDict['ellipse_PA']=-99
                        objDict['ellipse_A']=-99
                        objDict['ellipse_B']=-99
                        objDict['ellipse_x0']=-99
                        objDict['ellipse_y0']=-99
                        objDict['ellipse_e']=-99
                # Add to catalog (masks now applied before we get here so no need to check)
                catalog.append(objDict)
                idNumCount=idNumCount+1     
        
        # From here on, catalogs should be astropy Table objects...
        if len(catalog) > 0:
            catalog=catalogs.catalogListToTab(catalog)
            if makeDS9Regions == True:
                catalogs.catalog2DS9(catalog, imageDict[key]['filteredMap'].replace(".fits", ".reg"))
        imageDict[key]['catalog']=catalog

#------------------------------------------------------------------------------------------------------------
def getSNValues(imageDict, SNMap = 'file', useInterpolator = True, invertMap = False, prefix = '', 
                template = None):
    """Measures SNR values in maps at catalog positions.
    
    Set invertMap == True, to do a test for estimating the spurious source fraction
    
    Use prefix to add prefix before SNR in catalog (e.g., fixed_SNR)
    
    Use template to select a particular filtered map (i.e., label of mapFilter given in .par file, e.g.,
    a filter used to measure properties for all objects at fixed scale).
    
    """
            
    print(">>> Getting %sSNR values ..." % (prefix))
    
    # Do search on each filtered map separately
    for key in imageDict['mapKeys']:
        
        print("... searching %s ..." % (key))
        
        if template == None:
            templateKey=key
        else:
            templateKey=None
            for k in imageDict['mapKeys']:
                # Match the reference filter name and the tileName
                if k.split("#")[0] == template and k.split("#")[-1] == key.split("#")[-1]:
                    templateKey=k
            if templateKey == None:
                raise Exception("didn't find templateKey")
        
        if SNMap == 'file':
            img=pyfits.open(imageDict[templateKey]['SNMap'])
            wcs=astWCS.WCS(imageDict[templateKey]['SNMap'])
            data=img[0].data
        elif SNMap == 'array':
            data=imageDict[templateKey]['SNMap']
            wcs=imageDict[templateKey]['wcs']
        else:
            raise Exception("Didn't understand SNMap value '%s'" % (str(SNMap)))

        # This is for checking contamination
        if invertMap == True:
            data=data*-1
            
        # In the past this had problems - Simone using happily in his fork though, so now optional
        if useInterpolator == True:
            mapInterpolator=interpolate.RectBivariateSpline(np.arange(data.shape[0]), 
                                                            np.arange(data.shape[1]), 
                                                            data, kx = 3, ky = 3)
        else:
            print("... not using sub-pixel interpolation for SNRs ...")
        
        if len(imageDict[key]['catalog']) > 0:
            imageDict[key]['catalog'].add_column(atpy.Column(np.zeros(len(imageDict[key]['catalog'])), prefix+'SNR'))
        for obj in imageDict[key]['catalog']:
            x, y=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
            if int(x) > 0 and int(x) < data.shape[1] and int(y) > 0 and int(y) < data.shape[0]:
                if useInterpolator == True:
                    obj[prefix+'SNR']=mapInterpolator(y, x)[0][0]
                else:
                    obj[prefix+'SNR']=data[int(round(y)), int(round(x))] # read directly off of S/N map
       
#------------------------------------------------------------------------------------------------------------
def measureFluxes(imageDict, photFilter, diagnosticsDir, useInterpolator = True, unfilteredMapsDict = None):
    """Add flux measurements to each catalog pointed to in the imageDict. Measured in 'outputUnits' 
    specified in the filter definition in the .par file (and written to the image header as 'BUNIT').
    
    Maps should have been normalised before this stage to make the value of the filtered map at the object
    position equal to the flux that is wanted (yc, uK or Jy/beam). This includes taking out the effect of
    the beam.
    
    Use 'photFilter' to choose which filtered map in which to make flux measurements at fixed scale
    (fixed_delta_T_c, fixed_y_c etc.). Set to None if you don't want these.

    """
    
    print(">>> Doing photometry ...")

    # Adds fixed_SNR values to catalogs for all maps
    if photFilter is not None:
        getSNValues(imageDict, SNMap = 'file', prefix = 'fixed_', template = photFilter, 
                    useInterpolator = useInterpolator)

    for key in imageDict['mapKeys']:
        
        print("--> Map: %s ..." % (imageDict[key]['filteredMap']))
        catalog=imageDict[key]['catalog']        
        
        img=pyfits.open(imageDict[key]['filteredMap'])
        wcs=astWCS.WCS(imageDict[key]['filteredMap'])
        mapData=img[0].data

        mapUnits=wcs.header['BUNIT']                # could be 'yc' or 'Jy/beam'
        
        # These keywords would only be written if output units 'uK' (point source finding)
        # This also relies on the beam solid angle being in the header of the beam .txt file
        if 'BEAMNSR' in wcs.header.keys():
            beamSolidAngle_nsr=wcs.header['BEAMNSR']
            obsFreqGHz=wcs.header['FREQGHZ']
            reportJyFluxes=True
        else:
            reportJyFluxes=False
        
        # In the past this had problems - Simone using happily in his fork though, so now optional
        if useInterpolator == True:
            mapInterpolator=interpolate.RectBivariateSpline(np.arange(mapData.shape[0]), 
                                                            np.arange(mapData.shape[1]), 
                                                            mapData, kx = 3, ky = 3)
        else:
            print("... not using sub-pixel interpolation for fluxes ...")
            mapInterpolator=None
        
        # Add fixed filter scale maps
        mapDataList=[mapData]
        interpolatorList=[mapInterpolator]
        prefixList=['']
        tileName=key.split("#")[-1]
        if photFilter is not None:
            photImg=pyfits.open(imageDict[photFilter+"#"+tileName]['filteredMap'])
            photMapData=photImg[0].data
            mapDataList.append(photMapData)
            prefixList.append('fixed_')
            if useInterpolator == True:
                photMapInterpolator=interpolate.RectBivariateSpline(np.arange(photMapData.shape[0]), 
                                                                    np.arange(photMapData.shape[1]), 
                                                                    photMapData, kx = 3, ky = 3)
            else:
                photMapInterpolator=None
            interpolatorList.append(photMapInterpolator)

        keysToAdd=['deltaT_c', 'err_deltaT_c']
        if mapUnits == 'yc':
            keysToAdd=keysToAdd+['y_c', 'err_y_c']
        elif mapUnits == 'uK':
            keysToAdd=keysToAdd+['fluxJy', 'err_fluxJy']        
        for prefix in prefixList:
            for k in keysToAdd:
                if len(catalog) > 0:
                    catalog.add_column(atpy.Column(np.zeros(len(catalog)), prefix+k))

        for obj in catalog:
                        
            if key == obj['template']:
                x, y=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
                for data, prefix, interpolator in zip(mapDataList, prefixList, interpolatorList):
                    # NOTE: We might want to avoid 2d interpolation here because that was found not to be robust elsewhere
                    # 2018: Simone seems to now be using this happily, so now optional
                    if useInterpolator == True:
                        mapValue=interpolator(y, x)[0][0]
                    else:
                        mapValue=data[int(round(y)), int(round(x))]
                    # NOTE: remember, all normalisation should be done when constructing the filtered maps, i.e., not here!
                    if mapUnits == 'yc':
                        yc=mapValue
                        deltaTc=maps.convertToDeltaT(yc, obsFrequencyGHz = 148.0)
                        obj[prefix+'y_c']=yc/1e-4                            # So that same units as H13 in output catalogs
                        obj[prefix+'err_y_c']=obj[prefix+'y_c']/obj[prefix+'SNR']
                        obj[prefix+'deltaT_c']=deltaTc
                        obj[prefix+'err_deltaT_c']=abs(deltaTc/obj[prefix+'SNR'])
                    elif mapUnits == 'Y500':
                        print("add photometry.measureFluxes() for Y500")
                        IPython.embed()
                        sys.exit()
                    elif mapUnits == 'uK':
                        # For this, we want deltaTc to be source amplitude
                        deltaTc=mapValue
                        obj[prefix+'deltaT_c']=deltaTc
                        obj[prefix+'err_deltaT_c']=deltaTc/obj[prefix+'SNR']                        
                        if reportJyFluxes == True:
                            obj[prefix+"fluxJy"]=deltaTToJySr(obj[prefix+'deltaT_c'], obsFreqGHz)*beamSolidAngle_nsr*1.e-9
                            obj[prefix+"err_fluxJy"]=deltaTToJySr(obj[prefix+'err_deltaT_c'], obsFreqGHz)*beamSolidAngle_nsr*1.e-9

#------------------------------------------------------------------------------------------------------------
def addFreqWeightsToCatalog(imageDict, photFilter, diagnosticsDir):
    """Add relative weighting by frequency for each object in the optimal catalog, extracted from the 
    data cube saved under diagnosticsDir (this is made by makeSZMap in RealSpaceMatchedFilter). This is
    needed for multi-frequency cluster finding / analysis - for, e.g., weighting relativistic corrections to
    y0~ (which was estimated from the inverse variance weighted average of y0~ from each frequency map).
    
    NOTE: this is only applied for the reference filter (pointed to by photFilter)
    
    """

    if photFilter is None:
        return None
    
    catalog=imageDict['optimalCatalog']
    for tileName in imageDict['tileNames']:
        label=photFilter+"#"+tileName
        freqWeightMapFileName=diagnosticsDir+os.path.sep+"freqRelativeWeights_%s.fits" % (label)
        if os.path.exists(freqWeightMapFileName) == True:
            img=pyfits.open(freqWeightMapFileName)
            freqCube=img[0].data
            wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
            obsFreqGHzDict={}
            for key in wcs.header.keys():
                if key.find("FREQGHZ") != -1:
                    obsFreqGHzDict[int(key.split("FREQGHZ")[0])]=wcs.header[key]
            for i in range(len(obsFreqGHzDict.keys())):
                freqStr=("%.1f" % (obsFreqGHzDict[i])).replace(".", "p")    # MongoDB doesn't like '.' in key names
                keyToAdd='fixed_y_c_weight_%sGHz' % (freqStr)
                if keyToAdd not in catalog.keys():
                    catalog.add_column(atpy.Column(np.zeros(len(catalog)), keyToAdd))
            for row in catalog:
                if row['template'].split("#")[-1] == tileName:
                    x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
                    x=int(round(x))
                    y=int(round(y))
                    freqWeights=freqCube[:, y, x]
                    for i in range(len(obsFreqGHzDict.keys())):
                        freqStr=("%.1f" % (obsFreqGHzDict[i])).replace(".", "p")    # MongoDB doesn't like '.' in key names
                        row['fixed_y_c_weight_%sGHz' % (freqStr)]=freqWeights[i]

#------------------------------------------------------------------------------------------------------------
def deltaTToJySr(temp, obsFreqGHz):
    """Convert delta T (uK) to Jy/sr at the given frequency in GHz
    
    """
    
    # We should enforce consistency of constants used in various places at some point...
    kB=1.380658e-16
    h=6.6260755e-27
    c=29979245800.
    nu=obsFreqGHz*1.e9
    T0=2.726
    x=h*nu/(kB*T0)
    cNu=2*(kB*T0)**3/(h**2*c**2)*x**4/(4*(np.sinh(x/2.))**2)
    cNu*=1e23
    
    return temp*cNu*1e-6/T0

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
def makeAnnulus(innerScalePix, outerScalePix):
    """Makes the footprint of an annulus for feeding into the rank filter
    
    Returns annulus, numValues
    
    """

    # Make the footprint for the rank filter
    #bckInnerScalePix=beamSigmaPix*3#int(round((bckInnerRadiusArcmin/60.0)/yPixScale))
    #bckOuterScalePix=beamSigmaPix*5#int(round((bckOuterRadiusArcmin/60.0)/yPixScale))
    bckInnerScalePix=int(round(innerScalePix))
    bckOuterScalePix=int(round(outerScalePix))
    annulus=np.zeros([bckOuterScalePix*2, bckOuterScalePix*2])
    xRange=np.array([np.arange(annulus.shape[1])-bckOuterScalePix]*annulus.shape[0])
    yRange=np.array([np.arange(annulus.shape[0])-bckOuterScalePix]*annulus.shape[1]).transpose()
    rRange=np.sqrt(xRange**2+yRange**2)
    annulus[np.logical_and(np.greater(rRange, bckInnerScalePix),
                                np.less(rRange, bckOuterScalePix))]=1.0
    #numValues=annulus.flatten().nonzero()[0].shape[0]
    
    return annulus.astype('int64')

