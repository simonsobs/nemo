"""

Creates survey masks: 
    1. with no holes (for tiling)
    2. with holes (extended sources + dust + misc artifacts + sources to 10 mJy)

The holes in mask 2 here have no effect on the noise term in the filter, so don't affect recovered object
positions or properties. They are simply just holes in the map where clusters will not be detected.

We're using catalogs to mask point sources etc. for the purposes of constructing the noise term in the 
matched filter. That's not done here (see the maskPointSourcesFromCatalog parameter in nemo .yml files).

This script relies on .reg files and catalogs to define the survey region and mask holes.
    
"""

import os
import sys
import astropy.table as atpy
import astropy.io.fits as pyfits
from astLib import *
import pylab as plt
import numpy as np
from scipy import ndimage
#from PIL import Image
#from PIL import ImageDraw
import mahotas.polygon
import glob
import nemoCython
from nemo import catalogs
import time
import IPython
plt.matplotlib.interactive(True)

#-------------------------------------------------------------------------------------------------------------
def addArtifactsToMask(regFile, surveyMask, wcs):
    """Given a .reg file of circular regions, make holes in surveyMask...
    
    """
    
    RAs=[]
    decs=[]
    rArcmin=[]
    with open(regFile) as inFile:
        lines=inFile.readlines()
        for line in lines:
            if line[0] != "#" and line.find("circle") != -1:
                bits=line.split("(")[-1].split('")')[0].split(",") 
                RAs.append(float(bits[0]))
                decs.append(float(bits[1]))
                rArcmin.append(float(bits[2])/60.)
    
    tab=atpy.Table()
    tab.add_column(atpy.Column(RAs, 'RADeg'))
    tab.add_column(atpy.Column(decs, 'decDeg'))
    tab.add_column(atpy.Column(rArcmin, 'rArcmin'))   
    rArcminMap=np.ones(surveyMask.shape, dtype = float)*1e6
    count=0
    for row in tab:
        count=count+1
        print(count)
        x, y=wcs.wcs2pix(row['RADeg'], row['decDeg'])
        rArcminMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(rArcminMap, wcs, 
                                                                       row['RADeg'], row['decDeg'], 
                                                                       row['rArcmin']/60.)
        rArcminMap=rArcminMap*60
        surveyMask[rArcminMap < row['rArcmin']]=0
    
    return surveyMask

#-------------------------------------------------------------------------------------------------------------
def makeSurveyMask(surveyMaskRegFileName, shape, wcs):
    """Makes a survey mask from a DS9 .reg file.
    
    Returns survey mask (2d array)
    
    """
    
    with open(surveyMaskRegFileName, "r") as inFile:
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
                    IPython.embed()
                    sys.exit()
                x, y=wcs.wcs2pix(RADeg, decDeg)
                polyPoints.append((int(round(y)), int(round(x))))
            polyList.append(polyPoints)
    surveyMask=np.zeros(shape, dtype = int)
    for polyPoints in polyList:
        mahotas.polygon.fill_polygon(polyPoints, surveyMask)
        
    return surveyMask

#-------------------------------------------------------------------------------------------------------------
def maskObject(mapToMask, degreesMap, wcs, RADeg, decDeg, maskRadiusDeg):
    """Set circular region in mapToMask to zero. Done in place.
    
    """
    degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, RADeg, decDeg, 
                                                                   maskRadiusDeg)
    mapToMask[degreesMap < maskRadiusArcsec/3600.0]=0

#-------------------------------------------------------------------------------------------------------------
# Main

# Options
# NOTE: Can have two different survey geometry masks - mask 1 here doesn't have galactic latitude cut
# We're doing this so that tiles remain the same across versions (could be changed for next release)
mask1RegFileName="AdvACTSurveyMask_v7_S18.reg"                  # Defines survey geometry / tiling (no holes)
mask2RegFileName="AdvACTSurveyMask_v7_galLatCut_S18.reg"        # Defines survey geometry (with holes)
addDust=True                                                    # Add in dust from P353 map
addArtifacts=False # was True                                   # Mask anything in artifactRegFileName
addExtended=True                                                # Apply Simone's extended objects mask
add10mJy=False # was True                                       # Post-processing masking of 10 mJy sources
useThreads=True
artifactsRegFileName="all_artifacts_20200405.reg"
sourceCatalogFileName="PS_S18d_202003_optimalCatalog.fits"

# We can't automatically mask 10 mJy sources... because some clusters have rings that are picked up
# So, we will remove objects from the source catalog that are within some radius of bright clusters
# NOTE: These values are based on not wanting to zap ACT-CL J0417.5-1154 (SNR ~ 25, fixed_yc = 3.3e-04)
clusterCatalogFileName="S18d_202003_optimalCatalog.fits"        
ycBrightCut=1.0
ycBrightExclusionRadiusArcmin=9.0

# New format for output file names
mask1FileName=mask1RegFileName.replace(".reg", ".fits")    # Survey mask with no holes written here
mask2FileName=mask2RegFileName.replace(".reg", ".fits")    # The mask with holes
if addDust == True:
    mask2FileName=mask2FileName.replace(".fits", "-dust.fits")
if addArtifacts == True:
    mask2FileName=mask2FileName.replace(".fits", "-artifacts.fits")
if addExtended == True:
    mask2FileName=mask2FileName.replace(".fits", "-extended.fits")
if add10mJy == True:
    mask2FileName=mask2FileName.replace(".fits", "-post10mJy.fits")

# NOTE: Now in Mar2020 directory
with pyfits.open("act_s08_s18_cmb_f150_daynight_map.fits") as img:
    wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
    shape=img[0].data[0].shape

# Survey mask geometry defined using ds9
# NOTE: now can handle multiple polygons
# mahotas is needed for > 2 Gb images
print(">>> Making mask 1: survey mask for tiling (no holes) ...")
surveyMask=makeSurveyMask(mask1RegFileName, shape, wcs)
astImages.saveFITS(mask1FileName, surveyMask, wcs)

# Drill optional holes
print(">>> Making mask 2: survey mask with holes ...")
surveyMask=makeSurveyMask(mask2RegFileName, shape, wcs)
if addExtended == True:
    print("... adding extended sources ...")
    with pyfits.open("mask_extsrc_gt4arcmin_20190416.fits") as img:
        extMask=1-img[0].data
    surveyMask=surveyMask-extMask
    
if addArtifacts == True:
    print("... adding artifacts from file = %s ..." % (artifactsRegFileName))
    surveyMask=addArtifactsToMask(artifactsRegFileName, surveyMask, wcs)
    
if addDust == True:
    print("... adding dusty regions ...")
    with pyfits.open("down4_P353_AdvACT.fits") as img:
        dust=np.array(np.greater(img[0].data, 0.004), dtype = int)
        dustWCS=astWCS.WCS(img[0].header, mode = 'pyfits')
    for y in range(surveyMask.shape[0]):
        print("... %d/%d ..." % (y+1, surveyMask.shape[0]))
        xs=np.arange(surveyMask.shape[1])
        coords=np.array(wcs.pix2wcs(xs, [y]*xs.shape[0]))
        RAs=coords[:, 0]
        decs=coords[:, 1]
        dustPix=np.array(np.round(dustWCS.wcs2pix(RAs, decs)), dtype = int)
        dustXs=dustPix[:, 0]
        dustYs=dustPix[:, 1]
        mask=np.logical_and(np.greater_equal(dustXs, 0), np.less(dustXs, dust.shape[1]))
        dustMask=np.array((dust[dustYs, dustXs] == 1), dtype = int)
        surveyMask[y, mask]=surveyMask[y, mask]-dustMask[mask]

if add10mJy == True:
    print("... adding 10 mJy sources from file = %s ..." % (sourceCatalogFileName))
    tab=atpy.Table().read(sourceCatalogFileName)
    tab=tab[tab['fluxJy']*1000 > 10]
    # Have to remove spurious sources (rings) that are too close to bright clusters
    clTab=atpy.Table().read(clusterCatalogFileName)
    clTab=clTab[clTab['fixed_y_c'] > ycBrightCut]
    for row in clTab:
        rDeg=astCoords.calcAngSepDeg(row['RADeg'], row['decDeg'], tab['RADeg'].data, tab['decDeg'].data)
        rGoodMask=(rDeg > ycBrightExclusionRadiusArcmin/60.)
        tab=tab[rGoodMask]
    # Now do the masking
    maskRadiusArcsec=320
    if useThreads == True:
        #mask2FileName=mask2FileName.replace(".fits", "_threaded.fits") # If testing
        import multiprocessing
        from concurrent.futures import ThreadPoolExecutor
        from itertools import repeat
        print("... num threads = %d ..." % (multiprocessing.cpu_count()))
        RAs=tab['RADeg'].data
        decs=tab['decDeg'].data
        degreesMap=np.ones(shape, dtype = float)*1e6   # Updated in place (fast)
        with ThreadPoolExecutor(max_workers = multiprocessing.cpu_count()) as executor:
            executor.map(maskObject, repeat(surveyMask), repeat(degreesMap), repeat(wcs), RAs, decs, repeat(maskRadiusArcsec/3600))
    else:
        degreesMap=np.ones(shape, dtype = float)*1e6   # Updated in place (fast)
        count=0        
        for row in tab:
            print("... %d/%d ..." % (count, len(tab)))
            count=count+1
            degreesMap, xBounds, yBounds=nemoCython.makeDegreesDistanceMap(degreesMap, wcs, 
                                                                        row['RADeg'], row['decDeg'], 
                                                                        maskRadiusArcsec/3600.0)
            surveyMask[degreesMap < maskRadiusArcsec/3600.0]=0
    
# Don't want -ve area (holes where sources outside survey mask)
surveyMask[np.less(surveyMask, 0)]=0
astImages.saveFITS(mask2FileName, surveyMask, wcs)

