"""

Automated tiling using given survey mask - aims for 7.5 x 5 deg tiles (RA, dec degrees).

Writes tiles.reg file (for debugging), tiles.yml (for pasting into nemo .yml)

This should be incorporated into nemo.maps.makeTileDeck

NOTE: Assumes CAR pixelization (will break if not true)

"""

import os
import sys
import numpy as np
import astropy.io.fits as pyfits
from astLib import *
from scipy import ndimage
import IPython

# Will not get exactly these
targetTileHeight=5.0
targetTileWidth=7.5

with pyfits.open("surveyMask.fits.gz") as img:
    segMap=img[0].data
    wcs=astWCS.WCS(img[0].header, mode = 'pyfits')

segMap, numObjects=ndimage.label(np.greater(segMap, 0))
fieldIDs=np.arange(1, numObjects+1)

tileList=[]
for f in fieldIDs:
    
    ys, xs=np.where(segMap == f)
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
        assert(tileWidth < 10)
    
        for j in range(numCols):
            RALeft=RAMax-j*tileWidth
            RARight=RAMax-(j+1)*tileWidth
            if RALeft < 0:
                RALeft=RALeft+360
            if RARight < 0:
                RARight=RARight+360
            # HACK: Edge-of-map handling
            if RARight < 180.01 and RALeft < 195 and RALeft > 180.01:
                RARight=180.01
            tileList.append({'tileName': '%d_%d_%d' % (f, i, j), 
                             'RADecSection': [RARight, RALeft, decBottom, decTop]})

# For now, just write to .yml that we can copy + paste into the nemo config .yml
with open("tiles.yml", "w") as outFile:
    outFile.write("tileDefLabel: 'tiles_SO'\n")
    outFile.write("tileDefinitions:\n")
    for tileDict in tileList:
        outFile.write("- {tileName: '%s', RADecSection: %s}\n" % (tileDict['tileName'], 
                                                                  tileDict['RADecSection']))

# Make DS9 .reg file (pixel coords) - this is from maps.makeTileDeck
tileNames=[]
coordsList=[]
for tileDict in tileList:
    ra0, ra1, dec0, dec1=tileDict['RADecSection']
    x0, y0=wcs.wcs2pix(ra0, dec0)
    x1, y1=wcs.wcs2pix(ra1, dec1)
    xMin=min([x0, x1])
    xMax=max([x0, x1])
    yMin=min([y0, y1])
    yMax=max([y0, y1])
    coordsList.append([xMin, xMax, yMin, yMax])
    tileNames.append(tileDict['tileName'])   
                
with open("tiles.reg", "w") as outFile:
    outFile.write("# Region file format: DS9 version 4.1\n")
    outFile.write('global color=blue dashlist=8 3 width=1 font="helvetica 10 normal roman" select=1 highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 include=1 source=1\n')
    outFile.write("image\n")
    for c, name in zip(coordsList, tileNames):
        outFile.write('polygon(%d, %d, %d, %d, %d, %d, %d, %d) # text="%s"\n' % (c[0], c[2], c[0], c[3], c[1], c[3], c[1], c[2], name))
                
    
    
