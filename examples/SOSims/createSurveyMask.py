"""

Just creates the survey mask image from a DS9 .reg file.
    
"""

import os
import sys
import argparse
import astropy.io.fits as pyfits
from astLib import *
import numpy as np
import mahotas.polygon
import IPython

#-------------------------------------------------------------------------------------------------------------
# Main
parser=argparse.ArgumentParser()
parser.add_argument("headerTemplate")
parser.add_argument("regionFile")
args=parser.parse_args()

header=pyfits.header.Header.fromtextfile(args.headerTemplate)
wcs=astWCS.WCS(header, mode = 'pyfits')
shape=(header['NAXIS2'], header['NAXIS1'])

outFileName=args.regionFile.replace(".reg", ".fits")

with open(args.regionFile) as inFile:
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

# New method - need for > 2 Gb images
surveyMask=np.zeros(shape, dtype = int)
for polyPoints in polyList:
    mahotas.polygon.fill_polygon(polyPoints, surveyMask)
astImages.saveFITS(outFileName, surveyMask, wcs)

