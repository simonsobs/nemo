"""Save only the T part of Sigurd's 3d map

"""

from astLib import *
import glob
import pyfits

fileNames=glob.glob("*.fits")
for mapFileName in fileNames:
    if mapFileName.find("hits") == -1:
        img=pyfits.open(mapFileName)
        d=img[0].data[0]
        wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
        astImages.saveFITS("TOnly_"+mapFileName, d, wcs)
