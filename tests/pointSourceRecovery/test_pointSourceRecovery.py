"""

This script takes the E-D56 map, injects a bunch of point sources into it, runs nemo, and checks that the
average recovered source amplitude is unbiased. We also check the accuracy of source position recovery.

NOTE: This will check for and download the E-D56 maps (using wget) if not found.

"""

import os
import numpy as np
from pixell import utils, pointsrcs, enmap
from astLib import *
import astropy.io.fits as pyfits
import astropy.table as atpy
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import pylab as plt
import IPython

#------------------------------------------------------------------------------------------------------------
# Test threholds
FLUX_RATIO_THRESH=0.01         # Require median recovered source amplitude to be within 1% of input

#------------------------------------------------------------------------------------------------------------
def makeSimSourceCatalog(mapData, wcs, numSources = 1000):
    """Generates a source catalog with a reasonable flux distribution. Sources will be randomly placed
    in the footprint of the provided map.
    
    Args:
        mapData: A 2d array. Sources coordinates will be drawn from random non-zero locations in this.
        wcs: An astWCS object
        numSources: The number of sources in the catalog
    
    Returns:
        An astropy Table object
    
    """
    
    # Random fluxes (in delta T uK)
    # This is chosen to vaguely match that seen in deep56 f090
    I=np.random.lognormal(np.log(600), 1.1, numSources)
    
    # Random positions
    ys, xs=np.where(mapData != 0)
    ys=ys+np.random.uniform(0, 1, len(ys))
    xs=xs+np.random.uniform(0, 1, len(xs))
    indices=np.random.randint(0, len(ys), numSources)
    coords=wcs.pix2wcs(xs[indices], ys[indices])
    coords=np.array(coords)
    
    tab=atpy.Table()
    tab.add_column(atpy.Column(coords[:, 0], "ra"))
    tab.add_column(atpy.Column(coords[:, 1], "dec"))
    tab.add_column(atpy.Column(I, "I"))
    
    return tab
    
#------------------------------------------------------------------------------------------------------------
# Main

beamFileName="../../examples/equD56/profiles_ACT/profile_AR1_2009_pixwin_130224.txt"

# Download the E-D56 map if not found
inMapFileName="weightedMap_4.fits"
if os.path.exists(inMapFileName) == False:
    print(">>> Downloading E-D56 map ...")
    os.system("wget https://www.acru.ukzn.ac.za/~mjh/equD56Maps.tar.gz")
    os.system("tar -zxvf equD56Maps.tar.gz")

# Inject sources
outMapFileName="sourceInjectedMap.fits"
nsigma=5.0
smul=1.0

img=pyfits.open(inMapFileName)
mapData=img[0].data
if mapData.ndim > 2:
    mapData=mapData[0]
wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
imap=enmap.read_map(inMapFileName)

tab=makeSimSourceCatalog(mapData, wcs)
outCatalogFileName="inputSourcesCatalog.fits"
if os.path.exists(outCatalogFileName) == True:
    os.remove(outCatalogFileName)
tab.write(outCatalogFileName)
srcs=pointsrcs.tab2recarray(tab)

beam=pointsrcs.read_beam(beamFileName)
beam[0]*=utils.degree
srcparam=pointsrcs.src2param(srcs)

omap=pointsrcs.sim_srcs(imap.shape, imap.wcs, srcparam, beam, imap, smul=smul, nsigma=nsigma, pixwin=True)
enmap.write_map(outMapFileName, omap)

# Run nemo and match with input catalog
os.system("nemo PSTest_E-D56.yml")

inTab=atpy.Table().read(outCatalogFileName)
outTab=atpy.Table().read("PSTest_E-D56"+os.path.sep+"PSTest_E-D56_optimalCatalog.fits")
cat1=SkyCoord(ra = inTab['ra'].data, dec = inTab['dec'].data, unit = 'deg')
radiusArcmin=2.5    # For filtering out implausible matches
xMatchRadiusDeg=radiusArcmin/60.
cat2=SkyCoord(ra = outTab['RADeg'].data, dec = outTab['decDeg'].data, unit = 'deg')
xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
mask=np.less(rDeg.value, xMatchRadiusDeg)

# Source amplitude recovery check
print(">>> Amplitude recovery test:")
inT=inTab['I'][mask]
diffT=inTab['I'][mask]-outTab[xIndices]['deltaT_c'][mask]
ratioT=inTab['I'][mask]/outTab[xIndices]['deltaT_c'][mask]
SNRs=outTab[xIndices]['SNR'][mask]
print("... median amplitude in / amplitude out = %.6f" % (np.median(ratioT)))
if abs(1.0-np.median(ratioT)) > FLUX_RATIO_THRESH:
    plt.plot(inT, ratioT, 'r.')
    plt.xlabel("input deltaT (uK)")
    plt.ylabel("input / output deltaT")
    plt.savefig("amplitudeRecovery.png")
    raise Exception("... FAILED")
else:
    print("... PASSED")

# Position recovery check

# Map subtraction / residual check... make template of recovered sources and subtract from sourceInjectedMap.fits

# Clean up
os.system("rm -r PSTest_E-D56")
os.system("rm sourceInjectedMap.fits inputSourcesCatalog.fits")
