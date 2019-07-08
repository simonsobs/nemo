"""

Makes redshifts.fits catalog (table of names, redshifts matched against nemo output catalog), to allow 
nemoMass to run.

"""

import os
import sys
import argparse
import astropy.table as atpy
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import numpy as np
import IPython

parser=argparse.ArgumentParser()
parser.add_argument("nemoCatalog")
parser.add_argument("redshiftCatalog")
parser.add_argument("-O", "--output-catalog", type=str, default="redshifts.fits")
args=parser.parse_args()

xMatchRadiusArcmin=1.0  # Conservative

tab=atpy.Table().read(args.nemoCatalog)
xTab=atpy.Table().read(args.redshiftCatalog)
cat1=SkyCoord(ra = tab['RADeg'], dec = tab['decDeg'], unit = 'deg')
xMatchRadiusDeg=xMatchRadiusArcmin/60.
cat2=SkyCoord(ra = xTab['RADeg'].data, dec = xTab['decDeg'].data, unit = 'deg')
xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
matchTab=xTab[xIndices]

outTab=atpy.Table()
outTab.add_column(tab['name'])
outTab.add_column(tab['RADeg'])
outTab.add_column(tab['decDeg'])
outTab.add_column(atpy.Column(matchTab['z'].data, 'redshift'))
if 'zErr' in matchTab.keys():
    outTab.add_column(atpy.Column(matchTab['zErr'].data, 'redshiftErr'))
else:
    outTab.add_column(atpy.Column(np.zeros(len(tab)), 'redshiftErr'))

mask=np.less(rDeg.value, xMatchRadiusDeg)
outTab=outTab[mask]

outTab.write(args.output_catalog, overwrite = True)
