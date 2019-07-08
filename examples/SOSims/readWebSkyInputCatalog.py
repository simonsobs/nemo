"""

Reads in the halos.pksc file and outputs .fits and .reg format catalogs for 
the 1 million most massive halos (M200m > 1e14 MSun). Cuts down to SO dec 
range.

"""
import os
import sys
import healpy as hp
import astropy.table as atpy
from nemo import catalogs
from cosmology import *
import numpy as np
import IPython

rho = 2.775e11*omegam*h**2 # Msun/Mpc^3

f=open('../WebSky/halos.pksc', 'rb')
N=np.fromfile(f,count=3,dtype=np.int32)[0]

# only take first five entries for testing (there are ~8e8 halos total...)
# comment the following line to read in all halos
N = 1000000

catalog=np.fromfile(f,count=N*10,dtype=np.float32)
catalog=np.reshape(catalog,(N,10))

x  = catalog[:,0];  y = catalog[:,1];  z = catalog[:,2] # Mpc (comoving)
vx = catalog[:,3]; vy = catalog[:,4]; vz = catalog[:,5] # km/sec
R  = catalog[:,6] # Mpc

# convert to mass, comoving distance, radial velocity, redshfit, RA and DEc
M        = 4*np.pi/3.*rho*R**3        # Msun (M200m)
chi      = np.sqrt(x**2+y**2+z**2)    # Mpc
vrad     = (x*vx + y*vy + z*vz) / chi # km/sec
redshift = zofchi(chi)      

theta, phi  = hp.vec2ang(np.column_stack((x,y,z))) # in radians

decDeg=-1*(np.degrees(theta)-90) # Because HEALPix
RADeg=np.degrees(phi)

names=[]
for ra, dec in zip(RADeg, decDeg):
    names.append(catalogs.makeACTName(ra, dec, prefix= 'MOCK-CL'))

outFileName="halos.fits"
tab=atpy.Table()
tab.add_column(atpy.Column(names, 'name'))
tab.add_column(atpy.Column(RADeg, 'RADeg'))
tab.add_column(atpy.Column(decDeg, 'decDeg'))
tab.add_column(atpy.Column(M, 'M200m'))
tab.add_column(atpy.Column(redshift, "z"))
tab=tab[np.where(tab['decDeg'] < 30)]
tab=tab[np.where(tab['decDeg'] > -70)]
tab.write(outFileName, overwrite = True)

catalogs.catalog2DS9(tab, "halos.reg")

#IPython.embed()
#sys.exit()

### e.g. project to a map, matching the websky orientations
#nside = 1024
#map   = np.zeros((hp.nside2npix(nside)))

#pix = hp.vec2pix(nside, x, y, z)
#pix = hp.ang2pix(nside, theta, phi) does the same

#weight = 1. #1 for number density, array of size(x) for arbitrary
#np.add.at(map, pix, weight)
