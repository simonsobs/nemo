import numpy as np, argparse
from enlib import enmap, utils

enmap.extent_model.append('subgrid')

parser = argparse.ArgumentParser()
parser.add_argument("ifile")
parser.add_argument("ofile")
parser.add_argument("-N",             type=float, default=0.06)
parser.add_argument("-b", "--beta",   type=float, default=-2.33)
parser.add_argument("-a", "--alim",   type=float, default=1e-5)
parser.add_argument("-B", "--beam",   type=float, default=1.3)
parser.add_argument("-n", "--numsrc", type=int,   default=0)
args = parser.parse_args()

N     = args.N
beta  = args.beta
alim  = args.alim
sigma = args.beam*utils.fwhm*utils.arcmin
ntot  = args.numsrc

map   = enmap.read_map(args.ifile)
A     = map.area()/(utils.degree**2.0)

# Source number density is N/sqdeg/mK * (amp/1mK)**beta.
# Can first draw total number of sources in patch, then
# draw each source's amplitude and position.
# Total number is -N/(beta+1) * a0**(beta+1), and is
# poisson distributed, with that number being lambda.
lam  = -N*A/(beta+1)*alim**(beta+1)
if (not ntot):
   ntot = np.random.poisson(lam)

# Each source amplitude is power law distributed between
# a0 and infinity. P(a>a') = (a'/alim)**(beta+1),
# which has inverse a' = P**(1/(beta+1))*alim
amps = np.random.uniform(0,1,ntot)**(1/(beta+1))*alim
print np.min(amps),  np.max(amps), ntot

# Flat sky approximation
pix  = np.random.randint(0,map.size,ntot)

map[:] = np.bincount(pix, amps, minlength=map.size).reshape(map.shape)
# Apply beam
map = enmap.smooth_gauss(map, sigma) * (2*np.pi*sigma**2)/map.pixsize()
map = map*1.e3
enmap.write_map(args.ofile, map)

pix_x = np.mod(pix,map.shape[1])
pix_y = (pix - pix_x)/map.shape[1]
pix2 = np.array([pix_y,pix_x])

coord = enmap.pix2sky(map.shape,map.wcs,pix2)/utils.degree
to_write = np.array([coord[1,:],coord[0,:],amps*1.e3]).T
np.savetxt('Simulated_catalog.txt', to_write,delimiter='	')

