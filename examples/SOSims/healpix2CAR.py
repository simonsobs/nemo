"""

Based on the tile-c script by Mat M.

Here we just want plain-old healpix -> CAR projection, no Galactic coord transform involved,
and no polarization either.

We've added an option to add extra white noise - this is a quick/easy (but suboptimal) way
to avoid ringing in filtered maps (may go away when we move to NSIDE = 8192 sims).

"""

import argparse, os, sys, time
import astropy.io.fits as pyfits
from astropy.wcs import WCS
import IPython
parser = argparse.ArgumentParser()
parser.add_argument("ifiles", nargs="+")
parser.add_argument("template")
parser.add_argument("-L", '--label', type=str, default="CAR") 
parser.add_argument('-W', '--add-white-noise', type=float, default=0)
parser.add_argument("-l", "--lmax",    type=int, default=0)
parser.add_argument("-v", "--verbose", action="store_true")
parser.add_argument("-O", "--outputs", type=str, default="map,ivar")
parser.add_argument('-n', '--npol', type=int, default=1)
args = parser.parse_args()
import numpy as np, healpy
from pixell import enmap, utils, curvedsky, coordinates, mpi, memory

unit  = 1
#euler = np.array([57.06793215,  62.87115487, -167.14056929])*utils.degree
dtype = np.float32
ctype = np.result_type(dtype,0j)
rstep = 100
nside = 0
npol = 1
verbose = args.verbose
outputs = args.outputs.split(",")
if len(outputs) == 0:
    print("No outputs selected - nothing to do")
    sys.exit(0)
t0 = time.time()

def progress(msg):
    if verbose:
        print("%6.2f %6.2f %6.2f %s" % ((time.time()-t0)/60, memory.current()/1024.**3, memory.max()/1024.**3, msg))

comm = mpi.COMM_WORLD

header=pyfits.header.Header.fromtextfile(args.template)
wcs=WCS(header)
shape=(header['NAXIS2'], header['NAXIS1'])
#shape = shape[-2:]
progress("Allocating output map %s %s" % (str((1,)+shape), str(dtype)))
if "map" in outputs:
    for ifile in args.ifiles[comm.rank::comm.size]:
        name = os.path.basename(ifile)
        runit = unit
        #if "545" in name: 
        #	runit *= factor_545
        #	npol = 1
        #elif "857" in name: 
        #	runit *= factor_857
        #	npol = 1
        #else:
        #	npol = 3
        npol=args.npol
        fields = range(npol)
        omap = enmap.zeros((npol,)+shape, wcs, dtype)
        progress("%s TQU read" % name)
        imap  = np.array(healpy.read_map(ifile, fields)).astype(dtype)
        nside = healpy.npix2nside(imap.shape[-1])
        progress("%s TQU mask" % name)
        imap[healpy.mask_bad(imap)] = 0
        progress("%s TQU scale" % name)
        imap *= runit
        nside = healpy.npix2nside(imap.shape[-1])
        lmax  = args.lmax or 3*nside
        progress("%s TQU alm2map" % name)
        alm   = curvedsky.map2alm_healpix(imap, lmax=lmax)
        del imap
        # work around healpix bug
        #progress("%s TQU rotate_alm" % name)
        #alm   = alm.astype(np.complex128,copy=False)
        #healpy.rotate_alm(alm, euler[0], euler[1], euler[2])
        #alm   = alm.astype(ctype,copy=False)
        progress("%s TQU map2alm" % name)
        curvedsky.alm2map_cyl(alm, omap)
        del alm
        ofile = ifile.replace(".fits", "_%s.fits" % (args.label))
        if args.add_white_noise > 0:
            progress("%s adding extra white noise with sigma = %.1f per pixel" % (name, args.add_white_noise))
            omap[0]=omap[0]+np.random.normal(0, args.add_white_noise, omap[0].shape)
        progress("%s TQU write %s" % (name, ofile))
        enmap.write_map(ofile, omap[0])
