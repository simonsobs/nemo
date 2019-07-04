"""

Combines components from SO sims into maps for each channel of interest.

This is first done in healpix - CAR projection is done separately.

"""

import os
import sys
import healpy
import numpy as np
import IPython

bands=['la145', 'la093', 'la225', 'la280']
bands.sort()
healpixSimsDir="../SO_mbs"
refBand='la145' # file names in components below must contain this
components={'cmb': "201905_extragalactic/4096/cmb/0000/simonsobs_cmb_uKCMB_la145_nside4096_0000.fits",
            'tsz': "201905_extragalactic/4096/tsz/0000/simonsobs_tsz_uKCMB_la145_nside4096_0000.fits",
            'cib': "201905_extragalactic/4096/cib/0000/simonsobs_cib_uKCMB_la145_nside4096_0000.fits",
            'noise': "201906_noise_no_lowell/4096/noise/0000/simonsobs_noise_uKCMB_la145_nside4096_0000.fits"}

# At the moment, including the noise can add 'diamonds' pattern (from healpix?) when CAR projected?
componentsToUse=['cmb', 'tsz', 'cib']
componentsToUse=componentsToUse+['noise']

for b in bands:
    print("Making combined map for %s" % (b))
    combinedMap=None
    for c in componentsToUse:
        fileName=healpixSimsDir+os.path.sep+components[c].replace(refBand, b)
        print("   reading map %s" % (fileName))
        if combinedMap is None:
            combinedMap=healpy.read_map(fileName)
        else:
            combinedMap=combinedMap+healpy.read_map(fileName)
    outFileName='TOnly_%s.fits' % (b)
    healpy.write_map(outFileName, combinedMap, overwrite = True) 
