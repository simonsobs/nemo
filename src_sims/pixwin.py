import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from enlib import utils


## These are the official files from Marius
ell2, bl2 = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_tform_151107_v1_2014_pa2_jitmap_deep56.txt', usecols=[0,1], unpack=True)
ell, bl = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_tform_151107_v1_2014_pa2_jitter_deep56.txt', usecols=[0,1], unpack=True)

#Pixel-window function
wl     = np.sinc(0.5*ell/(2.*180.*60.))


plt.close()
plt.plot(ell, bl2/bl)
plt.plot(ell, wl)
plt.xscale('log')
plt.savefig('pixwin.png')
plt.yscale('log')
plt.savefig('pixwin2.png')

plt.close()
plt.plot(ell, (bl2/bl/wl)-1.)
plt.xscale('log')
plt.savefig('pixwin_ratio.png')

