import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from enlib import utils
from scipy.special import legendre

## These are the official files from Marius
t, bt   = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_profile_151107_v1_2014_pa2_jitter_deep56.txt', unpack=True)

t2, bt2 = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_profile_151107_v1_2014_pa2_jitmap_deep56.txt', unpack=True)
ell, bl = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_tform_151107_v1_2014_pa2_jitmap_deep56.txt', usecols=[0,1], unpack=True)


## We need to apply the Legendre transformation
bt2n = bt2*0.

ct = np.cos(t2/180.*np.pi)
for l in ell:
   print l
   p = legendre(l)
   bt2n += bl[l]*(2.*l+1.)/4./np.pi*p(ct)

plt.close()
plt.plot(t2, bt2)
plt.plot(t2, bt2n)
plt.yscale('log')
plt.xlim(0,max(t2))
plt.ylim(1e-7,1)
plt.savefig('window_test_legendre.png')
