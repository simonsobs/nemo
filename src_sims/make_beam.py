import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import numpy as np
from enlib import utils

t, bt = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_profile_151107_v1_2014_pa2_jitter_deep56.txt', unpack=True)
ell, bl = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_tform_151107_v1_2014_pa2_jitter_deep56.txt', usecols=[0,1], unpack=True)

sigma = 1.3*utils.fwhm*utils.arcmin
bl2 = np.exp(-(ell**2.)*(sigma**2.))
#bl2 = np.exp(-(ell**2.)*(sigma**2.)/2.) #After Sigurd's fix

bt2 = np.exp(-((t*utils.degree)**2.)/(4.*(sigma**2.)))
#bt2 = np.exp(-((t*utils.degree)**2.)/(2.*(sigma**2.))) #After Sigurd's fix

plt.close()
plt.plot(t,bt)
plt.plot(t,bt2)
plt.savefig('profile.png')

plt.close()
plt.plot(ell,bl)
plt.plot(ell,bl2)
#plt.yscale('log')
plt.savefig('tform.png')

np.savetxt('profile_gauss_1.3arcmin_v2.txt', np.array([t,bt2]).T, fmt=['%0.6f','%1.5e'], delimiter='  ')
np.savetxt('tform_gauss_1.3arcmin_v2.txt', np.array([ell,bl2]).T, fmt=['%5d','%1.8e'], delimiter='  ')
