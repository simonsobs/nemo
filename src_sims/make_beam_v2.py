import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np
from enlib import utils
from scipy import interpolate


## These are the official files from Marius
t, bt   = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_profile_151107_v1_2014_pa2_jitter_deep56.txt', unpack=True)

t2, bt2 = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_profile_151107_v1_2014_pa2_jitmap_deep56.txt', unpack=True)
ell, bl = np.loadtxt('/data/saiola/depot/Beams/151107_v1/beam_tform_151107_v1_2014_pa2_jitmap_deep56.txt', usecols=[0,1], unpack=True)


#Mulpoles from the rFFT of the real-space beam transform
#Let's try something that is a bit of a hack and increase the resolution in ell
dt      = (180./max(ell))/180.*np.pi #np.median(np.diff(t2))/180.*np.pi
t2n     = np.arange(0,50,dt/np.pi*180.)
ell_tmp = np.fft.rfftfreq(t2n.size,dt)*2.*np.pi
print max(ell_tmp), max(ell)

#Setting up the beam
#sigma  = 1.3*utils.fwhm*utils.arcmin
#bl     = np.exp(-(ell**2.)*(sigma**2.))
#bl     = np.exp(-(ell**2.)*(sigma**2.)/2.) #After Sigurd's fix


#Pixel-window function
#wl     = np.sinc(0.5*ell/(2.*180.*60.))
bl2     = bl*1.#wl


## We try to interpolate the logs to avoid to zeropad
f       = interpolate.splrep(np.log(ell[1:]), np.log(bl2[1:]), s=0)
bl2n    = interpolate.splev(np.log(ell_tmp[1:]), f, der=0)
bl2n    = np.exp(bl2n)
bl2n    = np.insert(bl2n, 0, 1)
bt2n    = np.fft.irfft(bl2n, n=t2n.size)
#bt2n   = bt2n/bt2n[0]

plt.close()
plt.plot(ell, bl)
plt.plot(ell_tmp, bl2n)
plt.xscale('log')
plt.savefig('window.png')
plt.yscale('log')
plt.savefig('window2.png')


plt.close()
plt.plot(t2, bt2)
plt.plot(t2n, bt2n)
plt.yscale('log')
plt.xlim(0,max(t))
plt.ylim(1e-7,1)
plt.savefig('window_test.png')
