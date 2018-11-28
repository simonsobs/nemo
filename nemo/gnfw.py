"""
Perform line-of-sight integral of GNFW model.

Matthew Hasselfield's code, slightly updated to include P0, c500
(makes life a bit easier for some applications of this).

"""

import numpy as np
from scipy.optimize import fmin
#import warnings
#warnings.filterwarnings("error")

# Below, params is a dictionary that specifies the GNFW powers as well
# as the integration control parameters.  The GNFW parameters should
# be keyed with 'alpha', 'beta', 'gamma'.
#
# The integration control parameters are:
#   'tol'   - the approximate precision in the answer.
#   'npts'  - number of points over which to evaulate the Riemann sum;
#             a couple hundred is enough, unless tol is smaller than 1e-6.
#
# Even at tol=1e-7, the difference between npts~100 and npts~1000 is
# less than 1e-8.

_default_params = {
        'P0': 8.403,
        'c500': 1.177,
        'gamma': 0.3081,
        'alpha': 1.0510,
        'beta':  5.4905,
        'tol': 1e-7,
        'npts': 100,
        }

def func(x, params):
    """
    The GNFW radial profile - now with added P0, c500.
    """
    G, A, B, c500, P0 = params['gamma'], params['alpha'], params['beta'], params['c500'], params['P0']
    prof=np.zeros(x.shape)
    mask=np.greater(x, 0)
    prof[mask]=P0*((x[mask]*c500)**-G * (1+(x[mask]*c500)**A)**((G-B)/A))
    #prof[x == 0]=np.inf
    return prof

def xfunc(x, b, params):
    """
    The log-scaled integrand for GNFW cylindrical integration along
    line-of-sight coordinate x, with impact parameter b.
    """
    x = np.array(x)
    if x.ndim == 0:
        return xfunc([x], b, params)[0]
    r = (x**2 + b**2)**.5
    y = x*func(r, params)
    y[x==0] = 0
    return y

def integrated(b, params = _default_params):
    """
    Return the line of sight integral of the GNFW profile at impact
    parameter b.  Since the GNFW is a smoothly varying power-law from
    0+ to +infty, we make a change of variables to
        u = ln x
    and perform the integral by summing over equally spaced logarithmic
    bins in x.
    """
    # Get GNFW vars as well as integration parameters
    G, A, B = params['gamma'], params['alpha'], params['beta']
    TH, N = params.get('tol',1e-6), params.get('npts', 200)
    # The function ( x * y(r) ) has a local maximum around b or so
    x_max = fmin(lambda x: -xfunc(x,b,params), b, disp=0)[0]
    y_max = xfunc(x_max, b, params)
    # The wings of x * y(r) fall off exponentially (in log x).  This makes
    # the truncation error easy to estimate (or even to include)
    x_lo = (y_max * TH)**(1/(1-G))
    x_hi = (y_max * TH)**(1/(1-B))
    # Take log-spaced bins
    u_lo, u_hi = np.log(x_lo), np.log(x_hi)
    du = (u_hi-u_lo) / N
    x = np.exp(np.arange(u_lo, u_hi, du))
    # Sum
    I1 = np.sum(du*xfunc(x,b,params))
    # Wing (under-)estimate
    x_hi = np.exp(u_hi)
    I2 = x_lo**(1-G)/(1-G) + x_hi**(1-B)/(1-B)
    return I1 + I2

# Test
#if __name__ == '__main__':
    #from pylab import semilogy, show, plot, subplot
    ## Use default GNFW params
    #params = _default_params.copy()
    ## Test a few different integration parameters
    #subplot(211)
    #bs = np.arange(0., 2., 1e-2)
    #zzs = []
    #for N in [1e2, 1e3]:
        #params['npts'] = N
        #for tol in [ 1e-7]:
            #params['tol'] = tol
            #zs = [integrated(b, params) for b in bs]
            #semilogy(bs, zs)
            #zzs.append(zs)
    # Show fractional difference of each curve from the last curve
    #zzs = np.array(zzs)
    #subplot(212)
    #for z in zzs[:-1]:
        #d = 1 - z/zzs[-1]
        #plot(bs, d)
    #show()
    
    #ipshell()
    #sys.exit()
