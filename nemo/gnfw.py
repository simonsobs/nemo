"""
Routines to perform line-of-sight integrals of the GNFW model. 

This code was originally written by Matthew Hasselfield, and has been slightly extended.

It now includes much faster routines based on those in Pixell by Sigurd Naess
(but here we use Arnaud+2010 style notation)

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
    """The GNFW radial profile.
    
    Args:
        x (:obj:`np.ndarray`): Radial coordinate.
        params (:obj:`dict`): Dictionary with keys `alpha`, `beta`, `gamma`, `c500`, and `P0` that defines
            the GNFW profile shape.
    
    Returns:
        Profile (1d :obj:`np.ndarray`).
        
    """
    G, A, B, c500, P0 = params['gamma'], params['alpha'], params['beta'], params['c500'], params['P0']
    prof=np.zeros(x.shape)
    mask=np.greater(x, 0)
    prof[mask]=P0*((x[mask]*c500)**-G * (1+(x[mask]*c500)**A)**((G-B)/A))
    #prof[x == 0]=np.inf
    return prof

def xfunc(x, b, params):
    """The log-scaled integrand for GNFW cylindrical integration along the line-of-sight coordinate `x`,
    with impact parameter `b`.
    
    Args:
        x (:obj:`np.ndarray`): Line-of-sight coordinate.
        b (:obj:`float`): Impact parameter.
        params (:obj:`dict`): Dictionary with keys `alpha`, `beta`, `gamma`, `c500`, and `P0` that defines
            the GNFW profile shape.
    
    Returns:
        Value of integrand.
        
    """
    x = np.array(x)
    if x.ndim == 0:
        return xfunc([x], b, params)[0]
    r = (x**2 + b**2)**.5
    y = x*func(r, params)
    y[x==0] = 0
    return y

def integrated(b, params = _default_params):
    """Returns the line of sight integral of the GNFW profile at impact parameter `b`.
    
    Args:
        b (:obj:`float`): Impact parameter.
        params (:obj:`dict`): Dictionary with keys `alpha`, `beta`, `gamma`, `c500`, and `P0` that defines
            the GNFW profile shape.
    
    Returns:
        Line of sight integral at given impact parameter.
        
    """
    # Since the GNFW is a smoothly varying power-law from
    # 0+ to +infty, we make a change of variables to
    #     u = ln x
    # and perform the integral by summing over equally spaced logarithmic
    # bins in x.
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

#------------------------------------------------------------------------------------------------------------
# Faster routines, based on Sigurd's code in pixell - but here we use A10-style GNFW function notation
# Docstrings not updated yet!

def gnfw(x, c, alpha, beta, gamma, P0 = 1):
    """The GNFW radial profile (in the style of Hasselfield+2013, Arnaud+2010).

    Args:
        x (:obj:`np.ndarray`): Radial coordinate.
        c (:obj:`float`): Concentration parameter (e.g. c500).
        alpha (:obj:`float`): GNFW shape parameter (see Arnaud+2010).
        beta (:obj:`float`): GNFW shape parameter (see Arnaud+2010).
        gamma (:obj:`float`): GNFW shape parameter (see Arnaud+2010).
        P0 (:obj:`float`): Pressure normalization (see Arnaud+2010).

    Returns:
        Profile (1d :obj:`np.ndarray`).
    """

    return P0*((x*c)**-gamma * (1+(x*c)**alpha)**((gamma-beta)/alpha))

def tsz_profile_raw(x, c = 1.177, alpha = 1.0510, beta = 5.4905, gamma = 0.3081):
    return gnfw(x, c = c, alpha = alpha, beta = beta, gamma = gamma)

_tsz_profile_los_cache = {}
def tsz_profile_los(x, c = 1.177, alpha = 1.0510, beta = 5.4905, gamma = 0.3081, zmax=1e5, npoint=100, x1=1e-8, x2=1e4, _a=8, cache=None):
    """Fast, highly accurate approximate version of tsz_profile_los_exact. Interpolates the exact
    function in log-log space, and caches the interpolator. With the default settings,
    it's accurate to better than 1e-5 up to at least x = 10000, and building the
    interpolator takes about 25 ms. After that, each evaluation takes 50-100 ns per
    data point. This makes it about 10000x faster than tsz_profile_los_exact.
    See tsz_profile_raw for the units."""
    from scipy import interpolate
    # Cache the fit parameters.
    if cache is None: global _tsz_profile_los_cache
    else: _tsz_profile_los_cache = {}
    key = (c, alpha, beta, gamma, zmax, npoint, _a, x1, x2)
    if key not in _tsz_profile_los_cache:
        xp = np.linspace(np.log(x1),np.log(x2),npoint)
        yp = np.log(tsz_profile_los_exact(np.exp(xp), c = c, alpha=alpha, beta=beta, gamma=gamma, zmax=zmax, _a=_a))
        _tsz_profile_los_cache[key] = (interpolate.interp1d(xp, yp, "cubic"), x1, x2, yp[0], yp[-1], (yp[-2]-yp[-1])/(xp[-2]-xp[-1]))
    spline, xmin, xmax, vleft, vright, slope = _tsz_profile_los_cache[key]
    # Split into 3 cases: x<xmin, x inside and x > xmax.
    x     = np.asfarray(x)
    left  = x<xmin
    right = x>xmax
    inside= (~left) & (~right)
    return np.piecewise(x, [inside, left, right], [
        lambda x: np.exp(spline(np.log(x))),
        lambda x: np.exp(vleft),
        lambda x: np.exp(vright + (np.log(x)-np.log(xmax))*slope),
    ])

def tsz_profile_los_exact(x, c = 1.177, alpha = 1.0510, beta = 5.4905, gamma = 0.3081, zmax=1e5, _a=8):
    """Line-of-sight integral of the cluster_pressure_profile. See tsz_profile_raw
    for the meaning of the arguments. Slow due to the use
    of quad and the manual looping this requires. Takes about 1 ms per data point.
    The argument _a controls a change of variable used to improve the speed and
    accuracy of the integral, and should not be necessary to change from the default
    value of 8.

    See tsz_profile_raw for the units and how to scale it to something physical.
    Without scaling, the profile has a peak of about 0.5 and a FWHM of about
    0.12 with the default parameters.

    Instead of using this function directly, consider using
    tsz_profile_los instead. It's 10000x faster and more than accurate enough.
    """
    from scipy.integrate import quad
    x     = np.asarray(x)
    xflat = x.reshape(-1)
    # We have int f(x) dx, but would be easier to integrate
    # int y**a f(y) dy. So we want y**a dy = dx => 1/(a+1)*y**(a+1) = x
    # => y = (x*(a+1))**(1/(a+1))
    def yfun(x): return (x*(_a+1))**(1/(_a+1))
    def xfun(y): return y**(_a+1)/(_a+1)
    res    = 2*np.array([quad(lambda y: y**_a*tsz_profile_raw((xfun(y)**2+x1**2)**0.5, c=c, alpha=alpha, beta=beta, gamma=gamma), 0, yfun(zmax))[0] for x1 in xflat])
    res   = res.reshape(x.shape)
    return res

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
