"""

This module contains global plot settings. 

For any routine that makes a plot using matplotlib, call plotSettings.update_rcParams() first.

"""

import pylab as plt
import matplotlib as mpl
from cycler import cycler

#------------------------------------------------------------------------------------------------------------
def update_rcParams(dict={}):
    """Based on Cristobal's preferred settings. Updates matplotlib rcParams in place.
    
    """
    default = {}
    for tick in ('xtick', 'ytick'):
        default['{0}.major.size'.format(tick)] = 8
        default['{0}.minor.size'.format(tick)] = 4
        default['{0}.major.width'.format(tick)] = 2
        default['{0}.minor.width'.format(tick)] = 2
        default['{0}.labelsize'.format(tick)] = 20
        default['{0}.direction'.format(tick)] = 'in'
    default['xtick.top'] = True
    default['ytick.right'] = True
    default['axes.linewidth'] = 2
    default['axes.labelsize'] = 22
    default['font.size'] = 22
    default['font.family']='sans-serif'
    default['legend.fontsize'] = 18
    default['lines.linewidth'] = 2

    for key in default:
        plt.rcParams[key] = default[key]
    # if any parameters are specified, overwrite anything previously
    # defined
    for key in dict:
        plt.rcParams[key] = dict[key]

    # From https://github.com/mhasself/rg_friendly
    plt.rcParams['axes.prop_cycle']=cycler(color=['#2424f0','#df6f0e','#3cc03c','#d62728','#b467bd','#ac866b','#e397d9','#9f9f9f','#ecdd72','#77becf'])

# Let's just do this whenever it's imported
update_rcParams()
