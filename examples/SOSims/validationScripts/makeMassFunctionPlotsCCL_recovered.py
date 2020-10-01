"""

Plot the mass function in z bins.

Range adjusted to drop the last bin, which is more incomplete in the sense that it may not cover that
full mass bin (whereas all other bins are guaranteed to by definition).

"""

import os
import sys
import astropy.table as atpy
import astropy.io.fits as pyfits
import IPython
import numpy as np
from nemo import plotSettings, completeness, signals
import pylab as plt
from scipy import stats
from astLib import *
import pyccl as ccl
from colossus.lss import mass_function

#------------------------------------------------------------------------------------------------------------
# Options
SNRCut=4.0
selFnDir="../MFMF_SOSim_3freq_tiles/selFn"
footprintLabel=None
massCol='M200m'
zBinEdges=[0.2, 0.5, 0.9, 1.2]
zMin=min(zBinEdges)
zMax=max(zBinEdges)
log10MBinEdges=np.linspace(13.8, 15.5, 18)
symbs=['D', 's', 'o']
# Handling different mass definitions
if massCol == 'M500c':
    delta=500
    rhoType="critical"
elif massCol == 'M200m':
    delta=200
    rhoType="matter"
else:
    raise Exception("Unsupported massCol - should be M500c or M200m")
deltaLabel="%d%s" % (delta, rhoType[0])
log10MBinCentres=(log10MBinEdges[1:]+log10MBinEdges[:-1])/2

# Set up Websky cosmology
H0, Om0, Ob0, sigma_8, ns = 68.0, 0.31, 0.049, 0.81, 0.965
selFn=completeness.SelFn(selFnDir, SNRCut, footprintLabel = footprintLabel, zStep = 0.02,
                         delta = delta, rhoType = rhoType)
scalingRelationDict=selFn.scalingRelationDict
selFn.update(H0, Om0, Ob0, sigma_8, ns, scalingRelationDict = scalingRelationDict)
print("Total area = %.3f square degrees" % (selFn.totalAreaDeg2))

# Load Nemo catalog
tab=atpy.Table().read('../MFMF_SOSim_3freq_tiles/MFMF_SOSim_3freq_tiles_M500.fits')
tab.rename_column("M500", "M500c")

# All the analysis first ------------------------------------------------------------------------------------
# WARNING: We're using halo catalogs, so disabled completeness correction
results={}
predMz=selFn.mockSurvey.clusterCount
for i in range(len(zBinEdges)-1):
    zMin=zBinEdges[i]
    zMax=zBinEdges[i+1]
    label='%.1f < z < %.1f' % (zMin, zMax)
    fSky=selFn.mockSurvey.areaDeg2/(4*np.pi*(180/np.pi)**2)
    shellVolumeMpc3=fSky*(selFn.mockSurvey._comovingVolume(zMax)-selFn.mockSurvey._comovingVolume(zMin))
    zMask=np.logical_and(selFn.mockSurvey.z >= zMin, selFn.mockSurvey.z < zMax)
    countsByMass=predMz[zMask, :].sum(axis = 0)
    
    predCounts=np.zeros(len(log10MBinEdges)-1)
    predNumDensity=np.zeros(len(log10MBinEdges)-1)
    obsCounts=np.zeros(len(log10MBinEdges)-1)
    obsCountsErr=np.zeros(len(log10MBinEdges)-1)
    obsNumDensity=np.zeros(len(log10MBinEdges)-1)
    obsNumDensityErr=np.zeros(len(log10MBinEdges)-1)
    complCorr=np.zeros(len(log10MBinEdges)-1) # Holds average completeness in each mass bin
    
    h=H0/100.
    binTab=tab[np.logical_and(tab['redshift'] >= zMin, tab['redshift'] < zMax)]
    obsLog10Ms=np.log10(binTab[massCol]*1e14)
    for j in range(len(log10MBinEdges)-1):
        mMin=log10MBinEdges[j]
        mMax=log10MBinEdges[j+1]
        mMask=np.logical_and(selFn.mockSurvey.log10M >= mMin, selFn.mockSurvey.log10M < mMax)
        predCounts[j]=countsByMass[mMask].sum()
        obsMask=np.logical_and(obsLog10Ms >= mMin, obsLog10Ms < mMax)
        obsCounts[j]=obsMask.sum()
        obsCountsErr[j]=np.sqrt(obsCounts[j])
        predNumDensity[j]=predCounts[j]/shellVolumeMpc3
        obsNumDensity[j]=obsCounts[j]/shellVolumeMpc3
        complCorr[j]=selFn.compMz[zMask, :].mean(axis = 0)[mMask].mean()
    validMask=(obsCounts > 0)
    fracErr=obsCountsErr[validMask]/obsCounts[validMask]
    results[label]={'log10MBinCentres': log10MBinCentres[validMask],
                    'predCounts': predCounts[validMask],
                    'obsCounts': obsCounts[validMask],
                    'obsCountsErr': obsCountsErr[validMask],
                    'predNumDensity': predNumDensity[validMask],
                    'obsNumDensity': obsNumDensity[validMask],
                    'obsNumDensityErr': fracErr*obsNumDensity[validMask],
                    # Completeness corrected
                    'corr_obsCounts': obsCounts[validMask]/complCorr[validMask],
                    'corr_obsCountsErr': fracErr*(obsCounts[validMask]/complCorr[validMask]),
                    'corr_obsNumDensity': obsNumDensity[validMask]/complCorr[validMask],
                    'corr_obsNumDensityErr': fracErr*(obsNumDensity[validMask]/complCorr[validMask]),      
                    }

# Counts comparison plot (just N as a function of mass) -----------------------------------------------------
plotSettings.update_rcParams()
plt.figure(figsize=(9,6.5))
ax=plt.axes([0.15, 0.12, 0.84, 0.85])
for key, symb in zip(results.keys(), symbs):
    plotLog10MBinCentres=results[key]['log10MBinCentres']
    pred=results[key]['predCounts']
    obs=results[key]['obsCounts']
    obsErr=results[key]['obsCountsErr']
    corr_obs=results[key]['corr_obsCounts']
    corr_obsErr=results[key]['corr_obsCountsErr']
    plt.errorbar(plotLog10MBinCentres, obs, yerr = obsErr, color = 'none', markeredgecolor = 'k',
                elinewidth = 3, fmt = symb, ms = 6, zorder = 900)
    plt.errorbar(plotLog10MBinCentres, corr_obs, yerr = corr_obsErr,
                elinewidth = 3, fmt = symb, ms = 6, zorder = 900, label = key)
    plt.plot(plotLog10MBinCentres, pred, 'k-')
plt.semilogy()
plt.ylim(0.1, 5e5)
plt.xlim(14.0, log10MBinEdges.max())
plt.xlabel("log$_{10}$($M_{\\rm %s}$ / $M_{\odot}$)" % (deltaLabel))
plt.ylabel("$N$")
plt.legend()
plt.savefig("Recovered_%s_counts.png" % (massCol))
plt.close()

# Counts per unit volume (N per Mpc^3) ----------------------------------------------------------------------
plotSettings.update_rcParams()
plt.figure(figsize=(9,6.5))
ax=plt.axes([0.15, 0.12, 0.84, 0.85])
for key, symb in zip(results.keys(), symbs):
    plotLog10MBinCentres=results[key]['log10MBinCentres']
    pred=results[key]['predNumDensity']
    obs=results[key]['obsNumDensity']
    obsErr=results[key]['obsNumDensityErr']
    corr_obs=results[key]['corr_obsNumDensity']
    corr_obsErr=results[key]['corr_obsNumDensityErr']
    plt.errorbar(plotLog10MBinCentres, obs, yerr = obsErr, color = 'none', markeredgecolor = 'k',
                elinewidth = 3, fmt = symb, ms = 6, zorder = 900)
    plt.errorbar(plotLog10MBinCentres, corr_obs, yerr = corr_obsErr,
                elinewidth = 3, fmt = symb, ms = 6, zorder = 900, label = key)
    plt.plot(plotLog10MBinCentres, pred, 'k-')
plt.semilogy()
#plt.ylim(0.1, 5e5)
plt.xlim(14.0, log10MBinEdges.max())
plt.xlabel("log$_{10}$($M_{\\rm %s}$ / $M_{\odot}$)" % (deltaLabel))
plt.ylabel("$N$ (Mpc$^{-3}$)")
plt.legend()
plt.savefig("Recovered_%s_numDensity.png" % (massCol))
plt.close()

IPython.embed()
sys.exit()
