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

# Cut to just the halos in the survey mask
cutTabFileName="halosInMask.fits"
if os.path.exists(cutTabFileName) == False:
    print("Cutting halos catalog to the survey mask")
    tab=atpy.Table().read('../halos.fits')
    checkMask=selFn.checkCoordsInAreaMask(tab['RADeg'], tab['decDeg'])
    tab=tab[checkMask]
    tab.write(cutTabFileName, overwrite = True)
print("Reading %s" % (cutTabFileName))
tab=atpy.Table().read(cutTabFileName)

# On-the-fly mass conversion as quick with CCL
if massCol == "M500c":
    print("Converting M200m to M500c")
    M500c=[]
    count=0
    M200mDef=ccl.halos.MassDef200m(c_m='Bhattacharya13')
    M500cDef=ccl.halos.MassDef(500, "critical")
    M500c=[]
    count=0
    for row in tab:
        M500c.append(M200mDef.translate_mass(selFn.mockSurvey.cosmoModel, row['M200m'], 1/(1+row['z']), M500cDef))
    tab['M500c']=M500c

# Bit of preprocessing to make life easier
tab['fixed_SNR']=100.0
tab.rename_column('z', 'redshift')
tab[massCol]=tab[massCol]/1e14

## Example (not used here) - N(z) with M500c > 5e13 MSun - with selection function applied
#predMz=selFn.compMz*selFn.mockSurvey.clusterCount
#countsByRedshift=predMz[:, np.greater(selFn.mockSurvey.log10M, np.log10(5e13))].sum(axis = 1)

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
        #complCorr[j]=selFn.compMz[zMask, :].mean(axis = 0)[mMask].mean()
    validMask=(obsCounts > 0)
    results[label]={'log10MBinCentres': log10MBinCentres[validMask],
                    'predCounts': predCounts[validMask],
                    'obsCounts': obsCounts[validMask],
                    'obsCountsErr': obsCountsErr[validMask],
                    'predNumDensity': predNumDensity[validMask],
                    'obsNumDensity': obsNumDensity[validMask],
                    'obsNumDensityErr': (obsCountsErr[validMask]/obsCounts[validMask])*obsNumDensity[validMask]}

# Counts comparison plot (just N as a function of mass) -----------------------------------------------------
plotSettings.update_rcParams()
plt.figure(figsize=(9,6.5))
ax=plt.axes([0.15, 0.12, 0.84, 0.85])
for key in results.keys():
    plotLog10MBinCentres=results[key]['log10MBinCentres']
    pred=results[key]['predCounts']
    obs=results[key]['obsCounts']
    obsErr=results[key]['obsCountsErr']
    plt.errorbar(plotLog10MBinCentres, obs, yerr = obsErr,
                elinewidth = 3, fmt = 'D', ms = 6, zorder = 900, label = key)
    plt.plot(plotLog10MBinCentres, pred, 'k-')
plt.semilogy()
plt.ylim(0.1, 5e5)
plt.xlim(14.0, log10MBinEdges.max())
plt.xlabel("log$_{10}$($M^{\\rm true}_{\\rm %s}$ / $M_{\odot}$)" % (deltaLabel))
plt.ylabel("$N$")
plt.legend()
plt.savefig("%s_counts.png" % (massCol))
plt.close()

# Counts per unit volume (N per Mpc^3) ----------------------------------------------------------------------
plotSettings.update_rcParams()
plt.figure(figsize=(9,6.5))
ax=plt.axes([0.15, 0.12, 0.84, 0.85])
for key in results.keys():
    plotLog10MBinCentres=results[key]['log10MBinCentres']
    pred=results[key]['predNumDensity']
    obs=results[key]['obsNumDensity']
    obsErr=results[key]['obsNumDensityErr']
    plt.errorbar(plotLog10MBinCentres, obs, yerr = obsErr,
                elinewidth = 3, fmt = 'D', ms = 6, zorder = 900, label = key)
    plt.plot(plotLog10MBinCentres, pred, 'k-')
plt.semilogy()
#plt.ylim(0.1, 5e5)
plt.xlim(14.0, log10MBinEdges.max())
plt.xlabel("log$_{10}$($M^{\\rm true}_{\\rm %s}$ / $M_{\odot}$)" % (deltaLabel))
plt.ylabel("$N$ (Mpc$^{3}$)")
plt.legend()
plt.savefig("%s_numDensity.png" % (massCol))
plt.close()
