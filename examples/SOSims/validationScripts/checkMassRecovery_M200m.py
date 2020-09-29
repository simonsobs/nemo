"""

Fit the scaling relation in the sims

"""

import os
import sys
import numpy as np
import astropy.table as atpy
from nemo import catalogs, signals, plotSettings, MockSurvey
from astropy.cosmology import FlatLambdaCDM
from scipy import stats
import pylab as plt
import IPython

#------------------------------------------------------------------------------------------------------------
def calcMass(tab, massOptions, tckQFitDict, fRelWeightsDict, mockSurvey):
    """Calculates masses for cluster data in table. Because why not re-do on the fly when hippo busy?
    
    """
    
    label=mockSurvey.mdefLabel
    
    count=0
    for row in tab:
        count=count+1
        #print("... %d/%d; %s (%.3f +/- %.3f) ..." % (count, len(tab), row['name'], 
                                                              #row['redshift'], row['redshiftErr']))

        tileName=row['tileName']
        
        # Cuts on z, fixed_y_c for forced photometry mode (invalid objects will be listed but without a mass)
        if row['fixed_y_c'] > 0 and np.isnan(row['redshift']) == False:
            # Corrected for mass function steepness
            massDict=signals.calcMass(row['fixed_y_c']*1e-4, row['fixed_err_y_c']*1e-4, 
                                            row['redshift'], row['redshiftErr'],
                                            tenToA0 = massOptions['tenToA0'],
                                            B0 = massOptions['B0'], 
                                            Mpivot = massOptions['Mpivot'], 
                                            sigma_int = massOptions['sigma_int'],
                                            tckQFit = tckQFitDict[tileName], mockSurvey = mockSurvey, 
                                            applyMFDebiasCorrection = True,
                                            applyRelativisticCorrection = True,
                                            fRelWeightsDict = fRelWeightsDict[tileName])
            row['%s' % (label)]=massDict['%s' % (label)]
            row['%s_errPlus' % (label)]=massDict['%s_errPlus' % (label)]
            row['%s_errMinus' % (label)]=massDict['%s_errMinus' % (label)]

    return tab

#------------------------------------------------------------------------------------------------------------
# Main

# Websky cosmo - for on-the-fly redone masses
minMass=1e13
areaDeg2=700.   # doesn't matter
zMin=0.0
zMax=2.0
H0, Om0, Ob0, sigma8, ns = 68.0, 0.049+0.261, 0.049, 0.81, 0.965
TCMB=2.72548
cosmoModel=FlatLambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0, Tcmb0 = TCMB)
mockSurvey=MockSurvey.MockSurvey(minMass, areaDeg2, zMin, zMax, H0, Om0, Ob0, sigma8, ns,
                                 rhoType = 'matter', delta = 200)
massOptions={'tenToA0': 1.7e-05,
             'B0': 0.0,
             'Mpivot': 3.0e+14, 
             'sigma_int': 0.0}
tckQFitDict=signals.loadQ("../MFMF_SOSim_3freq_tiles/selFn/QFit.fits")
fRelWeightsDict=signals.loadFRelWeights("../MFMF_SOSim_3freq_tiles/selFn/fRelWeights.fits")

# Make combined table
mergedTabFileName="trueMasses_MFMF_SOSim_3freq_tiles_M500.fits"
if os.path.exists(mergedTabFileName) == False:
    halos=atpy.Table().read("../halos.fits")
    tab=atpy.Table().read("../MFMF_SOSim_3freq_tiles/MFMF_SOSim_3freq_tiles_M500.fits")
    tab=tab[tab['fixed_SNR'] > 6]

    tab, halos, rDeg=catalogs.crossMatch(tab, halos, radiusArcmin = 1.0)

    zs=halos['z']
    yc=tab['fixed_y_c']
    M200m=halos['M200m']
    M500c=[]
    count=0
    for m, z in zip(M200m, zs):
        count=count+1
        print(count, len(M200m))
        M500c.append(signals.convertM200mToM500c(m, z))
    M500c=np.array(M500c)
    M500c=M500c[:, 0]
    tab['true_M500']=M500c/1e14
    tab['true_M200']=M200m/1e14
    tab['redshift']=zs
    tab.write(mergedTabFileName, overwrite = True)

# Re-do masses on the fly
tab=atpy.Table().read(mergedTabFileName)

# Cut on mass and z to do the fit
MMin=3.0
zBinEdges=[0.2, 0.4, 0.6, 0.8, 1.0]
for i in range(len(zBinEdges)-1):
    zMin=zBinEdges[i]
    zMax=zBinEdges[i+1]

    fitTab=tab[tab['M200m'] > MMin]
    fitTab=fitTab[fitTab['redshift'] > zMin]
    fitTab=fitTab[fitTab['redshift'] < zMax]

    # NOTE: This is done in place anyway
    fitTab=calcMass(fitTab, massOptions, tckQFitDict, fRelWeightsDict, mockSurvey)

    y=fitTab['M200m']
    x=fitTab['true_M200']
    result=stats.linregress(x, y)
    sumSqRes=np.sum((x-y)**2)
    calibFactor=np.mean(fitTab['true_M200'])/np.mean(fitTab['M200m'])
    
    # Scaling relation plot
    plotSettings.update_rcParams()
    plt.figure(figsize=(9.5,9))
    ax=plt.axes([0.1, 0.1, 0.89, 0.89])
    ax.set_aspect('equal')        
    plotRange=np.linspace(1.0, 50.0, 100)
    plt.plot(x, y, '.')
    plt.plot(plotRange, plotRange, 'k-')
    plt.xlabel("$M^{\\rm true}_{\\rm 200m}$ (10$^{14}$ $M_{\odot}$)")
    plt.ylabel("$M_{\\rm 200m}$ (10$^{14}$ $M_{\odot}$)")
    plt.xlim(2, 50)
    plt.ylim(2, 50)
    plt.loglog()
    plt.title("%.1f < z < %.1f" % (zMin, zMax))
    plt.savefig("massRecoveryM200m_%.1f_%.1f.png" % (zMin, zMax))
    plt.close()

    print("%.1f < z < %.1f:" % (zMin, zMax))
    print("    calibFactor = ", calibFactor)
    print("    sumSqRes = ", sumSqRes)
