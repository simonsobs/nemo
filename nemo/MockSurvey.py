"""Playing with the halo mass function... updated for the latest hmf which uses astropy


"""

import os
import sys
import numpy as np
import IPython
import atpy
import pylab as plt
import subprocess
import hmf
import cPickle
import astropy.cosmology as cosmology
from scipy import interpolate
from astLib import *
plt.matplotlib.interactive(True)

class MockSurvey(object):
    
    def __init__(self, minMass, areaDeg2, zMin, zMax, H0, Om0, Ob0, sigma_8, enableDrawSample = True):
        """Initialise a MockSurvey object. This first calculates the probability of drawing a cluster of 
        given M500, z, assuming the Tinker mass function, and the given (generous) selection limits. 
        An additional selection function can be dialled in later when using drawSample.
        
        NOTE: We've hard coded everything to use M500 wrt critical density at this point.
        
        Set enableDrawSample == True to be able to use drawSample() function.
        
        """
        # NOTE: beware factors of h - mf.M has h-1 MSun, dndm has h^4 Mpc-3
        # Below, we think we have cancelled all of those out correctly
        # For the total number of clusters, we step through the mass function calculation in z
        # For the survey volume at each z, we use the volume of that z slice
        #minMass=2e14
        #zMin=0.2
        #zMax=1.4
        #sigma8=0.80

        #mf=hmf.MassFunction(z = 0.2, Mmin = 13., Mmax = 16., delta_wrt = 'crit', delta_h = 500., sigma_8 = 0.8, cosmo_model = cosmo)

        #zRange=np.linspace(zMin, zMax, 51)
        zRange=np.linspace(zMin, zMax, 201)
        areaSr=np.radians(np.sqrt(areaDeg2))**2
        
        cosmo=cosmology.FlatLambdaCDM(H0 = H0, Om0 = Om0, Ob0 = Ob0)
        
        self.minMass=minMass
        
        # It's much faster to generate one mass function and then update its parameters (e.g., z)
        mf=hmf.MassFunction(z = zRange[0], Mmin = 13., Mmax = 16., delta_wrt = 'crit', delta_h = 500.0,
                            sigma_8 = sigma_8, cosmo_model = cosmo)#, force_flat = True, cut_fit = False)

        
        h=H0/100.
        astCalc.H0=H0
        astCalc.OMEGA_M0=Om0
        astCalc.OMEGA_L0=1.-Om0

        numClusters=0
        totalVolume=0
        numClustersByRedshift=[]
        self.tckLog10MByRedshift=[]
        self.shellVolumeMpc3ByRedshift=[]
        self.PLog10MByRedshift=[]
        self.log10MByRedshift=[]
        for i in range(len(zRange)-1):
            zShellMin=zRange[i]
            zShellMax=zRange[i+1]
            zShellMid=(zShellMax+zShellMin)/2.
            dVcdzs=[]
            zShellRange=np.linspace(zShellMin, zShellMax, 10)
            for z in zShellRange:
                dVcdzs.append(astCalc.dVcdz(z))
            dVcdzs=np.array(dVcdzs)
            shellVolumeMpc3=(areaSr*np.trapz(dVcdzs, zShellRange))/h**3 # put in units h-3 Mpc3, because dn/dm h^4 Mpc-3 MSun-1 
            self.shellVolumeMpc3ByRedshift.append(shellVolumeMpc3)
            totalVolume=totalVolume+shellVolumeMpc3
            mf.update(z = zShellMid)
            mask=np.greater(mf.M, minMass/h)
            Nmodel=np.trapz(mf.dndm[mask], mf.M[mask])*shellVolumeMpc3
            numClustersByRedshift.append(Nmodel)
            numClusters=numClusters+Nmodel
            # For M500 probability once we've picked a z (see below)
            # NOTE: commented out below for speed in likelihood, doesn't help much though
            if enableDrawSample == True:
                log10M=np.log10(mf.M[mask]/h)
                fracByM=mf.dndm[mask]/np.sum(mf.dndm[mask])
                cumFracByM=np.cumsum(fracByM) # tried dtype = longdouble, but splrep doesn't use anyway
                cMask=np.where(np.gradient(cumFracByM) != 0.0)
                tck=interpolate.splrep(cumFracByM[cMask], log10M[cMask])
                self.tckLog10MByRedshift.append(tck)
                self.PLog10MByRedshift.append(fracByM/np.trapz(fracByM, log10M))
                self.log10MByRedshift.append(log10M)
                if np.isnan(tck[1]).sum() > 0:
                    print "nan in mass interpolation"
                    IPython.embed()
                    sys.exit()
        
        self.mf=mf
        self.areaSr=areaSr
        self.volumeMpc3=totalVolume
        self.numClusters=numClusters
        self.numClustersByRedshift=np.array(numClustersByRedshift)
        self.zBinEdges=zRange
        self.z=(zRange[:-1]+zRange[1:])/2.
        
        # If we have a flat, fixed mass cut, we can use this below for z probability
        if enableDrawSample == True:
            fracByRedshift=numClustersByRedshift/np.sum(numClustersByRedshift)
            self.tckRedshift=interpolate.splrep(np.cumsum(fracByRedshift), self.z)
        
        self.enableDrawSample=enableDrawSample
    
    
    def getPLog10M(self, z):
        """Returns P(log10M), log10M  at closest given z.
        
        """
        
        mask=np.where(abs(self.z-z) == abs(self.z-z).min())[0][0]
        PLog10M=self.PLog10MByRedshift[mask]
        log10M=self.log10MByRedshift[mask]
            
        return PLog10M, log10M
        
            
    def drawSample(self, dlog10M = 0.01, dz = 0.001):
        """Draws a cluster sample from the MockSurvey. Returns it as an atpy.Table object, with columns
        name, z, dz, M500, dM500.
        
        NOTE: units of M500 are 1e14 MSun
        
        """
        
        if self.enableDrawSample == False:
            raise Exception, "enableDrawSample == False"
        
        zRolls=np.random.uniform(0, 1, self.numClusters)
        mRolls=np.random.uniform(0, 1, self.numClusters)
        
        zs=interpolate.splev(zRolls, self.tckRedshift)
        ms=[]
        for z, mRoll in zip(zs, mRolls):
            zMidIndex=np.where(abs(self.z-z) == abs(self.z-z).min())[0][0]
            ms.append(interpolate.splev(mRoll, self.tckLog10MByRedshift[zMidIndex]))
        ms=np.array(ms)
        msErr=(np.power(10, ms+dlog10M)-np.power(10, ms))/1e14
        ms=np.power(10, ms)/1e14
        
        tab=atpy.Table()
        tab.table_name="mock sample"
        tab.add_column(atpy.Column(zs, "z"))
        tab.add_column(atpy.Column([dz]*len(tab), "dz"))
        tab.add_column(atpy.Column(ms, "M500"))
        tab.add_column(atpy.Column(msErr, "dM500"))
        
        return tab
        
