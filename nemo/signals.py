"""

This module contains routines for modeling cluster and source signals.

"""

from pixell import enmap
import astropy.wcs as enwcs
import astropy.io.fits as pyfits
from astLib import *
from scipy import ndimage
from scipy import interpolate
from scipy import stats
import time
import astropy.table as atpy
from . import maps
from . import catalogs
from . import photometry
from . import gnfw
from . import plotSettings
from . import signals
import numpy as np
import numpy.fft as fft
import os
import math
import pylab as plt
import pickle
import sys
import operator
import pyximport; pyximport.install()
import nemoCython
import nemo
import glob
import shutil
import yaml
import IPython
np.random.seed()

#------------------------------------------------------------------------------------------------------------
# Global constants (we could move others here but then need to give chunky obvious names, not just e.g. h)
TCMB=2.726

#------------------------------------------------------------------------------------------------------------
def fSZ(obsFrequencyGHz):
    """Returns the frequency dependence of the (non-relativistic) Sunyaev-Zel'dovich effect.
    
    """

    h=6.63e-34
    kB=1.38e-23
    sigmaT=6.6524586e-29
    me=9.11e-31
    c=3e8
    x=(h*obsFrequencyGHz*1e9)/(kB*TCMB)
    fSZ=x*((np.exp(x)+1)/(np.exp(x)-1))-4.0
    
    return fSZ
    
#------------------------------------------------------------------------------------------------------------
def calcR500Mpc(z, M500):
    """Given z, M500 (in MSun), returns R500 in Mpc, with respect to critical density.
    
    """

    if type(M500) == str:
        raise Exception("M500 is a string - check M500MSun in your .yml config file: use, e.g., 1.0e+14 (not 1e14 or 1e+14)")

    Ez=astCalc.Ez(z)    # h(z) in Arnaud speak
    Hz=astCalc.Ez(z)*astCalc.H0  
    G=4.301e-9  # in MSun-1 km2 s-2 Mpc
    criticalDensity=(3*np.power(Hz, 2))/(8*np.pi*G)
    R500Mpc=np.power((3*M500)/(4*np.pi*500*criticalDensity), 1.0/3.0)
        
    return R500Mpc

#------------------------------------------------------------------------------------------------------------
def calcTheta500Arcmin(z, M500):
    """Given z, M500 (in MSun), returns angular size equivalent to R500, with respect to critical density.
    
    """
    
    R500Mpc=calcR500Mpc(z, M500)
    theta500Arcmin=np.degrees(np.arctan(R500Mpc/astCalc.da(z)))*60.0
    
    return theta500Arcmin
    
#------------------------------------------------------------------------------------------------------------
def makeArnaudModelProfile(z, M500, obsFreqGHz, GNFWParams = 'default'):
    """Given z, M500 (in MSun), returns dictionary containing Arnaud model profile (well, knots from spline 
    fit, 'tckP' - assumes you want to interpolate onto an array with units of degrees) and parameters 
    (particularly 'y0', 'theta500Arcmin').
    
    Use GNFWParams to specify a different shape. If GNFWParams = 'default', then the default parameters as listed
    in gnfw.py are used, i.e., 
    
    GNFWParams = {'gamma': 0.3081, 'alpha': 1.0510, 'beta': 5.4905, 'tol': 1e-7, 'npts': 100}
    
    Otherwise, give a dictionary that specifies the wanted values. This would usually be specified as
    GNFWParams in the filter params in the nemo .par file (see the example .par files).
    
    Used by ArnaudModelFilter
    
    """

    if GNFWParams == 'default':
        GNFWParams=gnfw._default_params
    
    # Adjust tol for speed vs. range of b covered
    bRange=np.linspace(0, 30, 1000)
    cylPProfile=[]
    tol=1e-6
    for i in range(len(bRange)):
        b=bRange[i]
        cylPProfile.append(gnfw.integrated(b, params = GNFWParams))
        if i > 0 and abs(cylPProfile[i] - cylPProfile[i-1]) < tol:
            break
    cylPProfile=np.array(cylPProfile)
    bRange=bRange[:i+1]
    
    # Normalise to 1 at centre
    cylPProfile=cylPProfile/cylPProfile.max()

    # Calculate R500Mpc, theta500Arcmin corresponding to given mass and redshift
    theta500Arcmin=calcTheta500Arcmin(z, M500)
    
    # Map between b and angular coordinates
    # NOTE: c500 now taken into account in gnfw.py
    thetaDegRange=bRange*(theta500Arcmin/60.)
    tckP=interpolate.splrep(thetaDegRange, cylPProfile)
    
    return {'tckP': tckP, 'theta500Arcmin': theta500Arcmin, 'rDeg': thetaDegRange}

#------------------------------------------------------------------------------------------------------------
def makeBeamModelSignalMap(degreesMap, wcs, beamFileName):
    """Makes a 2d signal only map containing the given beam.
        
    Returns signalMap (2d array), inputSignalProperties
    
    """
    
    # Load Matthew's beam profile
    beamData=np.loadtxt(beamFileName).transpose()
    profile1d=beamData[1]
    rArcmin=beamData[0]*60.0
    
    # Turn 1d profile into 2d
    rRadians=np.radians(rArcmin/60.0)
    r2p=interpolate.interp1d(rRadians, profile1d, bounds_error=False, fill_value=0.0)
    profile2d=r2p(np.radians(degreesMap))
    signalMap=profile2d

    # This used to contain info for undoing the pixel window function
    # But that's now done in mapFilters.filterMaps
    inputSignalProperties={}
                
    return signalMap, inputSignalProperties
    
#------------------------------------------------------------------------------------------------------------
def makeArnaudModelSignalMap(z, M500, obsFreqGHz, degreesMap, wcs, beamFileName, GNFWParams = 'default',
                             deltaT0 = None, maxSizeDeg = 15.0, convolveWithBeam = True):
    """Makes a 2d signal only map containing an Arnaud model cluster. Units of M500 are MSun.
    
    degreesMap is a 2d array containing radial distance from the centre - the output map will have the same
    dimensions and pixel scale (see nemoCython.makeDegreesDistanceMap).
    
    Use GNFWParams to specify a different shape. If GNFWParams = 'default', then the default parameters as listed
    in gnfw.py are used, i.e., 
    
    GNFWParams = {'gamma': 0.3081, 'alpha': 1.0510, 'beta': 5.4905, 'tol': 1e-7, 'npts': 100}
    
    Otherwise, give a dictionary that specifies the wanted values. This would usually be specified as
    GNFWParams in the filter params in the nemo .par file (see the example .par files).
    
    deltaT0 specifies the amplitude of the input signal (in map units, e.g., uK) - this is only needed if this
    routine is being used to inject sources (completely arbitrary for making filter kernels).
    
    maxSizeDeg is used to limit the region over which the beam convolution is done, for speed.
    
    If convolveWithBeam == False, no beam convolution is done (can be quicker to just do that over the whole
    source injected map rather than per object).
    
    Returns the map (2d array) and a dictionary containing the properties of the inserted cluster model.
    
    """

    # Making the 1d profile itself is the slowest part (~1 sec)
    signalDict=makeArnaudModelProfile(z, M500, obsFreqGHz, GNFWParams = GNFWParams)
    tckP=signalDict['tckP']
    theta500Arcmin=signalDict['theta500Arcmin']
    
    # Make cluster map (unit-normalised profile)
    rDeg=np.linspace(0.0, 1.0, 5000)
    profile1d=interpolate.splev(rDeg, tckP)
    if deltaT0 is not None:
        profile1d=profile1d*deltaT0
        y0=maps.convertToY(deltaT0, obsFrequencyGHz = obsFreqGHz)
    else:
        y0=1.0
    rRadians=np.radians(rDeg)
    radiansMap=np.radians(degreesMap)
    r2p=interpolate.interp1d(rRadians, profile1d, bounds_error=False, fill_value=0.0)
    profile2d=r2p(radiansMap)
    
    # Load beam profile and interpolate onto signal profile coords
    if convolveWithBeam == True:
        
        # This is slow...
        #signalMap=maps.convolveMapWithBeam(profile2d, wcs, beamFileName, maxDistDegrees = 1.0)        
        
        # This is fast...
        beamData=np.loadtxt(beamFileName).transpose()
        profile1d_beam=beamData[1]
        rDeg_beam=beamData[0]
        tck_beam=interpolate.splrep(rDeg_beam, profile1d_beam)
        profile1d_beam=interpolate.splev(rDeg, tck_beam, ext = 1) # ext = 1 sets out-of-range values to 0 rather than extrapolating beam
        
        ys, xs=np.where(degreesMap < maxSizeDeg)
        yMin=ys.min()
        yMax=ys.max()+1
        xMin=xs.min()
        xMax=xs.max()+1

        r2p_beam=interpolate.interp1d(rRadians, profile1d_beam, bounds_error=False, fill_value=0.0)
        profile2d_beam=r2p_beam(radiansMap)
        #smoothedProfile2d_1=fft.fftshift(fft.ifft2(fft.fft2(profile2d)*fft.fft2(profile2d_beam))).real
        smoothedProfile2d=np.zeros(degreesMap.shape)
        smoothedProfile2d[yMin:yMax, xMin:xMax]=fft.fftshift(fft.ifft2(fft.fft2(profile2d[yMin:yMax, xMin:xMax])*fft.fft2(profile2d_beam[yMin:yMax, xMin:xMax]))).real
        normFactor=profile2d.sum()/smoothedProfile2d.sum()
        smoothedProfile2d=smoothedProfile2d*normFactor
        signalMap=smoothedProfile2d   
        
    else:
        signalMap=profile2d
    
    # For sanity checking later, let's dump the input properties of the model into a dictionary
    inputSignalProperties={'deltaT0': deltaT0, 'y0': y0, 'theta500Arcmin': theta500Arcmin, 
                           'obsFreqGHz': obsFreqGHz}
    
    return signalMap, inputSignalProperties
    
#------------------------------------------------------------------------------------------------------------
def fitQ(config):
    """Calculates Q on a grid, and then fits (theta, Q) with a spline, saving a plot in the diagnosticDir
    and the (theta, Q) array under selFnDir.
    
    Needs a startUp.NemoConfig object.
        
    Use GNFWParams (in parDict) to specify a different shape.
    
    The calculation will be done in parallel, if MPIEnabled = True, and comm and rank are given, and extNames 
    is different for each rank. This is only needed for the first run of this routine by the nemoMass script.
    
    If extNames == [], then we figure out what the extNames are from the contents of the filteredMapsDir.
    
    NOTE: We're assuming that beamFileName is given under parDict['unfilteredMaps'].
    
    """

    if config.parDict['GNFWParams'] == 'default':
        GNFWParams=gnfw._default_params
    else:
        GNFWParams=config.parDict['GNFWParams']
        
    # Spin through the filter kernels
    photFilterLabel=config.parDict['photometryOptions']['photFilter']
    filterList=config.parDict['mapFilters']
    for f in filterList:
        if f['label'] == photFilterLabel:
            ref=f
            
    # M, z ranges for Q calc
    # NOTE: ref filter that sets scale we compare to must ALWAYS come first
    # To safely (numerically, at least) apply Q at z ~ 0.01, we need to go to theta500 ~ 500 arcmin (< 10 deg)
    MRange=[ref['params']['M500MSun']]
    zRange=[ref['params']['z']]
    #MRange=MRange+np.logspace(13.5, 16.8, 4).tolist()
    #zRange=zRange+[0.01, 0.1, 0.2, 0.3, 0.6, 0.9, 1.5, 2.0]
    # This should cover theta500Arcmin range fairly evenly without using too crazy masses
    theta500Arcmin_wanted=np.logspace(np.log10(1e-1), np.log10(500), 50)
    zRange_wanted=np.zeros(50)
    zRange_wanted[np.less(theta500Arcmin_wanted, 3.0)]=2.0
    zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 3.0), np.less(theta500Arcmin_wanted, 6.0))]=1.0
    zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 6.0), np.less(theta500Arcmin_wanted, 10.0))]=0.5
    zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 10.0), np.less(theta500Arcmin_wanted, 20.0))]=0.1
    zRange_wanted[np.logical_and(np.greater(theta500Arcmin_wanted, 20.0), np.less(theta500Arcmin_wanted, 30.0))]=0.05
    zRange_wanted[np.greater(theta500Arcmin_wanted, 30.0)]=0.01
    MRange_wanted=[]
    for theta500Arcmin, z in zip(theta500Arcmin_wanted, zRange_wanted):
        Ez=astCalc.Ez(z)
        Hz=astCalc.Ez(z)*astCalc.H0  
        G=4.301e-9  # in MSun-1 km2 s-2 Mpc
        criticalDensity=(3*np.power(Hz, 2))/(8*np.pi*G)
        R500Mpc=np.tan(np.radians(theta500Arcmin/60.0))*astCalc.da(z)
        M500=(4/3.0)*np.pi*np.power(R500Mpc, 3)*500*criticalDensity
        MRange_wanted.append(M500)
    MRange=MRange+MRange_wanted
    zRange=zRange+zRange_wanted.tolist()
    # Just for checking coverage
    #theta500Arcmin=[]
    #for M in MRange:
        #for z in zRange:
            #theta500Arcmin.append(np.degrees(np.arctan(calcR500Mpc(z, M)/astCalc.da(z)))*60.0)
    
    # Q calc - results for all tiles stored in one file
    outFileName=config.selFnDir+os.path.sep+"QFit.pickle"
    rank_QTabDict={}
    if os.path.exists(outFileName) == False:
        
        print(">>> Fitting for Q ...")
        
        # Do each tile in turn
        # Since our multi-freq filter adjusts pixel-by-pixel for frequencies available, we need to also account for that at some point...
        for extName in config.extNames:        
            
            print("... %s ..." % (extName))
            
            # Some faffing to get map pixel scale        
            with pyfits.open(config.filteredMapsDir+os.path.sep+photFilterLabel+"#%s_SNMap.fits" % (extName)) as img:
                wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
                RADeg, decDeg=wcs.getCentreWCSCoords()
                clipDict=astImages.clipImageSectionWCS(img[0].data, wcs, RADeg, decDeg, 15.0)
                wcs=clipDict['wcs']
            
            # NOTE: The block below is now only getting theta500 values, which could do elsewhere
            # So this can probably be deleted eventually...
            # Input signal maps to which we will apply filter(s)
            # We don't actually care about freq here because we deal with pressure profile itself
            # And Q measurement is relative 
            theta500Arcmin=[]
            signalMapDict={}
            signalMap=np.zeros(clipDict['data'].shape)
            degreesMap=nemoCython.makeDegreesDistanceMap(signalMap, wcs, RADeg, decDeg, 15.0)
            for z, M500MSun in zip(zRange, MRange):
                key='%.2f_%.2f' % (z, np.log10(M500MSun))
                modelDict=makeArnaudModelProfile(z, M500MSun, 148.0, GNFWParams = GNFWParams)   
                rDeg=modelDict['rDeg']
                profile1d=interpolate.splev(rDeg, modelDict['tckP'])
                r2p=interpolate.interp1d(rDeg, profile1d, bounds_error=False, fill_value=0.0)
                signalMap=r2p(degreesMap)
                # NOTE: missing a beam convolution here?
                signalMapDict[key]=signalMap
                theta500Arcmin.append(modelDict['theta500Arcmin'])
            theta500Arcmin=np.array(theta500Arcmin)

            # Set-up the beams and kernels
            # NOTE: adjusted for tileDeck files, and we have the RA, dec footprint info to deal with too
            beamsDict={}
            for mapDict in config.parDict['unfilteredMaps']:
                obsFreqGHz=mapDict['obsFreqGHz']
                beamsDict[obsFreqGHz]=mapDict['beamFileName']
            kernsDict={}            
            with pyfits.open(config.diagnosticsDir+os.path.sep+"kern2d_%s#%s.fits" % (photFilterLabel, extName)) as kernImg:
                kern2d=kernImg[0].data
                kernsDict={'kern2d': kern2d, 'header': kernImg[0].header}
            
            # Filter maps with the ref kernel
            # NOTE: keep only unique values of Q, theta500Arcmin (or interpolation routines will fail)
            Q=[]
            QTheta500Arcmin=[]
            count=0
            for z, M500MSun in zip(zRange, MRange):
                key='%.2f_%.2f' % (z, np.log10(M500MSun))
                signalMaps=[]
                y0=2e-4
                tol=1e-6
                for obsFreqGHz in list(beamsDict.keys()):
                    deltaT0=maps.convertToDeltaT(y0, obsFreqGHz)
                    # NOTE: Q is to adjust for mismatched filter shape - should this have beam in it? This does
                    signalMap, inputProperties=signals.makeArnaudModelSignalMap(z, M500MSun, obsFreqGHz, 
                                                                               degreesMap, wcs, 
                                                                               beamsDict[obsFreqGHz], 
                                                                               deltaT0 = deltaT0,
                                                                               convolveWithBeam = False)
                    if abs(y0-inputProperties['y0']) > tol:
                        raise Exception("y0 mismatch between input and output returned by makeSignalTemplateMap")
                    signalMaps.append(signalMap)
                signalMaps=np.array(signalMaps)
                peakFilteredSignals=[]
                if 'BCKSCALE' in kernsDict['header'].keys() and kernsDict['header']['BCKSCALE'] > 0.0:
                    filteredSignal=[]
                    for i in range(signalMaps.shape[0]):
                        filteredSignal.append(maps.subtractBackground(signalMaps[i], wcs, RADeg = RADeg,
                                                                      decDeg = decDeg,
                                                                      smoothScaleDeg = kernsDict['header']['BCKSCALE']/60.))
                    filteredSignal=np.array(filteredSignal)
                else:
                    filteredSignal=np.zeros(signalMaps.shape)+signalMaps
                for i in range(kern2d.shape[0]):
                    filteredSignal[i]=ndimage.convolve(filteredSignal[i], kern2d[i])
                filteredSignal=filteredSignal.sum(axis = 0)
                peakFilteredSignal=filteredSignal.max()*kernsDict['header']['SIGNORM']
                if peakFilteredSignal not in Q and theta500Arcmin[count] > 1:   # NOTE: this is to avoid problems at theta500 < beam size
                    Q.append(peakFilteredSignal)      
                    QTheta500Arcmin.append(theta500Arcmin[count])
                count=count+1
            Q=np.array(Q)
            Q=Q/Q[0]
                
            # Sort and do spline fit... save .fits table of theta, Q
            QTab=atpy.Table()
            QTab.add_column(atpy.Column(Q, 'Q'))
            QTab.add_column(atpy.Column(QTheta500Arcmin, 'theta500Arcmin'))
            QTab.sort('theta500Arcmin')
            rank_QTabDict[extName]=QTab
                       
            # Fit with spline
            tck=interpolate.splrep(QTab['theta500Arcmin'], QTab['Q'])
            
            # Plot
            plotSettings.update_rcParams()
            #fontSize=18.0
            #fontDict={'size': fontSize, 'family': 'serif'}
            plt.figure(figsize=(9,6.5))
            ax=plt.axes([0.10, 0.11, 0.88, 0.88])
            #plt.tick_params(axis='both', which='major', labelsize=15)
            #plt.tick_params(axis='both', which='minor', labelsize=15)       
            thetaArr=np.linspace(0, 500, 100000)
            #plt.plot(thetaArr, np.poly1d(coeffs)(thetaArr), 'k-')
            plt.plot(thetaArr, interpolate.splev(thetaArr, tck), 'k-')
            plt.plot(QTheta500Arcmin, Q, 'D', ms = 8)
            #plt.plot(thetaArr, simsTools.calcQ_H13(thetaArr), 'b--')
            #plt.xlim(0, 9)
            plt.ylim(0, Q.max()*1.05)
            #plt.xlim(0, thetaArr.max())
            plt.xlim(0.1, 500)
            plt.semilogx()
            plt.xlabel("$\\theta_{\\rm 500c}$ (arcmin)")
            plt.ylabel("$Q$ ($M_{\\rm 500c}$, $z$)")
            plt.savefig(config.diagnosticsDir+os.path.sep+"QFit_%s.pdf" % (extName))
            plt.close()
        
        # Gather and save all the Q fits
        if config.MPIEnabled == True:
            gathered_QTabDicts=config.comm.gather(rank_QTabDict, root = 0)
            if config.rank != 0:
                assert gathered_QTabDicts is None
                print("... MPI rank %d finished ..." % (config.rank))
                sys.exit()
            else:
                print("... gathering QTabDicts ...")
                QTabDict={}
                for tabDict in gathered_QTabDicts:
                    for key in tabDict:
                        QTabDict[key]=tabDict[key]
        else:
            QTabDict=rank_QTabDict
                    
        with open(outFileName, "wb") as pickleFile:
            pickler=pickle.Pickler(pickleFile)
            pickler.dump(QTabDict)
            
    else:
        
        if config.MPIEnabled == True and config.rank != 0:
            sys.exit()
        
        print(">>> Loading previously cached Q fit ...")
        return loadQ(outFileName)
        
    tckDict={}
    for key in QTabDict:
        tckDict[key]=interpolate.splrep(QTabDict[key]['theta500Arcmin'], QTabDict[key]['Q'])
    
    return tckDict

#------------------------------------------------------------------------------------------------------------
def loadQ(pickledQFileName):
    """Loads tckQFitDict from given path.
    
    """
    
    with open(pickledQFileName, "rb") as pickleFile:
        unpickler=pickle.Unpickler(pickleFile)
        QTabDict=unpickler.load()
    
    tckDict={}
    for key in QTabDict:
        tckDict[key]=interpolate.splrep(QTabDict[key]['theta500Arcmin'], QTabDict[key]['Q'])
    
    return tckDict

#------------------------------------------------------------------------------------------------------------
def calcQ(theta500Arcmin, tck):
    """Returns Q, given theta500Arcmin, and a set of spline fit knots for (theta, Q).
    
    """
    
    #Q=np.poly1d(coeffs)(theta500Arcmin)
    Q=interpolate.splev(theta500Arcmin, tck)
    
    return Q

#------------------------------------------------------------------------------------------------------------
def calcWeightedFRel(z, M500, fRelWeightsDict):
    """Return fRel for the given (z, M500), weighted by frequency according to fRelWeightsDict
    
    """
    
    fRels=[]
    freqWeights=[]
    for obsFreqGHz in fRelWeightsDict.keys():
        if fRelWeightsDict[obsFreqGHz] > 0:
            fRels.append(calcFRel(z, M500, obsFreqGHz = obsFreqGHz))
            freqWeights.append(fRelWeightsDict[obsFreqGHz])
    fRel=np.average(fRels, weights = freqWeights)
    
    return fRel
    
#------------------------------------------------------------------------------------------------------------
def calcFRel(z, M500, obsFreqGHz = 148.0):
    """Calculates relativistic correction to SZ effect at specified frequency, given z, M500 in MSun.
       
    This assumes the Arnaud et al. (2005) M-T relation, and applies formulae of Itoh et al. (1998)
    
    As for H13, we return fRel = 1 + delta_SZE (see also Marriage et al. 2011)

    """
    
    # NOTE: we should define constants somewhere else...
    h=6.63e-34
    kB=1.38e-23
    sigmaT=6.6524586e-29
    me=9.11e-31
    e=1.6e-19
    c=3e8
    TCMB=2.726
    
    # Using Arnaud et al. (2005) M-T to get temperature
    A=3.84e14
    B=1.71
    TkeV=5.*np.power(((astCalc.Ez(z)*M500)/A), 1/B)
    TKelvin=TkeV*((1000*e)/kB)

    # Itoh et al. (1998) eqns. 2.25 - 2.30
    thetae=(kB*TKelvin)/(me*c**2)
    X=(h*obsFreqGHz*1e9)/(kB*TCMB)
    Xtw=X*(np.cosh(X/2.)/np.sinh(X/2.))
    Stw=X/np.sinh(X/2.)

    Y0=-4+Xtw

    Y1=-10. + (47/2.)*Xtw - (42/5.)*Xtw**2 + (7/10.)*Xtw**3 + np.power(Stw, 2)*(-(21/5.) + (7/5.)*Xtw)

    Y2=-(15/2.) +  (1023/8.)*Xtw - (868/5.)*Xtw**2 + (329/5.)*Xtw**3 - (44/5.)*Xtw**4 + (11/30.)*Xtw**5 \
        + np.power(Stw, 2)*(-(434/5.) + (658/5.)*Xtw - (242/5.)*Xtw**2 + (143/30.)*Xtw**3) \
        + np.power(Stw, 4)*(-(44/5.) + (187/60.)*Xtw)

    Y3=(15/2.) + (2505/8.)*Xtw - (7098/5.)*Xtw**2 + (14253/10.)*Xtw**3 - (18594/35.)*Xtw**4 + (12059/140.)*Xtw**5 - (128/21.)*Xtw**6 + (16/105.)*Xtw**7 \
        + np.power(Stw, 2)*(-(7098/10.) + (14253/5.)*Xtw - (102267/35.)*Xtw**2 + (156767/140.)*Xtw**3 - (1216/7.)*Xtw**4 + (64/7.)*Xtw**5) \
        + np.power(Stw, 4)*(-(18594/35.) + (205003/280.)*Xtw - (1920/7.)*Xtw**2 + (1024/35.)*Xtw**3) \
        + np.power(Stw, 6)*(-(544/21.) + (992/105.)*Xtw)

    Y4=-(135/32.) + (30375/128.)*Xtw - (62391/10.)*Xtw**2 + (614727/40.)*Xtw**3 - (124389/10.)*Xtw**4 \
        + (355703/80.)*Xtw**5 - (16568/21.)*Xtw**6 + (7516/105.)*Xtw**7 - (22/7.)*Xtw**8 + (11/210.)*Xtw**9 \
        + np.power(Stw, 2)*(-(62391/20.) + (614727/20.)*Xtw - (1368279/20.)*Xtw**2 + (4624139/80.)*Xtw**3 - (157396/7.)*Xtw**4 \
        + (30064/7.)*Xtw**5 - (2717/7.)*Xtw**6 + (2761/210.)*Xtw**7) \
        + np.power(Stw, 4)*(-(124389/10.) + (6046951/160.)*Xtw - (248520/7.)*Xtw**2 + (481024/35.)*Xtw**3 - (15972/7.)*Xtw**4 + (18689/140.)*Xtw**5) \
        + np.power(Stw, 6)*(-(70414/21.) + (465992/105.)*Xtw - (11792/7.)*Xtw**2 + (19778/105.)*Xtw**3) \
        + np.power(Stw, 8)*(-(682/7.) + (7601/210.)*Xtw)

    deltaSZE=((X**3)/(np.exp(X)-1)) * ((thetae*X*np.exp(X))/(np.exp(X)-1)) * (Y0 + Y1*thetae + Y2*thetae**2 + Y3*thetae**3 + Y4*thetae**4)

    fRel=1+deltaSZE
    
    return fRel

#------------------------------------------------------------------------------------------------------------
def getM500FromP(P, log10M, calcErrors = True):
    """Returns M500 as the maximum likelihood value from given P(log10M) distribution, together with 
    1-sigma error bars (M500, -M500Err, +M500 err).
    
    """

    # Find max likelihood and integrate to get error bars
    tckP=interpolate.splrep(log10M, P)
    fineLog10M=np.linspace(log10M.min(), log10M.max(), 10000)
    fineP=interpolate.splev(fineLog10M, tckP)
    fineP=fineP/np.trapz(fineP, fineLog10M)
    index=np.argmax(fineP)
    
    clusterLogM500=fineLog10M[index]
    clusterM500=np.power(10, clusterLogM500)/1e14

    if calcErrors == True:
        for n in range(fineP.shape[0]):
            minIndex=index-n
            maxIndex=index+n
            if minIndex < 0 or maxIndex > fineP.shape[0]:
                # This shouldn't happen; if it does, probably y0 is in the wrong units
                # Previously we threw an exception here, but we can't if using this for forced photometry
                print("WARNING: outside M500 range - check y0 units or for problem at cluster location in map (if not in forced photometry mode)")
                clusterM500MinusErr=0.
                clusterM500PlusErr=0.
                break
            p=np.trapz(fineP[minIndex:maxIndex], fineLog10M[minIndex:maxIndex])
            if p >= 0.6827:
                clusterLogM500Min=fineLog10M[minIndex]
                clusterLogM500Max=fineLog10M[maxIndex]
                clusterM500MinusErr=(np.power(10, clusterLogM500)-np.power(10, clusterLogM500Min))/1e14
                clusterM500PlusErr=(np.power(10, clusterLogM500Max)-np.power(10, clusterLogM500))/1e14
                break        
    else:
        clusterM500MinusErr=0.
        clusterM500PlusErr=0.
    
    return clusterM500, clusterM500MinusErr, clusterM500PlusErr

#------------------------------------------------------------------------------------------------------------
def y0FromLogM500(log10M500, z, tckQFit, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, sigma_int = 0.2,
                  H0 = None, OmegaM0 = None, OmegaL0 = None, fRelWeightsDict = {148.0: 1.0}):
    """Predict y0~ given logM500 (in MSun) and redshift. Default scaling relation parameters are A10 (as in
    H13).
    
    Use H0, OmegaM0, OmegaL0 to specify different cosmological parameters to the default stored under
    astCalc.H0, astCalc.OMEGA_M0, astCalc.OMEGA_L0 (used for E(z) and angular size calculation). If these
    are set to None, the defaults will be used.
    
    fRelWeightsDict is used to account for the relativistic correction when y0~ has been constructed
    from multi-frequency maps. Weights should sum to 1.0; keys are observed frequency in GHz.
    
    Returns y0~, theta500Arcmin, Q
    
    """

    if type(Mpivot) == str:
        raise Exception("Mpivot is a string - check Mpivot in your .yml config file: use, e.g., 3.0e+14 (not 3e14 or 3e+14)")
    
    # Change cosmology for this call, if required... then put it back afterwards (just in case)
    if H0 != None:
        oldH0=astCalc.H0
        astCalc.H0=H0
    if OmegaM0 != None:
        oldOmegaM0=astCalc.OMEGA_M0
        astCalc.OMEGA_M0=OmegaM0
    if OmegaL0 != None:
        oldOmegaL0=astCalc.OMEGA_L0
        astCalc.OMEGA_L0=OmegaL0
    
    # Filtering/detection was performed with a fixed fiducial cosmology... so we don't need to recalculate Q
    # We just need to recalculate theta500Arcmin and E(z) only
    M500=np.power(10, log10M500)
    theta500Arcmin=calcTheta500Arcmin(z, M500)
    Q=calcQ(theta500Arcmin, tckQFit)
    
    # Relativistic correction: now a little more complicated, to account for fact y0~ maps are weighted sum
    # of individual frequency maps, and relativistic correction size varies with frequency
    fRels=[]
    freqWeights=[]
    for obsFreqGHz in fRelWeightsDict.keys():
        fRels.append(calcFRel(z, M500, obsFreqGHz = obsFreqGHz))
        freqWeights.append(fRelWeightsDict[obsFreqGHz])
    fRel=np.average(fRels, weights = freqWeights)
    
    # UPP relation according to H13
    # NOTE: m in H13 is M/Mpivot
    # NOTE: this goes negative for crazy masses where the Q polynomial fit goes -ve, so ignore those
    y0pred=tenToA0*np.power(astCalc.Ez(z), 2)*np.power(M500/Mpivot, 1+B0)*Q*fRel
            
    # Restore cosmology if we changed it for this call
    if H0 != None:
        astCalc.H0=oldH0
    if OmegaM0 != None:
        astCalc.OMEGA_M0=oldOmegaM0
    if OmegaL0 != None:
        astCalc.OMEGA_L0=oldOmegaL0
        
    return y0pred, theta500Arcmin, Q
            
#------------------------------------------------------------------------------------------------------------
def calcM500Fromy0(y0, y0Err, z, zErr, tckQFit, mockSurvey, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, 
                   sigma_int = 0.2, applyMFDebiasCorrection = True, calcErrors = True,
                   fRelWeightsDict = {148.0: 1.0}):
    """Returns M500 +/- errors in units of 10^14 MSun, calculated assuming a y0 - M relation (default values
    assume UPP scaling relation from Arnaud et al. 2010), taking into account the steepness of the mass
    function. The approach followed is described in H13, Section 3.2.
    
    Here, mockSurvey is a MockSurvey object. We're using this to handle the halo mass function calculations
    (in turn using the hmf module). Supplying mockSurvey is no longer optional (and handles setting the 
    cosmology anyway when initialised or updated).
    
    tckQFit is a set of spline knots, as returned by fitQ.
    
    If applyMFDebiasCorrection == True, apply correction that accounts for steepness of mass function.
    
    If calcErrors == False, error bars are not calculated, they are just set to zero.
    
    fRelWeightsDict is used to account for the relativistic correction when y0~ has been constructed
    from multi-frequency maps. Weights should sum to 1.0; keys are observed frequency in GHz.
        
    Returns dictionary with keys M500, M500_errPlus, M500_errMinus
    
    """
    
    if y0 < 0:
        raise Exception('y0 cannot be negative')
    if y0 > 1e-2:
        raise Exception('y0 is suspiciously large - probably you need to multiply by 1e-4')
            
    P=calcPM500(y0, y0Err, z, zErr, tckQFit, mockSurvey, tenToA0 = tenToA0, B0 = B0, Mpivot = Mpivot, 
                sigma_int = sigma_int, applyMFDebiasCorrection = applyMFDebiasCorrection,
                fRelWeightsDict = fRelWeightsDict)
    
    M500, errM500Minus, errM500Plus=getM500FromP(P, mockSurvey.log10M, calcErrors = calcErrors)
    
    return {'M500': M500, 'M500_errPlus': errM500Plus, 'M500_errMinus': errM500Minus}

#------------------------------------------------------------------------------------------------------------
def calcPM500(y0, y0Err, z, zErr, tckQFit, mockSurvey, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, sigma_int = 0.2, 
             applyMFDebiasCorrection = True, fRelWeightsDict = {148.0: 1.0}, return2D = False):
    """Calculates P(M500) assuming a y0 - M relation (default values assume UPP scaling relation from Arnaud 
    et al. 2010), taking into account the steepness of the mass function. The approach followed is described 
    in H13, Section 3.2. The binning for P(M500) is set according to the given mockSurvey, as are the assumed
    cosmological parameters.
    
    This routine is used by calcM500Fromy0.
    
    If return2D == True, returns a grid of same dimensions / binning as mockSurvey.z, mockSurvey.log10M,
    normalised such that the sum of the values is 1.
    
    """
    
    # For marginalising over photo-z errors (we assume +/-5 sigma is accurate enough)
    if zErr > 0:
        zMin=z-zErr*5
        zMax=z+zErr*5
        zMask=np.logical_and(np.greater_equal(mockSurvey.z, zMin), np.less(mockSurvey.z, zMax))
        zRange=mockSurvey.z[zMask]
        #if zMin <= 0:
            #zMin=1e-3
        #zRange=np.arange(zMin, zMax, 0.005)
        Pz=np.exp(-np.power(z-zRange, 2)/(2*(np.power(zErr, 2))))
        Pz=Pz/np.trapz(Pz, zRange)
    else:
        zRange=[z]
        Pz=np.ones(len(zRange))

    log_y0=np.log(y0)
    log_y0Err=y0Err/y0
        
    # NOTE: Swap below if want to use bigger log10M range...
    log10Ms=mockSurvey.log10M
    #log10MStep=mockSurvey.log10M[1]-mockSurvey.log10M[0]
    #log10Ms=np.arange(-100.0, 100.0, log10MStep)
        
    PArr=[]
    for k in range(len(zRange)):
        
        zk=zRange[k]        
        mockSurvey_zIndex=np.argmin(abs(mockSurvey.z-zk))

        theta500s=interpolate.splev(log10Ms, mockSurvey.theta500Splines[mockSurvey_zIndex], ext = 3)
        Qs=interpolate.splev(theta500s, tckQFit, ext = 3)
        fRels=interpolate.splev(log10Ms, mockSurvey.fRelSplines[mockSurvey_zIndex], ext = 3)   
        fRels[np.less_equal(fRels, 0)]=1e-4   # For extreme masses (> 10^16 MSun) at high-z, this can dip -ve
        y0pred=tenToA0*np.power(mockSurvey.Ez[mockSurvey_zIndex], 2)*np.power(np.power(10, log10Ms)/Mpivot, 1+B0)*Qs*fRels
        if np.less(y0pred, 0).sum() > 0:
            # This generally means we wandered out of where Q is defined (e.g., beyond mockSurvey log10M limits)
            # Or fRel can dip -ve for extreme mass at high-z (can happen with large Om0)
            raise Exception("some predicted y0 values -ve")
        log_y0pred=np.log(y0pred)

        Py0GivenM=np.exp(-np.power(log_y0-log_y0pred, 2)/(2*(np.power(log_y0Err, 2)+np.power(sigma_int, 2))))
        Py0GivenM=Py0GivenM/np.trapz(Py0GivenM, log10Ms)

        # Mass function de-bias
        if applyMFDebiasCorrection == True:
            PLog10M=mockSurvey.getPLog10M(zk)
            PLog10M=PLog10M/np.trapz(PLog10M, log10Ms)
        else:
            PLog10M=1.0
        
        P=Py0GivenM*PLog10M*Pz[k]
        PArr.append(P)
        
    # 2D PArr is what we would want to project onto (M, z) grid
    PArr=np.array(PArr)
        
    # Marginalised over z uncertainty
    P=np.sum(PArr, axis = 0)
    P=P/np.trapz(P, log10Ms)
    
    # Reshape to (M, z) grid - use this if use different log10M range to mockSurvey
    #tck=interpolate.splrep(log10Ms, P)
    #P=interpolate.splev(mockSurvey.log10M, tck, ext = 1)
        
    if return2D == True:
        P2D=np.zeros(mockSurvey.clusterCount.shape)
        if zErr == 0:
            P2D[np.argmin(abs(mockSurvey.z-z))]=PArr
        else:
            P2D[zMask]=PArr
        P=P2D/P2D.sum()
        #astImages.saveFITS("test.fits", P.transpose(), None)
        
    return P

#------------------------------------------------------------------------------------------------------------
# Mass conversion routines

# For getting x(f) - see Hu & Kravtsov
x=np.linspace(1e-3, 10, 1000)
fx=(x**3)*(np.log(1+1./x)-np.power(1+x, -1))
XF_TCK=interpolate.splrep(fx, x)
FX_TCK=interpolate.splrep(x, fx)

#------------------------------------------------------------------------------------------------------------
def gz(zIn, zMax = 1000, dz = 0.1):
    """Calculates linear growth factor at redshift z. Use Dz if you want normalised to D(z) = 1.0 at z = 0.
    
    See http://www.astronomy.ohio-state.edu/~dhw/A873/notes8.pdf for some notes on this.
    
    """
    
    zRange=np.arange(zIn, zMax, dz)
    HzPrime=[]
    for zPrime in zRange:
        HzPrime.append(astCalc.Ez(zPrime)*astCalc.H0)
    HzPrime=np.array(HzPrime)
    gz=astCalc.Ez(zIn)*np.trapz((np.gradient(zRange)*(1+zRange)) / np.power(HzPrime, 3), zRange)
    
    return gz

#------------------------------------------------------------------------------------------------------------
def calcDz(zIn):
    """Calculate linear growth factor, normalised to D(z) = 1.0 at z = 0.
    
    """
    return gz(zIn)/gz(0.0)

#------------------------------------------------------------------------------------------------------------
def criticalDensity(z):
    """Returns the critical density at the given z.
    
    """
    
    G=4.301e-9  # in MSun-1 km2 s-2 Mpc, see Robotham GAMA groups paper
    Hz=astCalc.H0*astCalc.Ez(z)
    rho_crit=((3*np.power(Hz, 2))/(8*np.pi*G))
    
    return rho_crit

#------------------------------------------------------------------------------------------------------------
def meanDensity(z):
    """Returns the mean density at the given z.
    
    """
    
    rho_mean=astCalc.OmegaMz(z)*criticalDensity(z)
    
    return rho_mean  

#------------------------------------------------------------------------------------------------------------
def convertM200mToM500c(M200m, z):
    """Returns M500c (MSun), R500c (Mpc) for given M200m and redshift. Uses the Bhattacharya et al. c-M
    relation: http://adsabs.harvard.edu/abs/2013ApJ...766...32B

    See also Hu & Kravtsov: http://iopscience.iop.org/article/10.1086/345846/pdf
    
    """
        
    # c-M relation for full cluster sample
    Dz=calcDz(z)    # <--- this is the slow part. 3 seconds!
    nu200m=(1./Dz)*(1.12*np.power(M200m / (5e13 * np.power(astCalc.H0/100., -1)), 0.3)+0.53)
    c200m=np.power(Dz, 1.15)*9.0*np.power(nu200m, -0.29)
    
    rho_crit=criticalDensity(z)
    rho_mean=meanDensity(z)
    R200m=np.power((3*M200m)/(4*np.pi*200*rho_mean), 1/3.)

    rs=R200m/c200m
    
    f_rsOverR500c=((500*rho_crit) / (200*rho_mean)) * interpolate.splev(1./c200m, FX_TCK)
    x_rsOverR500c=interpolate.splev(f_rsOverR500c, XF_TCK)
    R500c=rs/x_rsOverR500c

    M500c=(4/3.0)*np.pi*R500c**3*(500*rho_crit)
        
    return M500c, R500c

#------------------------------------------------------------------------------------------------------------
def convertM500cToM200m(M500c, z):
    """Returns M200m given M500c
    
    """
    
    tolerance=1e-5
    scaleFactor=3.0
    ratio=1e6
    count=0
    while abs(1.0-ratio) > tolerance:
        testM500c, testR500c=convertM200mToM500c(scaleFactor*M500c, z)
        ratio=M500c/testM500c
        scaleFactor=scaleFactor*ratio
        count=count+1
        if count > 10:
            raise Exception("M500c -> M200m conversion didn't converge quickly enough")
        
    M200m=scaleFactor*M500c
    
    return M200m
