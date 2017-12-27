# -*- coding: utf-8 -*-
"""This module contains routines for comparing measured fluxes to input sims.

"""

import pyfits
from astLib import *
from scipy import ndimage
from scipy import interpolate
from scipy import stats
import time
import astropy.table as atpy
import mapTools
from flipper import liteMap
from flipper import fftTools
import catalogTools
import photometry
import gnfw
import plotSettings
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
import actDict
import glob
import IPython
np.random.seed()

#------------------------------------------------------------------------------------------------------------
def parseInputSimCatalog(fileName, wcs):
    """Parses input simulation catalog (i.e. as produced by Hy and friends, well actually they look like 
    they've been fiddled a bit by Ryan), returns dictionary list. Only returns objects that are found within
    the map boundaries defined by the wcs.
        
    """
    
    inFile=file(fileName, "r")
    lines=inFile.readlines()
    inFile.close()
    
    catalog=[]
    idNum=0
    for line in lines:
        if len(line) > 3 and line[0] != "#":
            objDict={}
            bits=line.split()
            idNum=idNum+1
            objDict['id']=idNum
            objDict['z']=float(bits[2])
            objDict['RADeg']=float(bits[0])
            objDict['decDeg']=float(bits[1])
            objDict['name']=catalogTools.makeACTName(objDict['RADeg'], objDict['decDeg'])
            objDict['Mvir']=float(bits[3]) 
            objDict['RvirMpc']=float(bits[4])
            objDict['fluxFromCatalog_arcmin2']=float(bits[5])
            if wcs.coordsAreInImage(objDict['RADeg'], objDict['decDeg']) == True:
                catalog.append(objDict)
    
    # Check this works okay
    catalogTools.catalog2DS9(catalog, 'inputSimObjects.reg')
    
    return catalog

#------------------------------------------------------------------------------------------------------------
def matchAgainstSimCatalog(catalog, inputSimCatalog, simFluxKey = 'fixedApertureFluxFromInputSimMap_arcmin2'):
    """Matches the given catalog against the given input sim catalog, adding the input sim flux, mass etc.
    info for all matched objects to the catalog. Needed for feeding the catalog into routines that estimate
    completeness, purity.
    
    """
    
    # For faster cross matching
    sRAs=[]
    sDecs=[]
    for s in inputSimCatalog:
        sRAs.append(s['RADeg'])
        sDecs.append(s['decDeg'])
    sRAs=np.array(sRAs)
    sDecs=np.array(sDecs)
    
    for m in catalog:
        mra=m['RADeg']
        mdec=m['decDeg']
        rMin=1e6
        bestMatch=None
                        
        # Faster matching - best matches here by definition only actually had sim fluxes measured
        # if they passed the Y, z limit cuts. simCatalog still contains every object in map area.
        rs=astCoords.calcAngSepDeg(mra, mdec, sRAs, sDecs)
        rMin=rs.min()
        rMinIndex=np.equal(rs, rMin).nonzero()[0][0]
        bestMatch=inputSimCatalog[rMinIndex]
        if bestMatch != None and rMin < catalogTools.XMATCH_RADIUS_DEG:
            # We want to track all matches we detected even if we don't necessarily want to compare flux
            m['inputSim_Mvir']=bestMatch['Mvir']
            m['inputSim_z']=bestMatch['z']
            m['inputSim_RvirMpc']=bestMatch['RvirMpc']
            if simFluxKey in bestMatch.keys():
                m['inputSim_flux_arcmin2']=bestMatch[simFluxKey]
            else:
                m['inputSim_flux_arcmin2']=None
        else:
            m['inputSim_flux_arcmin2']=None
            m['inputSim_z']=None
            m['inputSim_RvirMpc']=None
            m['inputSim_Mvir']=None
               
#------------------------------------------------------------------------------------------------------------
def compareToInputSimCatalog(imageDict, inputSimCatFileName, inputSimMapFileName, photometryOptions, 
                             outDir = ".", YLimit = None, zRange = None, convertFromJySrToMicroK = False,
                             obsFreqGHz = 148, clusterProfilesToPlot = []):
    """This is a wrapper for calling all the stuff we might want to check against the
    input sims, e.g. purity, completeness etc. Saves out results under outDir.
    
    Can use clusterProfilesToPlot to list names of objects to make plots of e.g. cumulative flux vs. radius
    for each filtered map.
    
    """
    
    print ">>> Checking fluxes against input sim ..."        
    if os.path.exists(outDir) == False:
        os.makedirs(outDir)
    
    signalTemplatesList=[]
    for k in imageDict.keys():
        if k != "mergedCatalog" and k != "optimalCatalog":
            signalTemplatesList.append(k)
    
    # Need a data map WCS for working out which objects in the sim catalog are in our data map
    dataMapWCS=astWCS.WCS(imageDict[signalTemplatesList[0]]['ycFilteredMap'])
    
    # Get sim catalog
    simCatalog=getInputSimApertureFluxes(inputSimCatFileName, inputSimMapFileName, dataMapWCS, \
                                         apertureRadiusArcmin = photometryOptions["apertureRadiusArcmin"], \
                                         YLimit = YLimit, zRange = zRange, \
                                         convertFromJySrToMicroK = convertFromJySrToMicroK, obsFreqGHz = obsFreqGHz)
    simFluxKey='fixedApertureFluxFromInputSimMap_arcmin2'

    # Checking fluxes in maps/catalog agree - they do, so this should no longer be necessary
    #simCatalog=checkInputSimCatalogFluxes(inputSimCatFileName, inputSimMapFileName, YLimit = YLimit, \
                            #zRange = zRange, convertFromJySrToMicroK = convertFromJySrToMicroK, \
                            #obsFreqGHz = obsFreqGHz)
    
    # We will want to make plots of optimal template, plus see individual tenplate results at the same time
    # For both the flux comparison and completeness/purity plots
    catalogsDict={}
    catalogsDict['optimal']=imageDict['optimalCatalog']
    for st in signalTemplatesList:
        catalogsDict[st]=imageDict[st]['catalog']
    
    # Match all catalogs against input sim
    matchAgainstSimCatalog(imageDict['optimalCatalog'], simCatalog)
    for key in catalogsDict.keys():
        matchAgainstSimCatalog(catalogsDict[key], simCatalog)

    # Flux comparison plot
    fig=plt.figure(num=4, figsize=(10,8))
    fig.canvas.set_window_title('Measured vs. Input Flux Comparison')
    topPlot=plt.axes([0.12, 0.35, 0.8, 0.6])
    bottomPlot=plt.axes([0.12, 0.1, 0.8, 0.23])
    keysList=catalogsDict.keys()
    keysList.reverse()  # so optimal gets plotted first
    for key in keysList:
        catalog=catalogsDict[key]
        names=[]
        simFlux=[]
        myFlux=[]
        myFluxErr=[]
        for obj in catalog:
            if 'inputSim_flux_arcmin2' in obj.keys() and obj['inputSim_flux_arcmin2'] != None:
                names.append(obj['name'])
                myFlux.append(obj['flux_arcmin2'])
                myFluxErr.append(obj['fluxErr_arcmin2'])
                simFlux.append(obj['inputSim_flux_arcmin2'])
        simFlux=np.array(simFlux)
        myFlux=np.array(myFlux)
        myFluxErr=np.array(myFluxErr)
        if len(simFlux) > 0:
            if key == 'optimal':
                plt.axes(topPlot)
                plt.errorbar(simFlux, myFlux, yerr=myFluxErr, fmt='ko', label = 'Optimal S/N template')
                plt.axes(bottomPlot)
                plt.errorbar(simFlux, simFlux-myFlux, yerr=myFluxErr, fmt='ko', label = 'Optimal S/N template')
            else:
                plt.axes(topPlot)
                plt.errorbar(simFlux, myFlux, yerr=myFluxErr, fmt='.', label = key)
                plt.axes(bottomPlot)
                plt.errorbar(simFlux, simFlux-myFlux, yerr=myFluxErr, fmt='.', label = key)

    plt.axes(bottomPlot)
    plotRange=np.linspace(0, 0.02, 10)
    plt.plot(plotRange, [0]*len(plotRange), 'k--')
    plt.ylim(-0.0015, 0.0015)
    plt.xlim(0.0, 0.0015)
    plt.ylabel('$\Delta$Y (arcmin$^2$)')
    plt.xlabel('Input Sim Catalog Y (arcmin$^2$)')
    
    plt.axes(topPlot)
    plotRange=np.linspace(0, 0.02, 10)
    plt.plot(plotRange, plotRange, 'k--')
    plt.ylabel('Measured Y (arcmin$^2$)')
    plt.xticks([], [])
    plt.xticks(bottomPlot.get_xticks())
    plt.xlim(0.0, 0.0015)
    plt.ylim(0.0, 0.0015)
    #legAxes=plt.axes([0.65, 0.3, 0.2, 0.4], frameon=False)
    #plt.xticks([], [])
    #plt.yticks([], [])
    #plt.legend(loc="center")
    plt.legend(loc="best", prop=plt.matplotlib.font_manager.FontProperties(family='sans-serif', size=10))
    plt.savefig(outDir+os.path.sep+"catalogFluxesVsMeasuredFluxes.png")

    #print "In flux recovery bit"
    #ipshell()
    #sys.exit()
    
    # Simple flux recovery stats
    # These are in units of 1-sigma error bars
    # So sigmaResidualSigma = 3 would mean that the standard deviation of my fluxes is 3 error bars
    # and medianResidualSigma = 0.3 means that we're reasonably unbiased (within 0.3 error bars)
    medianResidualSigma=np.median((simFlux-myFlux)/myFluxErr)
    sigmaResidualSigma=np.std(((simFlux-myFlux)/myFluxErr))
    print ">>> Flux recovery stats:"
    print "... median residual = %.3f error bars" % (medianResidualSigma)
    print "... stdev residual = %.3f error bars" % (sigmaResidualSigma)
    
    # Now we want to do completeness and purity as function of S/N and mass
    fig=plt.figure(num=5, figsize=(10,8))
    fig.canvas.set_window_title('Completeness and Purity')
    plt.subplots_adjust(left=0.1, bottom=0.07, right=0.95, top=0.95, wspace=0.02, hspace=0.27)
    compAxes=plt.subplot(311)
    trueDetAxes=plt.subplot(312)
    purityAxes=plt.subplot(313)
    keysList.reverse()   # in this case, we want optimal plotted last

    # Completeness
    simMasses=[]
    for obj in simCatalog:
        simMasses.append(obj['Mvir'])  
    simMasses=np.array(simMasses)
    simMasses.sort()
    simNum=np.ones(simMasses.shape[0])
    
    for key in keysList:
        catalog=catalogsDict[key]
        
        # Completeness
        detected=np.zeros(simMasses.shape[0])
        for obj in catalog:
            if obj['inputSim_Mvir'] != None:
                index=np.equal(simMasses, obj['inputSim_Mvir']).nonzero()[0][0]
                detected[index]=detected[index]+1    
        cumTotal=1+simNum.sum()-simNum.cumsum()
        cumDetected=1+detected.sum()-detected.cumsum()
        plt.axes(compAxes)
        if key == 'optimal':
            plt.plot(simMasses/1e14, cumDetected/cumTotal, 'k', lw=2, label = 'Optimal S/N template')
        else:
            plt.plot(simMasses/1e14, cumDetected/cumTotal, '--', label = key)
        
        # True detections
        realObjs=[]
        for obj in catalog:
            if obj['inputSim_Mvir'] != None:
                realObjs.append([obj['SNR'], 1])
            else:
                realObjs.append([obj['SNR'], 0])
        realObjs=sorted(realObjs, key=operator.itemgetter(0))
        realObjs=np.array(realObjs).transpose()
        plt.axes(trueDetAxes)
        if key == 'optimal':
            # yes, +1 is essential
            plt.plot(realObjs[0], 1+realObjs[1].sum()-realObjs[1].cumsum(), 'k', lw=2, label = 'Optimal S/N template')
        else:
            plt.plot(realObjs[0], 1+realObjs[1].sum()-realObjs[1].cumsum(), '--', label = key)

        # Purity
        purity=(1+realObjs[1].sum()-realObjs[1].cumsum())/ \
               (realObjs.shape[1]-np.arange(realObjs[1].shape[0], dtype=float))
        plt.axes(purityAxes)
        if key == 'optimal':
            plt.plot(realObjs[0], purity, 'k', lw=2, label = 'Optimal S/N template')
        else:
            plt.plot(realObjs[0], purity, '--', label = key)
            
    # Completeness
    plt.axes(compAxes)
    plt.xlabel("Mass ($\\times 10^{14}$ M$_\odot$)")
    plt.ylabel("Completeness > Mass")
    plt.ylim(0, 1.1)
    
    # True detections
    plt.axes(trueDetAxes)
    plt.xlabel("S/N")
    plt.ylabel("True Detections > S/N")
    plt.legend(loc="upper right", prop=plt.matplotlib.font_manager.FontProperties(family='sans-serif', size=10))
    
    # Purity
    plt.axes(purityAxes)
    plt.xlabel("S/N")
    plt.ylabel("Purity > S/N")
    plt.ylim(0, 1.1)
    
    plt.savefig(outDir+os.path.sep+"completenessAndPurity.png")
        
#------------------------------------------------------------------------------------------------------------
def getInputSimApertureFluxes(inputSimCatFileName, inputSimMapFileName, dataMapWCS, 
                              apertureRadiusArcmin = 3.0, YLimit = None, zRange = None, 
                              convertFromJySrToMicroK = False, obsFreqGHz = 148, saveAsY = True):
    """This parses the input sim catalog file, adding every object in it that falls within the map
    to a list of objects. For objects that additionally pass the given Y, z cuts, their fluxes
    are measured directly from the input sim map through the specified circular aperture.
    
    dataMapWCS is needed to work out which input sim catalog objects are actually in the map. Speeds things
    up!
    
    Note that YLimit, zLimit are on quantities in the input catalog, which are measured within virial radius.
    
    """
    
    # We may want these again, so pickle the results
    if os.path.exists("nemoCache") == False:
        os.makedirs("nemoCache")
        
    pickleFileName="nemoCache/inputSimFluxes.pickled"
    if os.path.exists(pickleFileName) == True:
        print ">>> Loading previously measured input sim map fluxes ..."
        pickleFile=file(pickleFileName, "r")
        unpickler=pickle.Unpickler(pickleFile)
        inputCatalog=unpickler.load()
    else:
        
        print ">>> Loading noiseless input sim catalog and getting fluxes from input sim map  ... "
        
        inputCatalog=parseInputSimCatalog(inputSimCatFileName, dataMapWCS)

        # Measure fluxes in actual input map from directly adding up flux within given aperture           
        img=pyfits.open(inputSimMapFileName)
        wcs=astWCS.WCS(inputSimMapFileName)
        data=img[0].data
        if convertFromJySrToMicroK == True: # from Jy/sr
            if obsFreqGHz == 148:
                print ">>> Converting from Jy/sr to micro kelvins assuming obsFreqGHz = %.0f" % (obsFreqGHz)
                data=(data/1.072480e+09)*2.726*1e6
            else:
                raise Exception, "no code added to support conversion to uK from Jy/sr for freq = %.0f GHz" % (obsFreqGHz)
        data=mapTools.convertToY(data, obsFrequencyGHz = obsFreqGHz)
        if saveAsY == True:
            astImages.saveFITS("yc_inputSimMap.fits", data, wcs)
            
        count=0
        for obj in inputCatalog:            
            # Progress report
            count=count+1
            tenPercent=len(inputCatalog)/10
            for j in range(0,11):
                if count == j*tenPercent:
                    print "... "+str(j*10)+"% complete ..."
            # Only measure flux if this object passes Y, z limit cuts (otherwise this will take forever!)
            if obj['z'] > zRange[0] and obj['z'] < zRange[1] and obj['fluxFromCatalog_arcmin2'] > YLimit:
                x0, y0=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
                ra0, dec0=[obj['RADeg'], obj['decDeg']]
                ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
                xLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
                yLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
                localDegPerPix=astCoords.calcAngSepDeg(obj['RADeg'], obj['decDeg'], ra1, dec1) 
                RvirDeg=math.degrees(math.atan(obj['RvirMpc']/astCalc.da(obj['z'])))
                RvirPix=RvirDeg/localDegPerPix                
                flux_arcmin2=photometry.objectFluxInAperture(obj, apertureRadiusArcmin, data, wcs)
                obj['fixedApertureFluxFromInputSimMap_arcmin2']=flux_arcmin2
            
        # Pickle results for speed
        pickleFile=file(pickleFileName, "w")
        pickler=pickle.Pickler(pickleFile)
        pickler.dump(inputCatalog)
                    
    return inputCatalog    
    
#------------------------------------------------------------------------------------------------------------
def checkInputSimCatalogFluxes(inputSimCatFileName, inputSimMapFileName, YLimit = None, zRange = None,
                               outDir = ".", convertFromJySrToMicroK = False, obsFreqGHz = 148):
    """This measures fluxes directly from the input sim map, and compares to catalog, doing a least squares 
    fit. Saves out plot under 'inputSimChecks'.
        
    """

    print ">>> Getting fluxes from noiseless input sim catalog ... "
    
    inputCatalog=parseInputSimCatalog(inputSimCatFileName, YLimit = YLimit, zRange = zRange)

    # Handy DS9 .reg file
    catalog2DS9(inputCatalog, outDir+os.path.sep+"inputSimCatalog.reg", idKeyToUse = 'name', \
                    color = "yellow")

    # Measure fluxes in actual input map from directly adding up flux within given Rvir    
    img=pyfits.open(inputSimMapFileName)
    wcs=astWCS.WCS(inputSimMapFileName)
    data=img[0].data
    if convertFromJySrToMicroK == True: # from Jy/sr
        if obsFreqGHz == 148:
            print ">>> Converting from Jy/sr to micro kelvins assuming obsFreqGHz = %.0f" % (obsFreqGHz)
            data=(data/1.072480e+09)*2.726*1e6
        else:
            raise Exception, "no code added to support conversion to uK from Jy/sr for freq = %.0f GHz" % (obsFreqGHz)
    data=mapTools.convertToY(data, obsFrequencyGHz = obsFreqGHz)

    for obj in inputCatalog:
        
        x0, y0=wcs.wcs2pix(obj['RADeg'], obj['decDeg'])
        ra0, dec0=[obj['RADeg'], obj['decDeg']]
        ra1, dec1=wcs.pix2wcs(x0+1, y0+1)
        xLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
        yLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
        localDegPerPix=astCoords.calcAngSepDeg(obj['RADeg'], obj['decDeg'], ra1, dec1)   
        arcmin2PerPix=xLocalDegPerPix*yLocalDegPerPix*60.0**2
        RvirDeg=math.degrees(math.atan(obj['RvirMpc']/astCalc.da(obj['z'])))
        RvirPix=RvirDeg/localDegPerPix
        clip=astImages.clipImageSectionPix(data, x0, y0, int(round(RvirPix*3)))
        
        # This automatically chucks out objects too close to edge of map
        if clip.shape[1] == clip.shape[0]:
            obj['RvirPix']=RvirPix
            obj['RvirArcmin']=RvirDeg*60.0
            x=np.arange(-clip.shape[1]/2, clip.shape[1]/2, dtype=float)*localDegPerPix
            y=np.arange(-clip.shape[0]/2, clip.shape[0]/2, dtype=float)*localDegPerPix
            rDeg=np.sqrt(x**2+y**2)            
            insideApertureMask=np.less(rDeg, RvirDeg)
            obj['fluxFromInputSimMap_arcmin2']=clip[insideApertureMask].sum()*arcmin2PerPix
    
    # Fit for not understood offset between directly measured input map fluxes and the catalog
    fluxesFromMap=[]
    fluxesFromCatalog=[]
    for obj in inputCatalog:
        if 'fluxFromInputSimMap_arcmin2' in obj.keys():
            fluxesFromMap.append(obj['fluxFromInputSimMap_arcmin2'])
            fluxesFromCatalog.append(obj['fluxFromCatalog_arcmin2'])
    fluxesFromMap=np.array(fluxesFromMap)
    fluxesFromCatalog=np.array(fluxesFromCatalog)
    
    # Clipped OLS fit
    res=np.zeros(fluxesFromMap.shape)
    sigma=1e6
    for i in range(10):
        fitData=[]
        for m, c, r in zip(fluxesFromMap, fluxesFromCatalog, res):
            if abs(r) < 2.0*sigma:
                fitData.append([c, m])
        fit=astStats.OLSFit(fitData)
        res=(fluxesFromCatalog*fit['slope']+fit['intercept'])-fluxesFromMap
        sigma=np.std(res)

    # Plot    
    fig=plt.figure(num=6, figsize=(8,8))
    fig.canvas.set_window_title('Noiseless Input Sim Catalog vs. Map Flux Comparison')
    fitRange=np.arange(0, 0.02, 0.001)
    fitLine=fit['slope']*fitRange+fit['intercept']
    plt.plot(fluxesFromCatalog, fluxesFromMap, 'ro')
    plt.plot(fitRange, fitLine, 'b--', label='fit = %.6f*x + %.6f' % (fit['slope'], fit['intercept']))
    #plt.plot(np.arange(0, 0.02, 0.001), np.arange(0, 0.02, 0.001), 'b--')
    plt.xlabel("Y input sim catalog (arcmin$^2$)")
    plt.ylabel("Y input sim map (arcmin$^2$)")
    plt.xlim(0, 0.02)
    plt.ylim(0, 0.02)
    plt.legend()
    plt.savefig(outDir+os.path.sep+"noiselessInputSim_catalogVsMapFluxes.png")
        
    # Put corrected fluxes into our input catalog
    for obj in inputCatalog:
        if 'fluxFromInputSimMap_arcmin2' in obj.keys():
            obj['correctedFlux_arcmin2']=fit['slope']*obj['fluxFromCatalog_arcmin2']+fit['intercept']

    return inputCatalog

#-------------------------------------------------------------------------------------------------------------
def fakeSourceSims(fakeSourceSimOptions, unfilteredMapsDictList, filtersList, detectionThresholdSigma, 
                   detectionMinObjPix, rejectBorderPix, minSNToInclude, photometryOptions, diagnosticsDir):
    """For each of the unfilteredMaps, inserts fake clusters with a range of different Ys, sizes etc.. and
    runs the source finding and photometry over them. Makes plots of the completeness as a function of Y, 
    size and plots of flux recovery (input fake source flux vs. recovered flux). All output is stored in
    a subdir under diagnosticsDir/
    
    """
    
    print ">>> Running completeness & flux recovery sims ..."
    
    outDir=diagnosticsDir+os.path.sep+"fakeSourceSims"
    if os.path.exists(outDir) == False:
        os.makedirs(outDir)
    
    # Make a bunch of fake maps to filter
    # We'll randomly populate each fake map with sources drawn from the list of allowed scales, deltaTs
    # rather than doing a separate run with each set of source parameters
    numRuns=fakeSourceSimOptions['numRuns']
    sourcesPerRun=fakeSourceSimOptions['sourcesPerRun']
    scalesArcminList=fakeSourceSimOptions['scalesArcminList']
    deltaTList=fakeSourceSimOptions['deltaTList']
    
    # Selection function stuff - how this works:
    #
    # selFnDict -> profileType keys -> propertiesToTrack keys -> recovered, total per run keys
    #
    # e.g. selFnDict -> 'betaProfile' -> 'deltaT'
    #                                 -> 'scaleArcmin'
    selFnDict={}
    propertiesToTrack=['deltaT', 'scaleArcmin']
    propertyPlotLabels=['$\Delta T_c$ ($\mu$K)', '$\\theta_c$ (arcmin)']
    propertyValuesLists=[deltaTList, scalesArcminList]
    
    for profileType in fakeSourceSimOptions['profilesList']:
        
        if profileType == 'betaModel':
            insertSourceIntoMap=insertBetaModelIntoMap
        elif profileType == 'arnaudModel':
            insertSourceIntoMap=insertArnaudModelIntoMap
        elif profileType == 'projectedNFWModel':
            insertSourceIntoMap=insertProjectedNFWModelIntoMap
        else:
            raise Exception, "didn't understand profileType"
        
        # Set up selection function storage
        if profileType not in selFnDict.keys():
            selFnDict[profileType]={}
            for prop, valuesList in zip(propertiesToTrack, propertyValuesLists):
                selFnDict[profileType][prop]={'recoveredByRun': np.zeros([numRuns, len(valuesList)], dtype=float), 
                                              'totalByRun': np.zeros([numRuns, len(valuesList)], dtype=float),
                                              'valuesList': valuesList}        
            
        inputYArcmin2=[]
        recoveredYArcmin2=[]
        inputDeltaTc=[]
        recoveredDeltaTc=[]
        for run in range(numRuns):
            
            t0=time.time()
            print "--> Run: %d" % (run+1)
            
            label='%s_run_%d' % (profileType, run+1)
            outFileName=outDir+os.path.sep+"fakeMap_%s.fits" % (label)
            
            # Load in map, do minimal pre-processing (trimming etc.). Assuming here 148 GHz only
            fakeMapDict={}
            fakeMapDict['RADecSection']=unfilteredMapsDictList[0]['RADecSection']
            fakeMapDict['obsFreqGHz']=unfilteredMapsDictList[0]['obsFreqGHz']
            fakeMapDict['units']=unfilteredMapsDictList[0]['units']
            fakeMapDict['beamFWHMArcmin']=unfilteredMapsDictList[0]['beamFWHMArcmin']
            fakeMapDict['pointSourceRemoval']=None
            fakeMapDict['mapFileName']=unfilteredMapsDictList[0]['mapFileName']
            fakeMapDict['weightsFileName']=unfilteredMapsDictList[0]['weightsFileName']
            mapTools.preprocessMapDict(fakeMapDict, diagnosticsDir = diagnosticsDir)
            
            # Stops cython complaining
            fakeMapDict['data']=np.array(fakeMapDict['data'], dtype=np.float64) 
            
            # Uncomment out below if want to do noiseless sanity checks on e.g. recovered Y etc.
            fakeMapDict['data']=np.zeros(fakeMapDict['data'].shape, dtype=np.float64) 
             
            # Generate fake source catalog
            fakeInputCatalog=[]
            xMin=0+rejectBorderPix*2
            xMax=fakeMapDict['data'].shape[1]-rejectBorderPix*2
            yMin=0+rejectBorderPix*2
            yMax=fakeMapDict['data'].shape[0]-rejectBorderPix*2
            xs=np.random.randint(xMin, xMax, sourcesPerRun)
            ys=np.random.randint(yMin, yMax, sourcesPerRun)
            wcsCoords=fakeMapDict['wcs'].pix2wcs(xs, ys)
            wcsCoords=np.array(wcsCoords)
            RAs=wcsCoords[:, 0]
            decs=wcsCoords[:, 1]
            for k in range(sourcesPerRun):
                objDict={}
                objDict['name']="fake_"+catalogTools.makeACTName(RAs[k], decs[k]).replace(" ", "_")
                objDict['RADeg']=RAs[k]
                objDict['decDeg']=decs[k]
                objDict['x']=xs[k]
                objDict['y']=ys[k]
                fakeInputCatalog.append(objDict)
            
            # Insert fake catalog sources into map, add source properties to catalog (integrated Y etc.)
            for obj in fakeInputCatalog:
                deltaT=deltaTList[np.random.randint(len(deltaTList))]
                scaleArcmin=scalesArcminList[np.random.randint(len(scalesArcminList))]
                obj['deltaT']=deltaT
                obj['scaleArcmin']=scaleArcmin
                fakeMapDict['data'], sourceProperties=insertSourceIntoMap(obj, deltaT, scaleArcmin, 
                                                                          fakeMapDict, photometryOptions)
                for key in sourceProperties.keys():
                    obj[key]=sourceProperties[key]
            
            # Complete map pre-processing - i.e. do same background subtraction, point source removal
            maskFileName=diagnosticsDir+os.path.sep+"psMask_%d.fits" % (fakeMapDict['obsFreqGHz'])
            maskingType=unfilteredMapsDictList[0]['pointSourceRemoval']['masking']
            fakeMapDict['data']=mapTools.applyPointSourceMask(maskFileName, fakeMapDict['data'], 
                                                              fakeMapDict['wcs'], mask = maskingType)
                                                              
            # Save here before backgroundSubtraction, in case we want to pass on to Matthew
            astImages.saveFITS(outFileName, fakeMapDict['data'], fakeMapDict['wcs'])
            fakeMapDict['mapFileName']=outFileName
            keysToWrite=['name', 'RADeg', 'decDeg', 'deltaT', 'scaleArcmin', 'YArcmin2']
            keyFormats=['%s', '%.6f', '%.6f', '%.3f', '%.1f', '%.10f']
            extraHeaderText="# Profile = %s, YArcmin2 measured in %.1f arcmin radius circular aperture\n" \
                             % (profileType, photometryOptions['apertureRadiusArcmin'])
            catalogTools.writeCatalog(fakeInputCatalog, outFileName.replace(".fits", ".csv"), keysToWrite, 
                                      keyFormats, [], extraHeaderText = extraHeaderText)
            catalogTools.catalog2DS9(fakeInputCatalog, outFileName.replace(".fits", ".reg"))

            # Ok, now do backgroundSubtraction
            if unfilteredMapsDictList[0]['backgroundSubtraction'] == True:
                fakeMapDict['data']=mapTools.subtractBackground(fakeMapDict['data'], fakeMapDict['wcs'])
            
            fakeInputMapsDictList=[fakeMapDict]
        
            # Filter maps, detect objects, measure fluxes, merge catalogs - in same way as in nemo script
            if os.path.exists(outDir+os.path.sep+"filteredMaps") == True:
                os.system("rm -r %s" % (outDir+os.path.sep+"filteredMaps"))
                os.system("rm -r %s" % (outDir+os.path.sep+"diagnostics"))
            imageDict=mapTools.filterMaps(fakeInputMapsDictList, filtersList, rootOutDir = outDir)
            photometry.findObjects(imageDict, threshold = detectionThresholdSigma, 
                                   minObjPix = detectionMinObjPix, rejectBorder = rejectBorderPix)            
            photometry.measureFluxes(imageDict, photometryOptions, diagnosticsDir, 
                                     unfilteredMapsDict = fakeInputMapsDictList)
            catalogTools.mergeCatalogs(imageDict)
            catalogTools.makeOptimalCatalog(imageDict, minSNToInclude)
                        
            # Match output catalog against input to get recovered fraction
            outputCatalog=imageDict['optimalCatalog']
            simpleCatalogMatch(fakeInputCatalog, outputCatalog, fakeMapDict['beamFWHMArcmin']*2)
             
            # Recovered fraction as fn. of scale
            for obj in fakeInputCatalog:
                for prop in propertiesToTrack:
                    valuesList=selFnDict[profileType][prop]['valuesList']
                    index=valuesList.index(obj[prop])
                    selFnDict[profileType][prop]['totalByRun'][run, index]+=1
                    if obj['recovered'] == True:
                        selFnDict[profileType][prop]['recoveredByRun'][run, index]+=1
                        if prop == propertiesToTrack[0]:    # avoid double counting
                            inputYArcmin2.append(obj['YArcmin2'])
                            recoveredYArcmin2.append(obj['recoveredMatch']['flux_arcmin2'])
                            inputDeltaTc.append(obj['deltaT'])
                            recoveredDeltaTc.append(obj['recoveredMatch']['deltaT_c'])
        
            t1=time.time()
            print "... time taken for fake source insertion & recovery run = %.3f sec" % (t1-t0)
            
        # Bring all results from sim runs together for a given profileType
        for prop in propertiesToTrack:
            fraction=selFnDict[profileType][prop]['recoveredByRun']/selFnDict[profileType][prop]['totalByRun']
            mean=fraction.mean(axis=0)
            stderr=fraction.std(axis=0)/np.sqrt(fraction.shape[0])
            selFnDict[profileType][prop]['meanRecoveredFraction']=mean
            selFnDict[profileType][prop]['stderrRecoveredFraction']=stderr            
        
        # Plot for each property
        plt.close()
        for prop, plotLabel in zip(propertiesToTrack, propertyPlotLabels):
            fraction=selFnDict[profileType][prop]['meanRecoveredFraction']
            errFraction=selFnDict[profileType][prop]['stderrRecoveredFraction']
            valuesList=selFnDict[profileType][prop]['valuesList']
            plt.plot(valuesList, fraction, 'b-')
            plt.errorbar(valuesList, fraction, yerr=errFraction, fmt='r.')
            plt.xlabel(plotLabel)
            plt.ylabel("Recovered Fraction")
            plt.ylim(0, 1.05)
            plt.title("Fake Source Profile = %s" % (profileType))
            outFileName=outDir+os.path.sep+"recoveredFraction_%s_%s.png" % (profileType, prop)
            plt.savefig(outFileName)
            plt.close()
        
        # Y recovery plot - residual
        inputYArcmin2=np.array(inputYArcmin2)
        recoveredYArcmin2=np.array(recoveredYArcmin2)
        diff=inputYArcmin2-recoveredYArcmin2
        norm=1e-3
        labelString="mean(Y$_{in}$-Y$_{out}$) = %.3f $\\times$ 10$^{-3}$ arcmin$^2$\n$\sigma$(Y$_{in}$-Y$_{out}$) = %.3f $\\times$ 10$^{-3}$ arcmin$^2$" % ((diff/norm).mean(), (diff/norm).std())
        plt.figure(figsize=(10, 8))
        plt.plot(inputYArcmin2, diff, 'r.')
        plt.plot(np.linspace(inputYArcmin2.min()-1, inputYArcmin2.max()+1, 3), [0]*3, 'k--')
        plt.xlabel("input Y(r<%d') (arcmin$^2$)" % (photometryOptions['apertureRadiusArcmin']))
        plt.ylabel("input Y(r<%d') - recovered Y(r<%d') (arcmin$^2$)" % (photometryOptions['apertureRadiusArcmin'], photometryOptions['apertureRadiusArcmin']))
        ax=plt.gca()
        plt.text(0.03, 0.03, labelString, transform = ax.transAxes, fontdict = {"size": 14, "linespacing" : 1.5})
        plt.title("Fake Source Profile = %s" % (profileType))
        plt.xlim(0, inputYArcmin2.max()*1.1)
        plt.ylim(-inputYArcmin2.max(), inputYArcmin2.max())
        outFileName=outDir+os.path.sep+"YRecovery_residual_%s.png" % (profileType)
        plt.savefig(outFileName)
        plt.close()
        
        # Y recovery plot - correlation
        plt.figure(figsize=(10, 8))
        oneToOneRange=np.linspace(inputYArcmin2.min()*0.5, inputYArcmin2.max()+1, 5)
        plt.plot(inputYArcmin2, recoveredYArcmin2, 'r.')
        plt.plot(oneToOneRange, oneToOneRange, 'k--')
        plt.loglog()
        plt.xlabel("input Y(r<%d') (arcmin$^2$)" % (photometryOptions['apertureRadiusArcmin']))
        plt.ylabel("recovered Y(r<%d') (arcmin$^2$)" % (photometryOptions['apertureRadiusArcmin']))
        plt.title("Fake Source Profile = %s" % (profileType))
        plt.xlim(inputYArcmin2.min()*0.5, inputYArcmin2.max()*1.5)
        plt.ylim(inputYArcmin2.min()*0.5, inputYArcmin2.max()*1.5)
        outFileName=outDir+os.path.sep+"YRecovery_correlation_%s.png" % (profileType)
        plt.savefig(outFileName)
        plt.close()

        # Delta T recovery plot - residual
        inputDeltaTc=np.array(inputDeltaTc)
        recoveredDeltaTc=np.array(recoveredDeltaTc)
        diff=inputDeltaTc-recoveredDeltaTc
        labelString="mean($\Delta$T$_{c(in)}$-$\Delta$T$_{c(out)}$) = %.3f $\mu$K\n$\sigma$($\Delta$T$_{c(in)}$-$\Delta$T$_{c(out)}$) = %.3f $\mu$K" % ((diff).mean(), (diff).std())
        plt.figure(figsize=(10, 8))
        plt.plot(inputDeltaTc, diff, 'r.')
        plt.plot(np.linspace(inputDeltaTc.min()-1000, inputDeltaTc.max()+1000, 3), [0]*3, 'k--')
        plt.xlabel("input $\Delta$T$_c$ ($\mu$K)")
        plt.ylabel("input $\Delta$T$_c$ - recovered $\Delta$T$_c$ ($\mu$K)")
        ax=plt.gca()
        plt.text(0.03, 0.03, labelString, transform = ax.transAxes, fontdict = {"size": 14, "linespacing" : 1.5})
        plt.title("Fake Source Profile = %s" % (profileType))
        plt.xlim(inputDeltaTc.min()*1.1, inputDeltaTc.max()*1.1)
        plt.ylim(-(abs(inputDeltaTc).max()*1.1), abs(inputDeltaTc).max()*1.1)
        outFileName=outDir+os.path.sep+"DeltaTc_residual_%s.png" % (profileType)
        plt.savefig(outFileName)
        plt.close()        
        
        # Also, work out fudge factor for Y recovery here
        # Insert this into photometry.measureApertureFluxes - obviously we have to turn this off there first
        # if we want to refit this.
        mask=np.greater(recoveredYArcmin2, 0)
        #slope, intercept, blah1, blah2, blah3=stats.linregress(recoveredYArcmin2[mask], inputYArcmin2[mask])
        slope, intercept, blah1, blah2, blah3=stats.linregress(np.log10(recoveredYArcmin2[mask]), np.log10(inputYArcmin2[mask]))
        #rec2=np.power(10.0, np.log10(recoveredYArcmin2)*slope+intercept)
        
        # Save Y data and fit for this profileType
        outFileName=outDir+os.path.sep+"YRecovery_%s.npz" % (profileType)
        np.savez(outFileName, inputYArcmin2, recoveredYArcmin2)
        
    # Save sel. fn. as pickle
    # May want something to convert this into more readable format
    pickleFileName=outDir+os.path.sep+"selFnDict.pickle"
    if os.path.exists(pickleFileName) == True:
        os.remove(pickleFileName)
    pickleFile=file(pickleFileName, "w")
    pickler=pickle.Pickler(pickleFile)
    pickler.dump(selFnDict)
    pickleFile.close()
    
    print "Done, check plots etc., what's going on with recovered Y"    
    ipshell()
    sys.exit()
            
#-------------------------------------------------------------------------------------------------------------
def insertBetaModelIntoMap(objDict, deltaT, scaleArcmin, mapDict, photometryOptions):
    """Inserts a beta model into the map. Adds source properties to the object.
    
    Returns updated map data with source inserted, and sourceProperties dictionary with intergrated Y etc.
    for the objDict.

    """
    
    sourceProperties={} # store things like integrated Y in here
    
    rDegMap=nemoCython.makeDegreesDistanceMap(mapDict['data'], mapDict['wcs'], objDict['RADeg'], 
                                              objDict['decDeg'], (10*scaleArcmin)/60.0)

    # beta fixed, for now
    beta=0.86
    rArcmin=np.linspace(0.0, 60.0, 5000)
    smoothArcmin=scaleArcmin
    profile1d=(1.0+(rArcmin/scaleArcmin)**2)**((1.0-3.0*beta)/2)
    
    # Scale to size of input central decrement, before beam smoothing
    profile1d=profile1d*deltaT
        
    # Apply beam as Gaussian filter to profile
    #beamSigma=mapDict['beamFWHMArcmin']/np.sqrt(8.0*np.log(2.0))            
    #beamSigmaPix=beamSigma/(rArcmin[1]-rArcmin[0])
    #profile1d=ndimage.gaussian_filter1d(profile1d, beamSigmaPix)
    
    # Truncate beyond 10 times core radius
    mask=np.greater(rArcmin, 10.0*scaleArcmin)
    profile1d[mask]=0.0
    
    # Turn 1d profile into 2d
    rDeg=rArcmin/60.0
    r2p=interpolate.interp1d(rDeg, profile1d, bounds_error=False, fill_value=0.0)
    profile2d=np.zeros(rDegMap.shape)
    mask=np.less(rDegMap, 1000)
    profile2d[mask]=r2p(rDegMap[mask])
    
    mapDict['data']=mapDict['data']+profile2d
    
    # What is the integrated Y within the aperture we're using for measuring Ys?
    mask=np.less(rDegMap, photometryOptions['apertureRadiusArcmin']/60.0)    
    sumPix=mapTools.convertToY(profile2d[mask], obsFrequencyGHz = mapDict['obsFreqGHz']).sum()
    ra0=objDict['RADeg']
    dec0=objDict['decDeg']
    x, y=mapDict['wcs'].wcs2pix(ra0, dec0)
    ra1, dec1=mapDict['wcs'].pix2wcs(x+1, y+1)    
    xLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
    yLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
    arcmin2PerPix=xLocalDegPerPix*yLocalDegPerPix*60.0**2
    YArcmin2=sumPix*arcmin2PerPix
    sourceProperties['YArcmin2']=YArcmin2
    
    return [mapDict['data'], sourceProperties]

#-------------------------------------------------------------------------------------------------------------
def insertProjectedNFWModelIntoMap(objDict, deltaT, scaleArcmin, mapDict, photometryOptions):
    """Inserts a projected 2d NFW profile (see Koester et al. 2007 and references therein)
    
    """
    
    rs=0.3 # 0.3 Mpc is sensible if R200 = 1.5 Mpc, c = 5, plus Giodini gets that. Koester uses 150 kpc
    r=np.linspace(0.0, 10.0, 4000)
    x=r/rs
    fx=np.zeros(x.shape)
    mask=np.greater(x, 1)
    fx[mask]=1-(2.0/np.sqrt(x[mask]**2-1))*np.arctan(np.sqrt((x[mask]-1)/(x[mask]+1)))
    mask=np.less(x, 1)
    fx[mask]=1-(2.0/np.sqrt(1-x[mask]**2))*np.arctanh(np.sqrt((1-x[mask])/(x[mask]+1)))
    mask=np.equal(x, 1)
    fx[mask]=0
    mask=np.greater(x, 20)
    fx[mask]=0
    sigmax=np.zeros(r.shape)
    mask=np.greater(r, rs)
    sigmax[mask]=(2*rs*fx[mask])/(x[mask]**2-1)   # ignoring rho_s, which is arbitrary
    
    # Fit power law for extrapolating in centre (NFW profile undefined here)
    mask=np.logical_and(np.greater(r, rs), np.less(r, rs*3))
    deg=1
    p=np.polyfit(np.log10(r[mask]), np.log10(sigmax[mask]), 1)
    fittedFunc=np.power(10, p[0]*np.log10(r)+p[1])
    sigmax[np.less(r, rs)]=fittedFunc[np.less(r, rs)]
    sigmax[0]=sigmax[1] # centre is still infinite
    sigmax=sigmax/sigmax.max()
    tckSigmax=interpolate.splrep(r, sigmax) # Note we defined interpolator in terms or r NOT x here!

    #plt.plot(r, fittedFunc, 'k--')
    
    ipshell()
    sys.exit()
    
#-------------------------------------------------------------------------------------------------------------
def insertArnaudModelIntoMap(objDict, deltaT, scaleArcmin, mapDict, photometryOptions):
    """Inserts an Arnaud (GNFW) model into the map. Adds source properties to the object.
    
    Returns updated map data with source inserted, and sourceProperties dictionary with intergrated Y etc.
    for the objDict.

    """
    
    sourceProperties={} # store things like integrated Y in here
    
    rDegMap=nemoCython.makeDegreesDistanceMap(mapDict['data'], mapDict['wcs'], objDict['RADeg'], 
                                              objDict['decDeg'], (10*scaleArcmin)/60.0)
   
    # The easy way - use Matthew's saved Arnaud GNFW profile and scale it according to scale arcmin
    print "Arnaud model!"
    ipshell()
    sys.exit()
 
    r=np.linspace(0.0001, 3, 1000)    # in Mpc
    r500=1.0                        # in Mpc
    
    x=r/r500
        
    P0=8.403
    c500=1.177
    gamma=0.3081
    alpha=1.0510
    beta=5.4905
    
    # dimensionlessP _is_ just the gNFW profile
    dimensionlessP=P0/(np.power(c500*x, gamma)*np.power((1+np.power(c500*x, alpha)), (beta-gamma)/alpha))
        
    # Things get physical here
    z=0.3   # redshift
    M500=5e14
    
    alphaP=0.12
    alphaPPrime=0.10-(alphaP+0.10)*(((x/0.5)**3)/(1.+(x/0.5)**3))
    
    Pr=1.65e-3*astCalc.Ez(z)*np.power((M500/3e14), 2.0/3.0+alphaP+alphaPPrime)*dimensionlessP
    
    # Turn from radial profile into projected radial cylindrical profile
    tck_Pr=interpolate.splrep(r, Pr)
    rCyn=np.zeros(r.shape)+r
    profCyn=[]
    dr=rCyn[1]-rCyn[0]
    for rc in rCyn:
        dimensionlessP=P0/(np.power(c500*x, gamma)*np.power((1+np.power(c500*x, alpha)), (beta-gamma)/alpha))

        dr = r[1] - r[0]
    y0 = array([profile.get((_x0**2+r**2)**0.5).sum() for _x0 in x0]) / dr
    
    # Convert to angular coords, do cylindrical integral
    rDeg=np.degrees(np.arctan(r/astCalc.da(z)))
    r2p=interpolate.interp1d(rDeg, Pr, bounds_error=False, fill_value=0.0)
    profile2d=np.zeros(rDegMap.shape)
    mask=np.less(rDegMap, 1000)
    profile2d[mask]=r2p(rDegMap[mask])
    
    
#-------------------------------------------------------------------------------------------------------------
def simpleCatalogMatch(primary, secondary, matchRadiusArcmin):
    """Simple catalog matching, for finding which fake objects we recovered
    
    Adds 'recovered' key to primary catalog, in place.
    
    """
        
    xMatchRadiusDeg=matchRadiusArcmin/60.0
    sobjRAs=[]
    sobjDecs=[]
    for sobj in secondary:
        sobjRAs.append(sobj['RADeg'])
        sobjDecs.append(sobj['decDeg'])
    sobjRAs=np.array(sobjRAs)
    sobjDecs=np.array(sobjDecs)
    for pobj in primary:   
        pobj['recovered']=False
        rMin=1e6
        match=None
        rs=astCoords.calcAngSepDeg(pobj['RADeg'], pobj['decDeg'], sobjRAs, sobjDecs)
        rMin=rs.min()
        rMinIndex=np.equal(rs, rMin).nonzero()[0][0]
        if rMin < xMatchRadiusDeg:
            match=secondary[rMinIndex]
        if match != None:
            pobj['recovered']=True
            pobj['recoveredMatch']=match

#-------------------------------------------------------------------------------------------------------------
def generateCMBOnlySim(parDict):
    """Generates a signal-only CMB sky, convolved with the beam.
    
    Returns beam-smoothed map (2d array) and wcs
    
    """
    
    # CAMB power spec with Planck 2015 parameters (ish)
    tab=atpy.Table().read(nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"planck_lensedCls.dat", format = 'ascii')
    ell=tab['L']
    DEll=tab['TT']
    Cl_uK2=(DEll*2*np.pi)/(ell*(ell+1))
    
    # NOTE: assuming the first map is the only one, used for getting dimensions... this will break if this isn't the case 
    # See here for flipper docs: http://www.hep.anl.gov/sdas/flipperDocumentation/map_manipulations.html
    img=pyfits.open(parDict['unfilteredMaps'][0]['mapFileName'])
    wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
    d=np.zeros(img[0].data.shape, dtype = float)
    grfMap=liteMap.liteMapFromDataAndWCS(d, wcs)
    grfMap.fillWithGaussianRandomField(ell, Cl_uK2)
    
    # Convolve with beam - NOTE: beam map only 4 x HWHM on a side (otherwise convolution very slow)
    beamData=np.loadtxt(parDict['unfilteredMaps'][0]['beamFileName']).transpose()
    profile1d=beamData[1]
    rArcmin=beamData[0]*60.0
    HWHMDeg=rArcmin[np.where(profile1d <= 0.5)[0][0]]/60.
    radiusDeg=HWHMDeg*4
    RADeg, decDeg=grfMap.wcs.getCentreWCSCoords()
    clip=astImages.clipUsingRADecCoords(grfMap.data, grfMap.wcs, RADeg+radiusDeg, RADeg-radiusDeg, decDeg-radiusDeg, decDeg+radiusDeg)
    RADeg, decDeg=clip['wcs'].getCentreWCSCoords()
    degreesMap=nemoCython.makeDegreesDistanceMap(clip['data'], clip['wcs'], RADeg, decDeg, 0.5)
    beamSignalMap, inputSignalProperties=makeBeamModelSignalMap(parDict['unfilteredMaps'][0]['obsFreqGHz'], degreesMap, 
                                                                clip['wcs'], parDict['unfilteredMaps'][0]['beamFileName'])
    beamSmoothedMap=ndimage.convolve(grfMap.data, beamSignalMap.data/beamSignalMap.data.sum()) 
    
    return beamSmoothedMap, grfMap.wcs

#-------------------------------------------------------------------------------------------------------------
def estimateContaminationFromSkySim(imageDict, parDictFileName, numSkySims, diagnosticsDir):
    """Run the whole filtering set up again, on sky simulation, with noise, that we generate here.
    
    Turns out we need to run many realisations, as results can vary by a lot.
    
    We will want to combine the output from this with the inverted maps test (which is quicker and easier).
    
    Writes a plot and a .fits table to the diagnostics dir.
    
    Runs over both SNR and fixed_SNR values.
    
    Returns a dictionaries containing the results
    
    """
    
    resultsList=[]
    for i in range(numSkySims):
        
        print "... sky sim %d/%d ..." % (i+1, numSkySims)
        simParDict=actDict.ACTDict()
        simParDict.read_from_file(parDictFileName)
            
        simFileName=diagnosticsDir+os.path.sep+"signalOnlyCMBSim.fits"
        print "... generating CMB-only sim ..."
        simMap, wcs=generateCMBOnlySim(simParDict)
        #astImages.saveFITS(simFileName, simMap, wcs)
        
        # Add noise...
        # Rather than using the 'noiseBoostFactor', and issuing a warning about inv var only weight maps, we could
        # just solve for the noiseBoostFactor needed to make the RMS map match the actual one
        print "... adding noise (WARNING: for this to work, the weight map must be an inverse variance map) ..."
        iVarImg=pyfits.open(simParDict['unfilteredMaps'][0]['weightsFileName'])
        iVarMap=iVarImg[0].data
        footprint=np.where(iVarMap != 0)
        blank=np.where(iVarMap == 0)
        noiseMap=np.zeros(iVarMap.shape)
        noiseMap[footprint]=np.random.normal(0, np.sqrt(1/iVarMap[footprint]))
        noiseMap=noiseMap*simParDict['noiseBoostFactor'] # test if reasonable: RMSMap of real data and sim should agree
        simMap[blank]=0.
        simMap=simMap+noiseMap
        astImages.saveFITS(simFileName.replace(".fits", "_withNoise.fits"), simMap, wcs)

        # After we get out sky sim map into this, we want to do the filtering again...
        # NOTE: redirect rootOutDir here as well
        simParDict['unfilteredMaps'][0]['mapFileName']=simFileName.replace(".fits", "_withNoise.fits")
        rootOutDir=diagnosticsDir+os.path.sep+"skySim"
        if os.path.exists(rootOutDir) == True:
            mapFileNames=glob.glob(rootOutDir+os.path.sep+"filteredMaps"+os.path.sep+"*.fits")
            for m in mapFileNames:
                os.remove(m)
        simImageDict=mapTools.filterMaps(simParDict['unfilteredMaps'], simParDict['mapFilters'], rootOutDir = rootOutDir)
            
        # Below here is same as inverted maps right now....
        # If we have makeDS9Regions = True here, we overwrite the existing .reg files from when we ran on the non-inverted maps
        photometry.findObjects(simImageDict, threshold = simParDict['thresholdSigma'], minObjPix = simParDict['minObjPix'],
                            rejectBorder = simParDict['rejectBorder'], diagnosticsDir = diagnosticsDir,
                            invertMap = False, makeDS9Regions = True)    

        # For fixed filter scale
        # Adds fixed_SNR values to catalogs for all maps
        photometryOptions=simParDict['photometryOptions']
        if 'photFilter' in photometryOptions.keys():
            photFilter=photometryOptions['photFilter']
        else:
            photFilter=None
        if photFilter != None:
            photometry.getSNValues(simImageDict, SNMap = 'file', prefix = 'fixed_', template = photFilter, invertMap = False)   
            SNRKeys=['SNR', 'fixed_SNR']
        else:
            SNRKeys=['SNR']
        
        minSNToIncludeInOptimalCatalog=simParDict['catalogCuts']
        catalogTools.mergeCatalogs(simImageDict)
        catalogTools.makeOptimalCatalog(simImageDict, minSNToIncludeInOptimalCatalog)
            
        contaminTabDict=estimateContamination(simImageDict, imageDict, SNRKeys, 'skySim', diagnosticsDir)
        resultsList.append(contaminTabDict)
    
    # Average results
    avContaminTabDict={}
    for k in resultsList[0].keys():
        avContaminTabDict[k]=atpy.Table()
        for kk in resultsList[0][k].keys():
            avContaminTabDict[k].add_column(atpy.Column(np.zeros(len(resultsList[0][k])), kk))
            for i in range(len(resultsList)):
                avContaminTabDict[k][kk]=avContaminTabDict[k][kk]+resultsList[i][k][kk]
            avContaminTabDict[k][kk]=avContaminTabDict[k][kk]/float(len(resultsList))

    for k in avContaminTabDict.keys():
        fitsOutFileName=diagnosticsDir+os.path.sep+"%s_contaminationEstimate.fits" % (k)
        if os.path.exists(fitsOutFileName) == True:
            os.remove(fitsOutFileName)
        contaminTab=avContaminTabDict[k]
        contaminTab.write(fitsOutFileName)
        
    return avContaminTabDict

#-------------------------------------------------------------------------------------------------------------
def estimateContaminationFromInvertedMaps(imageDict, thresholdSigma, minObjPix, rejectBorder, 
                                          minSNToIncludeInOptimalCatalog, photometryOptions, diagnosticsDir, findCenterOfMass = True):
    """Run the whole filtering set up again, on inverted maps.
    
    Writes a DS9. reg file, which contains only the highest SNR contaminants (since these
    are most likely to be associated with artefacts in the map - e.g., point source masking).
    
    Writes a plot and a .fits table to the diagnostics dir.
    
    Runs over both SNR and fixed_SNR values.
    
    Returns a dictionaries containing the results
    
    """
    
    invertedDict={}
    ignoreKeys=['optimalCatalog', 'mergedCatalog']
    for key in imageDict:
        if key not in ignoreKeys:
            invertedDict[key]=imageDict[key]
    
    # If we have makeDS9Regions = True here, we overwrite the existing .reg files from when we ran on the non-inverted maps
    photometry.findObjects(invertedDict, threshold = thresholdSigma, minObjPix = minObjPix,
                           rejectBorder = rejectBorder, diagnosticsDir = diagnosticsDir,
                           invertMap = True, makeDS9Regions = False, findCenterOfMass = findCenterOfMass)    

    # For fixed filter scale
    # Adds fixed_SNR values to catalogs for all maps
    if 'photFilter' in photometryOptions.keys():
        photFilter=photometryOptions['photFilter']
    else:
        photFilter=None
    if photFilter != None:
        photometry.getSNValues(invertedDict, SNMap = 'file', prefix = 'fixed_', template = photFilter, invertMap = True)   
        SNRKeys=['SNR', 'fixed_SNR']
    else:
        SNRKeys=['SNR']
        
    catalogTools.mergeCatalogs(invertedDict)
    catalogTools.makeOptimalCatalog(invertedDict, minSNToIncludeInOptimalCatalog)
    
    contaminTabDict=estimateContamination(invertedDict, imageDict, SNRKeys, 'invertedMap', diagnosticsDir)

    for k in contaminTabDict.keys():
        fitsOutFileName=diagnosticsDir+os.path.sep+"%s_contaminationEstimate.fits" % (k)
        if os.path.exists(fitsOutFileName) == True:
            os.remove(fitsOutFileName)
        contaminTab=contaminTabDict[k]
        contaminTab.write(fitsOutFileName)
        
    return contaminTabDict

#------------------------------------------------------------------------------------------------------------
def plotContamination(contaminTabDict, diagnosticsDir):
    """Makes contamination rate plots, output stored under diagnosticsDir
    
    """

    plotSettings.update_rcParams()

    for k in contaminTabDict.keys():
        if k.find('fixed') != -1:
            SNRKey="fixed_SNR"
            SNRLabel="SNR$_{\\rm 2.4}$"
        else:
            SNRKey="SNR"
            SNRLabel="SNR"
        binEdges=contaminTabDict[k][SNRKey]
        cumContamination=contaminTabDict[k]['cumContamination']
        plt.figure(figsize=(9,6.5))
        ax=plt.axes([0.10, 0.11, 0.87, 0.87])  
        plt.plot(binEdges, cumContamination, 'k-')# % (l))#, label = legl)
        plt.xlabel("%s" % (SNRLabel))#, fontdict = fontDict)
        plt.ylabel("Contamination fraction > %s" % (SNRLabel))#, fontdict = fontDict)
        allLabels=['4.0', '', '', '', '', '5.0', '', '', '', '', '6.0', '', '', '', '', '7.0', '', '', '', '', '8.0']
        allTicks=np.arange(4.0, 8.2, 0.2)
        plt.xticks(allTicks, allLabels)
        plt.xlim(4, 8)
        #plt.xlim(binMin, 10.01)#binMax)
        plt.ylim(-0.05, 0.6)
        #plt.legend()
        plt.savefig(diagnosticsDir+os.path.sep+"%s_contaminationEstimate.pdf" % (k))
        plt.close()  
        
#------------------------------------------------------------------------------------------------------------
def estimateContamination(contamSimDict, imageDict, SNRKeys, label, diagnosticsDir):
    """Performs the actual contamination estimate, makes output under diagnosticsDir.
        
    Use label to set a prefix for output (plots / .fits tables), e.g., label = "skySim"
    
    """
    
    invertedDict=contamSimDict
    contaminTabDict={}
    for SNRKey in SNRKeys:
        #catalogTools.catalog2DS9(invertedDict['optimalCatalog'], rootOutDir+os.path.sep+"skySimCatalog_%s_gtr_5.reg" % (SNRKey), 
                                 #constraintsList = ['%s > 5' % (SNRKey)])
        
        invertedSNRs=[]
        for obj in invertedDict['optimalCatalog']:
            invertedSNRs.append(obj[SNRKey])
        invertedSNRs=np.array(invertedSNRs)
        invertedSNRs.sort()
        numInverted=np.arange(len(invertedSNRs))+1
        
        candidateSNRs=[]
        for obj in imageDict['optimalCatalog']:
            candidateSNRs.append(obj[SNRKey])
        candidateSNRs=np.array(candidateSNRs)
        candidateSNRs.sort()
        numCandidates=np.arange(len(candidateSNRs))+1
        
        binMin=4.0
        binMax=20.0
        binStep=0.2
        binEdges=np.linspace(binMin, binMax, (binMax-binMin)/binStep+1)
        binCentres=(binEdges+binStep/2.0)[:-1]
        candidateSNRHist=np.histogram(candidateSNRs, bins = binEdges)
        invertedSNRHist=np.histogram(invertedSNRs, bins = binEdges)    
        
        cumSumCandidates=[]
        cumSumInverted=[]
        for i in range(binCentres.shape[0]):
            cumSumCandidates.append(candidateSNRHist[0][i:].sum())
            cumSumInverted.append(invertedSNRHist[0][i:].sum())
        cumSumCandidates=np.array(cumSumCandidates, dtype = float)
        cumSumInverted=np.array(cumSumInverted, dtype = float)
        
        # Plot cumulative contamination estimate (this makes more sense than plotting purity, since we don't know
        # that from what we're doing here, strictly speaking)
        cumContamination=np.zeros(cumSumCandidates.shape)
        mask=np.greater(cumSumCandidates, 0)
        cumContamination[mask]=cumSumInverted[mask]/cumSumCandidates[mask]
        
        # Remember, this is all cumulative (> SNR, so lower bin edges)
        contaminDict={}
        contaminDict['%s' % (SNRKey)]=binEdges[:-1]
        contaminDict['cumSumCandidates']=cumSumCandidates
        contaminDict['cumSumInverted']=cumSumInverted
        contaminDict['cumContamination']=cumContamination       
        
        # Convert to .fits table
        contaminTab=atpy.Table()
        for key in contaminDict.keys():
            contaminTab.add_column(atpy.Column(contaminDict[key], key))
        
        contaminTabDict['%s_%s' % (label, SNRKey)]=contaminTab
        
    return contaminTabDict

#------------------------------------------------------------------------------------------------------------
def calcR500Mpc(z, M500):
    """Given z, M500 (in MSun), returns R500 in Mpc, with respect to critical density.
    
    """
    
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
def makeArnaudModelProfile(z, M500, obsFreqGHz):
    """Given z, M500 (in MSun), returns dictionary containing Arnaud model profile (well, knots from spline 
    fit, 'tckP' - assumes you want to interpolate onto an array with units of degrees) and parameters 
    (particularly 'y0', 'theta500Arcmin').
    
    Used by ArnaudModelFilter
    
    """

    bRange=np.linspace(0, 30, 1000)
    cylPProfile=[]
    for b in bRange:
        cylPProfile.append(gnfw.integrated(b))
    cylPProfile=np.array(cylPProfile)
    
    # Normalise to 1 at centre
    cylPProfile=cylPProfile/cylPProfile.max()

    # Calculate R500Mpc, theta500Arcmin corresponding to given mass and redshift
    theta500Arcmin=calcTheta500Arcmin(z, M500)
    
    # Map between b and angular coordinates for random model
    # Note we fix c500 here to Arnaud value, we could leave it free
    c500=1.177
    thetaDegRange=(bRange*(theta500Arcmin/60.0))/c500
    tckP=interpolate.splrep(thetaDegRange, cylPProfile)
    
    # Get Y500 from M500 according to Arnaud et al. (eq. 25, cylindrical relation)
    # NOTE: Although sr is better units for comparing to Arnaud plots etc., arcmin2 is easier for rescaling below
    arnaudY500_arcmin2=calcY500FromM500_Arnaud(M500, z, units = 'arcmin2')

    # The above is the number we want, now normalise the profile to get that inside R500 and what deltaT0 is
    fidDeltaT0=-500.0
    yProfile=mapTools.convertToY(fidDeltaT0*cylPProfile, obsFrequencyGHz = obsFreqGHz)
    tcky=interpolate.splrep(thetaDegRange, yProfile)
    fineDegRange=np.linspace(0, theta500Arcmin/60.0, 1000)
    fineyProfile=interpolate.splev(fineDegRange, tcky)    
    YArcmin2=np.trapz(fineyProfile*np.pi*2*fineDegRange*60, fineDegRange*60)
    norm=arnaudY500_arcmin2/YArcmin2
    deltaT0=fidDeltaT0*norm
    y0=mapTools.convertToY(deltaT0, obsFrequencyGHz = obsFreqGHz)

    return {'tckP': tckP, 'y0': y0, 'deltaT0': deltaT0, 'theta500Arcmin': theta500Arcmin, 
            'Y500Arcmin2': arnaudY500_arcmin2, 'rDeg': thetaDegRange}

#------------------------------------------------------------------------------------------------------------
def calcY500FromM500_Arnaud(M500, z, units = 'sr'):
    """Calculate Y500 (in arcmin2) given M500, z following eq. 25 in Arnaud et al.
    
    Units can be 'sr' or 'arcmin2' 
    
    NOTE: this is the Y_cyl relation, so offset from Y_sph data (if e.g., looking at Fig. 10 of Arnaud) 
        
    """
    
    Ez=astCalc.Ez(z)    # h(z) in Arnaud speak
    Hz=astCalc.Ez(z)*astCalc.H0  
    
    # Now we need to adopt a scaling between Y500 (cylindrical for us) and mass
    # Let's go with Arnaud et al. Section 6.3 (eq. 25)
    alpha=1.78
    logBx=-4.665    # for 1 x R500 - read Arnaud more carefully!
    constantsTimesYSZ=np.power(10, logBx)*np.power(M500/3e14, alpha)
    YSZ_R500_sr=constantsTimesYSZ/(np.power(Ez, -2.0/3.0)*np.power(astCalc.da(z), 2))   # in steradians
    srToArcmin2=np.power(np.radians(1.0/60.0), 2)
    YSZ_R500_arcmin2=YSZ_R500_sr/srToArcmin2
        
    if units == 'sr':
        return YSZ_R500_sr
    elif units == 'arcmin2':
        return YSZ_R500_arcmin2
    else:
        raise Exception, "didn't understand units"

#------------------------------------------------------------------------------------------------------------
def makeBeamModelSignalMap(obsFreqGHz, degreesMap, wcs, beamFileName):
    """Makes a 2d signal only map containing the given beam.
    
    Returns signalMap, inputSignalProperties
    
    """
    
    # Load Matthew's beam profile
    beamData=np.loadtxt(beamFileName).transpose()
    profile1d=beamData[1]
    rArcmin=beamData[0]*60.0
                
    # Turn 1d profile into 2d
    rRadians=np.radians(rArcmin/60.0)
    r2p=interpolate.interp1d(rRadians, profile1d, bounds_error=False, fill_value=0.0)
    profile2d=r2p(np.radians(degreesMap))
    signalMap=liteMap.liteMapFromDataAndWCS(profile2d, wcs)
            
    # The ratio by which beam smoothing biases the intrinsic deltaT0
    beamDecrementBias=1.0#deltaT0/profile1d[0]  # assuming rDeg[0] is at 0 # 

    # Pixel window function
    inputSignalProperties={}
    fineWCS=wcs.copy()
    fineWCS.header['CDELT1']=fineWCS.header['CDELT1']*0.01
    fineWCS.header['CDELT2']=fineWCS.header['CDELT2']*0.01
    fineWCS.updateFromHeader()
    cRA, cDec=fineWCS.getCentreWCSCoords()
    degXMap, degYMap=nemoCython.makeXYDegreesDistanceMaps(np.zeros([fineWCS.header['NAXIS2'], fineWCS.header['NAXIS1']]), 
                                                            fineWCS, cRA, cDec, 1.0)
    degRMap=nemoCython.makeDegreesDistanceMap(np.zeros([fineWCS.header['NAXIS2'], fineWCS.header['NAXIS1']]), 
                                                fineWCS, cRA, cDec, 1.0)
    fineScaleProfile2d=r2p(np.radians(degRMap))
    mask=np.logical_and(np.less(abs(degXMap), wcs.getXPixelSizeDeg()/2.), 
                        np.less(abs(degYMap), wcs.getYPixelSizeDeg()/2.))
    inputSignalProperties['pixWindowFactor']=fineScaleProfile2d[mask].mean()
    
    return signalMap, inputSignalProperties
    
#------------------------------------------------------------------------------------------------------------
def makeArnaudModelSignalMap(z, M500, obsFreqGHz, degreesMap, wcs, beamFileName):
    """Makes a 2d signal only map containing an Arnaud model cluster. Units of M500 are MSun.
    
    degreesMap is a 2d array containing radial distance from the centre - the output map will have the same
    dimensions and pixel scale (see nemoCython.makeDegreesDistanceMap).
    
    Returns the map (2d array) and a dictionary containing the properties of the inserted cluster model.
    
    """
        
    # Broken out the Arnaud model code from here into simsTools
    signalDict=makeArnaudModelProfile(z, M500, obsFreqGHz)
    tckP=signalDict['tckP']
    y0=signalDict['y0']
    theta500Arcmin=signalDict['theta500Arcmin']
    deltaT0=signalDict['deltaT0']
    arnaudY500_arcmin2=signalDict['Y500Arcmin2']
    
    # Setup 1d profile
    rDeg=np.linspace(0.0, 1.0, 5000)
    profile1d=deltaT0*interpolate.splev(rDeg, tckP)
            
    # Apply beam to profile
    # NOTE: Do not disable this
    # Load Matthew's beam profile and interpolate onto signal profile coords
    beamData=np.loadtxt(beamFileName).transpose()
    profile1d_beam=beamData[1]
    rDeg_beam=beamData[0]
    tck_beam=interpolate.splrep(rDeg_beam, profile1d_beam)
    profile1d_beam=interpolate.splev(rDeg, tck_beam)
                        
    # Turn 1d profiles into 2d and convolve signal with beam
    # Convolving just redistributes the total signal to different spatial scales
    # So sum should be the same after the convolution - this is how we normalise below
    # (this ignores some power shifted beyond the boundary of the map if convolution kernel doesn't have
    # compact support)
    rRadians=np.radians(rDeg)
    radiansMap=np.radians(degreesMap)
    r2p=interpolate.interp1d(rRadians, profile1d, bounds_error=False, fill_value=0.0)
    profile2d=r2p(radiansMap)
    r2p_beam=interpolate.interp1d(rRadians, profile1d_beam, bounds_error=False, fill_value=0.0)
    profile2d_beam=r2p_beam(radiansMap)
    smoothedProfile2d=fft.fftshift(fft.ifft2(fft.fft2(profile2d)*fft.fft2(profile2d_beam))).real
    normFactor=profile2d.sum()/smoothedProfile2d.sum()
    smoothedProfile2d=smoothedProfile2d*normFactor
    signalMap=liteMap.liteMapFromDataAndWCS(smoothedProfile2d, wcs)        
    
    # Check profile2d integrates to give Arnaud value
    # Need solid angle map for this
    # NOTE: this does indeed recover the input Y500 IF we turn the beam smoothing off
    # With beam smoothing, get less than Arnaud value because the smoothing shifts some signal beyond R500 cut off
    # This makes sense
    # NOTE: to do area-type map scaling, we do need the area mask also
    pixAreaMapArcmin2=mapTools.getPixelAreaArcmin2Map(signalMap.data, wcs)
    R500Radians=np.radians(theta500Arcmin/60.0)
    mask=np.less(radiansMap, R500Radians)
    #YRec=mapTools.convertToY(np.sum(profile2d[mask]*pixAreaMapArcmin2[mask]), obsFrequencyGHz = mapObsFreqGHz)
    YRec=mapTools.convertToY(np.sum(smoothedProfile2d[mask]*pixAreaMapArcmin2[mask]), obsFrequencyGHz = obsFreqGHz)
    
    # Correction factor for signal smeared beyond R500 if cluster really does follow Arnaud profile
    # We would multiply integrated Ys by this to correct for this bias - if we were actually able to measure
    # on map within some radius
    # We could do this for 2' radius aperture or something also if we wanted to
    Y500BeamCorrection=arnaudY500_arcmin2/YRec
    
    # Beam decrement bias: in this case, now in 2d we don't have resolution below pixel size
    # But that's easy to fix: finer resolution grid interpolating, don't do full map size
    # 1D convolution for beam decrement bias
    # NOTE: changed how this is defined, i.e., now multiply by this to get true y0
    symBeam=np.zeros(profile1d_beam.shape[0]*2)
    symBeam[:profile1d_beam.shape[0]]=profile1d_beam[::-1]
    symBeam[profile1d_beam.shape[0]:]=profile1d_beam
    symProfile=np.zeros(profile1d.shape[0]*2)
    symProfile[:profile1d.shape[0]]=profile1d[::-1]
    symProfile[profile1d.shape[0]:]=profile1d
    smoothedProfile1d=fft.fftshift(fft.ifft(fft.fft(symProfile)*fft.fft(symBeam))).real
    normFactor=symProfile.sum()/smoothedProfile1d.sum()
    smoothedProfile1d=normFactor*smoothedProfile1d
    beamDecrementBias=abs(profile1d).max()/abs(smoothedProfile1d).max()               

    # For sanity checking later, let's dump the input properties of the model into a dictionary
    # Then we can apply whatever filter we like later and check that we're recovering these
    # These are BEFORE beam smoothing in this case, i.e., just from the input Arnaud model
    inputSignalProperties={'deltaT0': deltaT0, 'y0': y0, 'theta500Arcmin': theta500Arcmin, 
                            'Y500Arcmin2': arnaudY500_arcmin2, 'obsFreqGHz': obsFreqGHz}

    return signalMap, inputSignalProperties
    
#------------------------------------------------------------------------------------------------------------
def fitQ(parDict, diagnosticsDir, filteredMapsDir):
    """Calculates Q on a grid, and then fits (theta, Q) with a spline, saving a plot and the (theta, Q) array
    as a table in the diagnostics dir.
    
    This can be generalised, but for now is hard coded to use the Arnaud model.
    
    NOTE: This is also currently signal frequency only - we're assuming that beamFileName is given under
    parDict['unfilteredMaps'].
    
    """
    
    outFileName=diagnosticsDir+os.path.sep+"QFit.fits"
    
    if os.path.exists(outFileName) == False:
        print ">>> Fitting for Q ..."
        
        # Spin through the filter kernels
        photFilterLabel=parDict['photometryOptions']['photFilter']
        filterList=parDict['mapFilters']
        for f in filterList:
            if f['label'] == photFilterLabel:
                ref=f
    
        # Need the beam - this assumes we are single frequency only
        beamFileName=parDict['unfilteredMaps'][0]['beamFileName']
        
        # Some faffing to get map pixel scale
        img=pyfits.open(filteredMapsDir+os.path.sep+photFilterLabel+"_SNMap.fits")
        wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
        RADeg, decDeg=wcs.getCentreWCSCoords()
        clipDict=astImages.clipImageSectionWCS(img[0].data, wcs, RADeg, decDeg, 1.0)
        wcs=clipDict['wcs']

        # Ref kernel
        kernImg=pyfits.open(diagnosticsDir+os.path.sep+"kern2d_%s.fits" % (photFilterLabel))
        kern2d=kernImg[0].data

        # Blank map - to make Q calc more accurate, use finer pixel scale
        signalMap=np.zeros(clipDict['data'].shape)
        degreesMap=nemoCython.makeDegreesDistanceMap(signalMap, wcs, RADeg, decDeg, 1.0)

        # M, z ranges for Q calc
        # NOTE: ref filter that sets scale we compare to must ALWAYS come first
        MRange=[ref['params']['M500MSun']]
        zRange=[ref['params']['z']]
        MRange=MRange+np.logspace(13.5, 16.0, 10).tolist()
        zRange=zRange+np.arange(0.1, 1.7, 0.2).tolist()

        # Make signal only maps and filter them with the ref kernel
        t0=time.time()
        Q=[]
        theta500Arcmin=[]
        for z in zRange:
            for M500MSun in MRange:
                # We don't actually care about freq here because we deal with pressure profile itself
                # And Q measurement is relative 
                modelDict=makeArnaudModelProfile(z, M500MSun, 148.0)   
                rDeg=modelDict['rDeg']
                profile1d=interpolate.splev(rDeg, modelDict['tckP'])
                r2p=interpolate.interp1d(rDeg, profile1d, bounds_error=False, fill_value=0.0)
                signalMap=r2p(degreesMap)
                # Apply high-pass then convolve (high-pass effect seems negligible)
                signalMap=mapTools.subtractBackground(signalMap, wcs, smoothScaleDeg = kernImg[0].header['BCKSCALE']/60.)
                filteredSignal=ndimage.convolve(signalMap, kern2d) 
                Q.append(filteredSignal.max())
                theta500Arcmin.append(modelDict['theta500Arcmin'])
        Q=np.array(Q)
        Q=Q/Q[0]
        theta500Arcmin=np.array(theta500Arcmin)
        t1=time.time()
        
        # Sort and do spline fit... save .fits table of theta, Q
        QTab=atpy.Table()
        QTab.add_column(atpy.Column(Q, 'Q'))
        QTab.add_column(atpy.Column(theta500Arcmin, 'theta500Arcmin'))
        QTab.sort('theta500Arcmin')
        if os.path.exists(outFileName) == True:
            os.remove(outFileName)
        QTab.write(outFileName)
        
        # Fit with spline
        #coeffs=np.polyfit(theta500Arcmin, Q, 16)
        tck=interpolate.splrep(QTab['theta500Arcmin'], QTab['Q'])
        
        # Plot
        plotSettings.update_rcParams()
        #fontSize=18.0
        #fontDict={'size': fontSize, 'family': 'serif'}
        plt.figure(figsize=(9,6.5))
        ax=plt.axes([0.10, 0.11, 0.88, 0.88])
        #plt.tick_params(axis='both', which='major', labelsize=15)
        #plt.tick_params(axis='both', which='minor', labelsize=15)       
        thetaArr=np.linspace(0, 30, 300)
        #plt.plot(thetaArr, np.poly1d(coeffs)(thetaArr), 'k-')
        plt.plot(thetaArr, interpolate.splev(thetaArr, tck), 'k-')
        plt.plot(theta500Arcmin, Q, 'D', ms = 8)
        #plt.plot(thetaArr, simsTools.calcQ_H13(thetaArr), 'b--')
        #plt.xlim(0, 9)
        plt.ylim(0, Q.max()*1.05)
        plt.xlim(0, thetaArr.max())
        plt.xlabel("$\\theta_{\\rm 500c}$ (arcmin)")
        plt.ylabel("$Q$ ($M_{\\rm 500c}$, $z$)")
        plt.savefig(diagnosticsDir+os.path.sep+"QFit.pdf")
        plt.close()
    
    else:
        
        print ">>> Loading previously cached Q fit ..."
        #coeffs=np.load(outFileName)
        QTab=atpy.Table().read(outFileName)
        tck=interpolate.splrep(QTab['theta500Arcmin'], QTab['Q'])
    
    return tck

#------------------------------------------------------------------------------------------------------------
def calcQ(theta500Arcmin, tck):
    """Returns Q, given theta500Arcmin, and a set of spline fit knots for (theta, Q).
    
    """
    
    #Q=np.poly1d(coeffs)(theta500Arcmin)
    Q=interpolate.splev(theta500Arcmin, tck)
    
    return Q
    
#------------------------------------------------------------------------------------------------------------
def getQCoeffsH13():
    """Returns result of a polynomial fit to (theta, Q) given in H13
    
    NOTE: this fit is not accurate for theta500Arcmin > 5.9'
    
    """
    
    coeffs=np.array([-1.12816214e-04, 2.37255951e-03, -1.62915564e-02, 9.87118185e-03, 3.50817483e-01, 
                     -1.38970056e-01])
    
    return coeffs

#------------------------------------------------------------------------------------------------------------
def calcFRel(z, M500):
    """Calculates relativistic correction to SZ effect, given z, M500 in MSun, as per H13. 
    
    This assumes the Arnaud et al. (2005) M-T relation, and applies formulae of Itoh et al. (1998)
    
    See H13 Section 2.2.
    
    """
    
    #kB=1.38e-23
    #me=9.11e-31
    #e=1.6e-19
    #c=3e8
    
    ## Get T in keV from Arnaud et al. 2005
    #A=3.84e14
    #B=1.71
    #TkeV=5.*np.power(((astCalc.Ez(z)*M500)/A), 1/B)
    #T=(((TkeV*1000)*e)/kB)
    
    #t=(kB*T)/(me*c**2)

    # H13
    m=M500/3e14
    t=-0.00848*np.power(m*astCalc.Ez(z), -0.585)
    fRel=1+3.79*t-28.2*t*t
    
    return fRel

#------------------------------------------------------------------------------------------------------------
def getM500FromP(P, log10M, calcErrors = True):
    """Returns M500 as the maximum likelihood value from given P(log10M) distribution, together with 
    1-sigma error bars (M500, -M500Err, +M500 err).
    
    """

    # Find max likelihood and integrate to get error bars
    tckP=interpolate.splrep(log10M, P)
    fineLog10M=np.linspace(log10M.min(), log10M.max(), 1e5)
    fineP=interpolate.splev(fineLog10M, tckP)
    fineP=fineP/np.trapz(fineP, fineLog10M)
    try:
        index=np.where(fineP == fineP.max())[0][0]
    except:
        print "argh"
        IPython.embed()
        sys.exit()
    
    clusterLogM500=fineLog10M[index]
    clusterM500=np.power(10, clusterLogM500)/1e14

    if calcErrors == True:
        for n in range(fineP.shape[0]):
            minIndex=index-n
            maxIndex=index+n
            if minIndex < 0 or maxIndex > fineP.shape[0]:
                # This shouldn't happen; if it does, probably y0 is in the wrong units
                print "WARNING: outside M500 range"
                clusterLogM500=None
                break            
            p=np.trapz(fineP[minIndex:maxIndex], fineLog10M[minIndex:maxIndex])
            if p >= 0.6827:
                clusterLogM500Min=fineLog10M[minIndex]
                clusterLogM500Max=fineLog10M[maxIndex]
                break        
        clusterM500MinusErr=(np.power(10, clusterLogM500)-np.power(10, clusterLogM500Min))/1e14
        clusterM500PlusErr=(np.power(10, clusterLogM500Max)-np.power(10, clusterLogM500))/1e14
    else:
        clusterM500MinusErr=0.
        clusterM500PlusErr=0.
    
    return clusterM500, clusterM500MinusErr, clusterM500PlusErr

#------------------------------------------------------------------------------------------------------------
def y0FromLogM500(log10M500, z, tckQFit, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, sigma_int = 0.2):
    """Predict y0~ given logM500 (in MSun) and redshift. Default scaling relation parameters are A10 (as in
    H13).
    
    Returns y0~, theta500Arcmin, Q
    
    """
    
    M500=np.power(10, log10M500)
    theta500Arcmin=calcTheta500Arcmin(z, M500)
    Q=calcQ(theta500Arcmin, tckQFit)
    # UPP relation according to H13
    # NOTE: m in H13 is M/Mpivot
    # NOTE: this goes negative for crazy masses where the Q polynomial fit goes -ve, so ignore those
    y0pred=tenToA0*np.power(astCalc.Ez(z), 2)*np.power(M500/Mpivot, 1+B0)*Q*calcFRel(z, M500)
    
    return y0pred, theta500Arcmin, Q
            
#------------------------------------------------------------------------------------------------------------
def calcM500Fromy0(y0, y0Err, z, zErr, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, sigma_int = 0.2, 
                   tckQFit = None, mockSurvey = None, applyMFDebiasCorrection = True, calcErrors = True):
    """Returns M500 +/- errors in units of 10^14 MSun, calculated assuming a y0 - M relation (default values
    assume UPP scaling relation from Arnaud et al. 2010), taking into account the steepness of the mass
    function. The approach followed is described in H13, Section 3.2.
    
    Here, mockSurvey is a MockSurvey object. We're using this to handle the halo mass function calculations
    (in turn using the hmf module).
    
    tckQFit is a set of spline knots, as returned by fitQ.
    
    If applyMFDebiasCorrection == True, apply correction that accounts for steepness of mass function.
    
    If calcErrors == False, error bars are not calculated, they are just set to zero.
    
    """
    
    if y0 < 0:
        raise Exception, 'y0 cannot be negative'
    
    if mockSurvey == None and applyMFDebiasCorrection == True:
        raise Exception, 'MockSurvey object must be supplied for the mass function shape de-bias correction to work'
    
    try:
        log10M=mockSurvey.log10M
    except:
        log10M=np.linspace(13., 16., 300)
        
    # For marginalising over photo-z errors
    if zErr > 0:
        zRange=np.linspace(0, 2.0, 401)
        Pz=np.exp(-np.power(z-zRange, 2)/(2*(np.power(zErr, 2))))
        Pz=Pz/np.trapz(Pz, zRange)
    else:
        zRange=[z]
        Pz=np.ones(len(zRange))

    # M500
    Py0GivenM=[]
    QArr=[]
    theta500ArcminArr=[]
    for log10M500 in log10M:
        lnPy0=np.zeros(len(zRange))
        for i in range(len(zRange)):
            zi=zRange[i]
            # UPP relation according to H13
            # NOTE: m in H13 is M/Mpivot
            y0pred, theta500Arcmin, Q=y0FromLogM500(log10M500, zi, tckQFit, tenToA0 = tenToA0, B0 = B0, Mpivot = Mpivot, sigma_int = sigma_int)
            theta500ArcminArr.append(theta500Arcmin)
            QArr.append(Q)
            if y0pred > 0:
                log_y0=np.log(y0)
                log_y0Err=np.log(y0+y0Err)-log_y0
                log_y0pred=np.log(y0pred)
                lnprob=-np.power(log_y0-log_y0pred, 2)/(2*(np.power(log_y0Err, 2)+np.power(sigma_int, 2)))
            else:
                lnprob=-np.inf
            lnPy0[i]=lnprob
        Py0GivenM.append(np.sum(np.exp(lnPy0)*Pz))
    Py0GivenM=np.array(Py0GivenM)
    if np.any(np.isnan(Py0GivenM)) == True:
        print "nan"
        IPython.embed()
        sys.exit()
    
    # Normalise
    Py0GivenM=Py0GivenM/np.trapz(Py0GivenM, log10M)
    if applyMFDebiasCorrection == True and mockSurvey != None:
        PLog10M=mockSurvey.getPLog10M(z)
        PLog10M=PLog10M/np.trapz(PLog10M, log10M)
        M500, errM500Minus, errM500Plus=getM500FromP(Py0GivenM*PLog10M, log10M, calcErrors = calcErrors)
    else:
        M500, errM500Minus, errM500Plus=0.0, 0.0, 0.0

    # M500 without de-biasing for mass function shape (this gives the ~15% offset compared to Planck)
    M500Uncorr, errM500UncorrMinus, errM500UncorrPlus=getM500FromP(Py0GivenM, log10M, calcErrors = calcErrors)
    
    if M500Uncorr == 0:
        print "M500 fail"
        IPython.embed()
        sys.exit()

    #print "Add Q"
    #IPython.embed()
    #sys.exit()
    
    return {'M500': M500, 'M500_errPlus': errM500Plus, 'M500_errMinus': errM500Minus,
            'M500Uncorr': M500Uncorr, 'M500Uncorr_errPlus': errM500UncorrPlus, 
            'M500Uncorr_errMinus': errM500UncorrMinus}

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
            raise Exception, "M500c -> M200m conversion didn't converge quickly enough"
        
    M200m=scaleFactor*M500c
    
    return M200m
