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
import catalogTools
import photometry
import gnfw
import numpy as np
import os
import math
import pylab as plt
import pickle
import sys
import operator
import pyximport; pyximport.install()
import nemoCython
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
def estimateContaminationFromInvertedMaps(imageDict, thresholdSigma, minObjPix, rejectBorder, 
                                          minSNToIncludeInOptimalCatalog, photometryOptions, diagnosticsDir = None):
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
                           invertMap = True, makeDS9Regions = False)    

    # For fixed filter scale
    # Adds fixed_SNR values to catalogs for all maps
    if 'photFilter' in photometryOptions.keys():
        photFilter=photometryOptions['photFilter']
    else:
        photFilter=None
    if photFilter != None:
        photometry.getSNValues(invertedDict, SNMap = 'file', prefix = 'fixed_', template = photFilter, invertMap = True)   
    
    catalogTools.mergeCatalogs(invertedDict)
    catalogTools.makeOptimalCatalog(invertedDict, minSNToIncludeInOptimalCatalog)
    
    # Do everything else with both SNR and fixed_SNR
    contaminDictList=[]
    for SNRKey in ['SNR', 'fixed_SNR']:
        catalogTools.catalog2DS9(invertedDict['optimalCatalog'], diagnosticsDir+os.path.sep+"invertedMapsCatalog_%s_gtr_5.reg" % (SNRKey), 
                                constraintsList = ['%s > 5' % (SNRKey)])
        
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
        binStep=0.1
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

        xtickLabels=[]
        xtickValues=[]
        fmodMajorTicks=np.fmod(binEdges, 5)
        fmodMinorTicks=np.fmod(binEdges, 1)
        for i in range(len(binEdges)):
            if fmodMinorTicks[i] == 0:
                xtickValues.append(binEdges[i])
                if fmodMajorTicks[i] == 0:
                    xtickLabels.append('%d' % (binEdges[i]))
                else:
                    xtickLabels.append('')
                
        # Plot cumulative detections > SNR for both inverted map catalog and actual catalog
        plt.plot(binEdges[:-1], cumSumInverted, 'r-', label = 'inverted maps')
        plt.plot(binEdges[:-1], cumSumCandidates, 'b-', label = 'candidates')
        plt.xlabel("%s" % (SNRKey))
        plt.ylabel("Number > %s" % (SNRKey))
        plt.semilogx()
        plt.xticks(xtickValues, xtickLabels)
        plt.xlim(binMin, binMax)
        plt.legend()
        plt.savefig(diagnosticsDir+os.path.sep+"cumulative_%s.png" % (SNRKey))
        plt.close()
        
        # Plot cumulative contamination estimate (this makes more sense than plotting purity, since we don't know
        # that from what we're doing here, strictly speaking)
        cumContamination=cumSumInverted/cumSumCandidates
        cumContamination[np.isnan(cumContamination)]=0.0
        plt.plot(binEdges[:-1], cumContamination, 'k-')
        plt.xlabel("%s" % (SNRKey))
        plt.ylabel("Estimated contamination > %s" % (SNRKey))
        plt.semilogx()
        plt.xticks(xtickValues, xtickLabels)
        plt.xlim(binMin, binMax)
        plt.ylim(0, 1)
        plt.savefig(diagnosticsDir+os.path.sep+"contaminationEstimate_%s.png" % (SNRKey))
        plt.close()    
        
        # Remember, this is all cumulative (> SNR, so lower bin edges)
        contaminDict={}
        contaminDict['%s' % (SNRKey)]=binEdges[:-1]
        contaminDict['cumSumCandidates']=cumSumCandidates
        contaminDict['cumSumInverted']=cumSumInverted
        contaminDict['cumContamination']=cumContamination       
        
        # Wite a .fits table
        contaminTab=atpy.Table()
        for key in contaminDict.keys():
            contaminTab.add_column(atpy.Column(contaminDict[key], key))
        fitsOutFileName=diagnosticsDir+os.path.sep+"contaminationEstimate_%s.fits" % (SNRKey)
        if os.path.exists(fitsOutFileName) == True:
            os.remove(fitsOutFileName)
        contaminTab.write(fitsOutFileName)
        
        contaminDictList.append(contaminDict)
        
    return contaminDictList

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
def fitQ(parDict, diagnosticsDir, filteredMapsDir):
    """Calculates Q on a grid, and then fits (theta, Q) with a polynomial, saving a plot and the coeffs
    array in the diagnostics dir.
    
    This can be generalised, but for now is hard coded to use the Arnaud model.
    
    """
    
    outFileName=diagnosticsDir+os.path.sep+"QCoeffs.npy"
    
    if os.path.exists(outFileName) == False:
        print ">>> Fitting for Q ..."
        
        # Spin through the filter kernels
        photFilterLabel=parDict['photometryOptions']['photFilter']
        filterList=parDict['mapFilters']
        for f in filterList:
            if f['label'] == photFilterLabel:
                ref=f
                
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
        MRange=MRange+np.logspace(13.7, 15.3, 5).tolist()
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
                #signalNorm=modelDict['y0']/interpolate.splev(0., modelDict['tckP'])
                #signalMap=signalNorm*signalMap
                # Apply high-pass then convolve (high-pass effect seems negligible)
                signalMap=mapTools.subtractBackground(signalMap, wcs, smoothScaleDeg = kernImg[0].header['BCKSCALE']/60.)
                filteredSignal=ndimage.convolve(signalMap, kern2d) 
                Q.append(filteredSignal.max())
                theta500Arcmin.append(modelDict['theta500Arcmin'])
        Q=np.array(Q)
        Q=Q/Q[0]
        theta500Arcmin=np.array(theta500Arcmin)
        t1=time.time()

        # Fit with polynomial
        coeffs=np.polyfit(theta500Arcmin, Q, 12)

        # Plot
        plt.figure(figsize=(12, 8))
        plt.plot(theta500Arcmin, Q, 'r.')
        thetaArr=np.linspace(0, 10, 100)
        plt.plot(thetaArr, np.poly1d(coeffs)(thetaArr), 'k-')
        #plt.plot(thetaArr, simsTools.calcQ_H13(thetaArr), 'b--')
        plt.xlim(0, 10)
        plt.xlabel("$\\theta_{500}$ (arcmin)")
        plt.ylabel("Q")
        plt.savefig(diagnosticsDir+os.path.sep+"QFit.png")
        plt.close()

        # Save polynomial fit coeffs
        np.save(outFileName, coeffs)
    
    else:
        
        print ">>> Loading previously cached Q fit ..."
        coeffs=np.load(outFileName)
    
    return coeffs

#------------------------------------------------------------------------------------------------------------
def calcQ(theta500Arcmin, coeffs):
    """Returns Q, given theta500Arcmin, and a set of polynomial fit coefficients to (theta, Q).
    
    """
    
    Q=np.poly1d(coeffs)(theta500Arcmin)
    
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
def calcM500Fromy0(y0, y0Err, z, mockSurvey, tenToA0 = 4.95e-5, B0 = 0.08, Mpivot = 3e14, sigma_int = 0.2, 
                   QFitCoeffs = getQCoeffsH13()):
    """Returns M500 +/- errors in units of 10^14 MSun, calculated assuming a y0 - M relation (default values
    assume UPP scaling relation from Arnaud et al. 2010), taking into account the steepness of the mass
    function. The approach followed is described in H13, Section 3.2.
    
    Here, mockSurvey is a MockSurvey object. We're using this to handle the halo mass function calculations
    (in turn using the hmf module).
    
    QFitCoeffs are from a polynomial fit to (theta, Q), with theta given in arcmin.
    
    NOTE: not added marginalisation over z uncertainty yet
    
    """
    
    PLog10M, log10M=mockSurvey.getPLog10M(z)
    Py0GivenM=[]
    for log10M500 in log10M:

        M500=np.power(10, log10M500)
        theta500Arcmin=calcTheta500Arcmin(z, M500)
            
        # UPP relation according to H13
        # NOTE: m in H13 is M/Mpivot
        y0pred=tenToA0*np.power(astCalc.Ez(z), 2)*np.power(M500/Mpivot, 1+B0)*calcQ(theta500Arcmin, QFitCoeffs)*calcFRel(z, M500)
        log_y0=np.log(y0)
        log_y0Err=np.log(y0+y0Err)-log_y0
        log_y0pred=np.log(y0pred)
        lnprob=-np.power(log_y0-log_y0pred, 2)/(2*(np.power(log_y0Err, 2)+np.power(sigma_int, 2)))
        if np.isnan(lnprob) == False:
            Py0GivenM.append(np.exp(lnprob))
        else:
            Py0GivenM.append(0)

    Py0GivenM=np.array(Py0GivenM)

    # Find max likelihood and integrate to get error bars
    P=Py0GivenM*PLog10M
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
    for n in range(fineP.shape[0]):
        minIndex=index-n
        maxIndex=index+n
        if minIndex < 0 or maxIndex > fineP.shape[0]:
            # This shouldn't happen; if it does, probably y0 is in the wrong units
            print "outside M500 range"
            IPython.embed()
            sys.exit()
            clusterLogM500=None
            break            
        p=np.trapz(fineP[minIndex:maxIndex], fineLog10M[minIndex:maxIndex])
        if p >= 0.6827:
            clusterLogM500=fineLog10M[index]
            clusterLogM500Min=fineLog10M[minIndex]
            clusterLogM500Max=fineLog10M[maxIndex]
            break
        
    if np.any(clusterLogM500) != None:
        clusterM500=np.power(10, clusterLogM500)/1e14
        clusterM500MinusErr=(np.power(10, clusterLogM500)-np.power(10, clusterLogM500Min))/1e14
        clusterM500PlusErr=(np.power(10, clusterLogM500Max)-np.power(10, clusterLogM500))/1e14
    else:
        clusterM500=0.
        clusterM500MinusErr=0.
        clusterM500PlusErr=0.
    
    #print "return Q, QErr, theta500Arcmin, theta500ArcminErr"
    #IPython.embed()
    #sys.exit()
        
    return {'M500': clusterM500, 'M500_errPlus': clusterM500PlusErr, 'M500_errMinus': clusterM500MinusErr}

#------------------------------------------------------------------------------------------------------------
