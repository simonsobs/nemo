"""

Tools used by nemoSelFn

"""

import os
import sys
import resource
import glob
import numpy as np
import pylab as plt
import astropy.table as atpy
from astLib import *
from scipy import stats
from scipy import interpolate
from scipy.interpolate import InterpolatedUnivariateSpline as _spline
from scipy import ndimage
from scipy import optimize
from nemo import simsTools
from nemo import mapTools
from nemo import MockSurvey
from nemo import plotSettings
from collections import OrderedDict
import colorcet
import types
import pickle
import astropy.io.fits as pyfits
import time
import IPython
plt.matplotlib.interactive(False)

# If want to catch warnings as errors...
#import warnings
#warnings.filterwarnings('error')

#------------------------------------------------------------------------------------------------------------
def loadAreaMask(extName, diagnosticsDir):
    """Loads the survey area mask (i.e., after edge-trimming and point source masking, produced by nemo).
    
    Returns map array, wcs
    
    """
    
    areaImg=pyfits.open(diagnosticsDir+os.path.sep+"areaMask#%s.fits" % (extName))
    areaMap=areaImg[0].data
    wcs=astWCS.WCS(areaImg[0].header, mode = 'pyfits')
    areaImg.close()
    
    return areaMap, wcs

#------------------------------------------------------------------------------------------------------------
def getTileTotalAreaDeg2(extName, diagnosticsDir, masksList = [], footprintLabel = None):
    """Returns total area of the tile pointed at by extName (taking into account survey mask and point
    source masking).
    
    A list of other survey masks (e.g., from optical surveys like KiDS, HSC) can be given in masksList. 
    These should be file names for .fits images where 1 defines valid survey area, 0 otherwise. If given,
    this routine will return the area of intersection between the extra masks and the SZ survey.
            
    """
    
    areaMap, wcs=loadAreaMask(extName, diagnosticsDir)
    areaMapSqDeg=(mapTools.getPixelAreaArcmin2Map(areaMap, wcs)*areaMap)/(60**2)
    totalAreaDeg2=areaMapSqDeg.sum()
    
    if footprintLabel != None:  
        intersectMask=makeIntersectionMask(extName, diagnosticsDir, footprintLabel, masksList = masksList)
        totalAreaDeg2=(areaMapSqDeg*intersectMask).sum()        
        
    return totalAreaDeg2

#------------------------------------------------------------------------------------------------------------
def makeIntersectionMask(extName, diagnosticsDir, label, masksList = []):
    """Creates intersection mask between mask files given in masksList.
    
    Calculating the intersection is slow (~30 sec per tile per extra mask for KiDS), so intersection masks 
    are cached; label is used in output file names (e.g., diagnosticsDir/intersect_label#extName.fits)
    
    Can optionally be called without extraMasksList, IF the intersection mask has already been created and
    cached.
    
    NOTE: We assume masks have dec aligned with y direction for speed.
    
    Returns intersectionMask as array (1 = valid area, 0 = otherwise)
    
    """
        
    areaMap, wcs=loadAreaMask(extName, diagnosticsDir)
    RAMin, RAMax, decMin, decMax=wcs.getImageMinMaxWCSCoords()
    
    intersectFileName=diagnosticsDir+os.path.sep+"intersect_%s#%s.fits.gz" % (label, extName)
    if os.path.exists(intersectFileName) == True:
        intersectImg=pyfits.open(intersectFileName)
        intersectMask=intersectImg[0].data
        intersectImg.close()
    else:
        if masksList == []:
            raise Exception("didn't find previously cached intersection mask but makeIntersectionMask called with empty masksList")
        print("... creating %s intersection mask (%s) ..." % (label, extName)) 
        intersectMask=np.zeros(areaMap.shape)
        for fileName in masksList:
            maskImg=pyfits.open(fileName)
            for hdu in maskImg:
                if type(hdu) == pyfits.ImageHDU:
                    break
            maskWCS=astWCS.WCS(hdu.header, mode = 'pyfits')
            maskData=hdu.data
            # For speed: assuming mask is aligned with dec in y direction and no rotation
            blah, yMin=maskWCS.wcs2pix(0.0, decMin)
            blah, yMax=maskWCS.wcs2pix(0.0, decMax)
            yMin=int(round(yMin))
            yMax=int(round(yMax))
            if yMin < 0:
                yMin=0
            if yMax >= maskData.shape[0]:
                yMax=maskData.shape[0]-1
            ys, xs=np.where(maskData[yMin:yMax] == 1)
            if len(ys) > 0:
                RADec=maskWCS.pix2wcs(xs, ys+yMin)
                RADec=np.array(RADec)
                areaMapCoords=wcs.wcs2pix(RADec[:, 0], RADec[:, 1])
                areaMapCoords=np.array(areaMapCoords)
                areaMapCoords=np.array(np.round(areaMapCoords), dtype = int)
                xMask=np.logical_and(np.greater_equal(areaMapCoords[:, 0], 0), np.less(areaMapCoords[:, 0], areaMap.shape[1]-1))
                yMask=np.logical_and(np.greater_equal(areaMapCoords[:, 1], 0), np.less(areaMapCoords[:, 1], areaMap.shape[0]-1))
                for coord in areaMapCoords[np.logical_and(xMask, yMask)]:
                    intersectMask[coord[1], coord[0]]=1
        astImages.saveFITS(intersectFileName, intersectMask, wcs)
    
    return intersectMask

#------------------------------------------------------------------------------------------------------------
def getRMSTab(extName, photFilterLabel, diagnosticsDir, footprintLabel = None):
    """Makes a table containing fraction of map area in tile pointed to by extName against RMS values
    (so this compresses the information in the RMS maps). The first time this is run takes ~200 sec (for a
    1000 sq deg tile), but the result is cached.
    
    Can optionally take extra masks for specifying e.g. HSC footprint. Here, we assume these have already
    been made by makeIntersectionMask, and we can load them from the cache, identifying them through
    footprintLabel
        
    Returns RMSTab
    
    """
    
    # This can probably be sped up, but takes ~200 sec for a ~1000 sq deg tile, so we cache
    RMSTabFileName=diagnosticsDir+os.path.sep+"RMSTab_%s.fits" % (extName)
    if footprintLabel != None:
        RMSTabFileName=RMSTabFileName.replace(".fits", "_%s.fits" % (footprintLabel))
    if os.path.exists(RMSTabFileName) == False:
        print(("... making %s ..." % (RMSTabFileName)))
        RMSImg=pyfits.open(diagnosticsDir+os.path.sep+"RMSMap_Arnaud_M2e14_z0p4#%s.fits" % (extName))
        RMSMap=RMSImg[0].data

        areaMap, wcs=loadAreaMask(extName, diagnosticsDir)
        areaMapSqDeg=(mapTools.getPixelAreaArcmin2Map(areaMap, wcs)*areaMap)/(60**2)
        
        if footprintLabel != None:  
            intersectMask=makeIntersectionMask(extName, diagnosticsDir, footprintLabel)
            areaMapSqDeg=areaMapSqDeg*intersectMask        
            RMSMap=RMSMap*intersectMask

        RMSValues=np.unique(RMSMap[np.nonzero(RMSMap)])

        totalAreaDeg2=areaMapSqDeg.sum()
        
        fracArea=np.zeros(len(RMSValues))
        for i in range(len(RMSValues)):
            fracArea[i]=areaMapSqDeg[np.equal(RMSMap, RMSValues[i])].sum()
        RMSTab=atpy.Table()
        RMSTab.add_column(atpy.Column(fracArea, 'fracArea'))
        RMSTab.add_column(atpy.Column(RMSValues, 'y0RMS'))
        RMSTab.write(RMSTabFileName)
    else:
        #print("... reading %s ..." % (RMSTabFileName))
        RMSTab=atpy.Table().read(RMSTabFileName)
    
    return RMSTab

#------------------------------------------------------------------------------------------------------------
def calcTileWeightedAverageNoise(extName, photFilterLabel, diagnosticsDir, footprintLabel = None):
    """Returns weighted average noise value in the tile.
    
    """

    RMSTab=getRMSTab(extName, photFilterLabel, diagnosticsDir, footprintLabel = footprintLabel)
    RMSValues=np.array(RMSTab['y0RMS'])
    fracArea=np.array(RMSTab['fracArea'])
    tileRMSValue=np.average(RMSValues, weights = fracArea)

    return tileRMSValue

#------------------------------------------------------------------------------------------------------------
def completenessByFootprint(selFnCollection, mockSurvey, diagnosticsDir, additionalLabel = ""):
    """Write out average (M, z) grid for all tiles (extNames) given in selFnCollection (a dictionary with 
    keys corresponding to footprints: 'full' is the entire survey), weighted by fraction of total survey area
    within the footprint. We also produce a bunch of other stats and plots to do with completeness versus 
    redshift.
    
    Output is written to files named e.g. diagnosticsDir/MzCompleteness_label.npz, where label is the 
    footprint (key in selFnCollection); 'full' is the default (survey-wide average).
    
    additionalLabel is optional and will be added to the output filenames (use for e.g., tagging with the
    calcCompleteness method used).
    
    """
    
    zBinEdges=np.arange(0.05, 2.1, 0.1)
    zBinCentres=(zBinEdges[:-1]+zBinEdges[1:])/2.
        
    for footprintLabel in selFnCollection.keys():
        print(">>> Survey-averaged results inside footprint: %s ..." % (footprintLabel))
        selFnDictList=selFnCollection[footprintLabel]
        tileAreas=[]
        compMzCube=[]
        completeness=[]
        for selFnDict in selFnDictList:
            tileAreas.append(selFnDict['tileAreaDeg2'])
            massLimit_90Complete=calcMassLimit(0.9, selFnDict['compMz'], mockSurvey, zBinEdges = zBinEdges)
            completeness.append(massLimit_90Complete)
            compMzCube.append(selFnDict['compMz'])
        tileAreas=np.array(tileAreas)
        completeness=np.array(completeness)
        if np.sum(tileAreas) == 0:
            print("... no overlapping area with %s ..." % (footprintLabel))
            continue
        fracArea=tileAreas/np.sum(tileAreas)
        compMzCube=np.array(compMzCube)
        compMz_surveyAverage=np.average(compMzCube, axis = 0, weights = fracArea)

        outFileName=diagnosticsDir+os.path.sep+"MzCompleteness_%s%s.npz" % (footprintLabel, additionalLabel)
        np.savez(outFileName, z = mockSurvey.z, log10M500c = mockSurvey.log10M, 
                 M500Completeness = compMz_surveyAverage)
        
        makeMzCompletenessPlot(compMz_surveyAverage, mockSurvey.log10M, mockSurvey.z, footprintLabel, 
                               diagnosticsDir+os.path.sep+"MzCompleteness_%s%s.pdf" % (footprintLabel, additionalLabel))

        # 90% mass completeness limit and plots
        massLimit_90Complete=np.average(completeness, axis = 0, weights = fracArea)  # agrees with full mass limit map
        makeMassLimitVRedshiftPlot(massLimit_90Complete, zBinCentres, diagnosticsDir+os.path.sep+"completeness90Percent_%s%s.pdf" % (footprintLabel, additionalLabel))
        zMask=np.logical_and(zBinCentres >= 0.2, zBinCentres < 1.0)
        averageMassLimit_90Complete=np.average(massLimit_90Complete[zMask])
        print("... total survey area (after masking) = %.3f sq deg" % (np.sum(tileAreas)))
        print("... survey-averaged 90%% mass completeness limit (z = 0.5) = %.3f x 10^14 MSun" % (massLimit_90Complete[np.argmin(abs(zBinCentres-0.5))]))
        print("... survey-averaged 90%% mass completeness limit (0.2 < z < 1.0) = %.3f x 10^14 MSun" % (averageMassLimit_90Complete))

#------------------------------------------------------------------------------------------------------------
def makeMzCompletenessPlot(compMz, log10M, z, title, outFileName):
    """Makes a (M, z) plot. Here, compMz is a 2d array, and log10M and z are arrays corresponding to the axes.
    
    """
    
    plotSettings.update_rcParams()
    plt.figure(figsize=(9.5,6.5))
    ax=plt.axes([0.11, 0.11, 0.87, 0.80])

    plt.imshow((compMz*100).transpose(), cmap = colorcet.m_rainbow, origin = 'lower', aspect = 0.8)
    
    y_tck=interpolate.splrep(log10M, np.arange(log10M.shape[0]))
    plot_log10M=np.linspace(13.5, 15.5, 9)
    coords_log10M=interpolate.splev(plot_log10M, y_tck)
    labels_log10M=[]
    for lm in plot_log10M:
        labels_log10M.append("%.2f" % (lm))
    plt.yticks(interpolate.splev(plot_log10M, y_tck), labels_log10M)
    plt.ylim(coords_log10M.min(), coords_log10M.max())
    plt.ylabel("log$_{10}$ ($M / M_{\odot}$)")
    
    x_tck=interpolate.splrep(z, np.arange(z.shape[0]))
    plot_z=np.linspace(0.0, 2.0, 11)
    coords_z=interpolate.splev(plot_z, x_tck)
    labels_z=[]
    for lz in plot_z:
        labels_z.append("%.1f" % (lz))
    plt.xticks(interpolate.splev(plot_z, x_tck), labels_z)
    plt.xlim(coords_z.min(), coords_z.max())
    plt.xlabel("$z$")
    
    plt.colorbar(pad = 0.03)
    cbLabel="Completeness (%)" 
    plt.figtext(0.96, 0.52, cbLabel, ha="center", va="center", family = "sans-serif", rotation = "vertical")

    plt.title(title)
    plt.savefig(outFileName)

#------------------------------------------------------------------------------------------------------------
def calcMassLimit(completenessFraction, compMz, mockSurvey, zBinEdges = []):
    """Given a completeness (log10M, z) grid as made by calcCompleteness, return the mass limit (units of 
    10^14 MSun) as a function of redshift. By defualt, the same binning used in the given mockSurvey object -
    this can be overridden by giving zBinEdges.
    
    """
    
    massLimit=np.power(10, mockSurvey.log10M[np.argmin(abs(compMz-completenessFraction), axis = 1)])/1e14
    
    if zBinEdges != []:
        binnedMassLimit=np.zeros(len(zBinEdges)-1)
        for i in range(len(zBinEdges)-1):
            binnedMassLimit[i]=np.average(massLimit[np.logical_and(mockSurvey.z > zBinEdges[i], mockSurvey.z <= zBinEdges[i+1])])
        massLimit=binnedMassLimit
    
    return massLimit
    
#------------------------------------------------------------------------------------------------------------
def calcCompleteness(y0Noise, SNRCut, extName, mockSurvey, scalingRelationDict, tckQFitDict,
                     plotFileName = None, z = None, method = "fast", numDraws = 2000000, numIterations = 100):
    """Calculate completeness as a function of (log10 M500, z) on the mockSurvey grid for assumed fiducial
    cosmology and scaling relation, at the given SNRCut and noise level. Intrinsic scatter in the scaling
    relation is taken into account.
    
    If plotFileName is given, writes 90% completeness plot.
    
    Two methods for doing the calculation are available:
    - "fast"        : using calcPM500 to fill-out the (log10 M, z) grid
    - "Monte Carlo" : using samples drawn from mockSurvey to fill out the (log10 M, z) grid
    
    Adjust numDraws (per iteration) and numIterations to set the noise level of estimate when using the
    "Monte Carlo" method. For the "fast" method, numDraws and numIterations are ignored.
    
    The calculation can be done at a single fixed z, if given.

    Returns 2d array of (log10 M500, z) completeness
    
    """
        
    if z != None:
        zIndex=np.argmin(abs(mockSurvey.z-z))
        zRange=mockSurvey.z[zIndex:zIndex+1]
    else:
        zRange=mockSurvey.z
        
    if method == "Monte Carlo":
        # Monte-carlo sims approach: slow - but can use to verify the other approach below
        halfBinWidth=(mockSurvey.log10M[1]-mockSurvey.log10M[0])/2.0
        binEdges_log10M=(mockSurvey.log10M-halfBinWidth).tolist()+[np.max(mockSurvey.log10M)+halfBinWidth]
        halfBinWidth=(mockSurvey.z[1]-mockSurvey.z[0])/2.0
        binEdges_z=(zRange-halfBinWidth).tolist()+[np.max(zRange)+halfBinWidth]
        allMz=np.zeros([mockSurvey.clusterCount.shape[1], mockSurvey.clusterCount.shape[0]])
        detMz=np.zeros([mockSurvey.clusterCount.shape[1], mockSurvey.clusterCount.shape[0]])
        for i in range(numIterations):
            tab=mockSurvey.drawSample(y0Noise, scalingRelationDict, tckQFitDict, extName = extName, 
                                    SNRLimit = SNRCut, applySNRCut = False, z = z, numDraws = numDraws)
            allMz=allMz+np.histogram2d(np.log10(tab['true_M500']*1e14), tab['redshift'], [binEdges_log10M, binEdges_z])[0]
            detMask=np.greater(tab['fixed_y_c']*1e-4, y0Noise*SNRCut)
            detMz=detMz+np.histogram2d(np.log10(tab['true_M500'][detMask]*1e14), tab['redshift'][detMask], [binEdges_log10M, binEdges_z])[0]
        mask=np.not_equal(allMz, 0)
        compMz=np.ones(detMz.shape)
        compMz[mask]=detMz[mask]/allMz[mask]
        compMz=compMz.transpose()
        #astImages.saveFITS("test_compMz_MC_5000.fits", compMz.transpose(), None)

    elif method == "fast":
        # This is like the old approach... 50% completeness is by definition where P(M500) is centred
        tenToA0, B0, Mpivot, sigma_int=[scalingRelationDict['tenToA0'], scalingRelationDict['B0'], 
                                        scalingRelationDict['Mpivot'], scalingRelationDict['sigma_int']]
        compMz=np.zeros([mockSurvey.log10M.shape[0], mockSurvey.z.shape[0]])
        for k in range(len(mockSurvey.z)):
            zk=mockSurvey.z[k]
            P=simsTools.calcPM500(y0Noise*SNRCut, y0Noise, zk, 0.0, tckQFitDict[extName], mockSurvey, 
                                tenToA0 = tenToA0, B0 = B0, Mpivot = Mpivot, sigma_int = sigma_int, 
                                applyMFDebiasCorrection = True, fRelWeightsDict = {148.0: 1.0}, 
                                return2D=False)
            cumP=np.zeros(mockSurvey.log10M.shape)
            for i in range(len(mockSurvey.log10M)):
                cumP[i]=np.trapz(P[:i], mockSurvey.log10M[:i])
            compMz[:, k]=cumP
        compMz=compMz.transpose()
        #astImages.saveFITS("test_compMz_fastApproach_withoutMFDebiasCorr_fixedy0ErrInPM500_convScatterInPM500.fits", compMz.transpose(), None)
    
    else:
        raise Exception("calcCompleteness only has 'fast' and 'Monte Carlo' methods available")
            
    if plotFileName != None:
        # Calculate 90% completeness as function of z
        zBinEdges=np.arange(0.05, 2.1, 0.1)
        zBinCentres=(zBinEdges[:-1]+zBinEdges[1:])/2.
        massLimit_90Complete=calcMassLimit(0.9, compMz, mockSurvey, zBinEdges = zBinEdges)
        zMask=np.logical_and(zBinCentres >= 0.2, zBinCentres < 1.0)
        averageMassLimit_90Complete=np.average(massLimit_90Complete[zMask])
        makeMassLimitVRedshiftPlot(massLimit_90Complete, zBinCentres, plotFileName, 
                                   title = "%s: $M_{\\rm 500c}$ / $10^{14}$ M$_{\odot}$ > %.2f (0.2 < $z$ < 1)" % (extName, averageMassLimit_90Complete)) 
            
    return compMz
      
#------------------------------------------------------------------------------------------------------------
def makeMassLimitMap(SNRCut, z, extName, photFilterLabel, mockSurvey, scalingRelationDict, tckQFitDict, 
                     diagnosticsDir):
    """Makes a map of 90% mass completeness (for now, this fraction is fixed).
    
    NOTE: The map here is "downsampled" (i.e., binned) in terms of noise resolution - okay for display 
    purposes, may not be for analysis using the map (if anyone should decide to). Much quicker and saves
    on memory issues (the y0 noise estimates have uncertainties themselves anyway).
    
    """
    
    # Get the stuff we need...
    RMSImg=pyfits.open(diagnosticsDir+os.path.sep+"RMSMap_%s#%s.fits" % (photFilterLabel, extName))
    RMSMap=RMSImg[0].data
    wcs=astWCS.WCS(RMSImg[0].header, mode = 'pyfits')
    RMSTab=getRMSTab(extName, photFilterLabel, diagnosticsDir)
    
    # Downsampling in noise resolution
    stepSize=0.001*1e-4
    binEdges=np.arange(RMSTab['y0RMS'].min(), RMSTab['y0RMS'].max()+stepSize, stepSize)
    y0Binned=[]
    fracAreaBinned=[]
    binMins=[]
    binMaxs=[]
    for i in range(len(binEdges)-1):
        mask=np.logical_and(RMSTab['y0RMS'] > binEdges[i], RMSTab['y0RMS'] <= binEdges[i+1])
        if mask.sum() > 0:
            y0Binned.append(np.average(RMSTab['y0RMS'][mask]))
            fracAreaBinned.append(np.sum(RMSTab['y0RMS'][mask]))
            binMins.append(binEdges[i])
            binMaxs.append(binEdges[i+1])
            
    # Fill in blocks in map for each RMS value
    outFileName=diagnosticsDir+os.path.sep+"massLimitMap_z%s#%s.fits" % (str(z).replace(".", "p"), extName)
    if os.path.exists(outFileName) == False:
        massLimMap=np.zeros(RMSMap.shape)
        count=0
        t0=time.time()
        # New
        for y0Noise, binMin, binMax in zip(y0Binned, binMins, binMaxs):
            count=count+1
            print(("... %d/%d (%.3e) ..." % (count, len(y0Binned), y0Noise)))
            compMz=calcCompleteness(y0Noise, SNRCut, extName, mockSurvey, scalingRelationDict, tckQFitDict,
                                    z = z)
            mapMask=np.logical_and(RMSMap > binMin, RMSMap <= binMax)
            massLimMap[mapMask]=mockSurvey.log10M[np.argmin(abs(compMz-0.9))]
            # Does this stop memory leak?
            #massLimMap[RMSMap == y0Noise]=mockSurvey.log10M[np.argmin(abs(calcCompleteness(y0Noise, SNRCut, extName, mockSurvey, scalingRelationDict, tckQFitDict, z = z)-0.9))]
        t1=time.time()        
        massLimMap=np.power(10, massLimMap)/1e14
        astImages.saveFITS(outFileName, massLimMap, wcs)
        # Old 
        #for y0Noise in RMSTab['y0RMS']:
            #count=count+1
            #print(("... %d/%d (%.3e) ..." % (count, len(RMSTab), y0Noise)))
            ##compMz=calcCompleteness(y0Noise, SNRCut, extName, mockSurvey, scalingRelationDict, tckQFitDict,
                                    ##z = z)
            ## Does this stop memory leak?
            #massLimMap[RMSMap == y0Noise]=mockSurvey.log10M[np.argmin(abs(calcCompleteness(y0Noise, SNRCut, extName, mockSurvey, scalingRelationDict, tckQFitDict, z = z)-0.9))]
        #t1=time.time()        
        #massLimMap=np.power(10, massLimMap)/1e14
        #astImages.saveFITS(outFileName, massLimMap, wcs)
    
    # This is sanity checking survey average completeness... 
    # Summary: agrees within 0.2% on average, but some tiles out by up to 3%
    #massImg=pyfits.open(outFileName)
    #massLimMap=massImg[0].data
    #areaImg=pyfits.open(diagnosticsDir+os.path.sep+"areaMask#%s.fits" % (extName))
    #areaMap=areaImg[0].data
    #wcs=astWCS.WCS(areaImg[0].header, mode = 'pyfits')
    #areaMapSqDeg=(mapTools.getPixelAreaArcmin2Map(areaMap, wcs)*areaMap)/(60**2)
    #fromMassLimMap=np.average(massLimMap, weights = areaMapSqDeg/areaMapSqDeg.sum())
    #y0Noise=calcTileWeightedAverageNoise(extName, photFilterLabel, diagnosticsDir)
    #fitTab=calcCompleteness(y0Noise, SNRCut, extName, mockSurvey, parDict['massOptions'], tckQFitDict, diagnosticsDir,
                            #fitTabFileName = diagnosticsDir+os.path.sep+"selFn_fitTab#%s.fits" % (extName))
    #fromSurveyAverage=np.power(10, fitTab['log10MLimit_90%'][np.where(fitTab['z'] == 0.5)])/1e14
    #print("... extName = %s: limit from mass limit map / survey average = %.3f ..." % (extName, fromMassLimMap/fromSurveyAverage))

#------------------------------------------------------------------------------------------------------------
def makeMassLimitVRedshiftPlot(massLimit_90Complete, zRange, outFileName, title = None):
    """Write a plot of 90%-completeness mass limit versus z, adding a spline interpolation.
    
    """
    
    plotSettings.update_rcParams()
    plt.figure(figsize=(9,6.5))
    if title == None:
        ax=plt.axes([0.10, 0.11, 0.87, 0.86])
    else:
        ax=plt.axes([0.10, 0.11, 0.87, 0.80])
    tck=interpolate.splrep(zRange, massLimit_90Complete)
    plotRange=np.linspace(0, 2, 100)
    plt.plot(plotRange, interpolate.splev(plotRange, tck), 'k-')
    plt.plot(zRange, massLimit_90Complete, 'D', ms = 8)
    plt.xlabel("$z$")
    plt.ylim(0.5, 8)
    plt.xticks(np.arange(0, 2.2, 0.2))
    plt.xlim(0, 2)
    labelStr="$M_{\\rm 500c}$ (10$^{14}$ M$_{\odot}$) [90% complete]"
    plt.ylabel(labelStr)
    if title != None:
        plt.title(title)
    plt.savefig(outFileName)
    plt.close()   
    
#------------------------------------------------------------------------------------------------------------
def cumulativeAreaMassLimitPlot(z, diagnosticsDir):
    """Make a cumulative plot of 90%-completeness mass limit versus survey area, at the given redshift.
    
    """

    fileList=glob.glob(diagnosticsDir+os.path.sep+"massLimitMap_z%s#*.fits" % (str(z).replace(".", "p")))
    extNames=[]
    for f in fileList:
        extNames.append(f.split("#")[-1].split(".fits")[0])
    extNames.sort()
    
    # NOTE: We can avoid this if we add 'areaDeg2' column to selFn_mapFitTab_*.fits tables
    allLimits=[]
    allAreas=[]
    for extName in extNames:
        massLimImg=pyfits.open(diagnosticsDir+os.path.sep+"massLimitMap_z%s#%s.fits" % (str(z).replace(".", "p"), extName))
        massLimMap=massLimImg[0].data
        areaImg=pyfits.open(diagnosticsDir+os.path.sep+"areaMask#%s.fits" % (extName))
        areaMap=areaImg[0].data
        wcs=astWCS.WCS(areaImg[0].header, mode = 'pyfits')
        areaMapSqDeg=(mapTools.getPixelAreaArcmin2Map(areaMap, wcs)*areaMap)/(60**2)
        limits=np.unique(massLimMap).tolist()
        areas=[]
        for l in limits:
            areas.append(areaMapSqDeg[np.where(massLimMap == l)].sum())
        allLimits=allLimits+limits
        allAreas=allAreas+areas
    tab=atpy.Table()
    tab.add_column(atpy.Column(allLimits, 'MLim'))
    tab.add_column(atpy.Column(allAreas, 'areaDeg2'))
    tab.sort('MLim')

    plotSettings.update_rcParams()

    # Full survey plot
    plt.figure(figsize=(9,6.5))
    ax=plt.axes([0.155, 0.12, 0.82, 0.86])
    plt.minorticks_on()
    plt.plot(tab['MLim'], np.cumsum(tab['areaDeg2']), 'k-')
    plt.ylabel("survey area < $M_{\\rm 500c}$ limit (deg$^2$)")
    plt.xlabel("$M_{\\rm 500c}$ (10$^{14}$ M$_{\odot}$) [90% complete]")
    labelStr="total survey area = %.0f deg$^2$" % (np.cumsum(tab['areaDeg2']).max())
    plt.ylim(0.1, 1.01*np.cumsum(tab['areaDeg2']).max())
    plt.xlim(1.5, 10)
    plt.figtext(0.2, 0.9, labelStr, ha="left", va="center")
    plt.savefig(diagnosticsDir+os.path.sep+"cumulativeArea_massLimit_z%s.pdf" % (str(z).replace(".", "p")))
    plt.close()
    
    # Deepest 20%
    totalAreaDeg2=tab['areaDeg2'].sum()
    deepTab=tab[np.where(np.cumsum(tab['areaDeg2']) < 0.2 * totalAreaDeg2)]
    plt.figure(figsize=(9,6.5))
    ax=plt.axes([0.155, 0.12, 0.82, 0.86])
    plt.minorticks_on()
    plt.plot(deepTab['MLim'], np.cumsum(deepTab['areaDeg2']), 'k-')
    plt.ylabel("survey area < $M_{\\rm 500c}$ limit (deg$^2$)")
    plt.xlabel("$M_{\\rm 500c}$ (10$^{14}$ M$_{\odot}$) [90% complete]")
    labelStr="area of deepest 20%% = %.0f deg$^2$" % (np.cumsum(deepTab['areaDeg2']).max())
    plt.ylim(0.1, 1.01*np.cumsum(deepTab['areaDeg2']).max())
    plt.xlim(1.5, deepTab['MLim'].max())
    plt.figtext(0.2, 0.9, labelStr, ha="left", va="center")
    plt.savefig(diagnosticsDir+os.path.sep+"cumulativeArea_massLimit_z%s_deepest20Percent.pdf" % (str(z).replace(".", "p")))
    plt.close()
        
#------------------------------------------------------------------------------------------------------------
def makeFullSurveyMassLimitMapPlot(z, diagnosticsDir):
    """Reprojects tile mass limit maps onto the full map pixelisation, then makes a plot and saves a
    .fits image.
    
    NOTE: we hardcoded the pixelisation in here for now...
    
    """

    # Making the full res reprojected map takes ~1650 sec
    outFileName=diagnosticsDir+os.path.sep+"reproj_massLimitMap_z%s.fits" % (str(z).replace(".", "p"))
    if os.path.exists(outFileName) == False:

        print(">>> Making reprojected full survey mass limit map:")
        
        # Downsampled
        h=pyfits.Header()
        h['NAXIS']=2
        h['NAXIS1']=10800
        h['NAXIS2']=2580
        h['WCSAXES']=2
        h['CRPIX1']=5400.25
        h['CRPIX2']=1890.25
        h['CDELT1']=-0.0333333333333332
        h['CDELT2']=0.0333333333333332
        h['CUNIT1']='deg'
        h['CUNIT2']='deg'
        h['CTYPE1']='RA---CAR'
        h['CTYPE2']='DEC--CAR'
        h['CRVAL1']=0.0
        h['CRVAL2']=0.0
        h['LONPOLE']=0.0
        h['LATPOLE']=90.0
        h['RADESYS']='ICRS'
        # Full res
        #h['NAXIS']=2
        #h['NAXIS1']=43200
        #h['NAXIS2']=10320
        #h['WCSAXES']=2
        #h['CRPIX1']=21601.0
        #h['CRPIX2']=7561.0
        #h['CDELT1']=-0.0083333333333333
        #h['CDELT2']=0.0083333333333333
        #h['CUNIT1']='deg'
        #h['CUNIT2']='deg'
        #h['CTYPE1']='RA---CAR'
        #h['CTYPE2']='DEC--CAR'
        #h['CRVAL1']=0.0
        #h['CRVAL2']=0.0
        #h['LONPOLE']=0.0
        #h['LATPOLE']=90.0
        #h['RADESYS']='ICRS'
        wcs=astWCS.WCS(h, mode = 'pyfits')
        
        reproj=np.zeros([wcs.header['NAXIS2'], wcs.header['NAXIS1']])
        sumPix=np.zeros([wcs.header['NAXIS2'], wcs.header['NAXIS1']])
        
        fileList=glob.glob(diagnosticsDir+os.path.sep+"massLimitMap_z%s#*.fits" % (str(z).replace(".", "p")))
        extNames=[]
        for f in fileList:
            extNames.append(f.split("#")[-1].split(".fits")[0])
        extNames.sort()
        
        t0=time.time()
        for extName in extNames:
            print("... %s ..." % (extName))
            img=pyfits.open(diagnosticsDir+os.path.sep+"massLimitMap_z%s#%s.fits" % (str(z).replace(".", "p"), extName))
            data=img[0].data
            imgWCS=astWCS.WCS(img[0].header, mode = 'pyfits')
            areaImg=pyfits.open(diagnosticsDir+os.path.sep+"areaMask#%s.fits" % (extName))
            areaMap=areaImg[0].data
            for y in range(data.shape[0]):
                for x in range(data.shape[1]):
                    outRADeg, outDecDeg=imgWCS.pix2wcs(x, y)
                    inX, inY=wcs.wcs2pix(outRADeg, outDecDeg)
                    # Once, this returned infinity...
                    try:
                        inX=int(round(inX))
                        inY=int(round(inY))
                    except:
                        continue
                    # handle the overlap regions, which were zeroed
                    if areaMap[y, x] > 0:   
                        reproj[inY, inX]=reproj[inY, inX]+data[y, x]
                        sumPix[inY, inX]=sumPix[inY, inX]+1.0
        t1=time.time()
        #print(t1-t0)
        astImages.saveFITS(outFileName, reproj/sumPix, wcs)
    
    # Make plot
    img=pyfits.open(outFileName)
    reproj=np.nan_to_num(img[0].data)
    #reproj=np.ma.masked_where(reproj == np.nan, reproj)
    reproj=np.ma.masked_where(reproj <1e-6, reproj)
    wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
    plotSettings.update_rcParams()
    fontSize=20.0
    figSize=(16, 5.7)
    axesLabels="sexagesimal"
    axes=[0.08,0.15,0.91,0.88]
    cutLevels=[2, 9]
    colorMapName=colorcet.m_rainbow
    fig=plt.figure(figsize = figSize)
    p=astPlots.ImagePlot(reproj, wcs, cutLevels = cutLevels, title = None, axes = axes, 
                         axesLabels = axesLabels, colorMapName = colorMapName, axesFontFamily = 'sans-serif', 
                         RATickSteps = {'deg': 30.0, 'unit': 'h'}, decTickSteps = {'deg': 20.0, 'unit': 'd'},
                         axesFontSize = fontSize)
    cbLabel="$M_{\\rm 500c}$ (10$^{14}$ M$_{\odot}$) [90% complete]"
    cbShrink=0.7
    cbAspect=40
    cb=plt.colorbar(p.axes.images[0], ax = p.axes, orientation="horizontal", fraction = 0.05, pad = 0.18, 
                    shrink = cbShrink, aspect = cbAspect)
    plt.figtext(0.53, 0.04, cbLabel, size = 20, ha="center", va="center", fontsize = fontSize, family = "sans-serif")
    plt.savefig(outFileName.replace(".fits", ".pdf"))
    plt.close()
