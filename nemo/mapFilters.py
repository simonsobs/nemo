# -*- coding: utf-8 -*-
"""Filter classes are defined in this module.

There are two main base classes of filter - WienerFilter and MatchedFilter. These are both derived from the
overall base class for everything: MapFilter.

There are also base classes corresponding to filters with different signal templates: Beta, Profle, Gaussian.

The actual filters that can be used are derived from these, e.g.:

BetaModelMatchedFilter
ProfileMatchedFilter
GaussianMatchedFilter
ArnaudModelMatchedFilter
BetaModelWienerFilter
ProfileWienerFilter
GaussianWienerFilter
ArnaudModelWienerFilter

"""

import math
from flipper import liteMap
from flipper import fftTools
from astLib import *
import numpy as np
from numpy import fft
import pylab as plt
import os
from scipy import interpolate
from scipy import ndimage
import mapTools
import simsTools
import photometry
import catalogTools
import simsTools
import plotSettings
import gnfw
import pyfits
import copy
import sys
import glob
import itertools
import pyximport; pyximport.install()
import nemoCython
import time
import IPython

#------------------------------------------------------------------------------------------------------------
class MapFilter:
    """Generic map filter base class. Defines common interface.
    
    """
    def __init__(self, label, unfilteredMapsDictList, paramsDict, writeFilter = False, forceRebuild = False,
                 outDir = ".", diagnosticsDir = None):
        """Initialises a MapFilter. unfilteredMapsDictList describes the input maps, paramsDict describes
        the filter options. The convention is that single frequency only filters (e.g. GaussianWienerFilter)
        only operate on the first map in the unfilteredMapDictList.
        
        label is used to store the filter to save calculating it again. Filters are stored under filters/
        dir. If a filter already exists under this dir with the given label, it is loaded rather than
        recalculated. To force recalculating the filter, set forceRebuild == True.
                        
        """
        
        self.label=label
        self.params=paramsDict
        
        # Prepare all the unfilteredMaps (in terms of cutting sections, adding noise etc.)
        for mapDict in unfilteredMapsDictList:           
            mapDict=mapTools.preprocessMapDict(mapDict, diagnosticsDir = diagnosticsDir)
        self.wcs=mapDict['wcs']

        # Sanity check that all maps are the same dimensions
        shape=unfilteredMapsDictList[0]['data'].shape
        for mapDict in unfilteredMapsDictList:
            if mapDict['data'].shape != shape:
                raise Exception, "Maps at different frequencies have different dimensions!"
        
        self.unfilteredMapsDictList=unfilteredMapsDictList
                                
        # Set up storage if necessary, build this filter if not already stored
        self.diagnosticsDir=diagnosticsDir
        if os.path.exists(outDir) == False:
            os.makedirs(outDir)
        self.filterFileName=self.diagnosticsDir+os.path.sep+label+".fits"
               
        
    def makeRadiansMap(self):
        """Makes a map of distance in radians from centre, based on dimensions of sciMap.data.
        
        """
        
        mapDict=self.unfilteredMapsDictList[0]
        
        x0=mapDict['data'].shape[1]/2
        y0=mapDict['data'].shape[0]/2
        ra0, dec0=self.wcs.pix2wcs(x0, y0)
        ra1, dec1=self.wcs.pix2wcs(x0+1, y0+1)
        self.degPerPixX=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
        self.degPerPixY=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
        
        # Real space map og angular distance from centre in radians, used in making filters and beam
        xRadRange=np.array([np.arange(-mapDict['data'].shape[1]/2, mapDict['data'].shape[1]/2, \
                                    dtype=np.float64)*np.radians(self.degPerPixX)]*mapDict['data'].shape[0])
        yRadRange=np.array([np.arange(-mapDict['data'].shape[0]/2, mapDict['data'].shape[0]/2, \
                                    dtype=np.float64)*np.radians(self.degPerPixY)]*mapDict['data'].shape[1]).transpose()
        rRadRange=np.sqrt(xRadRange**2+yRadRange**2)
        self.radiansMap=rRadRange
        
        
    def buildAndApply(self):
        """Builds and applies the filter to the unfiltered map(s). Returns a dictionary, containing
        keys 'data', 'wcs', 'weights', 'obsFreqGHz'. If this routine converts to yc in place, the latter is set
        to 'yc'.
    
        """
        
        raise Exception, "Called a base filter class without a buildAndApply() function implemented."
        return None
    
    
    def makeForegroundsPower(self):
        """Makes a flipper power2d object from given foregrounds template file.
        This is assumed from Ryan, and is given in units l, dT (micro K).
        
        This should be called after setUpFFTStuff() has already been called
        Returns a flipper power2D object
        
        """
                
        inFile=file(self.params['foregroundsTemplate'], "r") 
        #inFile=file("/home/matty/Astro_Software/flipper-0.1.0/params/wmap_3year_toKevin_lensed_v3_lensedCls.dat", "r")
        lines=inFile.readlines()
        inFile.close()
        l=[]
        llClsMicroK2=[] # l(l+1)Cl/2pi, in uK^2
        for line in lines:
            if line[0] != "#" and len(line) > 3:
                bits=line.split()
                l.append(float(bits[0]))
                llClsMicroK2.append(float(bits[1]))
        l=np.array(l)
        
        # Square here because Ryan had them square rooted - if using the CMB power spec from flipper, take this off
        llClsMicroK2=np.array(llClsMicroK2)**2 
        
        # Note we square back up below, as we do actually want a power spectrum, if we're running from Ryan's
        # stuff
        ll=np.ravel(self.fSciMap.modLMap)
        fgPowerMap=np.zeros(self.fSciMap.kMap.shape)
        interpolator=interpolate.interp1d(l, llClsMicroK2, fill_value=0.0, bounds_error=False)
        kk=interpolator(ll)
        fgPowerMap=np.reshape(kk, [self.fSciMap.Ny, self.fSciMap.Nx])
        fgPower=fftTools.powerFromFFT(self.fSciMap)
        
        # Convert to raw Cls
        fgPower.powerMap=(fgPowerMap*np.pi*2)/(fgPower.modLMap*(fgPower.modLMap+1))
        fgPower.powerMap=np.nan_to_num(fgPower.powerMap) # division by l(l+1) above causes nan where l == 0
        
        # Check if can recover 1d power spectrum from 2d
        #results=fgPower.binInAnnuli("/home/matty/Astro_Software/flipper-0.1.0/params/BIN_100_200")
        #plt.plot(results[2], (results[2]*(results[2]+1)*results[3])/(np.pi*2))
        #plt.plot(l, llClsMicroK2)
        
        return fgPower
               
        
    def loadFilter(self, filterFileName):
        """Loads in a previously saved filter
        
        """
        img=pyfits.open(filterFileName)
        self.G=img[0].data
        
    
    def makeBeamPowerSpectrum(self):
        """Makes beam image power spectrum. This is normalised.
        
        """
        
        beamFWHMRad=np.radians(1.4/60.0)
        beamSigmaRad=beamFWHMRad/np.sqrt(8.0*np.log(2.0))       
        beam=np.exp(-(self.radiansMap**2)/(2*beamSigmaRad**2))            
        self.beamPowerSpectrum=fft.fft2(beam)
        self.beamPowerSpectrum=self.beamPowerSpectrum/self.beamPowerSpectrum.max()
        self.beamPowerSpectrum=self.beamPowerSpectrum*self.beamPowerSpectrum.conj()
    
    
    def besselTransformSmoothFunction(self, x):
        """Smoothing function used by the Bessel transform routine.
        
        """
        
        q=0.3
        smoothFn=np.zeros(x.shape)
        indexLow=np.less_equal(x, 0.5)
        indexHigh=np.greater_equal(x, 1.0)
        indexMed=np.logical_and(np.greater(x, 0.5), np.less(x, 1.0))

        smoothFn[indexLow]=1.0
        smoothFn[indexHigh]=0.0
        smoothFn[indexMed]=(1.+np.tanh(q/(x*(indexMed)-0.5)+q/(x*(indexMed))))/2.0
        
        return smoothFn
        

    def besselTransform(self, profile, lDim, thetaMax, oversamplingFactor):
        """Does the Bessel transform of the given profile in real space (i.e. function of theta in radians).
        
        thetaMax in radians too.
        
        Returns [bessel transformed array, l range]
        
        """
        
        # NOTE: IDL replicate(0, l_dim) => np.zeros(l_dim)
        #       IDL indgen(n_dim) => np.arange(0, n_dim)
        # x -> theta
        thetaDim=profile.shape[0]
        lMax=(lDim*(np.pi/2.0)/oversamplingFactor)/thetaMax
        thetaArray=(np.arange(0, thetaDim, dtype=float))/(thetaDim)*thetaMax
        lArray=(np.arange(0, lDim, dtype=float))/(lDim)*lMax
        thetaUnitArray=(np.arange(0, thetaDim, dtype=float))/(thetaDim)
        lUnitArray=(np.arange(0, lDim, dtype=float))/(lDim)
        dTheta=thetaMax/float(thetaDim)

        # Kavi's way
        #besselTransformedArray=np.zeros(lDim, dtype=float)
        #lSmoothed=self.besselTransformSmoothFunction(lUnitArray)
        #thetaSmoothed=self.besselTransformSmoothFunction(thetaUnitArray)              
        #for i in range(lDim):
            #print "... %d/%d ..." % (i, lDim)
            #btSum=0.0
            #for j in range(thetaDim):
                #z=lArray[i]*thetaArray[j]
                #btSum=btSum+special.j0(z)*thetaArray[j]*thetaSmoothed[j]*dTheta*profile[j]*lSmoothed[i]
            #besselTransformedArray[i]=btSum

        # Fast way, looks identical to me
        besselTransformedArray=np.zeros(lDim, dtype=float)
        lSmoothed=self.besselTransformSmoothFunction(lUnitArray)
        thetaSmoothed=self.besselTransformSmoothFunction(thetaUnitArray)  
        for i in range(lDim):
            tenPercent=lArray.shape[0]/10
            for j in range(0,11):
                if i == j*tenPercent:
                    print "... "+str(j*10)+"% complete ..."
            besselTransformedArray[i]=integrate.simps(profile*special.j0(lArray[i]*thetaArray)*thetaArray*thetaSmoothed, thetaArray)

        return [2*np.pi*besselTransformedArray, lArray]
                

    def noisePowerFrom4WayMaps(self, weights, obsFreqGHz):
        """This makes a noise power spectrum from combinations of the 4-way difference maps, and adds on
        noise from the CMB, point sources etc. from l-space templates.
        
        """
                
        noiseParams=self.params['noiseParams']       
        
        powerPath=self.diagnosticsDir+os.path.sep+"4WayNoise_power_%s.fits" % (int(obsFreqGHz))
        if os.path.exists(powerPath) == False:
            ar1MapFileNames=glob.glob(noiseParams['4WayMapsPattern_%s' % int(obsFreqGHz)])
            ar1MapFileNames.sort()
            
            # Run all null map combinations, take average
            # Corrected from previous version... scaling was all wrong
            combinations=[[0, 1, 2, 3], [0, 2, 1, 3], [0, 3, 1, 2]]
            sumNullMap=None
            countNullMaps=0
            for c in combinations:
                mapFileNames=(np.array(ar1MapFileNames)[c]).tolist()
                countNullMaps=countNullMaps+1
                nullMap=None
                splitCount=0
                for f in mapFileNames:
                    img=pyfits.open(f)
                    wcs=astWCS.WCS(f)
                    if sumNullMap == None:
                        sumNullMap=np.zeros(img[0].data.shape)
                    if nullMap == None:
                        nullMap=img[0].data
                    elif splitCount < 2:
                        nullMap=nullMap+img[0].data
                    else:
                        nullMap=nullMap-img[0].data
                    splitCount=splitCount+1
                nullMap=nullMap/8.0   # scales to same noise level as in full map
                sumNullMap=sumNullMap+nullMap
            nullMap=sumNullMap/countNullMaps

            # Corrected from previous version... scaling was all wrong
            #nullMap=None
            #splitCount=0
            #for f in ar1MapFileNames:
                #img=pyfits.open(f)
                #wcs=astWCS.WCS(f)
                #if nullMap == None:
                    #nullMap=img[0].data
                #elif splitCount < 2:
                    #nullMap=nullMap+img[0].data
                #else:
                    #nullMap=nullMap-img[0].data
                #splitCount=splitCount+1
            #nullMap=nullMap/8.0   # scales to same noise level as in full map
            
            # Trim, subtract background if needed
            if 'RADecSection' in self.unfilteredMapsDictList[0].keys() and self.unfilteredMapsDictList[0]['RADecSection'] != None:
                clipSection=True
            else:
                clipSection=False
            if 'backgroundSubtraction' in self.unfilteredMapsDictList[0].keys() and self.unfilteredMapsDictList[0]['backgroundSubtraction'] == True:
                bckSub=True      
            else:
                bckSub=False
            if clipSection == True:
                RAMin, RAMax, decMin, decMax=self.unfilteredMapsDictList[0]['RADecSection']
                clip=mapTools.clipUsingRADecCoords(nullMap, wcs, RAMin, RAMax, decMin, decMax)
                mapData=clip['data']
                mapWCS=clip['wcs']
            if bckSub == True:
                mapData=mapTools.subtractBackground(mapData, mapWCS)                
                
            # Make some noise - have to treat same way as data
            lm=liteMap.liteMapFromDataAndWCS(mapData*np.sqrt(weights/weights.max()), mapWCS)
            apodlm=lm.createGaussianApodization(pad=10, kern=5)
            apod=apodlm.data
            plm=liteMap.liteMapFromDataAndWCS(mapData*apod, mapWCS) # Note: apply apodization here
            power=fftTools.powerFromLiteMap(plm)
            powerMap=power.powerMap
            
            # Set noise to infinite where we have horizontal striping (although this doesn't look to be doing anything
            lxMap=np.array([power.lx]*power.ly.shape[0])
            mask=np.less(abs(lxMap), 100)
            powerMap[mask]=1e30
            #lyMap=np.array([power.ly]*power.lx.shape[1]).transpose()
            #mask=np.less(abs(lyMap), 100)
            #powerMap[mask]=1e30
            
            # Or just admit we're not going to cope with any really large scale power... 
            mask=np.less(abs(power.modLMap), 1000)
            powerMap[mask]=1e30
            
            astImages.saveFITS(powerPath, fft.fftshift(powerMap), None)
        else:
            img=pyfits.open(powerPath)
            powerMap=fft.fftshift(img[0].data)
                    
        # Now add in noise from CMB power spec and point sources ...
        foregroundsPower=self.makeForegroundsPower(obsFreqGHz)
        

        
        # Combine it all together - has to be returned as a flipper Power2D object
        NP=foregroundsPower
        try:
            NP.powerMap=NP.powerMap+powerMap
        except:
            print "Hmm? 4WayNoise"
            ipshell()
            sys.exit()
                
        return NP
        

    def makeForegroundsPower(self, obsFreqGHz):
        """Returns a Power2D object with foregrounds power from the CMB (using the WMAP spectrum included in
        flipper) plus the ACT point source contribution measured in Fowler et al. 2010, hardcoded. Need
        to add code for 220 GHz point source contribution.
                
        """
                
        # CMB from flipper WMAP spectrum file
        inFileName=os.environ['FLIPPER_DIR']+os.path.sep+"params"+os.path.sep+"wmap_3year_toKevin_lensed_v3_lensedCls.dat"
        inFile=file(inFileName, "r")
        lines=inFile.readlines()
        inFile.close()
        l=[]
        llClsMicroK2=[] # l(l+1)Cl/2pi, in uK^2
        for line in lines:
            if line[0] != "#" and len(line) > 3:
                bits=line.split()
                l.append(float(bits[0]))
                llClsMicroK2.append(float(bits[1]))
        l=np.array(l)
        llClsMicroK2=np.array(llClsMicroK2)
        
        # Convert to Cl
        Cl=(llClsMicroK2*2*np.pi)/(l*(l+1))
        
        # Add in point source power from Fowler et al.
        if obsFreqGHz == 148:
            Ap=11.2 # uK^2
            PSB=Ap*np.power(l/3000, 2)
            PSCl=(PSB*np.pi*2)/(l*(l+1))
        else:
            print "WARNING: need to add code for point source contribution at other frequencies than 148 GHz"
            PSCl=np.zeros(Cl.shape)
        
        # We should roll off the sources power spectrum in l-space with the beam - actually, this does nothing
        #fwhmDeg=1.4/60.0
        #eightLnTwo=8.0*math.log(2.0) # note: sqrt(8ln2)*sigma=fwhm 
        #beam=np.exp(-1.0*l*(l+1)*math.radians(fwhmDeg)**2/eightLnTwo)
        #beam=(beam*np.pi*2)/(l*(l+1))
        #PSCl=PSCl*beam
        
        # FFT stuff
        templm=liteMap.liteMapFromDataAndWCS(np.ones(self.unfilteredMapsDictList[0]['data'].shape), self.wcs)
        fTempMap=fftTools.fftFromLiteMap(templm)
        
        # Add point source power to CMB, make 2d
        ll=np.ravel(fTempMap.modLMap)
        fgPowerMap=np.zeros(fTempMap.kMap.shape)
        interpolator=interpolate.interp1d(l, Cl+PSCl, fill_value=0.0, bounds_error=False)
        kk=interpolator(ll)
        fgPowerMap=np.reshape(kk, [fTempMap.Ny, fTempMap.Nx])
        fgPower=fftTools.powerFromFFT(fTempMap)
        fgPower.powerMap=np.nan_to_num(fgPower.powerMap)
                
        return fgPower
        
        
    def writeFilteredMap(self, outFileName):
        """Writes out the filtered map as a fits file.
        
        """
        astImages.saveFITS(outFileName, self.filteredMap, self.wcs)


    def makeRealSpaceFilterProfile(self):
        """Makes a 1d real-space profile of the filter.
        
        Returns profile, arcminRange
        """
        wcs=self.wcs
        
        realSpace=fft.ifft2(self.G).real
        realSpace=fft.fftshift(realSpace)
        
        x0=realSpace.shape[1]/2
        y0=realSpace.shape[0]/2
        prof=realSpace[y0, x0:]/realSpace[y0, x0:].max() # normalise for plot
        arcminRange=np.arange(0, prof.shape[0])*self.degPerPixX*60.0
        
        return prof, arcminRange
    
        
    def saveRealSpaceFilterProfile(self):
        """Saves a real-space profile of the filter out as a .png plot and 2d .fits image.
        
        """
        
        prof, arcminRange=self.makeRealSpaceFilterProfile()
        
        # Measure characteristic FWHM
        tck=interpolate.splrep(arcminRange, prof)
        FWHMArcmin=interpolate.splev([0.5], tck)*2 # *2 because otherwise would be half width
        
        fig=plt.figure(num=2, figsize=(8,8))
        fig.canvas.set_window_title('Filter Profile in Real Space')
        plt.clf()
        plt.title("Filter Profile %s" % (self.label))
        plt.ylabel("Amplitude")
        plt.xlabel("$\\theta$ (arcmin)")
        plt.plot(arcminRange, prof)
        plt.xlim(0, 30.0)
        plt.ylim(prof.min(), prof.max()*1.1)
        plt.savefig(self.diagnosticsDir+os.path.sep+"realSpaceProfile1d_"+self.label+".png")
        plt.close()
        
        # Save 2D realspace filter image too
        #astImages.saveFITS(self.diagnosticsDir+os.path.sep+"realSpaceProfile2d_"+self.label+".fits", \
                                #realSpace, wcs)
    

#------------------------------------------------------------------------------------------------------------
class WienerFilter(MapFilter):
    """Generic class for Wiener filtering. All derivatives of this class must provide a 
    makeSignalTemplateMap() that returns a liteMap. This class includes the routines for various ways of 
    handling the noise that goes into the denominator (e.g. l-space templates).
    
    """
    
    def buildAndApply(self):
        
        print ">>> Building filter %s ..." % (self.label)

        if 'iterations' not in self.params.keys():
            self.params['iterations']=1
        iterations=self.params['iterations']
        clusterCatalog=[]
        for i in range(iterations):
            
            #print "... iteration %d ..." % (i) 
            
            # Filter both maps with the same scale filter
            filteredMaps={}
            for mapDict in self.unfilteredMapsDictList:   
                
                # Iterative filtering stuff currently disabled
                # Mask out clusters so they don't go into the noise component of the filter
                #maskedDict=mapTools.maskOutSources(mapDict['data'], self.wcs, clusterCatalog, 
                                                   #radiusArcmin = self.params['maskingRadiusArcmin'], 
                                                   #mask = 'whiteNoise')
                #maskedData=maskedDict['data']
                #if self.diagnosticsDir != None:
                    #outPath=self.diagnosticsDir+os.path.sep+"iteration_%d_%s_mask.fits" % (i, self.label)
                    #astImages.saveFITS(outPath, maskedDict['mask'], self.wcs)
                
                maskedData=mapDict['data']
                
                # FFT - note using flipper 0.1.3 Gaussian apodization
                # The apodization here is only used if we're estimation noise power directly from the map itself
                # otherwise, the FFT here is just so we have the modThetaMap etc. which we would use if
                # using e.g. beta model templates where we're setting our own normalisation, rather than
                # using e.g. the profile models.
                lm=liteMap.liteMapFromDataAndWCS(maskedData, self.wcs) 
                lm.data=lm.data*np.sqrt(mapDict['weights']/mapDict['weights'].max()) # this appears to make no difference
                apodlm=lm.createGaussianApodization()
                lm.data=lm.data*apodlm.data
                fMaskedMap=fftTools.fftFromLiteMap(lm)
                fMaskedMap.modThetaMap=180.0/(fMaskedMap.modLMap+1)       
                                
                # Make noise power spectrum
                if self.params['noiseParams']['method'] == 'dataMap':
                    NP=fftTools.powerFromFFT(fMaskedMap)
                elif self.params['noiseParams']['method'] == '4WayMaps':
                    NP=self.noisePowerFrom4WayMaps(self.weights, mapDict['obsFreqGHz'])
                else:
                    raise Exception, "noise method must be either 'dataMap' or '4WayMaps'"
                
                # Make signal power spectrum
                signalMapDict=self.makeSignalTemplateMap(mapDict['beamFWHMArcmin'], mapDict['obsFreqGHz'])
                signalMap=signalMapDict['signalMap']
                rInnerDeg=signalMapDict['normInnerDeg']
                rOuterDeg=signalMapDict['normOuterDeg']
                powerScaleFactor=signalMapDict['powerScaleFactor']
                SP=fftTools.powerFromLiteMap(signalMap)
                
                # Scale the signal power such that it is equal to the noise power within given annulus, 
                # if signal template has this sort of scaling applied (e.g. beta model case)
                # If it doesn't, we'd better have supplied powerScaleFactor (e.g. profile case)
                if rInnerDeg != None and rOuterDeg != None:
                    mask=np.logical_and(np.greater(fMaskedMap.modThetaMap, rInnerDeg), \
                                        np.less(fMaskedMap.modThetaMap, rOuterDeg))
                    matchNoise=np.sum(NP.powerMap[mask])/np.sum(SP.powerMap[mask])
                    SP.powerMap=SP.powerMap*matchNoise
                elif powerScaleFactor != None:
                    SP.powerMap=SP.powerMap*powerScaleFactor
                
                # Make Wiener filter
                #filt=1.0/(SP.powerMap+NP.powerMap)
                #med=np.median(filt)
                #filt[np.where(filt>10*med)]=med
                #kernelSize=(5,5)
                #filt=ndimage.gaussian_filter(filt, kernelSize)
                #self.G=SP.powerMap*filt
                self.G=SP.powerMap/(SP.powerMap+NP.powerMap)
                self.G=np.nan_to_num(self.G)
                #astImages.saveFITS(self.filterFileName, fft.fftshift(self.G), None)

                # Apply the filter - note we apply weights here now like Toby
                # We also window the map in the same way as noise was windowed
                weightedMap=mapDict['data']*np.sqrt(mapDict['weights']/mapDict['weights'].max())
                apodlm=lm.createGaussianApodization(pad=20, kern=10)
                apodlm=apodlm.data
                weightedMap=weightedMap*apodlm
                fMap=fftTools.fftFromLiteMap(liteMap.liteMapFromDataAndWCS(weightedMap, self.wcs))
                filteredMaps['%d' % int(mapDict['obsFreqGHz'])]=np.real(fft.ifft2(fMap.kMap[:,:]*self.G[:,:]))

                # Check and correct for bias - apply filter to the signal template - do we recover the same flux?
                peakCoords=np.where(abs(signalMap.data) == abs(signalMap.data).max())
                yPeak=peakCoords[0][0]
                xPeak=peakCoords[1][0]
                fSignalMap=fftTools.fftFromLiteMap(signalMap)
                filteredSignalMap=np.real(fft.ifft2(fSignalMap.kMap[:,:]*self.G[:,:]))                                     
                amplitudeIn=signalMap.data[yPeak, xPeak]
                amplitudeOut=filteredSignalMap[yPeak, xPeak].real
                filteredMaps['%d' % int(mapDict['obsFreqGHz'])]=filteredMaps['%d' % int(mapDict['obsFreqGHz'])]*(amplitudeIn/amplitudeOut)
                
            # Linearly combine filtered maps and convert to yc if asked to do so
            if 'mapCombination' in self.params.keys():
                combinedMap=np.zeros(filteredMaps[filteredMaps.keys()[0]].shape)
                for key in filteredMaps.keys():
                    combinedMap=combinedMap+self.params['mapCombination'][key]*filteredMaps[key]
                combinedMap=mapTools.convertToY(combinedMap, self.params['mapCombination']['rootFreqGHz'])
                combinedObsFreqGHz='yc'
                #print "So - did the map combination make sense?"
                #ipshell()
                #sys.exit()
            else:
                # If no linear map combination given, assume we only want the first item in the filtered maps list
                combinedObsFreqGHz=self.unfilteredMapsDictList[0]['obsFreqGHz']
                combinedMap=filteredMaps['%d' % combinedObsFreqGHz]
            
            # Now that we're applying the weights in here, we may as well write out the S/N map here too
            # So we measure the RMS of the whole filtered map
            goodAreaMask=np.greater_equal(apodlm, 1.0) # don't want the apodized edges of the map to bias this
            mapMean=np.mean(combinedMap[goodAreaMask])
            mapRMS=np.std(combinedMap[goodAreaMask])
            sigmaClip=3.0
            for i in range(10):
                mask=np.less(abs(combinedMap), abs(mapMean+sigmaClip*mapRMS))
                mask=np.logical_and(goodAreaMask, mask)
                mapMean=np.mean(combinedMap[mask])
                mapRMS=np.std(combinedMap[mask])
            SNMap=combinedMap/mapRMS
            
            # Save the filtered map at this iteration so we can see if we're making a difference
            #if self.diagnosticsDir != None:
                #outPath=self.diagnosticsDir+os.path.sep+"iteration_%d_%s.fits" % (i, self.label)
                #astImages.saveFITS(outPath, combinedMap, self.wcs)
                
            # Now find clusters so we can mask them on the next iteration
            # Also save DS9 .reg file here so we can see if there are differences in what we detect
            #sigmaImage=mapTools.makeSigmaMap(combinedMap, mapDict['weights'])
            #SNMap=mapTools.makeSNMap(combinedMap, sigmaImage)
            #imageDict={'iterFilter': {}}
            #imageDict['iterFilter']['SNMap']=SNMap
            #imageDict['iterFilter']['wcs']=self.wcs
            #photometry.findObjects(imageDict, SNMap = 'array', 
                                    #threshold = self.params['detParams']['thresholdSigma'],
                                    #minObjPix = self.params['detParams']['minObjPix'], 
                                    #rejectBorder = 10, makeDS9Regions = False, 
                                    #writeSegmentationMap = False)            
            #clusterCatalog=imageDict['iterFilter']['catalog']
            #if self.diagnosticsDir != None:
                #outRegFileName=self.diagnosticsDir+os.path.sep+"iteration_%d_%s.reg" % (i, self.label)
                #catalogTools.catalog2DS9(clusterCatalog, outRegFileName)
            
            # Save filter profile in real space
            self.saveRealSpaceFilterProfile()
            
            # Add stuff to header for doing photometry
            # Hacked here - flux inside 2' radius from signal template, needs to be in yc
            # Big assumption here - need to think how this works for e.g. beta profile etc.
            # Maybe we can multiply the map through by something so that peak value = 2' value?
            ycSignalData=mapTools.convertToY(signalMap.data, obsFrequencyGHz = 148)
            ra0, dec0=self.wcs.pix2wcs(xPeak, yPeak)
            ra1, dec1=self.wcs.pix2wcs(xPeak+1, yPeak+1)    
            xLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
            yLocalDegPerPix=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
            xPix=np.array([np.arange(0, ycSignalData.shape[1], dtype=float)]*ycSignalData.shape[0])-xPeak
            yPix=(np.array([np.arange(0, ycSignalData.shape[0], dtype=float)]*ycSignalData.shape[1])-yPeak).transpose()
            xDeg=xPix*xLocalDegPerPix
            yDeg=yPix*yLocalDegPerPix
            rDegMap=np.sqrt(xDeg**2+yDeg**2) 
            apertureMask=np.less(rDegMap, 2.0/60.0)  # 2' radius
            signalFluxInAperture=np.sum(ycSignalData[apertureMask])
            arcmin2PerPix=xLocalDegPerPix*yLocalDegPerPix*60.0**2
            signalFluxInAperture=signalFluxInAperture*arcmin2PerPix # ?
            self.wcs.header.update('SAPFLUX', signalFluxInAperture)
            self.wcs.header.update('SPEAK', ycSignalData[yPeak, xPeak])
            self.wcs.header.update('SRADIUS', 2.0/60.0)
            self.wcs.updateFromHeader()
            #if 'gauss' not in self.label:
                #print "So, can we recover the input yc from the peak in the filtered map?"
                #ipshell()
                #sys.exit()
            
        return {'data': combinedMap, 'wcs': self.wcs, 'obsFreqGHz': combinedObsFreqGHz,
                'SNMap': SNMap}
                
#------------------------------------------------------------------------------------------------------------
class MatchedFilter(MapFilter):
    """Yes, it doesn't make that much sense deriving this from Wiener filter, but this is a quick and
    dirty test for now and I want the noise stuff.
    
    """

    def buildAndApply(self):
        
        print ">>> Building filter %s ..." % (self.label)

        self.makeRadiansMap()

        if 'iterations' not in self.params.keys():
            self.params['iterations']=1
        iterations=self.params['iterations']
        clusterCatalog=[]
        for i in range(iterations):
            
            #print "... iteration %d ..." % (i) 
            
            # Filter both maps with the same scale filter
            filteredMaps={}
            for mapDict in self.unfilteredMapsDictList:   
                
                # Iterative filtering stuff currently disabled
                # Mask out clusters so they don't go into the noise component of the filter
                #maskedDict=mapTools.maskOutSources(mapDict['data'], self.wcs, clusterCatalog, 
                                                   #radiusArcmin = self.params['maskingRadiusArcmin'], 
                                                   #mask = 'whiteNoise')
                #maskedData=maskedDict['data']
                #if self.diagnosticsDir != None:
                    #outPath=self.diagnosticsDir+os.path.sep+"iteration_%d_%s_mask.fits" % (i, self.label)
                    #astImages.saveFITS(outPath, maskedDict['mask'], self.wcs)
                
                maskedData=mapDict['data']
                
                # FFT - note using flipper 0.1.3 Gaussian apodization
                # The apodization here is only used if we're estimating noise power directly from the map itself
                # otherwise, the FFT here is just so we have the modThetaMap etc. which we would use if
                # using e.g. beta model templates where we're setting our own normalisation, rather than
                # using e.g. the profile models.
                lm=liteMap.liteMapFromDataAndWCS(maskedData, self.wcs)                               
                lm.data=lm.data*np.sqrt(mapDict['weights']/mapDict['weights'].max()) # this appears to make no difference
                apodlm=lm.createGaussianApodization(pad=10, kern=5)               
                lm.data=lm.data*apodlm.data
                fMaskedMap=fftTools.fftFromLiteMap(lm)
                fMaskedMap.modThetaMap=180.0/(fMaskedMap.modLMap+1)       
                                
                # Make noise power spectrum
                if self.params['noiseParams']['method'] == 'dataMap':
                    NP=fftTools.powerFromFFT(fMaskedMap)
                elif self.params['noiseParams']['method'] == '4WayMaps':
                    NP=self.noisePowerFrom4WayMaps(mapDict['weights'], mapDict['obsFreqGHz'])
                else:
                    raise Exception, "noise method must be either 'dataMap' or '4WayMaps'"
                
                # FFT of signal
                signalMapDict=self.makeSignalTemplateMap(mapDict['beamFileName'], mapDict['obsFreqGHz'])
                signalMap=signalMapDict['signalMap']
                fftSignal=fftTools.fftFromLiteMap(signalMap)
                
                # Toby style -note smoothing noise is essential!
                filt=1.0/NP.powerMap      
                med=np.median(filt)
                filt[np.where(filt>10*med)]=med
                kernelSize=(5,5)
                filt=ndimage.gaussian_filter(filt, kernelSize)
                filt=filt*abs(fftSignal.kMap)
                cov=filt*np.abs(fftSignal.kMap)**2
                integral=cov.sum()/signalMap.Nx/signalMap.Ny
                filt=filt/integral
                self.G=filt
                
                # NOTE: Do not disable this weighting
                weightedMap=mapDict['data']*np.sqrt(mapDict['weights']/mapDict['weights'].max())
                apodlm=lm.createGaussianApodization(pad=10, kern=5)
                weightedMap=weightedMap*apodlm.data
                fMap=fftTools.fftFromLiteMap(liteMap.liteMapFromDataAndWCS(weightedMap, self.wcs))
                
                # Check and correct for bias - apply filter to the signal template - do we recover the same flux?
                # Rescale the filtered map so that this is taken care of (better way than this, but this
                # should be pretty good)
                peakCoords=np.where(abs(signalMap.data) == abs(signalMap.data).max())
                yPeak=peakCoords[0][0]
                xPeak=peakCoords[1][0]
                fSignalMap=fftTools.fftFromLiteMap(signalMap)
                filteredSignalMap=np.real(fft.ifft2(fSignalMap.kMap*self.G))                          
                
                #---
                # Use the signal map we made using MatchedFilter to figure out how much it has been rolled off by:
                # 1. The high pass filter (bck sub step)
                # 2. The matched filter itself (includes beam)
                signalProperties=signalMapDict['inputSignalProperties']
                if self.params['outputUnits'] == 'yc':
                    # Normalise such that peak value in filtered map == y0, taking out the effect of the beam
                    filteredSignalMap=mapTools.convertToY(filteredSignalMap, obsFrequencyGHz = signalProperties['obsFreqGHz'])
                    signalNorm=signalProperties['y0']/filteredSignalMap.max()
                elif self.params['outputUnits'] == 'Y500':
                    # Normalise such that peak value in filtered map == Y500 for the input SZ cluster model
                    # We can get this from yc and the info in signalProperties, so we don't really need this
                    print "implement signal norm for %s" % (self.params['outputUnits'])
                    IPython.embed()
                    sys.exit()
                elif self.params['outputUnits'] == 'uK':
                    # Normalise such that peak value in filtered map == peak value of source in uK
                    signalNorm=1.0/filteredSignalMap.max()
                elif self.params['outputUnits'] == 'Jy/beam':
                    # Normalise such that peak value in filtered map == flux density of the source in Jy/beam
                    print "implement signal norm for %s" % (self.params['outputUnits'])
                    IPython.embed()
                    sys.exit()
                else:
                    raise Exception, "didn't understand 'outputUnits' given in the .par file"

                #---
                # Filter actual map
                filteredMaps['%d' % int(mapDict['obsFreqGHz'])]=np.real(fft.ifft2(self.G*fMap.kMap))*signalNorm
                
                # Apply the filter to the noiseless signal only sim - for checking y recovery later
                # NOTE: This will want fixing up to work in multi-frequency mode...
                #lmSim=liteMap.liteMapFromDataAndWCS(mapDict['simData'], self.wcs)                               
                #lmSim.data=lmSim.data*apodlm.data
                #fSimMap=fftTools.fftFromLiteMap(lmSim)
                #filteredSimMap=np.real(fft.ifft2(self.G*fSimMap.kMap*filterNormFactor))
                #filteredSimMap=mapTools.convertToY(filteredSimMap, self.params['mapCombination']['rootFreqGHz'])  
                # Or... add the signal only sim to the real map
                if 'simData' in mapDict.keys():
                    lmSim=liteMap.liteMapFromDataAndWCS(mapDict['simData']+mapDict['data'], self.wcs)                               
                    lmSim.data=lmSim.data*apodlm.data
                    fSimMap=fftTools.fftFromLiteMap(lmSim)
                    filteredSimMap=np.real(fft.ifft2(self.G*fSimMap.kMap))*signalNorm
                    filteredSimMap=mapTools.convertToY(filteredSimMap, self.params['mapCombination']['rootFreqGHz'])  
                else:
                    filteredSimMap=None
                
            # Linearly combine filtered maps and convert to yc if asked to do so
            if 'mapCombination' in self.params.keys():
                combinedMap=np.zeros(filteredMaps[filteredMaps.keys()[0]].shape)
                for key in filteredMaps.keys():
                    combinedMap=combinedMap+self.params['mapCombination'][key]*filteredMaps[key]
            else:
                # If no linear map combination given, assume we only want the first item in the filtered maps list
                combinedObsFreqGHz=self.unfilteredMapsDictList[0]['obsFreqGHz']
                combinedMap=filteredMaps['%d' % combinedObsFreqGHz]

            # Convert to whatever output units we want:
            # Jy/sr (should go to Jy/beam eventually) for sources
            # yc for SZ clusters
            if 'outputUnits' in self.params.keys():
                if self.params['outputUnits'] == 'yc':
                    combinedMap=mapTools.convertToY(combinedMap, self.params['mapCombination']['rootFreqGHz'])
                    combinedObsFreqGHz='yc'
                    mapUnits='yc'
                elif self.params['outputUnits'] == 'uK':
                    combinedObsFreqGHz=self.params['mapCombination']['rootFreqGHz']
                    mapUnits='uK'
                elif self.params['outputUnits'] == 'Jy/beam':
                    print "Jy/beam here"
                    combinedObsFreqGHz=self.params['mapCombination']['rootFreqGHz']
                    mapUnits='Jy/beam'
                    IPython.embed()
                    sys.exit()
                else:
                    raise Exception, 'need to specify "outputUnits" ("yc", "uK", or "Jy/beam") in filter params'
                
            # We'll be saving this shortly...
            apodMask=np.zeros(apodlm.data.shape)
            apodMask[np.greater(apodlm.data, 0.999999)]=1.0 # not sure why == 1 doesn't work
            
            # Now that we're applying the weights in here, we may as well write out the S/N map here too
            # So we measure the RMS of the whole filtered map
            goodAreaMask=np.greater_equal(apodMask, 1.0) # don't want the apodized edges of the map to bias this
            mapMean=np.mean(combinedMap[goodAreaMask])
            mapRMS=np.std(combinedMap[goodAreaMask])
            sigmaClip=3.0
            for i in range(10):
                mask=np.less(abs(combinedMap), abs(mapMean+sigmaClip*mapRMS))
                mask=np.logical_and(goodAreaMask, mask)
                mapMean=np.mean(combinedMap[mask])
                mapRMS=np.std(combinedMap[mask])
            SNMap=combinedMap/mapRMS
            
            # NOTE: Applying weights in this way disabled Feb 2014
            # Undo the weighting in the yc map, otherwise we underestimate decrement sizes!
            #combinedMap=combinedMap/np.sqrt(mapDict['weights']/mapDict['weights'].max())
            #combinedMap=np.nan_to_num(combinedMap)
            
            maskFileName=self.diagnosticsDir+os.path.sep+"areaMask.fits"
            if os.path.exists(maskFileName) == False:
                astImages.saveFITS(maskFileName, apodMask, mapDict['wcs'])
                    
            # Blank out the apodised area
            combinedMap=combinedMap*apodMask
            SNMap=SNMap*apodMask
                                        
            # Save filter profile in real space
            self.saveRealSpaceFilterProfile()        
            
        return {'data': combinedMap, 'simData': filteredSimMap, 'wcs': self.wcs, 
                'obsFreqGHz': combinedObsFreqGHz, 'SNMap': SNMap, 'signalMap': signalMap.data, 
                'mapUnits': mapUnits, 'inputSignalProperties': signalMapDict['inputSignalProperties']}
            
#------------------------------------------------------------------------------------------------------------
class RealSpaceMatchedFilter(MapFilter):
    """Makes a matched-filter kernel using the noise properties of a specified region of the map (e.g.
    the deepest part) in Fourier space, which is converted into real space and truncated such that the 
    kernel is small enough in footprint to be applied by directly convolving with the map in real space 
    in a short amount of time. 
    
    First though we high-pass filter the maps on a similar scale using a Mexican hat (difference of Gaussian
    smoothed images), to get rid of large scale noise from CMB and atmosphere.
    
    Tne S/N map is constructed, measuring locally on a scale 3 x that of the filter kernel, with 3-sigma
    clipping applied. This allows us to take into account how the noise level varies across the map.
    
    A rank filter is used to zap noisy regions at the edge of the map, where the RMS values are
    not accurate.
    
    Optionally, a 'surveyMask' can be applied (e.g., for point source masking).
    
    A map of the final area searched for clusters called 'areaMask.fits' is written in the diagnostics/ 
    folder.
        
    """

    def buildKernel(self, mapDict, RADecSection, RADeg = 'centre', decDeg = 'centre'):
        """Builds the real space kernel itself. 
        
        RADeg, decDeg are used for figuring out pixel scales for background subtraction
        
        Returns kern2d, signalNorm, bckSubScaleArcmin, signalProperties
        
        """
        
        wcs=mapDict['wcs']
            
        # Build the matched-filter kernel in a small section of the map
        # Apply the same difference of Gaussians high pass filter here
        # NOTE: we could merge 'bckSubScaleArcmin' and 'maxArcmin' keys here!
        kernelMaxArcmin=self.params['noiseParams']['kernelMaxArcmin']
        #mapDict['bckSubScaleArcmin']=maxArcmin
        keysWanted=['mapFileName', 'weightsFileName', 'obsFreqGHz', 'units', 'beamFileName', 'addNoise', 
                    'pointSourceRemoval']
        kernelUnfilteredMapsDict={}
        for k in keysWanted:
            if k in mapDict.keys():
                kernelUnfilteredMapsDict[k]=mapDict[k]
        kernelUnfilteredMapsDict['RADecSection']=RADecSection
        kernelUnfilteredMapsDictList=[kernelUnfilteredMapsDict]
        kernelLabel="realSpaceKernel_%s" % (self.label)
        matchedFilterDir=self.diagnosticsDir+os.path.sep+kernelLabel
        if os.path.exists(matchedFilterDir) == False:
            os.makedirs(matchedFilterDir)
        if os.path.exists(matchedFilterDir+os.path.sep+'diagnostics') == False:
            os.makedirs(matchedFilterDir+os.path.sep+'diagnostics')
        matchedFilterClass=eval(self.params['noiseParams']['matchedFilterClass'])
        matchedFilter=matchedFilterClass(kernelLabel, kernelUnfilteredMapsDictList, self.params, 
                                            outDir = matchedFilterDir, 
                                            diagnosticsDir = matchedFilterDir+os.path.sep+'diagnostics')
        filteredMapDict=matchedFilter.buildAndApply()
                
        # Turn the matched filter into a smaller real space convolution kernel
        # This means we have to roll off the kernel to 0 at some radius
        # This is set by maxArcmin in the .par file
        prof, arcminRange=matchedFilter.makeRealSpaceFilterProfile()
        rIndex=np.where(arcminRange > kernelMaxArcmin)[0][0]
        # Alternatively, roll off to zero after the second zero crossing
        # NOTE: now setting in the .par file, uncomment below to switch back
        #segProf=ndimage.label(np.greater(prof, 0))[0]
        #rIndex=np.where(segProf == 2)[0][0]
        #maxArcmin=arcminRange[rIndex]

        # Make 2d kernel
        mask=np.less(arcminRange, kernelMaxArcmin)
        rRadians=np.radians(arcminRange/60.)
        r2p=interpolate.interp1d(rRadians[mask], prof[mask], bounds_error=False, fill_value=0.0)
        profile2d=r2p(matchedFilter.radiansMap)
        y, x=np.where(profile2d == profile2d.max())
        y=y[0]
        x=x[0]
        yMin=y-rIndex
        yMax=y+rIndex
        xMin=x-rIndex
        xMax=x+rIndex
        if (yMax-yMin) % 2 == 0:
            yMin=yMin+1
        if (xMax-xMin) % 2 == 0:
            xMin=xMin+1
        kern2d=profile2d[yMin:yMax, xMin:xMax]
        kern2dRadiansMap=matchedFilter.radiansMap[yMin:yMax, xMin:xMax]

        # This is what to high pass filter on
        bckSubScaleArcmin=arcminRange[prof == prof.min()][0]
        
        # Use the signal map we made using MatchedFilter to figure out how much it has been rolled off by:
        # 1. The high pass filter (bck sub step)
        # 2. The matched filter itself (includes beam)
        signalMap=filteredMapDict['signalMap']      # Note that this has had the beam applied already
        signalProperties=filteredMapDict['inputSignalProperties']
        # Should add an applyKernel function to do all this
        if self.params['bckSub'] == True:
            filteredSignal=mapTools.subtractBackground(signalMap, wcs, RADeg = RADeg, decDeg = decDeg,
                                                       smoothScaleDeg = bckSubScaleArcmin/60.0)
        else:
            filteredSignal=np.zeros(signalMap.shape)+signalMap
        filteredSignal=ndimage.convolve(filteredSignal, kern2d) 
        if self.params['outputUnits'] == 'yc':
            # Normalise such that peak value in filtered map == y0, taking out the effect of the beam
            filteredSignal=mapTools.convertToY(filteredSignal, obsFrequencyGHz = signalProperties['obsFreqGHz'])
            signalNorm=signalProperties['y0']/filteredSignal.max()
        elif self.params['outputUnits'] == 'Y500':
            # Normalise such that peak value in filtered map == Y500 for the input SZ cluster model
            # We can get this from yc and the info in signalProperties, so we don't really need this
            print "implement signal norm for %s" % (self.params['outputUnits'])
            IPython.embed()
            sys.exit()
        elif self.params['outputUnits'] == 'uK':
            # Normalise such that peak value in filtered map == peak value of source in uK
            # We take out the effect of the pixel window function here - this assumes we're working with sources
            # NOTE: for amplitude, the mappers want the amplitude of the delta function, not the beam
            signalNorm=1.0/(signalProperties['pixWindowFactor']*filteredSignal.max())
        elif self.params['outputUnits'] == 'Jy/beam':
            # Normalise such that peak value in filtered map == flux density of the source in Jy/beam
            print "implement signal norm for %s" % (self.params['outputUnits'])
            IPython.embed()
            sys.exit()
        else:
            raise Exception, "didn't understand 'outputUnits' given in the .par file"

        # Save 2d kernel - we need this (at least for the photometry ref scale) to calc Q later
        # Add bckSubScaleArcmin to the header
        kernWCS=wcs.copy()
        if self.params['bckSub'] == True:
            kernWCS.header['BCKSCALE']=bckSubScaleArcmin
        kernWCS.header['SIGNORM']=signalNorm
        RADecLabel=str(RADecSection).replace(",", "_").replace(" ", "").replace("[", "").replace("]", "").replace(".", "p")
        astImages.saveFITS(self.diagnosticsDir+os.path.sep+"kern2d_%s_%s.fits" % (self.label, RADecLabel), kern2d, kernWCS)
        
        # Filter profile plot   
        # Save the stuff we plot first, in case we want to make a plot with multiple filters on later
        np.savez(self.diagnosticsDir+os.path.sep+"filterProf1D_%s.npz" % (self.label), 
                    arcminRange = arcminRange, prof = prof, mask = mask, bckSubScaleArcmin = bckSubScaleArcmin)
        plotSettings.update_rcParams()
        #fontDict={'size': 18, 'family': 'serif'}
        plt.figure(figsize=(9,6.5))
        ax=plt.axes([0.13, 0.12, 0.86, 0.86])
        #plt.tick_params(axis='both', which='major', labelsize=15)
        #plt.tick_params(axis='both', which='minor', labelsize=15)
        tck=interpolate.splrep(arcminRange[mask], prof[mask])
        plotRange=np.linspace(0, arcminRange[mask].max(), 1000)
        #plt.plot(arcminRange[mask], prof[mask])
        plt.plot(plotRange, interpolate.splev(plotRange, tck), 'k-')
        plt.xlabel("$\\theta$ (arcmin)")
        plt.ylabel("Amplitude")
        #plt.title(self.label)
        #plt.plot(arcminRange[mask], [0]*len(arcminRange[mask]), 'k--')
        plt.xlim(0, arcminRange[mask].max())
        plt.plot([bckSubScaleArcmin]*3, np.linspace(-1, 1.2, 3), 'k--')
        plt.ylim(-0.2, 1.2)
        plt.savefig(self.diagnosticsDir+os.path.sep+"filterPlot1D_%s.pdf" % (self.label))
        plt.close()
        
        return kern2d, signalNorm, bckSubScaleArcmin, signalProperties
            
            
    def buildAndApply(self):
        
        print ">>> Building filter %s ..." % (self.label)
        
        filteredMaps={}
        for mapDict in self.unfilteredMapsDictList:   

            mapData=mapDict['data']
            wcs=mapDict['wcs']
            surveyMask=mapDict['surveyMask']
            psMask=mapDict['psMask']
            
            # Make kernels at different decs (can add RA needed later if necessary)
            # Copes with CAR distortion at large | dec |
            # NOTE: the way this is done currently means we should pick something contiguous in dec direction at fixed RA
            RAMin, RAMax, decMin, decMax=wcs.getImageMinMaxWCSCoords()
            if self.params['noiseParams']['RADecSection'][2] == 'numDecSteps':
                numDecSteps=float(self.params['noiseParams']['RADecSection'][3])
                decEdges=np.linspace(decMin, decMax, numDecSteps+1)
                RADecSectionDictList=[]
                for i in range(len(decEdges)-1):
                    RADecSectionDict={'RADecSection': [self.params['noiseParams']['RADecSection'][0],
                                                       self.params['noiseParams']['RADecSection'][1],
                                                       decEdges[i], decEdges[i+1]],
                                      'applyDecMin': decEdges[i], 'applyDecMax': decEdges[i+1],
                                      'applyRAMin': RAMin, 'applyRAMax': RAMax}
                    RADecSectionDictList.append(RADecSectionDict)
            else:
                RADecSectionDictList=[{'RADecSection': self.params['noiseParams']['RADecSection'],
                                       'applyDecMin': decMin, 'applyDecMax': decMax,
                                       'applyRAMin': RAMin, 'applyRAMax': RAMax}]
            
            # Building the filter in different regions (e.g., dec strips to handle CAR distortion)
            # NOTE: RADecSectionDictList is something we could pickle instead of saving individual kern2d files etc.?
            for RADecSectionDict in RADecSectionDictList:
                
                # We need this for the background subtraction x, y pixel scales
                applyDecCentre=(RADecSectionDict['applyDecMax']+RADecSectionDict['applyDecMin'])/2.
                applyRACentre=(RADecSectionDict['applyRAMax']+RADecSectionDict['applyRAMin'])/2.
                
                # For now we'll just restrict the y-range (dec) on which we do filtering... can add RA blocks later
                x, yMin=wcs.wcs2pix(applyRACentre, RADecSectionDict['applyDecMin'])
                yMin=int(round(yMin))
                x, yMax=wcs.wcs2pix(applyRACentre, RADecSectionDict['applyDecMax'])
                yMax=int(round(yMax))
                
                # Build the matched-filter kernel in a small section of the map
                kern2d, signalNorm, bckSubScaleArcmin, signalProperties=self.buildKernel(mapDict, RADecSectionDict['RADecSection'],
                                                                                         RADeg = applyRACentre, decDeg = applyDecCentre)
                RADecSectionDict['kern2d']=kern2d
                RADecSectionDict['signalNorm']=signalNorm
                RADecSectionDict['bckSubScaleArcmin']=bckSubScaleArcmin
                RADecSectionDict['signalProperties']=signalProperties
                RADecSectionDict['applyRACentre']=applyRACentre
                RADecSectionDict['applyDecCentre']=applyDecCentre
                RADecSectionDict['yMin']=yMin
                RADecSectionDict['yMax']=yMax
                if yMax+RADecSectionDict['kern2d'].shape[0] < mapData.shape[0]:
                    yOverlap=RADecSectionDict['kern2d'].shape[0]
                else:
                    yOverlap=0
                RADecSectionDict['yOverlap']=yOverlap
            
            # Apply the high pass filter - subtract background on larger scales using difference of Gaussians  
            for RADecSectionDict in RADecSectionDictList:
                yMin=RADecSectionDict['yMin']
                yMax=RADecSectionDict['yMax']
                yOverlap=RADecSectionDict['yOverlap']
                buff=np.zeros(mapData[yMax:yMax+yOverlap].shape)+mapData[yMax:yMax+yOverlap]
                if self.params['bckSub'] == True:
                    mapData[yMin:yMax+yOverlap, :]=mapTools.subtractBackground(mapData[yMin:yMax+yOverlap, :], wcs, 
                                                                      RADeg = RADecSectionDict['applyRACentre'], 
                                                                      decDeg = RADecSectionDict['applyDecCentre'],
                                                                      smoothScaleDeg = RADecSectionDict['bckSubScaleArcmin']/60.)
                    mapData[yMax:yMax+yOverlap]=buff
            if 'saveHighPassMap' in self.params['noiseParams'] and self.params['noiseParams']['saveHighPassMap'] == True:
                bckSubFileName=self.diagnosticsDir+os.path.sep+"bckSub_%s.fits" % (self.label)
                #astImages.saveFITS(bckSubFileName, bckSubData, mapDict['wcs'])
                astImages.saveFITS(bckSubFileName, mapData, mapDict['wcs'])
            
            # Apply the kernel
            backupData=np.zeros(mapData.shape)+mapData
            for RADecSectionDict in RADecSectionDictList:
                yMin=RADecSectionDict['yMin']
                yMax=RADecSectionDict['yMax']
                yOverlap=RADecSectionDict['yOverlap']
                buff=np.zeros(mapData[yMax:yMax+yOverlap].shape)+mapData[yMax:yMax+yOverlap]
                t0=time.time()
                print "... convolving map with kernel [%d:%d] ..." % (yMin, yMax)
                mapData[yMin:yMax+yOverlap, :]=ndimage.convolve(mapData[yMin:yMax+yOverlap, :], RADecSectionDict['kern2d'])   
                filtBuff=np.zeros(mapData[yMax:yMax+yOverlap].shape)+mapData[yMax:yMax+yOverlap]
                mapData[yMax:yMax+yOverlap]=buff
                t1=time.time()
                print "... took %.3f sec ..." % (t1-t0)
                # Apply the normalisation
                mapData[yMin:yMax, :]=mapData[yMin:yMax, :]*RADecSectionDict['signalNorm']
            
            #filteredMaps['%d' % int(mapDict['obsFreqGHz'])]=filteredMap
            filteredMaps['%d' % int(mapDict['obsFreqGHz'])]=mapData
        
        # Linearly combine filtered maps and convert to yc if asked to do so
        # NOTE: Here, we ARE converting to yc always... we can change this somewhere for point source folks (use Jy/beam)
        if 'mapCombination' in self.params.keys():
            combinedMap=np.zeros(filteredMaps[filteredMaps.keys()[0]].shape)
            for key in filteredMaps.keys():
                combinedMap=combinedMap+self.params['mapCombination'][key]*filteredMaps[key]
        else:
            # If no linear map combination given, assume we only want the first item in the filtered maps list
            combinedObsFreqGHz=self.unfilteredMapsDictList[0]['obsFreqGHz']
            combinedMap=filteredMaps['%d' % combinedObsFreqGHz]
        
        # We don't need the filteredMaps dict from here, so zap it to save memory
        del filteredMaps
        
        # Convert to whatever output units we want:
        # Jy/sr (should go to Jy/beam eventually) for sources
        # yc for SZ clusters
        if 'outputUnits' in self.params.keys():
            if self.params['outputUnits'] == 'yc':
                combinedMap=mapTools.convertToY(combinedMap, self.params['mapCombination']['rootFreqGHz'])
                combinedObsFreqGHz='yc'
                mapUnits='yc'
            elif self.params['outputUnits'] == 'uK':
                combinedObsFreqGHz=self.params['mapCombination']['rootFreqGHz']
                mapUnits='uK'
            elif self.params['outputUnits'] == 'Jy/beam':
                print "Jy/beam here"
                combinedObsFreqGHz=self.params['mapCombination']['rootFreqGHz']
                mapUnits='Jy/beam'
                IPython.embed()
                sys.exit()
            else:
                raise Exception, 'need to specify "outputUnits" ("yc", "uK", or "Jy/beam") in filter params'

        # Cython-based SNMap calc, which uses an annulus mask, but doesn't do clipping
        # Takes ~460 sec on equD56
        #t0=time.time()
        #annulus=photometry.makeAnnulus(rIndex, rIndex*2)
        #SNMap=nemoCython.makeLocalSNMap(combinedMap, annulus)
        #t1=time.time()
        
        # Make SN map by local RMS measurements on a grid, over the whole filtered map
        # We're not making a separate noise map here any more to save memory
        print "... making SN map ..."
        gridSize=int(round((self.params['noiseParams']['noiseGridArcmin']/60.)/wcs.getPixelSizeDeg()))
        #gridSize=rIndex*3
        overlapPix=gridSize/2
        numXChunks=combinedMap.shape[1]/gridSize
        numYChunks=combinedMap.shape[0]/gridSize
        yChunks=np.linspace(0, combinedMap.shape[0], numYChunks+1, dtype = int)
        xChunks=np.linspace(0, combinedMap.shape[1], numXChunks+1, dtype = int)
        SNMap=np.zeros(combinedMap.shape)
        apodMask=np.not_equal(combinedMap, 0)
        # We could make below behaviour default if match photFilter? Would need to see photFilter though...
        if 'saveRMSMap' in self.params['noiseParams'] and self.params['noiseParams']['saveRMSMap'] == True:
            RMSMap=np.zeros(combinedMap.shape)
        t0=time.time()
        for i in range(len(yChunks)-1):
            for k in range(len(xChunks)-1):
                y0=yChunks[i]-overlapPix
                y1=yChunks[i+1]+overlapPix
                x0=xChunks[k]-overlapPix
                x1=xChunks[k+1]+overlapPix
                if y0 < 0:
                    y0=0
                if y1 > combinedMap.shape[0]:
                    y1=combinedMap.shape[0]
                if x0 < 0:
                    x0=0
                if x1 > combinedMap.shape[1]:
                    x1=combinedMap.shape[1]
                chunkValues=combinedMap[y0:y1, x0:x1]
                # 3-sigma clipped stdev - 12 sec
                # Tried biweight scale version, 10 x slower than this
                if np.not_equal(chunkValues, 0).sum() != 0:
                    goodAreaMask=np.greater_equal(apodMask[y0:y1, x0:x1], 1.0)
                    chunkMean=np.mean(chunkValues[goodAreaMask])
                    chunkRMS=np.std(chunkValues[goodAreaMask])
                    sigmaClip=3.0
                    for c in range(10):
                        mask=np.less(abs(chunkValues), abs(chunkMean+sigmaClip*chunkRMS))
                        mask=np.logical_and(goodAreaMask, mask)
                        chunkMean=np.mean(chunkValues[mask])
                        chunkRMS=np.std(chunkValues[mask])
                else:
                    chunkRMS=0.
                if chunkRMS > 0:
                    SNMap[y0:y1, x0:x1]=combinedMap[y0:y1, x0:x1]/chunkRMS
                    if 'saveRMSMap' in self.params['noiseParams'] and self.params['noiseParams']['saveRMSMap'] == True:
                        RMSMap[y0:y1, x0:x1]=chunkRMS
        t1=time.time()
        if 'saveRMSMap' in self.params['noiseParams'] and self.params['noiseParams']['saveRMSMap'] == True:
            RMSFileName=self.diagnosticsDir+os.path.sep+"RMSMap_%s.fits" % (self.label)
            astImages.saveFITS(RMSFileName, RMSMap, mapDict['wcs'])
        #print "... took %.3f sec ..." % (t1-t0)
        
        # Below is global RMS, for comparison
        #apodMask=np.not_equal(mapData, 0)
        #goodAreaMask=np.greater_equal(apodMask, 1.0) # don't want the apodized edges of the map to bias this
        #mapMean=np.mean(combinedMap[goodAreaMask])
        #mapRMS=np.std(combinedMap[goodAreaMask])
        #sigmaClip=3.0
        #for i in range(10):
            #mask=np.less(abs(combinedMap), abs(mapMean+sigmaClip*mapRMS))
            #mask=np.logical_and(goodAreaMask, mask)
            #mapMean=np.mean(combinedMap[mask])
            #mapRMS=np.std(combinedMap[mask])

        # Use rank filter to zap edges where RMS will be artificially low - we use a bit of a buffer here
        # Fold point source mask into survey mask here
        edgeCheck=ndimage.rank_filter(abs(mapData), 0, size = (int(round(gridSize*5.0)), int(round(gridSize*5.0))))
        edgeCheck=np.array(np.greater(edgeCheck, 0), dtype = float)
        combinedMap=combinedMap*edgeCheck
        apodMask=np.not_equal(combinedMap, 0)
        surveyMask=edgeCheck*surveyMask*psMask
        del edgeCheck

        # Apply final survey mask to signal-to-noise map
        # NOTE: need to avoid NaNs in here, otherwise map interpolation for e.g. S/N will fail later on
        SNMap=SNMap*surveyMask
        SNMap[np.isnan(SNMap)]=0.

        maskFileName=self.diagnosticsDir+os.path.sep+"areaMask.fits"
        if os.path.exists(maskFileName) == False:
            astImages.saveFITS(maskFileName, np.array(surveyMask, dtype = int), mapDict['wcs'])
                
        return {'data': combinedMap, 'simData': None, 'wcs': self.wcs, 'obsFreqGHz': combinedObsFreqGHz,
                'SNMap': SNMap, 'signalMap': kern2d, 'mapUnits': mapUnits,
                'inputSignalProperties': signalProperties}
            
#------------------------------------------------------------------------------------------------------------
class GaussianFilter(MapFilter):
    """Base class for filters using user specified Gaussian profile.
    
    """
    
    def makeSignalTemplateMap(self, beamFWHMArcmin, mapObsFreqGHz):
        """Makes a Gaussian model signal template map.
        
        Returns dictionary of {'signalMap', 'normInnerDeg':, 'normOuterDeg'}. The latter two keys are
        used to define an annulus in which the signal power is scaled to match the noise power.
        
        """
        
        # Setup Gaussian filter profile
        sigmaArcmin=beamFWHMArcmin/np.sqrt(8.0*np.log(2.0))            
        rArcmin=np.linspace(0.0, 30.0, 5000)
        profile1d=np.exp(-((rArcmin**2)/(2*sigmaArcmin**2)))
                    
        # Turn 1d profile into 2d
        rRadians=np.radians(rArcmin/60.0)
        r2p=interpolate.interp1d(rRadians, profile1d, bounds_error=False, fill_value=0.0)
        profile2d=r2p(self.radiansMap)
        signalMap=liteMap.liteMapFromDataAndWCS(profile2d, self.wcs)
        
        # Angular scales within which we will scale signal to match the noise
        rInnerDeg=(sigmaArcmin/60.0)*5.0
        rOuterDeg=(sigmaArcmin/60.0)*10.0
        
        # The ratio by which beam smoothing biases the intrinsic deltaT0
        beamDecrementBias=1.0#deltaT0/profile1d[0]  # assuming rDeg[0] is at 0 # 
        
        return {'signalMap': signalMap, 'normInnerDeg': rInnerDeg, 'normOuterDeg': rOuterDeg, 
                'powerScaleFactor': None, 'beamDecrementBias': beamDecrementBias}

#------------------------------------------------------------------------------------------------------------
class BeamFilter(MapFilter):
    """Base class for filters using beam profile files in Matthew + Kavi's format.
    
    NOTE: We haven't actually implemented beamDecrementBias and signalAreaSum in here.
    
    """
    
    def makeSignalTemplateMap(self, beamFileName, mapObsFreqGHz):
        """Makes a Gaussian model signal template map.
        
        Returns dictionary of {'signalMap', 'normInnerDeg':, 'normOuterDeg'}. The latter two keys are
        used to define an annulus in which the signal power is scaled to match the noise power.
        
        """
        
        signalMap, inputSignalProperties=simsTools.makeBeamModelSignalOnlyMap(obsFreqGHz, np.degrees(self.radiansMap),
                                                                              self.wcs, 
                                                                              beamFileName)

        return {'signalMap': signalMap, 'normInnerDeg': None, 'normOuterDeg': None, 
                'powerScaleFactor': None, 'beamDecrementBias': 1.0, 'signalAreaSum': 1.0,
                'inputSignalProperties': inputSignalProperties}
    
#------------------------------------------------------------------------------------------------------------
class BetaModelFilter(MapFilter):
    """Base class for filters using beta model profile.
    
    """
    
    def makeSignalTemplateMap(self, beamFWHMArcmin, mapObsFreqGHz):
        """Makes a beta model signal template map.
        
        Returns dictionary of {'signalMap', 'normInnerDeg':, 'normOuterDeg'}. The latter two keys are
        used to define an annulus in which the signal power is scaled to match the noise power.
        
        """
        
        # Setup Beta model profile
        rArcmin=np.linspace(0.0, 60.0, 5000)
        smoothArcmin=self.params['coreRadiusArcmin']
        profile1d=(1.0+(rArcmin/self.params['coreRadiusArcmin'])**2)**((1.0-3.0*self.params['beta'])/2)
        mask=np.greater(rArcmin, 5.0*self.params['coreRadiusArcmin'])
        profile1d[mask]*=np.exp(-(rArcmin[mask]-5.0*self.params['coreRadiusArcmin'])/5.0) # same as Matthew
        
        # Apply beam as Gaussian filter to profile
        beamSigma=beamFWHMArcmin/np.sqrt(8.0*np.log(2.0))            
        beamSigmaPix=beamSigma/(rArcmin[1]-rArcmin[0])
        profile1d=ndimage.gaussian_filter1d(profile1d, beamSigmaPix)

        # Truncate beyond 5 times core radius
        mask=np.greater(rArcmin, 5.0*self.params['coreRadiusArcmin'])
        profile1d[mask]=0.0
        
        # Turn 1d profile into 2d
        rRadians=np.radians(rArcmin/60.0)
        r2p=interpolate.interp1d(rRadians, profile1d, bounds_error=False, fill_value=0.0)
        profile2d=r2p(self.radiansMap)
        signalMap=liteMap.liteMapFromDataAndWCS(profile2d, self.wcs)        
        
        # Angular scales within which we will scale signal to match the noise
        rInnerDeg=(5.0*self.params['coreRadiusArcmin']/60.0)
        rOuterDeg=(rInnerDeg+5.0/60.0)
                
        return {'signalMap': signalMap, 'normInnerDeg': rInnerDeg, 'normOuterDeg': rOuterDeg,
                'powerScaleFactor': None}

#------------------------------------------------------------------------------------------------------------
class ArnaudModelFilter(MapFilter):
    """Base class for filters using the GNFW profile as described in Arnaud et al. (2010).
    
    """
    
    def makeSignalTemplateMap(self, beamFileName, mapObsFreqGHz):
        """Makes a model signal template map.
        
        Returns dictionary of {'signalMap', 'inputSignalProperties'}
        
        """
        
        signalMap, modelDict=simsTools.makeArnaudModelSignalMap(self.params['z'], self.params['M500MSun'], 
                                                                mapObsFreqGHz, np.degrees(self.radiansMap),
                                                                self.wcs, beamFileName)
        
        return {'signalMap': signalMap, 'inputSignalProperties': modelDict}
        
#------------------------------------------------------------------------------------------------------------
class ProfileFilter(MapFilter):
    """Base class for filters using arbitrary signal profiles stored in text files as per Ryan's code for
    generating templates based on Bode et al. models.
    
    """
    
    def makeSignalTemplateMap(self, beamFWHMArcmin, mapObsFreqGHz):
        """Makes a signal map from file containing given temperature profile with columns arcmin, profile.
        Note the units of profile are Kelvin
        In this case (unlike beta model or Gaussian cases), we use the model y_c to set the normalization
        of the filter.
        
        Returns dictionary of {'signalMap', 'normInnerDeg':, 'normOuterDeg'}. The latter two keys are
        used to define an annulus in which the signal power is scaled to match the noise power.
        
        """
        
        # Load in temperature profile template in format r (arcmin), delta T CMB (Kelvin)
        inFile=file(self.params['profileFile'], "r")
        lines=inFile.readlines()
        inFile.close()
        rArcmin=[]
        profile=[]
        for line in lines:
            if line[0] != "#" and len(line) > 3:
                bits=line.split()
                rArcmin.append(np.float64(bits[0]))
                profile.append(np.float64(bits[1]))
        rRadians=np.radians(np.array(rArcmin, dtype=np.float64)/60.0)
        profile=np.array(profile, dtype=np.float64)*1e6   # to micro K
                
        # The profile might not be given with a constant linear bin size, so interpolate here as well
        r2prof=interpolate.interp1d(rRadians, profile, bounds_error=False, fill_value=0.0)
        rRadiansLinear=np.linspace(rRadians[0], rRadians[-1], 10000)
        profileLinear=r2prof(rRadiansLinear)
        
        # Apply beam as Gaussian filter to profile
        # this is probably not the best way to handle the beam params
        beamFWHMRad=np.radians(beamFWHMArcmin/60.0)
        beamSigmaRad=beamFWHMRad/np.sqrt(8.0*np.log(2.0))            
        beamSigmaPix=beamSigmaRad/(rRadiansLinear[1]-rRadiansLinear[0])
        profileLinear=ndimage.gaussian_filter1d(profileLinear, beamSigmaPix)
        
        # Convert from dT the profile template is made in to dT at the map frequency
        ycProfile=mapTools.convertToY(profileLinear, obsFrequencyGHz = self.params['profileObsFreqGHz'])
        profileLinear=mapTools.convertToDeltaT(ycProfile, obsFrequencyGHz = mapObsFreqGHz)
        
        # Log FWHM of profile as plot
        arcmin2prof=interpolate.interp1d(np.degrees(rRadiansLinear)*60.0, abs(profileLinear), bounds_error=False, fill_value=0.0)
        plotRangeArcmin=np.linspace(0, 10, 10000)
        plotProfile=arcmin2prof(plotRangeArcmin)/abs(profileLinear).max()
        plotProfile[0]=1.0
        diff=abs(plotProfile-0.5)
        FWHMIndex=diff.tolist().index(diff.min())
        FWHMArcmin=plotRangeArcmin[FWHMIndex]*2
        fig=plt.figure(num=2, figsize=(8,8))
        fig.canvas.set_window_title('Signal Profile in Real Space')
        plt.clf()
        plt.title("Signal Profile - %s - FWHM = %.3f arcmin" % (self.label, FWHMArcmin))
        plt.ylabel("Amplitude")
        plt.xlabel("$\\theta$ (arcmin)")
        plt.plot(plotRangeArcmin, plotProfile)
        plt.xlim(0, 10)
        plt.ylim(0, 1.05)
        plt.savefig(self.diagnosticsDir+os.path.sep+"SignalProfile1d_"+self.label+".png")

        # Turn 1d profile into 2d
        r2prof=interpolate.interp1d(rRadiansLinear, profileLinear, bounds_error=False, fill_value=0.0)
        profile2d=r2prof(self.radiansMap)
        profile2d[np.where((abs(self.radiansMap) == 0.0))]=profileLinear[0]
        signalMap=liteMap.liteMapFromDataAndWCS(profile2d, self.wcs)

        # Scale the signal map according to map area, so that when we later take the signal power spectrum,
        # it is normalized.
        height=signalMap.data.nonzero()[0].max()-signalMap.data.nonzero()[0].min()
        width=signalMap.data.nonzero()[1].max()-signalMap.data.nonzero()[1].min()
        signalArea=float(signalMap.data[np.where(signalMap.data != 0)].shape[0])
        mapArea=float(signalMap.Nx*signalMap.Ny)
        scaleFactor=mapArea/signalArea
        
        rInnerDeg=None
        rOuterDeg=None
                        
        return {'signalMap': signalMap, 'normInnerDeg': rInnerDeg, 'normOuterDeg': rOuterDeg, 
                'powerScaleFactor': scaleFactor}
        
#------------------------------------------------------------------------------------------------------------
# Definitions of actual filters that can be used
class BetaModelMatchedFilter(MatchedFilter, BetaModelFilter):
    pass
class ProfileMatchedFilter(MatchedFilter, ProfileFilter):
    pass
class GaussianMatchedFilter(MatchedFilter, GaussianFilter):
    pass
class ArnaudModelMatchedFilter(MatchedFilter, ArnaudModelFilter):
    pass
class BeamMatchedFilter(MatchedFilter, BeamFilter):
    pass

class ArnaudModelRealSpaceMatchedFilter(RealSpaceMatchedFilter, ArnaudModelFilter):
    pass
class BeamRealSpaceMatchedFilter(RealSpaceMatchedFilter, BeamFilter):
    pass

class BetaModelWienerFilter(WienerFilter, BetaModelFilter):
    pass
class ProfileWienerFilter(WienerFilter, ProfileFilter):
    pass
class GaussianWienerFilter(WienerFilter, GaussianFilter):
    pass
class ArnaudModelWienerFilter(WienerFilter, ArnaudModelFilter):
    pass

#------------------------------------------------------------------------------------------------------------
def combinations(iterable, r):
    """This is in itertools in python 2.6.
    
    """
    
    # combinations('ABCD', 2) --> AB AC AD BC BD CD
    # combinations(range(4), 3) --> 012 013 023 123
    pool = tuple(iterable)
    n = len(pool)
    if r > n:
        return
    indices = range(r)
    yield tuple(pool[i] for i in indices)
    while True:
        for i in reversed(range(r)):
            if indices[i] != i + n - r:
                break
        else:
            return
        indices[i] += 1
        for j in range(i+1, r):
            indices[j] = indices[j-1] + 1
        yield tuple(pool[i] for i in indices)
        
