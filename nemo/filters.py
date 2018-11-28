"""

Filter classes are defined in this module, together with the `filterMaps` function that uses them.

There are two main classes of filter: `MatchedFilter` and `RealSpaceMatchedFilter`. 
`RealSpaceMatcherFilter` is the preferred one to use.

New base classes can be derived from the overall base class: `MapFilter`.

There are also base classes corresponding to filters with different signal templates: 
e.g., `BeamFilter`, `ArnaudModelFilter`.

The actual filters that can be used are derived from these, e.g.:

* `BeamMatchedFilter`
* `ArnaudModelMatchedFilter`
* `BeamRealSpaceMatchedFilter`
* `ArnaudModelRealSpaceMatchedFilter`

etc.

"""

import math
from pixell import enmap
from pixell import fft as enfft
from pixell import powspec
import astropy.wcs as enwcs
from astLib import *
import numpy as np
from numpy import fft
import pylab as plt
import os
from scipy import interpolate
from scipy import ndimage
import astropy.io.fits as pyfits
import astropy.stats as apyStats
import copy
import sys
import glob
import itertools
import pyximport; pyximport.install()
import nemoCython
import nemo
from . import maps
from . import signals
from . import photometry
from . import catalogs
from . import plotSettings
from . import gnfw
import astropy.table as atpy
import time
import IPython

#-------------------------------------------------------------------------------------------------------------
def filterMaps(unfilteredMapsDictList, filtersList, extNames = ['PRIMARY'], rootOutDir = ".", verbose = True):
    """Build and applies filters to the unfiltered maps(s). The output is a filtered map in yc or uK (this
    can be set with outputUnits in the config file). All filter operations are done in the filter objects, 
    even if multifrequency (a change from previous behaviour).
   
    Filtered maps are written to rootOutDir/filteredMaps
    Filters, if stored, are written to rootOutDir/filters
    
    Returns a dictionary containing a map of filtered maps to keys in filterDict. We'll use this dictionary
    for keeping track of catalogs etc. subsequently.
    
    """
    
    # Storage, in case it doesn't already exist
    filteredMapsDir=rootOutDir+os.path.sep+"filteredMaps"
    diagnosticsDir=rootOutDir+os.path.sep+"diagnostics"
    dirList=[filteredMapsDir, diagnosticsDir]
    for d in dirList:
        if os.path.exists(d) == False:
            os.makedirs(d)
            
    # Dictionary to keep track of images we're going to make
    imageDict={}
    
    # For handling tileDeck style .fits files
    imageDict['extNames']=extNames
    
    # Since we're putting stuff like extNames in the top level, let's keep a separate list of mapDicts
    imageDict['mapKeys']=[]
    
    # Make filtered maps for each filter
    if verbose == True: print(">>> Making filtered maps and S/N maps ...")
    for f in filtersList:
        
        # Iterate over all extensions (for tileDeck files)...
        for extName in extNames:
            
            print("--> extName = %s ..." % (extName))
            label=f['label']+"#"+extName
            
            filteredMapFileName=filteredMapsDir+os.path.sep+"%s_filteredMap.fits"  % (label)
            SNMapFileName=filteredMapsDir+os.path.sep+"%s_SNMap.fits" % (label)
            signalMapFileName=diagnosticsDir+os.path.sep+"%s_signalMap.fits" % (label)
            #transferFnFileName=filteredMapsDir+os.path.sep+"%s_transferFunction.fits" % (f['label'])

            if os.path.exists(filteredMapFileName) == False:
                
                print("... making filtered map %s ..." % (label)) 
                filterClass=eval('%s' % (f['class']))
                filterObj=filterClass(label, unfilteredMapsDictList, f['params'], \
                                      extName = extName, 
                                      diagnosticsDir = diagnosticsDir)
                filteredMapDict=filterObj.buildAndApply()
                    
                # Keywords we need for photometry later
                filteredMapDict['wcs'].header['BUNIT']=filteredMapDict['mapUnits']
                if 'beamSolidAngle_nsr' in filteredMapDict.keys() and filteredMapDict['beamSolidAngle_nsr'] > 0:
                    filteredMapDict['wcs'].header['BEAMNSR']=filteredMapDict['beamSolidAngle_nsr']
                    filteredMapDict['wcs'].header['FREQGHZ']=filteredMapDict['obsFreqGHz']
                filteredMapDict['wcs'].updateFromHeader()
                                           
                # Undo pixel window function using Sigurd's FFT method (takes into account variable pixel scale etc.)
                # We only need to do this for maps of signal (cancels in S/N map)
                # We do this once because it does take some time... and then we can forget about if e.g. stacking or doing forced photometry later
                print("... undoing pixel window function ...")
                mask=np.equal(filteredMapDict['data'], 0)
                filteredMapDict['data']=enmap.apply_window(filteredMapDict['data'], pow=-1.0)
                filteredMapDict['data'][mask]=0 # just in case we rely elsewhere on zero == no data

                # Write maps
                astImages.saveFITS(filteredMapFileName, filteredMapDict['data'], filteredMapDict['wcs'])
                astImages.saveFITS(SNMapFileName, filteredMapDict['SNMap'], filteredMapDict['wcs'])            

            else:
                print("... filtered map %s already made ..." % (label)) 
            
            # Add file names to imageDict
            if label not in imageDict:
                imageDict[label]={}
            imageDict[label]['filteredMap']=filteredMapFileName
            imageDict[label]['SNMap']=SNMapFileName
            imageDict[label]['signalMap']=signalMapFileName
            
            # Track e.g. reference filter scale with this key
            imageDict[label]['template']=f['label']
            
            # Track which keys have filtered maps that we might want to iterate over
            imageDict['mapKeys'].append(label)
            
            # May be handy to keep track of for plotting etc. later
            imageDict[label]['unfilteredMapsDictList']=unfilteredMapsDictList  
            
    return imageDict

#------------------------------------------------------------------------------------------------------------
class MapFilter(object):
    """Generic map filter base class. Defines common interface.
    
    """
    def __init__(self, label, unfilteredMapsDictList, paramsDict, extName = 'PRIMARY', writeFilter = False, 
                 forceRebuild = False, diagnosticsDir = None):
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
        # NOTE: we're now copying the input unfilteredMapsDictList, for supporting multi-ext tileDeck files
        self.unfilteredMapsDictList=[]
        for mapDict in unfilteredMapsDictList:           
            mapDict=maps.preprocessMapDict(mapDict.copy(), extName = extName, diagnosticsDir = diagnosticsDir)
            self.unfilteredMapsDictList.append(mapDict)
        self.wcs=mapDict['wcs']
        
        # Get beam solid angle info (units: nanosteradians)... we'll need for fluxes in Jy later
        self.beamSolidAnglesDict={}
        for mapDict in self.unfilteredMapsDictList:    
            beamFileName=mapDict['beamFileName']
            with open(beamFileName, "r") as inFile:
                lines=inFile.readlines()
                foundLine=False
                for line in lines:
                    if line.find("solid angle") != -1:
                        foundLine=True
                        break
                if foundLine == True:
                    bits=line.split("=")
                    solidAngle_nsr=float(bits[1].split("nsr")[0])
                else:
                    solidAngle_nsr=0.0
            self.beamSolidAnglesDict[mapDict['obsFreqGHz']]=solidAngle_nsr
        
        # For pixell / enmap
        # NOTE: enki maps can have an additional axis, which we don't want
        enheader=self.wcs.header
        if 'NAXIS3' in enheader.keys():
            del enheader['NAXIS3']
        enheader['NAXIS']=2
        self.enwcs=enwcs.WCS(enheader)

        # We could make this adjustable... added after switch to pixell
        self.apodPix=20
        
        # Sanity check that all maps are the same dimensions
        shape=self.unfilteredMapsDictList[0]['data'].shape
        for mapDict in self.unfilteredMapsDictList:
            if mapDict['data'].shape != shape:
                raise Exception("Maps at different frequencies have different dimensions!")
        
        # This is used by routines that make signal templates
        self.makeRadiansMap()
                                        
        # Set up storage if necessary, build this filter if not already stored
        self.diagnosticsDir=diagnosticsDir
        self.filterFileName=self.diagnosticsDir+os.path.sep+label+".fits"
               
        
    def makeRadiansMap(self):
        """Makes a map of distance in radians from centre, for the map being filtered.
        
        """
        
        mapDict=self.unfilteredMapsDictList[0]
        
        # NOTE: int conversion for python3
        x0=int(mapDict['data'].shape[1]/2)
        y0=int(mapDict['data'].shape[0]/2)
        ra0, dec0=self.wcs.pix2wcs(x0, y0)
        ra1, dec1=self.wcs.pix2wcs(x0+1, y0+1)
        self.degPerPixX=astCoords.calcAngSepDeg(ra0, dec0, ra1, dec0)
        self.degPerPixY=astCoords.calcAngSepDeg(ra0, dec0, ra0, dec1)
        
        # Real space map og angular distance from centre in radians, used in making filters and beam
        # NOTE: floor and int conversion added for python3
        xRadRange=np.array([np.arange(int(np.floor(-mapDict['data'].shape[1]/2)), int(mapDict['data'].shape[1]/2), \
                                    dtype=np.float64)*np.radians(self.degPerPixX)]*mapDict['data'].shape[0])
        yRadRange=np.array([np.arange(int(np.floor(-mapDict['data'].shape[0]/2)), int(mapDict['data'].shape[0]/2), \
                                    dtype=np.float64)*np.radians(self.degPerPixY)]*mapDict['data'].shape[1]).transpose()
        rRadRange=np.sqrt(xRadRange**2+yRadRange**2)
        self.radiansMap=rRadRange
        
        
    def buildAndApply(self):
        """Builds and applies the filter to the unfiltered map(s). Returns a dictionary, containing
        keys 'data', 'wcs', 'weights', 'obsFreqGHz'. If this routine converts to yc in place, the latter is set
        to 'yc'.
    
        """
        
        raise Exception("Called a base filter class without a buildAndApply() function implemented.")
        return None
    
        
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
                    print("... "+str(j*10)+"% complete ...")
            besselTransformedArray[i]=integrate.simps(profile*special.j0(lArray[i]*thetaArray)*thetaArray*thetaSmoothed, thetaArray)

        return [2*np.pi*besselTransformedArray, lArray]
                

    def makeForegroundsPower(self, obsFreqGHz, whiteNoiseLevel):
        """Returns a Power2D object with foregrounds power from the CMB and white noise only.
        
        White noise level is per pixel.
                       
        """
                    
        # CAMB power spec with Planck 2015 parameters (ish)        
        ps=powspec.read_spectrum(nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"planck_lensedCls.dat", scale = True)
        randMap=enmap.rand_map(self.unfilteredMapsDictList[0]['data'].shape, self.enwcs, ps)
        fMaskedData=enmap.fft(randMap)#(enmap.apod(randMap, self.apodPix))
        fgPower=np.real(fMaskedData*fMaskedData.conj())
        
        # Add some white noise to avoid ringing
        #whiteNoise=np.random.normal(0, whiteNoiseLevel, fgPowerMap.shape)
        #lm=liteMap.liteMapFromDataAndWCS(whiteNoise, self.wcs)
        #whiteNoisePower=fftTools.powerFromLiteMap(lm)
        #fgPower.powerMap=fgPower.powerMap+whiteNoisePower.powerMap
                
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

        # Changed for python3 - we may want to check this...
        # i.e., force to always be int after /2 without this
        x0=int(realSpace.shape[1]/2)
        y0=int(realSpace.shape[0]/2)

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
class MatchedFilter(MapFilter):
    """Matched filter...
    
    """

    def buildAndApply(self):
        
        print(">>> Building filter %s ..." % (self.label))
                    
        # Filter both maps with the same scale filter
        filteredMaps={}
        for mapDict in self.unfilteredMapsDictList:   
                        
            maskedData=mapDict['data']
            
            # Replaced flipper with pixell
            # NOTE: modLMap, modThetaMap are not exactly the same as flipper gives... doesn't seem to be python3 thing
            # Normalisation is different to flipper, but doesn't matter if we are consistent throughout
            fMaskedData=enmap.fft(enmap.apod(maskedData, self.apodPix))
            try:
                modThetaMap=180.0/(enmap.modlmap(fMaskedData.shape, self.enwcs)+1)
            except:
                print("shape business")
                IPython.embed()
                sys.exit()
                            
            # Make noise power spectrum
            if self.params['noiseParams']['method'] == 'dataMap':
                print("... taking noise power for filter from map ...")
                NP=np.real(fMaskedData*fMaskedData.conj())
            elif self.params['noiseParams']['method'] == 'CMBOnly':
                print("... taking noise power for filter from model CMB power spectrum ...")
                highPassMap=maps.subtractBackground(maskedData, self.wcs, 0.25/60.0)
                whiteNoiseLevel=np.std(highPassMap)
                NP=self.makeForegroundsPower(mapDict['obsFreqGHz'], whiteNoiseLevel)
            elif self.params['noiseParams']['method'] == 'max(dataMap,CMB)':
                print("... taking noise power from max(map noise power, model CMB power spectrum) ...")
                NP=np.real(fMaskedData*fMaskedData.conj())
                NPCMB=self.makeForegroundsPower(mapDict['obsFreqGHz'], 0.0)
                NP=np.maximum.reduce([NP, NPCMB])
            else:
                raise Exception("noise method must be either 'dataMap', '4WayMaps', 'CMBOnly', or 'max(dataMap,CMB)'")
            
            # Save plot of 2d noise power (to easily compare different methods / frequencies)
            # NOTE: output path here will break if we're not using RealSpaceMatchedFilter?
            plotSettings.update_rcParams()
            plt.figure(figsize=(9.5,9.5))
            ax=plt.axes([0.05, 0.05, 0.9, 0.9])
            plotData=np.log10(fft.fftshift(NP))
            #cutImage=astImages.intensityCutImage(plotData, ['relative', 99.5])
            cutImage=astImages.intensityCutImage(plotData, [plotData.min(), plotData.max()])
            plt.imshow(cutImage['image'], cmap = 'gray', norm = cutImage['norm'], origin = 'lower')
            plt.xticks([], [])
            plt.yticks([], [])
            plt.savefig(self.diagnosticsDir+os.path.sep+self.label.replace("realSpaceKernel", "noisePower%.1f" % (mapDict['obsFreqGHz']))+".png")
            plt.close()
            
            # FFT of signal
            signalMapDict=self.makeSignalTemplateMap(mapDict['beamFileName'], mapDict['obsFreqGHz'])
            signalMap=signalMapDict['signalMap']
            fftSignal=enmap.fft(signalMap)
            
            # Toby style -note smoothing noise is essential!
            filt=1.0/NP
            med=np.median(filt)
            filt[np.where(filt>10*med)]=med
            kernelSize=(5,5)
            filt=ndimage.gaussian_filter(filt, kernelSize)
            filt=filt*abs(fftSignal)
            cov=filt*np.abs(fftSignal)**2
            integral=cov.sum()/float(signalMap.shape[0])/float(signalMap.shape[1])
            filt=filt/integral
            self.G=filt
            
            # NOTE: Do not disable this weighting
            weightedMap=mapDict['data']*np.sqrt(mapDict['weights']/mapDict['weights'].max())
            weightedMap=enmap.apod(weightedMap, self.apodPix)
            fMap=enmap.fft(weightedMap)
            
            # Check and correct for bias - apply filter to the signal template - do we recover the same flux?
            # Rescale the filtered map so that this is taken care of (better way than this, but this
            # should be pretty good)
            peakCoords=np.where(abs(signalMap) == abs(signalMap).max())
            yPeak=peakCoords[0][0]
            xPeak=peakCoords[1][0]
            fSignalMap=enmap.fft(signalMap)
            filteredSignalMap=np.real(enmap.ifft(fSignalMap*self.G))

            # Use the signal map we made using MatchedFilter to figure out how much it has been rolled off by the filter
            # NOTE: this is being done the same way in RealSpaceMatchedFilter - add a function instead...
            signalProperties=signalMapDict['inputSignalProperties']
            if self.params['outputUnits'] == 'yc':
                # Normalise such that peak value in filtered map == y0, taking out the effect of the beam
                filteredSignalMap=maps.convertToY(filteredSignalMap, obsFrequencyGHz = signalProperties['obsFreqGHz'])
                signalNorm=signalProperties['y0']/filteredSignalMap.max()
            elif self.params['outputUnits'] == 'Y500':
                # Normalise such that peak value in filtered map == Y500 for the input SZ cluster model
                # We can get this from yc and the info in signalProperties, so we don't really need this
                print("implement signal norm for %s" % (self.params['outputUnits']))
                IPython.embed()
                sys.exit()
            elif self.params['outputUnits'] == 'uK':
                # Normalise such that peak value in filtered map == peak value of source in uK
                signalNorm=1.0/filteredSignalMap.max()
            elif self.params['outputUnits'] == 'Jy/beam':
                # Normalise such that peak value in filtered map == flux density of the source in Jy/beam
                print("implement signal norm for %s" % (self.params['outputUnits']))
                IPython.embed()
                sys.exit()
            else:
                raise Exception("didn't understand 'outputUnits' given in the .par file")

            filteredMaps['%d' % int(mapDict['obsFreqGHz'])]=np.real(enmap.ifft(fMap*self.G, normalize = False))*signalNorm
                
            # Linearly combine filtered maps (optional)
            if 'mapCombination' in list(self.params.keys()):
                combinedMap=np.zeros(filteredMaps[list(filteredMaps.keys())[0]].shape)
                for key in list(filteredMaps.keys()):
                    combinedMap=combinedMap+self.params['mapCombination'][key]*filteredMaps[key]
                combinedObsFreqGHz=self.params['mapCombination']['rootFreqGHz']
            else:
                # If no linear map combination given, assume we only want the first item in the filtered maps list
                combinedObsFreqGHz=self.unfilteredMapsDictList[0]['obsFreqGHz']
                combinedMap=filteredMaps['%d' % combinedObsFreqGHz]

            # Convert to whatever output units we want:
            # Jy/sr (should go to Jy/beam eventually) for sources
            # yc for SZ clusters
            if 'outputUnits' in list(self.params.keys()):
                if self.params['outputUnits'] == 'yc':
                    combinedMap=maps.convertToY(combinedMap, combinedObsFreqGHz)
                    combinedObsFreqGHz='yc'
                    mapUnits='yc'
                elif self.params['outputUnits'] == 'uK':
                    mapUnits='uK'
                elif self.params['outputUnits'] == 'Jy/beam':
                    print("Jy/beam here")
                    mapUnits='Jy/beam'
                    IPython.embed()
                    sys.exit()
                else:
                    raise Exception('need to specify "outputUnits" ("yc", "uK", or "Jy/beam") in filter params')
                
            # We'll be saving this shortly...
            apodMap=enmap.apod(np.ones(maskedData.shape), self.apodPix)
            apodMask=np.zeros(apodMap.shape)
            apodMask[np.greater(apodMap, 0.999999)]=1.0
            
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
            
            maskFileName=self.diagnosticsDir+os.path.sep+"areaMask#%s.fits" % (mapDict['extName'])
            if os.path.exists(maskFileName) == False:
                astImages.saveFITS(maskFileName, apodMask, mapDict['wcs'])
                    
            # Blank out the apodised area
            combinedMap=combinedMap*apodMask
            SNMap=SNMap*apodMask
                                        
            # Save filter profile in real space
            self.saveRealSpaceFilterProfile()        
            
        return {'data': combinedMap, 'wcs': self.wcs, 'obsFreqGHz': combinedObsFreqGHz, 'SNMap': SNMap, 
                'signalMap': signalMap, 'mapUnits': mapUnits, 
                'inputSignalProperties': signalMapDict['inputSignalProperties']}
            
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

    def loadKernel(self, kern2DFileName):
        """Loads a previously cached kernel.
        
        Returns kern2d, signalNorm, bckSubScaleArcmin
        
        """
        
        print("... loading previously cached kernel %s ..." % (kern2DFileName))
        img=pyfits.open(kern2DFileName)
        kern2d=img[0].data
        signalNorm=img[0].header['SIGNORM']
        bckSubScaleArcmin=img[0].header['BCKSCALE']

        return kern2d, signalNorm, bckSubScaleArcmin

        
    def buildKernel(self, mapDict, RADecSection, RADeg = 'centre', decDeg = 'centre'):
        """Builds the real space kernel itself. 
        
        RADeg, decDeg are used for figuring out pixel scales for background subtraction
        
        Returns kern2d, signalNorm, bckSubScaleArcmin
        
        """
        
        # We will cache the kernel... we also may want to re-load it later (for e.g. sky sim contamination tests)
        freqLabel=int(mapDict['obsFreqGHz'])
        kern2DFileName=self.diagnosticsDir+os.path.sep+"kern2d_%s_%d.fits" % (self.label, freqLabel)
        if os.path.exists(kern2DFileName) == True:
            return self.loadKernel(kern2DFileName)
        
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
            if k in list(mapDict.keys()):
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
                                         extName = mapDict['extName'],
                                         diagnosticsDir = matchedFilterDir+os.path.sep+'diagnostics')
        filteredMapDict=matchedFilter.buildAndApply()
                
        # Turn the matched filter into a smaller real space convolution kernel
        # This means we have to roll off the kernel to 0 at some radius
        # This is set by maxArcmin in the .par file
        prof, arcminRange=matchedFilter.makeRealSpaceFilterProfile()
        rIndex=np.where(arcminRange > kernelMaxArcmin)[0][0]
        mask=np.less(arcminRange, kernelMaxArcmin)

        # Alternatively, roll off to zero after the second zero crossing
        # NOTE: now setting in the .par file, uncomment below to switch back
        #segProf=ndimage.label(np.greater(prof, 0))[0]
        #rIndex=np.where(segProf == 2)[0][0]
        #maxArcmin=arcminRange[rIndex]

        # Kernel can be either fully 2d, or be azimuthally averaged... in the ACTPol E-D56 paper, we used the latter
        if self.params['noiseParams']['symmetrize'] == False:
            profile2d=fft.ifft2(matchedFilter.G).real
            profile2d=fft.fftshift(profile2d)
        else:
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
        if 'bckSubScaleArcmin' in self.params.keys():
            bckSubScaleArcmin=self.params['bckSubScaleArcmin']
        else:
            bckSubScaleArcmin=arcminRange[prof == prof.min()][0]
        
        # Use the signal map we made using MatchedFilter to figure out how much it has been rolled off by:
        # 1. The high pass filter (bck sub step)
        # 2. The matched filter itself (includes beam)
        signalMap=filteredMapDict['signalMap']      # Note that this has had the beam applied already
        signalProperties=filteredMapDict['inputSignalProperties']
        # Should add an applyKernel function to do all this
        if self.params['bckSub'] == True:
            filteredSignal=maps.subtractBackground(signalMap, wcs, RADeg = RADeg, decDeg = decDeg,
                                                       smoothScaleDeg = bckSubScaleArcmin/60.0)
        else:
            filteredSignal=np.zeros(signalMap.shape)+signalMap
        filteredSignal=ndimage.convolve(filteredSignal, kern2d) 
        if self.params['outputUnits'] == 'yc':
            # Normalise such that peak value in filtered map == y0, taking out the effect of the beam
            filteredSignal=maps.convertToY(filteredSignal, obsFrequencyGHz = signalProperties['obsFreqGHz'])
            signalNorm=signalProperties['y0']/filteredSignal.max()
        elif self.params['outputUnits'] == 'Y500':
            # Normalise such that peak value in filtered map == Y500 for the input SZ cluster model
            # We can get this from yc and the info in signalProperties, so we don't really need this
            print("implement signal norm for %s" % (self.params['outputUnits']))
            IPython.embed()
            sys.exit()
        elif self.params['outputUnits'] == 'uK':
            # Normalise such that peak value in filtered map == peak value of source in uK
            # We take out the effect of the pixel window function here - this assumes we're working with sources
            # NOTE: for amplitude, the mappers want the amplitude of the delta function, not the beam
            # NOTE: we're now undoing the pixel window function in filterMaps (at the top)
            #signalNorm=1.0/(signalProperties['pixWindowFactor']*filteredSignal.max())
            signalNorm=1.0/filteredSignal.max()
        elif self.params['outputUnits'] == 'Jy/beam':
            # Normalise such that peak value in filtered map == flux density of the source in Jy/beam
            print("implement signal norm for %s" % (self.params['outputUnits']))
            IPython.embed()
            sys.exit()
        else:
            raise Exception("didn't understand 'outputUnits' given in the .par file")
        
        # Save 2d kernel - we need this (at least for the photometry ref scale) to calc Q later
        # Add bckSubScaleArcmin to the header
        kernWCS=wcs.copy()
        if self.params['bckSub'] == True:
            kernWCS.header['BCKSCALE']=bckSubScaleArcmin
        kernWCS.header['SIGNORM']=signalNorm
        RADecLabel=str(RADecSection).replace(",", "_").replace(" ", "").replace("[", "").replace("]", "").replace(".", "p")
        astImages.saveFITS(kern2DFileName, kern2d, kernWCS)
        
        # Filter profile plot   
        # Save the stuff we plot first, in case we want to make a plot with multiple filters on later
        np.savez(self.diagnosticsDir+os.path.sep+"filterProf1D_%s_%d.npz" % (self.label, freqLabel), 
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
        plt.savefig(self.diagnosticsDir+os.path.sep+"filterPlot1D_%s_%d.pdf" % (self.label, freqLabel))
        plt.close()
        
        return kern2d, signalNorm, bckSubScaleArcmin


    def makeNoiseMap(self, mapData):
        """Estimate the noise map using local RMS measurements on a grid, over the whole filtered map.
        
        """

        #print "... making SN map ..."
        gridSize=int(round((self.params['noiseParams']['noiseGridArcmin']/60.)/self.wcs.getPixelSizeDeg()))
        #gridSize=rIndex*3
        overlapPix=int(gridSize/2)
        numXChunks=mapData.shape[1]/gridSize
        numYChunks=mapData.shape[0]/gridSize
        yChunks=np.linspace(0, mapData.shape[0], numYChunks+1, dtype = int)
        xChunks=np.linspace(0, mapData.shape[1], numXChunks+1, dtype = int)
        #SNMap=np.zeros(mapData.shape)
        apodMask=np.not_equal(mapData, 0)
        # We could make below behaviour default if match photFilter? Would need to see photFilter though...
        #if 'saveRMSMap' in self.params['noiseParams'] and self.params['noiseParams']['saveRMSMap'] == True:
        RMSMap=np.zeros(mapData.shape)
        for i in range(len(yChunks)-1):
            for k in range(len(xChunks)-1):
                y0=yChunks[i]-overlapPix
                y1=yChunks[i+1]+overlapPix
                x0=xChunks[k]-overlapPix
                x1=xChunks[k+1]+overlapPix
                if y0 < 0:
                    y0=0
                if y1 > mapData.shape[0]:
                    y1=mapData.shape[0]
                if x0 < 0:
                    x0=0
                if x1 > mapData.shape[1]:
                    x1=mapData.shape[1]
                chunkValues=mapData[y0:y1, x0:x1]

                goodAreaMask=np.greater_equal(apodMask[y0:y1, x0:x1], 1.0)
                
                if 'RMSEstimator' in self.params['noiseParams'].keys() and self.params['noiseParams']['RMSEstimator'] == 'biweight':
                    if goodAreaMask.sum() >= 10:
                        # Astropy version is faster but gives identical results
                        chunkRMS=apyStats.biweight_scale(chunkValues[goodAreaMask], c = 9.0, modify_sample_size = True)
                        #chunkRMS=astStats.biweightScale(chunkValues[goodAreaMask], 6.0)
                    else:
                        chunkRMS=0.
                else:
                    # Default: 3-sigma clipped stdev
                    if np.not_equal(chunkValues, 0).sum() != 0:
                        goodAreaMask=np.greater_equal(apodMask[y0:y1, x0:x1], 1.0)
                        chunkMean=np.mean(chunkValues[goodAreaMask])
                        chunkRMS=np.std(chunkValues[goodAreaMask])
                        sigmaClip=3.0
                        for c in range(10):
                            mask=np.less(abs(chunkValues), abs(chunkMean+sigmaClip*chunkRMS))
                            mask=np.logical_and(goodAreaMask, mask)
                            if mask.sum() > 0:
                                chunkMean=np.mean(chunkValues[mask])
                                chunkRMS=np.std(chunkValues[mask])
                    else:
                        chunkRMS=0.
                
                if chunkRMS > 0:
                    RMSMap[y0:y1, x0:x1]=chunkRMS
                
        return RMSMap
    
    
    def _makeCorrMap(self, dataCube, outFileName):
        """This is a test: not generalized, just for 90 + 150 GHz for now. Records off-diagonal correlation coeff
        measured in the same cells as makeNoiseMap uses. Saves as a map in diagnosticsDir.
        
        """
                
        #print "... making SN map ..."
        gridSize=int(round((self.params['noiseParams']['noiseGridArcmin']/60.)/self.wcs.getPixelSizeDeg()))
        #gridSize=rIndex*3
        overlapPix=int(gridSize/2)
        numXChunks=dataCube[0].shape[1]/gridSize
        numYChunks=dataCube[0].shape[0]/gridSize
        yChunks=np.linspace(0, dataCube[0].shape[0], numYChunks+1, dtype = int)
        xChunks=np.linspace(0, dataCube[0].shape[1], numXChunks+1, dtype = int)
        #SNMap=np.zeros(mapData.shape)
        apodMask=np.logical_and(np.not_equal(dataCube[0], 0), np.not_equal(dataCube[1], 0))
        # We could make below behaviour default if match photFilter? Would need to see photFilter though...
        #if 'saveRMSMap' in self.params['noiseParams'] and self.params['noiseParams']['saveRMSMap'] == True:
        covMap=np.zeros(dataCube[0].shape)  # We're assuming 90, 150 GHz only, so we only need to collect one of the off diagonal terms
        for i in range(len(yChunks)-1):
            for k in range(len(xChunks)-1):
                y0=yChunks[i]-overlapPix
                y1=yChunks[i+1]+overlapPix
                x0=xChunks[k]-overlapPix
                x1=xChunks[k+1]+overlapPix
                if y0 < 0:
                    y0=0
                if y1 > dataCube[0].shape[0]:
                    y1=dataCube[0].shape[0]
                if x0 < 0:
                    x0=0
                if x1 > dataCube[0].shape[1]:
                    x1=dataCube[0].shape[1]
                chunkValues=dataCube[:, y0:y1, x0:x1]
                goodAreaMask=np.greater_equal(apodMask[y0:y1, x0:x1], 1.0)
                cov=np.corrcoef(chunkValues[0][goodAreaMask], chunkValues[1][goodAreaMask])
                covMap[y0:y1, x0:x1]=cov[0, 1]
        astImages.saveFITS(outFileName, covMap, self.wcs)
    
    
    def makeSZMap(self, filteredMaps):
        """Combines maps at multiple frequencies (or just one) to make a Compton yc map, using inverse
        variance weighting of the single frequency maps.
        
        Note:
            This routine will write a FITS datacube containing the weights used to the `diagnostics/`
            directory, if `saveFreqWeightMap` is set to True in `self.params`.
            
        Args:
            filteredMaps: A dictionary containing maps. Each key gives the corresponding map frequency 
                in GHz.
        
        Returns:
            yc map, noise map, signal-to-noise map
        
        """
        
        dataCube=[]
        noiseCube=[]
        obsFreqs=[]
        for key in list(filteredMaps.keys()):
            obsFreqs.append(float(key))
            mapData=filteredMaps[key]
            mapData=maps.convertToY(mapData, float(key))
            RMSMap=self.makeNoiseMap(mapData)
            dataCube.append(mapData)
            noiseCube.append(RMSMap)
        dataCube=np.array(dataCube)
        noiseCube=np.array(noiseCube)
        
        # Test: how correlated are the filtered maps? This is currently hardcoded for two frequencies...
        if dataCube.ndim == 3 and 'saveCorrMap' in self.params.keys() and self.params['saveCorrMap'] == True:
            self._makeCorrMap(dataCube, self.diagnosticsDir+os.path.sep+"corrMap_signal.fits")
            self._makeCorrMap(noiseCube, self.diagnosticsDir+os.path.sep+"corrMap_noiseLevel.fits")

        # Zap very low values as rounding errors (so have zero weight in combination)
        noiseCube[np.less(noiseCube, 1e-7)]=0.
                
        # Combining - inverse variance weighted average
        if dataCube.shape[0] > 1:
            invVar=1./noiseCube**2
            invVar[np.isinf(invVar)]=0.
            # Use mask to flag pixels with zero everywhere
            zeroMask=np.equal(np.sum(invVar, axis = 0), 0)
            for z in range(invVar.shape[0]):
                invVar[z][zeroMask]=1.
            # Make relative weights cube from inv var (for convenience really)
            # Since we used inverse variance for making a combined y0~ estimate,
            # we need these to weight fRel in mass estimates (see signals.calcM500Fromy0)
            relWeightsCube=np.zeros(dataCube.shape)
            for i in range(invVar.shape[0]):
                relWeightsCube[i]=invVar[i]/np.sum(invVar, axis = 0)
            weightedMap=np.average(dataCube, weights = invVar, axis = 0)
            weightedNoise=np.sqrt(1./np.sum(invVar, axis = 0))
            weightedNoise[np.isinf(weightedNoise)]=0.
        else:
            weightedMap=dataCube[0]
            weightedNoise=noiseCube[0]
            relWeightsCube=np.ones(dataCube.shape)
        
        weightedSNMap=np.zeros(weightedMap.shape)
        mask=np.greater(weightedNoise, 0)
        weightedSNMap[mask]=weightedMap[mask]/weightedNoise[mask]
        
        # Write out relative weights cube as .fits - handy not just for fRel, but also for coverage checking
        # i.e., we can put this info into object catalogs - indicate where map was 148 GHz only, for example
        # We store frequency info in FITS header as 0FREQGHZ 1FREQGHZ etc.
        if 'saveFreqWeightMap' in self.params and self.params['saveFreqWeightMap'] == True:
            cubeWCS=self.wcs.copy()
            count=0
            for obsFreqGHz in obsFreqs:
                cubeWCS.header['%dFREQGHZ' % (count)]=obsFreqGHz
                count=count+1
            outFileName=self.diagnosticsDir+os.path.sep+"freqRelativeWeights_%s.fits" % (self.label)
            astImages.saveFITS(outFileName, relWeightsCube, cubeWCS)
        
        return weightedMap, weightedNoise, weightedSNMap
                
            
    def buildAndApply(self):
        
        print(">>> Building filter %s ..." % (self.label))
        
        # We don't want to get rid of this option for the regular matched filter...
        # BUT we don't want it here... because of the way we're implementing the multi-frequency filter for SZ searches
        if 'mapCombination' in list(self.params.keys()):
            raise Exception("mapCombination key should not be given in .yml file for RealSpaceMatchedFilter")

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
            if self.params['noiseParams']['RADecSection'] == 'tileNoiseRegions':
                RADecSectionDictList=[{'RADecSection': [wcs.header['NRAMIN'], wcs.header['NRAMAX'], 
                                                        wcs.header['NDEMIN'], wcs.header['NDEMAX']],
                                       'applyDecMin': decMin, 'applyDecMax': decMax,
                                       'applyRAMin': RAMin, 'applyRAMax': RAMax}]
                print("... taking noise from tileDeck image header: %s ..." % (RADecSectionDictList[0]['RADecSection']))
                
            elif self.params['noiseParams']['RADecSection'] == 'auto':
                
                if self.params['noiseParams']['method'] != 'CMBOnly':
                    raise Exception("'auto' option for RADecSection disabled for all methods EXCEPT CMBOnly (needs fixing)")

                # NOTE: hardcoded max size... if we're using CMBOnly, we don't care (no real noise) 
                cRADeg, cDecDeg=wcs.getCentreWCSCoords()
                maxSizeDeg=5.0
                kernelBuildRAMin, kernelBuildRAMax, kernelBuildDecMin, kernelBuildDecMax=astCoords.calcRADecSearchBox(cRADeg, cDecDeg, maxSizeDeg)
                
                RADecSectionDictList=[{'RADecSection': [kernelBuildRAMin, kernelBuildRAMax, 
                                                        kernelBuildDecMin, kernelBuildDecMax],
                                       'applyDecMin': decMin, 'applyDecMax': decMax,
                                       'applyRAMin': RAMin, 'applyRAMax': RAMax}]
                #print "... taking noise from %s ..." % (RADecSectionDictList[0]['RADecSection'])
                
            elif self.params['noiseParams']['RADecSection'][2] == 'numDecSteps':
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
                kern2d, signalNorm, bckSubScaleArcmin=self.buildKernel(mapDict, RADecSectionDict['RADecSection'],
                                                                       RADeg = applyRACentre, decDeg = applyDecCentre)
                RADecSectionDict['kern2d']=kern2d
                RADecSectionDict['signalNorm']=signalNorm
                RADecSectionDict['bckSubScaleArcmin']=bckSubScaleArcmin
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
                    mapData[yMin:yMax+yOverlap, :]=maps.subtractBackground(mapData[yMin:yMax+yOverlap, :], wcs, 
                                                                      RADeg = RADecSectionDict['applyRACentre'], 
                                                                      decDeg = RADecSectionDict['applyDecCentre'],
                                                                      smoothScaleDeg = RADecSectionDict['bckSubScaleArcmin']/60.)
                    mapData[yMax:yMax+yOverlap]=buff
            if 'saveHighPassMap' in self.params and self.params['saveHighPassMap'] == True:
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
                print("... convolving map with kernel [%d:%d] ..." % (yMin, yMax))
                mapData[yMin:yMax+yOverlap, :]=ndimage.convolve(mapData[yMin:yMax+yOverlap, :], RADecSectionDict['kern2d'])   
                filtBuff=np.zeros(mapData[yMax:yMax+yOverlap].shape)+mapData[yMax:yMax+yOverlap]
                mapData[yMax:yMax+yOverlap]=buff
                t1=time.time()
                print("... took %.3f sec ..." % (t1-t0))
                # Apply the normalisation
                mapData[yMin:yMax, :]=mapData[yMin:yMax, :]*RADecSectionDict['signalNorm']
            
            # Apply the point source mask here (before noise estimates etc.)
            mapData=mapData*psMask
                        
            #filteredMaps['%d' % int(mapDict['obsFreqGHz'])]=filteredMap
            filteredMaps['%d' % int(mapDict['obsFreqGHz'])]=mapData
        
        # We don't need the filteredMaps dict from here, so zap it to save memory
        #del filteredMaps
        
        # Convert to whatever output units we want:
        # Jy/sr (should go to Jy/beam eventually) for sources
        # yc for SZ clusters
        if 'outputUnits' in list(self.params.keys()):
            if self.params['outputUnits'] == 'yc':
                # Output of below is already in yc
                combinedMap, RMSMap, SNMap=self.makeSZMap(filteredMaps)
                combinedObsFreqGHz='yc'
                mapUnits='yc'
                beamSolidAngle_nsr=0.0   # not used...
            elif self.params['outputUnits'] == 'uK':
                if len(list(filteredMaps.keys())) > 1:
                    raise Exception("multi-frequency filtering not currently supported for outputUnits 'uK' (point source finding)")
                keyFreq=list(filteredMaps.keys())[0]
                combinedMap=filteredMaps[keyFreq]
                combinedObsFreqGHz=float(keyFreq)
                RMSMap=self.makeNoiseMap(combinedMap)
                validMask=np.greater(RMSMap, 0)
                SNMap=np.zeros(combinedMap.shape)+combinedMap
                SNMap[validMask]=SNMap[validMask]/RMSMap[validMask]
                #SNMap[np.isinf(SNMap)]=0.
                mapUnits='uK'
                beamSolidAngle_nsr=self.beamSolidAnglesDict[combinedObsFreqGHz]
            else:
                raise Exception('need to specify "outputUnits" ("yc" or "uK") in filter params')

        # Use rank filter to zap edges where RMS will be artificially low - we use a bit of a buffer here
        # NOTE: Now point source mask is applied above, we fill the holes back in here when finding edges
        if 'edgeTrimArcmin' in self.params.keys():
            trimSizePix=int(round((self.params['edgeTrimArcmin']/60.)/self.wcs.getPixelSizeDeg()))
        else:
            gridSize=int(round((self.params['noiseParams']['noiseGridArcmin']/60.)/self.wcs.getPixelSizeDeg()))
            trimSizePix=int(round(gridSize*3.0))
        edgeCheck=ndimage.rank_filter(abs(combinedMap+(1-psMask)), 0, size = (trimSizePix, trimSizePix))
        edgeCheck=np.array(np.greater(edgeCheck, 0), dtype = float)
        combinedMap=combinedMap*edgeCheck
        apodMask=np.not_equal(combinedMap, 0)
        surveyMask=edgeCheck*surveyMask*psMask
        del edgeCheck

        # Apply final survey mask to signal-to-noise map and RMS map
        # NOTE: need to avoid NaNs in here, otherwise map interpolation for e.g. S/N will fail later on
        SNMap=SNMap*surveyMask
        SNMap[np.isnan(SNMap)]=0.
        RMSMap=RMSMap*surveyMask

        maskFileName=self.diagnosticsDir+os.path.sep+"areaMask#%s.fits" % (mapDict['extName'])
        if os.path.exists(maskFileName) == False:
            astImages.saveFITS(maskFileName, np.array(surveyMask, dtype = int), mapDict['wcs'])
        
        if 'saveRMSMap' in self.params and self.params['saveRMSMap'] == True:
            RMSFileName=self.diagnosticsDir+os.path.sep+"RMSMap_%s.fits" % (self.label)
            astImages.saveFITS(RMSFileName, RMSMap, mapDict['wcs'])

        return {'data': combinedMap, 'wcs': self.wcs, 'obsFreqGHz': combinedObsFreqGHz,
                'SNMap': SNMap, 'mapUnits': mapUnits, 'beamSolidAngle_nsr': beamSolidAngle_nsr}
                
#------------------------------------------------------------------------------------------------------------
class BeamFilter(MapFilter):
    """Base class for filters using beam profile files in Matthew + Kavi's format.
        
    """
    
    def makeSignalTemplateMap(self, beamFileName, mapObsFreqGHz):
        """Makes a beam model signal template map.
        
        """
        
        signalMap, inputSignalProperties=signals.makeBeamModelSignalMap(mapObsFreqGHz, np.degrees(self.radiansMap),
                                                                          self.wcs, 
                                                                          beamFileName)

        return {'signalMap': signalMap, 'inputSignalProperties': inputSignalProperties}
    
#------------------------------------------------------------------------------------------------------------
class ArnaudModelFilter(MapFilter):
    """Base class for filters using the GNFW profile as described in Arnaud et al. (2010).
    
    """
    
    def makeSignalTemplateMap(self, beamFileName, mapObsFreqGHz):
        """Makes a model signal template map.
        
        Returns dictionary of {'signalMap', 'inputSignalProperties'}
        
        """
        
        signalMap, modelDict=signals.makeArnaudModelSignalMap(self.params['z'], self.params['M500MSun'], 
                                                                mapObsFreqGHz, np.degrees(self.radiansMap),
                                                                self.wcs, beamFileName, 
                                                                GNFWParams = self.params['GNFWParams'])
        
        return {'signalMap': signalMap, 'inputSignalProperties': modelDict}
                
#------------------------------------------------------------------------------------------------------------
# Definitions of actual filters that can be used
class ArnaudModelMatchedFilter(MatchedFilter, ArnaudModelFilter):
    pass
class BeamMatchedFilter(MatchedFilter, BeamFilter):
    pass

class ArnaudModelRealSpaceMatchedFilter(RealSpaceMatchedFilter, ArnaudModelFilter):
    pass
class BeamRealSpaceMatchedFilter(RealSpaceMatchedFilter, BeamFilter):
    pass
