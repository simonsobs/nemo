"""

Filter classes are defined in this module, together with the `filterMaps` function that uses them.

There are two main classes of filter: `MatchedFilter` and `RealSpaceMatchedFilter`.

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
def filterMaps(unfilteredMapsDictList, filtersList, tileNames = ['PRIMARY'], rootOutDir = ".", verbose = True,
               undoPixelWindow = True):
    """Build and applies filters to the unfiltered maps(s). The output is a filtered map in yc or uK (this
    can be set with outputUnits in the config file). All filter operations are done in the filter objects, 
    even if multifrequency (a change from previous behaviour).
   
    Filtered maps are written to rootOutDir/filteredMaps
    Filters, if stored, are written to rootOutDir/filters
    
    Returns a dictionary containing a map of filtered maps to keys in filterDict. We'll use this dictionary
    for keeping track of catalogs etc. subsequently.
    
    """
    
    # Storage, in case it doesn't already exist (these are created by NemoConfig in startUp.py anyway)
    filteredMapsDir=rootOutDir+os.path.sep+"filteredMaps"
    diagnosticsDir=rootOutDir+os.path.sep+"diagnostics"
    selFnDir=rootOutDir+os.path.sep+"selFn"
    dirList=[filteredMapsDir, diagnosticsDir, selFnDir]
    for d in dirList:
        if os.path.exists(d) == False:
            os.makedirs(d)
            
    # Dictionary to keep track of images we're going to make
    imageDict={}
    
    # For handling tileDeck style .fits files
    imageDict['tileNames']=tileNames
    
    # Since we're putting stuff like tileNames in the top level, let's keep a separate list of mapDicts
    imageDict['mapKeys']=[]
    
    # Make filtered maps for each filter
    if verbose == True: print(">>> Making filtered maps and S/N maps ...")
    for f in filtersList:
        
        # Iterate over all extensions (for tileDeck files)...
        for tileName in tileNames:
            
            print("--> tileName = %s ..." % (tileName))
            
            # This is the label tracked in imageDict, catalog template column
            # Should NOT be fed in as filterClass label
            label=f['label']+"#"+tileName
            
            filteredMapFileName=filteredMapsDir+os.path.sep+"%s_filteredMap.fits"  % (label)
            SNMapFileName=filteredMapsDir+os.path.sep+"%s_SNMap.fits" % (label)
            signalMapFileName=diagnosticsDir+os.path.sep+"%s_signalMap.fits" % (label)

            if os.path.exists(filteredMapFileName) == False:
                
                print("... making filtered map %s ..." % (label)) 
                filterClass=eval('%s' % (f['class']))
                filterObj=filterClass(f['label'], unfilteredMapsDictList, f['params'], \
                                      tileName = tileName, 
                                      diagnosticsDir = diagnosticsDir,
                                      selFnDir = selFnDir)
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
                if undoPixelWindow == True:
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
    def __init__(self, label, unfilteredMapsDictList, paramsDict, tileName = 'PRIMARY', writeFilter = False, 
                 forceRebuild = False, diagnosticsDir = None, selFnDir = None):
        """Initialises a MapFilter. unfilteredMapsDictList describes the input maps, paramsDict describes
        the filter options. The convention is that single frequency only filters (e.g. GaussianWienerFilter)
        only operate on the first map in the unfilteredMapDictList.
        
        label is used to store the filter to save calculating it again. Filters are stored under filters/
        dir. If a filter already exists under this dir with the given label, it is loaded rather than
        recalculated. To force recalculating the filter, set forceRebuild == True.
                        
        """
        
        self.label=label
        self.params=paramsDict
        
        # Set up storage if necessary, build this filter if not already stored
        self.diagnosticsDir=diagnosticsDir
        self.selFnDir=selFnDir
        self.tileName=tileName
        self.filterFileName=self.diagnosticsDir+os.path.sep+"filter_%s#%s.fits" % (self.label, self.tileName)
        
        # Prepare all the unfilteredMaps (in terms of cutting sections, masks etc.)
        # NOTE: we're now copying the input unfilteredMapsDictList, for supporting multi-ext tileDeck files
        self.unfilteredMapsDictList=[]
        for mapDict in unfilteredMapsDictList:           
            mapDict=maps.preprocessMapDict(mapDict.copy(), tileName = tileName, diagnosticsDir = diagnosticsDir)
            self.unfilteredMapsDictList.append(mapDict)
        self.wcs=mapDict['wcs']
        self.shape=mapDict['data'].shape
                                
        # Get beam solid angle info (units: nanosteradians)... we'll need for fluxes in Jy later
        self.beamSolidAnglesDict={}
        for mapDict in self.unfilteredMapsDictList:    
            if 'solidAngle_nsr' in mapDict.keys():
                solidAngle_nsr=mapDict['solidAngle_nsr']
            else:
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
        
        # This is a default - will get modified (calculated) by buildAndApply
        self.signalNorm=1.0
        self.fRelWeights={}
               
        
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
        
        # Real space map of angular distance from centre in radians, used in making filters and beam
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
    
    
    def makeForegroundsPower(self):
        """Returns a Power2D object with foregrounds power from the CMB (this could be extended to add other
        foregrounds if necessary).
                       
        """

        # CAMB power spec with Planck 2015 parameters (ish)    
        tab=atpy.Table().read(nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"planck_lensedCls.dat", format = 'ascii')
        tab['TT']=(tab['TT']*2*np.pi)/(tab['L']*(tab['L']+1))
        lmap=enmap.modlmap(self.unfilteredMapsDictList[0]['data'].shape, self.enwcs)
        l2p=interpolate.interp1d(tab['L'], tab['TT'], bounds_error=False, fill_value=0.0)
        fgPower=l2p(lmap)*lmap.shape[0]*lmap.shape[1]
                            
        return fgPower
        

    def makeRealSpaceFilterProfile(self):
        """Makes a 1d real-space profile of the filter. This is normalised to 1 at the frequency of the 
        first map given in the unfiltered maps list.
        
        Returns profile, arcminRange
        
        """
        realSpace=fft.ifft2(self.filt).real
        realSpace=fft.fftshift(realSpace)

        # Changed for python3 - we may want to check this...
        # i.e., force to always be int after /2 without this
        x0=int(realSpace.shape[2]/2)
        y0=int(realSpace.shape[1]/2)

        # Arbitrarily normalise to 1 at the maximum (in whichever frequency that occurs)
        normFactor=abs(realSpace[:, y0, x0:]).max() 
        prof=realSpace[:, y0, x0:]/normFactor
        #prof=realSpace[0, y0, x0:]/realSpace[0, y0, x0:].max() # normalise for plot
        arcminRange=np.arange(0, prof.shape[1])*self.degPerPixX*60.0
                
        return prof, arcminRange
    
        
    def saveRealSpaceFilterProfile(self):
        """Saves a real-space profile of the filter out as a .png plot and 2d .fits image.
        
        """
        
        prof, arcminRange=self.makeRealSpaceFilterProfile()
        
        # Measure characteristic FWHM
        #tck=interpolate.splrep(arcminRange, prof)
        #FWHMArcmin=interpolate.splev([0.5], tck)*2 # *2 because otherwise would be half width
        
        fig=plt.figure(figsize=(8,8))
        ax=plt.axes([0.14, 0.11, 0.85, 0.86])
        #fig.canvas.set_window_title('Filter Profile in Real Space')
        #plt.title("Filter Profile %s" % (self.label))
        plt.ylabel("Amplitude")
        plt.xlabel("$\\theta$ (arcmin)")
        for row, mapDict in zip(prof, self.unfilteredMapsDictList):
            if mapDict['obsFreqGHz'] is not None:
                label = '%d GHz' % (mapDict['obsFreqGHz'])
            elif mapDict['units'] == 'yc':
                label = 'yc'
            plt.plot(arcminRange, row, label = label)
        plt.xlim(0, 30.0)
        plt.ylim(prof.min(), prof.max()*1.1)
        plt.legend()
        plt.savefig(self.diagnosticsDir+os.path.sep+"realSpaceProfile1d_"+self.label+"#"+self.tileName+".png")
        plt.close()
        
        # Save 2D realspace filter image too
        #astImages.saveFITS(self.diagnosticsDir+os.path.sep+"realSpaceProfile2d_"+self.label+".fits", \
                                #realSpace, wcs)
      

    def makeNoiseMap(self, mapData):
        """Estimate the noise map using local RMS measurements on a grid, over the whole filtered map.
        
        """
        
        # Grid method
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
                elif 'RMSEstimator' in self.params['noiseParams'].keys() and self.params['noiseParams']['RMSEstimator'] == 'percentile':
                    chunkRMS=np.percentile(abs(chunkValues[goodAreaMask]), 68.3)
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
    
    
    def loadFRelWeights(self):
        """Reads frequency weights used for relativistic corrections from filter header.
        
        """
        with pyfits.open(self.filterFileName) as img:
            self.fRelWeights={}
            for i in range(1, 10):
                if 'RW%d_GHZ' % (i) in img[0].header.keys():
                    freqGHz=img[0].header['RW%d_GHZ' % (i)]
                    self.fRelWeights[freqGHz]=img[0].header['RW%d' % (i)]
    
#------------------------------------------------------------------------------------------------------------
class MatchedFilter(MapFilter):
    """Multi-frequency matched filter...
    
    """

    def buildAndApply(self):
        
        fMapsToFilter=[]
        for mapDict in self.unfilteredMapsDictList:
            fMapsToFilter.append(enmap.fft(enmap.apod(mapDict['data'], self.apodPix)))
        fMapsToFilter=np.array(fMapsToFilter)
        
        # NOTE: We've tidied up the config file, so we don't have to feed in surveyMask and psMask like this
        # (see startUp.parseConfig)
        surveyMask=self.unfilteredMapsDictList[0]['surveyMask']
        psMask=self.unfilteredMapsDictList[0]['psMask']
            
        if os.path.exists(self.filterFileName) == False:
            print(">>> Building filter %s#%s ..." % (self.label, self.tileName))
                        
            fMapsForNoise=[]
            for mapDict in self.unfilteredMapsDictList: 
                d=mapDict['data']
                fMapsForNoise.append(enmap.fft(enmap.apod(d, self.apodPix)))
            fMapsForNoise=np.array(fMapsForNoise)
        
            # Smoothing noise here is essential
            kernelSize=(3,3)
            noiseCov=[]
            for i in range(len(self.unfilteredMapsDictList)):
                iMap=self.unfilteredMapsDictList[i]
                row=[]
                for j in range(len(self.unfilteredMapsDictList)):
                    jMap=self.unfilteredMapsDictList[j]
                    if self.params['noiseParams']['method'] == 'dataMap':
                        NP=np.real(fMapsForNoise[i]*fMapsForNoise[j].conj())
                    elif self.params['noiseParams']['method'] == 'max(dataMap,CMB)':
                        NP=np.real(fMapsForNoise[i]*fMapsForNoise[j].conj())
                        NPCMB=self.makeForegroundsPower() # This needs a beam convolution adding
                        NP=np.maximum.reduce([NP, NPCMB])
                    else:
                        raise Exception("Other noise models not yet re-implemented")
                    NP=ndimage.gaussian_filter(NP, kernelSize)
                    #astImages.saveFITS("test_%d_%d_%s.fits" % (i,  j, self.tileName), fft.fftshift(NP), None)
                    row.append(NP)
                noiseCov.append(row)
            del fMapsForNoise
            noiseCov=np.array(noiseCov)
            
            # Signal frequency weighting
            w=[]
            for mapDict in self.unfilteredMapsDictList:
                if mapDict['units'] != 'yc':
                    if self.params['outputUnits'] == 'yc':
                        w.append(signals.fSZ(mapDict['obsFreqGHz']))
                    elif self.params['outputUnits'] == 'uK':
                        # This is where variable spectral weighting for sources would be added...
                        w.append(1.0)
                    else:
                        raise Exception('need to specify "outputUnits" ("yc" or "uK") in filter params')
                else:
                    w.append(1.0)   # For tILe-C: there should only be one map if input units are yc anyway...
            w=np.array(w)
                
            # Make FFTs of unit-normalised signal templates for each band
            signalMapsList=[]
            fSignalsArr=[]
            for mapDict in self.unfilteredMapsDictList:
                signalMap=self.makeSignalTemplateMap(mapDict['beamFileName'])
                fSignal=enmap.fft(signalMap)
                signalMapsList.append(signalMap)
                fSignalsArr.append(fSignal)
            fSignalsArr=np.array(fSignalsArr)
                    
            # Build the filter itself
            self.filt=np.zeros([len(self.unfilteredMapsDictList), self.shape[0], self.shape[1]])
            for y in range(0, self.shape[0]):
                for x in range(0, self.shape[1]):
                    try:
                        self.filt[:, y, x]=np.dot(np.linalg.inv(noiseCov[:, :, y, x]), w*abs(fSignalsArr[:, y, x])) 
                    except:
                        continue
            del fSignalsArr
            del noiseCov
                        
            # Use a map with known input signal to figure out how much it has been rolled off by     
            if self.params['outputUnits'] == 'yc':
                # Normalise such that peak value in filtered map == y0, taking out the effect of the beam
                signalMaps=[]
                fSignalMaps=[]
                y0=2e-4
                for mapDict in self.unfilteredMapsDictList:
                    if mapDict['units'] == 'yc':    # For handling tILe-C maps
                        signalMap=self.makeSignalTemplateMap(mapDict['beamFileName'], amplitude = y0)
                    else:                           # The normal case
                        deltaT0=maps.convertToDeltaT(y0, mapDict['obsFreqGHz'])
                        signalMap=self.makeSignalTemplateMap(mapDict['beamFileName'], 
                                                             amplitude = deltaT0)
                    signalMaps.append(signalMap)
                    fSignal=enmap.fft(signalMap)
                    fSignalMaps.append(fSignal)
                signalMaps=np.array(signalMaps)
                fSignalMaps=np.array(fSignalMaps)        
                filteredSignal=self.applyFilter(fSignalMaps)
                self.signalNorm=y0/filteredSignal.max()

                # For relativistic corrections (see signals module)
                totalSignal=filteredSignal.flatten()[np.argmax(filteredSignal)]
                filteredSignalCube=np.real(enmap.ifft(fSignalMaps*self.filt, normalize = False))
                self.fRelWeights={}
                for filteredSignalPlane, mapDict in zip(filteredSignalCube, self.unfilteredMapsDictList):
                    freqGHz=mapDict['obsFreqGHz']
                    fRelWeight=filteredSignalPlane.flatten()[np.argmax(filteredSignal)]/totalSignal
                    self.fRelWeights[freqGHz]=fRelWeight
                del fSignalMaps
                self.signalNorm=y0/filteredSignal.max()
            elif self.params['outputUnits'] == 'uK':
                if len(self.unfilteredMapsDictList) > 1:
                    raise Exception("multi-frequency filtering not currently supported for outputUnits 'uK' (point source finding)")
                combinedObsFreqGHz=float(list(self.beamSolidAnglesDict.keys())[0])  # Make less clunky...
                signalMaps=[]
                fSignalMaps=[]
                for mapDict in self.unfilteredMapsDictList:
                    signalMap=self.makeSignalTemplateMap(mapDict['beamFileName'])
                    signalMaps.append(signalMap)
                    fSignal=enmap.fft(signalMap)
                    fSignalMaps.append(fSignal)
                signalMaps=np.array(signalMaps)
                fSignalMaps=np.array(fSignalMaps)        
                filteredSignal=self.applyFilter(fSignalMaps)
                self.signalNorm=1.0/filteredSignal.max()
                del fSignalMaps
            else:
                raise Exception('need to specify "outputUnits" ("yc" or "uK") in filter params')
        else:
            print(">>> Loading cached filter %s#%s ..." % (self.label, self.tileName))
            self.loadFilter()
        
        # Apply filter
        filteredMap=self.applyFilter(fMapsToFilter)
        del fMapsToFilter
                
        # Units etc.
        if self.params['outputUnits'] == 'yc':
            mapUnits='yc'
            combinedObsFreqGHz='yc'
            beamSolidAngle_nsr=0.0   # not used for clusters...
        elif self.params['outputUnits'] == 'uK':
            if len(self.unfilteredMapsDictList) > 1:
                raise Exception("multi-frequency filtering not currently supported for outputUnits 'uK' (point source finding)")
            combinedObsFreqGHz=float(list(self.beamSolidAnglesDict.keys())[0])  # Make less clunky...
            mapUnits='uK'
            beamSolidAngle_nsr=self.beamSolidAnglesDict[combinedObsFreqGHz]
        else:
            raise Exception('need to specify "outputUnits" ("yc" or "uK") in filter params')
        
        # Apply the point source mask here (before noise estimates etc.)        
        filteredMap=filteredMap*psMask
                
        # Make noise and S/N maps
        RMSMap=self.makeNoiseMap(filteredMap)
        validMask=np.greater(RMSMap, 0)
        SNMap=np.zeros(filteredMap.shape)+filteredMap
        SNMap[validMask]=SNMap[validMask]/RMSMap[validMask]
        
        # If we did a good job of subtracting / filling mask holes, we could apply mask after noise estimates
        #filteredMap=filteredMap*psMask
        #SNMap=SNMap*psMask
        
        # Use rank filter to zap edges where RMS will be artificially low - we use a bit of a buffer here
        # NOTE: Now point source mask is applied above, we fill the holes back in here when finding edges
        if 'edgeTrimArcmin' in self.params.keys() and self.params['edgeTrimArcmin'] > 0:
            trimSizePix=int(round((self.params['edgeTrimArcmin']/60.)/self.wcs.getPixelSizeDeg()))
        else:
            gridSize=int(round((self.params['noiseParams']['noiseGridArcmin']/60.)/self.wcs.getPixelSizeDeg()))
            trimSizePix=int(round(gridSize*3.0))
        edgeCheck=ndimage.rank_filter(abs(filteredMap+(1-psMask)), 0, size = (trimSizePix, trimSizePix))
        edgeCheck=np.array(np.greater(edgeCheck, 0), dtype = float)
        filteredMap=filteredMap*edgeCheck
        apodMask=np.not_equal(filteredMap, 0)
        surveyMask=edgeCheck*surveyMask*psMask
        del edgeCheck

        # Apply final survey mask to signal-to-noise map and RMS map
        # NOTE: need to avoid NaNs in here, otherwise map interpolation for e.g. S/N will fail later on
        SNMap=SNMap*surveyMask
        SNMap[np.isnan(SNMap)]=0.
        RMSMap=RMSMap*surveyMask

        maskFileName=self.selFnDir+os.path.sep+"areaMask#%s.fits.gz" % (self.tileName)
        if os.path.exists(maskFileName) == False:
            astImages.saveFITS(maskFileName, np.array(surveyMask, dtype = int), self.wcs)
        
        if 'saveRMSMap' in self.params and self.params['saveRMSMap'] == True:
            RMSFileName=self.selFnDir+os.path.sep+"RMSMap_%s#%s.fits.gz" % (self.label, self.tileName)
            astImages.saveFITS(RMSFileName, RMSMap, self.wcs)

        try:
            self.saveRealSpaceFilterProfile()   
        except:
            raise Exception("Error saving real space filter profile for filter '%s', tile '%s'" % (self.label, self.tileName))
        
        if 'saveFilter' in self.params and self.params['saveFilter'] == True:
            img=pyfits.PrimaryHDU()                                                                                                                                                              
            img.header['SIGNORM']=self.signalNorm
            count=0
            for key in list(self.fRelWeights.keys()):
                count=count+1
                img.header['RW%d_GHZ' % (count)]=key
                img.header['RW%d' % (count)]=self.fRelWeights[key]
            img.data=self.filt                                                                                                                                                                      
            img.writeto(self.filterFileName, overwrite = True) 
            
        # NOTE: What to do about frequency here? Generalise for non-SZ
        return {'data': filteredMap, 'wcs': self.wcs, 'obsFreqGHz': combinedObsFreqGHz, 'SNMap': SNMap, 
                'mapUnits': mapUnits, 'beamSolidAngle_nsr': beamSolidAngle_nsr}


    def loadFilter(self):
        """Loads in a previously saved filter.
        
        """
        with pyfits.open(self.filterFileName) as img:
            self.filt=img[0].data
            self.signalNorm=img[0].header['SIGNORM']
        self.loadFRelWeights()
        

    def reshapeFilter(self, shape):
        """Use linear interpolation to transform the filter to the given shape.
        
        Returns:
            Reshaped filter (2d numpy array)
        
        """
        
        # If we feed in a 2d shape, make sure we add an axis to keep generalised to multi-frequency
        if len(shape) == 2:
            shape=[self.filt.shape[0], shape[0], shape[1]]
        assert(len(shape) == 3)

        lx, ly=enmap.laxes(self.unfilteredMapsDictList[0]['data'].shape, self.enwcs)
        lxToX=interpolate.interp1d(lx, np.arange(lx.shape[0]), fill_value = 'extrapolate')
        lyToY=interpolate.interp1d(ly, np.arange(ly.shape[0]), fill_value = 'extrapolate')
        lxOut, lyOut=enmap.laxes([shape[1], shape[2]], self.enwcs)
        xOut=lxToX(lxOut)  
        yOut=lyToY(lyOut)
        reshapedFilt=np.zeros(shape)
        for i in range(self.filt.shape[0]):
            filtInterp=interpolate.interp2d(np.arange(ly.shape[0]), np.arange(lx.shape[0]), self.filt[i])
            reshapedFilt[i]=filtInterp(yOut, xOut)
        
        return reshapedFilt
    
            
    def applyFilter(self, mapDataToFilter):
        """Apply the filter to the given map data (must be a cube - with each plane corresponding to a 
        frequency). If the map data is not complex, it will be Fourier transformed. If the map data 
        is not the same shape as the filter, the filter will be interpolated to match.
        
        An optional additional high-pass filter can be applied if 'bckSub' and 'bckSubScaleArcmin' are
        given in self.params.
                
        Returns:
            Filtered map (2d numpy array)
        
        """
        
        # NOTE: need to check appropriate signalNorm after reshaping
        if mapDataToFilter.shape == self.filt.shape:
            filt=self.filt
        else:
            filt=self.reshapeFilter(mapDataToFilter.shape)
        
        if np.any(np.iscomplex(mapDataToFilter)) == True:
            fMapsToFilter=mapDataToFilter
        else:
            fMapsToFilter=enmap.fft(enmap.apod(mapDataToFilter, self.apodPix))

        filteredMap=np.real(enmap.ifft(fMapsToFilter*filt, normalize = False)).sum(axis = 0)

        # Optional additional high-pass filter
        if 'bckSub' in self.params.keys() and 'bckSubScaleArcmin' in self.params.keys() and self.params['bckSub'] == True:
            filteredMap=maps.subtractBackground(filteredMap, self.wcs, smoothScaleDeg = self.params['bckSubScaleArcmin']/60.)
        
        filteredMap=filteredMap*self.signalNorm
        
        return filteredMap
        
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
    
    A map of the final area searched for clusters called 'areaMask.fits.gz' is written in the selFn/ 
    folder.
        
    """       
    
    def loadFilter(self):
        """Loads a previously cached kernel.
        
        Returns kern2d, signalNorm, bckSubScaleArcmin
        
        """        
        with pyfits.open(self.filterFileName) as img:
            kern2d=img[0].data
            signalNorm=img[0].header['SIGNORM']
            if 'BCKSCALE' in img[0].header.keys():
                bckSubScaleArcmin=img[0].header['BCKSCALE']
            else:
                bckSubScaleArcmin=0
            self.applyRACentre=img[0].header['APP_RA']
            self.applyDecCentre=img[0].header['APP_DEC']
        self.kern2d=kern2d
        self.signalNorm
        self.bckSubScaleArcmin=bckSubScaleArcmin

        
    def buildKernel(self, RADecSection, RADeg = 'centre', decDeg = 'centre'):
        """Builds the real space kernel itself. 
        
        RADeg, decDeg are used for figuring out pixel scales for background subtraction
        
        Returns kern2d, signalNorm, bckSubScaleArcmin
        
        """
        
        if os.path.exists(self.filterFileName) == True:
            return self.loadFilter()
        
        wcs=self.wcs

        # Build the matched-filter kernel in a small section of the map
        # Apply the same difference of Gaussians high pass filter here
        # NOTE: we could merge 'bckSubScaleArcmin' and 'maxArcmin' keys here!
        #mapDict['bckSubScaleArcmin']=maxArcmin
        keysWanted=['mapFileName', 'weightsFileName', 'obsFreqGHz', 'units', 'beamFileName', 'addNoise', 
                    'pointSourceRemoval', 'weightsType']
        kernelUnfilteredMapsDictList=[]
        for mapDict in self.unfilteredMapsDictList:
            kernelUnfilteredMapsDict={}
            for k in keysWanted:
                if k in list(mapDict.keys()):
                    kernelUnfilteredMapsDict[k]=mapDict[k]
            kernelUnfilteredMapsDict['RADecSection']=RADecSection
            kernelUnfilteredMapsDictList.append(kernelUnfilteredMapsDict)
        kernelLabel="realSpaceKernel_%s" % (self.label)
        matchedFilterDir=self.diagnosticsDir+os.path.sep+kernelLabel+"#"+self.tileName
        diagnosticsDir=matchedFilterDir+os.path.sep+'diagnostics'
        selFnDir=matchedFilterDir+os.path.sep+'selFn'
        for d in [matchedFilterDir, diagnosticsDir, selFnDir]:
            if os.path.exists(d) == False:
                os.makedirs(d)
        matchedFilterClass=eval(self.params['noiseParams']['matchedFilterClass'])
        matchedFilter=matchedFilterClass(kernelLabel, kernelUnfilteredMapsDictList, self.params, 
                                         tileName = mapDict['tileName'],
                                         diagnosticsDir = matchedFilterDir+os.path.sep+'diagnostics',
                                         selFnDir = matchedFilterDir+os.path.sep+'selFn')
        filteredMapDict=matchedFilter.buildAndApply()   
        
        # Turn the matched filter into a smaller real space convolution kernel
        # This means we have to roll off the kernel to 0 at some radius
        # This is set by maxArcmin in the config file
        kernelMaxArcmin=self.params['noiseParams']['kernelMaxArcmin']
        prof, arcminRange=matchedFilter.makeRealSpaceFilterProfile()
        rIndex=np.where(arcminRange > kernelMaxArcmin)[0][0]
        mask=np.less(arcminRange, kernelMaxArcmin)

        # Kernel can be either fully 2d, or be azimuthally averaged... in the ACTPol E-D56 paper, we used the latter
        if 'symmetrize' in self.params['noiseParams'].keys() and self.params['noiseParams']['symmetrize'] == True:
            rRadians=np.radians(arcminRange/60.)
            profile2d=[]
            for i in range(prof.shape[0]):
                r2p=interpolate.interp1d(rRadians[mask], prof[i, mask], bounds_error=False, fill_value=0.0)
                profile2d.append(r2p(matchedFilter.radiansMap))
            profile2d=np.array(profile2d)
        else:
            profile2d=fft.ifft2(matchedFilter.filt).real
            profile2d=fft.fftshift(profile2d)
            
        # z is not needed here - just because we switched to multi-freq throughout
        z, y, x=np.where(abs(profile2d) == abs(profile2d).max()) 
        #y, x=np.where(profile2d[0] == profile2d[0].max())
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
        self.kern2d=profile2d[:, yMin:yMax, xMin:xMax]
        kern2dRadiansMap=matchedFilter.radiansMap[yMin:yMax, xMin:xMax]

        # This is what to high pass filter on
        if 'bckSubScaleArcmin' in self.params.keys():
            self.bckSubScaleArcmin=self.params['bckSubScaleArcmin']
        else:
            # This is here just so we can reproduce old results where this was done automatically
            # Now we have to fiddle about a bit to check the sign of the filter 
            # (depends on spectral response)
            # We're only considering the first frequency given in a set of multi-freq maps
            if np.greater(prof[0, 0], 0) == True:
                func=np.min
            else:
                func=np.max
            self.bckSubScaleArcmin=arcminRange[prof[0] == func(prof[0])][0]
            
        # Use a map with known input signal to figure out how much it has been rolled off by:
        # 1. The high pass filter (bck sub step)
        # 2. The matched filter itself (includes beam)       
        # NOTE: This is SZ-specific again and needs generalising
        signalMaps=[]
        for mapDict in self.unfilteredMapsDictList:
            # This has been made more complicated because of tILe-C
            if self.params['outputUnits'] == 'yc':
                y0=2e-4
                if mapDict['obsFreqGHz'] is not None:   
                    # Normal case
                    deltaT0=maps.convertToDeltaT(y0, mapDict['obsFreqGHz'])
                    signalMap=self.makeSignalTemplateMap(mapDict['beamFileName'],
                                                         amplitude = deltaT0)
                else:
                    # tILe-C case
                    signalMap=self.makeSignalTemplateMap(mapDict['beamFileName'], amplitude = y0)
            elif self.params['outputUnits'] == 'uK':
                signalMap=self.makeSignalTemplateMap(mapDict['beamFileName'])
            else:
                raise Exception('need to specify "outputUnits" ("yc" or "uK") in filter params')
            signalMaps.append(signalMap)
        signalMaps=np.array(signalMaps)
        
        filteredSignal=self.applyFilter(signalMaps, calcFRelWeights = True)
        if self.params['outputUnits'] == 'yc':
            # Normalise such that peak value in filtered map == y0, taking out the effect of the beam
            self.signalNorm=y0/filteredSignal.max()            
        elif self.params['outputUnits'] == 'uK':
            self.signalNorm=1.0/filteredSignal.max()
        else:
            raise Exception('need to specify "outputUnits" ("yc" or "uK") in filter params')
        
        # Save 2d kernel - we need this (at least for the photometry ref scale) to calc Q later
        # Add bckSubScaleArcmin to the header
        kernWCS=wcs.copy()
        if self.params['bckSub'] == True:
            kernWCS.header['BCKSCALE']=self.bckSubScaleArcmin
        kernWCS.header['SIGNORM']=self.signalNorm
        kernWCS.header['APP_RA']=self.applyRACentre
        kernWCS.header['APP_DEC']=self.applyDecCentre
        count=0
        for key in list(self.fRelWeights.keys()):
            count=count+1
            kernWCS.header['RW%d_GHZ' % (count)]=key
            kernWCS.header['RW%d' % (count)]=self.fRelWeights[key]
        astImages.saveFITS(self.filterFileName, self.kern2d, kernWCS)
        
        # Filter profile plot   
        # Save the stuff we plot first, in case we want to make a plot with multiple filters on later
        np.savez(self.diagnosticsDir+os.path.sep+"filterProf1D_%s#%s.npz" % (self.label, self.tileName), 
                 arcminRange = arcminRange, prof = prof, mask = mask, 
                 bckSubScaleArcmin = self.bckSubScaleArcmin)
        plotSettings.update_rcParams()
        plt.figure(figsize=(9,6.5))
        ax=plt.axes([0.13, 0.12, 0.86, 0.86])
        #plt.tick_params(axis='both', which='major', labelsize=15)
        #plt.tick_params(axis='both', which='minor', labelsize=15)
        for row, mapDict in zip(prof, self.unfilteredMapsDictList):
            tck=interpolate.splrep(arcminRange[mask], row[mask])
            plotRange=np.linspace(0, arcminRange[mask].max(), 1000)
            #plt.plot(arcminRange[mask], prof[mask])
            if mapDict['obsFreqGHz'] is not None:
                label = '%d GHz' % (mapDict['obsFreqGHz'])
            elif mapDict['units'] == 'yc':
                label = 'yc'
            plt.plot(plotRange, interpolate.splev(plotRange, tck), '-', label = label)
        plt.xlabel("$\\theta$ (arcmin)")
        plt.ylabel("Amplitude")
        plt.legend()
        #plt.title(self.label)
        #plt.plot(arcminRange[mask], [0]*len(arcminRange[mask]), 'k--')
        plt.xlim(0, arcminRange[mask].max())
        if self.params['bckSub'] == True:
            plt.plot([self.bckSubScaleArcmin]*3, np.linspace(-1.2, 1.2, 3), 'k--')
        plt.ylim(-1.2, 0.2)
        plt.savefig(self.diagnosticsDir+os.path.sep+"filterPlot1D_%s#%s.pdf" % (self.label, self.tileName))
        plt.close()

            
    def buildAndApply(self):

        surveyMask=self.unfilteredMapsDictList[0]['surveyMask']
        psMask=self.unfilteredMapsDictList[0]['psMask']
            
        if os.path.exists(self.filterFileName) == False:
            print(">>> Building filter %s#%s ..." % (self.label, self.tileName))
            
            # Noise region to use
            RAMin, RAMax, decMin, decMax=self.wcs.getImageMinMaxWCSCoords()
            if self.params['noiseParams']['RADecSection'] == 'tileNoiseRegions':
                RADecSection=[self.wcs.header['NRAMIN'], self.wcs.header['NRAMAX'], 
                            self.wcs.header['NDEMIN'], self.wcs.header['NDEMAX']]
                print("... taking noise region from tileDeck image header: %s ..." % (RADecSection))
            else:
                RADecSection=self.params['noiseParams']['RADecSection']
            self.applyDecCentre=(decMax+decMin)/2.
            self.applyRACentre=(RAMax+RAMin)/2.
            
            # Build kernel   
            self.buildKernel(RADecSection, RADeg = self.applyRACentre, decDeg = self.applyDecCentre)
        else:
            print(">>> Loading cached filter %s#%s ..." % (self.label, self.tileName))
            self.loadFilter()

        # Apply kernel        
        t0=time.time()
        mapDataToFilter=[]
        for mapDict in self.unfilteredMapsDictList:
            mapDataToFilter.append(mapDict['data'])
        mapDataToFilter=np.array(mapDataToFilter)
        print("... convolving map with kernel ...")
        filteredMap=self.applyFilter(mapDataToFilter)
        t1=time.time()
        print("... took %.3f sec ..." % (t1-t0))
                                                
        # Apply the point source mask here (before noise estimates etc.)
        filteredMap=filteredMap*psMask
        
        # Make noise and S/N maps
        RMSMap=self.makeNoiseMap(filteredMap)
        validMask=np.greater(RMSMap, 0)
        SNMap=np.zeros(filteredMap.shape)+filteredMap
        SNMap[validMask]=SNMap[validMask]/RMSMap[validMask]
        
        # Units etc.
        if self.params['outputUnits'] == 'yc':
            mapUnits='yc'
            combinedObsFreqGHz='yc'
            beamSolidAngle_nsr=0.0   # not used for clusters...
        elif self.params['outputUnits'] == 'uK':
            if len(self.unfilteredMapsDictList) > 1:
                raise Exception("multi-frequency filtering not currently supported for outputUnits 'uK' (point source finding)")
            combinedObsFreqGHz=float(list(self.beamSolidAnglesDict.keys())[0])  # Make less clunky...
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
        if trimSizePix > 0:
            edgeCheck=ndimage.rank_filter(abs(filteredMap+(1-psMask)), 0, size = (trimSizePix, trimSizePix))
        else:
            edgeCheck=np.ones(filteredMap.shape)
        edgeCheck=np.array(np.greater(edgeCheck, 0), dtype = float)
        filteredMap=filteredMap*edgeCheck
        apodMask=np.not_equal(filteredMap, 0)
        surveyMask=edgeCheck*surveyMask*psMask
        del edgeCheck

        # Apply final survey mask to signal-to-noise map and RMS map
        # NOTE: need to avoid NaNs in here, otherwise map interpolation for e.g. S/N will fail later on
        SNMap=SNMap*surveyMask
        SNMap[np.isnan(SNMap)]=0.
        RMSMap=RMSMap*surveyMask

        maskFileName=self.selFnDir+os.path.sep+"areaMask#%s.fits.gz" % (self.tileName)
        if os.path.exists(maskFileName) == False:
            astImages.saveFITS(maskFileName, np.array(surveyMask, dtype = int), self.wcs)
        
        if 'saveRMSMap' in self.params and self.params['saveRMSMap'] == True:
            RMSFileName=self.selFnDir+os.path.sep+"RMSMap_%s#%s.fits.gz" % (self.label, self.tileName)
            astImages.saveFITS(RMSFileName, RMSMap, self.wcs)

        return {'data': filteredMap, 'wcs': self.wcs, 'obsFreqGHz': combinedObsFreqGHz,
                'SNMap': SNMap, 'mapUnits': mapUnits, 'beamSolidAngle_nsr': beamSolidAngle_nsr}

    
    def applyFilter(self, mapDataToFilter, calcFRelWeights = False):
        """Apply the kernel to the given map data (must be a cube - with each plane corresponding to a 
        frequency).
        
        NOTE: calcFRelWeights = True should ONLY be set if this routine is being applied to an ideal
        signal map (when calculating signalNorm). This operation is being done in here to save 
        complications / time from doing unnecessary convolution operations elsewhere.
        
        Returns:
            Filtered map (2d numpy array)
        
        """
        
        # Apply the high pass filter - subtract background on larger scales using difference of Gaussians 
        filteredMap=np.zeros(mapDataToFilter.shape)
        if self.params['bckSub'] == True and self.bckSubScaleArcmin > 0:
            for i in range(mapDataToFilter.shape[0]):
                filteredMap[i]=maps.subtractBackground(mapDataToFilter[i], self.wcs, 
                                                       RADeg = self.applyRACentre,
                                                       decDeg = self.applyDecCentre,
                                                       smoothScaleDeg = self.bckSubScaleArcmin/60.)
        else:
            filteredMap=filteredMap+mapDataToFilter
        
        # Apply the kernel
        for i in range(filteredMap.shape[0]):
            filteredMap[i]=ndimage.convolve(filteredMap[i], self.kern2d[i])
        
        # For relativistic corrections (see signals module)
        if calcFRelWeights == True:
            self.fRelWeights={}
            maxIndex=np.argmax(filteredMap.sum(axis = 0))
            totalSignal=filteredMap.sum(axis = 0).flatten()[maxIndex]
            for filteredSignalPlane, mapDict in zip(filteredMap, self.unfilteredMapsDictList):
                freqGHz=mapDict['obsFreqGHz']
                fRelWeight=filteredSignalPlane.flatten()[maxIndex]/totalSignal
                self.fRelWeights[freqGHz]=fRelWeight
                
        filteredMap=filteredMap.sum(axis = 0)
        
        # Apply the normalisation
        filteredMap=filteredMap*self.signalNorm
        
        return filteredMap
    
        
#------------------------------------------------------------------------------------------------------------
class BeamFilter(MapFilter):
    """Base class for filters using beam profile files in Matthew + Kavi's format.
        
    """
    
    def makeSignalTemplateMap(self, beamFileName, amplitude = None):
        """Makes a beam model signal template map.
        
        """
        
        signalMap=signals.makeBeamModelSignalMap(np.degrees(self.radiansMap),
                                                            self.wcs, 
                                                            beamFileName,
                                                            amplitude = amplitude)

        return signalMap
    
#------------------------------------------------------------------------------------------------------------
class ArnaudModelFilter(MapFilter):
    """Base class for filters using the GNFW profile as described in Arnaud et al. (2010).
    
    """
    
    def makeSignalTemplateMap(self, beamFileName, amplitude = None):
        """Makes a model signal template map.
                
        """
        
        signalMap=signals.makeArnaudModelSignalMap(self.params['z'], self.params['M500MSun'], 
                                                   np.degrees(self.radiansMap),
                                                   self.wcs, beamFileName, 
                                                   GNFWParams = self.params['GNFWParams'],
                                                   amplitude = amplitude)
        
        return signalMap
                
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
