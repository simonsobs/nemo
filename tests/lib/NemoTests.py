"""

Library for running nemo tests using Robot Framework

"""

import os
import sys
import subprocess
import shutil
import numpy as np
from pixell import utils, pointsrcs, enmap
from astLib import *
import astropy.io.fits as pyfits
import astropy.table as atpy
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import pylab as plt

#------------------------------------------------------------------------------------------------------------
class NemoTests(object):

    def __init__(self):
        """Set-up for equD56 cluster tests. This downloads the needed maps if they are not found.
        
        """

        self._status = ''
        
        # Needed for inject_sources
        self.beamFileName="../examples/equD56/profiles_ACT/profile_AR1_2009_pixwin_130224.txt"

        cacheDir="testsCache"
        if os.path.exists(cacheDir) == False:
            os.makedirs(cacheDir)
        thisDir=os.getcwd()
        
        self.plotsDir="plots"
        if os.path.exists(self.plotsDir) == False:
            os.makedirs(self.plotsDir)

        # Download the E-D56 map if not found
        self.inMapFileName=cacheDir+os.path.sep+"weightedMap_4.fits"
        if os.path.exists(self.inMapFileName) == False:
            print(">>> Downloading E-D56 map ...")
            os.chdir(cacheDir)
            os.system("wget https://www.acru.ukzn.ac.za/~mjh/equD56Maps.tar.gz")
            os.system("tar -zxvf equD56Maps.tar.gz")
            os.remove("equD56Maps.tar.gz")
            os.chdir(thisDir)


    def set_config(self, configFileName):
        self.configFileName=configFileName
        
        
    def run_nemo(self):
        self._run_command(["nemo", self.configFileName])


    def run_nemo_mass(self, catalogFileName = None):
        args=['nemoMass', self.configFileName]
        if catalogFileName != None:
            args=args+['-c', catalogFileName]
        self._run_command(args)


    def run_nemo_selfn(self):
        self._run_command(["nemoSelFn", self.configFileName])


    def run_nemo_mock(self):
        self._run_command(["nemoMock", self.configFileName])
        

    def inject_sources(self, outMapFileName = "sourceInjectedMap.fits", 
                       catalogFileName = "inputSourcesCatalog.fits", numSources = 1000):
        """Injects a bunch of point sources into a map, using a vaguely plausible amplitude distribution
        based on the sources found in the deep56 map.
        
        """
        
        nsigma=5.0
        smul=1.0

        img=pyfits.open(self.inMapFileName)
        mapData=img[0].data
        if mapData.ndim > 2:
            mapData=mapData[0]
        wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
        imap=enmap.read_map(self.inMapFileName)

        # Random fluxes (in delta T uK)
        # This is chosen to vaguely match that seen in deep56 f090
        I=np.random.lognormal(np.log(600), 1.1, numSources)
        
        # Random positions
        ys, xs=np.where(mapData != 0)
        ys=ys+np.random.uniform(0, 1, len(ys))
        xs=xs+np.random.uniform(0, 1, len(xs))
        indices=np.random.randint(0, len(ys), numSources)
        coords=wcs.pix2wcs(xs[indices], ys[indices])
        coords=np.array(coords)
        
        tab=atpy.Table()
        tab.add_column(atpy.Column(coords[:, 0], "ra"))
        tab.add_column(atpy.Column(coords[:, 1], "dec"))
        tab.add_column(atpy.Column(I, "I"))
        
        if os.path.exists(catalogFileName) == True:
            os.remove(catalogFileName)
        tab.write(catalogFileName)
        srcs=pointsrcs.tab2recarray(tab)

        beam=pointsrcs.read_beam(self.beamFileName)
        beam[0]*=utils.degree
        srcparam=pointsrcs.src2param(srcs)

        omap=pointsrcs.sim_srcs(imap.shape, imap.wcs, srcparam, beam, imap, smul=smul, nsigma=nsigma, pixwin=True)
        enmap.write_map(outMapFileName, omap)
        
        
    def cross_match(self, inCatalogFileName, outCatalogFileName, radiusArcmin = 2.5):
        """Cross matches input and output source catalogs.
        
        """
        inTab=atpy.Table().read(inCatalogFileName)
        outTab=atpy.Table().read(outCatalogFileName)
        RAKey1, decKey1=self.getRADecKeys(inTab)
        RAKey2, decKey2=self.getRADecKeys(outTab)
        cat1=SkyCoord(ra = inTab[RAKey1].data, dec = inTab[decKey1].data, unit = 'deg')
        xMatchRadiusDeg=radiusArcmin/60.
        cat2=SkyCoord(ra = outTab[RAKey2].data, dec = outTab[decKey2].data, unit = 'deg')
        xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
        mask=np.less(rDeg.value, xMatchRadiusDeg)  
        matched_outTab=outTab[xIndices]
        inTab=inTab[mask]
        matched_outTab=matched_outTab[mask]
        rDeg=rDeg.value[mask]
        self.inTab=inTab
        self.outTab=matched_outTab
        self.rDeg=rDeg
    
    
    def getRADecKeys(self, tab):
        """Returns the column names in which RA, dec coords are stored, after trying a couple of variations.
        
        """
        RAKeysToTry=['ra', 'RADeg']
        decKeysToTry=['dec', 'decDeg']
        RAKey, decKey=None, None
        for key in RAKeysToTry:
            if key in tab.keys():
                RAKey=key
                break
        for key in decKeysToTry:
            if key in tab.keys():
                decKey=key
                break
        if RAKey == None or decKey == None:
            raise Exception("couldn't identify RA, dec columns in the supplied table")
        
        return RAKey, decKey
        
        
    def check_recovered_ratio(self, inKey, outKey, tolerance = 0.01, SNRCut = 4,
                              SNRKey = 'fixed_SNR', plotFileName = None):
        """Cross match the input and output catalogs pointed to by the respective file names, and calculate 
        the median of the ratio of the columns pointed to by inKey, outKey. If tolerance (defined as
        1 - median ratio) is exceeded, the test is failed. SNRCut is applied to SNRKey in the output catalog.
        
        """
        inTab=self.inTab
        outTab=self.outTab
        ratio=inTab[inKey]/outTab[outKey]
        medRatio=np.median(ratio[outTab[SNRKey] > SNRCut])
        label="median input %s / output %s = %.3f (%s > %.1f)" % (inKey, outKey, medRatio, SNRKey, SNRCut)
        print("... %s" % (label))            
        if plotFileName != None:
            plt.plot(outTab[SNRKey], ratio, 'r.')
            plt.xlabel(SNRKey)
            plt.ylabel("input %s / output %s" % (inKey, outKey))
            plt.plot(np.linspace(0, 16, 3), [1]*3, 'k--')
            plt.xlim(SNRCut, 16)
            plt.title(label)
            plt.ylim(0.5, 1.5)
            plt.savefig(plotFileName)
        plt.close()
        if abs(1.0-medRatio) > tolerance:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
            
            
    def check_recovered_positions(self, toleranceArcsec = 12.0, SNRKey = 'SNR', SNRMax = 10.0, 
                                  plotFileName = None):
        """Blah
        
        """
        inTab=self.inTab
        outTab=self.outTab
        rDeg=self.rDeg
        SNRs=outTab[SNRKey]
        medOffsetArcmin=np.median(rDeg)*60
        SNRMask=np.less(SNRs, SNRMax)
        medOffsetArcmin_SNR10=np.median(rDeg[SNRMask])*60#[SNRMask])*60
        print('... median recovered position offset = %.2f" (full sample)' % (medOffsetArcmin*60))
        label='median recovered position offset = %.2f" (SNR < 10)' % (medOffsetArcmin_SNR10*60)
        print("... %s" % (label))
        if plotFileName != None:
            plt.plot(SNRs, rDeg*3600., 'r.')
            plt.semilogx()
            plt.xlabel("SNR")
            plt.ylabel('Position offset (")')
            plt.title(label)
            plt.savefig(plotFileName)
            plt.close()
        if medOffsetArcmin_SNR10*60 > toleranceArcsec:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
            
            
    def status_should_be(self, expected_status):
        if expected_status != self._status:
            raise AssertionError("Expected status to be '%s' but was '%s'."
                                 % (expected_status, self._status))


    def _run_command(self, args):
        process = subprocess.run(args, universal_newlines=True, stdout=subprocess.PIPE,
                                 stderr=subprocess.STDOUT)
        #self._status = process.communicate()[0].strip()

