"""

Library for running tests using Robot Framework

"""

import os
import sys
import subprocess
import numpy as np
from pixell import utils, pointsrcs, enmap
from astLib import *
import astropy.io.fits as pyfits
import astropy.table as atpy
from astropy.coordinates import SkyCoord
from astropy.coordinates import match_coordinates_sky
import pylab as plt

#------------------------------------------------------------------------------------------------------------
# Test threholds
FLUX_RATIO_THRESH=0.01          # Require median recovered source amplitude to be within 1% of input
POSREC_SNR10_THRESH=0.2         # Require median recovered source position within 0.2' of input for SNR < 10 

#------------------------------------------------------------------------------------------------------------
class PointSourceTests(object):

    def __init__(self):
        """Set-up for point source recovery tests. This downloads the needed maps if they are not found.
        
        """

        self._status = ''
        
        self.beamFileName="../examples/equD56/profiles_ACT/profile_AR1_2009_pixwin_130224.txt"

        cacheDir="testsCache"
        if os.path.exists(cacheDir) == False:
            os.makedirs(cacheDir)
        thisDir=os.getcwd()

        # Download the E-D56 map if not found
        self.inMapFileName=cacheDir+os.path.sep+"weightedMap_4.fits"
        if os.path.exists(self.inMapFileName) == False:
            print(">>> Downloading E-D56 map ...")
            os.chdir(cacheDir)
            os.system("wget https://www.acru.ukzn.ac.za/~mjh/equD56Maps.tar.gz")
            os.system("tar -zxvf equD56Maps.tar.gz")
            os.remove("equD56Maps.tar.gz")
            os.chdir(thisDir)

        self.outMapFileName="sourceInjectedMap.fits"            
        self.outCatalogFileName="inputSourcesCatalog.fits"
        self.configFileName="configs"+os.path.sep+"PSTest_E-D56.yml"

    
    def inject_sources(self, numSources = 1000):
        
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
        
        if os.path.exists(self.outCatalogFileName) == True:
            os.remove(self.outCatalogFileName)
        tab.write(self.outCatalogFileName)
        srcs=pointsrcs.tab2recarray(tab)

        beam=pointsrcs.read_beam(self.beamFileName)
        beam[0]*=utils.degree
        srcparam=pointsrcs.src2param(srcs)

        omap=pointsrcs.sim_srcs(imap.shape, imap.wcs, srcparam, beam, imap, smul=smul, nsigma=nsigma, pixwin=True)
        enmap.write_map(self.outMapFileName, omap)
        
        
    def run_nemo(self):
        os.system("nemo %s" % (self.configFileName))

    
    def cross_match_with_input(self):
        radiusArcmin=2.5    # Prior for filtering out implausible matches
        inTab=atpy.Table().read(self.outCatalogFileName)
        outTab=atpy.Table().read("configs"+os.path.sep+"PSTest_E-D56"+os.path.sep+"PSTest_E-D56_optimalCatalog.fits")
        cat1=SkyCoord(ra = inTab['ra'].data, dec = inTab['dec'].data, unit = 'deg')
        xMatchRadiusDeg=radiusArcmin/60.
        cat2=SkyCoord(ra = outTab['RADeg'].data, dec = outTab['decDeg'].data, unit = 'deg')
        xIndices, rDeg, sep3d = match_coordinates_sky(cat1, cat2, nthneighbor = 1)
        mask=np.less(rDeg.value, xMatchRadiusDeg)        
        return inTab, outTab, xIndices, mask, rDeg
        
        
    def check_recovered_amplitudes(self):
        # Source amplitude recovery check
        inTab, outTab, xIndices, mask, rDeg=self.cross_match_with_input()
        print(">>> Amplitude recovery test:")
        inT=inTab['I'][mask]
        diffT=inTab['I'][mask]-outTab[xIndices]['deltaT_c'][mask]
        ratioT=inTab['I'][mask]/outTab[xIndices]['deltaT_c'][mask]
        SNRs=outTab[xIndices]['SNR'][mask]
        print("... median amplitude in / amplitude out = %.6f" % (np.median(ratioT)))
        plt.plot(inT, ratioT, 'r.')
        plt.xlabel("input deltaT (uK)")
        plt.ylabel("input / output deltaT")
        plt.savefig("amplitudeRecovery.png")
        plt.close()
        if abs(1.0-np.median(ratioT)) > FLUX_RATIO_THRESH:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
            
            
    def check_recovered_positions(self):
        # Position recovery check
        inTab, outTab, xIndices, mask, rDeg=self.cross_match_with_input()
        print(">>> Position recovery test:")
        SNRs=outTab[xIndices]['SNR'][mask]
        medOffsetArcmin=np.median(rDeg[mask].data)*60
        SNRMask=np.less(SNRs, 10.0)
        medOffsetArcmin_SNR10=np.median(rDeg[mask][SNRMask].data)*60
        print("... median recovered position offset = %.6f arcmin (full sample)" % (medOffsetArcmin))
        print("... median recovered position offset = %.6f arcmin (SNR < 10)" % (medOffsetArcmin_SNR10))
        plt.plot(SNRs, rDeg[mask]*60., 'r.')
        plt.semilogx()
        plt.xlabel("SNR")
        plt.ylabel("Position offset (arcmin)")
        plt.savefig("positionRecovery.png")
        plt.close()
        if medOffsetArcmin_SNR10 > POSREC_SNR10_THRESH:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
            
            
    def status_should_be(self, expected_status):
        if expected_status != self._status:
            raise AssertionError("Expected status to be '%s' but was '%s'."
                                 % (expected_status, self._status))


    def _run_command(self, command, *args):
        command = [sys.executable, self._sut_path, command] + list(args)
        process = subprocess.Popen(command, universal_newlines=True, stdout=subprocess.PIPE,
                                   stderr=subprocess.STDOUT)
        self._status = process.communicate()[0].strip()
