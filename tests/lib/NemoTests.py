"""

Library for running nemo tests using Robot Framework

"""

import os
import sys
import subprocess
import shutil
import numpy as np
from nemo import catalogs, maps, plotSettings
from astLib import *
import astropy.io.fits as pyfits
import astropy.table as atpy
from astropy.coordinates import SkyCoord, match_coordinates_sky
import pylab as plt

plotSettings.update_rcParams()
plotTitleSize=14

#------------------------------------------------------------------------------------------------------------
class NemoTests(object):

    def __init__(self):
        """Basic set-up for running tests. For those that need datasets to be downloaded, a setup_<name> 
        routine should be added (see, e.g., setup_equD56).
        
        """

        self._status = ''

        self.cacheDir="testsCache"
        if os.path.exists(self.cacheDir) == False:
            os.makedirs(self.cacheDir)
        
        # Why don't we just run all tests from the cache dir?
        self.runDir=os.path.abspath(self.cacheDir)
            
        self.plotsDir="plots"
        if os.path.exists(self.plotsDir) == False:
            os.makedirs(self.plotsDir)
    
    
    def setup_equD56(self):
        """Set-up for tests that use E-D56 maps - downloads them if not found.
        
        """
    
        # Needed for inject_sources
        self.beamFileName=self.cacheDir+os.path.sep+"profiles_ACT/profile_AR1_2009_pixwin_130224.txt"
        
        thisDir=os.getcwd()
        
        # Download the E-D56 map if not found
        self.inMapFileName=self.cacheDir+os.path.sep+"weightedMap_4.fits"
        if os.path.exists(self.inMapFileName) == False:
            print(">>> Downloading E-D56 map ...")
            os.chdir(self.cacheDir)
            os.system("wget https://astro.ukzn.ac.za/~mjh/equD56Maps.tar.gz")
            os.system("tar -zxvf equD56Maps.tar.gz")
            os.remove("equD56Maps.tar.gz")
            os.chdir(thisDir)
        
        # This one is actually in the distribution as it stands but anyway...
        if os.path.exists(self.cacheDir+os.path.sep+"E-D56Clusters.fits") == False:
            os.chdir(self.cacheDir)
            os.system("wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/actpol_2016_lensing/E-D56Clusters.fits")
            os.chdir(thisDir)
            
        # Copy files that are in the git/source distribution (saves re-organising)
        if os.path.exists(self.cacheDir+os.path.sep+"surveyMask.fits.gz") == False:
            os.system("cp %s %s/" % ("../examples/ACT-DR3-clusters/surveyMask.fits.gz", self.cacheDir))
        if os.path.exists(self.cacheDir+os.path.sep+"pointSourceMask.fits.gz") == False:
            os.system("cp %s %s/" % ("../examples/ACT-DR3-clusters/pointSourceMask.fits.gz", self.cacheDir))
        if os.path.exists(self.cacheDir+os.path.sep+"ACTPol_redshifts.fits") == False:
            os.system("cp %s %s/" % ("../examples/ACT-DR3-clusters/ACTPol_redshifts.fits", self.cacheDir))
        if os.path.exists(self.cacheDir+os.path.sep+"profiles_ACT") == False:
            os.system("cp -R %s %s/" % ("../examples/ACT-DR3-clusters/profiles_ACT", self.cacheDir))
            
        
    def setup_south2008(self):
        """Set-up for tests that use southern 2008 ACT data - downloads needed files from LAMBDA if not 
        found.
        
        """
        
        thisDir=os.getcwd()
        self.inMapFileName=self.cacheDir+os.path.sep+"ACT_148_south_season_2_1way_v3_summed.fits"
        if os.path.exists(self.inMapFileName) == False:
            print(">>> Downloading South 2008 data ...")
            os.chdir(self.cacheDir)
            os.system("wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/data2013/Maps/AR1/South/ACT_148_south_season_2_1way_v3_summed.fits")
            os.system("wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/data2013/Weights/AR1/South/ACT_148_south_season_2_1way_hits_v3.fits")
            os.system("wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/data2013/Beams/profiles/profile_AR1_2008_pixwin_130224.txt")
            os.system("wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/Cluster_src/Ptsrc_cat/act_source_catalog_AR1_2008.txt")
            os.chdir(thisDir)
            
        # Need to convert published catalog such that comparison routines work
        # NOTE: loading the table this way breaks the name column but we soldier on...
        tabFileName=self.cacheDir+os.path.sep+"act_source_catalog_AR1_2008.fits"
        if os.path.exists(tabFileName) == False:
            tab=atpy.Table().read(self.cacheDir+os.path.sep+"act_source_catalog_AR1_2008.txt", format = 'ascii')
            tab.rename_column("col2", "name")
            tab.rename_column("col3", "RADeg")
            tab.rename_column("col4", "decDeg")
            tab.rename_column("col6", "fluxJy")
            tab['fluxJy']=tab['fluxJy']/1000.0
            tab.write(tabFileName)


    def setup_quickstart(self):
        """Set-up for tests that use the files from the quickstart tutorial - downloads them if not found.
        
        """
    
        # Maps and beams
        thisDir=os.getcwd()        
        self.inMapFileName=self.cacheDir+os.path.sep+"maps"+os.path.sep+"f150_1_10_8_map.fits"
        if os.path.exists(self.inMapFileName) == False:
            print(">>> Downloading quickstart files ...")
            os.chdir(self.cacheDir)
            os.system("wget 'https://www.dropbox.com/scl/fi/5sjapfshk8g3kxgy4zc56/nemo-quickstart-maps.tar.gz?rlkey=swpdttpgvkffbgr9pcmebslw2&dl=0' -O nemo-quickstart-maps.tar.gz")
            os.system("tar -zxvf nemo-quickstart-maps.tar.gz")
            os.remove("nemo-quickstart-maps.tar.gz")
            os.chdir(thisDir)
        
        # DR5 release catalog for comparing signals
        comparisonCatalog=self.cacheDir+os.path.sep+"DR5_cluster-catalog_v1.1.fits"
        if os.path.exists(comparisonCatalog) == False:
            os.chdir(self.cacheDir)
            os.system("wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/DR5_cluster-catalog_v1.1.fits")
            os.chdir(thisDir)
        
        # Generate mask file(s) for simming
        headerFileNames=["configs/smallTestSurveyMaskHeader.txt", "configs/MPITestSurveyMaskHeader.txt"]
        maskFileNames=["smallTestSurveyMask.fits", "MPITestSurveyMask.fits"]
        for headerFileName, maskFileName in zip(headerFileNames, maskFileNames):
            header=pyfits.Header().fromtextfile(headerFileName)
            wcs=astWCS.WCS(header, mode = 'pyfits')
            d=np.ones([wcs.header['NAXIS2'], wcs.header['NAXIS1']], dtype = int)
            maps.saveFITS(self.cacheDir+os.path.sep+maskFileName, d, wcs, compressionType = 'PLIO_1')
        
        # Map/frequency-related defaults
        self.bandsDict={'f150': {'beam': "maps/s16_pa2_f150_nohwp_night_beam_profile_jitter.txt",
                                 'freq': str(149.6)},
                        'f090': {'beam': "maps/s16_pa3_f090_nohwp_night_beam_profile_jitter.txt",
                                 'freq': str(97.8)}
                       }
        self.sizesDict={'small': {'mask': "smallTestSurveyMask.fits"},
                        'large': {'mask': "MPITestSurveyMask.fits"}
                       }
    

    def set_config(self, configFileName):
        """Set the config file to be used by all nemo scripts executed by tests.
        Path given here can be relative to the current working directory (it gets 
        turned into an absolute path if so).
        
        """
        self.configFileName=os.path.abspath(configFileName)
        self.selFnDir=configFileName.replace(".yml", "")+os.path.sep+"selFn"
        self.mocksDir=configFileName.replace(".yml", "")+os.path.sep+"mocks"
        
        
    def run_nemo(self, forcedPhotometryCatalog = None):
        args=['nemo', self.configFileName]
        if forcedPhotometryCatalog is not None:
            args=args+['-f', forcedPhotometryCatalog]
        self._run_command(args)


    def run_nemo_injection_test(self):
        self._run_command(["nemo", self.configFileName, "-I"])


    def run_parallel_nemo(self):
        self._run_command(["mpiexec", "-np", "4", "nemo", self.configFileName, "-M"])


    def run_nemo_mass(self, catalogFileName = None, forcedPhotometry = False):
        args=['nemoMass', self.configFileName]
        if catalogFileName is not None:
            args=args+['-c', catalogFileName]
        if forcedPhotometry == True:
            args=args+['-F']
        self._run_command(args)


    def make_sim(self, catalogFileName = None, noiseLevel = None, seed = None, band = None, size = None):
        args=['nemoModel', catalogFileName, self.sizesDict[size]['mask'], self.bandsDict[band]['beam'], 
              "sim_%s.fits" % (band), '-f', self.bandsDict[band]['freq'], 
              '-C', '-N', noiseLevel, '-S', seed]
        self._run_command(args)


    def make_parallel_sim(self, catalogFileName = None, noiseLevel = None, seed = None, band = None, size = None, profile = None):
        args=['mpiexec', '-np', '4', 'nemoModel', catalogFileName, self.sizesDict[size]['mask'], self.bandsDict[band]['beam'],
              "sim_%s.fits" % (band), '-f', self.bandsDict[band]['freq'],
              '-C', '-N', noiseLevel, '-S', seed, '-M']
        if profile is not None:
            args=args+['-p', profile]
        self._run_command(args)


    def make_signal_free_sim(self, noiseLevel = None, seed = None, band = None, size = None):
        args=['nemoModel', "pointsources-0", self.sizesDict[size]['mask'], self.bandsDict[band]['beam'],
              "signal_free_%s.fits" % (band), '-f', self.bandsDict[band]['freq'],
              '-C', '-N', noiseLevel, '-S', seed]
        self._run_command(args)
        

    def make_signal_only_sim(self, catalogFileName = None, band = None, size = None):
        args=['nemoModel', catalogFileName, self.sizesDict[size]['mask'], self.bandsDict[band]['beam'], 
              "signal_model_only_%s.fits" % (band), '-f', self.bandsDict[band]['freq']]
        self._run_command(args)
        
        
    def make_parallel_signal_only_sim(self, catalogFileName = None, band = None, size = None, profile = None):
        args=['mpiexec', '-np', '4', 'nemoModel', catalogFileName, self.sizesDict[size]['mask'], self.bandsDict[band]['beam'],
              "signal_model_only_%s.fits" % (band), '-f', self.bandsDict[band]['freq'], '-M']
        if profile is not None:
            args=args+['-p', profile]
        self._run_command(args)
    

    def run_nemo_selfn(self, configFileName = None):
        self._run_command(["nemoSelFn", configFileName])


    def run_nemo_mock(self, selFnDir = None, mocksDir = None):
        self._run_command(["nemoMock", selFnDir, mocksDir])
        

    def inject_sources_using_nemo(self, outMapFileName = "sourceInjectedMap.fits", 
                       catalogFileName = "inputSourcesCatalog.fits", numSources = 1000):
        """Injects a bunch of point sources into a map, using a vaguely plausible amplitude distribution
        based on the sources found in the deep56 map.
        
        """
        
        with pyfits.open(self.inMapFileName) as img:
            mapData=img[0].data
            if mapData.ndim > 2:
                mapData=mapData[0]
            wcs=astWCS.WCS(img[0].header, mode = 'pyfits')
        simCat=catalogs.generateRandomSourcesCatalog(mapData, wcs, numSources)
        simCat.write(catalogFileName, overwrite = True)
        #catalogs.catalog2DS9(simCat, catalogFileName.replace(".fits", ".reg"))

        modelMap=maps.makeModelImage(mapData.shape, wcs, simCat, self.beamFileName)
        astImages.saveFITS(outMapFileName, modelMap+mapData, wcs)
        
        
    def cross_match(self, inCatalogFileName, outCatalogFileName, radiusArcmin = 1.0):
        """Cross matches input and output source catalogs.
        
        """
        inTab=atpy.Table().read(inCatalogFileName)
        outTab=atpy.Table().read(outCatalogFileName)
        self.inTab, self.outTab, self.rDeg=catalogs.crossMatch(inTab, outTab, 
                                                               radiusArcmin = radiusArcmin)
    
    
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
        
        
    def check_recovered_ratio(self, inKey, outKey, toleranceSigma = 1.0, expectedRatio = 1.0, SNRCut = 4,
                              SNRKey = 'fixed_SNR', errInKey = None, errOutKey = None,
                              plotLabel = None, plotsDir = "plots"):
        """Catalogs must have been cross matched before this can be run.
        Calculate the ratio of the columns pointed to by inKey, outKey. If tolerance (defined as
        1 - average ratio) is exceeded, the test is failed. SNRCut is applied to SNRKey in the output catalog.
        
        """

        inTab=self.inTab
        outTab=self.outTab
        mask=outTab[SNRKey] > SNRCut
        x=inTab[inKey]
        y=outTab[outKey]
        meanRatio=np.mean(y[mask])/np.mean(x[mask])
        bsRatios=[]
        for i in range(5000):
            indices=np.random.randint(0, len(x[mask]), len(x[mask]))
            bsX=x[mask][indices]
            bsY=y[mask][indices]
            bsRatios.append(np.mean(bsY)/np.mean(bsX))
        meanRatioErr=np.percentile(abs(bsRatios-meanRatio), 68.3)

        label="<input %s>/<output %s> = %.3f ± %.3f (%s > %.1f)" % (inKey, outKey, meanRatio, meanRatioErr, SNRKey, SNRCut)
        print("%s" % (label))

        if plotLabel is not None:
            plotMin=0.1*min([inTab[inKey].min(), outTab[outKey].min()])
            plotMax=10*max([inTab[inKey].max(), outTab[outKey].max()])
            plotRange=np.linspace(plotMin, plotMax, 100)
            plt.figure(figsize = (10, 8))
            if errInKey is not None and errOutKey is not None:
                plt.errorbar(inTab[inKey], outTab[outKey], yerr = outTab[errOutKey],
                     xerr = inTab[errInKey], elinewidth = 1.5, ecolor = '#AAAAFF',
                     fmt = 'D', ms = 6, label = None)
            else:
                plt.plot(inTab[inKey], outTab[outKey], 'D')
            plt.xlabel("input "+inKey)
            plt.ylabel("output "+outKey)
            plt.plot(plotRange, plotRange, 'k--')
            plt.xlim(plotMin, plotMax)
            plt.ylim(plotMin, plotMax)
            plt.loglog()
            plt.title(label, fontdict = {'size': plotTitleSize})
            plt.savefig(plotsDir+os.path.sep+plotLabel+"_XvY.png")
            plt.close()
        if abs((expectedRatio-meanRatio)/meanRatioErr) > toleranceSigma:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
            
            
    def check_recovered_positions(self, toleranceArcsec = 12.0, SNRKey = 'SNR', SNRMax = 10.0, 
                                  plotLabel = None, plotsDir = "plots"):
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
        if plotLabel is not None:
            plt.figure(figsize = (10, 8))
            plt.plot(SNRs, rDeg*3600., 'r.')
            plt.semilogx()
            plt.xlabel("SNR")
            plt.ylabel('Position offset (")')
            plt.title(label, fontdict = {'size': plotTitleSize})
            plt.savefig(plotsDir+os.path.sep+plotLabel+"_posRec.png")
            plt.close()
        if medOffsetArcmin_SNR10*60 > toleranceArcsec:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
    
    
    def subtract_maps(self, map1FileName, map2FileName, outMapFileName):
        with pyfits.open(map1FileName) as img1:
            d1=img1[0].data
            wcs=astWCS.WCS(img1[0].header, mode = 'pyfits')
        with pyfits.open(map2FileName) as img2:
            d2=img2[0].data
        maps.saveFITS(outMapFileName, d1-d2, wcs)


    def check_map_sigma(self, map1FileName, expectedSigma, tol = 1.0):
        with pyfits.open(map1FileName) as img1:
            d1=img1[0].data
        print("args:", map1FileName, expectedSigma)
        diff=abs(np.std(d1.flatten()) - float(expectedSigma))
        print("diff from expected sigma", diff)
        if diff > tol:
            self._status="FAILED"
        else:
            self._status="SUCCESS"
    
            
    def status_should_be(self, expected_status):
        if expected_status != self._status:
            raise AssertionError("Expected status to be '%s' but was '%s'."
                                 % (expected_status, self._status))


    def _run_command(self, args):
        thisDir=os.getcwd()
        os.chdir(self.runDir)
        print(args)
        process=subprocess.run(args, universal_newlines=True, stdout=subprocess.PIPE,
                               stderr=subprocess.STDOUT)
        if process.returncode != 0:
             print(process.stdout)
             raise AssertionError("Return code of '%s' is non-zero." % (str(args)))
        os.chdir(thisDir)
