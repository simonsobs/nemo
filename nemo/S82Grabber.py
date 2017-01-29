# -*- coding: utf-8 -*-
"""This module provides a class for getting the S82 FITS images and mosaicing them (uses IRAF currently, should
change that) given RA, dec coords.

"""

import urllib
import glob
import os
import sys
import numpy
import pyfits
import time
import operator
try:
    from pyraf import iraf
except:
    pass
from astLib import *
import nemo
import gzip
import IPython

#-------------------------------------------------------------------------------------------------------------
class S82Grabber:
    
    def __init__(self, cacheDir = "default"):
        """Initialises a S82Grabber object. If the image array cannot be found, then it is created by scraping
        the S82 webpages.
        
        If cacheDir == 'default', FITS images are downloaded to $HOME/.S82Grabber/cache, otherwise they are 
        stored under the given cacheDir path
        
        """
        
        self.verbose=True   # show status messages about e.g. downloading images if True
                
        S82ArrayPath=nemo.__path__[0]+os.path.sep+"data"+os.path.sep+"S82ImageMapArray.npy"

        dataDir=os.environ['HOME']+os.path.sep+".S82Grabber"
        if os.path.exists(dataDir) == False:
            os.makedirs(dataDir)
        
        # We'll download and store SDSS images under this dir - we allow path to be changed if needed
        if cacheDir == "default":
            cacheDir=dataDir+os.path.sep+"cache"
            if os.path.exists(cacheDir) == False:
                os.makedirs(cacheDir)
        self.cacheDir=cacheDir
        
        self.columns=['imageNumber', 'RADeg', 'decDeg', 'run', 'rerun', 'camcol', 'field']
        
        if os.path.exists(S82ArrayPath) == False:
            print "How did this happen?"
            ipshell()
            sys.exit()
            runs=[100006, 200006]
            reruns=[2]
            camcols=range(1,7)
            fields=range(62,801)

            ras=[]
            decs=[]
            imageMapArray=[]
            imageCounter=0
            for run in runs:
                print "> run = %d" % (run)
                for rerun in reruns:
                    for camcol in camcols:
                        print ">> camcol = %d" % (camcol)
                        for field in fields:
                            print ">>> field = %d" % (field)
                            imageCounter=imageCounter+1
                            response=urllib.urlretrieve("http://das.sdss.org/www/cgi-bin/field?RUN=%d&RERUN=%d&CAMCOL=%d&FIELD=%d" % (run, rerun, camcol, field))
                            inFile=file(response[0], "r")
                            lines=inFile.readlines()
                            inFile.close()
                            for line in lines:
                                if line.find("Coordinates of center:") != -1:
                                    coordBits=line[line.find("Coordinates of center: ")+len("Coordinates of center: "):line.find("</a>")].split(",")
                                    ra=float(coordBits[0])
                                    dec=float(coordBits[1])
                                    imageMapArray.append([imageCounter, ra, dec, run, rerun, camcol, field])
                                    time.sleep(0.5)

            self.imageMapArray=numpy.array(imageMapArray)
            numpy.save(S82ArrayPath, imageMapArray)

        else:
            
            self.imageMapArray=numpy.load(S82ArrayPath)
    
    
    def getNearestFieldsList(self, RADeg, decDeg, number = 5):
        """Fetches info on a given number of fields closest to given RA, dec coords, sorted by radial 
        distance. An assumption is made here that will only hold on the celestial equator for speed!
        
        If the nearest fields are > 0.5 degrees away, we return an empty list
        
        """
        
        diffRA=self.imageMapArray[:, self.columns.index('RADeg')]-RADeg
        diffDec=self.imageMapArray[:, self.columns.index('decDeg')]-decDeg
        diff=(abs(diffRA)+abs(diffDec))/2.0
        
        sortable=numpy.array([self.imageMapArray[:, 0], diff])
        sortable=sortable.transpose()
        sortable=sorted(sortable.tolist(), key=operator.itemgetter(1))
        sortable=numpy.array(sortable)
        
        imageNumbers=[]
        for i in range(number):
            if sortable[i, 1] < 0.3:    # don't bother with anything more than 0.5 degrees away
                imageNumbers.append(int(sortable[i, 0]))
        
        fieldsDictList=[]
        for i in imageNumbers:
            row=self.imageMapArray[numpy.where(self.imageMapArray[:, 0] == i)][0]
            fieldDict={}
            for c in self.columns:
                fieldDict[c]=row[self.columns.index(c)]
            fieldsDictList.append(fieldDict)
            
        return fieldsDictList
    
    
    def fetchImages(self, fieldDict, bands = ['g', 'r', 'i']):
        """Fetches S82 .fits images for given bands using info in the fieldDict. Use getNearestFieldsList to
        make a fieldDict at a given RA, dec.
        
        """
        
        for band in bands:
            urlDict=self.makeImageURLAndPath(fieldDict, band)
            gotFile=False
            count=0
            while gotFile == False:
                count=count+1
                # SDSS images are approx 10 Mb
                if os.path.exists(urlDict['cachedPath']) == True and os.path.getsize(urlDict['cachedPath'])/1e6 > 11: 
                    gotFile=True
                if gotFile == False:
                    if os.path.exists(urlDict['cachedPath']) == True:
                        os.remove(urlDict['cachedPath'])
                    if self.verbose == True:
                        print "... S82Grabber: downloading image from %s" % (urlDict['url'])
                    # Use wget here because I think the file scanning thing built into the proxy looks just 
                    # like the proxy taking a long time to respond, which causes urlretrieve to time out
                    # sometimes (?)
                    topDir=os.getcwd()
                    os.chdir(os.path.split(urlDict['cachedPath'])[0])
                    os.system("wget %s" % (urlDict['url']))
                    os.chdir(topDir)
                    #response=urllib.urlretrieve(urlDict['url'], urlDict['cachedPath'])
                if count > 5:
                    print "What the hell happened?"
                    ipshell()
                    sys.exit()
                
    def makeImageURLAndPath(self, fieldDict, band):
        """Given a fieldDict, return as a dict the URL on the DAS and the local cache file name where the 
        image is/will be stored.
        
        """
        run=fieldDict['run']
        rerun=fieldDict['rerun']
        camcol=fieldDict['camcol']
        field=fieldDict['field']
        url="http://das.sdss.org/imaging/%d/%d/corr/%d/fpC-%d-%s%d-%04d.fit.gz" % \
                (run, rerun, camcol, run, band, camcol, field)
        outPath=self.cacheDir+os.path.sep+url.split("/")[-1]
        
        return {'url': url, 'cachedPath': outPath}
        

    def listToString(self, inList, prefix = "", delimiter=","):
        """Converts a list into a comma delimited string, can add prefixes if needed.
        
        """
        outString=""
        for i in inList:
            if outString != "":
                outString=outString+delimiter
            outString=outString+prefix+str(i)
        return outString
    
    
    def mosaicFields(self, fieldsDictList, outFileName, band = 'r'):
        """Mosaic the fields specified by the fieldsDictList and band together to make a huge image, written
        to outFileName as a .fits file.
        
        """
        
        # First, copy images here and decompress
        fileList=[]
        for fieldDict in fieldsDictList:
            urlDict=self.makeImageURLAndPath(fieldDict, band)
            shortPath=os.path.split(urlDict['cachedPath'])[-1]
            if os.path.exists(urlDict['cachedPath']) == False:
                raise Exception, "can't find locally cached image %s" % (urlDict['cachedPath'])
            if os.path.exists(shortPath.replace(".gz", "")) == False:
                os.system("cp %s ." % (urlDict['cachedPath']))
                # Apparently, it's not possible to download gzipped files over http at Nottingham without
                # some security scanner program decompressing them before they reach me. It doesn't remove
                # the .gz suffix either. Charming.
                gzippedFile=gzip.GzipFile(shortPath)
                try:
                    gzippedFile.read()
                    os.system("gunzip %s" % (shortPath))
                except IOError:
                    os.system("mv %s %s" % (shortPath, shortPath.replace(".gz", "")))
            fileList.append(shortPath.replace(".gz", ""))
        
        if os.path.exists(outFileName) == True:
            os.remove(outFileName)
        if len(fileList) > 0:
            iraf.imcombine(input=self.listToString(fileList), output=outFileName, offsets='wcs', 
                           combine='average')
        
        # Clean up images
        for f in fileList:
            os.remove(f)        
        
