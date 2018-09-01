"""

An example of how to use the nemo selection function

"""

import os
import sys
import numpy as np
import pylab as plt
import astropy.table as atpy
from astLib import *
from scipy import stats
from scipy import interpolate
from nemo import SelFn
from nemo import selFnTools
import time
import IPython

#------------------------------------------------------------------------------------------------------------
def printNumClusters(H0, Om0, Ob0, sigma_8, scalingRelationDict = None):
    """Example of how to update selection function.
    
    """
    
    # This function has to be called every time a parameter value is changed
    selFn.update(H0, Om0, Ob0, sigma_8, scalingRelationDict = scalingRelationDict)
    
    print(">>> H0 = %.2f km/s/Mpc Om0 = %.2f Ob0 = %.2f sigma_8 = %.2f tenToA0 = %.3e sigma_int = %.2f" \
        % (selFn.mockSurvey.H0, selFn.mockSurvey.Om0, selFn.mockSurvey.Ob0, selFn.mockSurvey.sigma_8, 
           selFn.scalingRelationDict['tenToA0'], selFn.scalingRelationDict['sigma_int']))
    numClusters=selFn.mockSurvey.calcNumClustersExpected(selFn = selFn.compMz)
    print("... number of clusters expected = %.2f (in %.2f square degrees) ..." %  (numClusters, selFn.totalAreaDeg2))
    
#------------------------------------------------------------------------------------------------------------
# Main
if len(sys.argv) < 2:
    print("Run: % selFnExample.py <.yml config file>")
else:
    
    parDictFileName=sys.argv[1]
    
    SNRCut=5.0
    
    print(">>> Setting up SNR > %.2f selection function ..." % (SNRCut))
    selFn=SelFn.SelFn(parDictFileName, SNRCut)
    
    # If we want to play with scaling relation also
    scalingRelationDict=selFn.scalingRelationDict
    
    # A complete iteration of changing cosmological parameters    
    t0=time.time()
    H0, Om0, Ob0, sigma_8 = 72.0, 0.27, 0.05, 0.74    
    printNumClusters(H0, Om0, Ob0, sigma_8, scalingRelationDict = scalingRelationDict)
    t1=time.time()
    print("... iteration took %.3f sec ..." % (t1-t0))
    selFnTools.makeMzCompletenessPlot(selFn.compMz, selFn.mockSurvey.log10M, selFn.mockSurvey.z, "test1", "test1_Mz.pdf")
    
    # Changing the scaling relation normalisation
    t0=time.time()
    scalingRelationDict['tenToA0']=3e-5
    printNumClusters(H0, Om0, Ob0, sigma_8, scalingRelationDict = scalingRelationDict)
    t1=time.time()
    print("... iteration took %.3f sec ..." % (t1-t0))
    selFnTools.makeMzCompletenessPlot(selFn.compMz, selFn.mockSurvey.log10M, selFn.mockSurvey.z, "test2", "test2_Mz.pdf")

    # Changing the scaling relation intrinsic scatter
    t0=time.time()
    scalingRelationDict['tenToA0']=4.95e-5
    scalingRelationDict['sigma_int']=0.001
    printNumClusters(H0, Om0, Ob0, sigma_8, scalingRelationDict = scalingRelationDict)
    t1=time.time()
    print("... iteration took %.3f sec ..." % (t1-t0))
    selFnTools.makeMzCompletenessPlot(selFn.compMz, selFn.mockSurvey.log10M, selFn.mockSurvey.z, "test3", "test3_Mz.pdf")    