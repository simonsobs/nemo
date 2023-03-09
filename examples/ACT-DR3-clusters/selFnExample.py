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
from nemo import completeness
import argparse
import time
import IPython

#------------------------------------------------------------------------------------------------------------
def printNumClusters(H0, Om0, Ob0, sigma_8, n_s, scalingRelationDict = None):
    """Example of how to update selection function.
    
    """
    
    # This function has to be called every time a parameter value is changed
    selFn.update(H0, Om0, Ob0, sigma_8, n_s, scalingRelationDict = scalingRelationDict)
    
    print(">>> H0 = %.2f km/s/Mpc Om0 = %.2f Ob0 = %.2f sigma_8 = %.2f tenToA0 = %.3e sigma_int = %.2f" \
        % (selFn.mockSurvey.H0, selFn.mockSurvey.Om0, selFn.mockSurvey.Ob0, selFn.mockSurvey.sigma8, 
           selFn.scalingRelationDict['tenToA0'], selFn.scalingRelationDict['sigma_int']))
    numClusters=selFn.mockSurvey.calcNumClustersExpected(compMz = selFn.compMz)
    print("... number of clusters expected = %.2f (in %.2f square degrees) ..." %  (numClusters, selFn.totalAreaDeg2))
    
#------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':

    parser=argparse.ArgumentParser("selFnExample.py")
    parser.add_argument("selFnDir", help="""Directory containing files needed for computing the selection 
                        function.""")
    parser.add_argument("-c", "--config", dest = "configFileName", help="""A .yml configuration file. If
                        this is not given, the config.yml file in selFnDir will be used.""",
                        default = None)
    parser.add_argument("-f", "--footprint", dest = "footprint", help="""Footprint to use, e.g., DES,
                        HSC, KiDS (default: full).""", default = None)
    parser.add_argument("-S", "--SNR-cut", dest = "SNRCut", help="""Use only clusters with fixed_SNR > 
                        this value.""", default = 5.0, type = float)
    parser.add_argument("-z", "--z-step", dest = "zStep", help="""Redshift bin width (default: 0.1).""", 
                        default = 0.02, type = float)
    args = parser.parse_args()
    
    configFileName=args.configFileName
    selFnDir=args.selFnDir
    SNRCut=args.SNRCut
    zStep=args.zStep
    footprintLabel=args.footprint
    
    print(">>> Setting up SNR > %.2f selection function ..." % (SNRCut))
    selFn=completeness.SelFn(selFnDir, SNRCut, configFileName = configFileName, 
                             footprint = footprintLabel, zStep = zStep)
    massLabel="M%d%s" % (selFn.mockSurvey.delta, selFn.mockSurvey.rhoType[0])

    # If we want to play with scaling relation also
    scalingRelationDict=selFn.scalingRelationDict
    
    # Default parameters    
    t0=time.time()
    H0, Om0, Ob0, sigma_8, n_s = 70.0, 0.30, 0.05, 0.80, 0.96
    printNumClusters(H0, Om0, Ob0, sigma_8, n_s, scalingRelationDict = scalingRelationDict)
    t1=time.time()
    print("... iteration took %.3f sec ..." % (t1-t0))
    completeness.makeMzCompletenessPlot(selFn.compMz, selFn.mockSurvey.log10M, selFn.mockSurvey.z, "test0", massLabel, "test0_Mz.pdf")
    
    # A complete iteration of changing cosmological parameters    
    t0=time.time()
    H0, Om0, Ob0, sigma_8, n_s = 72.0, 0.27, 0.05, 0.74, 1.00
    printNumClusters(H0, Om0, Ob0, sigma_8, n_s, scalingRelationDict = scalingRelationDict)
    t1=time.time()
    print("... iteration took %.3f sec ..." % (t1-t0))
    completeness.makeMzCompletenessPlot(selFn.compMz, selFn.mockSurvey.log10M, selFn.mockSurvey.z, "test1", massLabel, "test1_Mz.pdf")
    
    # Changing the scaling relation normalisation
    t0=time.time()
    scalingRelationDict['tenToA0']=3e-5
    printNumClusters(H0, Om0, Ob0, sigma_8, n_s, scalingRelationDict = scalingRelationDict)
    t1=time.time()
    print("... iteration took %.3f sec ..." % (t1-t0))
    completeness.makeMzCompletenessPlot(selFn.compMz, selFn.mockSurvey.log10M, selFn.mockSurvey.z, "test2", massLabel, "test2_Mz.pdf")

    # Changing the scaling relation intrinsic scatter
    t0=time.time()
    scalingRelationDict['tenToA0']=4.95e-5
    scalingRelationDict['sigma_int']=0.001
    printNumClusters(H0, Om0, Ob0, sigma_8, n_s, scalingRelationDict = scalingRelationDict)
    t1=time.time()
    print("... iteration took %.3f sec ..." % (t1-t0))
    completeness.makeMzCompletenessPlot(selFn.compMz, selFn.mockSurvey.log10M, selFn.mockSurvey.z, "test3", massLabel, "test3_Mz.pdf")
    
