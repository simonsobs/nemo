"""

N(z) check - for comparing across codes using different mass function implementations etc..

NOTE: We messed up when setting this task, so the mass limit is M500c > 5e13 MSun/h

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
if __name__ == '__main__':

    parser=argparse.ArgumentParser("NzCheck.py")
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
                        default = 0.1, type = float)
    args = parser.parse_args()
    
    configFileName=args.configFileName
    selFnDir=args.selFnDir
    SNRCut=args.SNRCut
    zStep=args.zStep
    footprintLabel=args.footprint
    
    print(">>> Setting up SNR > %.2f selection function ..." % (SNRCut))
    selFn=completeness.SelFn(selFnDir, SNRCut, configFileName = configFileName, 
                             footprintLabel = footprintLabel, zStep = zStep)

    scalingRelationDict=selFn.scalingRelationDict
    H0, Om0, Ob0, sigma8, ns = 70.0, 0.30, 0.05, 0.80, 0.95
    h=H0/100.0
    selFn.update(H0, Om0, Ob0, sigma8, ns, scalingRelationDict = scalingRelationDict)
    print("Total area = %.3f square degrees" % (selFn.totalAreaDeg2))
    
    # N(z) with M500c > 5e13 MSun - no selection function applied
    countsByRedshift=selFn.mockSurvey.clusterCount[:, np.greater(selFn.mockSurvey.log10M, np.log10(5e13/h))].sum(axis = 1)
    with open("NzCheck_noSelFn.csv", "w") as outFile:
        for i in range(len(selFn.mockSurvey.z)):
            outFile.write("%.1f <= z < %.1f\t%.3f\n" % (selFn.mockSurvey.zBinEdges[i], selFn.mockSurvey.zBinEdges[i+1], countsByRedshift[i]))
    
    # N(z) with M500c > 5e13 MSun - with S/N > 5 selection function applied
    predMz=selFn.compMz*selFn.mockSurvey.clusterCount
    countsByRedshift=predMz[:, np.greater(selFn.mockSurvey.log10M, np.log10(5e13/h))].sum(axis = 1)
    with open("NzCheck_withSelFn.csv", "w") as outFile:
        for i in range(len(selFn.mockSurvey.z)):
            outFile.write("%.1f <= z < %.1f\t%.3f\n" % (selFn.mockSurvey.zBinEdges[i], selFn.mockSurvey.zBinEdges[i+1], countsByRedshift[i]))

       
