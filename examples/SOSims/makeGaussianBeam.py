"""

Code to churn out a Gaussian beam in a format that nemo understands

"""

import os
import sys
import argparse
import numpy as np
import IPython

parser = argparse.ArgumentParser()
parser.add_argument("FWHMArcmin")
parser.add_argument("outFileName")
args = parser.parse_args()

FWHMArcmin=float(args.FWHMArcmin)
sigmaDeg=(FWHMArcmin/(2*np.sqrt(2*np.log(2))))/60. 

# Radial distance in degrees
xRange=np.linspace(0.0, 0.5, 1800)

# Gaussian
resp=np.exp(-np.power(xRange, 2)/(2*np.power(sigmaDeg, 2)))

# Solid angle in nsr
solidAngle=np.trapz(resp, np.pi*(np.radians(xRange))**2)*1e9

# Output format: plain text; 1st column degrees; 2nd column response (normalised to 1.0 at 0.00)
with open(args.outFileName, "w") as outFile:
    outFile.write("# solid angle = %.2f nsr\n" % (solidAngle))
    for x, r in zip(xRange, resp):
        outFile.write("%.6f %.6e\n" % (x, r))
