"""Numerical-equivalence check for the vectorized fast-completeness rewrite.

`SelFn.calcFastCompletenessInTile` / `_calcCompMzCube` were rewritten to vectorize over the RMS-noise
rows and to re-use the intrinsic-scatter kernel across S/N bins. This script reproduces the *original*
per-RMS-row algorithm (`referenceCompMzCube` below, a verbatim copy of the old code) and asserts that the
new implementation reproduces it to floating-point tolerance, across:

  - the intrinsic-scatter path with the optimization-bias correction + truncation (the real hot path),
  - the zero-scatter (pure erf) path,
  - the no-bias-model path,

and for several scaling-relation parameter sets and tiles.

Usage:
    python test_fastCompMz_equivalence.py [selFnDir] [nTiles]

Defaults to the in-repo quickstart selFn dir.
"""

import os, sys
import numpy as np
from nemo import completeness

defaultSelFnDir = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                               "..", "examples", "quickstart", "quickstart-clusters", "selFn")
selFnDir = sys.argv[1] if len(sys.argv) > 1 else defaultSelFnDir
nTiles = int(sys.argv[2]) if len(sys.argv) > 2 else 6

RTOL, ATOL = 1e-9, 1e-12


def referenceCompMzCube(selFn, y0Grid, RMSTab, snBins):
    """Verbatim re-implementation of the ORIGINAL _fastCompMz / calcFastCompletenessInTile loop,
    used only as the ground truth for the equivalence check."""
    cube = np.zeros((len(snBins), y0Grid.shape[0], y0Grid.shape[1]))
    for k in range(len(snBins)):
        minSN, maxSN = snBins[k]
        areaWeights = RMSTab['areaDeg2'] / RMSTab['areaDeg2'].sum()
        compMzTile = np.zeros(y0Grid.shape)
        for i in range(len(RMSTab)):
            if selFn.biasModel is not None:
                trueSNR = y0Grid / RMSTab['y0RMS'][i]
                corrFactors = selFn.biasModel['func'](trueSNR, selFn.biasModel['params'])
                if selFn.truncateDeltaSNR is not None:
                    corrFactors[trueSNR < minSN - selFn.truncateDeltaSNR] = 1.0
            else:
                corrFactors = np.ones(y0Grid.shape)
            if selFn.scalingRelationDict['sigma_int'] == 0:
                compMzTile = compMzTile + selFn._get_erf_diff(
                    (y0Grid * corrFactors) / RMSTab['y0RMS'][i], minSN, maxSN, selFn.SNRCut) * areaWeights[i]
            else:
                scatter = selFn.scalingRelationDict['sigma_int']
                lnyy = np.linspace(np.min(np.log(y0Grid)), np.max(np.log(y0Grid)), 44)
                yy0 = np.exp(lnyy)
                mu = np.log(y0Grid * corrFactors)
                fac = 1. / np.sqrt(2. * np.pi * scatter ** 2)
                arg = selFn._get_erf_diff(yy0 / RMSTab['y0RMS'][i], minSN, maxSN, minSN)
                cc = arg * areaWeights[i]
                arg0 = (lnyy[:, None, None] - mu) / (np.sqrt(2.) * scatter)
                args = fac * np.exp(-arg0 ** 2.) * cc[:, None, None]
                compMzTile = compMzTile + np.trapezoid(args, x=lnyy, axis=0)
        if selFn.maxTheta500Arcmin is not None:
            compMzTile = compMzTile * np.array(selFn._theta500Grid < selFn.maxTheta500Arcmin, dtype=float)
        cube[k] = compMzTile
    return cube


def buildSNBins(selFn):
    snBins = np.zeros((selFn.SNBinEdges.shape[0], 2))
    snBins[0] = [selFn.SNRCut, 1e5]
    snBins[1:, 0] = selFn.SNBinEdges[:-1]
    snBins[1:, 1] = selFn.SNBinEdges[1:]
    return snBins


def srd(tenToA0=2.25e-05, B0=0.08, sigma_int=0.2):
    return {'tenToA0': tenToA0, 'B0': B0, 'Mpivot': 3.0e+14, 'sigma_int': sigma_int,
            'Ez_gamma': 2.0, 'onePlusRedshift_power': 0.0, 'zpivot': 0.0}


print("Building SelFn from %s ..." % selFnDir)
SNBinEdges = np.logspace(np.log10(5.5), np.log10(50), 11)
selFn = completeness.SelFn(selFnDir, SNRCut=SNBinEdges[0], zStep=0.1, zMin=0.3, zMax=2.0,
                           massFunction='Tinker08', numMassBins=40, applyRelativisticCorrection=False,
                           rhoType='critical', delta=500, method='fast', QSource='fit',
                           footprint=None, maxFlags=None, downsampleRMS=8, setUpAreaMask=False,
                           biasModel={'func': completeness.optBiasPowerModelFunc, 'params': 2.1},
                           theoryCode='CCL', useAverageQ=False, massBinsTheory=200, zStepTheory=0.01,
                           SNBinEdges=SNBinEdges)
H0, Om0, Ob0, sigma8, ns = 67.62, 0.3116, 0.0492, 0.8149, 0.9709
tiles = list(selFn.tileNames)[:nTiles]
print("Checking %d tiles, %d S/N planes\n" % (len(tiles), selFn.SNBinEdges.shape[0]))

# (label, sets biasModel back on?, scalingRelationDict)
biasOn = {'func': completeness.optBiasPowerModelFunc, 'params': 2.1}
cases = [
    ("scatter>0, bias+truncation (hot path)", biasOn, 3.0,  srd(B0=0.08, sigma_int=0.20)),
    ("scatter>0, different params",           biasOn, 3.0,  srd(tenToA0=3.0e-05, B0=0.30, sigma_int=0.35)),
    ("scatter>0, no truncation",              biasOn, None, srd(B0=0.15, sigma_int=0.25)),
    ("scatter>0, no bias model",              None,   3.0,  srd(B0=0.10, sigma_int=0.30)),
    ("scatter==0, bias+truncation",           biasOn, 3.0,  srd(B0=0.12, sigma_int=0.0)),
    ("scatter==0, no bias model",             None,   3.0,  srd(B0=0.12, sigma_int=0.0)),
]

worstAbs, worstRel, allPass = 0.0, 0.0, True
for label, biasModel, truncDeltaSNR, scalingRelationDict in cases:
    selFn.biasModel = biasModel
    selFn.truncateDeltaSNR = truncDeltaSNR
    selFn.update(H0, Om0, Ob0, sigma8, ns, scalingRelationDict=scalingRelationDict)
    snBins = buildSNBins(selFn)
    caseAbs, caseRel, caseOK = 0.0, 0.0, True
    for tileName in tiles:
        y0Grid = selFn._makeSignalGrid(tileName=tileName)
        RMSTab = selFn.RMSDict[tileName]
        ref = referenceCompMzCube(selFn, y0Grid, RMSTab, snBins)
        new = selFn.calcFastCompletenessInTile(tileName, return_y0Grid=False)
        absDiff = np.abs(new - ref)
        relDiff = absDiff / np.abs(ref).clip(min=1e-30)
        caseAbs = max(caseAbs, absDiff.max())
        caseRel = max(caseRel, relDiff[ref != 0].max() if np.any(ref != 0) else 0.0)
        if not np.allclose(new, ref, rtol=RTOL, atol=ATOL):
            caseOK = False
    worstAbs, worstRel = max(worstAbs, caseAbs), max(worstRel, caseRel)
    allPass = allPass and caseOK
    print("  [%s] %-38s maxAbs=%.2e maxRel=%.2e" %
          ("PASS" if caseOK else "FAIL", label, caseAbs, caseRel))

print("\nWorst over all cases: maxAbs=%.2e maxRel=%.2e" % (worstAbs, worstRel))
if allPass:
    print("EQUIVALENCE OK (rtol=%.0e, atol=%.0e)" % (RTOL, ATOL))
    sys.exit(0)
else:
    print("EQUIVALENCE FAILED")
    sys.exit(1)
