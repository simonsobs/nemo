"""Check that objects near RA = +/-180 deg are painted into maps correctly.

:meth:`nemo.signals._paintSignalMap` paints each object into a postage stamp rather than over the whole
map, for speed. The stamp bounds used to be found by feeding an RA, dec range to
``astImages.clipUsingRADecCoords``, which folds RA at the +/-180 deg branch cut. For any object within
``maxSizeDeg`` of there (in a map whose RA axis runs 180 -> -180 deg, as ACT/SO maps do), that returned the
whole map *except* the stamp, so the object was silently not painted at all - nemoModel simply dropped
clusters and sources in a band around RA = 180 deg.

The bounds are now found in pixel coordinates instead (see :meth:`nemo.signals._getStampPixelBounds`),
taking into account that the angle spanned by a pixel in the RA direction varies with declination.

Objects right on the wrap were dropped even earlier than that, by the valid area test in
:meth:`nemo.catalogs.getCatalogWithinImage`, which now folds x into the map for maps that cover all RA.

This compares the stamp path against painting over the whole map, which pixell's object painter does
handle across the wrap, for point sources and clusters at a range of positions either side of the seam
and at a range of declinations.

It can be run standalone, or via the Robot Framework ``quick.robot`` suite (see the "Objects near RA =
180 deg are painted" test case and the ``Check RA wrap painting`` keyword).

Usage:
    python test_RAWrapPainting.py [beamFileName]

"""

import os, sys
import numpy as np
import astropy.table as atpy
from pixell import enmap, utils
from astLib import astWCS
from nemo import maps, signals

defaultBeamFileName=os.path.join(os.path.dirname(os.path.abspath(__file__)), "testsCache", "maps",
                                 "s16_pa3_f090_nohwp_night_beam_profile_jitter.txt")

# Positions to test - either side of the wrap, on it, plus controls well away from it.
# NOTE: RA = +/-180 deg exactly, and the sliver either side of it, are the awkward cases. There,
# wcs2pix reports an object at the far end of the map, or at x = -1e-10 rather than 0, so as well as the
# postage stamp bounds this exercises the RA wrap handling in catalogs.getCatalogWithinImage (without which
# such objects are dropped before they ever reach the painter).
testRADegs=[0.0, 90.0, 170.0, 179.0, 179.5, 179.9, 179.99, 180.0,
            -180.0, -179.99, -179.9, -179.5, -179.0, -170.0]

#------------------------------------------------------------------------------------------------------------
def makeTestGeometry(decMinDeg, decMaxDeg, resArcmin = 2.0):
    """Returns (shape, astWCS.WCS) for a full-RA CAR strip, i.e., with the RA axis running from
    180 -> -180 deg, so that the RA wrap falls at the map edge (as it does for ACT/SO maps).

    """

    box=np.array([[decMinDeg, 180.0], [decMaxDeg, -180.0]])*utils.degree
    shape, awcs=enmap.geometry(pos = box, res = resArcmin*utils.arcmin, proj = 'car')

    return shape, astWCS.WCS(awcs.to_header(), mode = 'pyfits')

#------------------------------------------------------------------------------------------------------------
def runRAWrapCheck(beamFileName = None, rtol = 1e-2):
    """Checks that sources and clusters either side of RA = +/-180 deg are painted with the same flux as
    they are when painted over the whole map (which is the slow path, but is wrap-safe).

    Args:
        beamFileName (:obj:`str`, optional): Path to a beam profile file.
        rtol (:obj:`float`, optional): Fractional tolerance on peak signal and total flux. The painter
            itself varies at the ~1e-4 level with sub-pixel position and stamp size, so this only needs to
            be tight enough to catch objects that are misplaced or missing.

    Returns:
        True if all checks pass.

    """

    if beamFileName is None:
        beamFileName=defaultBeamFileName
    beam=signals.BeamProfile(beamFileName = beamFileName)

    allPass=True

    # Point sources, at a range of declinations. The number of pixels spanned by maxSizeDeg in the RA
    # direction goes as 1/cos(dec), so the band affected by the bug was wider at high declination.
    for decMinDeg, decMaxDeg, decDegs in [(-10.0, 10.0, [0.0, 9.0]), (60.0, 88.0, [62.0, 85.0])]:
        shape, wcs=makeTestGeometry(decMinDeg, decMaxDeg)
        for decDeg in decDegs:
            for RADeg in testRADegs:
                tab=atpy.Table()
                tab['RADeg']=[RADeg]
                tab['decDeg']=[decDeg]
                tab['deltaT_c']=[1000.0]
                painted=maps.makeModelImage(shape, wcs, tab, beamFileName, applyPixelWindow = False)
                reference=signals.makeBeamModelSignalMap(shape, wcs, beam, RADeg = RADeg, decDeg = decDeg,
                                                         maxSizeDeg = 1.0, amplitude = 1000.0)
                label="source  RADeg = %7.2f  decDeg = %5.1f" % (RADeg, decDeg)
                allPass=_compare(painted, reference, label, rtol) and allPass
                del painted, reference

    # Clusters - here maxSizeDeg is set by the model, and so is much larger than for a point source
    shape, wcs=makeTestGeometry(-10.0, 10.0)
    z, M500=0.4, 2e14
    maxSizeDeg=10*signals.calcTheta500Arcmin(z, M500, signals.fiducialCosmoModel)/60
    for RADeg in testRADegs:
        for decDeg in [0.0, 9.0]:
            tab=atpy.Table()
            tab['RADeg']=[RADeg]
            tab['decDeg']=[decDeg]
            tab['y_c']=[5.0]
            tab['template']=['Arnaud_M2e14_z0p4']
            painted=maps.makeModelImage(shape, wcs, tab, beamFileName, obsFreqGHz = 150.0,
                                        applyPixelWindow = False, maxSizeDegMultiplier = 10)
            reference=signals.makeArnaudModelSignalMap(z, M500, shape, wcs, beam = beam, RADeg = RADeg,
                                                      decDeg = decDeg, amplitude = 5e-4,
                                                      maxSizeDeg = maxSizeDeg, convolveWithBeam = True,
                                                      obsFrequencyGHz = 150.0)
            label="cluster RADeg = %7.2f  decDeg = %5.1f" % (RADeg, decDeg)
            allPass=_compare(painted, reference, label, rtol) and allPass
            del painted, reference

    print("\nRA WRAP PAINTING OK (rtol = %.0e)" % (rtol) if allPass else "\nRA WRAP PAINTING FAILED")

    return allPass

#------------------------------------------------------------------------------------------------------------
def _compare(painted, reference, label, rtol):
    """Compares a model image painted via the postage stamp path with one painted over the whole map.

    """

    if painted is None:
        print("  [FAIL] %s - nothing painted (object dropped from the catalog)" % (label))
        return False
    refPeak=abs(reference).max()
    refFlux=reference.sum()
    peakDiff=abs(abs(painted).max()-refPeak)/refPeak
    fluxDiff=abs(painted.sum()-refFlux)/abs(refFlux)
    ok=peakDiff < rtol and fluxDiff < rtol
    print("  [%s] %s  peak = %9.3f (%.1e)  flux = %11.3f (%.1e)" % ("PASS" if ok else "FAIL", label,
                                                                    abs(painted).max(), peakDiff,
                                                                    painted.sum(), fluxDiff))

    return ok

#------------------------------------------------------------------------------------------------------------
if __name__ == '__main__':
    beamFileName=sys.argv[1] if len(sys.argv) > 1 else None
    sys.exit(0 if runRAWrapCheck(beamFileName) else 1)
