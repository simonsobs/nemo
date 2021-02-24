.. _ConfigReference:

=============================
Configuration File Parameters
=============================

This section describes all the options that can be set in the YAML-format configuration file.


Input Maps and Masks
====================

**unfilteredMaps** (list):

    This is a list of dictionaries, each of which contains information about a sky map
    at a particular observing frequency. All maps given in this list must have the same
    dimensions and World Coordinate System. Units should be ΔT (μK) with respect to the
    CMB. Each map dictionary must have the following keys:
        
        **mapFileName** (str):
        
            The path to a FITS-format map.
            
        **weightsFileName** (str):
        
            The path to a FITS-format inverse-variance map. This may be set to ``null`` if
            necessary.
            
        **obsFreqGHz** (float):
        
            The frequency of the map. For cluster searches, this should be the SZ-weighted
            effective frequency of the passband corresponding to the given map - this value
            is used to convert from ΔT (μK) to Compton y. For sources, this value is used
            to convert source amplitudes from ΔT (μK) to flux density (Jy).
            
        **units** (str):
            
            Map pixel intensity units. At present, only 'uK' is accepted, but in principle
            **Nemo** could be extended to use maps in surface brightness units (MJy/sr), for
            example.
            
        **beamFileName** (str):
            
            Path to a text file containing the beam profile, in the format used by ACT.
            The beam solid angle must be present in the header of this file, in order
            for :ref:`nemoCommand` to be able to output source flux densities in Jy.
            (e.g., the ACT beam profile files contain lines such as
            ``# solid angle = 188.17 +/- 4.84 nsr`` as part of the header).
    
    *Example:*
    
    .. code-block:: yaml

       unfilteredMaps:
           - {mapFileName: "maps/Jun2020/act_s08_s18_cmb_f150_daynight_map.fits",
              weightsFileName: "maps/Jun2020/act_s08_s18_cmb_f150_daynight_ivar.fits",
              obsFreqGHz: 149.6, units: 'uK',
              beamFileName: "maps/Jun2020/beams/s16_pa2_f150_nohwp_night_beam_profile_jitter.txt"}
           - {mapFileName: "maps/Jun2020/act_s08_s18_cmb_f090_daynight_map.fits",
              weightsFileName: "maps/Jun2020/act_s08_s18_cmb_f090_daynight_ivar.fits",
              obsFreqGHz: 97.8, units: 'uK',
              beamFileName: "maps/Jun2020/beams/s16_pa3_f090_nohwp_night_beam_profile_jitter.txt"}   


**surveyMask** (str):

    This is the path to a FITS format image, with the same dimensions and pixelization
    as the input map file names (see unfilteredMaps above), with values of 1 taken to
    indicate valid area, and values of 0 taken to indicate areas that should be ignored.
    Objects found within pixels with value 0 in this map will not be included in the
    output catalogs produced by :ref:`nemoCommand`. Maps which are compressed using the
    ``PLIO_1`` method (as implemented by `astropy.io.fits <https://docs.astropy.org/en/stable/io/fits/>`_)
    are supported.
    
    *Example:*
    
    .. code-block:: yaml
    
       surveyMask: "maps/Jun2020/AdvACTSurveyMask_v7_galLatCut_S18-dust-artifacts-extended-post10mJy.fits"


**maskPointSourcesFromCatalog** (list):

    This is a list, each element of which points to a file containing an object catalog
    (FITS format is guaranteed to work, but in principle any format that the
    `astropy.table <https://docs.astropy.org/en/stable/table/index.html>`_ module
    understands can be used). The catalog should contain object coordinates decimal
    degrees in columns named ``RADeg``, ``decDeg``, and either a column specifying the
    masking radius in arcmin to use (``rArcmin``), or shape information columns
    (as produced by running :ref:`nemoCommand` with ``measureShapes: True`` set).
    
    *Example:*
    
    .. code-block:: yaml
       
       maskPointSourcesFromCatalog:
           - "PSCatalog_rArcmin/PS_S18_f150_auto_rArcmin.fits"
           - "customPSMask_S18/customPSCatalog_S18.fits"
           

**noiseMaskCatalog** (str):
    
    This is the path to a **Nemo** object catalog (containing either sources or clusters).
    If this is given, a model image will be constructed from the catalog on-the-fly when
    running the :ref:`nemoCommand` command, and subtracted from the maps used to create
    the noise term in the matched filters. This mitigates potential bias and signal loss
    from using the map itself to construct the filter noise term.


Output Options
==============
outputDir
makeQuickLookMaps
stitchTiles

Tiling
======
makeTileDir
selFnFootprints
tileOverlapDeg
tileDefLabel
tileDefinitions
tileNoiseRegions
tileNameList

Object Detection and Photometry
===============================
thresholdSigma
forcedPhotometryCatalog
photFilter
measureShapes
useInterpolator
minObjPix
findCenterOfMass
removeRings
ringThresholdSigma
rejectBorder
objIdent
longNames
catalogCuts
twoPass
undoPixelWindow


Filters
=======
mapFilters
allFilters
GNFWParams

Cluster Mass Estimates
======================
fitQ
massOptions

Selection Function
==================
calcSelFn
selFnOptions
massLimitMaps


Mock Catalogs
=============
applyPoissonScatter
applyIntrinsicScatter
applyNoiseScatter
makeMockCatalogs
numMocksToMake
seed

Source Injection Simulations
============================
sourceInjectionTest
sourceInjectionIterations
sourceInjectionModels
sourcesPerTile

Other Diagnostics
=================

.. note::  The parameters listed in this section are depreciated and not currently enabled,
           i.e., they are ignored by the :ref:`nemoCommand` command.

numSkySims
estimateContaminationFromInvertedMaps
estimateContaminationFromSkySim
