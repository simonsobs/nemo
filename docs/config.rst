.. _ConfigReference:

=============================
Configuration File Parameters
=============================

This section describes all the options that can be set in the YAML-format configuration file.

.. _InputMaps:

Input Maps and Masks
====================

unfilteredMaps
^^^^^^^^^^^^^^

    This is a list of dictionaries, each of which contains information about a sky map
    at a particular observing frequency. All maps given in this list must have the same
    dimensions and World Coordinate System. Units should be ΔT (μK) with respect to the
    CMB. Each map dictionary must have the following keys:
        
        :mapFileName (str):
        
            The path to a FITS-format map.
            
        :weightsFileName (str):
        
            The path to a FITS-format inverse-variance map. This may be set to ``null`` if
            necessary.
            
        :obsFreqGHz (float):
        
            The frequency of the map. For cluster searches, this should be the SZ-weighted
            effective frequency of the passband corresponding to the given map - this value
            is used to convert from ΔT (μK) to Compton y. For sources, this value is used
            to convert source amplitudes from ΔT (μK) to flux density (Jy).
            
        :units (str):
            
            Map pixel intensity units. At present, only 'uK' is accepted, but in principle
            **Nemo** could be extended to use maps in surface brightness units (MJy/sr), for
            example.
            
        :beamFileName (str):
            
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


surveyMask
^^^^^^^^^^

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


maskPointSourcesFromCatalog
^^^^^^^^^^^^^^^^^^^^^^^^^^^

    This is a list, each element of which points to a file containing an object catalog
    (FITS format is guaranteed to work, but in principle any format that the
    `astropy.table <https://docs.astropy.org/en/stable/table/index.html>`_ module
    understands can be used). The catalog should contain object coordinates in decimal
    degrees in columns named ``RADeg``, ``decDeg``, and either a column specifying the
    masking radius in arcmin to use (``rArcmin``), or shape information columns
    (as produced by running :ref:`nemoCommand` with ``measureShapes: True`` set).
    
    *Example:*
    
    .. code-block:: yaml
       
       maskPointSourcesFromCatalog:
           - "PSCatalog_rArcmin/PS_S18_f150_auto_rArcmin.fits"
           - "customPSMask_S18/customPSCatalog_S18.fits"
           

noiseMaskCatalog
^^^^^^^^^^^^^^^^
    
    This is the path to a **Nemo** object catalog (containing either sources or clusters).
    If this is given, a model image will be constructed from the catalog on-the-fly when
    running the :ref:`nemoCommand` command, and subtracted from the maps used to create
    the noise term in the matched filters. This mitigates potential bias and signal loss
    from using the map itself to construct the filter noise term.

    *Example:*
    
    .. code-block:: yaml
    
       noiseMaskCatalog: "S18d_202006/S18d_202006_optimalCatalog.fits"


.. _OutputOptions:

Output Options
==============

outputDir
^^^^^^^^^
    
    Output from :ref:`nemoCommand` will be written into this directory. It is not
    necessary to include this parameter in any config file - by default,
    :ref:`nemoCommand` will set this to the name of the config file, minus the .yml
    extension.

    *Example:*
    
    .. code-block:: yaml
    
       outputDir: "nemoOutput"
    
    
stitchTiles
^^^^^^^^^^^
    
    If True, and :ref:`nemoCommand` is set to break a large map into tiles (see
    :ref:`Tiling` below), all of the output filtered tile maps will be combined together
    to make a single,
    monolithic map. These are placed in the ``nemoOutput/filteredMaps`` directory,
    and have the prefix ``stitched_``.

    *Example:*
    
    .. code-block:: yaml
    
       stitchTiles: True
       

makeQuickLookMaps
^^^^^^^^^^^^^^^^^
    
    If True, and :ref:`nemoCommand` is set to break a large map into tiles
    (see :ref:`Tiling` below), write
    monolithic, downsampled (by a factor of 4) output filtered maps. These maps
    have the prefix ``quickLook_``.
    
    .. note::  This feature is likely to be removed in the near future.

    *Example:*
    
    .. code-block:: yaml
    
       makeQuickLookMaps: True


.. _Tiling:

Tiling
======

makeTileDir
^^^^^^^^^^^
    
    If True, break the map into tiles according to the settings specified in
    `tileDefinitions`_. For each map and mask used by :ref:`nemoCommand`,
    a directory with prefix ``tileDir_`` followed by the map filename will be
    created (see also `tileDefLabel`_ below). The individual tiles cut 
    from each map can be found within, in subdirectories labeled according to
    the tile name. This is done such that
    on subsequent :ref:`nemoCommand` runs, each process (if running in parallel
    under MPI) can read its data from disk without any bottleneck, if using
    a cluster file system (e.g., `Lustre <https://www.lustre.org/>`_).

    .. note::  If you change the definitions of the tiles, you *must* delete
               the corresponding ``tileDir_*`` directories in order for your
               changes to take effect the next time that you run
               :ref:`nemoCommand`.

    *Example:*
    
    .. code-block:: yaml
    
       makeTileDir: True
       

tileDefLabel
^^^^^^^^^^^^    
    This is used to label ``tileDir`` directories, i.e., directories will be
    created with names in the format ``tileDir_tileDefLabel_<map file name>``.
    
    *Example:*
    
    .. code-block:: yaml
    
       tileDefLabel: 'auto'

    
tileDefinitions
^^^^^^^^^^^^^^^
    
    This sets how the map is broken into tiles, which is how **Nemo** handles
    parallel processing through MPI (dividing the tiles up as equally as
    possible between processors). There are two ways in which this can be done. 
    
    The easiest way to define tiles is using the automatic tiling feature. Here,
    a mask image must be supplied, and **Nemo** will use this mask to break maps
    down into tiles with sizes as close as possible to the user-supplied
    dimensions. To use this, a dictionary must be given containing the following
    fields:
        
        :mask (str):
        
            The path to a FITS-format mask used to define the area to be tiled.
            Pixels with value 1 denote valid area, and pixels with value 0
            indicate regions that are not of interest. Masks which are compressed
            using the ``PLIO_1`` method (as implemented by 
            `astropy.io.fits <https://docs.astropy.org/en/stable/io/fits/>`_)
            are supported.
            
        :targetTileWidth (float):
            
            The desired tile width, specified in degrees. The autotiling algorithm
            will create tiles with at least this minimum size, but will increase
            the size of the tiles if needed to ensure that all tiles at some
            declination have the same width.
            
        :targetTileHeight (float):
            
            The desired tile height, specified in degrees. The autotiling algorithm
            will create tiles with at least this minimum size, but will increase
            the size of the tiles if needed to ensure that all tiles are the
            same height.
    
    *Example:*
    
    .. code-block:: yaml
    
       tileDefinitions: {mask: 'maps/Jun2020/AdvACTSurveyMask_v7_S18.fits',
                         targetTileWidthDeg: 10.0, 
                         targetTileHeightDeg: 5.0}
    
    It is also possible to set the tiling explicitly by providing a list of
    dictionaries, each containing the following keys:
        
        :tileName (str):
            
            User-supplied name for the tile. This can be anything, but should be
            unique to each tile.
        
        :RADecSection (list):
            
            This defines the coordinates of the tile, in the order minimum RA,
            maximum RA, minimum declination, maximum declination (all given in
            decimal degrees).
    
    *Example:*
    
    .. code-block:: yaml
    
       tileDefinitions:
           - {'tileName': 'tile_1',
              'RADecSection': [354.8, 6.2, -33.4, -28.1]}
           - {'tileName': 'tile_2',
              'RADecSection': [343.5, 354.8, -33.4, -28.1]}
           
    .. note:: Tiling in **Nemo** is handled such that each pixel within a tile
              is uniquely mapped to a corresponding pixel in the monolithic,
              input map.


tileOverlapDeg
^^^^^^^^^^^^^^

    This sets the overlap region (specified in degrees) between tiles, i.e., a
    buffer zone is added to all of the defined tiles (see `tileDefinitions`_), 
    increasing their width and height by the given amount. This overlap region
    is accounted for in all processing by **Nemo**, and is only used
    to ensure that filtered maps are constructed using valid data all the way
    to the edge of each tile, without any apodization within the tile area
    itself.
    
    *Example:*
    
    .. code-block:: yaml
    
       tileOverlapDeg: 1.0


tileNameList
^^^^^^^^^^^^

    If given, only the tiles named in this list will be processed. This is
    particularly useful for testing purposes.

    *Example:*
    
    .. code-block:: yaml
    
       tileNameList:
           - '1_10_7'
           - '1_10_8'      # contains J2327 (next to a source)
           - '1_11_7'


tileNoiseRegions
^^^^^^^^^^^^^^^^

    To be added - this is only used by the ``RealSpaceMatchedFilter`` method.
    

Object Detection and Photometry
===============================

thresholdSigma
^^^^^^^^^^^^^^

    This sets the minimum signal-to-noise (S/N) level for object detection. Only
    objects with S/N greater than this threshold will be included in the output
    catalog.

    .. note:: This setting is ignored if running in forced photometry mode
              (see `forcedPhotometryCatalog`_).

    *Example:*
    
    .. code-block:: yaml
    
       thresholdSigma: 4.0
 

minObjPix
^^^^^^^^^

    This sets the number of pixels above `thresholdSigma`_ that an object must
    have in order to be included in the output catalog. Higher values will be less
    susceptible to spurious noise fluctuations being detected as objects, at the
    expense of completeness.

    .. note:: This setting is ignored if running in forced photometry mode
              (see `forcedPhotometryCatalog`_).

    *Example:*
    
    .. code-block:: yaml
    
       minObjPix: 1.0


findCenterOfMass
^^^^^^^^^^^^^^^^

    If True, object properties such as position and amplitude are reported for
    the center-of-mass of the pixels above `thresholdSigma`_. Otherwise, the
    maximum pixel value is used.

    *Example:*
    
    .. code-block:: yaml
    
       findCenterOfMass: True
    

useInterpolator
^^^^^^^^^^^^^^^

    If True, subpixel interpolation will be performed using a bicubic spline.
    This affects reported source and cluster amplitudes and S/N values.
    
    *Example:*

    .. code-block:: yaml

       useInterpolator: True


rejectBorder
^^^^^^^^^^^^

    If set, a border with this width in pixels will be removed from the survey
    area and the output survey area mask will be adjusted accordingly. This can
    be used to remove "junk" at map edges, although it is better to supply a
    `surveyMask`_ that makes this unnecessary.
    
    *Example:*

    .. code-block:: yaml
    
       rejectBorder: 0


objIdent
^^^^^^^^

    Identification string that will be used as a prefix for object names reported
    in the output catalog. The rest of the name strings will be constructed from
    the object coordinates following the IAU convention, in the format
    JHHMM.m+/-DDMM, with coordinates being truncated, not rounded.
    
    *Example:*

    .. code-block:: yaml    
    
       objIdent: 'ACT-CL'


longNames
^^^^^^^^^

    If True, object names in the output catalog (see `objIdent`_) will follow the
    format JHHMMSS.s+/-DDMMSS (this may be appropriate for source rather than
    cluster catalogs).

    *Example:*

    .. code-block:: yaml 
    
       longNames: False


measureShapes
^^^^^^^^^^^^^

    If True, object shapes are estimated from the 2nd-order moments of the pixels
    above `thresholdSigma`_, using the method employed by
    `Source Extractor <https://www.astromatic.net/software/sextractor>`_. 
    The shape parameters are recorded in columns with the prefix ``ellipse_`` in
    the output catalog. The ellipse major and minor axis lengths, position angle,
    center coordinates, and eccentricty are recorded.
    
    This feature can be used to help identify and mask extended sources in the
    maps.
    
    .. note:: All size information is given in pixel units, i.e., there is no
              correction applied here for distortion due to the map projection.
              Objects must contain a minimum of 9 pixels with S/N > `thresholdSigma`_
              in order for their shape measurements to be included in the catalog.

    *Example:*

    .. code-block:: yaml 
    
       measureShapes: True


removeRings
^^^^^^^^^^^

    If True, **Nemo** will attempt to identify spurous ring-features in the map
    that can occur in the vicinity of extremely bright sources (or even clusters
    in some circumstances). This is done by performing object segmentation down
    to the level set by `ringThresholdSigma`_, and testing whether the central
    coordinates of the object are found within significant pixels. If this is
    not the case, the object is deemed to be a ring-like feature and is excluded
    from the output catalog.

    *Example:*

    .. code-block:: yaml 
    
       removeRings: True
       

ringThresholdSigma
^^^^^^^^^^^^^^^^^^

    This sets the threshold used by the algorithm that finds and excludes
    spurious ring-like features (see `removeRings`_).

    *Example:*

    .. code-block:: yaml 
    
       ringThresholdSigma: 3


undoPixelWindow
^^^^^^^^^^^^^^^

    If True (the default), all signal amplitudes reported in the output catalog
    will be corrected for the pixel window effect.

    *Example:*

    .. code-block:: yaml 
    
       undoPixelWindow: True

       
forcedPhotometryCatalog
^^^^^^^^^^^^^^^^^^^^^^^
    
    Path to a file containing an object catalog. If this parameter is set,
    :ref:`nemoCommand` will perform forced photometry at the object coordinates,
    rather than detecting objects. The coordinates should be RA and declination
    in decimal degrees, and contained in columns named ``ra``, ``RA``, or ``RADeg``,
    and ``dec``, ``DEC``, ``decDeg``, or ``Dec``.
   
    .. note:: If this mode is used, `thresholdSigma`_ and `minObjPix`_ will be
              ignored. 
    
    .. note:: Forced photometry can be performed by both :ref:`nemoCommand`
              and :ref:`nemoMassCommand` (see :ref:`AdvancedForcedPhotometry`).
    
    *Example:*

    .. code-block:: yaml 
       
       forcedPhotometryCatalog: "redMaPPer.fits"


twoPass
^^^^^^^

    If True, the ``twoPass`` source detection pipeline is used. This runs an
    initial search for extremely bright sources in the first pass using a
    simple noise model, so that these sources can be subtracted when using the
    map itself as the filter noise model in the second pass. This allows
    source catalogs to be constructed from the maps with zero (or very little)
    masking.
    
    .. note::  At present, this method works only for detecting sources, and
               not clusters. An equivalent cluster pipeline can be constructed
               by using two :ref:`nemoCommand` runs, using the output of the
               first run as `noiseMaskCatalog`_ in the second run).

    *Example:*

    .. code-block:: yaml 
    
       twoPass: False


catalogCuts
^^^^^^^^^^^
    
    A list of constraints that can be used to filter output catalogs, where each
    item is a string of the form  "key < value", "key > value", etc., and "key"
    can be any column name in the output catalog. Note that the spaces between
    the key, operator (e.g. '<'), and value are essential.

    *Example:*

    .. code-block:: yaml 
    
       catalogCuts: ['fluxJy > 0.01']
       
    
Filters
=======
mapFilters
allFilters
GNFWParams

Cluster Mass Estimates
======================
photFilter
fitQ
massOptions

Selection Function
==================
calcSelFn
selFnOptions

**selFnFootprints** (list):
    
    Blah
    
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
