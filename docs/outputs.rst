.. _Output:

======
Output
======

This section describes the output produced by :ref:`nemoCommand`, and the main
files of interest.


Directory Structure
-------------------

The main :ref:`nemoCommand` command creates a folder in the current working directory
under which it stores all of its output. By default this directory has the same name
as the configuration file, but with the ``.yml`` extension removed. We will refer to
this as the ``output`` directory hereafter. 

After :ref:`nemoCommand` has run, within the ``output`` directory, you will find the
following::

    output/
        diagnostics/
        filteredMaps/
        selFn/
        output_optimalCatalog.csv
        output_optimalCatalog.fits
        output_optimalCatalog.reg

diagnostics/
    This directory contains output that is useful for debugging (mainly plots in ``.pdf`` and ``.png`` 
    format, but the filter kernels are also cached here).

filteredMaps/
    This directory contains the filtered maps, named according to the filter definitions given in the
    configuration file. If you have chosen to break up a large map into tiles for parallel processing
    (see :ref:`Tiling`), then this directory will contain subdirectories named after each tile. If
    ``stitchTiles: True`` is set in the config file (see :ref:`OutputOptions`), then you will find
    monolithic maps named ``stitched_*.fits`` in this directory. For more on the maps themselves,
    see :ref:`FilteredMaps` below.

selFn/
    This directory contains all the necessary information for calculating the completeness of a
    cluster search, given a model for the relationship between the SZ signal and mass, and an
    assumed halo mass function. For both point source and cluster searches, you will find the survey
    area masks and noise maps in this directory.

output_optimalCatalog.fits, output_optimalCatalog.csv
    These files contain the object catalogs. Both contain the same information, but in different
    formats. The CSV file is plain text, but despite the extension is tab-delimited. Largely this is
    provided for 'historical' reasons, as the FITS table format file is more convenient (it can be
    read by both `TopCat <http://www.star.bris.ac.uk/%7Embt/topcat/>`_ and
    `Astropy Table <https://docs.astropy.org/en/stable/table/index.html>`_ without any special handling).
    For information on the contents of these catalogs, see :ref:`Catalogs`.

output_optimalCatalog.reg
    This is a region file that can be loaded into `SAOImage DS9 <https://sites.google.com/cfa.harvard.edu/saoimageds9/home>`_
    to overplot the object locations on any FITS image.

.. _Catalogs:
    
Catalogs
--------

The catalogs produced by :ref:`nemoCommand` may contain the following columns:
    
name
    The name of the object.
    
RADeg
    The right ascension of the object in decimal degrees (J2000 typically, but this depends upon the
    WCS of the maps used with :ref:`nemoCommand`).

decDeg:
    The declination of the object in decimal degrees (J2000 typically, but this depends upon the
    WCS of the maps used with :ref:`nemoCommand`).

SNR:
    Signal-to-noise ratio, optimized over all filter scales.

fixed_SNR:
    Signal-to-noise ratio at the reference filter scale, chosen using ``photFilterLabel`` in the config file
    (see :ref:`ClusterMassEstimates`).
    
numSigPix:
    The number of pixels above ``thresholdSigma`` as set in the config file (see :ref:`Detection`).

template:
    Name of the matched filter template resulting in the highest SNR detection of this object.

tileName:
    The name of the tile in which the object was detected.
    
galacticLatDeg:
    Galactic latitude of the object in decimal degrees (can be used for pruning objects close to the
    Galactic plane, where contamination is naturally expected to be higher).

deltaT_c:
    The amplitude of the source, in ΔTemperature (μK) with respect to the CMB.

err_deltaT_c:
    The 1-sigma uncertainty on the amplitude of the source, in ΔTemperature (μK) with respect to the CMB.

fluxJy:
    The flux density of the source in Jy. Note this quantity is derived from ``deltaT_c``, and depends
    upon the beam solid angle and the assumed effective frequency for the map.

err_fluxJy:
    The 1-sigma uncertainty on the flux density of the source, in Jy.

y_c:
    Central Comptonization parameter (10\ :sup:`-4`) measured using the optimal matched filter template
    (i.e., the one that maximizes SNR).

err_y_c:
    The 1-sigma uncertainty on the central Comptonization parameter (10\ :sup:`-4`) measured using the
    optimal matched filter template (i.e., the one that maximizes SNR).

fixed_y_c:
    Central Comptonization parameter (10\ :sup:`-4`) measured at the reference filter scale, chosen
    using ``photFilterLabel`` in the config file (see :ref:`ClusterMassEstimates`). This quantity is
    called ỹ\ :sub:`0` in the ACT cluster papers. 

err_fixed_y_c:
    The 1-sigma uncertainty on the central Comptonization parameter (10\ :sup:`-4`) measured at the
    reference filter scale, chosen using ``photFilterLabel`` in the config file
    (see :ref:`ClusterMassEstimates`).

tileBoundarySplit:
    If ``True``, this object *may* have been de-blended across a tile boundary. At present, the only way
    to determine if this really is the case is by visual inspection of the maps and catalogs (it happened
    rarely when using the ACT DR5 maps).


.. _FilteredMaps:
    
Filtered Maps
-------------

The :ref:`nemoCommand` command produces two different types of filtered maps, assuming that you have
configured :ref:`nemoCommand` to write filtered maps to disk (see :ref:`Filters`):

\*_filteredMap.fits:
    These are maps in signal units, set by ``outputUnits`` in the config file (see :ref:`Filters`).
    For point source searches, these will be ΔTemperature (μK) with respect to the CMB
    (if ``outputUnits: uK`` is set in the config file). For cluster searches, these will be in terms
    of the central Comptonization parameter, y\ :sub:`0` (if ``outputUnits: yc`` is set in the config
    file).

\*_SNMap.fits:
    These are maps in signal-to-noise (S/N) units.
    
