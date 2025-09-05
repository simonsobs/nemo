Here we provide a guide to reproducing the ACT DR6 cluster search
data products, as described in the `ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_.

.. note::  The config files and scripts needed for this tutorial can be
           found in the `examples/ACT-DR6-clusters <https://github.com/simonsobs/nemo/tree/main/examples/ACT-DR6-clusters>`_
           directory of the **Nemo** source distribution.

.. note::  Fully up-to-date documentation might not
           be found in the **Nemo** source distribution
           corresponding to the latest stable release.


Downloading maps and other needed inputs
========================================

The first step is to download the co-added ACT+Planck maps used for
the ACT DR6 cluster search. These are described in
`Naess et al. (2025) <https://ui.adsabs.harvard.edu/abs/2025arXiv250314451N/abstract>`_,
and available on `LAMBDA <https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/>`_.
In addition to this, we require the beam transform files and some
masks. These can be downloaded and extracted by running a script:

.. code-block::

   sh FETCH_MAPS.sh

Each map is 5 Gb in size, so this may take a while.


Extracting the cluster candidate list
=====================================

**Nemo** breaks large maps into tiles, which can be processed
separately in parallel, before combining the output to make
the final catalogs and filtered maps. This requires the
:ref:`nemoCommand` command to run with the ``-M`` switch as below:

.. code-block::

   mpiexec nemo DR6ClusterSearch.yml -M

The **Nemo** configuration for the ACT DR6 cluster search is
found in ``DR6ClusterSearch.yml``. See :ref:`ConfigReference`
for a description of the available options and parameters.
This particular configuration finds point sources and clusters
in the ACT maps using a multi-pass approach, as described in
Section 2 of the
`ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_.
Point sources are subtracted from the maps before performing
the cluster search and affected pixels are flagged.

The file ``DR6ClusterSearch.pbs`` contains an example job script
for running **Nemo** on a compute cluster. The DR6 cluster search
takes 03h40m to complete, running on 144 CPU cores on Lengau at
the `CHPC <https://www.chpc.ac.za/>`_.

After this job is completed, the cluster candidate list,
filtered maps, and various other files will be found under the
``DR6ClusterSearch`` directory. For more information on what
these are, see :ref:`Output`.


Post-processing runs (footprint and flagging adjustments)
=========================================================

**Nemo** can be run again on existing output without re-running
the full map filtering and object detection pipeline, after adjusting
some of the options in the configuration file.

To produce information related to completeness for
different footprints that overlap with the cluster search area, run
:ref:`nemoCommand` with the ``-S`` switch as below:

.. code-block::

   nemo DR6ClusterSearch.yml -S

This does not currently benefit from being run in parallel, and
takes 1h12m to complete on one of the Lengau compute nodes at the
`CHPC <https://www.chpc.ac.za/>`_. The file
``DR6ClusterSearch-S.pbs`` contains a job script for running
:ref:`nemoCommand` in this way.

The ``selFnOptions`` and ``selFnFootprints`` keywords in the
``DR6ClusterSearch.yml`` configuration file control which regions
completeness statistics are calculated for. The information on
mass limits depends upon the assumed scaling relation parameters
set in ``massOptions`` (see :ref:`ConfigReference`).
After running :ref:`nemoCommand` with the ``-S`` switch,
plots corresponding to Figures 8, 9 and 10 of the
`ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_.
will be found
under the ``DR6ClusterSearch/diagnostics`` directory.

Post-processing flags can also be adjusted by re-running
:ref:`nemoCommand` after editing the ``postFlags`` keyword in the
``DR6ClusterSearch.yml`` configuration file. This allows masks
that are not related to object finding to be changed, such as the
dust mask, and regions flagged around extended objects (e.g. nearby
galaxies) and bright stars.


Estimating cluster masses
=========================

To infer cluster masses, a redshift catalog is needed. This is provided
for the ACT DR6 clusters in the file ``redshifts.fits``,
which is included in the ``ACTDR6ClustersExtras.tar.gz`` archive downloaded
by the ``FETCH_MAPS.sh`` script.

The :ref:`nemoMassCommand` tool can be used to estimate cluster masses:

.. code-block::

   mpiexec nemoMass DR6ClusterSearch.yml -M

This should take a few minutes to run, and will write a FITS-table
format catalog to ``DR6ClusterSearch/DR6ClusterSearch_mass.fits``.

The inferred masses depend on the assumed scaling relation
parameters, which can be edited in the ``massOptions`` section of the
``DR6ClusterSearch.yml`` config file.

Running ``nemoMass`` with the ``-I`` switch adds some additional columns
for inferred cluster SZ properties (e.g. ``inferred_y_c`` and
``inferred_Y500Arcmin2``; see :ref:`Catalogs`).


Making model cluster signal maps
================================

To make the model cluster signal maps as provided with the DR6 cluster
search data products, it is easiest to run:

.. code-block::

   mpiexec nemoMass DR6ClusterSearch.yml -M -I -m

This will produce cluster model maps for clusters with redshifts only,
writing them to the directory ``DR6ClusterSearch/clusterModelMaps/``.

To produce model maps from the candidate list, you can do something like
the following:

.. code-block::

   nemoModel DR6ClusterSearch/DR6ClusterSearch_optimalCatalog.fits masks/ExtendedSurveyMask_v3.fits beams/beam_f090_tform.txt candidatesModelMap_f090.fits -f 97.8 -m
   nemoModel DR6ClusterSearch/DR6ClusterSearch_optimalCatalog.fits masks/ExtendedSurveyMask_v3.fits beams/beam_f150_tform.txt candidatesModelMap_f150.fits -f 149.6 -m
   nemoModel DR6ClusterSearch/DR6ClusterSearch_optimalCatalog.fits masks/ExtendedSurveyMask_v3.fits beams/beam_f220_tform.txt candidatesModelMap_f220.fits -f 216.5 -m

These should take roughly 3 minutes to run per frequency map. The above
set-up will include objects that are flagged (``flags > 0``) in the
candidates list, as the whole candidate list is fed into
:ref:`nemoModelCommand` in this example.

With the default settings, the output maps may look boxier than you may like,
as cluster model images are painted inside a region of size
:math:`10 \times \theta`\ :sub:`500c`.
This can be eliminated by adding ``-x 40`` to the list of arguments given
to :ref:`nemoModelCommand`. This paints cluster model images out to
:math:`40 \times \theta`\ :sub:`500c`, at the cost of additional run time.

Note that running :ref:`nemoModelCommand` requires a minimum of 16 GB
of memory for DR6-sized maps.


Forced photometry
=================

To run forced photometry on eRASS1 cluster positions as described in
Section 5.2 of the
`ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_,
first download the appropriate eRASS1 cluster catalog from
`this page <https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/>`_,
and then run:

.. code-block::

   nemoMass DR6ClusterSearch.yml -c erass1cl_primary_v3.2.fits -F

This will create a catalog called ``erass1cl_primary_v3.2_mass.fits`` in
the current directory. This contains measurements extracted by forced
photometry at the positions given in the ``erass1cl_primary_v3.2.fits``
catalog, found in the ``fixed_SNR`` and ``fixed_y_c`` columns, in
addition to SZ mass estimates inferred using the scaling relation
set-up specified in ``DR6ClusterSearch.yml``.

A subset of the eRASS1 forced photometry catalog produced in the above
way is plotted in Figure 18 of the
`ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_.

The forced photometry mode can be used to do stacking analyses of
external catalogs by averaging the output SZ quantities in bins of
some other observable (e.g., optical richness). Other tools are
available for stacking in image space - see for example the
`DR6 y-map stacking notebook <https://github.com/ACTCollaboration/DR6_Notebooks/blob/main/ACT_DR6_ymap_stacking.ipynb>`_.


Running source injection simulations
====================================

Source injection simulations can be run using:

.. code-block::

   mpiexec nemo DR6ClusterSearch.yml -M -I

This was used to produce the results shown in Section 2.4 and Figure 7 of the
`ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_.

This takes 12h to complete, running on 144 CPU cores on Lengau at
the `CHPC <https://www.chpc.ac.za/>`_ using the set-up given in
``DR6ClusterSearch.yml``. See :ref:`SourceInjection` for documentation on
configuration file parameters that affect this mode.


Making random catalogs
======================

The random catalogs included in the DR6 cluster search products were
generated by running:

.. code-block::

   nemoMock DR6ClusterSearch/selFn mocks -N 100 -p 2.1 -C

This includes the optimization bias model as described in Appendix B
of the
`ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_.
See the documentation for :ref:`nemoMockCommand` for more
information on arguments for this command. Note that the output of
:ref:`nemoMockCommand` depends upon the assumed cosmology and scaling
relation parameters set in the configuration file.
