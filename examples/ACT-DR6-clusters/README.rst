.. note::  This is under construction. It will probably only be
           completed after a new **Nemo** release is made together
           with an updated release of the ACT DR6 cluster search
           products. So, fully up-to-date documentation might not
           be found in the **Nemo** source distribution
           corresponding to the latest stable release.

Here we provide a guide to reproducing the ACT DR6 cluster search
data products, as described in the `ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_.

.. note::  The config files and scripts needed for this tutorial can be
           found in the `examples/ACT-DR6-clusters <https://github.com/simonsobs/nemo/tree/main/examples/ACT-DR6-clusters>`_
           directory of the **Nemo** source distribution.


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

To get **Nemo** to produce information related to completeness for
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


Making model cluster signal maps
================================



Forced photometry
=================



Running source injection simulations
====================================

