Here we provide a guide to reproducing the 
`ACT DR5 cluster catalog <https://ui.adsabs.harvard.edu/abs/2020arXiv200911043H/abstract>`_.

.. note::  The config files and scripts needed for this tutorial can be
           found in the `examples/ACT-DR5-clusters <https://github.com/simonsobs/nemo/tree/main/examples/ACT-DR5-clusters>`_
           directory of the **Nemo** source distribution.

.. note::  The `released ACT DR5 cluster catalog v1.1 <https://lambda.gsfc.nasa.gov/product/act/actpol_dr5_szcluster_catalog_get.html>`_
           was produced with code that is closest (but not identical) to
           `Nemo v0.5.0 <https://github.com/simonsobs/nemo/releases/tag/v0.5.0>`_.
           Later versions of **Nemo** may produce slightly different output as the
           code is improved and bugs are fixed.


Downloading maps and other needed files
=======================================

The first step is to download the ACT DR5 maps, described in 
`Naess et al. (2020) <https://ui.adsabs.harvard.edu/abs/2020JCAP...12..046N/abstract>`_,
from `LAMBDA <https://lambda.gsfc.nasa.gov/product/act/actpol_prod_table.cfm>`_. In
addition to this, we require some mask files, point source catalogs, and the beam
profiles. These can all be downloaded and extracted by running a script:
    
.. code-block::

   sh FETCH_MAPS.sh
   
Each map is 5 Gb in size, so this may take a while. The archive named
``nemo-dr5-masking.tar.gz`` downloaded by this script contains the masks, point
source catalogs, beam profiles, and a catalog of cluster redshifts
(see `Inferring cluster masses`_ below).


Extracting the cluster candidate list
=====================================

The DR5 maps are large, and so **Nemo** breaks these maps into tiles, which can
be processed separately in parallel, before combining the output to make the
final catalog and filtered maps. A script (``DR5Clusters.slurm``) is provided
that will run :ref:`nemoCommand` in this mode using the
`Slurm Workload Manager <https://slurm.schedmd.com/overview.html>`_ :
    
.. code-block::
    
   #!/bin/sh
   #SBATCH --nodes=18
   #SBATCH --ntasks-per-node=16
   #SBATCH --mem=64000
   #SBATCH --time=03:59:00
   #SBATCH --output=DR5ClusterSearch.log
   #SBATCH --error=DR5ClusterSearch.err
   
   time mpiexec nemo DR5ClusterSearch.yml -M 

Note the  :ref:`nemoCommand` ``-M`` switch is needed to enable parallel
processing using MPI. The **Nemo** configuration for the ACT DR5 cluster
search can be found in ``DR5ClusterSearch.yml``.

After this job is completed, the cluster candidate list, filtered maps, and
various other files will be found under the ``DR5ClusterSearch`` directory.
For more information on what these are, see :ref:`Output`.


Inferring cluster masses
========================

To infer cluster masses, a redshift catalog is needed. This is provided
for the ACT DR5 clusters in the file ``redshifts_confirmed_v1.1.fits``,
which is included in the ``nemo-dr5-masking.tar.gz`` archive downloaded
by the ``FETCH_MAPS.sh`` script.

The :ref:`nemoMassCommand` tool can be used to estimate cluster masses:
    
.. code-block::

   mpiexec nemoMass DR5ClusterSearch.yml -M

This should take a few minutes to run, and will write a FITS-table
format catalog to ``DR5ClusterSearch/DR5ClusterSearch_mass.fits``.

The inferred masses depend on the assumed scaling relation parameters,
which can be edited in the ``DR5ClusterSearch.yml`` config file.
