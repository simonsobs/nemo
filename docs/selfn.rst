.. _SelFnPage:

============================================================
Tutorial: Using the Selection Function for Cluster Cosmology
============================================================

Nemo stores outputs related to cluster selection in the ``output/selFn`` directory
(see :ref:`Output`). Here we show how to use the selection function outputs produced
after running Nemo on the Simons Observatory Simulations (see :ref:`SOSimsPage`).
You can find the files referred to here in the 
`examples/SOSims <https://github.com/simonsobs/nemo/tree/master/examples/SOSims>`_
directory in the source code distribution.

Small Maps
----------

If you have not already run Nemo using the ``MFMF_SOSim_3freq_small.yml`` config
file (see :ref:`SOSimsPage`), you can download the catalogs and selection function
files using::

    wget https://acru.ukzn.ac.za/~mjh/SO/releases/MFMF_SOSim_3freq_small.tar.gz --user=user --password=password

.. note::
    
    You can find the needed username and password on the `SO Wiki <http://simonsobservatory.wikidot.com/awg:sz>`_. 

Unpacking this archive will create a directory called ``MFMF_SOSim_3freq_small`` that
contains the ``selFn`` directory within it, as well as the catalog files.

The ``selFnExample.py`` script shows how to use the selection function files with 
Nemo routines to calculate the completeness as a function of mass and redshift, given
a set of cosmological and scaling relation parameters. You can run it like this::

    python selFnExample.py MFMF_SOSim_3freq_small/selFn

This will output the number of clusters that are expected to be found in the 600 
square degrees simulated map used here, for various parameter combinations. It also
produces ``.pdf`` plots of the completeness level as a function of mass and redshift
(these are named ``test0_Mz.pdf`` etc.). The script accepts a number of other 
command-line arguments; use ``python selFnExample.py -h`` to see them.

The ``selFn`` directory contains all that you need to perform completeness calculations
and generate mock catalogs. For example, you can generate 20 mocks from the small SO
sim maps using::
    
    nemoMock MFMF_SOSim_3freq_small/selFn mocks_small -N 20

The `nemoCosmo <https://github.com/simonsobs/nemo/blob/master/bin/nemoCosmo>`_ script 
shows how the selection function information can be used in a cosmological likelihood
code. Note that this script is not an attempt at writing a full-blown publication-ready
likelihood code - it is used for testing/sanity checking the completeness calculations
and mock catalogs produced by Nemo. You will need to install 
`Cobaya <https://cobaya.readthedocs.io/en/latest/index.html>`_ to use ``nemoCosmo``.

The ``nemoCosmo`` script requires a catalog to run - you can obtain an archive full of mock 
catalogs generated from the ``MFMF_SOSim_3freq_small.yml`` configuration from 
`here <https://acru.ukzn.ac.za/~mjh/SO/releases/mocks_small.tar.gz>`_ 
(again, you will need the username and password `from the wiki <http://simonsobservatory.wikidot.com/awg:sz>`_).
These mock catalogs were generated using the ``nemoMock`` script, as described above.

You can run ``nemoCosmo`` like this::

    nemoCosmo mocks_small/mockCatalog_1.fits MFMF_SOSim_3freq_small/selFn -m 10
    
The ``-m`` switch is used to set the maximum number of samples (per process; set to 10 in this 
case for testing; the default value is 3000). You can use ``nemoCosmo -h`` to see the other available 
options. For this example, you will find output (chains) under the ``cosmo_mockCatalog_1_5.00/`` directory 
after ``nemoCosmo`` completes its run. You can analyse the chains with 
`GetDist <https://getdist.readthedocs.io/en/latest/>`_ (although that is pointless with only 10 samples).

The parameters to be estimated using ``Cobaya`` can be controlled using a ``.yml`` configuration file. You
can dump the default settings to a file named ``default.yml`` in the current working directory using the 
``-d`` switch::
    
    nemoCosmo mocks_small/mockCatalog_1.fits MFMF_SOSim_3freq_small/selFn -d

The contents of ``default.yml`` should look something like this::
    
    params:
    Ob0: 0.05
    tenToA0: 4.95e-05
    B0: 0.08
    Mpivot: 300000000000000.0
    sigma_int: 0.2
    H0:
        latex: H_0
        prior:
        dist: norm
        loc: 70.0
        scale: 4.0
        proposal: 5.0
    Om0:
        latex: \Omega_{\rm m0}
        prior:
        max: 0.5
        min: 0.1
        proposal: 0.05
        ref:
        dist: norm
        loc: 0.3
        scale: 0.1
    sigma8:
        latex: \sigma_8
        prior:
        max: 0.9
        min: 0.6
        proposal: 0.02
        ref:
        dist: norm
        loc: 0.8
        scale: 0.1
    sampler:
    mcmc:
        burn_in: 50
        max_samples: 3000
        max_tries: .inf

You can pass ``nemoCosmo`` your own Cobaya configuration file by using the ``-c`` switch.

Since ``Cobaya`` is used as the sampler, you can easily run ``nemoCosmo`` using MPI by doing, e.g.,::
    
    mpiexec nemoCosmo mocks_small/mockCatalog_1.fits MFMF_SOSim_3freq_small/selFn

The `examples/SOSims/slurm_cosmo.sh <https://github.com/simonsobs/nemo/blob/master/examples/SOSims/slurm_cosmo.sh>`_
file shows how to run ``nemoCosmo`` using the Slurm job scheduler. Running on 80 processors, 
it should take approximately two hours for the MCMC run to converge, when using the small mock catalogs 
made from the SO sims.

If you want to run on a "real" rather than mock catalog, the ``nemoCosmo`` script understands the 
format made by the ``nemoMass`` script. Note that although these catalogs contain mass estimates computed by 
Nemo for a fixed set of cosmological and scaling relation parameters, these are 
not used by ``nemoCosmo`` itself (it just needs the cluster positions, SZ-signal measurements, and
the redshifts). For this example, you can find the catalog here: 
``MFMF_SOSim_3freq_small/MFMF_SOSim_3freq_small_M500.fits``.

Tiled Maps
----------

Nemo breaks up large maps into tiles. At the moment, the time needed to run the selection 
function routines scales with the number of tiles - 540 in the case of the 
``MFMF_SOSim_3freq_tiles.yml`` configuration. So, while you can run the ``selFnExample.py``
and ``nemoCosmo`` scripts on the ``MFMF_SOSim_3freq_tiles.yml`` configuration, it will be
slower - this is being worked on.

You can find an archive containing the catalogs and selection function files for a full run
on the nominal SO survey area (~17,000 square degrees) 
`here <https://acru.ukzn.ac.za/~mjh/SO/releases/MFMF_SOSim_3freq_tiles.tar.gz>`_.

Note that the SO simulations have a different mass-scaling relation to the default 
Arnaud et al. (2010) / Universal Pressure Profile derived scaling relation assumed in the
example Nemo config files.
