.. _SelFnPage:

======================================
Tutorial: Using the Selection Function
======================================

Nemo stores outputs related to cluster selection in the ``output/selFn`` directory
(see :ref:`Output`). Here we show how to use the selection function outputs produced
after running Nemo on the Simons Observatory Simulations (see :ref:`SOSimsPage`).
You can find the files referred to here in the 
`examples/SOSims <https://github.com/simonsobs/nemo/tree/master/examples/SOSims>`_
directory in the source code distribution.

Small Maps
----------

If you have not already run Nemo using the ``MFMF_SOSim_3freq_small.yml`` config
file (see :ref:`SOSimsPage`), you can download the selection function files using::

    wget https://acru.ukzn.ac.za/~mjh/SO/selFn_MFMF_SOSim_3freq_small.tar.gz --user=user --password=password

.. note::
    
    You can find the needed username and password on the `SO Wiki <http://simonsobservatory.wikidot.com/awg:sz>`_. 

Unpack this archive in the same directory as where the ``MFMF_SOSim_3freq_small.yml`` 
config file is found.

The ``selFnExample.py`` script shows how to use the selection function files with 
Nemo routines to calculate the completeness as a function of mass and redshift, given
a set of cosmological and scaling relation parameters. You can run it like this::

    python selFnExample.py MFMF_SOSim_3freq_small.yml selFn_MFMF_SOSim_3freq_small

This will output the number of clusters that are expected to be found in the 600 
square degrees simulated map used here, for various parameter combinations. It also
produces ``.pdf`` plots of the completeness level as a function of mass and redshift
(these are named ``test0_Mz.pdf`` etc.).

The `nemoCosmo <https://github.com/simonsobs/nemo/blob/master/bin/nemoCosmo>`_ script 
shows how the selection function information can be used in a cosmological likelihood
code. Note that this script is not an attempt at writing a full-blown publication-ready
likelihood code - it is used for testing/sanity checking the completeness calculations
and mock catalogs produced by Nemo. 

The ``nemoCosmo`` script requires a catalog to run - 
you can obtain an archive full of mock catalogs generated from the ``MFMF_SOSim_3freq_small.yml`` 
configuration from `here <https://acru.ukzn.ac.za/~mjh/SO/mocks_MFMF_SOSim_3freq_small.tar.gz>`_ 
(again, you will need the username and password `from the wiki <http://simonsobservatory.wikidot.com/awg:sz>`_).
These mock catalogs were generated using the ``nemoMock`` script.

You can run ``nemoCosmo`` like this::

    nemoCosmo MFMF_SOSim_3freq_small.yml mockCatalog_1.fits selFn_MFMF_SOSim_3freq_small 100 10
    
The last two numbers here give the number of steps per walker (100 in this case) and the number
of walkers (10), used by the `emcee <https://emcee.readthedocs.io/en/stable/>`_ MCMC package. 
By default ``nemoCosmo`` will run in parallel using Python's multiprocessing - but other 
options are available, including a serial mode (use ``nemoCosmo -h`` to see the available options).
For this example, you will find output (a corner plot and the chains) under the 
``cosmo_mockCatalog_1_5.00/`` directory after ``nemoCosmo`` completes its run.

If you want to run on a "real" rather than mock catalog, the ``nemoCosmo`` script understands the 
format made by the ``nemoMass`` script. Note that although these catalogs contain mass estimates computed by 
Nemo for a fixed set of cosmological and scaling relation parameters, these are 
not used by ``nemoCosmo`` itself (it just needs the cluster positions, SZ-signal measurements, and
the redshifts). You can obtain the cluster catalog extracted by Nemo from the
SO sims map from `here <https://acru.ukzn.ac.za/~mjh/SO/MFMF_SOSim_3freq_small_M500.fits>`_.

Tiled Maps
----------

Nemo breaks up large maps into tiles. At the moment, the time needed to run the selection 
function routines scales with the number of tiles - 540 in the case of the 
``MFMF_SOSim_3freq_tiles.yml`` configuration. So, while you can run the ``selFnExample.py``
and ``nemoCosmo`` scripts on the ``MFMF_SOSim_3freq_tiles.yml`` configuration, it will be very
slow - this is being worked on.

Here are some links for output produced by Nemo using the ``MFMF_SOSim_3freq_tiles.yml``
configuration, if you would like to experiment with them:
    
    * `Selection function files <https://acru.ukzn.ac.za/~mjh/SO/selFn_MFMF_SOSim_3freq_tiles.tar.gz>`_
    * `Mock catalogs <https://acru.ukzn.ac.za/~mjh/SO/mocks_MFMF_SOSim_3freq_tiles.tar.gz>`_
    * `Cluster catalog extracted from the SO maps <https://acru.ukzn.ac.za/~mjh/SO/MFMF_SOSim_3freq_tiles_M500.fits>`_

Note that the SO simulations have a different mass-scaling relation to the default 
Arnaud et al. (2010) / Universal Pressure Profile derived scaling relation assumed in the
example Nemo config files.
