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
file (see :ref:`SOSimsPage`), you can download the catalogs and selection function
files using::

    wget https://acru.ukzn.ac.za/~mjh/SO/MFMF_SOSim_3freq_small.tar.gz --user=user --password=password

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

    nemoCosmo MFMF_SOSim_3freq_small/selFn/config.yml MFMF_SOSim_3freq_small/mocks/mockCatalog_1.fits MFMF_SOSim_3freq_small/selFn 100 10
    
The last two numbers here give the number of steps per walker (100 in this case for testing; 
2000 would be more reasonable) and the number of walkers (10), used by the 
`emcee <https://emcee.readthedocs.io/en/stable/>`_ MCMC package. By default ``nemoCosmo`` will
run in parallel using Python's multiprocessing - but other options are available, including a serial
mode (use ``nemoCosmo -h`` to see the available options). For this example, you will find output
(a corner plot and the chains) under the ``cosmo_mockCatalog_1_5.00/`` directory after
``nemoCosmo`` completes its run.

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
slow - this is being worked on.

You can find an archive containing the catalogs and selection function files for a full run
on the nominal SO survey area (~17,000 square degrees) 
`here <https://acru.ukzn.ac.za/~mjh/SO/MFMF_SOSim_3freq_tiles.tar.gz>`_.

Note that the SO simulations have a different mass-scaling relation to the default 
Arnaud et al. (2010) / Universal Pressure Profile derived scaling relation assumed in the
example Nemo config files.
