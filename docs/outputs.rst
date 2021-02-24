.. _Output:

======
Output
======

This section describes the output produced by Nemo, and the main files of interest.

Directory Structure
-------------------

The main Nemo script (launched with the ``nemo`` command) creates a folder (by default sharing the name
of the configuration file, but with ``.yml`` removed) in the current working directory under which it
stores all of its output. We will refer to this as the ``output`` directory hereafter. 

After ``nemo`` has run, within the ``output`` directory, you will find the following::

    output/
        diagnostics/
        filteredMaps/
        mocks/
        selFn/
        output_optimalCatalog.csv
        output_optimalCatalog.fits

diagnostics/
    This directory contains output that is useful for debugging (mainly plots in ``.pdf`` and ``.png`` 
    format, but the filter kernels are also cached here).

filteredMaps/
    This directory contains the filtered maps, named according to the filter definitions given in the
    configuration file. For more on the maps themselves, see :ref:`FilteredMaps` below.

mocks/
    If you run the ``nemoMock`` script, numbered mock catalogs will be written into this directory 
    (otherwise it will be empty).

selFn/
    This directory contains all the necessary information for calculating the completeness of the 
    survey, given a model for the relationship between the SZ signal and mass, and an assumed halo 
    mass function. For more information, see :ref:`SelFnPage` and/or :ref:`OutputSelFn` below.

.. _Catalogs:
    
Catalogs
--------

And region files.
M500 too.

.. _FilteredMaps:
    
Filtered Maps
-------------

.. _OutputSelFn:
    
Selection Function
------------------

Diagnostics
-----------
