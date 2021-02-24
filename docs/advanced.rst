.. _Scripts:

===============
Advanced Topics
===============

Here we describe how to use some selected features of **Nemo** that are not
necessarily covered in the :ref:`TutorialsIndex`.


.. _AdvancedForcedPhotometry:
    
Forced photometry
=================

Forced photometry is the extraction of flux or SZ-signal measurements at fixed prior positions - i.e.,
an external catalog is used to supply a list of object positions, instead of using Nemo to detect the
objects. Nemo offers two modes in which forced photometry is performed, depending on whether you are
interested in performing forced photometry on sources or clusters.

If you have already run the ``nemo`` command and produced an output directory containing filtered maps 
on which you wish to perform cluster photometry, the ``nemoMass`` command can be used. For example,::

    nemoMass S18d_202003.yml -c externalCatalogChecks/mergeMaDCoWS/MADCOWSUnion.fits

runs ``nemoMass`` using the configuration file ``S18d_202003.yml`` on the external catalog specified using
the ``-c`` switch (in this case, ``externalCatalogChecks/mergeMaDCoWS/MADCOWSUnion.fits``, which contains
the union of Tables 4, 5, and 6 in the Gonzalez et al. (2019) MaDCoWS galaxy cluster catalog paper).
This will extract the SZ signal at the location of each object in the catalog that lies in the 
filtered map. For objects that have redshift estimates in the catalog, the mass will be estimated 
using the scaling relation parameters given in the configuration file 
(``S18d_202003.yml`` in this example). By default, ``nemoMass`` will look for a column named 
``redshift``, but a few possible variations are detected, e.g., ``Photz``, as in this example. 
You can also specify the name of the redshift column using the ``-z`` switch.

The output will be written to a file with the same name as the
catalog specified by the ``-c`` switch, but with ``_M500.fits`` appended 
(``externalCatalogChecks/mergeMaDCoWS/MADCOWSUnion_M500.fits`` for the example above). One could 
perform a stacking analysis to look at, e.g., mass versus richness, by binning and averaging the
``fixed_y_c`` values in the output file and feeding the resulting catalog through the 
``nemoMass`` script.

Alternatively, forced photometry can be performed by the main ``nemo`` command. In this case, the path
to the catalog that provides the object prior positions is given by the ``forcedPhotometryCatalog`` key 
in the ``nemo`` configuration file. Forced photometry on sources of any kind can be done using this 
method.


Making model images
===================

The nemoModel script can make model images from a Nemo catalog. If run using MPI, it will break the map
into tiles to speed up processing, based on the given mask file. For example,::

    mpiexec nemoModel S18d_202003_test/S18d_202003_test_optimalCatalog.fits maps/Mar2020/AdvACTSurveyMask_v7_S18.fits Beams/190809/b20190809_s16_pa2_f150_nohwp_night_beam_profile_jitter_cmb.txt testModel_f150.fits -f 149.6 -M

will make a model image over the whole AdvACT footprint.


Using the selection function
----------------------------

Using the selection function...


Modifying GNFW parameters
-------------------------

How to do...


Source injection and flux recovery tests
----------------------------------------

Refer to examples...

