Here we present minimal examples of how to use **Nemo** to detect clusters and sources
in a small section cut from the f090, f150 ACT DR5 maps. These should take only a minute
or so to run.

.. note::  The config files needed for this tutorial can be found in the 
           `examples/quickstart <https://github.com/simonsobs/nemo/tree/main/examples/quickstart>`_
           directory of the **Nemo** source distribution.
           
First, download the needed maps and beams (about 20 Mb), and extract the archive:

.. code-block::

   wget https://astro.ukzn.ac.za/~mjh/nemo-quickstart-maps.tar.gz
   tar -zxvf nemo-quickstart-maps.tar.gz

This will extract the maps and beam files under a ``maps/`` folder, and place a point
source catalog (used for masking in the cluster search), in the current working
directory.


Running a cluster search
========================

The file ``quickstart-clusters.yml`` contains a **Nemo** configuration that will conduct
a multifrequency matched filter search for clusters using a single cluster signal template.
Run it using:
    
.. code-block::

   nemo quickstart-clusters.yml

After about a minute, you should find an output directory named ``quickstart-clusters`` has
been created. You can find the output cluster candidate list in the file
``quickstart-clusters/quickstart-clusters_optimalCatalog.fits`` (FITS table format). You
will find FITS format filtered maps (both in terms of signal-to-noise and central
Comptonization parameter) under the ``quickstart-clusters/filteredMaps/PRIMARY/`` directory.


Running a point source search
=============================

The file ``quickstart-sources.yml`` contains a **Nemo** configuration that will conduct
a single frequency (98 GHz) matched filter search for 5σ point sources, using the
"two pass" source detection pipeline. This runs an initial search for extremely bright
sources in the first pass using a simple noise model, so that these sources can be subtracted
when using the map itself as the filter noise model in the second pass.
Run it using:
    
.. code-block::

   nemo quickstart-sources.yml

After about a minute, you should find an output directory named ``quickstart-sources`` has
been created. Similarly to the example cluster run above, you can find the output source
list in the file ``quickstart-sources/quickstart-sources_optimalCatalog.fits``, and the
filtered maps under the ``quickstart-sources/filteredMaps/PRIMARY/`` directory.


Next steps
==========

* See `Tutorials <https://nemo-sz.readthedocs.io/en/latest/tutorials.html>`_ for a number
  of more advanced step-by-step guides.

* See `Configuration File Parameters <https://nemo-sz.readthedocs.io/en/latest/config.html>`_
  for information about the settings that control **Nemo**.

* Refer to `Output <https://nemo-sz.readthedocs.io/en/latest/outputs.html>`_ for a
  guide to the data products produced by **Nemo**.
