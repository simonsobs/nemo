.. note::  This is under construction.

Here we provide a guide to reproducing the ACT DR6 cluster search
data products, as described in the `ACT DR6 cluster catalog paper <https://ui.adsabs.harvard.edu/abs/2025arXiv250721459H/abstract>`_.

.. note::  The config files and scripts needed for this tutorial can be
           found in the `examples/ACT-DR6-clusters <https://github.com/simonsobs/nemo/tree/main/examples/ACT-DR6-clusters>`_
           directory of the **Nemo** source distribution.


Downloading maps and other needed inputs
========================================

The first step is to download the co-added ACT+Planck maps used for the ACT DR6 cluster search.
These are described in `Naess et al. (2025) <https://ui.adsabs.harvard.edu/abs/2025arXiv250314451N/abstract>`_,
and available on `LAMBDA <https://lambda.gsfc.nasa.gov/product/act/act_dr6.02/>`_. In
addition to this, we require the beam transform files and some masks. These can be downloaded and extracted by
running a script:

.. code-block::

   sh FETCH_MAPS.sh

Each map is 5 Gb in size, so this may take a while.
