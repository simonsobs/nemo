.. _Usage:

=============
Nemo Commands
=============

The **Nemo** package includes a number of command-line programs, each of which is described below.


.. _nemoCommand:
    
nemo
----

.. argparse::
   :filename: ../bin/nemo
   :func: makeParser
   :prog: nemo
   
   :program:`nemo` produces object catalogs and filtered maps using the parameter settings given
   in the YAML-format configuration file.


.. _nemoMassCommand:

nemoMass
--------

.. argparse::
   :filename: ../bin/nemoMass
   :func: makeParser
   :prog: nemoMass
   
   :program:`nemoMass` infers cluster masses based on their SZ signal as measured by
   :program:`nemo`. Cosmological and scaling relation parameters are set in the
   YAML-format configuration file.


.. _nemoMaskCommand:

nemoMask
--------

.. argparse::
   :filename: ../bin/nemoMask
   :func: makeParser
   :prog: nemoMask

   :program:`nemoMask` creates a mask image from a SAOImage DS9 region file.


.. _nemoMockCommand:

nemoMock
---------
   
.. argparse::
   :filename: ../bin/nemoMock
   :func: makeParser
   :prog: nemoMock
   
   :program:`nemoMock` generates random catalogs, using the selection function files
   produced by :program:`nemo`. Cosmological and scaling relation parameters are set
   in the YAML-format configuration file.


.. _nemoModelCommand:

nemoModel
---------

.. argparse::
   :filename: ../bin/nemoModel
   :func: makeParser
   :prog: nemoModel
   
   :program:`nemoModel` generates model images, i.e., signal-only maps containing
   either clusters or sources.


.. _nemoSpecCommand:
   
nemoSpec
--------

.. argparse::
   :filename: ../bin/nemoSpec
   :func: makeParser
   :prog: nemoSpec
   
   :program:`nemoSpec` extracts spectra (in ΔT with respect to the CMB)
   at locations given in the catalog, after first matching the angular resolution
   of the maps to the one with the worst resolution. Two modes are offered: 
   (i) compensated aperture photometry; and (ii) a matched filter. The latter
   was used in `Li et al. (2021) <https://ui.adsabs.harvard.edu/abs/2021arXiv210612467L/abstract>`_.

   
.. _nemoCatalogCheckCommand:

nemoCatalogCheck
----------------

.. argparse::
   :filename: ../bin/nemoCatalogCheck
   :func: makeParser
   :prog: nemoCatalogCheck
   
   :program:`nemoCatalogCheck` cross matches an external catalog against the
   output produced by :program:`nemo`, and reports which objects are detected in
   the :program:`nemo` catalog, which are missing, and which are outside the
   survey footprint.
