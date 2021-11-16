Nemo uses the `Robot Framework <http://robotframework.org/>`_ for tests.

.. note::  The files needed to run the test suite are found in the
           `tests <https://github.com/simonsobs/nemo/tree/main/tests>`_
           directory of the **Nemo** source distribution.

To run a set of tests that should take about 15 minutes:

.. code-block::

   robot quick.robot

Check the ``plots/`` directory for output from some tests.

The other test suites are in the process of being retired and/or
reformulated. They can take a long time to run on real data.

To run a single test suite do, e.g.,

.. code-block::

   robot point_sources.robot

To run a single test in a test suite do, e.g., 

.. code-block::

   robot -t "Recover published 2008 survey source fluxes - Fourier space filter" point_sources.robot
   robot -t "Recover published ACTPol cluster masses - Fourier space filter" clusters.robot

To run all of the tests:

.. code-block::

   robot -N "Nemo Tests" *.robot
