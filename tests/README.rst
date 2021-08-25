Nemo uses the `Robot Framework <http://robotframework.org/>`_ for tests.

To run a set of tests that should take only a few minutes, use,

.. code-block::

   robot quick.robot

The other test suites in are in the process of being replaced and/or
reformulated. These can take a long time to run on real data.

To run a single test suite do, e.g.,

.. code-block::

   robot point_sources.robot

To run a single test in a test suite do, e.g., 

.. code-block::

   robot -t "Recover published 2008 survey source fluxes - Fourier space filter" point_sources.robot

To run all of the tests:

.. code-block::

   robot -N "Nemo Tests" *.robot

or see the ``RUN_ALL.sh`` script.

Check the ``plots/`` directory for output from some tests.
