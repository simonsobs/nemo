Nemo uses the [Robot Framework](http://robotframework.org/) for tests. At the 
moment, these are end-to-end tests that take a while to run and use real data.

To run a single test suite do, e.g.,

```
robot point_sources.robot
```

To run a single test in a test suite do, e.g., 

```
robot -t "Recover published 2008 survey source fluxes - Fourier space filter" point_sources.robot
```

To run all of the tests:

```
robot -N "Nemo Tests" *.robot
```

or see the `RUN_ALL.sh` script.

Check the `plots/` directory for output from some tests.
