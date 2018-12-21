Nemo uses the [Robot Framework](http://robotframework.org/) for tests. At the 
moment, these are end-to-end tests that take a while to run and use real data.

To run a single test do, e.g.,

```
robot point_sources.robot
```

To run all of the tests, do:

```
robot -N "Nemo Tests" *.robot
```

or see the `RUN_ALL.sh` script.

Check the `plots/` directory for output from some tests.
