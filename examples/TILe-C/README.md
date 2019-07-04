The [examples/TILe-C](https://github.com/simonsobs/nemo/tree/master/examples/TILe-C)
directory contains a config file (`y_f090beam.yml`) that can be
used to run Nemo on y-maps produced by
[TILe-C](https://github.com/ACTCollaboration/tile-c).
You will need access to the `release_v0.2.3` D56 y-maps maps and the
90 GHz beam (the filename is given in the config file). You will also
need a survey mask that uses only the valid region of the map 
(to mitigate the effect of the apodization that was applied by
TILe-C). You can download that (204 kb) using:

```
wget https://acru.ukzn.ac.za/~mjh/tilec_deep56_surveyMask.fits.gz
```

After adjusting the paths in the config file to point to where you 
have kept the maps and beams, you can run Nemo using:

```
nemo y_f090beam.yml
```

This will produce filtered maps at 5 different scales (edit the config
file if you want to add more), and will take ~17 min to run. Other 
nemo scripts will work on the output produced here (e.g., `nemoMass`,
`nemoSelFn`, `nemoCatalogCheck` etc.).
