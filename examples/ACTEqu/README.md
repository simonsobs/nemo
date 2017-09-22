# Nemo Quick Start

This directory contains .par files for running over the ACT equatorial maps that are available on 
LAMBDA. You'll need the files named:

* `ACT_148_equ_season_3_1way_hits_v3.fits`  
* `ACT_148_equ_season_3_1way_v3_src_free.fits`

You can fetch them with:

`wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/data2013/Maps/AR1/Equa/ACT_148_equ_season_3_1way_v3_src_free.fits`
`wget https://lambda.gsfc.nasa.gov/data/suborbital/ACT/data2013/Weights/AR1/Equa/ACT_148_equ_season_3_1way_hits_v3.fits`

If you don't save these maps in this directory, you'll need to edit the .par files accordingly to 
point to wherever the maps live on your machine. In practice, one would make a weighted average of
all of the season / array maps and run on that to use the full available depth.

There are two examples .par files - one for clusters (`equClusters.par`), and another for point 
sources (`equPointSources.par`).

**Note that flux or y measurements extracted with the filters used here will not be normalized
correctly. If you require photometry, check out `examples/AdvACT/` and `examples/equD56/`, which
use the RealSpaceMatchedFilter.**

## Clusters (equClusters.par)

This example runs only a single cluster model, but you can add as many more as you like - they 
will be combined to make an "optimal" catalog, selecting the higher signal-to-noise result as 
the preferred one. It just takes longer.

For speed, this example runs over only a small section of the map, which you can adjust with the
`RADecSection` key in the .par file.

To run:

```
nemo equClusters.par
```

For the settings in this example, this should take ~20 sec on an Intel Core i5-3320M (2.60 GHz). 
Subsequent runs (if you fiddle with the code and decide to re-run, for example) will be quicker, 
if you don't modify the filters, as they are cached.

Output will be written to a directory called `equClusters`, in which you will find:

* `equClusters/diagnostics`           - various diagnostic plots
* `equClusters/filteredMaps`          - the filtered maps themselves, with .reg files for each model
* `equClusters/equ_optimal_catalog.*` - the combined, optimal catalog, in various formats 

The filtered map is in units of yc.

The file `Hasselfield2013_upp.fits` contains the Hasselfield et al. (2013) catalog, if you want to 
compare by cross matching using TopCat. Please note the different area and maps used in this quick
example when doing so.

Note that the output for photometric measurements (Y500 etc.) is not correct at the moment - this 
was being experimented with a while ago and has been left 'as is'.

## Point Sources (equPointSources.par)

This is the equivalent set-up to run a matched filter for point sources, using the beam profile
specified under `unfilteredMaps` in `equPointSources.par`.

To run:

```
nemo equPointSources.par
```

Output will be written in the `equPointSources` directory. Note that in this case, the filtered map is
in units of uK, and the same caveats about photometric measurements apply. One can change 
`outputUnits` in the .par file to Jy/beam, but the code to actually do this has not been added yet.
