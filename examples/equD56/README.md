# Re-creating the ACTPol two-season cluster catalog

Here is an example of how to re-create the cluster catalog presented 
in the [two-season ACTPol cluster catalog paper](http://adsabs.harvard.edu/abs/2017arXiv170905600H). 
This uses nemo's RealSpaceMatchedFilter method.

If you want to see an example .yml config file that breaks the map 
up into tiles and runs in parallel, see the [examples/AdvACT](../AdvACT/) 
directory instead.

## Step 1:

This directory already contains some of the things needed - e.g, the
survey and point source masks (`surveyMask.fits.gz` and 
`pointSourceMask.fits.gz` respectively), the beam profile used (under 
`profiles_ACT/`; from ACT actually, but ACTPol is not significantly 
different for cluster-finding purposes), and the .yml config file, from
which nemo reads its settings (see below). However, the maps are too 
big to include in the github repository, so we need to get them from 
elsewhere.

For the ACTPol two-season cluster catalog, we combined the ACT+ACTPol
maps (which are available on LAMBDA). You can download the combined 
map and weight files (403 Mb) using:

```
wget http://www.acru.ukzn.ac.za/~mjh/ACT/equD56Maps.tar.gz --user=act --password=atacamallama
```

Extract this in the current directory (i.e., the same place as the 
equD56.yml file). 

## Step 2:

Run `nemo` using:

```
nemo equD56.yml
```

You can check the settings by opening the .yml file in a text editor.
This may take ~20-30 minutes to run, because a lot of different filter
scales are used (for a "cosmological sample", you would only need to
run the `'Arnaud_M2e14_z0p4'` filter, and could comment out with # the
other dictionaries defining each other filter in the mapFilters 
list - this is already done in the `equD56_quick.yml` file that you can
find in this directory).

Output is written to the `equD56` directory. Here you will find 
catalogs (.fits tables, e.g., `equD56_optimalCatalog.fits`), DS9 region
(.reg files), images (both in terms of y0 and signal-to-noise; under
`filteredMaps/`), and a bunch of other stuff under `diagnostics/`.

## Step 3:

If you want to measure masses, a .fits table that includes the columns
'redshift' and 'redshiftErr' is needed. The `ACTPol_redshifts.fits` file
in this directory contains all redshifts that were assigned to 
cluster candidates on the web database. You can specify the redshift catalog 
used for mass estimates with the `massOptions` key in the .yml file. 
Since this is already filled in, you can run the mass estimation 
script with:

```
nemoMass equD56.yml
```

This will take ~10 minutes to run (initially), as it calculates the 
filter mismatch function Q (the result of this is cached, so subsequent
runs would be faster). The output is written as a .fits table, 
`equD56/equD56_M500.fits`.

If you're interested in measuring photo-zs for clusters detected with
`nemo`, check out `zCluster`: 

<https://github.com/ACTCollaboration/zCluster>

This can take the `equD56/equD56_optimalCatalog.fits` file as input.

## Step 4:

Steps 1-3 are all that are needed to recreate the two-season ACTPol
cluster catalog. You can check the output by cross-matching against
the `ACTPol_clusters.fits` catalog using e.g. TopCat.

If you wanted to run simulations to estimate completeness, you can use
the `nemoSelFn` script. You can run it with:

```
nemoSelFn equD56.yml
```

The output for this script is written in the `diagnostics/` directory,
and includes a plot of the 90% completeness limit, averaged over the
survey (in this case the E-D56 field), and an equivalent mass limit 
map, evaluated at z = 0.5 (a .fits image). It is making the latter
that takes up most of the time - this can be disabled by removing
the `massLimitMaps` key from the `selFnOptions` dictionary in the .yml
file. No doubt this part could be sped up, but at the time of writing,
this script takes ~40 minutes to run. As with `nemoMass`, some results
are cached, so repeat runs are much quicker.

Note that `nemoSelFn` has been completely rewritten and is different
to the version used for the ACTPol paper. Hence, the results are 
slightly different.
