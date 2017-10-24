# Re-creating the ACTPol two-season cluster catalog

Here is an example of how to re-create the cluster catalog presented 
in the [two-season ACTPol cluster catalog paper](http://adsabs.harvard.edu/abs/2017arXiv170905600H). 
This uses nemo's RealSpaceMatchedFilter method.


## Step 1:

This directory already contains some of the things needed - e.g, the
survey and point source masks (`surveyMask.fits.gz` and 
`pointSourceMask.fits.gz` respectively), the beam profile used (under 
`profiles_ACT/`; from ACT actually, but ACTPol is not significantly 
different for cluster-finding purposes), and the .par file, from which
nemo reads its settings (see below). However, the maps are too big to
include in the github repository, so we need to get them from 
elsewhere.

For the ACTPol two-season cluster catalog, we combined the ACT+ACTPol
maps (which are available on LAMBDA). You can download the combined 
map and weight files (403 Mb) using:

```
wget http://www.acru.ukzn.ac.za/~mjh/ACT/equD56Maps.tar.gz --user=act --password=atacamallama
```

Extract this in the current directory (i.e., the same place as the 
equD56.par file). 

## Step 2:

Run `nemo` using:

```
nemo equD56.par
```

You can check the settings by opening the .par file in a text editor.
This may take ~20-30 minutes to run, because a lot of different filter
scales are used (for a "cosmological sample", you would only need to
run the `'Arnaud_M2e14_z0p4'` filter, and could comment out with # the
other dictionaries defining each other filter in the mapFilters 
list).

Output is written to the equD56 directory. Here you will find 
catalogs (.fits tables, e.g., `equD56_optimalCatalog.fits`), DS9 region
(.reg files), images (both in terms of y and signal-to-noise; under
`filteredMaps/`), and a bunch of other stuff under `diagnostics/` (e.g.,
statistics (.fits tables) and plots of the contamination rate as a 
function of SNR, estimated by running the detection algorithm over 
inverted maps - you should recognise several of the plots here from
the paper draft linked to above).

## Step 3:

If you want to measure masses, a .fits table that includes the columns
'redshift' and 'redshiftErr' is needed. The `ACTPol_redshifts.fits` file
in this directory contains all redshifts that have been assigned to 
clusters on the web database. You can specify the redshift catalog 
used for mass estimates with the `massOptions` key in the .par file. 
Since this is already filled in, you can run the mass estimation 
script with:

```
nemoMass equD56.par
```

The output is written as a .fits table, `equD56/equD56_M500.fits`.

If you wanted to extract a catalog containing redshifts from the web 
database, this link shows the query needed:

<http://www.acru.ukzn.ac.za/actpol-sourcery/updateQueryParams?queryRADeg=0%3A360&queryDecDeg=-90%3A90&querySearchBoxArcmin=&queryOtherConstraints=sourceList+%3D+equD56-MJH+and+redshift+%3E+0&queryApply=Apply>

(you may need to enter 'act', 'atacamallama' to see the web database).
You could download the catalog (in .fits format) using the links at
the bottom of the web database table page.

Similarly, if you wanted everything flagged as a cluster (i.e., the
sample of 182 clusters in the two-season ACTPol catalog), you can use
this query:

<http://www.acru.ukzn.ac.za/actpol-sourcery/updateQueryParams?queryRADeg=0%3A360&queryDecDeg=-90%3A90&querySearchBoxArcmin=&queryOtherConstraints=sourceList+%3D+equD56-MJH+and+classification+%3D+%27cluster%27&queryApply=Apply>

This should match the `ACTPol_clusters.fits` file included in the 
current directory.

If you're interested in measuring photo-zs for clusters detected with
`nemo`, check out `zCluster`: 

<https://github.com/ACTCollaboration/zCluster>

This can take the `equD56/equD56_optimalCatalog.fits` file as input.

## Step 4:

Steps 1-3 are all that are needed to recreate the two-season ACTPol
cluster catalog. If you wanted to run simulations to estimate 
completeness, you can check out the nemoSelFn script. This runs in the
same way as the others (taking the .par file for input), and writes
output to `equD56/diagnostics`. With the settings in the equD56.par 
given in this directory, it only makes sense to run this on a cluster
using MPI. For applications of this, see the `nemo/MockSurvey.py` and 
`nemo/SelFn.py` modules, as well as the `bin/nemoSelFn` script itself.

