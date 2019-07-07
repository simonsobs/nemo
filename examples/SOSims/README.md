The [examples/SOSims](https://github.com/simonsobs/nemo/tree/master/examples/SOSims) 
directory contains scripts for working with
the simulated Simons Observatory maps made by the 
[Map-Based Simulation Pipeline Working Group](https://github.com/simonsobs/map_based_simulations).
More tools for handling the simulated maps may be available elsewhere,
but here you should be able to find enough to get you going with 
running Nemo on them.

# Making multi-component simulated maps

The map-based simulations are provided in HEALPix format, with one map 
for each component (CMB, CIB, tSZ, noise etc.) at each frequency. The 
`combineComponentMaps.py` script can be used to add these
together. To use this, you will first have to download the needed 
[simulation maps](https://github.com/simonsobs/map_based_simulations) and 
edit `healpixSimsDir` in the script to point to the directory where 
these can be found. You can then simply run this script with:

```
python combineComponentMaps.py
```

As provided, this script will make combined maps at four frequencies
(93, 145, 225, 280 GHz) featuring CMB, CIB, tSZ, and noise.

# Making CAR projection maps

At the moment, Nemo runs on FITS images that have a valid World 
Coordinate System (WCS) in the header. For Simons Observatory, these
use the plate carrée (CAR) projection. However, the map-based 
simulations are currently provided in HEALPix format. The 
`healpix2CAR.py` script (a modified version of a script by
Mathew Madhavacheril taken from the [TILe-C](https://github.com/ACTCollaboration/tile-c)
package) can be used to re-project these to CAR. Look at the 
`HEALPIX_TO_CAR.sh` script to see an example of its use.

As one of its inputs, the `healpix2CAR.py` takes a plain-text file 
containing the desired output image header. This defines both the sky 
area and the pixel scale. There are a few examples provided:

* `template_small.txt`: Corresponds to a subsection of the ACT D56 field. 
* `template_E-D56.txt`: Contains the full ACTPol E-D56 field footprint.
* `template_AdvACT.txt`: The pixelization used for AdvACT maps. To use Nemo on very large maps (as in `template_AdvACT.txt`), you will need to set-up the configuration file to break the map into tiles (see below).

If you look at the `HEALPIX_TO_CAR.sh`, you will see that it uses an option
to add a small amount of additional white noise to the output CAR maps. 
This is currently being used to avoid the effect of ringing in the CAR maps 
(which arises because the HEALPix maps are lower resolution than the output
CAR maps - when NSIDE = 8192 maps are available, this may not be a problem).
With the additional white noise, the matched-filtered maps made by Nemo from
these simulations look okay when viewed with, e.g., 
[DS9](http://ds9.si.edu/site/Home.html).

# Making beam files

Nemo understands beams given in a simple plain-text format, as used by
ACT. For the map-based simulations, Gaussian beams are used with FHWM
set according to the values given in Table 1 of the 
[Simons Observatory: Science Goals and Forecasts](https://ui.adsabs.harvard.edu/abs/2019JCAP...02..056A/abstract) 
paper. To generate the necessary beam files, run:
```
MAKE_BEAMS.sh
```
which in turn calls the `makeGaussianBeam.py` script.

# Running Nemo on small maps

After running the above scripts, you should now have a directory that contains
CAR-projected maps (e.g., `TOnly_la093_CAR.fits`) and beam files 
(e.g., `beam_gaussian_la093.txt`). Various configuration files are provided 
that you can use to produce matched-filtered maps and catalogs using Nemo - the
one that seems to work the best currently uses only three frequencies 
(93, 145, 225 GHz). You can run this with:

```
nemo MFMF_SOSim_3freq_small.yml
```

After this finishes, you will find matched-filtered maps under 
`MFMF_SOSim_3freq_small/filteredMaps/`, while 
`MFMF_SOSim_3freq_small/MFMF_SOSim_3freq_small_optimalCatalog.fits` contains the 
candidate cluster catalog. Note that this config file filters the maps at only 
one scale - you can filter at multiple scales by uncommenting the entries under
`mapFilters` in `MFMF_SOSim_3freq_small.yml`.

Configuration files are also supplied to filter using either two or four 
frequencies (at the time of writing, more ringing is seen when including the 
280 GHz map).

# Extracting the survey mass limit

The `nemoSelFn` script can be used to obtain estimates of the 90% mass 
completeness threshold at some signal-to-noise level (typically S/N > 5). This 
is subject to the assumed scaling relation parameters defined in the Nemo 
configuration file. You can run this with:
    
```
nemoSelFn MFMF_SOSim_3freq_small.yml
```

This script produces various plots related to the survey completeness as a 
function of redshift under `MFMF_SOSim_3freq_small/diagnostics/`, including a
mass-limit map.

The limits within the footprints of other surveys that intersect with the 
filtered maps can be obtained by supplying survey mask files in the 
`selFnFootprints` parameter in the configuration file (these are commented out
in `MFMF_SOSim_3freq_small.yml`). 

# Running Nemo on large maps

To do this, you will first need to output large area CAR maps using (for example)
the `healpix2CAR.py` script with `template_AdvACT.txt`. You will then need to 
use a config file that tells Nemo how to break the map up into tiles, such as
the supplied `MFMF_SOSims_3freq_tiles.yml` (the tile definitions can be found at
the bottom of this file). 

Nemo runs in parallel over MPI, distributing tiles to different processor cores. 
An automatic tiling algorithm will eventually be included in Nemo, but for now 
you can find a script that generates the tile definitions that have already been
included in `MFMF_SOSims_3freq_tiles.yml` in `autotiler.py`. This script makes 
use of a "survey mask", which defines the region that Nemo will search for 
clusters (or sources). 

The script `createSurveyMask.py` generates the survey mask FITS image, taking as
input a DS9 region file containing polygon-shaped regions (see 
`surveyMask.reg`), and a plain-text file that defines the image header
(`template_AdvACT.txt`). You can run this with:

```
python createSurveyMask.py template_AdvACT.txt surveyMask.reg
```

This writes the output survey mask to `surveyMask.fits.gz`. You will need this 
file to run Nemo using the included `MFMF_SOSims_3freq_tiles.yml` configuration.

You can run Nemo in parallel using, e.g.,
```
mpiexec --np $NUM_PROCESSES nemo MFMF_SOSim_3freq_tiles.yml -M
```
replacing `$NUM_PROCESSES` with the number of cores you want to run on. Nemo
will divide up the tiles between processors as evenly as it can. The file 
`slurm_nemo.sh` shows how Nemo can be run on a cluster that uses the 
[Slurm](https://slurm.schedmd.com/overview.html) job scheduler.

# Using the input simulation catalogs

Add stuff on how to use the WebSky catalogs here, running `nemoMass` with the given
redshifts there, etc.. - to do.
