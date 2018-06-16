# Running nemo on AdvACT maps

In this directory you will find example .yml config files that show how to 
run nemo in "tileDeck" mode - i.e., splitting a very large map into
tiles and then (optionally) running nemo in parallel using MPI. Given
that this is work in progress, no example maps or beam files are 
provided - check the ACTPol wiki for where to find such things.

For a quick tutorial on how to run `nemo`, see 
[examples/equD56/README.md](../equD56/README.md) first. The key difference 
between the AdvACT .yml files in this directory and the equD56 example are 
the "tileDeck" options. These are documented in the .yml files themselves, 
but the key parameters are:

* `makeTileDeck` - set to True for any of the other tileDeck options to be
  used (otherwise nemo will not split the map into tiles).
* `tileDefinitions` - a dictionary list that defines the names of the tiles, and
  the corresponding RA, dec range covered.
* `tileNoiseRegions` - a dictionary that defines the regions used for 
  estimating the noise in RealSpaceMatchedFilter. Hence, `RADecSection` in 
  `noiseParams` in the filter definitions is set to `tileNoiseRegions` for
  these to be used.
* `extNameList` - this can be used to select a particular set of tiles only,
  most useful for testing purposes. To see this in action, see
  `MF_AdvACT_multiScale_tileDeck_hybrid_quick.yml`.

Additionally, if `useMPI: True` is found in the .yml file, the `nemo` 
scripts should be launched with `mpirun`, e.g.,

```
mpirun --np 16 nemo MF_AdvACT_multiScale_tileDeck_hybrid.yml
```

The `*.sh` files in this directory are the batch files used for running 
the code on [hippo](https://www.acru.ukzn.ac.za/~hippo/) using `slurm`.
