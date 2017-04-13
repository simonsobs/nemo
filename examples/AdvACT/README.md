An example of using BeamRealSpaceMatchedFilter to find point sources.

Using Sigurd's maps: 

`~sigurdkn/data/actpol/maps/s4test3/cmb_huge`

Look at the .par files to see which maps are used there, and copy/link them into this directory. 
Note that since Sigurd's maps include polarisation, the script makeSigurd1D.py should be used to 
extract the temperature map only.

For testing out the code, run: 

```
nemo AdvACT_PointSources_quick.par
```

This is very quick (~30 sec), as it uses only a small part of the map. Under the 
`AdvACT_PointSources_quick` directory that is created, you'll find the catalogs and filtered maps in the usual format.

To run over the whole map, use:

```
nemo AdvACT_PointSources.par
```

This will start thrashing the disk on a machine with 8 Gb of RAM. Currently, it will still take a 
loooong time on a machine with lots of memory (maybe ~1 hour to do the filtering step) - the code 
needs to be parallelised to cope with the size of the AdvACT maps.