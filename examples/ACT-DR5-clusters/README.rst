Running the ACT DR5 Cluster Search
==================================

This is in the process of being updated for the DR5 release.

The `examples/AdvACT <https://github.com/simonsobs/nemo/tree/master/examples/AdvACT>`_
directory contains config files used for the AdvACT cluster 
search. There are two sets of files - one for night-time only
data taken up to 2016 (these have S16 in the file names) and
one for data taken up to 2018 (these have S18 in the file names).
In the examples that follow here, we will refer only to the S16 
files, but note that the set-up for both sets of maps is identical.

AdvACT has produced maps that cover most of the southern sky. To
search this large area efficiently, Nemo breaks up these maps into
tiles, which can be distributed to multiple processor cores using
MPI. The tiling procedure is controlled with the parameters listed
under 'tileDeck options' in each of the ``.yml`` config files. In the
current configuration, an automated procedure is used to divide the
maps up into 280 approximately 10 x 5 degree tiles.
The ``.sh`` files in this directory are batch job submission scripts
for the Slurm workload manager (as used on `Hippo <https://www.acru.ukzn.ac.za/~hippo/>`_\ ). 
The ``S16Pipeline.py`` and  ``S18Pipeline.py`` scripts show how to submit 
both the point source and cluster searches as dependent jobs.

For the AdvACT cluster search, we use the map itself as the noise
term in the matched filter constructed by Nemo (since the map itself
is not dominated by signal from clusters). However, powerful point
sources in the map can lead to ringing when the filter is applied.
So, these sources must be detected and masked before performing the
cluster search. In the current set-up, two passes of point source 
finding are performed. In the first pass (\ ``PS_S16_f150_auto_pass1.yml``\ ),
a catalog of 10-sigma sources is produced. These are then masked
(the holes are filled with a smoothed version of the map plus white
noise) and a second round of point source finding
(\ ``PS_S16_f150_auto_pass2.yml``\ ) is used to detect all sources down to
5-sigma. Both the point source searches are conducted on 150 GHz maps
only.

The cluster search (\ ``MFMF_S16_auto.yml``\ ) uses a multi-frequency matched
filter (e.g., Melin et al. 2006), modelling the cluster signal as a 
set of Universal Pressure Profiles with various different scales. 
Here, both sets of point source  catalogs are used to define a mask
before constructing the filter. Again, the holes in the mask are 
filled to ensure that the noise term in the matched filter is well 
behaved when Fourier transformed.
