# Validation scripts

This directory contains some scripts that can be used to sanity check
the WebSky sims data and Nemo. These can be run on the output produced 
by running Nemo using the config files in the parent directory.

The `massRecovery*.py` scripts plot the masses recovered by Nemo, in
comparison with the true halo masses from the WebSky catalog.

The `makeMassFunctionPlotsCCL*.py` scripts can be used to check that
the cluster counts routines in Nemo (implemented using CCL) match
what is expected from the WebSky simulations.
