Nemo is a map filtering and source detection and characterization pipeline, designed to find
galaxy clusters using the Sunyaev-Zel'dovich effect. It can also be used to find sources. 

* Documentation: <https://astro.ukzn.ac.za/~mjh/nemo/>
* License: Will be BSD 2-clause when released (this is a currently a private repository: for 
now Nemo cannot be used for non-ACT, non-SO approved projects; however, the aim is to make the 
code fully public by the time the first AdvACT cluster catalog is published).
* Authors: Matt Hilton, with contributions from Simone Aiola, David Alonso, Matthew Hasselfield, 
Toby Marriage, Sigurd Naess, and Cristóbal Sifón.

Nemo is *not* the pipeline used for [Hasselfield et al. (2013)](http://adsabs.harvard.edu/abs/2013JCAP...07..008H), 
but implements many of the ideas presented there, and should give similar results, given 
the same map (or at least it did in the past). *It is* the pipeline, that has been used for the
[two-season ACTPol cluster catalog paper](http://adsabs.harvard.edu/abs/2017arXiv170905600H),
and the [ACT DR5](https://ui.adsabs.harvard.edu/abs/2020arXiv200911043H/abstract) cluster catalog paper.

A slightly modified version of Matthew Hasselfield's `gnfw.py` code (used for modeling cluster
pressure profiles) is included in Nemo.

See [examples/equD56/](examples/equD56/) if you would like to re-create the 
cluster catalog produced for the E-D56 field; [examples/AdvACT/](examples/AdvACT/) if you want to
see example config files currently being used for the AdvACT cluster search; and 
[examples/pointSources/](examples/pointSources) for config files that can be used for detecting
point sources.

# Current status and previous internal releases

We now tag releases using the format vMAJOR.MINOR.PATCH. The first tagged release using this 
scheme is v0.1.0 (January 2020). The previously tagged internal releases of Nemo have been 
removed, except for the version used for the 
[two-season ACTPol cluster catalog paper](http://adsabs.harvard.edu/abs/2017arXiv170905600H) - this is now tagged as ACTPol2018.

# Software needed

Nemo itself is written in Python (3.6+), and requires the following additional modules to be installed 
(currently used versions are given in brackets, later versions also probably work):

* numpy (1.13.3)
* scipy (1.3.0)
* matplotlib (2.1.1)
* astLib (0.11.3)
* [pixell](https://github.com/simonsobs/pixell/) (0.6.3 or git version)
* Pillow (5.1.0)
* astropy (3.2.1)
* Cython (0.24.1)
* PyYAML (3.12)
* Colossus (1.2.9)
* [CCL](https://github.com/LSSTDESC/CCL) (2.1 or later)
* mpi4py (3.0.0)
* colorcet (1.0.0; https://github.com/bokeh/colorcet/releases)

All of the dependencies can be installed using `pip`, and should be installed automatically as needed
by the `setup.py` script (your mileage may vary in terms of how successful `pip` is at building
some of the external dependencies, depending on your set up).

# Installation

As root:
    
```
sudo python setup.py install
```

Alternatively, 

```
python setup.py install --user
```

will install `nemo` under `$HOME/.local` (on Ubuntu), and in some other default location on Mac.

You can also use the `--prefix` option, e.g.,
    
```
python setup.py install --prefix=$HOME/local
```

then add `$HOME/local/bin` to $PATH, and e.g., `$HOME/local/lib/python3.6/site-packages` to 
$PYTHONPATH (adjust the path according to your Python version number).

```
export PATH=$HOME/local/bin:$PATH    
export PYTHONPATH=$HOME/local/lib/python3.6/site-packages:$PYTHONPATH
```


# Running Nemo

Documentation is available at <https://astro.ukzn.ac.za/~mjh/nemo/>, including a number of
tutorials.

See [examples/equD56/README.md](examples/equD56/README.md) for a tutorial on how to re-create 
the ACTPol two-season cluster catalog (including mass estimates). 

See [examples/AdvACT/](examples/AdvACT/) for example .yml config files for the current AdvACT
cluster search. Refer to the comments in the .yml config files themselves for information on what
each parameter does.

# Related codes

Catalogs produced by Nemo can be fed into:

* [Sourcery](https://github.com/mattyowl/sourcery) - a package that creates and serves a live,
  editable web database with multi-wavelength imaging.
* [zCluster](https://github.com/ACTCollaboration/zCluster) - a code for estimating cluster
  photometric redshifts.

# Comments, bug reports, help, suggestions etc.

Please contact Matt Hilton <matt.hilton@mykolab.com>.
