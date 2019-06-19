Nemo is a map filtering and source detection and characterization pipeline, designed to find
galaxy clusters using the Sunyaev-Zel'dovich effect. It can also be used to find sources. 

* Documentation: https://acru.ukzn.ac.za/~mjh/nemo/
* License: Will be BSD 2-clause when released (this is a currently a private repository: for 
now nemo cannot be used for non-ACT, non-SO approved projects; however, the aim is to make the 
code fully public by the time the first AdvACT cluster catalog is published)

Nemo is *not* the pipeline used for [Hasselfield et al. (2013)](http://adsabs.harvard.edu/abs/2013JCAP...07..008H), 
but implements many of the ideas presented there, and should give similar results, given 
the same map (or at least it did in the past). *It is* the pipeline, that has been used for the
[two-season ACTPol cluster catalog paper](http://adsabs.harvard.edu/abs/2017arXiv170905600H).

A slightly modified version of Matthew Hasselfield's `gnfw.py` code (used for modeling cluster
pressure profiles) is included in nemo.

See [examples/equD56/](examples/equD56/) if you would like to re-create the 
cluster catalog produced for the E-D56 field; [examples/AdvACT/](examples/AdvACT/) if you want to
see example config files currently being used for the AdvACT cluster search; and 
[examples/pointSources/](examples/pointSources) for config files that can be used for detecting
point sources.

# Current status and previous internal releases

**Note**: Versions in the git respository after the last tagged release (v18.06) use a new config 
file format (.yml) - see the [examples/](examples/) directory. We only guarantee compatibility
with Python 3.x (and we're no longer testing on Python 2.7).

Listed below are tagged, internal releases of nemo:

* **v18.06:** Brought up-to-date with changes that were made in tileDeck branch (can run in 
  parallel using MPI); this is the version used for making S16 point source catalogs (by Sigurd) 
  and preliminary AdvACT cluster catalogs.
* **v17.10:** Version of the code at the time of submission of the ACTPol cluster catalog paper.

# Software needed

Nemo itself is written in Python (3.x), and requires the following additional modules to be installed 
(currently used versions are given in brackets, earlier and later versions also probably work):

* numpy (1.13.3)
* scipy (0.19.3)
* matplotlib (2.1.1)
* astLib (0.10.x + or git version: get it with `git clone http://git.code.sf.net/p/astlib/git astlib-git`)
* [pixell](https://github.com/simonsobs/pixell/) (0.5.3 or git version)
* Python Imaging Library (PIL or PILLOW - e.g., `sudo easy_install Pillow` on Mac)
* astropy (3.0.5)
* IPython (7.2.0)
* Cython (0.24.1)
* PyYAML (3.12)
* hmf (3.0.6)
* mpi4py (3.0.0)
* colorcet (1.0.0; https://github.com/bokeh/colorcet/releases)

All of the dependencies can be installed using `pip`.

# Installation

As root:
    
```
sudo python setup.py install
```

Or, in your home directory:
    
```
python setup.py install --prefix=$HOME/local
```

then add `$HOME/local/bin` to $PATH, and e.g., `$HOME/local/lib/python2.7/site-packages` to 
$PYTHONPATH (adjust the path according to your Python version number).

```
export PATH=$HOME/local/bin:$PATH    
export PYTHONPATH=$HOME/local/lib/python2.7/site-packages:$PYTHONPATH
```

Alternatively, 

```
python setup.py install --user
```

will install `nemo` under `$HOME/.local` (on Ubuntu), and in some other default location on Mac.

# Running nemo

See [examples/equD56/README.md](examples/equD56/README.md) for a tutorial on how to re-create 
the ACTPol two-season cluster catalog (including mass estimates). 
See [examples/AdvACT/](examples/AdvACT/) for example .yml config files for the current AdvACT
cluster search. Refer to the comments in the .yml config files themselves for information on what
each parameter does.

# Related codes

Catalogs produced by nemo can be fed into:

* [Sourcery](https://github.com/mattyowl/sourcery) - a package that creates and serves a live,
  editable web database with multi-wavelength imaging.
* [zCluster](https://github.com/ACTCollaboration/zCluster) - a code for estimating cluster
  photometric redshifts.

# Comments, bug reports, help, suggestions etc.

Please contact Matt Hilton <matt.hilton@mykolab.com>.
