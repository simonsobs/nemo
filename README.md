# nemo

Map filtering and SZ cluster detection and characterization pipeline. *Not* the pipeline used
for [Hasselfield et al. (2013)](http://adsabs.harvard.edu/abs/2013JCAP...07..008H), but should 
give similar results for signal-to-noise, given the same map (or at least it did in the past).

This code is posted "as is", in time for the Feb 2017 meeting, although it is in need of tidy up. 
There may be all sorts of weird things in here that don't make sense, wrong turns, and all the rest.
In particular, don't expect the source characteristics (Y500 etc.) besides signal-to-noise to 
correspond to actual reality. Note that the catalog matching and (static) web page generation
is disabled, as that's all now handled by the [sourcery](https://github.com/mattyowl/sourcery) 
package, which runs a live, editable web database.

Matthew Hasselfield's `gnfw.py` code (used for Arnaud profiles) is included in nemo, as is a modified 
version of `actDict.py` - thanks to whoever wrote that back in the day.

## Software needed

Nemo itself is written in python (2.7.x), and requires the following additional modules to be installed 
(current versions used by the author are given in brackets, earlier and later versions also probably work):

* pyfits (3.3)
* numpy (1.11.1)
* scipy (0.17.1)
* matplotlib (1.5.2)
* astLib (git version probably: get it with `git clone http://git.code.sf.net/p/astlib/git astlib-git`)
* flipper ([ACT collaboration git version](https://github.com/ACTCollaboration/flipper))
* astropy (1.3)
* IPython (2.4.1)
* Cython (0.24.1)

_Note:_ Switched from `atpy` to `astropy` for handling .fits tables (Feb 2017).

## Installation

As root:
    
```
sudo python setup.py install
```

Or, in your home directory:
    
```
python setup.py install --prefix=$HOME/local
```

Then add `$HOME/local/bin` to $PATH, and e.g., `$HOME/local/lib/python2.7/site-packages` to $PYTHONPATH.

```
export PATH=$HOME/local/bin:$PATH    
export PYTHONPATH=$HOME/local/lib/python2.7/site-packages:$PYTHONPATH
```

## Running nemo

See the README file under `examples/ACTEqu/` for a quick tutorial.

## Comments, bug reports, help, suggestions etc..

Please contact Matt Hilton <matt.hilton@mykolab.com>.

