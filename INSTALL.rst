**Nemo** is written in `Python <https://www.python.org/>`_ (3.6+), and requires the
following additional modules to be installed (currently used versions are given in
brackets, later versions also probably work):

* numpy (1.19)
* scipy (1.3.0)
* matplotlib (2.1.1)
* astLib (0.11.7)
* `pixell <https://github.com/simonsobs/pixell/>`_ (0.17 or later)
* Pillow (5.1.0)
* astropy (4.0)
* PyYAML (3.12)
* `CCL <https://github.com/LSSTDESC/CCL>`_ (3.0.0+)
* mpi4py (3.0.0)
* colorcet (1.0.0; https://github.com/bokeh/colorcet/releases)

The latest tagged version of **Nemo** can be installed using ``pip``:
    
.. code-block::

   pip install nemo-sz

Dependencies will be installed by ``pip``, except for ``pyccl`` and ``mpi4py``.

You may also install using the standard ``setup.py`` script, e.g., as root:

.. code-block::

   sudo python setup.py install

Alternatively, 

.. code-block::

   python setup.py install --user

will install ``nemo`` under ``$HOME/.local`` (on Ubuntu), and in some other default location on Mac.

You can also use the ``--prefix`` option, e.g.,

.. code-block::

   python setup.py install --prefix=$HOME/local

and then add ``$HOME/local/bin`` to $PATH, and e.g., ``$HOME/local/lib/python3.6/site-packages`` to 
$PYTHONPATH (adjust the path according to your Python version number).

.. code-block::

   export PATH=$HOME/local/bin:$PATH    
   export PYTHONPATH=$HOME/local/lib/python3.6/site-packages:$PYTHONPATH

If **Nemo** has installed correctly, then you should find its command line tools are available, for
example,

.. code-block::
   
   nemo -h
   
should display a helpful message about the command-line options for the main ``nemo`` command.
