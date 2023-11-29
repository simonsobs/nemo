.. image:: https://readthedocs.org/projects/nemo-sz/badge/?version=latest

**Nemo** is a millimeter-wave Sunyaev-Zel'dovich galaxy cluster and
compact source detection package, originally developed for the
`Atacama Cosmology Telescope <https://act.princeton.edu/>`_ project,
and capable of analyzing maps that will be produced by the
`Simons Observatory <https://simonsobservatory.org/>`_.

* **Documentation:** https://nemo-sz.readthedocs.io
* **License:** `BSD 3-Clause <https://github.com/simonsobs/nemo/blob/main/LICENSE>`_
* **Authors:** Matt Hilton, with contributions from Simone Aiola, David Alonso,
  Boris Bolliet, Matthew Hasselfield, Kevin Huffenberger, Toby Marriage, Niall MacCrann,
  Sigurd Naess, Cristóbal Sifón, and Íñigo Zubeldia (not all reflected on GitHub).
* **Installation:** ``pip install nemo-sz``
* **Support:** Please use the `GitHub issues page <https://github.com/simonsobs/nemo/issues>`_, 
  and/or contact `Matt Hilton <mailto:matt.hilton@mykolab.com>`_.

**Nemo** is written in `Python <https://www.python.org/>`_ and
provides `several modules <https://nemo-sz.readthedocs.io/en/latest/reference.html>`_ that
may be useful for analyzing ACT/SO data, in addition to the command-line programs provided
in the package.

**Nemo** is *not* the pipeline used for 
`Hasselfield et al. (2013) <http://adsabs.harvard.edu/abs/2013JCAP...07..008H>`_, but implements many
of the ideas presented there. *It is* the package that was used to produce the
`ACT DR3 cluster catalog <https://ui.adsabs.harvard.edu/abs/2018ApJS..235...20H/abstract>`_,
and the `ACT DR5 cluster catalog <https://ui.adsabs.harvard.edu/abs/2020arXiv200911043H/abstract>`_.
All ACT data products are available from `LAMBDA <https://lambda.gsfc.nasa.gov/product/act/actpol_prod_table.cfm>`_.

**Nemo** is under active development, and not all documentation or example config files are up to date
(some may contain references to files that are not yet publicly available). The package also contains
some experimental features that are not necessarily well tested.

If you need to acknowledge the use of **Nemo**, please cite
`Hilton et al. (2021) <https://ui.adsabs.harvard.edu/abs/2020arXiv200911043H/abstract>`_.
