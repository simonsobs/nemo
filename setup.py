# -*- coding: iso-8859-1 -*-
#
# nemo install script

import os
import glob
from setuptools import setup
from setuptools import Extension
from Cython.Distutils import build_ext
import numpy
import versioneer

cmdclass=versioneer.get_cmdclass()
cmdclass['build_ext']=build_ext

setup(name='nemo-sz',
      version=versioneer.get_version(),
      cmdclass=cmdclass,
      url="https://nemo-sz.readthedocs.io",
      author='Matt Hilton + Nemo Contributors',
      author_email='hiltonm@ukzn.ac.za',
      classifiers=['Development Status :: 4 - Beta',
                   'Environment :: Console',
                   'Intended Audience :: Science/Research',
                   'License :: OSI Approved :: BSD License',
                   'Natural Language :: English',
                   'Operating System :: POSIX',
                   'Programming Language :: Python',
                   'Topic :: Scientific/Engineering :: Astronomy'],
      description='Millimeter-wave galaxy cluster and source detection package.',
      long_description="""Millimeter-wave map filtering and Sunyaev-Zel'dovich galaxy cluster/source detection package. Originally developed for the Atacama Cosmology Telescope project.""",
      packages=['nemo'],
      package_data={'nemo': ['data/*']},
      scripts=['bin/nemo', 'bin/nemoMass', 'bin/nemoSelFn', 'bin/nemoMock', 'bin/nemoCatalogCheck', 'bin/nemoCosmo', 'bin/nemoModel', 'bin/nemoSpec'],
      ext_modules=[Extension("nemoCython", ["nemo/nemoCython.pyx"], include_dirs=[numpy.get_include()])],
      install_requires=["astropy >= 4.0",
                        "numpy >= 1.19",
                        "matplotlib >= 2.0",
                        "astLib >= 0.11.6",
                        "pixell >= 0.5",
                        "scipy >= 1.0",
                        "pillow",
                        "cython",
                        "PyYAML",
                        #"colossus",
                        #"pyccl >= 2.1",
                        #"mpi4py",
                        "colorcet >= 1.0",
                        "mahotas >= 1.4"]
)
