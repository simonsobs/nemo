# -*- coding: iso-8859-1 -*-
#
# nemo install script

import os
import glob
from setuptools import setup
from setuptools import Extension
#from distutils.core import setup
#from distutils.extension import Extension
from Cython.Distutils import build_ext
import numpy
#import popen2

setup(name='nemo',
      version="1.0.dev0",
      url="https://acru.ukzn.ac.za/~mjh/nemo",
      author='Matt Hilton',
      author_email='hiltonm@ukzn.ac.za',
      classifiers=[],
      description='SZ cluster detection tool.',
      long_description="""Map filtering and SZ cluster detection and characterization pipeline.""",
      packages=['nemo'],
      package_data={'nemo': ['data/*']},
      scripts=['bin/nemo', 'bin/nemoMass', 'bin/nemoSelFn', 'bin/nemoMock', 'bin/nemoCatalogCheck', 'bin/nemoCosmo'],
      cmdclass={'build_ext': build_ext},
      ext_modules=[Extension("nemoCython", ["nemo/nemoCython.pyx"], include_dirs=[numpy.get_include()])],
      install_requires=["astropy >= 3.2",
                        "numpy >= 1.10",
                        "matplotlib >= 2.0",
                        "astLib >= 0.10",
                        "pixell >= 0.5",
                        "scipy >= 1.0",
                        "pillow",
                        "cython",
                        "PyYAML",
                        "colossus",
                        "mpi4py",
                        "colorcet >= 1.0"]
)
