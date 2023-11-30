# -*- coding: iso-8859-1 -*-
#
# nemo install script

import os
import glob
from setuptools import setup
from setuptools import Extension
import numpy
import versioneer

cmdclass=versioneer.get_cmdclass()

setup(name='nemo-sz',
      version=versioneer.get_version(),
      cmdclass=cmdclass,
      url="https://nemo-sz.readthedocs.io",
      author='Matt Hilton + Nemo Contributors',
      author_email='matt.hilton@wits.ac.za',
      packages=['nemo'],
      package_data={'nemo': ['data/*']},
      scripts=['bin/nemo', 'bin/nemoMass', 'bin/nemoMock', 'bin/nemoCatalogCheck', 'bin/nemoMask', 'bin/nemoModel', 'bin/nemoSpec'],
)
