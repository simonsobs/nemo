# -*- coding: iso-8859-1 -*-
#
# nemo install script

import os
import glob
from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
import popen2

setup(name='nemo',
      version="git",
      url=None,
      author='Matt Hilton',
      author_email='hiltonm@ukzn.ac.za',
      classifiers=[],
      description='SZ cluster detection tool.',
      long_description="""Map filtering and SZ cluster detection and characterization pipeline.""",
      packages=['nemo'],
      package_data={'nemo': ['data/*']},
      scripts=['bin/nemo', 'bin/nemoFakeSources', 'bin/nemomontage', 'bin/nemostack-Planck'],
      cmdclass={'build_ext': build_ext},
      ext_modules=[Extension("nemoCython", ["nemo/nemoCython.pyx"])]
)
