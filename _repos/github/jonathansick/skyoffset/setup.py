#!/usr/bin/env python
# Licensed under a 3-clause BSD style license - see LICENSE.rst

import glob
import os
import sys

from setuptools import setup

# Set affiliated package-specific settings
PACKAGENAME = 'skyoffset'
DESCRIPTION = 'Toolkit for sky offset solving and mosaicing from ANDROIDS'
LONG_DESCRIPTION = ''
AUTHOR = 'Jonathan Sick'
AUTHOR_EMAIL = 'jonathansick@mac.com'
LICENSE = 'BSD'
URL = 'http://jonathansick.ca'

# VERSION should be PEP386 compatible (http://www.python.org/dev/peps/pep-0386)
VERSION = '0.1.dev'

# Indicates if this version is a release version
RELEASE = 'dev' not in VERSION

# Treat everything in scripts except README.rst as a script to be installed
scripts = [fname for fname in glob.glob(os.path.join('scripts', '*'))
           if os.path.basename(fname) != 'README.rst']

setup(name=PACKAGENAME,
      version=VERSION,
      description=DESCRIPTION,
      scripts=scripts,
      requires=['astropy', 'numpy', 'cython', 'matplotlib', 'shapely',
                'Polygon2', 'pymongo', 'pyregion'],
      install_requires=['astropy', 'cython'],
      provides=[PACKAGENAME],
      author=AUTHOR,
      author_email=AUTHOR_EMAIL,
      license=LICENSE,
      url=URL,
      long_description=LONG_DESCRIPTION,
      zip_safe=False)
