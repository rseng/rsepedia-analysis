#!/usr/bin/env python

from setuptools import setup, find_packages

with open("README.md", "r") as f:
    long_description = f.read()

# Get common version number (https://stackoverflow.com/a/7071358)
import re
VERSIONFILE="mhealpy/_version.py"
verstrline = open(VERSIONFILE, "rt").read()
VSRE = r"^__version__ = ['\"]([^'\"]*)['\"]"
mo = re.search(VSRE, verstrline, re.M)
if mo:
    verstr = mo.group(1)
else:
    raise RuntimeError("Unable to find version string in %s." % (VERSIONFILE,))

setup(name='mhealpy',
      version = verstr,
      author='Israel Martinez-Castellanos',
      author_email='imc@umd.edu',
      url='https://gitlab.com/burstcube/HealpixMap',
      packages = find_packages(include=["mhealpy", "mhealpy.*"]),
      install_requires = ['healpy'],
      description = "Object-oriented healpy wrapper supporting multi-resolution maps",
      long_description = long_description,
      long_description_content_type="text/markdown",
      )

