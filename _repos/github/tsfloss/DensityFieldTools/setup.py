import os
from setuptools import setup

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "DensityFieldTools",
    version = "0.0.1",
    author = "Thomas Flöss",
    author_email = "tsfloss@gmail.com",
    description = ("Tools for manipulating density fields and measuring power spectra and bispectra"),
    license = "MIT",
    url = "https://github.com/tsfloss/DensityFieldTools",
    packages=['DensityFieldTools'],
    install_requires=[
          'numpy',"numba","tqdm","ducc0",
      ],
)