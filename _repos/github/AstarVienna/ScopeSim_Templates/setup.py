#!/usr/bin/env python3
"""
ScopeSim_templates: Helper functions for making ScopeSim Source objects
=======================================================================

How to compile and put these on pip::

    $ python setup.py sdist bdist_wheel
    $ twine upload dist/*

"""

from datetime import datetime
from distutils.core import setup
from setuptools import find_packages

# Version number
with open('scopesim_templates/version.py') as f:
    __version__ = f.readline().split("'")[1]

with open("README.rst", "r", encoding='utf-8') as fh:
    long_description = fh.read()


def setup_package():
    setup(name='ScopeSim_Templates',
          version=__version__,
          description="On-sky source templates for ScopeSim",
          long_description=long_description,
          author="A* Vienna",
          author_email="astar.astro@univie.ac.at",
          url="https://github.com/AstarVienna/ScopeSim_Templates",
          package_dir={'scopesim_templates': 'scopesim_templates'},
          packages=find_packages(),
          include_package_data=True,
          install_requires=["numpy>=1.16",
                            "scipy>0.17",
                            "astropy>1.1.2",
                            "matplotlib>1.5.0",
                            "requests>2.0",
                            "pyyaml>5",
                            "synphot>0.1",
                            "scopesim",
                            "pyckles",
                            "spextra",
                            ],
          classifiers=["Programming Language :: Python :: 3",
                       "License :: OSI Approved :: MIT License",
                       "Operating System :: OS Independent",
                       "Intended Audience :: Science/Research",
                       "Topic :: Scientific/Engineering :: Astronomy", ]
          )


if __name__ == '__main__':
    setup_package()
