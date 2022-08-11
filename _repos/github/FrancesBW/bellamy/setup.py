#! /usr/bin/env python
"""
Setup for BELLAMY
"""
import os
import sys
from setuptools import setup

reqs=['numpy>=1.15.4', 'astropy>=2.0, <3','scipy>=1.1.0','argparse>=1.2.1', 'psutil>=5.4.8', 'matplotlib>=2.2.3']

setup(
    name="bellamy",
    version="0.1",
    author="Frances Buckland-Willis",
    description="BELLAMY: A cross-matching package for the cynical astronomer",
    url="https://github.com/FrancesBW/belamy",
    long_description=['README.md'],
    packages=['functions'],
    install_requires=reqs,
    scripts=['scripts/bellamy'],
    classifiers=['Programming Language :: Python :: 2.7']
)
