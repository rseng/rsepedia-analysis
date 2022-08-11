#!/usr/bin/env python

import re
import os
import sys
from setuptools import setup


setup(
    name="Cobra",
    version="1.0",
    author="Lindley Lentati",
    author_email="lindleylentati@gmail.com",
    url="https://gitlab.com/lindleylentati/Cobra",
    py_modules=["Cobra", "Candidate", "DatClass", "File"],
    description="Bayesian Pulsar Searching.",
    long_description=open("README.rst").read(),
    package_data={"": ["LICENSE"]},
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "License :: OSI Approved :: BSD License",
        "Intended Audience :: Developers",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Programming Language :: Python",
    ],
)
