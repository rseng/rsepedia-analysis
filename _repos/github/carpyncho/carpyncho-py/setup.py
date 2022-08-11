#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (c) 2020, Juan B Cabral
# License: BSD-3-Clause
#   Full Text: https://github.com/carpyncho/carpyncho-py/blob/master/LICENSE


# =============================================================================
# DOCS
# =============================================================================

"""This file is for distribute carpyncho

"""


# =============================================================================
# IMPORTS
# =============================================================================

import os
import pathlib

from ez_setup import use_setuptools

use_setuptools()

from setuptools import setup


# =============================================================================
# CONSTANTS
# =============================================================================

PATH = pathlib.Path(os.path.abspath(os.path.dirname(__file__)))

REQUIREMENTS = [
    "numpy",
    "requests",
    "attrs",
    "pandas",
    "diskcache",
    "tqdm",
    "typer",
    "humanize",
    # formats
    "fastparquet",
    "pyarrow",
    "openpyxl",
]

with open(PATH / "README.md") as fp:
    LONG_DESCRIPTION = fp.read()

DESCRIPTION = "Python client for Carpyncho VVV dataset collection."

with open(PATH / "carpyncho.py") as fp:
    VERSION = (
        [line for line in fp.readlines() if line.startswith("__version__")][0]
        .split("=", 1)[-1]
        .strip()
        .replace('"', "")
    )


# =============================================================================
# FUNCTIONS
# =============================================================================


def do_setup():
    setup(
        name="carpyncho",
        version=VERSION,
        description=DESCRIPTION,
        long_description=LONG_DESCRIPTION,
        long_description_content_type="text/markdown",
        author="J.B. Cabral",
        author_email="jbc.develop@gmail.com",
        url="https://carpyncho-py.readthedocs.io/",
        license="3 Clause BSD",
        keywords=["astronomy", "vvv", "catalogs"],
        classifiers=(
            "Development Status :: 4 - Beta",
            "Intended Audience :: Education",
            "Intended Audience :: Science/Research",
            "License :: OSI Approved :: BSD License",
            "Operating System :: OS Independent",
            "Programming Language :: Python",
            "Programming Language :: Python :: 3.8",
            "Programming Language :: Python :: 3.9",
            "Programming Language :: Python :: 3.10",
            "Programming Language :: Python :: Implementation :: CPython",
            "Topic :: Scientific/Engineering",
        ),
        entry_points={"console_scripts": ["carpyncho=carpyncho:main"]},
        py_modules=["carpyncho", "ez_setup"],
        install_requires=REQUIREMENTS,
        include_package_data=True,
    )


if __name__ == "__main__":
    do_setup()
