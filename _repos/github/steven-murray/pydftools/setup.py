#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""The setup script."""

from setuptools import setup, find_packages
import io
import os
import re

with open("README.rst") as readme_file:
    readme = readme_file.read()

with open("HISTORY.rst") as history_file:
    history = history_file.read()

requirements = [
    "scipy",
    "numpy>=1.6.2",
    "Click>=6.0",
    "attrs>=17.0",
    "cached_property",
    "chainconsumer",
    "matplotlib"
    # TODO: put package requirements here
]


def read(*names, **kwargs):
    with io.open(
        os.path.join(os.path.dirname(__file__), *names),
        encoding=kwargs.get("encoding", "utf8"),
    ) as fp:
        return fp.read()


def find_version(*file_paths):
    version_file = read(*file_paths)
    version_match = re.search(r"^__version__ = ['\"]([^'\"]*)['\"]", version_file, re.M)
    if version_match:
        return version_match.group(1)
    raise RuntimeError("Unable to find version string.")


setup_requirements = [
    "pytest-runner",
    # TODO(steven-murray): put setup requirements (distutils extensions, etc.) here
]

test_requirements = [
    "pytest",
    # TODO: put package test requirements here
]

setup(
    name="pydftools",
    version=find_version("pydftools", "__init__.py"),
    description="A pure-python port of the dftools R package.",
    long_description=readme + "\n\n" + history,
    author="Steven Murray",
    author_email="steven.murray@curtin.edu.au",
    url="https://github.com/steven-murray/pydftools",
    packages=find_packages(include=["pydftools"]),
    entry_points={"console_scripts": ["pydftools=pydftools.cli:main"]},
    include_package_data=True,
    install_requires=requirements,
    license="MIT license",
    zip_safe=False,
    keywords="pydftools",
    classifiers=[
        "Development Status :: 2 - Pre-Alpha",
        "Intended Audience :: Developers",
        "License :: OSI Approved :: MIT License",
        "Natural Language :: English",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.6",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.3",
        "Programming Language :: Python :: 3.4",
        "Programming Language :: Python :: 3.5",
    ],
    test_suite="tests",
    tests_require=test_requirements,
    setup_requires=setup_requirements,
)
