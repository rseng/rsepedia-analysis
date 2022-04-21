#!/usr/bin/env python
# -*- coding: utf-8 -*-
from os.path import exists, dirname, realpath
from setuptools import setup, find_packages
import sys

author = u"PiSCAT developers"
# authors in alphabetical order
authors = [
    "Houman Mirzaalian",
    "Reza Gholami",
]
description = 'Interferometric scattering microscopy with Python'
name = 'piscat'
year = "2021"


sys.path.insert(0, realpath(dirname(__file__))+"/"+name)
try:
    from piscat import version
except BaseException:
    version = "unknown"

setup(
    name=name,
    author=author,
    author_email='houman.mirzaalian-dastjerdi@mpl.mpg.de',
    project_urls={
    'Documentation': 'https://piscat.readthedocs.io/',
    'Source': 'https://github.com/SandoghdarLab/PiSCAT/',
    'Tracker': 'https://github.com/SandoghdarLab/PiSCAT/issues'
    },
    version=version,
    packages=find_packages(),
    package_dir={name: name},
    include_package_data=True,
    license="GPL v3",
    description=description,
    long_description=open('README.md').read() if exists('README.md') else '',
    long_description_content_type="text/markdown",
    install_requires=[
                        "numpy>=1.19.0",
                        "scipy>=0.14.0",
                        "matplotlib",
                        "opencv-python",
                        "pandas",
                        "PyQt6",
                        "Pyside6",
                        "pyqtgraph",
                        "PyQtWebEngine",
                        "scikit_learn",
                        "scikit-image",
                        "tqdm",
                        "joblib",
                        "psutil",
                        "ipywidgets",
                        "networkx",
                        "numba",
                        "trackpy",
                        "requests",
                        "wget",
                        "h5py",
                        "imageio-ffmpeg",
                        "astropy"
                      ],
    # not to be confused with definitions in pyproject.toml [build-system]
    setup_requires=["pytest-runner"],
    python_requires=">=3.8",
    tests_require=["pytest", "pytest-mock"],
    keywords=[],
    classifiers=['Operating System :: OS Independent',
                 'Programming Language :: Python :: 3',
                 'Topic :: Scientific/Engineering :: Visualization',
                 'Intended Audience :: Science/Research',
                 ],
    platforms=['ALL'],
)