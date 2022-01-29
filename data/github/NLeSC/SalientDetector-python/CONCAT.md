# Change Log
All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/). 

## 1.0.3 - 2016-09-28
## Changed
- Fix setup.py for PyPi installation

## 1.0.2 - 2016-09-28
## Changed
- Fix setup.py for PyPi installation

## 1.0.1 - 2016-09-28
## Changed
- Changed setup.py so that installation is completely possible through pip

## 1.0.0 - 2016-07-15
### Added
- The first release of the salient regions detectors software in Python. This implementation is (almost) equivalent to 
the original [MATLAB code](https://github.com/NLeSC/SalientDetector-matlab) and corresponds to the MATLAB's repo initial release. 
# Python software for image processing
[![Build Status](https://travis-ci.org/NLeSC/SalientDetector-python.svg?branch=master)](https://travis-ci.org/NLeSC/SalientDetector-python) [![Codacy Badge](https://api.codacy.com/project/badge/Grade/1c9f59fcbc6d48bbb35addc7d51e0bf1)](https://www.codacy.com/app/d-vankuppevelt/SalientDetector-python?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NLeSC/SalientDetector-python&amp;utm_campaign=Badge_Grade) [![Codacy Badge](https://api.codacy.com/project/badge/Coverage/1c9f59fcbc6d48bbb35addc7d51e0bf1)](https://www.codacy.com/app/d-vankuppevelt/SalientDetector-python?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=NLeSC/SalientDetector-python&amp;utm_campaign=Badge_Coverage) [![DOI](https://zenodo.org/badge/55982608.svg)](https://zenodo.org/badge/latestdoi/55982608)

This folder contains a Python  implementation of the Salient Region Detector code as part of the [image processing part of eStep](https://www.esciencecenter.nl/technology/expertise/computer-vision). The software conforms with the [eStep standards](https://github.com/NLeSC/estep-checklist).

The original MATLAB implementation can be found at [this repository](https://github.com/NLeSC/SalientDetector-matlab)

Documentation can be found on [Read the Docs](http://salientdetector-python.readthedocs.io/).

The repository contains the following sub-folders:

## Notebooks
Several iPython notebooks testing and illustrating major functionality.

## salientregions
The module for salient region detection functionality.

## tests
Unit tests for the code in salientregions.

# Installation
## Prerequisites
* Python 2.7 or 3.5
* pip (8.1.2)

## Installing the package (via pip)

`pip install salientregions`

Afterwards you can import it in Python:

```python
import salientregions as sr
```

## Installing the package manually
Clone or download repo. Install requirements from requirements.txt file:

`pip install -r requirements.txt`

Install the package itself:

`pip install .`

To perform tests:

`nosetests test`

## Manually install OpenCV
OpenCV should be installed via pip, but if you have an older version of pip, you might need to  install it manually. There is two ways to install OpenCV:
  * If you're using Conda, you can install OpenCV with the following command:

  `conda install -c https://conda.anaconda.org/menpo opencv3`

  * Otherwise, follow the instructions to [install OpenCV 3.1.0 manually](http://opencv-python-tutroals.readthedocs.org/en/latest/py_tutorials/py_setup/py_table_of_contents_setup/py_table_of_contents_setup.html#py-table-of-content-setup)


# Getting started
The source code documentation can be found [here](http://salientdetector-python.readthedocs.io/)

This code makes heavily use of the OpenCV library, so in order to understand how the code works, it helps to have a look at the [OpenCV Documentation](http://docs.opencv.org/3.1.0/).

## Images
In OpenCV, images are represented as numpy arrays. Grayscale images are represented by a 2-dimensional array. Color images have a third dimension for the color channel. The Salient Region Detector has a few simplifying assumptions:
* Color images have BGR channels
* Images are assumed to be 8-bit. This is also the case for binary images, so they only have values of 0 and 255.

## Detector object
The complete functionality of the salient region detectors are found in the Detector object. The SalientDetector implements DMSR detection, and MSSRDetector implements MSSR detection (see referred papers for more information about these algorithms).
An example of how to use the Detector can be found in [this iPython Notebook](https://github.com/NLeSC/SalientDetector-python/blob/master/Notebooks/DetectorExample.ipynb).

# Contributing
If you want to contribute to the code, please have a look at the [eStep standarts](https://github.com/NLeSC/estep-checklist).

We use numpy-style code documentation.

# References
Ranguelova, E.
[A Salient Region Detector for Structured Images.](http://ieeexplore.ieee.org/document/7945643/)
Proceedings of International Conference on Computer Systems and Applications (AICCSA), 2016, p.1–8,
DOI: 10.1109/AICCSA.2016.7945643 

##Salient Region Detectors software.

The prototype software is in Software/MATLAB/AffineRegions/Detectors
