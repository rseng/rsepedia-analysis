# Change Log
All notable changes to this project will be documented in this file.

The format is based on [Keep a Changelog](http://keepachangelog.com/)
and this project adheres to [Semantic Versioning](http://semver.org/).

## Unreleased

## [0.5.0] - 2020-07-07

### Fixed

* Compatible with Jupyter Lab 2.1.5 ([#32](https://github.com/eWaterCycle/jupyterlab_thredds/issues/32))

## [0.4.1] - 2019-09-02

### Fixed

* Added missing fetchers sub-package

## [0.4.0] - 2019-09-02

### Added

* ESGF free text search support ([#21](https://github.com/eWaterCycle/jupyterlab_thredds/issues/21))

### Fixed

* Compatible with Jupyter Lab 1.1.0 ([#28](https://github.com/eWaterCycle/jupyterlab_thredds/issues/28))

## [0.3.0] - 2018-12-05

### Added

* timeout config key to limit the crawling time

### Fixed

* THREDDS server v5 compound serviceType not expanded (#22)
* Crawling takes too long for large (>100 datasets) catalogs (#23)

### Changed

* Retrieve WMS layers in cell instead of while crawling (#23)
* Switched from thredds_crawler to siphon Python package (#22 + #23)
* Replaced workers config key with maxtasks

## [0.2.1] - 2018-11-29

### Fixed

* Skip misbehaving wms getcapabilities urls (#20)
* Row not disabled when open as not supported for row (#19)

### Changed

* Compatible with Jupyter Lab v0.35.4
* Pinned xarray dependency so ipyleaflet is satisfied
* Show crawling in progress
* Vertical scrollbar when datasets do not fit on screen

## [0.2.0] - 2018-09-03

### Changed

* Compatible with Jupyter Lab v0.34.7

### Fixed 

* Redirection loop

## [0.1.0] - 2018-05-04

Initial release
# jupyterlab_thredds

[![Build Status](https://travis-ci.org/eWaterCycle/jupyterlab_thredds.svg?branch=master)](https://travis-ci.org/eWaterCycle/jupyterlab_thredds)
[![SonarCloud Quality](https://sonarcloud.io/api/project_badges/measure?project=jupyterlab_thredds&metric=alert_status)](https://sonarcloud.io/dashboard?id=jupyterlab_thredds)
[![SonarCloud Coverage](https://sonarcloud.io/api/project_badges/measure?project=jupyterlab_thredds&metric=coverage)](https://sonarcloud.io/component_measures?id=jupyterlab_thredds&metric=coverage)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1241006.svg)](https://doi.org/10.5281/zenodo.1241006)

JupyterLab dataset browser for [THREDDS catalog](https://www.unidata.ucar.edu/software/thredds/v4.6/tds/catalog/index.html)

Can inject [iris](http://scitools.org.uk/iris/docs/latest/index.html)/[xarray](https://xarray.pydata.org)/[leaflet](https://github.com/jupyter-widgets/ipyleaflet) code cells into a Python notebook of a selected dataset to further process/visualize the dataset.

![screenshot](https://github.com/eWaterCycle/jupyterlab_thredds/blob/master/jupyterlab_thredds.gif "Screenshot")

## Prerequisites

* JupyterLab, `pip install jupyterlab`
* ipywidgets, `jupyter labextension install @jupyter-widgets/jupyterlab-manager`, requirement for ipyleaflet
* ipyleaflet, `jupyter labextension install jupyter-leaflet`, to load a WMS layer
* [iris](http://scitools.org.uk/iris/docs/latest/index.html), `conda install -c conda-forge iris`

## Installation

```bash
pip install jupyterlab_thredds
jupyter labextension install @ewatercycle/jupyterlab_thredds
```

## Usage

0. Start Jupyter lab with `jupyter lab`
1. In Jupyter lab open a notebook
2. Open the `THREDDS` tab on the left side.
3. Fill the catalog url
4. Press search button
5. Select how you would like to open the dataset, by default it uses [iris](http://scitools.org.uk/iris/docs/latest/index.html) Python package.
6. Press a dataset to insert code into a notebook

## Development

For a development install, do the following in the repository directory:

```bash
pip install -r requirements.txt
jlpm
jlpm build
jupyter labextension link .
jupyter serverextension enable --sys-prefix jupyterlab_thredds
```
(`jlpm` command is JupyterLab's pinned version of [yarn](https://yarnpkg.com/) that is installed with JupyterLab.)

To rebuild the package and the JupyterLab app:

```bash
jlpm build
jupyter lab build
```

Watch mode
```bash
# shell 1
jlpm watch
# shell 2
jupyter lab --ip=0.0.0.0 --no-browser --watch
```

## Release

To make a new release perform the following steps:
1. Update version in `package.json` and `jupyterlab_thredds/version.py`
2. Record changes in `CHANGELOG.md`
3. Make sure tests pass by running `jlpm test` and `pytest`
5. Commit and push all changes
6. Publish lab extension to npmjs with `jlpm build` and `jlpm publish --access=public`
7. Publish server extension to pypi with `python setup.py sdist bdist_wheel` and `twine upload dist/*`
8. Create GitHub release
9. Update DOI in `CITATION.cff`
