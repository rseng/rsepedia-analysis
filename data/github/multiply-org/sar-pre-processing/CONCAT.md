<img alt="MULTIPLY" align="right" src="https://raw.githubusercontent.com/multiply-org/sar-pre-processing/master/docs/images/multiply_multi_colour.png" />

# SenSARP

[![Build Status](https://www.travis-ci.com/multiply-org/sar-pre-processing.svg?branch=master)](https://travis-ci.com/McWhity/sar-pre-processing)
[![Documentation Status](https://readthedocs.org/projects/multiply-sar-pre-processing/badge/?version=master)](https://multiply-sar-pre-processing.readthedocs.io/en/master/?badge=master)

This repository contains the functionality SenSARP used within the MULTIPLY main platform.
The [SenSARP specific documentation](https://multiply-sar-pre-processing.readthedocs.io/en/master/) is hosted on ReadTheDocs. It is part of the [MULTIPLY core documentation](http://multiply.readthedocs.io/).
Please find the pdf version of the SenSARP documentation [here](https://multiply-sar-pre-processing.readthedocs.io/_/downloads/en/master/pdf/) and for the MULTIPLY platform [here](https://readthedocs.org/projects/multiply/downloads/pdf/latest/).
SenSARP is a pipeline to pre-process Sentinel-1 SLC data by using ESA SNAP Sentinel-1 Toolbox.
Expert users can adapt the default pre-processing chain to their needs and benefit from functionalities provided by SenSARP.

## Statement of need

Sentinel-1 satellites will provide continuous free available microwave remote sensing data of the entire globe at least until the end of 2030.
Furthermore, ESA is not only providing Sentinel satellite images (e.g. Sentinel-1, Sentinel-2, Sentinel-3) but they also developed free open source toolboxes (Sentinel-1, 2, 3 toolboxes) for scientific exploitation.
The toolboxes can be accessed and used via the Sentinel Application Platform (SNAP).
SNAP offers a graphical interface were expert users can develop different processing schemes and apply them on the satellite images.
Although, Sentinel-1 satellite data and a processing software are freely available, the usage of the data is mainly limited to expert users in the field of microwave remote sensing as different pre-processing steps need to be applied before using Sentinel-1 images.

SenSARP was developed to provide a push-button option to easily apply a rigid pre-processing pipeline with sensible defaults to a Sentinel-1 Level 1 SLC time series data as well as single Sentinel-1 Level 1 SLC images.
Thus, non-expert users in the field of pre-processing microwave data are able to use radiometric and geometric corrected sigma nought backscatter data for their specific applications.
Beside a rigid pre-processing pipeline, SenSARP provides filter options to retrieve only images of a specific year or images that contain a specific area of interest from a stack of downloaded Sentinel-1 data.
Furthermore, the default processing scheme of SenSARP can handle if an area of interest is contained in two tiles of the same swath (due to storage reasons data of one Sentinel-1 satellite swath is provided by ESA within different tiles).
Additionally, SenSARP checks if within a stack of Sentinel-1 images, one specific image was multiple processed by ESA and uses the newest.

For expert users, SenSARP provides the possibility to automate their pre-processing on a large scale by either modifying the default pre-processing scheme (modification of xml graph pre_processing_step1.xml) or create their own pre-processing scheme (create a new xml graph) with the graph builder of the SNAP software.
They can benefit from the filter options, the default pre-processing step 2 (co-registration of images) and the SenSARP functions to stack all processed and co-registered images within a netCDF file with additional image information e.g. satellite name, relative orbit and orbitdirection.

## Content of this repository

* `docs/` - The auto generated documentation
* `recipe/` Conda installation recipe
* `sar_pre_processing/` - The main sar pre processing software package
* `test/` - The test package.
* `AUTHORS.rst` - Author information.
* `CHANGES.md` - Package change log.
* `LICENSE.rst` - License of software in repository.
* `README.md` - Readme.
* `environmental.yml` - Requirements.
* `setup.py` - main build script, to be run with Python 3.6

## How to install

The first step is to clone the latest code and step into the check out directory::

    git clone https://github.com/multiply-org/sar-pre-processing.git
    cd sar-pre-processing

### Installation with Conda

Download and install [Anaconda](https://www.anaconda.com/products/individual) or [Miniconda](https://docs.conda.io/en/latest/miniconda.html). Anaconda/Miniconda installation instructions can be found [Anaconda](https://conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent)

To install all required modules, use

    conda env create --prefix ./env --file environment.yml
    conda activate ./env # activate the environment

To install SenSARP into an existing Python environment, use::

    python setup.py install

To install for development, use::

    python setup.py develop

### Installation with virtualenv and python

Install system requirements

    sudo apt install python3-pip python3-tk python3-virtualenv python3-venv virtualenv

Create a virtual environment

    virtualenv -p /usr/bin/python3 env
    source env/bin/activate # activate the environment
    pip install --upgrade pip setuptools # update pip and setuptools

To install SenSARP into an existing Python environment, use::

    python setup.py install

To install for development, use::

    python setup.py develop

GDAL package needs to be installed too

    sudo apt install gdal-bin libgdal-dev

    python -m pip install pygdal=="`gdal-config --version`.*"


### Further information

Please see the [environment file](environment.yml) for a list of all installed dependencies during the installation process.
Additionally, ESA's SNAP Sentinel-1 Toolbox (Version >8.0.3) has to be installed prerequisite. The Software can be downloaded [here](http://step.esa.int/main/download/snap-download/). To install the SNAP toolbox, open a terminal window and use

    bash esa-snap_sentinel_unix_8_0.sh

SenSARP uses only functionalities of the Sentinel-1 Toolbox.
Currently, only SNAP version 8.0 can be downloaded from the website.
To update SNAP to a version >8.0.3 please start the SNAP software.
You will be asked if you want to search for update.
Please search for updates and install all updates.
After the updates are installed, you need to restart SNAP to initialize the updates correctly.
SNAP Toolbox need libgfortran for specific operations but currently libgfortran is not installed during the installation process of SNAP, therefore you might use

    sudo apt-get install gfortran

## Usage

For usage checkout the [juypter notebook](https://nbviewer.jupyter.org/github/multiply-org/sar-pre-processing/tree/master/docs/notebooks/)

## Documentation

We use [Sphinx](http://www.sphinx-doc.org/en/stable/rest.html) to generate the documentation of SenSARP on [ReadTheDocs](https://multiply.readthedocs.io/). The SenSARP documentation is available [here](https://multiply-sar-pre-processing.readthedocs.io/en/latest/)

## Support, contributing and testing
Please contribute using [Github Flow](https://guides.github.com/introduction/flow/). Create a branch, add commits, and [open a pull request](https://github.com/multiply-org/sar-pre-processing/issues/new).

### Reporting bugs
If you find a bug in SenSARP, please open an new [issue](https://github.com/multiply-org/sar-pre-processing/issues/new) and tag it "bug".

### Suggesting enhancements
If you want to suggest a new feature or an improvement of a current feature, you can submit this on the [issue tracker](https://github.com/multiply-org/sar-pre-processing/issues/new) and tag it "enhancement".

### Testing

The package is currently tested for Python >= 3.6 on Unix-like systems.
To run unit tests, execute the following line from the root of the repository:

    pytest

## Authors

[Authors](AUTHORS.rst)

## License

This project is licensed under the GPLv3 License - see the [LICENSE.rst](LICENSE.rst) file for details.
# Changelog SenSARP module

## [0.1] - 2021-04-20

initial version

## [0.2] - 2022-01-18 (JOSS Paper)

### JOSS Paper

- extensive revision

### Documentation

- add installation guide using virtualenv and pip
- revise introduction
- add Statement of need
- add explanations of created files and used abbreviations
- name change from "MULTIPLY SAR pre-processing" to "SenSARP"
- add notebooks (use cases)
    - default_process_single_image.ipynb
    - default_process_time_series.ipynb
    - use_user_defined_graphs.ipynb
- remove notebook running_test_application.ipynb

### Software

- add more explanation to config file
- add more documentation to functions of class SARPreProcessor
- add functionality that to user defined xml-graphs can be used
- add functionality that a single file can be processed
- class NetcdfStackCreator was partly rewritten
- functionality of shell file (solve_projection_problem.sh) was included within python code
- bug fixes
---
title: 'SenSARP: A pipeline to pre-process Sentinel-1 SLC data by using ESA SNAP Sentinel-1 Toolbox'
tags:
  - Python
  - Sentinel-1
  - SAR
  - SNAP
authors:
  - name: Thomas Weiß
    orcid: 0000-0001-5278-8379
    affiliation: 1
  - name: Tonio Fincke
    affiliation: 2
affiliations:
 - name: Department of Geography, Ludwig-Maximilians-Universität München, 80333 Munich, Germany
   index: 1
 - name: Brockmann Consult GmbH, 21029 Hamburg, Germany;
   index: 2
date: 20 April 2021
bibliography: paper.bib
---

# Summary

The Sentinel-1 mission consists of two polar-orbiting satellites acquiring Synthetic Aperture Radar data (SAR) at C-band (frequency of 5.405 GHz) with a revisit time of 6 days.
The SAR data is distributed free of charge via the Copernicus Open Access Hub (https://scihub.copernicus.eu/) by European Space Agency (ESA) and the European Commission.
Large archives are also provided by Data and Information Access Services (DIAS) which serve the purpose to facilitate the access and use of Sentinel Data.
Due to the specific imaging geometry of the radar system, the acquired radar data contains different radiometric and geometric distortions.
The radiometric quality is affected by spreading loss effect, the non-uniform antenna pattern, possible gain changes, saturation, and speckle noise.
Geometric distortions such as foreshortening, layover or shadowing effects are based on the side looking radar acquisition system.
To account for these radiometric and geometric distortions, the Sentinel-1 Level 1 data has to be corrected radiometrically and geometrically before the data can be used for further analysis or within third party applications.
Therefore, either an automatic or manual pre-processing of Sentinel-1 images is needed.

# Statement of need

Sentinel-1 satellites will provide continuous free available microwave remote sensing data of the entire globe at least until the end of 2030.
Furthermore, ESA is not only providing Sentinel satellite images (e.g. Sentinel-1, Sentinel-2, Sentinel-3) but they also developed free open source toolboxes (Sentinel-1, 2, 3 toolboxes) for scientific exploitation.
The toolboxes can be accessed and used via the Sentinel Application Platform (SNAP).
SNAP offers a graphical interface were expert users can develop different processing schemes and apply them on the satellite images.
Although, Sentinel-1 satellite data and a processing software are freely available, the usage of the data is mainly limited to expert users in the field of microwave remote sensing as different pre-processing steps need to be applied before using Sentinel-1 images.

SenSARP was developed to provide a push-button option to easily apply a rigid pre-processing pipeline with sensible defaults to a Sentinel-1 Level 1 SLC time series data as well as single Sentinel-1 Level 1 SLC images.
Thus, non-expert users in the field of pre-processing microwave data are able to use radiometric and geometric corrected sigma nought backscatter data for their specific applications.
Beside a rigid pre-processing pipeline, SenSARP provides filter options to retrieve only images of a specific year or images that contain a specific area of interest from a stack of downloaded Sentinel-1 data.
Furthermore, the default processing scheme of SenSARP can handle if an area of interest is contained in two tiles of the same swath (due to storage reasons data of one Sentinel-1 satellite swath is provided by ESA within different tiles).
Additionally, SenSARP checks if within a stack of Sentinel-1 images, one specific image was multiple processed by ESA and uses the newest.

For expert users, SenSARP provides the possibility to automate their pre-processing on a large scale by either modifying the default pre-processing scheme (modification of xml graph pre_processing_step1.xml) or create their own pre-processing scheme (create a new xml graph) with the graph builder of the SNAP software.
They can benefit from the filter options, the default pre-processing step 2 (co-registration of images) and the SenSARP functions to stack all processed and co-registered images within a netCDF file with additional image information e.g. satellite name, relative orbit and orbit direction.

# Method

This Python package generates a file list of to be processed Sentinel-1 images (already downloaded and stored in a specific folder) based on different user defined criteria (specific year, area of interest).
Additionally, specific cases of repeatedly processed data are handled, as sometimes Sentinel-1 data were initially processed multiple times and stored under similar names on the Copernicus Open Access Hub. Also, cases where Sentinel-1 data within the user-defined area of interest might be stored in consecutive tiles are considered.

Based on the generated file list the default processing pipeline of the Python package applies a pre-processing chain to Sentinel-1 Single Look Complex (SLC) time series or single images to generate radiometrically and geometrically corrected sigma nought backscatter values.
Furthermore, if a time series is processed the images are co-registered and additional output files of multi-temporal speckle filtered data are generated.
In addition, a single speckle filter instead of a multi-temporal one is applied as well and the output will be stored as a separate layer.
To pre-process the images, the Python package uses the GPT (Graph Processing Tool) of SNAP to execute different operators provided by the Sentinel-1 Toolbox.
The Sentinel Toolbox is available for download at http://step.esa.int/, its source code is available in the senbox-org organization on GitHub.
Each of these operators performs a pre-processing step. The operators can be chained together to form a graph, which is used by the Python package to run on the Sentinel-1 data using the Graph Processing Framework (GPF). The graphs are stored in xml-files. Users may change the graphs by modifying the files directly or via the Sentinel Toolbox.
User Guides to show how the GPF can be used are provided here: https://senbox.atlassian.net/wiki/spaces/SNAP/pages/70503053/Processing.

After the pre-processing the resulting radiometrically and geometrically corrected images are stored for further usage within a NetCDF4 stack file.
The processing workflow was developed and optimized to use a Sentinel-1 time series of pre-processed sigma nought backscatter values to retrieve biophysical land surface parameters by the use of radiative transfer models.
The sigma nought backscatter values provided by the default workflow of SenSARP might be used in other applications like flood risk analysis, monitoring land cover changes or monitoring global food security but it has to be mentioned that different applications have different demands and therefore, slight adjustments of the default workflow might be required.
In the future, many more new products and operational third party services based on consistent Sentinel-1 time series might be developed.

# Applications

This Python package was developed within the Horizon 2020 project called MULTIscale SENTINEL land surface information retrieval Platform (MULTIPLY) (http://www.multiply-h2020.eu/, https://cordis.europa.eu/project/id/687320, https://multiply.obs-website.eu-de.otc.t-systems.com).
Furthermore, data processed by this package is used within Sentinel-Synergy-Study S3 project (https://www.researchgate.net/project/Sentinel-Synergy-Study-S3).
In addition, the Python code was used to process Sentinel-1 time series images for the detection and analysis of temporary flooded vegetation [@tsyganskaya_detection_2018; @tsyganskaya_flood_2019] and for the evaluation of different radiative transfer models for microwave backscatter estimation of wheat fields [@weis_evaluation_2020].

# Other available Python software packages using ESA's SNAP software to pre-process SAR data

The ESA's SNAP toolbox has been written in Java. For Python users the developers provide a Python interface called Snappy. However, the Snappy interface is lacking in terms of installation, processing performance and usability. Hence, the remote sensing community developed different wrappers (e.g. SenSARP, snapista or pyroSAR) to use SNAP processing functionalities by utilizing the SNAP Graph Processing Tool (GPT).

## snapista

Snapista (https://snap-contrib.github.io/snapista/index.html) targets mainly experts remote sensing users with Python programming skills.
It provides access to the processing operators of all toolboxes (e.g. Sentinel-1, Sentinel-2 or Sentinel-3) within SNAP.
Expert users can generate processing graphs and execute their generated graphs in a pure Pythonic way.
Guidelines about which processing steps are needed for different applications, or about which processing steps can or have to be combined, are not provided yet.
Establishing guidelines about how to process different satellite data for different applications is not an easy task to do and would exceed the goal of snapista as a Python wrapper for the SNAP software.
Summarizing, snapista provides access to all SNAP toolboxes (not just to Sentinel-1 Toolbox) via Python. But as it provides no default processing chains, snapista will be primarily usable by expert remote sensing users.
The advantage of snapista is the accessibility of processing operators for SAR and optical data.

## pyroSAR

PyroSAR (https://pyrosar.readthedocs.io/en/latest/index.html) is a Python library which provides a Python wrapper to SAR pre-processing software SNAP and GAMMA [@wegnuller_sentinel-1_2016; @werner_gamma_2000].
The library provides utilities to read and store metadata information of downloaded satellite data within a database.
Furthermore, pyroSAR provides access to processing operators of SNAP and GAMMA.
A default workflow with different user options is provided to process single or time-series Sentinel-1 images.
After executing the default processing workflow radiometric and geometric corrected gamma nought backscatter data are provided in Geotiff format [@truckenbrodt_pyrosar_2019].
The processed images can also be stored within an Open Data Cube.
For expert users which might want to use a different processing workflow pyroSAR provides an option to create SNAP xml-workflows and execute them with the GPT.
Summarizing, pyroSAR provides a similar push-button option to process Sentinel-1 data with a slightly different default workflow (pyroSAR: no temporal speckle filter, gamma nought backscatter output in Geotiff format) than SenSARP (SenSARP: temporal speckle filter, sigma nought backscatter output in netCDF format).
PyroSAR, as a more complex library than SenSARP, provides on the one hand more changeable parameters within the processing workflow but on the other hand the usability for non-expert users might be narrowed compared to SenSARP.
An advantage of SenSARP, especially for non-expert users, might be the provision of background information (theory/purpose) of the different pre-processing steps within the documentation.

# Acknowledgements

The project leading to this application has received funding from the European Union’s Horizon 2020 research and innovation program under Grant Agreement No. 687320.
We would like to thank Alexander Löw and Philip Marzahn for guiding discussions that lead to this publication.
We also would like to thank Thomas Ramsauer for discussions and suggestions.
<!-- for providing comments on the manuscript -->
The author also wishes to thank the reviewers and editors for their efforts and for their helpful comments to improve this paper and the software package.

# References
