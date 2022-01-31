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
* Thomas weiß <"weiss.thomas@lmu.de">
* Tonio Fincke <"tonio.fincke@brockmann-consult.de">
Explanation of config file
----------------------------------

.. code-block:: yaml

    # Sample config file
    #====================
    ## Necessary parameters
    #-----------------------
    ### Input folder with SAR data (zip format)
    input_folder: '/media/test/Desktop/data' # values: string

    ### Output folder to store the pre-processed data
    output_folder: '/media/test/Desktop/data' # values: string

    ### Location of SNAP's graph-processing-tool
    gpt: /home/tweiss/snap/bin/gpt

    ## Optional parameters
    #-----------------------
    ### Year of interest (only images of the specified year will be processed)
    #~~~~~~~~~~~~~~~~~~~~
    year: 2021 # values: integer

    ### Area of interest (only images containing the specified year will be processed)
    #~~~~~~~~~~~~~~~~~~~
    region:
      subset: 'yes' # values: 'no' or 'yes'
      ul:
        lat: 48.40 # values: float
        lon: 11.60 # values: float
      lr:
        lat: 48.10 # values: float
        lon: 11.90 # values: float

    ### Used parameters of multi-temporal speckle filter
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    speckle_filter:
      multi_temporal:
        apply: 'yes' # values: 'no' or 'yes'
        files: '5' # values: integer

    ### Used parameter of incidence normalization (default angle is 35°)
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    normalisation_angle: 35 # values: float

    ### Single file processing (if filelist contains only one file single_file option is automatically set to 'yes')
    #~~~~~~~~~~~~~~~~~~~~~~~~~
    single_file: 'yes' # values: 'no' or 'yes'

    ### Usage of user defined xml graphs
    #~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    #### Specification if user defined xml graph should be used
    use_user_definde_graphs: 'no' # values: 'no' or 'yes'

    #### Location of user defined xml files for processing
    xml_graph_path: /media/test/Desktop/sar-pre-processing/sar_pre_processing/user_defined_graphs

    ##### file names of user defined xml graphs
    pre_process_step1: example_of_expert_user.xml




.. _changes:
.. include:: ../CHANGES.md
Authors
=======
.. _authors:
.. include:: ../AUTHORS.rst
.. _Introduction:

Introduction
============
The Sentinel-1 mission consists of two polar-orbiting satellites acquiring Synthetic Aperture Radar data (SAR) at C-band (frequency of 5.405 GHz) with a revisit time of 6 days.
The SAR data is distributed free of charge via the Copernicus Open Access Hub (https://scihub.copernicus.eu/) by European Space Agency (ESA) and the European Commission.
Large archives are also provided by Data and Information Access Services (DIAS) which serve the purpose to facilitate the access and use of Sentinel Data.
Due to the specific imaging geometry of the radar system, the acquired radar data contains different radiometric and geometric distortions.
The radiometric quality is affected by spreading loss effect, the non-uniform antenna pattern, possible gain changes, saturation, and speckle noise.
Geometric distortions such as foreshortening, layover or shadowing effects are based on the side looking radar acquisition system.
To account for these radiometric and geometric distortions, the Sentinel-1 Level 1 data has to be corrected radiometrically and geometrically before the data can be used for further analysis or within third party applications.
Therefore, either an automatic or manual pre-processing of Sentinel-1 images is needed.

Statement of need
------------------
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


Getting Started
------------------
Please find instructions on how to download and install SenSARP in the :ref:`Installation` section.

Support, contributing and testing
----------------------------------
Please contribute using `Github Flow <https://guides.github.com/introduction/flow/>`_. Create a branch, add commits, and `open a pull request <https://github.com/multiply-org/sar-pre-processing/issues/new>`_.

Reporting bugs
~~~~~~~~~~~~~~~
If you find a bug in SenSARP, please open an new `issue <https://github.com/multiply-org/sar-pre-processing/issues/new>`_ and tag it "bug".

Suggesting enhancements
~~~~~~~~~~~~~~~
If you want to suggest a new feature or an improvement of a current feature, you can submit this on the `issue tracker <https://github.com/multiply-org/sar-pre-processing/issues/new>`_ and tag it "enhancement".

Testing
~~~~~~~~~~~~~~~
The package is currently tested for Python >= 3.6 on Unix-like systems.
To run unit tests, execute the following line from the root of the repository:

.. code:: bash

   pytest

.. _ProcessingChain:

Default Pre-Processing Chain
============================

The default pre-processing chain was developed to process Sentinel-1 SLC data to provide geometric and radiometric corrected Sigma nought backscatter values. Among other things the processed data can be used within different radiative transfer models (e.g. Integration Equation Model :cite:`fung_backscattering_1992`, Oh Model :cite:`oh_empirical_1992`, :cite:`yisok_oh_quantitative_2004`, Dubois Model :cite:`dubois_measuring_1995`, Water Cloud Model :cite:`attema_vegetation_1978`, Single Scattering Radiative Transfer :cite:`ulaby_microwave_2014`, :cite:`weis_evaluation_2020`) to retrieve different land surface and/or vegetation parameters.

Overview
--------
The different preprocessing steps are shown in :numref:`work_flow_step_1` and :numref:`work_flow_step_2`. Additionally, every processing step is explained in more detail in the following subsections. As it can be seen in :numref:`work_flow_step_1` and :numref:`work_flow_step_2` the preprocessing work-flow is split in two main parts. The preprocessing methods in :numref:`work_flow_step_1` can be applied separately for every image. Whereas the work-flow shown in :numref:`work_flow_step_2` needs several images which were preprocessed by the different steps presented in :numref:`work_flow_step_1`.

.. _work_flow_step_1:
.. figure:: images/work_flow_01.png
    :align: center
    :width: 80%

    Preprocessing chain showing processing steps to archive geometric and radiometric corrected Sentinel-1 data.

.. _work_flow_step_2:
.. figure:: images/work_flow_02.png
    :align: center
    :width: 80%

    Preprocessing chain showing processing steps to archive co-registered images which are multi-temporal speckle filtered

Sentinel-1 Level-1 SLC data
---------------------------
The preprocessing work-flow of :numref:`work_flow_step_1` is based on Sentinel-1 Level-1 SLC data. Among some other sources Sentinel-1 data can be downloaded from ESA's Copernicus Open Access Hub (`<https://scihub.copernicus.eu/>`_).


Sentinel-1 Level-1 SLC data are generated by the operational ESA Instrument Processing Facility (IPF). The SLC products are situated in slant range geometry. The slant range geometry is the natural radar one and is defined by the line-of-sight distance of the radar system to each reflecting object. The SLC product consists of focused SAR data in zero-Doppler orientation. Furthermore, for geo-referencing orbit and attitude information directly provided by the satellite are stored within the SLC product. Moreover the SAR data is corrected for errors caused by the well known azimuth bi-static delay, elevation antenna pattern and range spreading loss :cite:`sentinel-1_team_sentinel-1_2013`. In contrary to Level-1 Ground Range Detected (GRD) products SLC data preserve the real and imaginary part of the backscatter signal and contain therefore also the phase information :cite:`sentinel-1_team_sentinel-1_2013`. The IPF is generating SLC data for all available acquisition modes (StripMap (SM), Interferometric Wide (IW), Extra Wide (EW), and Wave (WV)) of the Sentinel-1 satellites. Further information about Sentinel-1 Level-1 products are gathered in ESA's Sentinel-1 User Handbook :cite:`sentinel-1_team_sentinel-1_2013` available at `<https://earth.esa.int/documents/247904/685163/Sentinel-1_User_Handbook>`_.


Precise orbit file
-------------------
Theory / Purpose
~~~~~~~~~~~~~~~~~~

During the acquisition of Sentinel-1 data the satellite position is recorded by a Global Navigation Satellite System (GNSS). To assure a fast delivery of Sentinel-1 products orbit information generated by an on-board navigation solution are stored within the Sentinel-1 Level-1 products. The orbit positions are later refined and made available as restituted or precise orbit files by the Copernicus Precise Orbit Determination (POD) Service. The POD products for Sentinel-1 data with given accuracy and availability after data acquisition are listed in :numref:`POD_table`.

.. _POD_table:
.. table:: Accuracy specification for Sentinel-1 POD products :cite:`sentinels_pod_team_sentinels_2016`
    :widths: auto

    +------------+--------------------------------------------+-------------+----------+
    |   Mission  | POD Product                                | Accuracy    | Latency  |
    +------------+--------------------------------------------+-------------+----------+
    |            | Restituted Orbit File                      | < 10 cm     | 3 hours  |
    |            +--------------------------------------------+-------------+----------+
    | Sentinel-1 | Precise Orbit Ephemerides (POE) Orbit file | < 5 cm      | 20 days  |
    |            +--------------------------------------------+-------------+----------+
    |            | Attitude Restituted Data                   | < 0.005 deg | 20 days  |
    +------------+--------------------------------------------+-------------+----------+

Precise orbit information can have a high influence on the quality of several preprocessing steps especially e.g. for the geo-referencing of the data. Therefore, it is always preferable to use the most accurate orbit information that is available.

Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
Since the preprocessing for the MULTIPLY project doesn't depend on near-real-time data the precise orbit file (available within 20 days) is used to update the orbit and velocity information within the Sentinel-1 SLC product. Therefore the operator "Apply Orbit Correction" of SNAP S1TBX toolbox is used.

Input:
    - Sentinel-1 SLC IW image (downloaded from Copernicus Open Access Hub)
    - Precise orbit file (automatic download by SNAP S1TBX)

Output:
    - Sentinel-1 SLC IW image with updated orbit information


Thermal noise removal
---------------------
Theory / Purpose
~~~~~~~~~~~~~~~~~~
Thermal noise is caused by the background energy of a SAR receiver and independent from the received signal power. Like some other noise factors thermal noise appears randomly over the entire image. But in contrary to quantization noise like speckle, which is connected to the signal power, thermal noise is hardly noticeable. Therefore, high impact of thermal noise on the quality of the data is especially given in areas like calm lakes, rivers and other with a low mean signal response detected by the SAR system. For the purpose of correction the IPF is calculating a thermal noise Look up Table (LUT) which is stored within the Sentinel-1 Level-1 product. More information about the calculation of the thermal noise for Sentinel-1 is given in :cite:`bourbigot_sentinel-1_2015`.


Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
The "Thermal Noise Removal" operator of SNAP S1TBX software is used to remove the thermal noise which is stored within a LUT within Sentinel-1 Level-1 products. Thermal noise removal can only applied on backscatter intensity therefore the phase information of the SLC data get lost.

Input:
    - Sentinel-1 SLC IW image with updated orbit information

Output:
    - Sentinel-1 SLC Intensity corrected by thermal noise

.. _radiometric_calibration:

Radiometric calibration
-------------------------
Theory / Purpose
~~~~~~~~~~~~~~~~~
Sentinel-1 Level-1 products are not radiometric corrected by default. However, for the quantitative use of SAR images a radiometric calibration of radar reflectivity (stored as Digital Numbers (DN) within Sentinel-1 Level-1 products) to physical units (radar backscatter) is essential. Otherwise a comparison of SAR images from different sensors or even the same sensor for different acquisition dates or different acquisition modes is not possible. To apply a radiometric calibration a Calibration Annotation Data Set (CADS) with four Look Up Tables (LUTs) are provided within the Sentinel-1 Level-1 products by Sentinel-1 Instrument Processing Facility (IPF). The four LUTs are used to convert DN to sigma naught, beta naught and gamma or vice versa. More information about the radiometric calibration is given in :cite:`miranda_radiometric_2015`.

Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
The "Radiometric Calibration" operator of SNAP S1TBX software is used to perform the conversion of DN to radar backscatter. In our case the output radar backscatter information is calibrated in Sigma naught.

Input:
    - Sentinel-1 SLC Intensity corrected by thermal noise

Output:
    - Sigma naught calibrated radar backscatter


TOPSAR Deburst
---------------
Theory / Purpose
~~~~~~~~~~~~~~~~~
Sentinel-1 Level-1 SLC images acquired in IW or EW swath mode consists of one image per swath and polarisation. IW products are made up of three swaths which means three images for single polarisation and six images for dual polarisation. EW products are made up of five swaths which means five images for single polarisation and ten images for dual polarisation. The sub-swath images consists of different bursts which are all processed as separate images. The different bursts are stored in one single image whereby each burst is separated by a black-filled demarcation :cite:`sentinel-1_team_sentinel-1_2013`. For the usage of Sentinel-1 Level-1 SLC data only one sub-swath can be extracted or several/all sub-swath can be combined to one image with fluent transitions between the sub-swaths. More detailed information are provided in :cite:`sentinel-1_team_sentinel-1_2013`, :cite:`daria_burst-mode_2007` and :cite:`de_zan_topsar_2006`.

Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
The "TOPSAR-Deburst" operator of SNAP S1TBX software is used to merge all sub-swath to retrieve one fluent image.

Input:
    - Sigma naught calibrated radar backscatter (with different sub-swath)

Output:
    - Sigma naught calibrated radar backscatter (with fluent transitions)


Geometric correction
---------------------
Theory / Purpose
~~~~~~~~~~~~~~~~~
An important part of the preprocessing chain is the geometric terrain correction. The geometric correction is a conversion of the Sentinel-1 SLC data from slant range geometry into a map coordinate system. Due to the acquisition geometry of the SAR different topographical distortions like foreshortening, layover or shadowing effects occur. The appropriate way to correct these distortions is the Range-Doppler approach. The method needs information about the topography (normally provided by a Digital Elevation Model (DEM)) as well as orbit and velocity information from the satellite (stored within Sentinel-1 SLC product) to correct the mentioned distortions and derive a precise geolocation for each pixel of the image.

Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
A geometric correction of the input data is performed by using the "Range Doppler Terrain Correction" method implement in SNAP's S1TBX software. Data from the Shuttle Radar Topography Mission (SRTM) with a resolution of 1-arc second (30 meters) is used for the necessary DEM.

Input:
    - Sigma naught calibrated radar backscatter (with fluent transitions)
    - SRTM data with 1-arc second resolution (automatic download by SNAP S1TBX)

Output:
    - Geometric corrected sigma naught calibrated radar backscatter (Map Projection WGS84)
    - Incidence angle from ellipsoid
    - Local incidence angle (based on SRTM)

Radiometric correction
---------------------------------------
Theory / Purpose
~~~~~~~~~~~~~~~~~
For the conversion of Sentinel-1 backscatter values to sigma or gamma naught, LUT's stored within the Sentinel-1 product are used (see :ref:`radiometric_calibration`). For the creation of the LUT's Sentinel-1 IPF is using an incidence angle of an ellipsoid inflated earth model :cite:`miranda_radiometric_2015`. Therefore, the local terrain variation within the image and their radiometric impact on the backscatter is considered insufficiently. A simple and widely used practice to consider the radiometric impact due to local terrain variations represents the approach to use the local incidence angle instead of the ellipsoid one :cite:`kellndorfer_toward_1998`. The radiometric corrected backscatter :math:`\sigma_{NORLIM}^{0}` used by kellndorfer_toward_1998 et al. :cite:`kellndorfer_toward_1998` can be calculated as

.. math::
    \sigma_{NORLIM}^{0} = \sigma_{Ell} \frac{sin \theta_{LIA}}{sin \theta_{Ell}}
    :label: kellndorfer_toward_1998

with :math:`\theta_{LIA}` as the local incidence angle, :math:`\theta_{Ell}` as the ellipsoid incidence angle used by IPF and the radar backscatter :math:`\sigma_{Ell}` calculated by using LUT's provided by IPF.

Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
Within the "Range Doppler Terrain Correction" method of SNAP's S1TBX software the radiometric normalisation approach of kellndorfer_toward_1998 et al. :cite:`kellndorfer_toward_1998` is implemented as a additional option. Unfortunately, the SNAP internal option can not be used with our kind of data. Therefore, normalisation after kellndorfer_toward_1998 et al :cite:`kellndorfer_toward_1998` is done by coding the equations within the "BandMath" operator of SNAP's S1TBX. The used local incidence angle is provided by the previous applied "Range Doppler Terrain Correction" and therefore the local incidence angle is based on the SRTM data.

Input:
    - Geometric corrected sigma naught calibrated radar backscatter (Map Projection WGS84)
    - Incidence angle from ellipsoid
    - Local incidence angle (based on SRTM)

Output:
    - Radiometric and geometric corrected sigma naught calibrated radar backscatter (Map Projection WGS84)

Backscatter normalisation
------------------------------------
Theory / Purpose
~~~~~~~~~~~~~~~~~
Beside the previously discussed geometric and radiometric distortions some other specific backscattering coefficient variations within the range direction of the image are caused by the image geometry of the SAR sensor. The backscattered energy of an illuminated area has not only a dependency on the area itself but also on the incidence angle. This means, backscatter values of a specific area with a small incidence angle return higher backscatter values then data of the same area acquired with a higher incidence angle. Incidence angle induced variations not only occur inside one image but also between images form different sensors as well as within one sensor through different acquisition geometries or different tracks or orbits. For a usage of Sentinel-1A and 1B time-series acquired with different orbits and/or different tracks and therefore most likly a high change between the incidence angles a backscatter normalisation is vital. A often and widely used technique to minimize backscatter variations caused by the incidence angle is the cosine correction :cite:`ulaby_microwave_1986`. The cosine correction is based on the Lambert's law for optics. Therefore, under the assumption that the backscattered energy in the upper hemisphere follows a cosine law and also the radiation variability has a cosine dependency, the received backscatter :math:`\sigma_{\theta_i}^{0}` and its dependency on the incidence angle can be written as

.. math::
    \sigma_{\theta_i}^{0} = \sigma_0^{0}cos^{n}(\theta_i)
    :label: cosine_1

with a weighting factor n and the incidence angle independent backscatter :math:`\sigma_{0}^{0}`.
With the cosine correction the backscatter of the Sentinel-1 products can therefore normalised to a reference angle :math:`\theta_{ref}` with

.. math::
    \sigma_{ref}^{0} = \frac{\sigma_{\theta_i}^{0}cos^{n}(\theta_{ref})}{cos^{n}_{\theta_i}}
    :label: cosine_2

Studies show that the weighting factor n is dependent on the roughness :cite:`ardila_angular_2010` and therefore the backscatter variations can vary with different land cover types. A schematic illustration of the backscatter variations considering the incidence angle is given in :numref:`wagner1999`.

.. _wagner1999:
.. figure:: images/wagner_1999.png
    :align: center
    :width: 60%

    Illustration of the backscatter variations considering the incidence angle dependency :cite:`wagner_study_1999`.


Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
The backscatter normlisation is applied by coding :eq:`cosine_2` in SNAP's S1TBX operator "BandMaths". As default a reference angle of 37,55° (average incidence angle for IW swath mode :cite:`bourbigot_sentinel-1_2015`) and a weighting factor of 2 (standard value) is specified. Through a configuration file the user can replace the default values for the reference angle and weighting factor to probably more suitable values of their specific applications.

Input:
    - Radiometric and geometric corrected sigma naught calibrated radar backscatter (Map Projection WGS84)
    - reference angle (default is 35°)
    - weighting factor (default is 2)

Output:
    - Radiometric and geometric corrected sigma naught calibrated radar backscatter values normalised to reference angle (Map Projection WGS84)


Co-registration
----------------
Theory / Purpose
~~~~~~~~~~~~~~~~~
For time-series analysis especially when applying a :ref:`multi_temporal_speckle_filter` the SAR image has to be co-registered. The co-registration is a method to get every image of the time-series on the same grid and also the pixel resolution.

Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
The co-registration as a requirement for the :ref:`multi_temporal_speckle_filter` is accomplished by the "Co-Registration" operator within SNAP's S1TBX. The "Co-Registration" operator in SNAP is defined as a completely automatic process. The operator consists of a stack creation (collocating master and slave image), a cross correlation (allignment between master ans slave image) and a warp (resamples pixels from the slave image to pixels of the master image).

Input:
    - Master image
    - Slave image(s)

Output:
    - Co-registered images

.. _multi_temporal_speckle_filter:

Multi-temporal speckle filter
-----------------------------
Theory / Purpose
~~~~~~~~~~~~~~~~~
A characteristic of images acquired by a SAR system is the visibility of random noise which look like "salt and pepper" within the image and is called speckle. The appearance of speckle is caused by the interferences of coherent echoes from individual scatterers within one pixel :cite:`woodhouse_introduction_2006`.The presence of speckle degrades the quality of the image and therefore it makes the interpretation of the SAR data more difficult. Over the years several approaches for speckle reduction were developed. They are mainly based on either multi-looking or filtering methods. Different filtering approaches like Frost, Lee etc. can be applied as a single or multi-temporal speckle filter. First findings with Sentinel-1 data show that a multi-temporal speckle filter provides better results in form of speckle reduction and resolution preservation then a single speckle filter. A major advantage for the usage of a multi-temporal speckle filter on Sentinel-1 data is the high temporal resolution availability. Nevertheless, more detailed studies on analysing the impact of different multi-temporal speckle filters on the retrieval of bio- and geophysical parameters from Sentinel-1 data are still lacking. Anyway, a usage of a multi-temporal filter significantly reduces the speckle and is therefore a essentinal part of our preprocessing chain.

Practical implementation
~~~~~~~~~~~~~~~~~~~~~~~~~
For the speckle reduction the "Multi-temporal Speckle Filter" operator within SNAP's S1TBX software is used. As default, 7 temporally consecutive images are used within the "Multi-temporal Speckle Filter" whereby the target image is temporally situated in the middle. The applied filter is a Lee filter with a spatial window size of 5x5 pixels, a sigma of 0.9, and a target window size of 3x3 pixels. The spatial averaging over pixel has a significant influence on spatial resolution information loss of the image. Therefore, the averaging pixel size might change during the project. If the image consists of two polarisations the filter is applied on each polarisation separately. The practical implementation in case of filter type, used polarisation, number of used images etc. may change with more experience of applying multi-temporal speckle filters and the occurring results.

Input:
    - x co-registered images (can be specified within configuration file)

Output:
    - speckle filtered images

Folder and data creation during pre-processing steps
----------------------------------------------------
Creation of
    - step1 (folder)
        - temporary results after applying processing steps shown in :numref:`work_flow_step_1`
    - step2 (folder)
        - co-registered images of step1
    - step3 (folder)
        - final results
    - foldername.nc (final netcdf stack file)
within specified output folder (config file)

Output layers of final netcdf stack file
-----------------------------------------
Output layer of default pre-processing chain
    - theta (local incidence angle)
    - sigma0_vv_single (single speckle filtered radiometric and geometric corrected sigma nought backscatter)
    - sigma0_vh_single (single speckle filtered radiometric and geometric corrected sigma nought backscatter)
    - sigma0_vv_multi (multi speckle filtered radiometric and geometric corrected sigma nought backscatter)
    - sigma0_vh_multi (multi speckle filtered radiometric and geometric corrected sigma nought backscatter)
    - sigma0_vv_norm_single (single speckle filtered radiometric and geometric corrected sigma nought backscatter normalized to a specific incidence angle)
    - sigma0_vh_norm_single (single speckle filtered radiometric and geometric corrected sigma nought backscatter normalized to a specific incidence angle)
    - sigma0_vv_norm_single (multi speckle filtered radiometric and geometric corrected sigma nought backscatter normalized to a specific incidence angle)
    - sigma0_vh_norm_single (multi speckle filtered radiometric and geometric corrected sigma nought backscatter normalized to a specific incidence angle)

Abbreviations and values within netcdf stack file
---------------------------------------------------------------------

Abbreviation within variable names
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    - theta = local incidence angle
    - sigma0 = radiometric and geometric corrected sigma nought backscatter
    - vv = VV polarization
    - vh = VH polarization
    - single = single speckle filter was applied
    - multi = multitemporal speckle filtered
    - norm = backscatter was normalized to a specific incidence angle


Values of specific variables
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
- orbitdirection
    - 0 = Ascending
    - 1 = Descending

- relorbit
    - number of relative orbit

- satellite
    - 0 = Sentinel-1 A
    - 1 = Sentinel-1 B

- name tags of processed data
    - theta (local incidence angle)
    - sigma0 (radiometric and geometric corrected sigma nought backscatter)
    - vv (polarization vv)
    - vh (polarization vh)
    - single (single speckle filtered)
    - multi (multi temporal speckle filtered)
    - norm (backscatter normalized to a specific incidence angle)

.. rubric:: References
.. bibliography:: references.bib
    :style: unsrt



.. _Technicaldocumantation:

Technical documentation
=======================

Sar-Pre-Processing
-------------------
.. automodule:: sar_pre_processing.sar_pre_processor
   :members:
   :undoc-members:
   :show-inheritance:

File list for Sar-Pre-Processing
---------------------------------
.. automodule:: sar_pre_processing.file_list_sar_pre_processing
   :members:
   :undoc-members:
   :show-inheritance:

netcdf-stack
---------------------------------
.. automodule:: sar_pre_processing.netcdf_stack
   :members:
   :undoc-members:
   :show-inheritance:

Attribute Dict
---------------------------------
.. automodule:: sar_pre_processing.attribute_dict
   :members:
   :undoc-members:
   :show-inheritance:
.. _Installation:

Installation for Linux (tested with Ubuntu 20.04)
==================================================
.. note::
    The SenSARP has been developed against Python 3.6.
    It cannot be guaranteed to work with previous Python versions.

The first step is to clone the latest code and step into the check out directory::

    git clone https://github.com/multiply-org/sar-pre-processing.git
    cd sar-pre-processing

Installation with Conda
------------------------
Download and install `Anaconda <https://www.anaconda.com/products/individual>`_ or `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. Anaconda/Miniconda installation instructions can be found `here <https://conda.io/projects/conda/en/latest/user-guide/install/linux.html#install-linux-silent>`_

To install all required modules, use::

    conda env create --prefix ./env --file environment.yml
    conda activate ./env # activate the environment

To install SenSARP into an existing Python environment, use::

    python setup.py install

To install for development, use::

    python setup.py develop

Installation with virtualenv and python
----------------------------------------
Install system requirements::

    sudo apt install python3-pip python3-tk python3-virtualenv python3-venv virtualenv

Create a virtual environment::

    virtualenv -p /usr/bin/python3 env
    source env/bin/activate # activate the environment
    pip install --upgrade pip setuptools # update pip and setuptools

To install SenSARP into an existing Python environment, use::

    python setup.py install

To install for development, use::

    python setup.py develop

GDAL package needs to be installed too::

    sudo apt install gdal-bin libgdal-dev

    python -m pip install pygdal=="`gdal-config --version`.*"


Further information
-------------------

.. literalinclude:: ./environment.yml

Please see the `environment file <https://github.com/multiply-org/sar-pre-processing/blob/master/environment.yml>`_ for a list of all installed dependencies during the installation process.
Additionally, ESA's SNAP Sentinel-1 Toolbox (Version >8.0.3) has to be installed prerequisite. The Software can be downloaded `here <http://step.esa.int/main/download/snap-download/>`_.
To install the SNAP toolbox, open a terminal window and use::

    bash esa-snap_sentinel_unix_8_0.sh

SenSARP uses only functionalities of the Sentinel-1 Toolbox.

.. note::
    Currently, only SNAP version 8.0 can be downloaded from the website.
    To update SNAP to a version >8.0.3 please start the SNAP software.
    You will be asked if you want to search for update.
    Please search for updates and install all updates.
    After the updates are installed, you need to restart SNAP to initialize the updates correctly.

SNAP Toolbox need libgfortran for specific operations but currently libgfortran is not installed during the installation process of SNAP, therefore you might use::

    sudo apt-get install gfortran
MULTIPLY - SenSARP
=============================

|buildstatus| |docstatus|

Content
---------

.. toctree::
    :maxdepth: 2

    01_Introduction
    02_Installation
    03_Usage
    04_Theory
    annex_technical
    Authors <authors>
    Changelog <changes>
    License

Credits
-------------
The project leading to this application has received funding from the
European Union’s Horizon 2020 research and innovation program
under grant agreement No 687320.

.. |buildstatus| image:: https://travis-ci.com/multiply-org/sar-pre-processing.svg?branch=getting_to_release
    :target: https://travis-ci.org/multiply-org/sar-pre-processing

.. |docstatus| image:: https://readthedocs.org/projects/multiply-sar-pre-processing/badge/?version=latest
    :target: https://multiply-sar-pre-processing.readthedocs.io/en/getting_to_release/?badge=getting_to_release
    :alt: Documentation Status
Usage
======

.. toctree::
    :maxdepth: 2

    ./notebooks/default_process_single_image.ipynb
    ./notebooks/default_process_time_series.ipynb
    ./notebooks/use_user_defined_graphs.ipynb
    031_config

