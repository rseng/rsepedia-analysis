<img src="https://github.com/NREL/OpenOA/blob/develop/Open%20OA%20Final%20Logos/Color/Open%20OA%20Color%20Transparent%20Background.png?raw=true" alt="OpenOA" width="300"/>

[![Binder Badge](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/NREL/OpenOA/main?filepath=examples) [![Gitter Badge](https://badges.gitter.im/NREL_OpenOA/community.svg)](https://gitter.im/NREL_OpenOA/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge) [![Journal of Open Source Software Badge](https://joss.theoj.org/papers/d635ef3c3784d49f6e81e07a0b35ff6b/status.svg)](https://joss.theoj.org/papers/d635ef3c3784d49f6e81e07a0b35ff6b)

[![Documentation Badge](https://readthedocs.org/projects/openoa/badge/?version=latest)](https://openoa.readthedocs.io) ![Tests Badge](https://github.com/NREL/OpenOA/workflows/Tests/badge.svg?branch=develop) [![Code Coverage Badge](https://codecov.io/gh/NREL/OpenOA/branch/develop/graph/badge.svg)](https://codecov.io/gh/NREL/OpenOA)

[![pre-commit](https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white)](https://github.com/pre-commit/pre-commit) [![Code style: black](https://img.shields.io/badge/code%20style-black-000000.svg)](https://github.com/psf/black) [![Imports: isort](https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336)](https://pycqa.github.io/isort/)

-----

This library provides a framework for working with large timeseries data from wind plants, such as SCADA.
Its development has been motivated by the WP3 Benchmarking (PRUF) project,
which aims to provide a reference implementation for plant-level performance assessment.

Analysis routines are grouped by purpose into methods,
and these methods in turn rely on more abstract toolkits.
In addition to the provided analysis methods,
anyone can write their own, which is intended to provide natural
growth of tools within this framework.

The library is written around Pandas Data Frames, utilizing a flexible backend
so that data loading, processing, and analysis could be performed using other libraries,
such as Dask and Spark, in the future.

If you would like to try out the code before installation or simply explore the possibilities, please see our examples on [Binder](https://mybinder.org/v2/gh/NREL/OpenOA/main?filepath=examples).

If you use this software in your work, please cite our JOSS article with the following BibTex:

```
@article{Perr-Sauer2021,
  doi = {10.21105/joss.02171},
  url = {https://doi.org/10.21105/joss.02171},
  year = {2021},
  publisher = {The Open Journal},
  volume = {6},
  number = {58},
  pages = {2171},
  author = {Jordan Perr-Sauer and Mike Optis and Jason M. Fields and Nicola Bodini and Joseph C.Y. Lee and Austin Todd and Eric Simley and Robert Hammond and Caleb Phillips and Monte Lunacek and Travis Kemper and Lindy Williams and Anna Craig and Nathan Agarwal and Shawn Sheng and John Meissner},
  title = {OpenOA: An Open-Source Codebase For Operational Analysis of Wind Farms},
  journal = {Journal of Open Source Software}
}
```

### Requirements

  * Python 3.6+ with pip.

We strongly recommend using the Anaconda Python distribution and creating a new conda environment for OpenOA. You can download Anaconda through [their website.](https://www.anaconda.com/products/individual)

After installing Anaconda, create and activate a new conda environment with the name "openoa-env":

```
conda create --name openoa-env python=3.8
conda activate openoa-env
```

### Installation

Clone the repository and install the library and its dependencies using pip:

```
git clone https://github.com/NREL/OpenOA.git
pip install ./OpenOA
```

You should now be able to import operational_analysis from the Python interpreter:

```
python
>>> import operational_analysis
```

#### Common Installation Issues:

- In Windows you may get an error regarding geos_c.dll. To fix this install Shapely using:

```
conda install Shapely
```

- In Windows, an ImportError regarding win32api can also occur. This can be resolved by fixing the version of pywin32 as follows:

```
pip install --upgrade pywin32==255
```

### Development

Development dependencies are provided through the develop extra flag in setup.py. Here, we install OpenOA, with development dependencies, in editable mode, and activate the pre-commit workflow (note: this second step must be done before committing any
changes):

```
pip install -e "./OpenOA[develop]"
pre-commit install
```

Occasionally, you will need to update the dependencies in the pre-commit workflow, which will provide an error when this needs to happen. When it does, this can normally be resolved with the below code, after which you can continue with your normal git workflow:
```
pre-commit autoupdate
git add .pre-commit-config.yaml
```

#### Example Notebooks and Data

The example data will be automaticaly extracted as needed by the tests. To manually extract the example data for use with the example notebooks, use the following command:

```
unzip examples/data/la_haute_borne.zip -d examples/data/la_haute_borne/
```

In addition, you will need to install the packages required for running the examples with the following command:

```
pip install -r ./OpenOA/examples/requirements.txt
```

The example notebooks are located in the `examples` directory. We suggest installing the Jupyter notebook server to run the notebooks interactively. The notebooks can also be viewed statically on [Read The Docs](http://openoa.readthedocs.io/).

```
jupyter notebook
```

#### Testing
Tests are written in the Python unittest framework and are runnable using pytest. There are two types of tests, unit tests (located in `test/unit`) run quickly and are automatically for every pull request to the OpenOA repository. Regression tests (located at `test/regression`) provide a comprehensive suite of scientific tests that may take a long time to run (up to 20 minutes on our machines). These tests should be run locally before submitting a pull request, and are run weekly on the develop and main branches.

To run all unit and regresison tests:
```
pytest
```

To run unit tests only:
```
pytest test/unit
```

To run all tests and generate a code coverage report
```
pytest --cov=operational_analysis
```

#### Documentation

Documentation is automatically built by, and visible through, [Read The Docs](http://openoa.readthedocs.io/).

You can build the documentation with [sphinx](http://www.sphinx-doc.org/en/stable/), but will need to ensure [Pandoc is installed](https://pandoc.org/installing.html) on your computer first:

```
cd sphinx
pip install -r requirements.txt
make html
```


### Contributors

Alphabetically:
Nathan Agarwal,
Nicola Bodini,
Anna Craig,
Jason Fields,
Rob Hammond,
Travis Kemper,
Joseph Lee,
Monte Lunacek,
John Meissner,
Mike Optis,
Jordan Perr-Sauer,
Sebastian Pfaffel,
Caleb Phillips,
Charlie Plumley,
Eliot Quon,
Sheungwen Sheng,
Eric Simley, and
Lindy Williams.
# Changelog
All notable changes to this project will be documented in this file. If you make a notable change to the project, please add a line describing the change to the "unreleased" section. The maintainers will make an effort to keep the [Github Releases](https://github.com/NREL/OpenOA/releases) page up to date with this changelog. The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/).

## [2.3 - 2022-01-18]
- Replaced hard-coded reanalysis dates in plant analysis with automatic valid date selection and added optional user-defined end date argument. Fixed bug in normalization to 30-day months.
- Toolkit added for downloading reanalysis data using the PlanetOS API
- Added hourly resolution to AEP calculation
- Added wind farm plotting function to pandas_plotting toolkit using the Bokeh library
- Split the QC methods into a more generic `WindToolKitQualityControlDiagnosticSuite` class and WTK-specific subclass: `WindToolKitQualityControlDiagnosticSuite`.
- Updated filter algorithms in AEP calculation, now with a proper outlier filter

## [2.2 - 2021-05-28]
- IAV incorporation in AEP calculation
- Set power to 0 for windspeeds above and below cutoff in IEC power curve function.
- Split unit tests from regression tests and updated CI pipeline to run the full regression tests weekly.
- Flake8 with Black code style implemented with git hook to run on commit
- Updated long-term loss calculations to weight by monthly/daily long-term gross energy
- Added wind turbine asset data to example ENGIE project
- Reduce amount of time it takes to run regression tests by decreasing number of monte carlo iterations. Reduce tolerance of float comparisons in plant analysis regression test. Linear regression on daily data is removed from test.
- Bugfixes, such as fixing an improper python version specifier in setup.py and replacing some straggling references to the master branch with main.

## [2.1 - 2021-02-17]
- Modify bootstrapping approach for period of record sampling. Data is now sampled with replacement, across 100% of the POR data.
- Cleaned up dependencies for JOSS review. Adding peer-reviewed JOSS paper.
- Add Binder button to Readme which makes running the example notebooks easier.
- Set maximum python version to 3.8, due to an install issue for dependency Shapely on Mac with Python 3.9.

## [2.0.1 - 2020-10-13]
- Replaced `GeoPandas` functionality with `pyproj` and `Shapely` for coordinate
reference system conversion and distance measurements.
- Moved and renamed tests and updated the documentation accordingly.

## [2.0.0 - 2020-08-11]
- Switch to [semantic versioning](https://semver.org) from this release forward.
- Efficiency improvements in AEP calculation
- Energy Yield Analysis (EYA) added to Operational Assessment (OA) Gap Analysis method
- Uncertainty quantification for electrical losses and longterm turbine gross energy
- Implemented open source Engie example data
- Complete update of example notebooks
- Switch to standard BSD-3 Clause license
- Automated quality control method to assist with data ingestion. Tools in this method include daylight savings time change detection and identification of the diurnal cycle.
- Add electrical losses method
- Method for estimating long-term turbine gross energy (excluding downtime and underperformance losses)
- CI pipeline using Github Actions includes regression testing with Pytest, code coverage reporting via CodeCov, packaging and distribution via Pypi, and automatic documentation using ReadTheDocs.

## [1.1] - 2019-01-29
- Python3 Support
- Addition of reanalysis schemas to the Sphinx documentation
- Easy import of EIA data using new module: Metadata_Fetch
- Updated contributing.md document
- Quality checks for reanalysis data
- Improved installation instructions
- Integration tests are now performed in CI
- Performed PEP8 linting

## [1.0] - 2018-12-06
- Refactor many analysis and toolkit modules to make them conform to a standard API (init, prepare, and run method).
- Timeseries Table is now an integrated component, no sparkplug-datastructures dependency
- Plant Level AEP method w/ Monte Carlo
- Turbine / Scada level toolkits: Filtering, Imputing, Met, Pandas Plotting, Timeseries, Unit Conversion
- Most toolkits and all methods are fully documented in Sphinx.
- Two example notebooks: Operational AEP Analysis and Turbine Analysis
- All toolkits except for Pandas Plotting have > 80% test coverage.

Contributing
============

## Issue Tracking

New feature requests, changes, enhancements, non-methodology features, and bug reports can be filed as new issues in the
[Github.com issue tracker](https://github.com/NREL/OpenOA/issues) at any time. Please be sure to fully describe the
issue.

For other issues, please email the OpenOA distribution list at `openoa@nrel.gov`.

## Repository

The OpenOA repository is hosted on Github, and located here: http://github.com/NREL/OpenOA

This repository is organized using a modified git-flow system. Branches are organized as follows:

- main: Stable release version. Must have good test coverage and may not have all the newest features.
- develop: Development branch which contains the newest features. Tests must pass, but code may be unstable.
- feature/xxx: Branch from develop, should reference a github issue number.

To work on a feature, please fork OpenOA first and then create a feature branch in your own fork.
Work out of this feature branch before submitting a pull request.
Be sure to periodically synchronize the upstream develop branch into your feature branch to avoid conflicts in the pull request.

When the feature branch is ready, make a pull request to NREL/OpenOA through the Github.com UI. When submitting a pull request, you will need to accept the Contributor License Agreement(CLA). [CLA Language](https://gist.github.com/Dynorat/118aaa0c8277be986c59c32029898faa)

[![CLA assistant](https://cla-assistant.io/readme/badge/NREL/OpenOA)](https://cla-assistant.io/NREL/OpenOA)


## Pull Request

Pull requests must be made for all changes.
Most pull requests should be made against the develop branch.
Only core developers should make pull requests to the main branch.
Pull requests must include updated documentation and pass all unit tests and integration tests.
In addition, code coverage should not be negatively affected by the pull request.

**Scope:** Encapsulate the changes of ideally one, or potentially a couple, issues.
It is greatly preferable to submit three small pull requests than it is to submit one large pull request.
Write a complete description of these changes in the pull request body.

**Tests:** Must pass all tests. Pull requests will be rejected if tests do not pass.
Tests are automatically run through Github Actions for any pull request or push to the main or develop branches.

**Documentation:** Include any relevant changes to inline documentation, as well as any changes to the RST files
located in /sphinx.

**Coverage:** The testing framework (described below) will generate a coverage report. Please ensure that your
work is fully covered.

**Changelog:** For pull requests that encapsulate a user-facing feature, or is significant to users of OpenOA for some other reason, please add a line to CHANGELOG.md in the [Unreleased] section.

## Coding Style

This code uses a ``pre-commit`` workflow where code styling and linting is taken care of when a user
commits their code. Specifically, this code utilizes ``black`` for automatic formatting (line length, quotation usage, hanging
lines, etc.), ``isort`` for automatic import sorting, and ``flake8`` for linting.

To activate the ``pre-commit`` workflow, the user must install the develop version as outlined in the
[Readme](https://github.com/NREL/OpenOA/tree/develop#Development), and run the following line:

```
pre-commit install
```

## Documentation Style

Documentation is written using RST, and is located both inline and within the /sphinx directory.
Any changes to the analysis methodology should be discussed there or offline. Once a methodology change is decided,
create new tickets in this repository towards implementing the change.

## Testing

All code should be paired with a corresponding unit or integration test.
OpenOA uses pytest and the built in unittest framework.
For instructions on running tests, please see the [Readme](https://github.com/NREL/OpenOA/tree/develop#Testing).

## Release Process

 - Bump version number and metadata in
   - operational_analysis/__init__.py
   - operational_analysis/setup.py
   - sphinx/config.py
 - Bump version numbers of any dependencies in
   - setup.py
   - requirements.txt
   - sphinx/requirements.txt
 - Update the changelog, removing the UNRELEASED section and converting it into a release heading.
 - Make a pull request into develop with these updates
   - Note: Ensure all tests pass and the documentation is building correctly prior to merging
 - Merge develop into main through the git command line
   ```
    git checkout main
    git merge develop
    git push
    ```
 - Tag the new release version:
   ```
    git tag -a v1.2.3 -m "Tag messgae for v1.2.3"
    git push origin v1.2.3
    ```
 - Deploying a Package to PyPi
    - The repository is equipped with a github action to build and publish new versions to PyPi. A maintainer can invoke this workflow by pushing a tag to the NREL/OpenOA reposiory with prefix "v", such as "v1.1.0".
    - The action is defined in `.github/workflows/tags-to-pypi.yml`.
<img src="https://github.com/NREL/OpenOA/blob/develop/Open%20OA%20Final%20Logos/Color/Open%20OA%20Color%20Transparent%20Background.png?raw=true" alt="OpenOA" width="300"/>

Project Goals
=============
OpenOA is an open source software project started at the
National Renewable Energy Laboratory (NREL) to standardize,
streamline, and improve data analysis workflows in the wind
energy industry.

### Impact

- Provide a reference implementation
- Develop and improve standards for
    - Data
    - Methods
    - Metrics
    - Status Codes
    - Failures

### Capabilities

- Performance
    - AEP
    - Benchmarking
        - Preconstruction EYA
        - Operations
        - Performance Changes
        - Gap Analysis
    - Evaluate / Diagnosis
        - Sub-optimal Operation
        - Performance Improvements
        - Plant vs Turbine
        - Anomaly Detection
    - Lost Energy Estimates
    - Portfolio Analysis
- Reliability
    - Fault Prediction
    - Anomaly detection

### Software Development

- Data Standardization
    - Digital Exchange File Format
    - OpenOA Data Model
- ML Methods
- Uncertainty Quantification
- Scalability and Speed
    - Parallel + Distributed computing
- Portability
- Packaging and ease of software installation
- Unit testing
- Reproducibility
---
title: 'OpenOA: An Open-Source Codebase For Operational Analysis of Wind Farms'
tags:
  - Python
  - wind energy
  - operational analysis
  - data analysis
  - standardization
authors:
  - name: Jordan Perr-Sauer
    orcid: 0000-0003-1571-1887
    affiliation: 1 # (Multiple affiliations must be quoted)
  - name: Mike Optis
    orcid: 0000-0001-5617-6134
    affiliation: 1
  - name: Jason M. Fields
    affiliation: 1
    orcid: 0000-0002-8781-6138
  - name: Nicola Bodini
    orcid: 0000-0002-2550-9853
    affiliation: 1
  - name: Joseph C.Y. Lee
    orcid: 0000-0003-1897-6290
    affiliation: 1
  - name: Austin Todd
    orcid: 0000-0002-1123-0982
    affiliation: 1
  - name: Eric Simley
    orcid: 0000-0002-1027-9848
    affiliation: 1
  - name: Robert Hammond
    orcid: 0000-0003-4476-6406
    affiliation: 1
  - name: Caleb Phillips
    affiliation: 1
  - name: Monte Lunacek
    orcid: 0000-0003-3755-224X
    affiliation: 1
  - name: Travis Kemper
    affiliation: 1
  - name: Lindy Williams
    affiliation: 1
  - name: Anna Craig
    affiliation: 1
  - name: Nathan Agarwal
    orcid: 0000-0002-2734-5514
    affiliation: 1
  - name: Shawn Sheng
    orcid: 0000-0003-0134-0907
    affiliation: 1
  - name: John Meissner
    affiliation: 1
affiliations:
 - name: National Renewable Energy Laboratory, Golden, CO, USA
   index: 1
date: 20 December 2019
bibliography: paper.bib
---

<!--
JOSS welcomes submissions from broadly diverse research areas. For this reason, we require that authors include in the paper some sentences that explain the software functionality and domain of use to a non-specialist reader. We also require that authors explain the research applications of the software. The paper should be between 250-1000 words.

Your paper should include:

A list of the authors of the software and their affiliations, using the correct format (see the example below).
A summary describing the high-level functionality and purpose of the software for a diverse, non-specialist audience.
A clear Statement of Need that illustrates the research purpose of the software.
A list of key references, including to other software addressing related needs.
Mention (if applicable) a representative set of past or ongoing research projects using the software and recent scholarly publications enabled by it.
Acknowledgement of any financial support.
As this short list shows, JOSS papers are only expected to contain a limited set of metadata (see example below), a Statement of Need, Summary, Acknowledgements, and References sections. You can look at an example accepted paper. Given this format, a “full length” paper is not permitted, and software documentation such as API (Application Programming Interface) functionality should not be in the paper and instead should be outlined in the software documentation.


Review Checklist:
Summary: Has a clear description of the high-level functionality and purpose of the software for a diverse, non-specialist audience been provided?
A statement of need: Do the authors clearly state what problems the software is designed to solve and who the target audience is?
State of the field: Do the authors describe how this software compares to other commonly-used packages?
Quality of writing: Is the paper well written (i.e., it does not require editing for structure, language, or writing quality)?
References: Is the list of references complete, and is everything cited appropriately that should be cited (e.g., papers, datasets, software)? Do references in the text use the proper citation syntax?

-->

# Summary

OpenOA is an open source framework for operational data analysis of wind energy plants, implemented in the Python programming language [@Rossum2009]. OpenOA provides a common data model, high level analysis workflows, and low-level convenience functions that engineers, analysts, and researchers in the wind energy industry can use to facilitate analytics workflows on operational data sets. OpenOA contains documentation, worked out examples in Jupyter notebooks, and a corresponding example dataset from the [Engie Renewable's La Haute Borne Dataset](https://opendata-renewables.engie.com/explore/dataset/d543716b-368d-4c53-8fb1-55addbe8d3ad/information).

Originally released to the public in 2018 [@osti_1478526], OpenOA is now actively developed through a public Github repository. With over 50 stars on Github, a dozen contributors, and an active issues forum, OpenOA is becoming a mature project that provides a high level interface to solve practical problems in the wind energy industry. OpenOA V2 is released as a Python package and is freely available under a business-friendly, open-source, BSD license. By committing to open source development, OpenOA hopes to facilitate reproducibility of research in this field, provide benchmark implementations of commonly performed data transformation and analysis tasks, and to serve as a conduit that delivers state-of-the-art analysis methods from researchers to practitioners.

Most users will interface with OpenOA through its analysis `methods` module. This includes Python classes which conform to a common interface (e.g., they implement `__init__`, `prepare`, and `run` methods). Version 2 of OpenOA implements three high level analysis methods for the calculation of: (1) Long term corrected annual energy production (AEP), (2) electrical losses, and (3) turbine level losses. Uncertainty quantification is achieved in each analysis using a Monte Carlo approach. A more detailed description of these analyses are provided [in the documentation](https://openoa.readthedocs.io). Low level functions that operate on Pandas series objects are organized in the `toolkit` module. These Python functions are written to be as generic as possible, and can be applied across multiple domains.

The OpenOA data model is implemented in the `types` module using a class called `PlantData`, which contains at least one Pandas data frame [@Mckinney2010]. These classes add convenience functions and a domain-specific schema based on the IEC 6400-25 standard [@iec25]. OpenOA is part of the ENTR alliance consortium, which envisions a complete software stack centered around an open source implementation of this standard.

OpenOA depends on scikit-learn [@Pedregosa2011] and numpy [@oliphant2006guide], with graphing functions implemented using matplotlib [@hunter2007matplotlib]. The OpenOA development team strives to use modern software development practices. Documentation is compiled from the source code and automatically published to [ReadTheDocs](https://openoa.readthedocs.io). We use Github actions to implement our continuous integration pipeline, including automated unit and regression tests, test coverage reporting via CodeCov, automated packaging and publication to the PyPI package index. We utilize a modified git-flow development workflow, with pull requests and issue tracking on Github driving the development.



# Statement of Need

OpenOA was created and is primarily developed by researchers at the National Renewable Energy Laboratory (NREL) through the Performance, Risk, Uncertainty, and Finance (PRUF) project. The PRUF team recognized the need to compute a long term corrected AEP (comparable to a 20-year estimate) from operational data as part of an industry-wide benchmarking study [@lunacek2018]. Due to access restrictions on the input data, open source publication of the code was necessary to foster trust in the results of the benchmarking study. Furthermore, after talking with our industry partners, it became clear that there was no industry standard method for computing a long term corrected AEP. Currently, participants in the wind energy industry who wish to compute metrics like AEP must rely on commercial secondary supervisory control and data acquisition (SCADA) software providers, or must develop their own methodologies internally. We know of no other open source software package to compute long term corrected AEP.

![Figure](openoa-joss-figure.png)
*Figure 1: A subset of graphical outputs from the OpenOA documentation. Clockwise from the top, (A) power curve with extreme values highlighted in red, (B) distribution of long term corrected AEP, (C) time-series of wind speed from multiple reanalysis products showing anomalously low wind speed for a highlighted period of record.*

Operational analysis involves obtaining time-series data from an industrial plant's SCADA system, performing quality control processes on these data, and computing metrics that provide insight into the performance charactertistics of a wind plant. Figure 1 contains some graphical outputs that are generated by OpenOA. Since its inception, OpenOA has been used in several published studies at NREL. An early version of the code was used in @Craig2018 to quantify the uncertainty in AEP estimates resulting from analyst choices. In @Bodini2020, it is used to calculate long-term operational AEP from over 470 wind plants in the US to assess correlation between uncertainty components. OpenOA will also be used in an upcoming technical report for the PRUF project's industry-wide benchmarking study.

# Acknowledgements
The authors would like to acknowledge that Jordan Perr-Sauer and Mike Optis have made an equal contribution to this work.
This work was authored by the National Renewable Energy Laboratory, operated by Alliance for Sustainable Energy, LLC, for the U.S. Department of Energy (DOE) under Contract No. DE-AC36-08GO28308.
Funding provided by the U.S. Department of Energy Office of Energy Efficiency and Renewable Energy Wind Energy Technologies Office, within the Atmosphere to Electrons research program.
The views expressed in the article do not necessarily represent the views of the DOE or the U.S. Government.
The U.S. Government retains and the publisher, by accepting the article for publication, acknowledges that the U.S. Government retains a nonexclusive, paid-up, irrevocable, worldwide license to publish or reproduce the published form of this work, or allow others to do so, for U.S. Government purposes.

# References
Please be sure to install the required packages for running the example notebooks by
executing the following command:

```
pip install -r requirements.txt
```
License Information for Example Data
=============
The example data provided in **la_haute_borne.zip** are based on data for the "La Haute Borne" wind farm published by ENGIE under the [Open License version 2.0](https://www.etalab.gouv.fr/wp-content/uploads/2017/04/ETALAB-Licence-Ouverte-v2.0.pdf) (English version of license [here](https://www.etalab.gouv.fr/wp-content/uploads/2018/11/open-licence.pdf)).

ENGIE - Original data downloaded from https://opendata-renewables.engie.com/explore/dataset/d543716b-368d-4c53-8fb1-55addbe8d3ad/information, updated on 9 October 2019.

More information about the original data is available at https://opendata-renewables.engie.com
.. _contributing:


.. mdinclude:: ../contributing.md
.. methods:

Analysis Methods
======================

Analysis Methods are a collection of modules and classes to accomplish high level analysis of plant data. They integrate
the more abstract toolkit modules with specific PlantData classes to create a reproducible process. Analysis modules may
respect the engine parameter of PlantData objects to support multi-platform, scalable analysis.

Plant Level Analysis
--------------------

.. automodule:: operational_analysis.methods.plant_analysis
    :members:
    :undoc-members:
    :show-inheritance:

Turbine Level Analysis
----------------------

.. automodule:: operational_analysis.methods.turbine_long_term_gross_energy
    :members:
    :undoc-members:
    :show-inheritance:

Electrical Losses Analysis
--------------------------

.. automodule:: operational_analysis.methods.electrical_losses
    :members:
    :undoc-members:
    :show-inheritance:
.. types:

Project Data
==================

New data imported into the operational_analysis toolkit can take advantage of the data structures in the `types` module.
The base class is operational_analysis.types.PlantData, which represents all known pieces of data about a given wind plant.


Schemas
-------

Operational Data
^^^^^^^^^^^^^^^^

PlantData.scada

==================== ====================
 Field Name           Data Type
==================== ====================
 time                 datetime64[ns]
 id                   string
 power_kw             float64
 windspeed_ms         float64
 winddirection_deg    float64
 status_label         string
 pitch_deg            float64
 temp_c               float64
==================== ====================

PlantData.meter

==================== ====================
 Field Name           Data Type
==================== ====================
 time                 datetime64[ns]
 power_kw             float64
 energy_kw            float64
==================== ====================


PlantData.tower

==================== ====================
 Field Name           Data Type
==================== ====================
 time                 datetime64[ns]
 id                   float64
==================== ====================

PlantData.curtail

==================== ====================
 Field Name           Data Type
==================== ====================
 time                 datetime64[ns]
 curtailment_pct      float64
 availability_pct     float64
 net_energy           float64
==================== ====================


PlantData.status

==================== ====================
 Field Name           Data Type
==================== ====================
 time                 datetime64[ns]
 id                   string
 status_id            int64
 status_code          int64
 status_text          string
==================== ====================


PlantData.asset

==================== ====================
 Field Name           Data Type
==================== ====================
 id                   string
 latitude             float64
 longitude            float64
 rated_power_kw       float64
 type                 string
==================== ====================


Reanalysis Products
^^^^^^^^^^^^^^^^^^^

Reanalysis products are included as Plant Data objects and, regardless
of data source, have a standardized set of field names and types (see below).
That said, the data sources are obviously different, as are the methods use to
calculate these standard fields from the raw datasets. These methods are described here.

PlantData.reanalysis.product["merra2"]

MERRA-2 data are based on the single-level diagnostic data available here:

https://disc.gsfc.nasa.gov/datasets/M2T1NXSLV_V5.12.4/summary?keywords=%22MERRA-2%22

Wind speed and direction are taken directly from the diagnostic 50-m u- and v-wind
fields provided in this dataset. Air density at 50m is calculated using temperature
and pressure estimations at 50m and the ideal gas law. Temperature at 50m is estimated by taking the 10-m
temperature data provided by this dataset and assuming a constant lapse rate of -9.8
degrees Celsius per vertical kilometer. Pressure at 50m is extrapolated from surface pressure
data provided in this dataset using the hypsometric equation.

==================== ====================
 Field Name           Data Type
==================== ====================
 time                 datetime64[ns]
 windspeed_ms         float64
 winddirection_deg    float64
 rho_kgm-3            float64
==================== ====================

PlantData.reanalysis.product["ncep2"]

NCEP-2 data are based on the single-level diagnostic data available here:

https://rda.ucar.edu/datasets/ds091.0/

Wind speed and direction are taken directly from the diagnostic 10-m u- and v-wind
fields provided in this dataset. Air density at 10m is calculated using temperature
and pressure estimations at 10m and the ideal gas law. Temperature at 10m is estimated by taking the 2-m
temperature data provided by this dataset and assuming a constant lapse rate of -9.8
degrees Celsius per vertical kilometer. Pressure at 10m is extrapolated from surface pressure
data provided in this dataset using the hypsometric equation.

==================== ====================
 Field Name           Data Type
==================== ====================
 time                 datetime64[ns]
 windspeed_ms         float64
 winddirection_deg    float64
 rho_kgm-3            float64
==================== ====================


PlantData.reanalysis.product["erai"]

ERA-interim data are based on the model-level data available here:

https://rda.ucar.edu/datasets/ds627.0/

Model levels are based on sigma coordinates (i.e. fractions of surface pressure). From this dataset, we
extract temperature, u-wind, and v-wind at the 58th model level, which is on average about 72m above ground level
(https://www.ecmwf.int/en/forecasts/documentation-and-support/60-model-levels). We also extract surface pressure
data. Air density at the 58th model level is calculated using temperature data extracted at that level and an estimation
of pressure at that level using the ideal gas law. Pressure at the 58th model level is extrapolated from surface pressure
data provided in this dataset using the hypsometric equation.

==================== ====================
 Field Name           Data Type
==================== ====================
 time                 datetime64[ns]
 windspeed_ms         float64
 winddirection_deg    float64
 rho_kgm-3            float64
==================== ====================



PlantData
----------------------------------------

.. automodule:: operational_analysis.types.plant
    :members:
    :no-undoc-members:
    :show-inheritance:


AssetData
----------------------------------------

.. automodule:: operational_analysis.types.asset
    :members:
    :no-undoc-members:
    :show-inheritance:


ReanalysisData
----------------------------------------

.. automodule:: operational_analysis.types.reanalysis
    :members:
    :no-undoc-members:
    :show-inheritance:
.. _credit:


Credit
******

Alphabetically:
Nathan Agarwal,
Nicola Bodini,
Anna Craig,
Jason Fields,
Rob Hammond,
Travis Kemper,
Joseph Lee,
Monte Lunacek,
John Meissner,
Mike Optis,
Jordan Perr-Sauer,
Sebastian Pfaffel,
Caleb Phillips,
Charlie Plumley,
Eliot Quon,
Sheungwen Sheng,
Eric Simley, and
Lindy Williams.
.. highlight:: rst


OpenOA Operational Analysis Framework
#####################################

|Binder Badge| |Gitter Badge| |Journal of Open Source Software Badge|

|Documentation Badge| |Tests Badge| |Code Coverage Badge|

|pre-commit| |Code style: black| |Imports: isort|

.. ::

    # with overline, for parts
    * with overline, for chapters
    =, for sections
    -, for subsections
    ^, for subsubsections
    ", for paragraphs

This library provides a generic framework for working with large timeseries data from wind plants. Its development
has been motivated by the WP3 Benchmarking (PRUF) project, which aims to provide a reference implementaiton for
plant-level performance assessment.

The implementation makes use of a flexible backend, so that data loading, processing, and analysis can be performed
locally (e.g., with Pandas dataframes), in a semi-distributed manner (e.g., with Dask dataframes), or in a fully
distributed matter (e.g., with Spark dataframes).

Data processing and ETL is handled by the PlantData class and by project-specific modules which implement subclasses.
These modules can be used to import, inspect, pre-process, and save the raw data from wind turbines, meters, met towers,
and reanalysis products such as Merra2.

Analysis routines are grouped by purpose into toolkits - which provide an abstract low level API for common
computations, and methods - which provide higher level wind industry specific API. In addition to these provided modules,
anyone can write their own, which is intended to provide natural growth of tools within this framework.

To interact with how each of these components of OpenOA are used, please visit our examples notebooks on
`Binder <https://mybinder.org/v2/gh/NREL/OpenOA/master?filepath=examples>`_, or view them statically on the
`examples page <examplesout.ipynb>`_.

If you use this software in your work, please cite our JOSS article with the following BibTex::

    @article{Perr-Sauer2021,
        doi = {10.21105/joss.02171},
        url = {https://doi.org/10.21105/joss.02171},
        year = {2021},
        publisher = {The Open Journal},
        volume = {6},
        number = {58},
        pages = {2171},
        author = {Jordan Perr-Sauer and Mike Optis and Jason M. Fields and Nicola Bodini and Joseph C.Y. Lee and Austin Todd and Eric Simley and Robert Hammond and Caleb Phillips and Monte Lunacek and Travis Kemper and Lindy Williams and Anna Craig and Nathan Agarwal and Shawn Sheng and John Meissner},
        title = {OpenOA: An Open-Source Codebase For Operational Analysis of Wind Farms},
        journal = {Journal of Open Source Software}
    }


.. toctree::
    :maxdepth: 2
    :caption: Contents:

    install.rst
    examplesout
    toolkits.rst
    methods.rst
    types.rst
    contributing.rst
    credit.rst


.. |Binder Badge| image:: https://mybinder.org/badge_logo.svg
   :target: https://mybinder.org/v2/gh/NREL/OpenOA/main?filepath=examples
.. |Gitter Badge| image:: https://badges.gitter.im/NREL_OpenOA/community.svg
   :target: https://gitter.im/NREL_OpenOA/community?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge
.. |Journal of Open Source Software Badge| image:: https://joss.theoj.org/papers/d635ef3c3784d49f6e81e07a0b35ff6b/status.svg
   :target: https://joss.theoj.org/papers/d635ef3c3784d49f6e81e07a0b35ff6b
.. |Documentation Badge| image:: https://readthedocs.org/projects/openoa/badge/?version=latest
   :target: https://openoa.readthedocs.io
.. |Tests Badge| image:: https://github.com/NREL/OpenOA/workflows/Tests/badge.svg?branch=develop
.. |Code Coverage Badge| image:: https://codecov.io/gh/NREL/OpenOA/branch/develop/graph/badge.svg
   :target: https://codecov.io/gh/NREL/OpenOA
.. |pre-commit| image:: https://img.shields.io/badge/pre--commit-enabled-brightgreen?logo=pre-commit&logoColor=white
   :target: https://github.com/pre-commit/pre-commit
.. |Code style: black| image:: https://img.shields.io/badge/code%20style-black-000000.svg
   :target: https://github.com/psf/black
.. |Imports: isort| image:: https://img.shields.io/badge/%20imports-isort-%231674b1?style=flat&labelColor=ef8336
   :target: https://pycqa.github.io/isort/
.. toolkits:

Toolkits
=============================

Toolkit modules provide abstract functions that operate on Pandas data frames and series. They can be used in isolation
as a useful library, or can be used through a higher level Analysis method. Toolkit modules are organized by function
and in general will only operate on data types from one particular backend. Currently, every toolkit function is
implemented with Pandas.

Filters
---------------------

.. automodule:: operational_analysis.toolkits.filters
    :members:
    :undoc-members:
    :show-inheritance:


Power Curve
---------------------

.. automodule:: operational_analysis.toolkits.power_curve
    :undoc-members:

.. automodule:: operational_analysis.toolkits.power_curve.functions
    :members:
    :undoc-members:
    :show-inheritance:


Imputing
---------------------

.. automodule:: operational_analysis.toolkits.imputing
    :members:
    :undoc-members:
    :show-inheritance:

Timeseries
---------------------

.. automodule:: operational_analysis.toolkits.timeseries
    :members:
    :undoc-members:
    :show-inheritance:

Met Data Processing
---------------------

.. automodule:: operational_analysis.toolkits.met_data_processing
    :members:
    :undoc-members:
    :show-inheritance:


Metadata Fetch
---------------------

.. automodule:: operational_analysis.toolkits.metadata_fetch
    :members:
    :undoc-members:
    :show-inheritance:


Unit Conversion
---------------------

.. automodule:: operational_analysis.toolkits.unit_conversion
    :members:
    :undoc-members:
    :show-inheritance:


Plotting
---------------------

.. automodule:: operational_analysis.toolkits.pandas_plotting
    :members:
    :undoc-members:
    :show-inheritance:


Reanalysis Downloading
----------------------

.. automodule:: operational_analysis.toolkits.reanalysis_downloading
    :members:
    :undoc-members:
    :show-inheritance:
.. _install::


.. ::

    # with overline, for parts
    * with overline, for chapters
    =, for sections
    -, for subsections
    ^, for subsubsections
    ", for paragraphs

Install
*******

.. mdinclude:: ../readme.md
