# Change Log

All notable changes to this project will be documented in this file.
This project adheres to [Semantic Versioning](http://semver.org/).
Formatted as described on [https://keepachangelog.com](https://keepachangelog.com/en/1.0.0/).

## [Unreleased]

## [1.1.4] (2022-01-14)

### Added

- 2020.1.3 version of wflow model ([#270](https://github.com/eWaterCycle/ewatercycle/issues/270))

### Changed

- Replace Cartesius section in system setup docs with Snellius ([#273](https://github.com/eWaterCycle/ewatercycle/issues/273))

### Fixed

- Test suite fails with fresh conda env ([#275](https://github.com/eWaterCycle/ewatercycle/issues/275))
- incompatible numpy typings ([#285](https://github.com/eWaterCycle/ewatercycle/issues/285))

## [1.1.3] (2021-10-18)

### Added

- 2020.1.2 version of wflow model ([#268](https://github.com/eWaterCycle/ewatercycle/pull/268))
- Document how to add a new version of a model ([#266](https://github.com/eWaterCycle/ewatercycle/pull/266))

## [1.1.2] (2021-09-29)

### Added

- Type information according to [PEP-0561](https://www.python.org/dev/peps/pep-0561/)
- Pre-commit hooks and black formatting ([#111](https://github.com/eWaterCycle/ewatercycle/issues/111))

### Changed

- Timeout for model setup set to 5 minutes ([#244](https://github.com/eWaterCycle/ewatercycle/issues/244))
- Use mamba for installation instructions ([#136](https://github.com/eWaterCycle/ewatercycle/issues/136))
- Use [version 1.2.0](https://github.com/citation-file-format/citation-file-format/releases/tag/1.2.0) of CITATION.cff format
- Moved package to src/ ([#228](https://github.com/eWaterCycle/ewatercycle/issues/228))

### Fixed

- Name particle in CITATION.cff ([#204](https://github.com/eWaterCycle/ewatercycle/issues/204))
- Build Sphinx locally with config file ([#169](https://github.com/eWaterCycle/ewatercycle/issues/169))
- Type errors in notebooks ([#262](https://github.com/eWaterCycle/ewatercycle/issues/262))
- Lisflood.finalize() ([#257](https://github.com/eWaterCycle/ewatercycle/issues/257))

## [1.1.1] (2021-08-10)

### Fixed

- Zenodo DOI

## [1.1.0] (2021-08-10)

### Added

- Column name argument to `get_grdc_data()` ([#190](https://github.com/eWaterCycle/ewatercycle/issues/190))
- Copy to clipboard button to documentation ([#216](https://github.com/eWaterCycle/ewatercycle/issues/216))

### Changed

- Compatible with ESMValTool 2.3 . Older versions (<2.3) of ESMValTool are no longer supported. ([#219](https://github.com/eWaterCycle/ewatercycle/issues/219))
- README, CONTRIBUTING, CHANGELOG reformated from RestructedText to Markdown ([#199](https://github.com/eWaterCycle/ewatercycle/issues/199))

### Fixed

- ParameterSet can be outside CFG['parametersets_dir'] ([#217](https://github.com/eWaterCycle/ewatercycle/issues/217))
- Link to nbviewer ([#202](https://github.com/eWaterCycle/ewatercycle/issues/202))
- Pinned esmpy as temporary workaround for single CPU affinity ([#234](https://github.com/eWaterCycle/ewatercycle/issues/234))

### Removed

- Unused esmvaltool_config field in CFG ([#152](https://github.com/eWaterCycle/ewatercycle/issues/152))

## [1.0.0] (2021-07-21)

### Added

- Documentation
  - Example notebooks
  - Setup guide
        ([\#120](https://github.com/eWaterCycle/ewatercycle/issues/120))
  - HPC cluster guide
- Forcing generation using [ESMValTool](https://www.esmvaltool.org/)
    ([\#28](https://github.com/eWaterCycle/ewatercycle/issues/28),
    [\#87](https://github.com/eWaterCycle/ewatercycle/issues/87),)
- Available parameter sets
    ([\#118](https://github.com/eWaterCycle/ewatercycle/issues/118))
- [PyMT](https://pymt.readthedocs.io/) inspired interface for
    following models
  - LISFLOOD
  - MARRMoT M01 and M14
  - PCR-GLOBWB
  - wflow
- Model methods to get and set values based on spatial coordinates
    ([\#53](https://github.com/eWaterCycle/ewatercycle/issues/53),
    [\#140](https://github.com/eWaterCycle/ewatercycle/issues/140))
- Model method to get value as a xarray dataset
    ([\#36](https://github.com/eWaterCycle/ewatercycle/issues/36))
- Containerized models using
    [grpc4bmi](https://github.com/eWaterCycle/grpc4bmi)
- Configuration files for system setup
- Hydrograph plotting
    ([\#54](https://github.com/eWaterCycle/ewatercycle/issues/54))
- Typings
- iso8601 time format
    ([\#90](https://github.com/eWaterCycle/ewatercycle/issues/90))

### Changed

- GRDC returns Pandas dataframe and metadata dict instead of xarray
    dataset
    ([\#109](https://github.com/eWaterCycle/ewatercycle/issues/109))

## [0.2.0] (2021-03-17)

### Added

- Observations from GRDC and USGS
- Empty Python project directory structure
- Added symlink based data files copier

[Unreleased]: https://github.com/eWaterCycle/ewatercycle/compare/1.1.4...HEAD
[1.1.4]: https://github.com/eWaterCycle/ewatercycle/compare/1.1.3...1.1.4
[1.1.3]: https://github.com/eWaterCycle/ewatercycle/compare/1.1.2...1.1.3
[1.1.2]: https://github.com/eWaterCycle/ewatercycle/compare/1.1.1...1.1.2
[1.1.1]: https://github.com/eWaterCycle/ewatercycle/compare/1.1.0...1.1.1
[1.1.0]: https://github.com/eWaterCycle/ewatercycle/compare/1.0.0...1.1.0
[1.0.0]: https://github.com/eWaterCycle/ewatercycle/compare/0.2.x-observation_data...1.0.0
[0.2.0]: https://github.com/eWaterCycle/ewatercycle/releases/tag/0.2.x-observation_data
# ewatercycle

![image](https://github.com/eWaterCycle/ewatercycle/raw/main/docs/examples/logo.png)

A Python package for running hydrological models.

[![image](https://github.com/eWaterCycle/ewatercycle/actions/workflows/ci.yml/badge.svg)](https://github.com/eWaterCycle/ewatercycle/actions/workflows/ci.yml)
[![image](https://sonarcloud.io/api/project_badges/measure?project=eWaterCycle_ewatercycle&metric=alert_status)](https://sonarcloud.io/dashboard?id=eWaterCycle_ewatercycle)
[![image](https://sonarcloud.io/api/project_badges/measure?project=eWaterCycle_ewatercycle&metric=coverage)](https://sonarcloud.io/component_measures?id=eWaterCycle_ewatercycle&metric=coverage)
[![Documentation Status](https://readthedocs.org/projects/ewatercycle/badge/?version=latest)](https://ewatercycle.readthedocs.io/en/latest/?badge=latest)
[![PyPI](https://img.shields.io/pypi/v/ewatercycle)](https://pypi.org/project/ewatercycle/)
[![image](https://img.shields.io/badge/fair--software.eu-%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8F%20%20%E2%97%8B-yellow)](https://fair-software.eu)
[![image](https://zenodo.org/badge/DOI/10.5281/zenodo.5119389.svg)](https://doi.org/10.5281/zenodo.5119389)
[![Research Software Directory Badge](https://img.shields.io/badge/rsd-ewatercycle-00a3e3.svg)](https://www.research-software.nl/software/ewatercycle)

The eWaterCycle package makes it easier to use hydrological models
without having intimate knowledge about how to install and run the
models.

- Uses container for running models in an isolated and portable way
    with [grpc4bmi](https://github.com/eWaterCycle/grpc4bmi)
- Generates rain and sunshine required for the model using
    [ESMValTool](https://www.esmvaltool.org/)
- Supports observation data from [GRDC or
    USGS](https://ewatercycle.readthedocs.io/en/latest/observations.html)
- Exposes [simple
    interface](https://ewatercycle.readthedocs.io/en/latest/examples/ewatercycle_api_notebook.html)
    to quickly get up and running

## Install

The ewatercycle package needs some geospatial non-python packages to
generate forcing data. It is preferred to create a Conda environment to
install those dependencies:

```shell
wget https://raw.githubusercontent.com/eWaterCycle/ewatercycle/main/environment.yml
conda install mamba -n base -c conda-forge -y
mamba env create --file environment.yml
conda activate ewatercycle
```

The ewatercycle package is installed with

```shell
pip install ewatercycle
```

Besides installing software you will need to create a configuration
file, download several data sets and get container images. See the
[system setup
chapter](https://ewatercycle.readthedocs.org/en/latest/system_setup.html)
for instructions.

## Usage

Example using the [Marrmot M14
(TOPMODEL)](https://github.com/wknoben/MARRMoT/blob/master/MARRMoT/Models/Model%20files/m_14_topmodel_7p_2s.m)
hydrological model on Merrimack catchment to generate forcing, run it
and produce a hydrograph.

```python
import pandas as pd
import ewatercycle.analysis
import ewatercycle.forcing
import ewatercycle.models
import ewatercycle.observation.grdc

forcing = ewatercycle.forcing.generate(
    target_model='marrmot',
    dataset='ERA5',
    start_time='2010-01-01T00:00:00Z',
    end_time='2010-12-31T00:00:00Z',
    shape='Merrimack/Merrimack.shp'
)

model = ewatercycle.models.MarrmotM14(version="2020.11", forcing=forcing)

cfg_file, cfg_dir = model.setup(
    threshold_flow_generation_evap_change=0.1,
    leakage_saturated_zone_flow_coefficient=0.99,
    zero_deficit_base_flow_speed=150.0,
    baseflow_coefficient=0.3,
    gamma_distribution_phi_parameter=1.8
)

model.initialize(cfg_file)

observations_df, station_info = ewatercycle.observation.grdc.get_grdc_data(
    station_id=4147380,
    start_time=model.start_time_as_isostr,
    end_time=model.end_time_as_isostr,
    column='observation',
)

simulated_discharge = []
timestamps = []
while (model.time < model.end_time):
    model.update()
    value = model.get_value('flux_out_Q')[0]
    # flux_out_Q unit conversion factor from mm/day to m3/s
    area = 13016500000.0  # from shapefile in m2
    conversion_mmday2m3s = 1 / (1000 * 24 * 60 * 60)
    simulated_discharge.append(value * area * conversion_mmday2m3s)
    timestamps.append(model.time_as_datetime.date())
simulated_discharge_df = pd.DataFrame({'simulated': simulated_discharge}, index=pd.to_datetime(timestamps))

ewatercycle.analysis.hydrograph(simulated_discharge_df.join(observations_df), reference='observation')

model.finalize()
```

More examples can be found in the
[documentation](https://ewatercycle.readthedocs.io).

## Contributing

If you want to contribute to the development of ewatercycle package,
have a look at the [contribution guidelines](CONTRIBUTING.md).

## License

Copyright (c) 2018, Netherlands eScience Center & Delft University of
Technology

Apache Software License 2.0
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our
project and our community a harassment-free experience for everyone,
regardless of age, body size, disability, ethnicity, gender identity and
expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

- Using welcoming and inclusive language
- Being respectful of differing viewpoints and experiences
- Gracefully accepting constructive criticism
- Focusing on what is best for the community
- Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

- The use of sexualized language or imagery and unwelcome sexual
    attention or advances
- Trolling, insulting/derogatory comments, and personal or political
    attacks
- Public or private harassment
- Publishing others\' private information, such as a physical or
    electronic address, without explicit permission
- Other conduct which could reasonably be considered inappropriate in
    a professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of
acceptable behavior and are expected to take appropriate and fair
corrective action in response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit,
or reject comments, commits, code, wiki edits, issues, and other
contributions that are not aligned to this Code of Conduct, or to ban
temporarily or permanently any contributor for other behaviors that they
deem inappropriate, threatening, offensive, or harmful.

## Scope

This Code of Conduct applies both within project spaces and in public
spaces when an individual is representing the project or its community.
Examples of representing a project or community include using an
official project e-mail address, posting via an official social media
account, or acting as an appointed representative at an online or
offline event. Representation of a project may be further defined and
clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may
be reported by contacting the project team at
<s.verhoeven@esciencecenter.nl>. All complaints will be reviewed and
investigated and will result in a response that is deemed necessary and
appropriate to the circumstances. The project team is obligated to
maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted
separately.

Project maintainers who do not follow or enforce the Code of Conduct in
good faith may face temporary or permanent repercussions as determined
by other members of the project\'s leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor
Covenant](https://www.contributor-covenant.org), version 1.4, available
at
<https://www.contributor-covenant.org/version/1/4/code-of-conduct.html>
# Contributing guidelines

We welcome any kind of contributions to our software, from simple
comment or question to a full fledged [pull
request](https://help.github.com/articles/about-pull-requests/). Please
read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. you have a question;
2. you think you may have found a bug (including unexpected behavior);
3. you want to make some kind of change to the code base (e.g. to fix a
    bug, to add a new feature, to update documentation).
4. you want to make a release

The sections below outline the steps in each case.

## You have a question

1. use the search functionality
    [here](https://github.com/eWaterCycle/ewatercycle/issues) to see if
    someone already filed the same issue;
2. if your issue search did not yield any relevant results, make a new
    issue;
3. apply the \"Question\" label; apply other labels when relevant.

## You think you may have found a bug

1. use the search functionality
    [here](https://github.com/eWaterCycle/ewatercycle/issues) to see if
    someone already filed the same issue;
2. if your issue search did not yield any relevant results, make a new
    issue, making sure to provide enough information to the rest of the
    community to understand the cause and context of the problem.
    Depending on the issue, you may want to include: - the [SHA
    hashcode](https://help.github.com/articles/autolinked-references-and-urls/#commit-shas)
    of the commit that is causing your problem; - some identifying
    information (name and version number) for dependencies you\'re
    using; - information about the operating system;
3. apply relevant labels to the newly created issue.

## You want to make some kind of change to the code base

1. (**important**) announce your plan to the rest of the community
    *before you start working*. This announcement should be in the form
    of a (new) issue;
2. (**important**) wait until some kind of consensus is reached about
    your idea being a good idea;
3. if needed, fork the repository to your own Github profile and create
    your own feature branch off of the latest main commit. While working
    on your feature branch, make sure to stay up to date with the main
    branch by pulling in changes, possibly from the \'upstream\'
    repository (follow the instructions
    [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/)
    and [here](https://help.github.com/articles/syncing-a-fork/));
4. install the package in editable mode and its dependencies with
    `pip3 install -e .[dev]`;
4. make sure pre commit hook is installed by running `pre-commit install`, causes linting and formatting to be applied during commit;
5. make sure the existing tests still work by running `pytest`;
6. make sure the existing documentation can still by generated without
    warnings by running `cd docs && make html`. [Pandoc](https://pandoc.org/) is required to generate docs, it can be installed with ``conda install -c conda-forge pandoc`` ;
7. add your own tests (if necessary);
8. update or expand the documentation; Please add [Google Style Python
    docstrings](https://google.github.io/styleguide/pyguide.html#38-comments-and-docstrings).
9. [push](http://rogerdudler.github.io/git-guide/) your feature branch
    to (your fork of) the ewatercycle repository on GitHub;
10. create the pull request, e.g. following the instructions
    [here](https://help.github.com/articles/creating-a-pull-request/).

In case you feel like you\'ve made a valuable contribution, but you
don\'t know how to write or run tests for it, or how to generate the
documentation: don\'t let this discourage you from making the pull
request; we can help you! Just go ahead and submit the pull request, but
keep in mind that you might be asked to append additional commits to
your pull request.

## You want to make a release

This section is for maintainers of the package.

1. Checkout ``HEAD`` of ``main`` branch with ``git checkout main`` and ``git pull``.
2. Determine what new version (major, minor or patch) to use. Package uses `semantic versioning <https://semver.org>`_.
3. Run ``bump2version <major|minor|patch>`` to update version in package files.
4. Update CHANGELOG.md with changes between current and new version.
5. Make sure pre-commit hooks are green for all files by running ``pre-commit run --all-files``.
6. Commit & push changes to GitHub.
7. Wait for [GitHub
    actions](https://github.com/eWaterCycle/ewatercycle/actions?query=branch%3Amain+)
    to be completed and green.

8. Create a [GitHub release](https://github.com/eWaterCycle/ewatercycle/releases/new)

    - Use version as title and tag version.
    - As description use intro text from README.md (to give context to
        Zenodo record) and changes from CHANGELOG.md

9. Create a PyPI release.

    1. Create distribution archives with `python3 -m build`.
    2. Upload archives to PyPI with `twine upload dist/*` (use your
        personal PyPI account).

10. Verify

    1. Has [new Zenodo
        record](https://zenodo.org/search?page=1&size=20&q=ewatercycle)
        been created?
    2. Has [stable
        ReadTheDocs](https://ewatercycle.readthedocs.io/en/stable/) been
        updated?
    3. Can new version be installed with pip using
        `pip3 install ewatercycle==<new version>`?

11. Celebrate
This is a small collection of example notebooks for the eWaterCycle package.

If GitHub fails to render these notebooks, you can also view them using the Jupyter nbviewer here: https://nbviewer.jupyter.org/github/eWaterCycle/ewatercycle/tree/main/docs/examples/
System setup
============

To use eWaterCycle package you need to setup the system with software
and data.

This chapter is for system administrators or Research Software Engineers who need to set up a system for the eWaterCycle platform.

The setup steps:

1.  `Conda environment <#conda-environment>`__
2.  `Install ewatercycle package <#install-ewatercycle-package>`__
3.  `Configure ESMValTool <#configure-ESMValTool>`__
4.  `Download climate data <#download-climate-data>`__
5.  `Install container engine <#install-container-engine>`__
6.  `Configure ewatercycle <#configure-ewatercycle>`__
7.  `Model container images <#model-container-images>`__
8.  `Download example parameter sets <#download-example-parameter-sets>`__
9.  `Prepare other parameter sets <#prepare-other-parameter-sets>`_
10. `Download example forcing <#download-example-forcing>`__
11. `Download observation data <#download-observation-data>`__

Conda environment
-----------------

The eWaterCycle Python package uses a lot of geospatial dependencies
which can be installed using `Conda <https://conda.io/>`__ package
management system.

Install Conda by using the `miniconda
installer <https://docs.conda.io/en/latest/miniconda.html>`__.

After conda is installed you can install the software dependencies with
a `conda environment
file <https://github.com/eWaterCycle/ewatercycle/blob/main/environment.yml>`__.

.. code:: shell

    wget https://raw.githubusercontent.com/eWaterCycle/ewatercycle/main/environment.yml
    conda install mamba -n base -c conda-forge -y
    mamba env create --file environment.yml
    conda activate ewatercycle

Do not forget that any terminal or Jupyter kernel should activate the conda environment before the eWaterCycle Python package can be used.

Install eWaterCycle package
---------------------------

The Python package can be installed using pip

.. code:: shell

    pip install ewatercycle


Configure ESMValTool
--------------------

ESMValTool is used to generate forcing (temperature, precipitation,
etc.) files from climate data for hydrological models. The
ESMValTool has been installed as a dependency of the package.

See https://docs.esmvaltool.org/en/latest/quickstart/configuration.html
how configure ESMValTool.

Download climate data
---------------------

The ERA5 and ERA-Interim data can be used to generate
forcings.

ERA5
~~~~

To download ERA5 data files you can use the
`era5cli <https://era5cli.readthedocs.io/>`__ tool.

.. code:: shell

    pip install era5cli

Follow `instructions <https://era5cli.readthedocs.io/en/stable/instructions.html>`_ to get access to data.

As an example, the hourly ERA5 data for the years 1990
and 1991 and for variables pr, psl, tas, taxmin, tasmax, tdps, uas,
vas, rsds, rsdt and fx orog are downloaded as:

.. code:: shell

    cd <ESMValTool ERA5 raw directory for example /projects/0/wtrcycle/comparison/rawobs/Tier3/ERA5/1>
    era5cli hourly --startyear 1990 --endyear 1991 --variables total_precipitation
    era5cli hourly --startyear 1990 --endyear 1991 --variables mean_sea_level_pressure
    era5cli hourly --startyear 1990 --endyear 1991 --variables 2m_temperature
    era5cli hourly --startyear 1990 --endyear 1991 --variables minimum_2m_temperature_since_previous_post_processing
    era5cli hourly --startyear 1990 --endyear 1991 --variables maximum_2m_temperature_since_previous_post_processing
    era5cli hourly --startyear 1990 --endyear 1991 --variables 2m_dewpoint_temperature
    era5cli hourly --startyear 1990 --endyear 1991 --variables 10m_u_component_of_wind
    era5cli hourly --startyear 1990 --endyear 1991 --variables 10m_v_component_of_wind
    era5cli hourly --startyear 1990 --endyear 1991 --variables surface_solar_radiation_downwards
    era5cli hourly --startyear 1990 --endyear 1991 --variables toa_incident_solar_radiation
    era5cli hourly --startyear 1990 --endyear 1991 --variables orography
    cd -

The hourly data needs need be converted to daily using a `ESMValTool recipe <https://docs.esmvaltool.org/en/latest/input.html#cmorization-as-a-fix>`_

.. code:: shell

    esmvaltool run cmorizers/recipe_era5.yml

ERA-Interim
~~~~~~~~~~~

ERA-Interim has been superseeded by ERA5, but could be useful for
reproduction studies and its smaller size. The ERA-Interim data files
can be downloaded at
https://www.ecmwf.int/en/forecasts/datasets/reanalysis-datasets/era-interim

Or you can use the `download_era_interim.py <https://github.com/ESMValGroup/ESMValTool/blob/main/esmvaltool/cmorizers/obs/download_scripts/download_era_interim.py>`_
script to download ERA-Interim data files. See first lines of script for documentation.
The files should be downloaded to the ESMValTool ERA-Interim raw directory for example ``/projects/0/wtrcycle/comparison/rawobs/Tier3/ERA-Interim``.

The ERA5-Interim raw data files need to be cmorized using `script <https://docs.esmvaltool.org/en/latest/input.html#using-a-cmorizer-script>`_:

.. code:: shell

    cmorize_obs -o ERA-Interim

Install container engine
------------------------

In eWaterCycle package, the hydrological models are run in containers
with engines like `Singularity <https://singularity.lbl.gov/>`__ or
`Docker <https://www.docker.com/>`__. At least Singularity or Docker
should be installed.

Installing a container engine requires root permission on the machine.

Singularity
~~~~~~~~~~~

Install Singularity using
`instructions <https://singularity.hpcng.org/user-docs/master/quick_start.html>`__.

Docker
~~~~~~

Install Docker using
`instructions <https://docs.docker.com/engine/install/>`__. Docker
should be configured so it can be `called without
sudo <https://docs.docker.com/engine/install/linux-postinstall/#manage-docker-as-a-non-root-user>`__

Configure eWaterCycle
---------------------

The eWaterCycle package simplifies the API by reading some of the
directories and settings from a configuration file.

The configuration can be set in Python with

.. code:: ipython3

    import logging
    logging.basicConfig(level=logging.INFO)
    import ewatercycle
    import ewatercycle.parameter_sets
    # Which container engine is used to run the hydrological models
    ewatercycle.CFG['container_engine'] = 'singularity'  # or 'docker'
    # If container_engine==singularity then where can the singularity images files (*.sif) be found.
    ewatercycle.CFG['singularity_dir'] = './singularity-images'
    # Directory in which output of model runs is stored. Each model run will generate a sub directory inside output_dir
    ewatercycle.CFG['output_dir'] = './'
    # Where can GRDC observation files (<station identifier>_Q_Day.Cmd.txt) be found.
    ewatercycle.CFG['grdc_location'] = './grdc-observations'
    # Where can parameters sets prepared by the system administator be found
    ewatercycle.CFG['parameterset_dir'] = './parameter-sets'
    # Where is the configuration saved or loaded from
    ewatercycle.CFG['ewatercycle_config'] = './ewatercycle.yaml'

and then written to disk with

.. code:: ipython3

    ewatercycle.CFG.save_to_file()

Later it can be loaded by using:

.. code:: ipython3

    ewatercycle.CFG.load_from_file('./ewatercycle.yaml')

To make the ewatercycle configuration load by default for current user
it should be copied to ``~/.config/ewatercycle/ewatercycle.yaml`` .

To make the ewatercycle configuration available to all users on the
system it should be copied to ``/etc/ewatercycle.yaml`` .

Configuration file for Snellius system
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users part of the eWaterCycle project can use the following configurations on the `Snellius system of
SURF <https://servicedesk.surfsara.nl/wiki/display/WIKI/Snellius>`_:

.. code:: yaml

   container_engine: singularity
   singularity_dir: /projects/0/wtrcycle/singularity-images
   output_dir: /scratch-shared/ewatercycle
   grdc_location:  /projects/0/wtrcycle/GRDC/GRDC_GCOSGTN-H_27_03_2019
   parameterset_dir: /projects/0/wtrcycle/parameter-sets

The `/scratch-shared/ewatercycle` output directory will be automatically removed if its content is older than 14 days.
If the output directory is missing it can be recreated with

.. code:: shell

    mkdir /scratch-shared/ewatercycle
    chgrp wtrcycle /scratch-shared/ewatercycle
    chmod 2770 /scratch-shared/ewatercycle

Configuration file for ewatecycle Jupyter machine
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Users can use the following configurations on systems constructed with eWaterCycle application on SURF Research
Cloud:

.. code:: yaml

   container_engine: singularity
   singularity_dir: /mnt/data/singularity-images
   output_dir: /scratch
   grdc_location: /mnt/data/GRDC
   parameterset_dir: /mnt/data/parameter-sets

Model container images
----------------------

As hydrological models run in containers, their container images should be
made available on the system.

The names of the images can be found in the ``ewatercycle.models.*``
classes.

Docker
~~~~~~

Docker images will be downloaded with ``docker pull``:

.. code:: shell

    docker pull ewatercycle/lisflood-grpc4bmi:20.10
    docker pull ewatercycle/marrmot-grpc4bmi:2020.11
    docker pull ewatercycle/pcrg-grpc4bmi:setters
    docker pull ewatercycle/wflow-grpc4bmi:2020.1.1
    docker pull ewatercycle/wflow-grpc4bmi:2020.1.2
    docker pull ewatercycle/wflow-grpc4bmi:2020.1.3

Singularity
~~~~~~~~~~~

Singularity images should be stored in configured directory
(``ewatercycle.CFG['singularity_dir']``) and can build from Docker with:

.. code:: shell

    cd {ewatercycle.CFG['singularity_dir']}
    singularity build ewatercycle-lisflood-grpc4bmi_20.10.sif docker://ewatercycle/lisflood-grpc4bmi:20.10
    singularity build ewatercycle-marrmot-grpc4bmi_2020.11.sif docker://ewatercycle/marrmot-grpc4bmi:2020.11
    singularity build ewatercycle-pcrg-grpc4bmi_setters.sif docker://ewatercycle/pcrg-grpc4bmi:setters
    singularity build ewatercycle-wflow-grpc4bmi_2020.1.1.sif docker://ewatercycle/wflow-grpc4bmi:2020.1.1
    singularity build ewatercycle-wflow-grpc4bmi_2020.1.2.sif docker://ewatercycle/wflow-grpc4bmi:2020.1.2
    singularity build ewatercycle-wflow-grpc4bmi_2020.1.3.sif docker://ewatercycle/wflow-grpc4bmi:2020.1.3
    cd -

Download example parameter sets
-------------------------------

To quickly run the models it is advised to setup a example parameter
sets for each model.

.. code:: ipython3

    ewatercycle.parameter_sets.download_example_parameter_sets()


.. parsed-literal::

    INFO:ewatercycle.parameter_sets._example:Downloading example parameter set wflow_rhine_sbm_nc to /home/verhoes/git/eWaterCycle/ewatercycle/docs/examples/parameter-sets/wflow_rhine_sbm_nc...
    INFO:ewatercycle.parameter_sets._example:Download complete.
    INFO:ewatercycle.parameter_sets._example:Adding parameterset wflow_rhine_sbm_nc to ewatercycle.CFG...
    INFO:ewatercycle.parameter_sets._example:Downloading example parameter set pcrglobwb_rhinemeuse_30min to /home/verhoes/git/eWaterCycle/ewatercycle/docs/examples/parameter-sets/pcrglobwb_rhinemeuse_30min...
    INFO:ewatercycle.parameter_sets._example:Download complete.
    INFO:ewatercycle.parameter_sets._example:Adding parameterset pcrglobwb_rhinemeuse_30min to ewatercycle.CFG...
    INFO:ewatercycle.parameter_sets._example:Downloading example parameter set lisflood_fraser to /home/verhoes/git/eWaterCycle/ewatercycle/docs/examples/parameter-sets/lisflood_fraser...
    INFO:ewatercycle.parameter_sets._example:Download complete.
    INFO:ewatercycle.parameter_sets._example:Adding parameterset lisflood_fraser to ewatercycle.CFG...
    INFO:ewatercycle.parameter_sets:3 example parameter sets were downloaded
    INFO:ewatercycle.config._config_object:Config written to /home/verhoes/git/eWaterCycle/ewatercycle/docs/examples/ewatercycle.yaml
    INFO:ewatercycle.parameter_sets:Saved parameter sets to configuration file /home/verhoes/git/eWaterCycle/ewatercycle/docs/examples/ewatercycle.yaml


Example parameter sets have been downloaded and added to the
configuration file.

.. code:: shell

    cat ./ewatercycle.yaml


.. parsed-literal::

    container_engine: null
    grdc_location: None
    output_dir: None
    parameter_sets:
      lisflood_fraser:
        config: lisflood_fraser/settings_lat_lon-Run.xml
        directory: lisflood_fraser
        doi: N/A
        supported_model_versions: !!set {'20.10': null}
        target_model: lisflood
      pcrglobwb_rhinemeuse_30min:
        config: pcrglobwb_rhinemeuse_30min/setup_natural_test.ini
        directory: pcrglobwb_rhinemeuse_30min
        doi: N/A
        supported_model_versions: !!set {setters: null}
        target_model: pcrglobwb
      wflow_rhine_sbm_nc:
        config: wflow_rhine_sbm_nc/wflow_sbm_NC.ini
        directory: wflow_rhine_sbm_nc
        doi: N/A
        supported_model_versions: !!set {2020.1.1: null}
        target_model: wflow
    parameterset_dir: /home/verhoes/git/eWaterCycle/ewatercycle/docs/examples/parameter-sets
    singularity_dir: None


.. code:: ipython3

    ewatercycle.parameter_sets.available_parameter_sets()


.. parsed-literal::

    ('lisflood_fraser', 'pcrglobwb_rhinemeuse_30min', 'wflow_rhine_sbm_nc')



.. code:: ipython3

    parameter_set = ewatercycle.parameter_sets.get_parameter_set('pcrglobwb_rhinemeuse_30min')
    print(parameter_set)


.. parsed-literal::

    Parameter set
    -------------
    name=pcrglobwb_rhinemeuse_30min
    directory=/home/verhoes/git/eWaterCycle/ewatercycle/docs/examples/parameter-sets/pcrglobwb_rhinemeuse_30min
    config=/home/verhoes/git/eWaterCycle/ewatercycle/docs/examples/parameter-sets/pcrglobwb_rhinemeuse_30min/setup_natural_test.ini
    doi=N/A
    target_model=pcrglobwb
    supported_model_versions={'setters'}

The ``parameter_set`` variable can be passed to a model class
constructor.

Prepare other parameter sets
----------------------------

The example parameter sets downloaded in the previous section are nice to show off the platform features but are a bit small.
To perform more advanced experiments, additional parameter sets are needed.
Users could use :py:class:`ewatercycle.parameter_sets.ParameterSet` to construct parameter sets themselves.
Or they can be made available via :py:func:`ewatercycle.parameter_sets.available_parameter_sets` and :py:func:`ewatercycle.parameter_sets.get_parameter_set` by extending the configuration file (ewatercycle.yaml).

A new parameter set should be added as a key/value pair in the ``parameter_sets`` map of the configuration file.
The key should be a unique string on the current system.
The value is a dictionary with the following items:

* directory: Location on disk where files of the parameter set are stored. If Path is relative then relative to :py:const:`ewatercycle.CFG['parameterset_dir']`.
* config: Model configuration file which uses files from directory. If Path is relative then relative to :py:const:`ewatercycle.CFG['parameterset_dir']`.
* doi: Persistent identifier of the parameter set. For example a DOI for a Zenodo record.
* target_model: Name of the model that parameter set can work with
* supported_model_versions: Set of model versions that are supported by this parameter set. If not set then parameter set will be supported by all versions of model

For example the parameter set for PCR-GLOBWB from https://doi.org/10.5281/zenodo.1045339 after downloading and unpacking to ``/data/pcrglobwb2_input/`` could be added with following config:

.. code:: yaml

    pcrglobwb_rhinemeuse_30min:
        directory: /data/pcrglobwb2_input/global_30min/
        config: /data/pcrglobwb2_input/global_30min/iniFileExample/setup_30min_non-natural.ini
        doi: https://doi.org/10.5281/zenodo.1045339
        target_model: pcrglobwb
        supported_model_versions: !!set {setters: null}


Download example forcing
------------------------

To be able to run the Marrmot example notebooks you need a forcing file.
You can use ``ewatercycle.forcing.generate()`` to make it or use an
already prepared `forcing
file <https://github.com/wknoben/MARRMoT/blob/master/BMI/Config/BMI_testcase_m01_BuffaloRiver_TN_USA.mat>`__.

.. code:: shell

    cd docs/examples
    wget https://github.com/wknoben/MARRMoT/raw/master/BMI/Config/BMI_testcase_m01_BuffaloRiver_TN_USA.mat
    cd -

Download observation data
-------------------------

Observation data is needed to calculate metrics of the model performance or plot a hydrograph . The
ewatercycle package can use `Global Runoff Data Centre
(GRDC) <https://www.bafg.de/GRDC>`__ or `U.S. Geological Survey Water
Services (USGS) <https://waterservices.usgs.gov/>`__ data.

The GRDC daily data files can be ordered at
https://www.bafg.de/GRDC/EN/02_srvcs/21_tmsrs/riverdischarge_node.html.

The GRDC files should be stored in ``ewatercycle.CFG['grdc_location']``
directory.

Examples
========

.. toctree::

  Generate forcing <examples/generate_forcing.ipynb>
  LISFLOOD <examples/lisflood.ipynb>
  MARRMoT M01 <examples/MarrmotM01.ipynb>
  MARRMoT M14 <examples/MarrmotM14.ipynb>
  PCRGlobWB <examples/pcrglobwb.ipynb>
  Wflow <examples/wflow.ipynb>
  Irrigation experiment <examples/Irrigation.ipynb>
Adding a model
==============

Integrating a new model into the eWaterCycle system involves the following steps:

* Create model as subclass of ``AbstractModel`` (``src/ewatercycle/models/abstract.py``)
* Import model in ``src/ewatercycle/models/__init__.py``
* Add ``src/ewatercycle/forcing/<model>.py``
* Register model in ``src/ewatercycle/forcing/__init__.py:FORCING_CLASSES``
* Add model to ``docs/conf.py``
* Write example notebook
* Write tests?
* If model needs custom parameter set class add it in ``src/ewatercycle/parameter_sets/_<model name>.py``
* Add example parameter set in ``src/ewatercycle/parameter_sets/__init__.py``
* Add container image to :doc:`system_setup`
* Add container image to infrastructure data preparation scripts_

We will expand this documentation in due time.

Adding a new version of a model
-------------------------------

A model can have different versions.
A model version in the eWaterCycle Python package corresponds to the tag of Docker image and the version in a Singularity container image filename.
The version of the container image should preferably be one of release versions of the model code. Alternativly the version could be the name of a feature branch or a date.

Also parameter sets can be specify which versions of a model they support.

To add a new version of a model involves the following steps:

Create container image
~~~~~~~~~~~~~~~~~~~~~~

* Create Docker container image named ``ewatercycle/<model>-grpc4bmi:<version>`` with `grpc4bmi server running as entrypoint <https://grpc4bmi.readthedocs.io/en/latest/container/building.html>`_
* Host Docker container image on `Docker Hub <https://hub.docker.com/u/ewatercycle>`_
* Create Singularity image from Docker with ``singularity build ./ewatercycle-<model>-grpc4bmi_<version>.sif docker://ewatercycle/<model>-grpc4bmi:<version>``

Add to Python package
~~~~~~~~~~~~~~~~~~~~~

* Add container image to :doc:`system_setup` page by editing ``docs/system_setup.rst``
* In ``src/ewatercycle/models/<model>.py``

  * add new version to ``available_versions`` class property.
  * to ``__init__()`` method add support for new version

* Optionally: Add new version to existing example parameter set or add new parameter set in ``src/ewatercycle/parameter_sets/_<model>.py:example_parameter_sets()``
* Add new version to supported parameter sets in local eWaterCycle config file (``/etc/ewatercycle.yaml`` and ``~/.config/ewatercycle/ewatercycle.yaml``)
* Test it out locally
* Create pull request and get it merged
* Create new release of Python package. Done by package maintainers

Add to platform
~~~~~~~~~~~~~~~

For platform developers and deployers.

* Add Singularity image to dCache shared folder ``ewcdcache:/singularity-images/<model>-grpc4bmi_<version>.sif``
* Add container image to infrastructure repository

  * data preparation scripts_
  * `config generation <https://github.com/eWaterCycle/infra/blob/main/roles/ewatercycle/templates/ewatercycle.yaml.j2>`_

* Install version/branch of eWaterCycle Python package with new model version on any running virtual machines
* Optionally: Add example parameter set to `explorer catalog <https://github.com/eWaterCycle/TerriaMap/blob/ewatercycle-v8/wwwroot/init/ewatercycle.json>`_. The forcing, parameter set and model image should be available on Jupyter server connected to explorer.

.. _scripts: https://github.com/eWaterCycle/infra/tree/main/roles/prep_shared_data
Migrate from HPC to Cluster (Snellius) guide
=============================================
The HPC node jupyter.ewatercycle.org can be used for small test experiments, to do actual work you will need to run your notebook/script on the cluster (Snellius). On Snellius the forcing data is already present and many users can run jobs at the same time without interfering each other.

Familiarize yourself with Linux by reading this simple guide:

- https://maker.pro/linux/tutorial/basic-linux-commands-for-beginners

*************************
Migration Preparation
*************************

**1. Create Github repository**

Start by creating a Github repository to store (only) your code by following these guides:

- https://docs.github.com/en/github/getting-started-with-github/set-up-git
- https://docs.github.com/en/github/getting-started-with-github/create-a-repo

**2. Create Conda environment.yml** (not required)

For ease of transfer it can be helpful to create a environment.yml file. This file contains a list of all the packages you use for running code. This is good practice because it allows users of your Github repository to quickly install the necessary package requirements.

- https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#create-env-file-manually

**3. Copy files from HPC to Snellius**

To copy files from the eWaterCycle HPC to Snellius the following command example can be used:

- ``scp -r {YourUserNameOnTheHPC}@jupyter.ewatercycle.org:/mnt/{YourUserNameOnTheHPC}/{PathToFolder}/ /home/{YourUserNameOnTheSnellius}/{PathToFolder}/``

When prompted, enter your eWaterCycle HPC password.

********************
Login to Snellius
********************

**1. VPN Connection**

Cluster computer hosting institutes have a strict policy on which IP-addresses are allowed to connect with the Cluster (Snellius). For this reason you need to first establish a VPN connection to your University or Research Institute that has a whitelisted IP-address.

**2. MobaXterm**

To connects with Snellius a SSH client is required. One such free client is MobaXterm and can be downloaded here: https://mobaxterm.mobatek.net/.

- After installation open the client and click on the session tab (top left), click on SSH, at remote host fill in "snellius.surf.nl", tick the specify username box, fill in your Snellius username and click OK (bottom). Fill in the snellius password when prompted.

**3. Login Node & Compute Node**

Once you are logged in you are on the login node. This node should not be used to run scripts as it is only a portal to communicate with the compute nodes running on the background (the actual computers). The compute nodes are where you will do the calculations. We communicate with compute nodes using Bash (.sh) scripts. This will be explained later.

**4. Home Directory & Scratch Directory**

When you login you are directed to your Home Directory:

- ``/home/{YourUserNameOnTheSnellius}/``

The Home Directory has slower diskspeeds than the Scratch Directory. The Scratch Directory needs to be created using the following commands:

- ``cd /scratch-shared/``
- ``mkdir {YourUserNameOnTheSnellius}``

You can now access the Scratch Directory at ``/scratch/shared/{YourUserNameOnTheSnellius}/``. Best practice is to modify your code such that it first copies all the required files (excluding code) to the Scratch Directory, followed by running the code, after completion copying the files back to the Home Directory, and cleaning up the Scratch Directory.

*************************
First Run preparations
*************************
**1. Clone Github repository**

Clone Github repository containing scripts using:

- ``git clone https://github.com/example_user/example_repo``


**2. Install MiniConda**

Go to home directory:

- ``cd /home/username/``

Download MiniConda:

- ``wget https://repo.anaconda.com/miniconda/Miniconda3-latest-Linux-x86_64.sh``

Install MiniConda:

- ``bash Miniconda3-latest-Linux-x86_64.sh``

**Restart the connection with Snellius**

- ``conda update conda``

**3. Create Conda environment**

Create a Conda enviroment and install required packages following the description:

https://docs.conda.io/projects/conda/en/latest/user-guide/tasks/manage-environments.html#creating-an-environment-with-commands

Make sure that Jupyter Lab is installed in the Conda environment:

- ``wget https://raw.githubusercontent.com/eWaterCycle/ewatercycle/main/environment.yml``
- ``conda install mamba -n base -c conda-forge -y``
- ``mamba env create --file environment.yml``
- ``conda activate ewatercycle``
- ``conda install -c conda-forge jupyterlab``

Install eWatercycle package:

- ``pip install ewatercycle``

**4. Create Singularity Container**

On Snellius, Docker requires root access and can therefore not be used. Singularity is similar to, and integrates well with Docker.         It also requires root access, but it is pre-installed on the compute nodes on Snellius.

The first step to run the model on a compute node is thus to use singularity to create a Singularity image (``.sif`` file) based on the Docker image. This is done with (note the ``srun`` command to access the compute node):

- ``srun -N 1 -t 40 -p short singularity build --disable-cache ewatercycle-wflow-grpc4bmi.sif docker://ewatercycle/wflow-grpc4bmi:latest``

This is an example for the wflow_sbm model, change to the correct Docker container:

-  ``docker://ewatercycle/{model}-grpc4bmi:{version}``

**5. Adjust code to run Singularity container**

Code should be adjusted to run Singularity instead of Docker following:
::

    from grpc4bmi.bmi_client_singularity import BmiClientSingularity

    model = BmiClientSingularity(image='ewatercycle-wflow-grpc4bmi.sif', input_dirs=[input_dir], work_dir=work_dir)
    ...

**6. Adjust code to use Scratch directory**

Before running the model copy the model instance to the scratch directory:

``/scratch-shared/{YourUsernameOnTheSnellius}/``

Run the model from this directory and copy the output back to the home directory:

``/home/{YourUsernameOnTheSnellius}/``

Cleanup files in the scratch directory.


**************************************
Submitting Jupyter Job on Cluster node
**************************************
Here we briefly explain general SBATCH parameters and how to launch a Jupyter Lab environment on Snellius. Start by opening a text editor on Snellius (e.g. ``nano``) or (easier) your local machine (e.g. notepad). Copy the following text inside your text editor, edit the Conda environment name, and save as **run_jupyter_on_snellius.sh** (make sure the extension is ``.sh``):
::

    #!/bin/bash

    # Serve a jupyter lab environment from a compute node on Snellius
    # usage: sbatch run_jupyter_on_compute_node.sh

    # SLURM settings
    #SBATCH -J jupyter_lab
    #SBATCH -t 09:00:00
    #SBATCH -N 1
    #SBATCH -p normal
    #SBATCH --output=slurm_%j.out
    #SBATCH --error=slurm_%j.out

    # Use an appropriate conda environment
    . ~/miniconda3/etc/profile.d/conda.sh
    conda activate {YourEnvironmentName}

    # Some security: stop script on error and undefined variables
    set -euo pipefail

    # Specify (random) port to serve the notebook
    port=8123
    host=$(hostname -s)

    # Print command to create ssh tunnel in log file
    echo -e "

    Command to create ssh tunnel (run from another terminal session on your local machine):
    ssh -L ${port}:${host}:${port} $(whoami)@snellius.surf.nl
    Below, jupyter will print a number of addresses at which the notebook is served.
    Due to the way the tunnel is set up, only the latter option will work.
    It's the one that looks like
    http://127.0.0.1:${port}/?token=<long_access_token_very_important_to_copy_as_well>
    Copy this address in your local browser and you're good to go

    Starting notebooks server
    **************************************************
    "

    # Start the jupyter lab session

    jupyter lab --no-browser --port ${port} --ip=${host}

**Explanation of SBATCH Parameters**

- ``#SBATCH -J jupyter_lab``

Here you can set the job name.

- ``#SBATCH -t 09:00:00``

Here you specify job runtime. On the Snellius we have a budget, each half hour cpu runtime costs 1 point on the budget. A Node consists of 24 cores meaning that the specified runtime (9 hours) costs 24*2*9 points on the budget.

- ``#SBATCH -N 1``

Specifies the amount of nodes used by the run, keep at default value of 1.

- ``#SBATCH -p normal``

Specifies the type of Node, keep at default value of "normal".

- ``#SBATCH --output=slurm_%j.out``

Specifies the location and name of the job log file.

- More information on SBATCH parameters can be found here: https://servicedesk.surfsara.nl/wiki/display/WIKI/Creating+and+running+jobs

**Specifying job runtime**

Good practice for calculating job runtime is by for example running a model first for 1 year, calculate the time it takes. Multiply it by the total amount of years for your study. Add a time buffer of around 10-20 percent.

- For example: 1 year takes 2 hours, total run is 10 years, 20 hours total, add time buffer, estimated runtime equals 22-24 hours.

**Running the bash (.sh) script**

Enter this command to run the bash script:

- ``sbatch run_jupyter_on_snellius.sh``

(If you get DOS and UNIX linebreak errors, run the following command:)

- ``dos2unix run_jupyter_on_snellius.sh``



**Job control**

To view which jobs are running you can enter:

- ``squeue -u {YourUserNameOnTheSnellius}``

To cancel a running job you can enter:

- ``scancel {jobID}``

More information on job control can be found here: https://userinfo.surfsara.nl/systems/lisa/user-guide/creating-and-running-jobs#interacting

=====================================
Launching Jupyter Lab on Cluster Node
=====================================

**1. Open Slurm output log file**

- Open slurm output log file by double clicking in the file browser or by using a text editor (``nano``) and read the output carefully.

**2. Create ssh tunnel between local machine and cluster**

To create a ssh connection between your local machine and the cluster you need to open a command prompt interface on your local machine. For example ``PowerShell`` or ``cmd`` on Windows.

- copy the line ``ssh -L ${port}:${host}:${port} $(whoami)@snellius.surf.nl`` from the slurm log file (not the bash script) into the command prompt and run.

**3. Connect through browser**

- Open a browser (e.g. Chrome) and go to the url: ``localhost:8123/lab``

**4. Enter the access token**

- Copy the access token from the slurm otput log file and paste in the browser at access token or password.

You have now succesfully launched a Jupyter Lab environment on a cluster node.
.. ewatercycle documentation master file, created by
   sphinx-quickstart on Thu Jun 21 11:07:11 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ewatercycle's documentation!
=======================================

The eWaterCycle Python package brings together many components from the
`eWaterCycle project <https://www.ewatercycle.org/>`_. An overall goal of this
project is to make hydrological modelling fully reproducible, open, and FAIR.

Modelled after `PyMT <https://pymt.readthedocs.io/en/latest/index.html>`_, it
enables interactively running a model from a Python environment like so:

.. code-block:: python

   from ewatercycle.models import Wflow
   model = Wflow(version="2020.1.1", parameterset=example_parameter_set, forcing=example_forcing)
   cfg_file, cfg_dir = model.setup(end_time="2020-01-01T00:00:00Z")
   model.initialize(cfg_file)

   output = []
   while model.time < model.end_time:
       model.update()
       discharge = model.get_value_at_coords("RiverRunoff", lat=[52.3], lon=[5.2])
       output.append(discharge)

To learn how to use the package, see the `User guide <user_guide.html>`_ and
`example pages <examples.html>`_.

Typically the eWaterCycle platform is deployed on a system that can be accessed
through the browser via JupyterHub, and comes preconfigured with readily
available parameter sets, meteorological forcing data, model images, etcetera.
This makes it possible for researchers to quickly run an experiment without the
hassle of installing a model or creating suitable input data. To learn more
about the system setup, read our `System setup <system_setup.html>`_ page.

In general eWaterCycle tries to strike a balance between making it easy to use
standard available elements of an experiment (datasets, models, analysis
algorithms), and supplying custom elements. This does mean that a simple usecase
sometimes requires slightly more lines of code than strictly nescessary, for the
sake of making it easy to adapt this code to more complex and/or custom
usecases.


Glossary
--------

To avoid miscommunication, here we define explicitly what we mean by some terms
that are commonly used throughout this documentation.

- **Experiment**: A notebook running one or more hydrological models and producing a scientific result.
- **Model**: Software implementation of an algorithm. Note this excludes data required for this model.
- **Forcing**: all time dependent data needed to run a model, and that is not impacted by the model.
- **Model Parameters**: fixed parameters (depth of river, land use, irrigation channels, dams). Considered constant during a model run.
- **Parameter Set**: File based collection of parameters for a certain model, resolution, and possibly area.
- **Model instance**: single running instance of a model, including all data required, and with a current state.


.. toctree::
  :maxdepth: 2
  :hidden:

  user_guide
  system_setup
  adding_models
  examples
  hpc_to_cluster
  observations
  API Reference <apidocs/ewatercycle.rst>
Observations
============

The eWaterCycle platform supports observations relevant for calibrating and validating models. We currently support USGS and GRDC river discharge observations.

USGS
----

The `U.S. Geological Survey Water Services <https://waterservices.usgs.gov/>`_ provides public discharge data for a large number of US based stations. In eWaterCycle we make use of the `USGS web service <https://waterservices.usgs.gov/rest/IV-Service.html>`_ to automatically retrieve this data.
The Discharge timestamp is corrected to the UTC timezone. Units are converted from cubic feet per second to cubic meter per second.

GRDC
----

The `Global Runoff Data Centre <https://www.bafg.de/GRDC/EN/Home/homepage_node.html>`_ provides discharge data for a large number of stations around the world. In eWaterCycle we support GRDC data. This is not downloaded automatically, but required to be present on the infrastructure where the eWaterCycle platform is deployed. By special permission from GRDC our own instance contains data from the ArcticHYCOS and GCOS/GTN-H, GTN-R projects.
