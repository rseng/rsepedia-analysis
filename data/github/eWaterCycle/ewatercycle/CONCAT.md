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
