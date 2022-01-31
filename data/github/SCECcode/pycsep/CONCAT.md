# v0.5.2 (01/25/2021)
## Change-log
Fixed failing build from matplotlib 3.5.0 release (#162)  
Updates to documentation and tutorials (#165)  
Added theory of tests to documentation (#171)  
Added notebooks folder for community recipes (#173)  
Pin numpy version to 1.25.1 to fix (#168) 

## Credits
William Savran (@wsavran)
Kirsty Bayliss (@kirstybayliss)

# v0.5.1 (11/15/2021)

## Change-log
Modified plot_spatial_dataset and plot_catalog to correctly handle Global plots (#150)  
Updated plot_customizations example in the docs (#150)  
Added DOI badge and conda downloads badge (#156)  
Added plotting args to catalog evaluation plots (#154)  
Add option to ignore colorbar in spatial dataset plots (#154)

## Credits
William Savran (@wsavran)  
Pablo Iturrieta (@pabloitu)

# v0.5.0 (11/03/2021)

## Change-log
- Removed normalization of rates on CL-Test (#117)
- Added function to compute bin-wise log-likelihood scores (#118)
- Properly use region to compute spatial counts and spatial magnitude counts in region class (#122)
- Fix for simplified 'fast' lat-lon ratio calculation (#125)
- Add feature to read forecasts with swapped lat/lon values in file (#130)
- Add 'percentile' argument for plotting Poisson evaluations (#131)
- Modify comparison plot to simultaneously plot T and W tests (#132)
- Add feature to trace region outline when plotting spatial data sets (#133)
- Better handling for magnitude ticks in plotting catalogs (#134)
- Refactor polygon to models module (#135)
- Added arguments to modify the fontsize for grid labels and basemap plots (#136)
- Added function to return the midpoints of the valid testing region (#137)
- Include regions when serializing catalogs to JSON (#138)
- Add support for spatial forecasts (#142)
- Upated CI workflows to reduce time required and fix intermittent OOM issues (#145)
- Adds function `get_event_counts` to catalog forecasts (#146)
- Updated README.md (#147)

## Credits
Jose Bayona (@bayonato89)
Marcus Hermann (@mherrmann3)
Pablo Iturrieta (@pabloitu)
Philip Maechling (@pjmaechling)


# v0.4.1 04/14/2021
    - Added 'fast' projection option for plotting spatial datasets (#110)
    - Fixed region border missing when plotted in various projections (#110)
    - Fixed bug where ascii catalog-based forecasts could be incorrectly loaded (#111)

# v0.4.0 03/24/2021 
    - Fixed issue in plot_poisson_consistency_test where one_sided_lower argument not coloring markers correctly
    - Added several plot configurations based on Cartopy 
      - Plotting spatial datasets with ESRI basemap
      - Plotting catalog
      - Plotting regions using outline of polygon
      - Added defaults to forecasts and catalogs
    - Added reader for gzipped UCERF3-ETAS forecasts
    - Updates for INGV readers
    - Fixed bug causing certain events to be placed into incorrrect bins
      
# v0.2 11/11/2020
  Added new catalog formats, support for masked forecast bins, and bug fixes, where applicable PR id are shown in parenthesis.

    - Fixed bug where filtering by catalog by lists did not remove all desired events (#37)
    - Included fast reader for Horus catalog (#39)
    - Modified INGV emrcmt reader (#40)
    - Fixed ndk reader and added unit tests (#44)
    - Fixed issue where magnitues were not correctly bound to gridded-forecast class (#46)
    - Fixed issue where forecasts did not work if lat/lon were not sorted (#47)
    - Fixed minor bug where catalog class did not implement inherited method (#52)
    - Gridded forecasts now parse flag from the ASCII file (#50)
    - Fixed issue where catalog did not filter properly using datetimes (#55)
    


# v0.1 10/08/2020
    Initial release to PyPI and conda-forge

    - Poisson evaluations for gridded forecasts
    - Likelihood-free evaluations for catalog-based forecasts
    - Catalog gridding and filtering
    - Plotting utilities
    - Forecast input and output
    - Documentation at docs.cseptesting.org

# v0.1-dev, 08/16/2018 -- Initial release.
# pyCSEP: Collaboratory for the Study of Earthquake Predictability
![](https://i.postimg.cc/Bb60rVQP/CSEP2-Logo-CMYK.png)
<p align=center>
    <a target="_blank" href="https://python.org" title="Python version"><img src="https://gist.githubusercontent.com/wsavran/efce311162c32460336a4f9892218532/raw/1b9c060efd1c6e52eb53f82d4249107417d6a5ec/pycsep_python_badge.svg">
    <a target="_blank" href="https://pypi.org/project/pycsep"><img src="https://anaconda.org/conda-forge/pycsep/badges/downloads.svg">
    <a target="_blank" href="https://github.com/SCECcode/pycsep/actions"><img src="https://github.com/SCECcode/pycsep/actions/workflows/python-app.yml/badge.svg">
    <a target="_blank" href="https://github.com/SCECcode/pycsep/actions"><img src="https://github.com/SCECcode/pycsep/actions/workflows/build-sphinx.yml/badge.svg">
    <a target="_blank" href="https://codecov.io/gh/SCECcode/pycsep"><img src="https://codecov.io/gh/SCECcode/pycsep/branch/master/graph/badge.svg?token=HTMKM29MAU">
    <a target="_blank" href="https://www.zenodo.org/badge/latestdoi/149362283"><img src="https://www.zenodo.org/badge/149362283.svg" alt="DOI"></a>
</p>

# Description:
The pyCSEP Toolkit helps earthquake forecast model developers evaluate their forecasts with the goal of understanding
earthquake predictability.

pyCSEP should:
1. Help modelers become familiar with formats, procedures, and evaluations used in CSEP Testing Centers.
2. Provide vetted software for model developers to use in their research.
3. Provide quantitative and visual tools to assess earthquake forecast quality.
4. Promote open-science ideas by ensuring transparency and availability of scientific code and results.
5. Curate benchmark models and data sets for modelers to conduct retrospective experiments of their forecasts.

# Table of Contents:
1. [Software Documentation](https://docs.cseptesting.org)
2. [Installation](#installation)
3. [Usage](#usage)
4. [Contributing](#contributing)
5. [Change Log](https://github.com/SCECcode/pycsep/blob/master/CHANGELOG.md)
6. [Credits](#credits)
7. [License](#license)

# Installation:
pyCSEP can be installed in several ways. It can be installed using conda or pip package managers or from the 
source code found in the pyCSEP github repo. Researchers interested in contributing to pyCSEP development should 
install pyCSEP from source code. pyCSEP depends on the following software packages. 
These which may be installed automatically, or manually, depending on the installation method used.
* Python 3.7 or later (https://python.org)
* NumPy 1.21.3 or later (https://numpy.org)
* SciPy 1.7.1 or later (https://scipy.org)
* pandas 1.3.4 or later (https://pandas.pydata.org)
* cartopy 0.20.0 or later (https://scitools.org.uk/cartopy/docs/latest)
* GEOS 3.7.2 or later (https://trac.osgeo.org/geos/)
* PROJ 8.0.0 or later (https://proj.org/)

Please see the [requirements file](https://github.com/SCECcode/pycsep/blob/master/requirements.yml) for a complete list 
of requirements. These are installed automatically when using the `conda` distribution.

Detailed pyCSEP [installation instructions](https://docs.cseptesting.org/getting_started/installing.html) can be found 
in the online pyCSEP documentation.

# Usage: 
Once installed, pyCSEP methods can be invoked from python code by importing package csep. pyCSEP provides objects and 
utilities related to several key concepts:
* Earthquake Catalogs
* Earthquake Forecasts
* Earthquake Forecast Evaluations
* Regions

An simple example to download and plot an earthquake catalog from the USGS ComCat:
<pre>
import csep
from csep.core import regions
from csep.utils import time_utils
start_time = time_utils.strptime_to_utc_datetime('2019-01-01 00:00:00.0')
end_time = time_utils.utc_now_datetime()
catalog = csep.query_comcat(start_time, end_time)
catalog.plot(show=True)
</pre>

Please see [pyCSEP Getting Started](https://docs.cseptesting.org/getting_started/core_concepts) documentation for more examples and tutorials.

# Software Support:
Software support for pyCSEP is provided by that Southern California Earthquake Center (SCEC) Research Computing Group. 
This group supports several research software distributions including UCVM. Users can report issues and feature requests 
using the pyCSEP github-based issue tracking link below. Developers will also respond to emails sent to the SCEC software contact listed below.
1. [pyCSEP Issues](https://github.com/SCECcode/pycsep/issues)
2. Email Contact: software [at] scec [dot] usc [dot] edu

# Contributing:
We welcome contributions to the pyCSEP Toolkit.  If you would like to contribute to this package, including software, tests, and documentation, 
please visit the [contribution guidelines](https://github.com/SCECcode/pycsep/blob/master/CONTRIBUTING.md) for guidelines on how to contribute to pyCSEP development.
pyCSEP contributors agree to abide by the code of conduct found in our [Code of Conduct](CODE_OF_CONDUCT.md) guidelines.

# Credits:
Development of pyCSEP is a group effort. A list of developers that have contributed to the PyCSEP Toolkit 
are listed in the [credits](CREDITS.md) file in this repository.

# License:
The pyCSEP software is distributed under the BSD 3-Clause open-source license. Please see the [license](LICENSE.txt) file for more information.
PyCSEP development is a group effort.

Contributors:
* William Savran, Southern California Earthquake Center
* Pablo Iturrieta, GFZ Potsdam
* Khawaja Asim, GFZ Potsdam
* Han Bao, University of California, Los Angeles
* Kirsty Bayliss, University of Edinburgh
* Jose Bayona, University of Bristol
* Thomas Beutin, GFZ Potsdam
* Marcus Hermann, University of Naples 'Frederico II'
* Edric Pauk, Southern California Earthquake Center
* Max Werner, University of Bristol
* Danijel Schorlemmner, GFZ Potsdam
* Philip Maechling, Southern California Earthquake Center

Thanks to everyone for all your contributions! 
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to make participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, sex characteristics, gender identity and expression,
level of experience, education, socio-economic status, nationality, personal
appearance, race, religion, or sexual identity and orientation.

## Our Standards

Examples of behavior that contributes to creating a positive environment
include:

* Using welcoming and inclusive language
* Being respectful of differing viewpoints and experiences
* Gracefully accepting constructive criticism
* Focusing on what is best for the community
* Showing empathy towards other community members

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery and unwelcome sexual attention or
  advances
* Trolling, insulting/derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or electronic
  address, without explicit permission
* Other conduct which could reasonably be considered inappropriate in a
  professional setting

## Our Responsibilities

Project maintainers are responsible for clarifying the standards of acceptable
behavior and are expected to take appropriate and fair corrective action in
response to any instances of unacceptable behavior.

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

## Scope

This Code of Conduct applies within all project spaces, and it also applies when
an individual is representing the project or its community in public spaces.
Examples of representing a project or community include using an official
project e-mail address, posting via an official social media account, or acting
as an appointed representative at an online or offline event. Representation of
a project may be further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at software@scec.org. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org/), version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq# Contributing to pyCSEP

This document provides an overview on how to contribute to pyCSEP. It will provide step-by-step instructions and hope to 
answer some questions.


## Getting Started

* Make sure you have an active GitHub account
* Download and install `git`
* Read git documentaion if you aren't familiar with `git`
* Install the **development version** of PyCSEP
* If you haven't worked with git Forks before, make sure to read the documentation linked here:
[some helping info](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/working-with-forks).

## Developer Installation

This shows you how to install a copy of the repository that you can use to create Pull Requests and sync with the upstream
repository. First, fork the repo on GitHub. It will now live at `https://github.com/<YOUR_GITHUB_USERNAME>/pycsep`. 
We recommend using `conda` to install the development environment.

    git clone https://github.com/<YOUR_GITHUB_USERNAME>/pycsep.git
    cd pycsep
    conda env create -f requirements.yml
    conda activate csep-dev
    pip install -e .[all]
    # add upstream repository
    git remote add upstream https://github.com/SCECCode/pycsep.git
    
Note: use the command `conda deactivate` to go back to your regular environment when you are done working with pyCSEP.

## Submitting a Pull Request

### Some notes for starting a pull request

Pull requests are how we use your changes to the code! Please submit them to us! Here's how:

1. Make a new branch. For features/additions base your new branch at `master`.
2. Make sure to add tests! Only pull requests for documentation, refactoring, or plotting features do not require a test.
3. Also, documentation must accompany new feature requests.
    - Note: We would really appreciate pull requests that help us improve documentation.
3. Make sure the tests pass. Run `./run_tests.sh` in the top-level directory of the repo.
4. Push your changes to your fork and submit a pull request. Make sure to set the branch to `pycsep:master`.
5. Wait for our review. There may be some suggested changes or improvements. Changes can be made after
the pull request has been opening by simply adding more commits to your branch.

Pull requests can be changed after they are opened, so create a pull request as early as possible.
This allows us to provide feedback during development and to answer any questions.

Also, if you find pyCSEP to be useful, but don't want to contribute to the code we highly encourage updates to the documentation!

Please make sure to set the correct branch for your pull request. Also, please do not include large files in your pull request.
If you feel that you need to add large files, such as a benchmark forecast, let us know and we can figure something out.

### Tips to get your pull request accepted quickly

1. Any new feature that contains calculations must contain unit-tests to ensure that the calculations are doing what you
expect. Some exceptions to this are documentation changes and new plotting features. 
2. Documentation should accompany any new feature additions into the package.
    * Plotting functions should provide a sphinx-gallery example, which can be found [here](https://github.com/SCECcode/pycsep/blob/master/examples/tutorials/catalog_filtering.py).
    * More complex features might require additional documentation. We will let you know upon seeing your pull request.
    * The documentation use sphinx which compiles reST. Some notes on that can be found [here](https://www.sphinx-doc.org/en/master/usage/quickstart.html).
3. pyCSEP uses pytest as a test runner. Add new tests to the `tests` folder in an existing file or new file starting matching `test_*.py`
4. New scientific capabilities that are not previously published should be presented to the CSEP science group as part of a 
science review. This will consist of a presentation that provides a scientific justification for the feature.
5. Code should follow the [pep8](https://pep8.org/) style-guide.
6. Functions should use [Google style docstrings](https://www.sphinx-doc.org/en/master/usage/extensions/example_google.html). These 
get compiled by Sphinx to become part of the documentation.
 
## Submitting an Issue

Please open an issue if you want to ask a question about PyCSEP.

* Please search through the past issues to see if your question or the bug has already been addressed
* Please apply the correct tag to your issue so others can search

If you want to submit a bug report, please provide the information below:
* pyCSEP version, Python version, and Platform (Linux, Windows, Mac OSX, etc)
* How did you install pyCSEP (pip, anaconda, from source...)
* Please provide a short, complete, and correct example that demonstrates the issue.
* If this broke in a recent update, please tell us when it used to work.

## Additional Resources
* [Working with Git Forks](https://help.github.com/en/github/collaborating-with-issues-and-pull-requests/working-with-forks)
* [Style Guide](http://google.github.io/styleguide/pyguide.html)
* [Docs or it doesn’t exist](https://lukeplant.me.uk/blog/posts/docs-or-it-doesnt-exist/)
* [Quickstart guide for Sphinx](https://www.sphinx-doc.org/en/master/usage/quickstart.html)
* [Pep8 style guide](https://pep8.org/)
* Performance Tips:
  * [Python](https://wiki.python.org/moin/PythonSpeed/PerformanceTips)
  * [NumPy and ctypes](https://scipy-cookbook.readthedocs.io/)
  * [SciPy](https://www.scipy.org/docs.html)
  * [NumPy Book](http://csc.ucdavis.edu/~chaos/courses/nlp/Software/NumPyBook.pdf)
# pyCSEP Pull Request Checklist

Please check out the [contributing guidelines](https://github.com/SCECcode/pycsep/blob/master/CONTRIBUTING.md) for some tips 
on making pull requests to pyCSEP. 

Fixes issue #(*please fill in or delete if not needed*).

## Type of change:

Please delete options that are not relevant.

- [ ] Bug fix (non-breaking change which fixes an issue)
- [ ] New feature (non-breaking change which adds functionality)
- [ ] Breaking change (fix or feature that would cause existing functionality to not work as expected)
- [ ] Documentation update (this pull request adds no new code)
- [ ] Unpublished science feature (This may require a science review)
- [ ] This change requires a documentation update

## Checklist:

- [ ] I have performed a self-review of my own code
- [ ] I have commented my code, particularly in hard-to-understand areas
- [ ] I have made corresponding changes to the documentation
- [ ] My changes generate no new warnings
- [ ] I have added tests that prove my fix is effective or that my feature works
- [ ] New and existing unit tests pass locally with my changes
# Helpful notebooks

This folder contains notebooks that might be helpful for providing information about pyCSEP and providing recipes for common
processing tasks. These are meant to supplement the [tutorials](https://docs.cseptesting.org/tutorials/).

1. Runnable notebook showing behavior of CSEP tests
pyCSEP: Tools for Earthquake Forecast Developers
================================================

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Getting Started

    getting_started/installing
    getting_started/core_concepts
    getting_started/theory

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Tutorials and Examples

    tutorials/catalog_filtering.rst
    tutorials/plot_gridded_forecast.rst
    tutorials/gridded_forecast_evaluation.rst
    tutorials/working_with_catalog_forecasts.rst
    tutorials/catalog_forecast_evaluation.rst
    tutorials/plot_customizations.rst

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: User Guide

    concepts/catalogs
    concepts/forecasts
    concepts/evaluations
    concepts/regions

.. toctree::
    :maxdepth: 2
    :hidden:
    :caption: Help & Reference

    reference/glossary
    reference/publications
    reference/roadmap
    reference/developer_notes
    reference/api_reference


*PyCSEP tools help earthquake forecast model developers evaluate their forecasts and provide the machinery to implement
experiments within CSEP testing centers.*

About
-----
The Collaboratory for the Study of Earthquake Predictability (CSEP) supports an international effort to conduct earthquake
forecasting experiments. CSEP supports these activities by developing the cyberinfrastructure necessary to run earthquake
forecasting experiments including the statistical framework required to evaluate probabilistic earthquake forecasts.

PyCSEP is a python library that provides tools for (1) evaluating probabilistic earthquake forecasts, (2) working
with earthquake catalogs in this context, and (3) creating visualizations. Official experiments that run in CSEP testing centers
will be implemented using the code provided by this package.

Project Goals
-------------
1. Help modelers become familiar with formats, procedures, and evaluations used in CSEP Testing Centers.
2. Provide vetted software for model developers to use in their research.
3. Provide quantitative and visual tools to assess earthquake forecast quality.
4. Promote open-science ideas by ensuring transparency and availability of scientific code and results.
5. Curate benchmark models and data sets for modelers to conduct retrospective experiments of their forecasts.

Contributing
------------
We highly encourage users of this package to get involved in the development process. Any contribution is helpful, even
suggestions on how to improve the package, or additions to the documentation (those are particularly welcome!). Check out
the `Contribution guidelines <https://github.com/SCECCode/csep2/blob/dev/CONTRIBUTING.md>`_ for a step by step on how to contribute to the project. If there are
any questions, please contact us!

Contacting Us
-------------
* For general discussion and bug reports please post issues on the `PyCSEP GitHub <https://github.com/SCECCode/csep2>`_.
* This project adheres to a `Code of Conduct <https://github.com/SCECCode/csep2/blob/dev/CODE_OF_CONDUCT.md>`_. By participating you agree to follow its terms.

List of Contributors
--------------------
* (Lead) William Savran ``wsavran [at] usc [dot] edu``
* Muhammad Asim Khawaja
* Thomas Beutin
* Pablo Iturrieta




.. _forecast-reference:

#########
Forecasts
#########

PyCSEP supports two types of earthquake forecasts that can be evaluated using the tools provided in this package.

1. Grid-based forecasts
2. Catalog-based forecasts

These forecast types and the PyCSEP objects used to represent them will be explained in detail in this document.

.. contents:: Table of Contents
    :local:
    :depth: 2

*****************
Gridded forecasts
*****************

Grid-based forecasts assume that earthquakes occur in independent and discrete space-time-magnitude bins. The occurrence
of these earthquakes are described only by their expected rates. This forecast format provides a general representation
of seismicity that can accommodate forecasts without explicit likelihood functions, such as those created using smoothed
seismicity models. Gridded forecasts can also be produced using simulation-based approaches like
epidemic-type aftershock sequence models.

Currently, grid-based forecasts define their spatial component using a 2D Cartesian (rectangular) grid, and
their magnitude bins using a 1D Cartesian (rectangular) grid. The last bin (largest magnitude) bin is assumed to
continue until infinity. Forecasts use latitude and longitude to define the bin edge of the spatial grid. Typical values
for the are 0.1° x 0.1° (lat x lon) and 0.1 ΔMw units. These choices are not strictly enforced and can defined
according the specifications of an experiment.

Working with gridded forecasts
##############################

PyCSEP provides the :class:`GriddedForecast<csep.core.forecasts.GriddedForecast>` class to handle working with
grid-based forecasts. Please see visit :ref:`this example<grid-forecast-evaluation>` for an end-to-end tutorial on
how to evaluate a grid-based earthquake forecast.

.. autosummary:: csep.core.forecasts.GriddedForecast

Default file format
--------------------

The default file format of a gridded-forecast is a tab delimited ASCII file with the following columns
(names are not included): ::

    LON_0 	LON_1 	LAT_0 	LAT_1 	DEPTH_0 DEPTH_1 MAG_0 	MAG_1 	RATE					FLAG
    -125.4	-125.3	40.1	40.2	0.0     30.0	4.95	5.05	5.8499099999999998e-04	1

Each row represents a single space-magnitude bin and the entire forecast file contains the rate for a specified
time-horizon. An example of a gridded forecast for the RELM testing region can be found
`here <https://github.com/SCECcode/csep2/blob/dev/csep/artifacts/ExampleForecasts/GriddedForecasts/helmstetter_et_al.hkj.aftershock-fromXML.dat>`_.


The coordinates (LON, LAT, DEPTH, MAG) describe the independent space-magnitude region of the forecast. The lower
coordinates are inclusive and the upper coordinates are exclusive. Rates are incremental within the magnitude range
defined by [MAG_0, MAG_1). The FLAG is a legacy value from CSEP testing centers that indicates whether a spatial cell should
be considered by the forecast. Currently, the implementation does not allow for individual space-magnitude cells to be
flagged. Thus, if a spatial cell is flagged then all corresponding magnitude cells are flagged.

.. note::
    PyCSEP only supports regions that have a thickness of one layer. In the future, we plan to support more complex regions
    including those that are defined using multiple depth regions. Multiple depth layers can be collapsed into a single
    layer by summing. This operations does reduce the resolution of the forecast.

Custom file format
------------------

The :meth:`GriddedForecast.from_custom<csep.core.forecasts.GriddedForecast.from_custom>` method allows you to provide
a function that can read custom formats. This can be helpful, because writing this function might be required to convert
the forecast into the appropriate format in the first place. This function has no requirements except that it returns the
expected data.

.. automethod:: csep.core.forecasts.GriddedForecast.from_custom


***********************
Catalog-based forecasts
***********************

Catalog-based earthquake forecasts are issued as collections of synthetic earthquake catalogs. Every synthetic catalog
represents a realization of the forecast that is representative the uncertainty present in the model that generated
the forecast. Unlike grid-based forecasts, catalog-based forecasts retain the space-magnitude dependency of the events
they are trying to model. A grid-based forecast can be easily computed from a catalog-based forecast by assuming a
space-magnitude region and counting events within each bin from each catalog in the forecast. There can be issues with
under sampling, especially for larger magnitude events.

Working with catalog-based forecasts
####################################

.. autosummary:: csep.core.forecasts.CatalogForecast

Please see visit :ref:`this<catalog-forecast-evaluation>` example for an end-to-end tutorial on how to evaluate a catalog-based
earthquake forecast. An example of a catalog-based forecast stored in the default PyCSEP format can be found
`here<https://github.com/SCECcode/pycsep/blob/dev/csep/artifacts/ExampleForecasts/CatalogForecasts/ucerf3-landers_1992-06-28T11-57-34-14.csv>`_.


The standard format for catalog-based forecasts a comma separated value ASCII format. This format was chosen to be
human-readable and easy to implement in all programming languages. Information about the format is shown below.

.. note::
    Custom formats can be supported by writing a custom function or sub-classing the
    :ref:`AbstractBaseCatalog<csep.core.forecasts.AbstractBaseCatalog>`.

The event format matches the follow specfication: ::

    LON, LAT, MAG, ORIGIN_TIME, DEPTH, CATALOG_ID, EVENT_ID
    -125.4, 40.1, 3.96, 1992-01-05T0:40:3.1, 8, 0, 0

Each row in the catalog corresponds to an event. The catalogs are expected to be placed into the same file and are
differentiated through their `catalog_id`. Catalogs with no events can be handled in a couple different ways intended to
save storage.

The events within a catalog should be sorted in time, and the *catalog_id* should be increasing sequentially. Breaks in
the *catalog_id* are interpreted as missing catalogs.

The following two examples show how you represent a forecast with 5 catalogs each containing zero events.

**1. Including all events (verbose)** ::

    LON, LAT, MAG, ORIGIN_TIME, DEPTH, CATALOG_ID, EVENT_ID
    ,,,,,0,
    ,,,,,1,
    ,,,,,2,
    ,,,,,3,
    ,,,,,4,

**2. Short-hand** ::

    LON, LAT, MAG, ORIGIN_TIME, DEPTH, CATALOG_ID, EVENT_ID
    ,,,,,4,

The following three example show how you could represent a forecast with 5 catalogs. Four of the catalogs contain zero events
and one catalog contains one event.

**3. Including all events (verbose)** ::

    LON, LAT, MAG, ORIGIN_TIME, DEPTH, CATALOG_ID, EVENT_ID
    ,,,,,0,
    ,,,,,1,
    ,,,,,2,
    ,,,,,3,
    -125.4, 40.1, 3.96, 1992-01-05T0:40:3.1, 8, 4, 0

**4. Short-hand** ::

    LON, LAT, MAG, ORIGIN_TIME, DEPTH, CATALOG_ID, EVENT_ID
    -125.4, 40.1, 3.96, 1992-01-05T0:40:3.1, 8, 4, 0

The simplest way to orient the file follow (3) in the case where some catalogs contain zero events. The zero oriented
catalog_id should be assigned to correspond with the total number of catalogs in the forecast. In the case where every catalog
contains zero forecasted events, you would specify the forecasting using (2). The *catalog_id* should be assigned to
correspond with the total number of catalogs in the forecast.

.. _plots-reference:

#####
Plots
#####

PyCSEP provides several functions to produce commonly used plots, such as an earthquake forecast or the evaluation catalog
or perhaps a combination of the two.

.. contents:: Table of Contents
    :local:
    :depth: 2

************
Introduction
************



**************
Plot arguments
**************

***************
Available plots
***************

.. _catalogs-reference:

########
Catalogs
########

PyCSEP provides routines for working with and manipulating earthquake catalogs for the purposes of evaluating earthquake
forecasting models.

If you are able to make use of these tools for other reasons, please let us know. We are especially interested in
including basic catalog statistics into this package. If you are interested in helping implement routines like
b-value estimation and catalog completeness that would be much appreciated.

.. contents:: Table of Contents
    :local:
    :depth: 2

************
Introduction
************

PyCSEP catalog basics
=====================

An earthquake catalog contains a collection of seismic events each defined by a set of attributes. PyCSEP implements
a simple event description that is suitable for evaluating earthquake forecasts. In this format, every seismic event is
defined by its location in space (longitude, latitude, and depth), magnitude, and origin time. In addition, each event
can have an optional ``event_id`` as a unique identifier.

PyCSEP provides :class:`csep.core.catalogs.CSEPCatalog` to represent an earthquake catalog. The essential event data are stored in a
`structured NumPy array <https://numpy.org/doc/stable/user/basics.rec.html>`_ with the following data type. ::

    dtype = numpy.dtype([('id', 'S256'),
                         ('origin_time', '<i8'),
                         ('latitude', '<f4'),
                         ('longitude', '<f4'),
                         ('depth', '<f4'),
                         ('magnitude', '<f4')])

Additional information can be associated with an event using the ``id`` field in the structured array in a class member called
``metadata``. The essential event data must be complete meaning that each event should have these attributes defined. The metadata
is more freeform and PyCSEP does not impose any restrictions on the way that event metadata is stored. Only that the metadata
for an event should be accessible using the event ``id``. An example of this could be ::

    catalog = csep.load_catalog('catalog_file.csv')
    event_metadata = catalog.metadata[event_id]

This would load a catalog stored in the PyCSEP .csv format. PyCSEP contains catalog readers for the following formats
(We are also looking to support other catalog formats. Please suggest some or better yet help us write the readers!):

1. CSEP ascii format
2. NDK format (used by the gCMT catalog)
3. INGV gCMT catalog
4. ZMAP format
5. pre-processed JMA format

PyCSEP supports the ability to easily define a custom reader function for a catalog format type that we don't currently support.
If you happen to implement a reader for a new catalog format please check out the
`contribution guidelines <https://github.com/SCECCode/csep2/blob/dev/CONTRIBUTING.md>`_ and make a pull request so
we can include this in the next release.

Catalog as Pandas dataframes
============================

You might be comfortable using Pandas dataframes to manipulate tabular data. PyCSEP provides some routines for accessing
catalogs as a :class:`pandas.DataFrame`. You can use ``df = catalog.to_dataframe(with_datetimes=True)`` to return the
DataFrame representation of the catalog. Using the ``catalog = CSEPCatalog.from_dataframe(df)`` you can return back to the
PyCSEP data model.

.. note::

    Going between a DataFrame and CSEPCatalog is a lossy transformation. It essentially only retains the essential event
    attributes that are defined by the ``dtype`` of the class.

****************
Loading catalogs
****************

Load catalogs from files
========================

You can easily load catalogs in the supported format above using :func:`csep.load_catalog`. This function provides a
top-level function to load catalogs that are currently supported by PyCSEP. You must specify the type of the catalog and
the format you want it to be loaded. The type of the catalog can be: ::

    catalog_type = ('ucerf3', 'csep-csv', 'zmap', 'jma-csv', 'ndk')
    catalog_format = ('csep', 'native')

The catalog type determines which reader :mod:`csep.utils.readers` will be used to load in the file. The default is the
``csep-csv`` type and the ``native`` format. The ``jma-csv`` format can be created using the ``./bin/deck2csv.pl``
Perl script.

.. note::
    The format is important for ``ucerf3`` catalogs, because those are stored as big endian binary numbers by default.
    If you are working with ``ucerf3-etas`` catalogs and would like to convert them into the CSEPCatalog format you can
    use the ``format='csep'`` option when loading in a catalog or catalogs.

Load catalogs from ComCat
=========================

PyCSEP provides top-level functions to load catalogs using ComCat. We incorporated the work done by Mike Hearne and
others from the U.S. Geological Survey into PyCSEP in an effort to reduce the dependencies of this project. The top-level
access to ComCat catalogs can be accessed from :func:`csep.query_comcat`. Some lower level functionality can be accessed
through the :mod:`csep.utils.comcat` module. All credit for this code goes to the U.S. Geological Survey.

:ref:`Here<tutorial-catalog-filtering>` is complete example of accessing the ComCat catalog.

Writing custom loader functions
===============================

You can easily add custom loader functions to import data formats that are not currently included with the PyCSEP tools.
Both :meth:`csep.core.catalogs.CSEPCatalog.load_catalog` and :func:`csep.load_catalog` support an optional argument
called ``loader`` to support these custom data formats.

In the simplest form the function should have the following stub: ::

    def my_custom_loader_function(filename):
        """ Custom loader function for catalog data.

        Args:
            filename (str): path to the file containing the path to the forecast

        Returns:
            eventlist: iterable of event data with the order:
                    (event_id, origin_time, latitude, longitude, depth, magnitude)
        """

        # imagine there is some logic to read in data from filename

        return eventlist

This function can then be passed to :func:`csep.load_catalog` or :meth:`CSEPCatalog.load_catalog<csep.core.catalogs.CSEPCatalog.load_catalog>`
with the ``loader`` keyword argument. The function should be passed as a first-class object like this: ::

    import csep
    my_custom_catalog = csep.load_catalog(filename, loader=my_custom_loader_function)

.. note::
    The origin_time is actually an integer time. We recommend to parse the timing information as a
    :class:`datetime.datetime` object and use the :func:`datetime_to_utc_epoch<csep.utils.time_utils.datetime_to_utc_epoch>`
    function to convert this to an integer time.

Notice, we did not actually call the function but we just passed it as a reference. These functions can also access
web-based catalogs like we implement with the :func:`csep.query_comcat` function. This function doesn't work with either
:func:`csep.load_catalog` or :meth:`CSEPCatalog.load_catalog<csep.core.catalogs.CSEPCatalog.load_catalog>`,
because these are intended for file-based catalogs. Instead, we can create the catalog object directly.
We would do that like this ::

    def my_custom_web_loader(...):
        """ Accesses catalog from online data source.

        There are no requirements on the arguments if you are creating the catalog directly from the class.

        Returns:
            eventlist: iterable of event data with the order:
                (event_id, origin_time, latitude, longitude, depth, magnitude)
        """

        # custom logic to access online catalog repository

        return eventlist

As you might notice, all loader functions are required to return an event-list. This event-list must be iterable and
contain the required event data.

.. note::
    The events in the eventlist should follow the form ::

        eventlist = my_custom_loader_function(...)

        event = eventlist[0]

        event[0] = event_id
        # see note above about using integer times
        event[1] = origin_time
        event[2] = latitude
        event[3] = longitude
        event[4] = depth
        event[5] = magnitude


Once you have a function that returns an eventlist, you can create the catalog object directly. This uses the
:class:`csep.core.catalogs.CSEPCatalog` as an example. ::

    import csep

    eventlist = my_custom_web_loader(...)
    catalog = csep.catalogs.CSEPCatalog(data=eventlist, **kwargs)

The **kwargs represents any other keyword argument that can be passed to
:class:`CSEPCatalog<csep.core.catalogs.CSEPCatalog>`. This could be the ``catalog_id`` or the
:class:`CartesianGrid2D<csep.core.regions.CartesianGrid2D>`.

Including custom event metadata
===============================

Catalogs can include additional metadata associated with each event. Right now, there are no direct applications for
event metadata. Nonetheless, it can be included with a catalog object.

The event metadata should be a dictionary where the keys are the ``event_id`` of the individual events. For example, ::

    event_id = 'my_dummy_id'
    metadata_dict = catalog.metadata[event_id]

Each event meta_data should be a JSON-serializable dictionary or a class that implements the to_dict() and from_dict() methods. This is
required to properly save the catalog files into JSON format and verify whether two catalogs are the same. You can see
the :meth:`to_dict<csep.core.catalogs.CSEPCatalog>` and :meth:`from_dict<csep.core.catalogs.CSEPCatalog>` methods for an
example of how these would work.

***************************
Accessing Event Information
***************************

In order to utilize the low-level acceleration from Numpy, most catalog operations are vectorized. The catalog classes
provide some getter methods to access the essential catalog data. These return arrays of :class:`numpy.ndarray` with the
``dtype`` defined by the class.

.. automodule:: csep.core.catalogs

The following functions return :class:`numpy.ndarrays<numpy.ndarray>` of the catalog information.

.. autosummary::

   CSEPCatalog.event_count
   CSEPCatalog.get_magnitudes
   CSEPCatalog.get_longitudes
   CSEPCatalog.get_latitudes
   CSEPCatalog.get_depths
   CSEPCatalog.get_epoch_times
   CSEPCatalog.get_datetimes
   CSEPCatalog.get_cumulative_number_of_events

The catalog data can be iterated through event-by-event using a standard for-loop. For example, we can do something
like ::

    for event in catalog.data:
        print(
            event['id'],
            event['origin_time'],
            event['latitude'],
            event['longitude'],
            event['depth'],
            event['magnitude']
        )

The keyword for the event tuple are defined by the ``dtype`` of the class. The keywords for
:class:`CSEPCatalog<csep.core.catalogs.CSEPCatalog>` are shown in the snippet directly above. For example, a quick and
dirty plot of the cumulative events over time can be made using the :mod:`matplotlib.pyplot` interface ::

    import csep
    import matplotlib.pyplot as plt

    # lets assume we already loaded in some catalog
    catalog = csep.load_catalog("my_catalog_path.csv")

    # quick and dirty plot
    fig, ax = plt.subplots()
    ax.plot(catalog.get_epoch_times(), catalog.get_cumulative_number_of_events())
    plt.show()

****************
Filtering events
****************

Most of the catalog files (or catalogs accessed via the web) contain more events that are desired for a given use case.
PyCSEP provides a few routines to help filter events out of the catalog. The following methods help to filter out
unwanted events from the catalog.

.. autosummary::

    CSEPCatalog.filter
    CSEPCatalog.filter_spatial
    CSEPCatalog.apply_mct

Filtering events by attribute
=============================

The function :meth:`CSEPCatalog.filter<csep.core.catalogs.CSEPCatalog.filter>` provides the ability to filter events
based on their essential attributes. This function works by parsing filtering strings and applying them using a logical
`and` operation. The catalog strings have the following format ``filter_string = f"{attribute} {operator} {value}"``. The
filter strings represent a statement that would evaluate as `True` after they are applied. For example, the statement
``catalog.filter('magnitude >= 2.5')`` would retain all events in the catalog greater-than-or-equal to magnitude 2.5.

The attributes are determined by the dtype of the catalog, therefore you can filter based on the ``(origin_time, latitude,
longitude, depth, and magnitude).`` Additionally, you can use the attribute ``datetime`` and provide a :class:`datetime.datetime`
object to filter events using the data type.

The filter function can accept a string or a list of filter statements. If the function is called without any arguments
the function looks to use the ``catalog.filters`` member. This can be provided during class instantiation or bound
to the class afterward. :ref:`Here<tutorial-catalog-filtering>` is complete example of how to filter a catalog
using the filtering strings.

Filtering events in space
=========================

You might want to supply a non-rectangular polygon that can be used to filter events in space. This is commonly done
to prepare an observed catalog for forecast evaluation. Right now, this can be accomplished by supplying a
:class:`region<csep.core.regions.CartesianGrid2D>` to the catalog or
:meth:`filter_spatial<csep.core.catalogs.CSEPCatalog.filter_spatial>`. There will be more information about using regions
in the :ref:`user-guide page<regions-reference>` page. The :ref:`catalog filtering<tutorial-catalog-filtering>` contains
a complete example of how to filter a catalog using a user defined aftershock region based on the M7.1 Ridgecrest
mainshock.

Time-dependent magnitude of completeness
========================================

Seismic networks have difficulty recording events immediately after a large event occurs, because the passing seismic waves
from the larger event become mixed with any potential smaller events. Usually when we evaluate an aftershock forecast, we should
account for this time-dependent magnitude of completeness. PyCSEP provides the
:ref:`Helmstetter et al., [2006]<helmstetter-2006>` implementation of the time-dependent magnitude completeness model.

This requires information about an event which can be supplied directly to :meth:`apply_mct<csep.core.catalogs.CSEPCatalog>`.
Additionally, PyCSEP provides access to the ComCat API using :func:`get_event_by_id<csep.utils.comcat.get_event_by_id>`.
An exmaple of this can be seen in the :ref:`filtering catalog tutorial<tutorial-catalog-filtering>`.


**************
Binning Events
**************

Another common task requires binning earthquakes by their spatial locations and magnitudes. This is routinely done when
evaluating earthquake forecasts. Like filtering a catalog in space, you need to provide some information about the region
that will be used for the binning. Please see the :ref:`user-guide page<regions-reference>` for more information about
regions.

.. note::
    We would like to make this functionality more user friendly. If you have suggestions or struggles, please open an issue
    on the GitHub page and we'd be happy to incorporate these ideas into the toolkit.

The following functions allow binning of catalogs using space-magnitude regions.

.. autosummary::

   CSEPCatalog.spatial_counts
   CSEPCatalog.magnitude_counts
   CSEPCatalog.spatial_magnitude_counts

These functions return :class:`numpy.ndarrays<numpy.ndarray>` containing the count of the events determined from the
catalogs. This example shows how to obtain magnitude counts from a catalog.
The index of the ndarray corresponds to the index of the associated space-magnitude region. For example, ::

    import csep
    import numpy

    catalog = csep.load_catalog("my_catalog_file")

    # returns bin edges [2.5, 2.6, ... , 7.5]
    bin_edges = numpy.arange(2.5, 7.55, 0.1)

    magnitude_counts = catalog.magnitude_counts(mag_bins=bin_edges)

In this example, ``magnitude_counts[0]`` is the number of events with 2.5 ≤ M < 2.6. All of the magnitude binning assumes
that the final bin extends to infinity, therefore ``magnitude_counts[-1]`` contains the number of events with
7.5 ≤ M < ∞... _regions-reference

#######
Regions
#######

.. automodule:: csep.utils.basic_types

PyCSEP includes commonly used CSEP testing regions and classes that facilitate working with gridded data sets. This
module is early in development and will be a focus of future development.

.. contents:: Table of Contents
    :local:
    :depth: 2

.. :currentmodule:: csep

.. automodule:: csep.core.regions

Practically speaking, earthquake forecasts, especially time-dependent forecasts, treat time differently than space and
magnitude. If we consider a family of monthly forecasts for the state of California for earthquakes with **M** 3.95+,
each of these forecasts would use the same space-magnitude region, even though the time periods are
different. Because the time horizon is an implicit property of the forecast, we do not explicitly consider time in the region
objects provided by PyCSEP. This module contains tools for working with gridded regions in both space and magnitude.

First, we will describe how the spatial regions are handled. Followed by magnitude regions, and how these two aspects
interact with one another.

**************
Region objects
**************

Currently, PyCSEP provides the :class:`CartesianGrid2D<csep.core.regions.CartesianGrid2D>` to handle binning catalogs
and defining regions for earthquake forecasting evaluations. We plan to expand this module in the future to include
more complex spatial regions.

2D Cartesian grids
##################

This section contains information about using 2D cartesian grids.

.. autosummary::

    CartesianGrid2D

.. note::
    We are planning to do some improvements to this module and to expand its capabilities. For example, we would like to
    handle non-regular grids such as a quad-tree. Also, a single Polygon should be able to act as the spatial component
    of the region. These additions will make this toolkit more useful for crafting bespoke experiments and for general
    catalog analysis. Feature requests are always welcome!

The :class:`CartesianGrid2D<csep.core.regions.CartesianGrid2D>` acts as a data structure that can associate a spatial
location (eg., lon and lat) with its corresponding spatial bin. This class is optimized to work with regular grids,
although they do not have to be complete (they can have holes) and they do not have to be rectangular (each row / column
can have a different starting coordinate).

The :class:`CartesianGrid2D<csep.core.regions.CartesianGrid2D>` maintains a list of
:class:`Polygon<csep.core.regions.Polygon>` objects that represent the individual spatial bins from the overall
region. The origin of each polygon is considered to be the lower-left corner (the minimum latitude and minimum longitude).

.. autosummary::

    CartesianGrid2D.num_nodes
    CartesianGrid2D.get_index_of
    CartesianGrid2D.get_location_of
    CartesianGrid2D.get_masked
    CartesianGrid2D.get_cartesian
    CartesianGrid2D.get_bbox
    CartesianGrid2D.midpoints
    CartesianGrid2D.origins
    CartesianGrid2D.from_origins


Creating spatial regions
########################

Here, we describe how the class works starting with the class constructors. ::

    @classmethod
    def from_origins(cls, origins, dh=None, magnitudes=None, name=None):
        """ Convenience function to create CartesianGrid2D from list of polygon origins """

For most applications, using the :meth:`from_origins<csep.core.regions.CartesianGrid2D.from_origins>` function will be
the easiest way to create a new spatial region. The method accepts a 2D :class:`numpy.ndarray` containing the x (lon) and y (lat)
origins of the spatial bin polygons. These should be the complete set of origins. The function will attempt to compute the
grid spacing by comparing the x and y values between adjacent origins. If this does not seem like a reliable approach
for your region, you can explicitly provide the grid spacing (dh) to this method.

When a :class:`CartesianGrid2D<csep.core.regions.CartesianGrid2D>` is created the following steps occur:

    1. Compute the bounding box containing all polygons (2D array)
    2. Create a map between the index of the 2D bounding box and the list of polygons of the region.
    3. Store a boolean flag indicating whether a given cell in the 2D array is valid or not

Once these mapping have been created, we can now associate an arbitrary (lon, lat) point with a spatial cell using the
mapping defined in (2). The :meth:`get_index_of<csep.core.regions.CartesianGrid2D.get_index_of>` accepts a list
of longitudes and latitudes and returns the index of the polygon they are associated with. For instance, this index can
now be used to access a data value stored in another data structure.

***************
Testing Regions
***************

CSEP has defined testing regions that can be used for earthquake forecasting experiments. The following functions in the
:mod:`csep.core.regions` module returns a :class:`CartesianGrid2D<csep.core.regions.CartesianGrid2D>` consistent with
these regions.

.. autosummary::

    california_relm_region
    italy_csep_region
    global_region

****************
Region Utilities
****************

PyCSEP also provides some utilities that can facilitate working with regions. As we expand this module, we will include
functions to accommodate different use-cases.

.. autosummary::

    magnitude_bins
    create_space_magnitude_region
    parse_csep_template
    increase_grid_resolution
    masked_region
    generate_aftershock_region.. _evaluation-reference:

.. automodule:: csep.core.poisson_evaluations

###########
Evaluations
###########

PyCSEP provides routines to evaluate both gridded and catalog-based earthquake forecasts. This page explains how to use
the forecast evaluation routines and also how to build "mock" forecast and catalog classes to accommodate different
custom forecasts and catalogs.

.. contents:: Table of Contents
    :local:
    :depth: 2

.. :currentmodule:: csep

****************************
Gridded-forecast evaluations
****************************

Grid-based earthquake forecasts assume earthquakes occur in discrete space-time-magnitude bins and their rate-of-
occurrence can be defined using a single number in each magnitude bin. Each space-time-magnitude bin is assumed to be
an independent Poisson random variable. Therefore, we use likelihood-based evaluation metrics to compare these
forecasts against observations.

PyCSEP provides two groups of evaluation metrics for grid-based earthquake forecasts. The first are known as
consistency tests and they verify whether a forecast in consistent with an observation. The second are comparative tests
that can be used to compare the performance of two (or more) competing forecasts.
PyCSEP implements the following evaluation routines for grid-based forecasts. These functions are intended to work with
:class:`GriddedForecasts<csep.core.forecasts.GriddedForecast>` and :class:`CSEPCatalogs`<csep.core.catalogs.CSEPCatalog>`.
Visit the :ref:`catalogs reference<catalogs-reference>` and the :ref:`forecasts reference<forecasts-reference>` to learn
more about to import your forecasts and catalogs into PyCSEP.

.. note::
    Grid-based forecast evaluations act directly on the forecasts and catalogs as they are supplied to the function.
    Any filtering of catalogs and/or scaling of forecasts must be done before calling the function.
    This must be done before evaluating the forecast and should be done consistently between all forecasts that are being
    compared.

See the :ref:`example<grid-forecast-evaluation>` for gridded forecast evaluation for an end-to-end walkthrough on how
to evaluate a gridded earthquake forecast.


Consistency tests
=================

.. autosummary::

   number_test
   magnitude_test
   spatial_test
   likelihood_test
   conditional_likelihood_test


Comparative tests
=================

.. autosummary::

   paired_t_test
   w_test

Publication references
======================

1. Number test (:ref:`Schorlemmer et al., 2007<schorlemmer-2007>`; :ref:`Zechar et al., 2010<zechar-2010>`)
2. Magnitude test (:ref:`Zechar et al., 2010<zechar-2010>`)
3. Spatial test (:ref:`Zechar et al., 2010<zechar-2010>`)
4. Likelihood test (:ref:`Schorlemmer et al., 2007<schorlemmer-2007>`; :ref:`Zechar et al., 2010<zechar-2010>`)
5. Conditional likelihood test (:ref:`Werner et al., 2011<werner-2011>`)
6. Paired t test (:ref:`Rhoades et al., 2011<rhoades-2011>`)
7. Wilcoxon signed-rank test (:ref:`Rhoades et al., 2011<rhoades-2011>`)

**********************************
Catalog-based forecast evaluations
**********************************

Catalog-based forecasts are issued as a family of stochastic event sets (synthetic earthquake catalogs) and can express
the full uncertainty of the forecasting model. Additionally, these forecasts retain the inter-event dependencies that
are lost when using discrete space-time-magnitude grids. This problem can impact the evaluation performance of
time-dependent forecasts like the epidemic type aftershock sequence model (ETAS).

In order to support generative or simulator-based models, we define a suite of consistency tests that compare forecasted
distributions against observations without the use of a parametric likelihood function. These evaluations take advantage
of the fact that the forecast and the observations are both earthquake catalogs. Therefore, we can compute identical
statistics from these catalogs and compare them against one another.

We provide four statistics that probe fundamental aspects of the earthquake forecasts. Please see
:ref:`Savran et al., 2020<savran-2020>` for a complete description of the individual tests. For the implementation
details please follow the links below and see the :ref:`example<catalog-forecast-evaluation>` for catalog-based
forecast evaluation for an end-to-end walk through.

.. automodule:: csep.core.catalog_evaluations

Consistency tests
=================

.. autosummary::

   number_test
   spatial_test
   magnitude_test
   pseudolikelihood_test
   calibration_test

Publication reference
=====================

1. Number test (:ref:`Savran et al., 2020<savran-2020>`)
2. Spatial test (:ref:`Savran et al., 2020<savran-2020>`)
3. Magnitude test (:ref:`Savran et al., 2020<savran-2020>`)
4. Pseudolikelihood test (:ref:`Savran et al., 2020<savran-2020>`)
5. Calibration test (:ref:`Savran et al., 2020<savran-2020>`)

****************************
Preparing evaluation catalog
****************************

The evaluations in PyCSEP do not implicitly filter the observed catalogs or modify the forecast data when called. For most
cases, the observation catalog should be filtered according to:
    1. Magnitude range of the forecast
    2. Spatial region of the forecast
    3. Start and end-time of the forecast

Once the observed catalog is filtered so it is consistent in space, time, and magnitude as the forecast, it can be used
to evaluate a forecast. A single evaluation catalog can be used to evaluate multiple forecasts so long as they all cover
the same space, time, and magnitude region.

*********************
Building mock classes
*********************

Python is a duck-typed language which means that it doesn't care what the object type is only that it has the methods or
functions that are expected when that object is used. This can come in handy if you want to use the evaluation methods, but
do not have a forecast that completely fits with the forecast classes (or catalog classes) provided by PyCSEP.

.. note::
    Something about great power and great responsibility... For the most reliable results, write a loader function that
    can ingest your forecast into the model provided by PyCSEP. Mock-classes can work, but should only be used in certain
    circumstances. In particular, they are very useful for writing software tests or to prototype features that can
    be added into the package.

This section will walk you through how to compare two forecasts using the :func:`paired_t_test<csep.core.poisson_evaluations>`
with mock forecast and catalog classes. This sounds much more complex than it really is, and it gives you the flexibility
to use your own formats and interact with the tools provided by PyCSEP.

.. warning::

    The simulation-based Poisson tests (magnitude_test, likelihood_test, conditional_likelihood_test, and spatial_test)
    are optimized to work with forecasts that contain equal-sized spatial bins. If your forecast uses variable sized spatial
    bins you will get incorrect results. If you are working with forecasts that have variable spatial bins, create an
    issue on GitHub because we'd like to implement this feature into the toolkit and we'd love your help.

If we look at the :func:`paired_t_test<csep.core.poisson_evaluations>` we see that it has the following code ::

    def paired_t_test(gridded_forecast1, gridded_forecast2, observed_catalog, alpha=0.05, scale=False):
        """ Computes the t-test for gridded earthquake forecasts.

        Args:
            gridded_forecast_1 (csep.core.forecasts.GriddedForecast): nd-array storing gridded rates, axis=-1 should be the magnitude column
            gridded_forecast_2 (csep.core.forecasts.GriddedForecast): nd-array storing gridded rates, axis=-1 should be the magnitude column
            observed_catalog (csep.core.catalogs.AbstractBaseCatalog): number of observed earthquakes, should be whole number and >= zero.
            alpha (float): tolerance level for the type-i error rate of the statistical test
            scale (bool): if true, scale forecasted rates down to a single day

        Returns:
            evaluation_result: csep.core.evaluations.EvaluationResult
        """

        # needs some pre-processing to put the forecasts in the context that is required for the t-test. this is different
        # for cumulative forecasts (eg, multiple time-horizons) and static file-based forecasts.
        target_event_rate_forecast1, n_fore1 = gridded_forecast1.target_event_rates(observed_catalog, scale=scale)
        target_event_rate_forecast2, n_fore2 = gridded_forecast2.target_event_rates(observed_catalog, scale=scale)

        # call the primative version operating on ndarray
        out = _t_test_ndarray(target_event_rate_forecast1, target_event_rate_forecast2, observed_catalog.event_count, n_fore1, n_fore2,
                              alpha=alpha)

        # prepare evaluation result object
        result = EvaluationResult()
        result.name = 'Paired T-Test'
        result.test_distribution = (out['ig_lower'], out['ig_upper'])
        result.observed_statistic = out['information_gain']
        result.quantile = (out['t_statistic'], out['t_critical'])
        result.sim_name = (gridded_forecast1.name, gridded_forecast2.name)
        result.obs_name = observed_catalog.name
        result.status = 'normal'
        result.min_mw = numpy.min(gridded_forecast1.magnitudes)

Notice that the function expects two forecast objects and one catalog object. The ``paired_t_test`` function calls a
method on the forecast objects named :meth:`target_event_rates<csep.core.forecasts.GriddedForecast.target_event_rates>`
that returns a tuple (:class:`numpy.ndarray`, float) consisting of the target event rates and the expected number of events
from the forecast.

.. note::
    The target event rate is the expected rate for an observed event in the observed catalog assuming that
    the forecast is true. For a simple example, if we forecast a rate of 0.3 events per year in some bin of a forecast,
    each event that occurs within that bin has a target event rate of 0.3 events per year. The expected number of events
    in the forecast can be determined by summing over all bins in the gridded forecast.

We can also see that the ``paired_t_test`` function uses the ``gridded_forecast1.name`` and calls the :func:`numpy.min`
on the ``gridded_forecast1.magnitudes``. Using this information, we can create a mock-class that implements these methods
that can be used by this function.

.. warning::
    If you are creating mock-classes to use with evaluation functions, make sure that you visit the corresponding
    documentation and source-code to make sure that your methods return values that are expected by the function. In
    this case, it expects the tuple (target_event_rates, expected_forecast_count). This will not always be the case.
    If you need help, please create an issue on the GitHub page.

Here we show an implementation of a mock forecast class that can work with the
:func:`paired_t_test<csep.core.poisson_evaluations.paired_t_test>` function. ::

    class MockForecast:

        def __init__(self, data=None, name='my-mock-forecast', magnitudes=(4.95)):

            # data is not necessary, but might be helpful for implementing target_event_rates(...)
            self.data = data
            self.name = name
            # this should be an array or list. it can be as simple as the default argument.
            self.magnitudes = magnitudes

        def target_event_rates(catalog, scale=None):
            """ Notice we added the dummy argument scale. This function stub should match what is called paired_t_test """

            # whatever custom logic you need to return these target event rates given your catalog can go here
            # of course, this should work with whatever catalog you decide to pass into this function

            # this returns the tuple that paired_t_test expects
            return (ndarray_of_target_event_rates, expected_number_of_events)

You'll notice that :func:`paired_t_test<csep.core.poisson_evaluations.paired_t_test>` expects a catalog class. Looking back
at the function definition we can see that it needs ``observed_catalog.event_count`` and ``observed_catalog.name``. Therefore
the mock class for the catalog would look something like this ::

    class MockCatalog:

        def __init__(self, event_count, data=None, name='my-mock-catalog'):

            # this is not necessary, but adding data might be helpful for implementing the
            # logic needed for the target_event_rates(...) function in the MockForecast class.
            self.data = data
            self.name = name
            self.event_count = event_count


Now using these two objects you can call the :func:`paired_t_test<csep.core.poisson_evaluations.paired_t_test>` directly
without having to modify any of the source code. ::

    # create your forecasts
    mock_forecast_1 = MockForecast(some_forecast_data1)
    mock_forecast_2 = MockForecast(some_forecast_data2)

    # lets assume that catalog_data is an array that contains the catalog data
    catalog = MockCatalog(len(catalog_data))

    # call the function using your classes
    eval_result = paired_t_test(mock_forecast_1, mock_forecast_2, catalog)

The only requirement for this approach is that you implement the methods on the class that the calling function expects.
You can add anything else that you need in order to make those functions work properly. This example is about
as simple as it gets.

.. note::

    If you want to use mock-forecasts and mock-catalogs for other evaluations. You can just add the additional methods
    that are needed onto the mock classes you have already built.
API Reference
=============

This contains a reference document to the PyCSEP API.

.. automodule:: csep

.. :currentmodule:: csep

Loading catalogs and forecasts
------------------------------

.. autosummary::
   :toctree: generated

   load_stochastic_event_sets
   load_catalog
   query_comcat
   load_gridded_forecast
   load_catalog_forecast

Catalogs
--------

.. :currentmodule:: csep.core.catalogs

.. automodule:: csep.core.catalogs


Catalog operations are defined using :class:`AbstractBaseCatalog` class.

.. autosummary::
   :toctree: generated

   AbstractBaseCatalog
   CSEPCatalog
   UCERF3Catalog

Catalog operations
------------------

Input and output operations for catalogs:

.. autosummary::
   :toctree: generated

   CSEPCatalog.to_dict
   CSEPCatalog.from_dict
   CSEPCatalog.to_dataframe
   CSEPCatalog.from_dataframe
   CSEPCatalog.write_json
   CSEPCatalog.load_json
   CSEPCatalog.load_catalog
   CSEPCatalog.write_ascii
   CSEPCatalog.load_ascii_catalogs
   CSEPCatalog.get_csep_format
   CSEPCatalog.plot

Accessing event information:

.. autosummary::
   :toctree: generated

   CSEPCatalog.event_count
   CSEPCatalog.get_magnitudes
   CSEPCatalog.get_longitudes
   CSEPCatalog.get_latitudes
   CSEPCatalog.get_depths
   CSEPCatalog.get_epoch_times
   CSEPCatalog.get_datetimes
   CSEPCatalog.get_cumulative_number_of_events

Filtering and binning:

.. autosummary::
   :toctree: generated

   CSEPCatalog.filter
   CSEPCatalog.filter_spatial
   CSEPCatalog.apply_mct
   CSEPCatalog.spatial_counts
   CSEPCatalog.magnitude_counts
   CSEPCatalog.spatial_magnitude_counts

Other utilities:

.. autosummary::
   :toctree: generated

   CSEPCatalog.update_catalog_stats
   CSEPCatalog.length_in_seconds
   CSEPCatalog.get_bvalue

.. currentmodule:: csep.core.forecasts
.. automodule:: csep.core.forecasts

Forecasts
---------

PyCSEP provides classes to interact with catalog and grid based Forecasts

.. autosummary::
   :toctree: generated

   GriddedForecast
   CatalogForecast

Gridded forecast methods:

.. autosummary::
   :toctree: generated

   GriddedForecast.data
   GriddedForecast.event_count
   GriddedForecast.sum
   GriddedForecast.magnitudes
   GriddedForecast.min_magnitude
   GriddedForecast.magnitude_counts
   GriddedForecast.spatial_counts
   GriddedForecast.get_latitudes
   GriddedForecast.get_longitudes
   GriddedForecast.get_magnitudes
   GriddedForecast.get_index_of
   GriddedForecast.get_magnitude_index
   GriddedForecast.load_ascii
   GriddedForecast.from_custom
   GriddedForecast.get_rates
   GriddedForecast.target_event_rates
   GriddedForecast.scale_to_test_date
   GriddedForecast.plot

Catalog forecast methods:

.. autosummary::
   :toctree: generated

   CatalogForecast.magnitudes
   CatalogForecast.min_magnitude
   CatalogForecast.spatial_counts
   CatalogForecast.magnitude_counts
   CatalogForecast.get_expected_rates
   CatalogForecast.get_dataframe
   CatalogForecast.write_ascii
   CatalogForecast.load_ascii

.. automodule:: csep.core.catalog_evaluations

Evaluations
-----------

PyCSEP provides implementations of evaluations for both catalog-based forecasts and grid-based forecasts.

Catalog-based forecast evaluations:

.. autosummary::
   :toctree: generated

   number_test
   spatial_test
   magnitude_test
   pseudolikelihood_test
   calibration_test

.. automodule:: csep.core.poisson_evaluations

Grid-based forecast evaluations:

.. autosummary::
   :toctree: generated

   number_test
   magnitude_test
   spatial_test
   likelihood_test
   conditional_likelihood_test
   paired_t_test
   w_test

.. automodule:: csep.core.regions

Regions
-------

PyCSEP includes commonly used CSEP testing regions and classes that facilitate working with gridded data sets. This
module is early in development and help is welcome here!

Region class(es):

.. autosummary::
    :toctree: generated

    CartesianGrid2D

Testing regions:

.. autosummary::
    :toctree: generated

    california_relm_region
    italy_csep_region
    global_region

Region utilities:

.. autosummary::
    :toctree: generated

    magnitude_bins
    create_space_magnitude_region
    parse_csep_template
    increase_grid_resolution
    masked_region
    generate_aftershock_region
    california_relm_region


Plotting
--------

.. automodule:: csep.utils.plots

General plotting:

.. autosummary::
   :toctree: generated

   plot_histogram
   plot_ecdf
   plot_basemap
   plot_spatial_dataset
   add_labels_for_publication

Plotting from catalogs:

.. autosummary::
   :toctree: generated

   plot_magnitude_versus_time
   plot_catalog

Plotting stochastic event sets and evaluations:

.. autosummary::
   :toctree: generated

   plot_cumulative_events_versus_time
   plot_magnitude_histogram
   plot_number_test
   plot_magnitude_test
   plot_distribution_test
   plot_likelihood_test
   plot_spatial_test
   plot_calibration_test

Plotting gridded forecasts and evaluations:

.. autosummary::
   :toctree: generated

   plot_spatial_dataset
   plot_comparison_test
   plot_poisson_consistency_test

.. automodule:: csep.utils.time_utils

Time Utilities
--------------

.. autosummary::
   :toctree: generated

   epoch_time_to_utc_datetime
   datetime_to_utc_epoch
   millis_to_days
   days_to_millis
   strptime_to_utc_epoch
   timedelta_from_years
   strptime_to_utc_datetime
   utc_now_datetime
   utc_now_epoch
   create_utc_datetime
   decimal_year

.. automodule:: csep.utils.comcat

Comcat Access
-------------

We integrated the code developed by Mike Hearne and others at the USGS to reduce the dependencies of this package. We plan
to move this to an external and optional dependency in the future.

.. autosummary::
   :toctree: generated

   search
   get_event_by_id

.. automodule:: csep.utils.calc

Calculation Utilities
---------------------

.. autosummary::
   :toctree: generated

   nearest_index
   find_nearest
   func_inverse
   discretize
   bin1d_vec

.. automodule:: csep.utils.stats

Statistics Utilities
--------------------

.. autosummary::
   :toctree: generated

   sup_dist
   sup_dist_na
   cumulative_square_diff
   binned_ecdf
   ecdf
   greater_equal_ecdf
   less_equal_ecdf
   min_or_none
   max_or_none
   get_quantiles
   poisson_log_likelihood
   poisson_joint_log_likelihood_ndarray
   poisson_inverse_cdf

.. automodule:: csep.utils.basic_types

Basic types
-----------

.. autosummary::
    :toctree: generated

    AdaptiveHistogram=====================
Terms and Definitions
=====================

The page contains terms and their definitions (and possible mathematical definitions) that are commonly used throughout the documentation
and CSEP literature.

.. contents:: Table of Contents
    :local:
    :depth: 2


.. _earthquake-catalog:

Earthquake catalog
------------------
List of earthquakes (either tectonic or non-tectonic) defined through their location in space, origin time of the event, and
their magnitude.


.. _earthquake_forecast:

Earthquake Forecast
-------------------
A probabilistic statement about the occurrence of seismicity that can include information about the magnitude and spatial
location. CSEP supports earthquake forecasts expressed as the expected rate of seismicity in disjoint space and magnitude bins
and as families of synthetic earthquake catalogs.

.. _stochastic-event-set:

Stochastic Event Set
--------------------
Collection of synthetic earthquakes (events) that are produced by an earthquake.
A *stochastic event set* consists of *N* events that represent a continuous representation of seismicity that can sample
the uncertainty present within in the forecasting model.

.. _time-dependent-forecast:

Time-dependent Forecast
-----------------------
The forecast changes over time using new information not available at the time the forecast was issued. For example,
epidemic-type aftershock sequence models (ETAS) models can utilize updated using newly observed seismicity to produce
new forecasts consistent with the model.

.. _time-independent-forecast:

Time-independent Forecast
-------------------------
The forecast does not change with time. Time-independent forecasts are generally used for long-term forecasts
needed for probabalistic seismic hazard analysis.
#######################
Referenced Publications
#######################



.. _helmstetter-2006:

Helmstetter, A., Y. Y. Kagan, and D. D. Jackson (2006). Comparison of short-term and time-independent earthquake
forecast models for southern California, *Bulletin of the Seismological Society of America* **96** 90-106.

.. _rhoades-2011:

Rhoades, D. A., D. Schorlemmer, M. C. Gerstenberger, A. Christophersen, J. D. Zechar, and M. Imoto (2011). Efficient
testing of earthquake forecasting models, *Acta Geophys* **59** 728-747.

.. _savran-2020:

Savran, W., M. J. Werner, W. Marzocchi, D. Rhoades, D. D. Jackson, K. R. Milner, E. H. Field, and A. J. Michael (2020).
Pseudoprospective evaluation of UCERF3-ETAS forecasts during the 2019 Ridgecrest Sequence,
*Bulletin of the Seismological Society of America*.

.. _schorlemmer-2007:

Schorlemmer, D., M. Gerstenberger, S. Wiemer, D. D. Jackson, and D. A. Rhoades (2007). Earthquake likelihood model
testing, *Seismological Research Letters* **78** 17-29.

.. _werner-2011:

Werner, M. J., A. Helmstetter, D. D. Jackson, and Y. Y. Kagan (2011). High-Resolution Long-Term and Short-Term
Earthquake Forecasts for California, *Bulletin of the Seismological Society of America* **101** 1630-1648.

.. _zechar-2010:

Zechar, J. D., M. C. Gerstenberger, and D. A. Rhoades (2010). Likelihood-Based Tests for Evaluating Space-Rate-Magnitude
Earthquake Forecasts, *Bulletin of the Seismological Society of America* **100** 1184-1195.






.. _roadmap:

###################
Development Roadmap
###################

This page contains expected changes for new releases of `pyCSEP`.
Last updated 3 November 2021.

v0.6.0
======

1. Include receiver operating characteristic (ROC) curve
2. Kagan I1 score
3. Add function to plot spatial log-likelihood scores
4. Add documentation section to explain maths of CSEP tests



Developer Notes
===============

Last updated: 25 January 2022

Creating a new release of pyCSEP
--------------------------------

These are the steps required to create a new release of pyCSEP. This requires a combination of updates to the repository
and Github. You will need to build the wheels for distribution on PyPI and upload them to GitHub to issue a release.
The final step involves uploading the tar-ball of the release to PyPI. CI tools provided by `conda-forge` will automatically
bump the version on `conda-forge`. Note: permissions are required to push new versions to PyPI.

1. Code changes
***************
1. Bump the version number in `_version.py <https://github.com/SCECcode/pycsep/tree/master/csep/_version.py>`_
2. Update `codemeta.json <https://github.com/SCECcode/pycsep/blob/master/codemeta.json>`_
3. Update `CHANGELOG.md <https://github.com/SCECcode/pycsep/blob/master/CHANGELOG.md>`_. Include links to Github pull requests if possible.
4. Update `CREDITS.md <https://github.com/SCECcode/pycsep/blob/master/CREDITS.md>`_ if required.
5. Update the version in `conf.py <https://github.com/SCECcode/pycsep/blob/master/docs/conf.py>`_.
6. Issue a pull request that contains these changes.
7. Merge pull request when all changes are merged into `master` and versions are correct.

2. Creating source distribution
*******************************

Issue these commands from the top-level directory of the project::

    python setup.py check

If that executes with no warnings or failures build the source distribution using the command::

    python setup.py sdist

This creates a folder called `dist` that contains a file called `pycsep-X.Y.Z.tar.gz`. This is the distribution
that will be uploaded to `PyPI`, `conda-forge`, and Github.

Upload to PyPI using `twine`. This requires permissions to push to the PyPI account::

    twine upload dist/pycsep-X.Y.Z.tar.gz

3. Create release on Github
***************************
1. Create a new `release <https://github.com/SCECcode/pycsep/releases>`_ on GitHub. This can be saved as a draft until it's ready.
2. Copy new updates information from `CHANGELOG.md <https://github.com/SCECcode/pycsep/blob/master/CHANGELOG.md>`_.
3. Upload tar-ball created from `setup.py`.
4. Publish release.Installing pyCSEP
=================

We are working on a ``conda-forge`` recipe and PyPI distribution.
If you plan on contributing to this package, visit the
`contribution guidelines <https://github.com/SCECcode/pycsep/blob/master/CONTRIBUTING.md>`_ for installation instructions.

.. note:: This package requires >=Python 3.7.

The easiest way to install PyCSEP is using ``conda``. It can also be installed using ``pip`` or built from source.

Using Conda
-----------
For most users, you can use ::

    conda install --channel conda-forge pycsep

Using Pip
---------

Before this installation will work, you must **first** install the following system dependencies. The remaining dependencies
should be installed by the installation script. To help manage dependency issues, we recommend using virtual environments
like `virtualenv`.

| Python 3.7 or later (https://python.org)
|
| NumPy 1.10 or later (https://numpy.org)
|     Python package for scientific computing and numerical calculations.
|
| GEOS 3.3.3 or later (https://trac.osgeo.org/geos/)
|     C++ library for processing geometry.
|
| PROJ 4.9.0 or later (https://proj4.org/)
|     Library for cartographic projections.

Example for Ubuntu: ::

    sudo apt-get install libproj-dev proj-data proj-bin
    sudo apt-get install libgeos-dev
    pip install --upgrade pip
    pip install numpy

Example for MacOS: ::

    brew install proj geos
    pip install --upgrade pip
    pip install numpy

Installing from Source
----------------------

Use this approach if you want the most up-to-date code. This creates an editable installation that can be synced with
the latest GitHub commit.

We recommend using virtual environments when installing python packages from source to avoid any dependency conflicts. We prefer
``conda`` as the package manager over ``pip``, because ``conda`` does a good job of handling binary distributions of packages
across multiple platforms. Also, we recommend using the ``miniconda`` installer, because it is lightweight and only includes
necessary pacakages like ``pip`` and ``zlib``.

Using Conda
***********

If you don't have ``conda`` on your machine, download and install `Miniconda <https://docs.conda.io/en/latest/miniconda.html>`_. ::

    git clone https://github.com/SCECcode/pycsep
    cd pycsep
    conda env create -f requirements.yml
    conda activate csep-dev
    # Installs in editor mode with all dependencies
    pip install -e .

Note: If you want to go back to your default environment use the command ``conda deactivate``.

Using Pip / Virtualenv
**********************

We highly recommend using Conda, because this tools helps to manage binary dependencies on Python packages. If you
must use `Virtualenv <https://packaging.python.org/guides/installing-using-pip-and-virtual-environments/>`_
follow these instructions: ::

    git clone https://github.com/SCECcode/pycsep
    cd pycsep
    python -m virtualenv venv
    source venv/bin/activate
    # Installs in editor mode dependencies are installed by conda
    pip install -e .[all]

Note: If you want to go back to your default environment use the command ``deactivate``.

Developers Installation
-----------------------

This shows you how to install a copy of the repository that you can use to create Pull Requests and sync with the upstream
repository. First, fork the repo on GitHub. It will now live at ``https://github.com/<YOUR_GITHUB_USERNAME>/pycsep``.
We recommend using ``conda`` to install the development environment. ::

    git clone https://github.com/<YOUR_GITHUB_USERNAME>/pycsep.git
    cd pycsep
    conda env create -f requirements.yml
    conda activate csep-dev
    pip install -e .
    # Allow sync with default repository
    git remote add upstream https://github.com/SCECCode/pycsep.git

Now you can pull from upstream using ``git pull upstream master`` to keep your copy of the repository in sync with the
latest commits.Theory of CSEP Tests
====================

This page describes the theory of each of the forecast tests
included in pyCSEP along with working code examples. You will find
information on the goals of each test, the theory behind the tests, how
the tests are applied in practice, and how forecasts are ‘scored’ given
the test results. Also, we include the code required to run in each test
and a description of how to interpret the test results.

.. code:: ipython3

    import csep
    from csep.core import (
        regions,
        catalog_evaluations,
        poisson_evaluations as poisson
    )
    from csep.utils import (
        datasets,
        time_utils,
        comcat,
        plots,
        readers
    )

    # Filters matplotlib warnings
    import warnings
    warnings.filterwarnings('ignore')

Grid-based Forecast Tests
-------------------------

These tests are designed for grid-based forecasts (e.g., Schorlemmer et
al., 2007), where expected rates are provided in discrete Poisson
space-magnitude cells covering the region of interest. The region
:math:`\boldsymbol{R}` is then the product of the spatial rate
:math:`\boldsymbol{S}` and the binned magnitude rate
:math:`\boldsymbol{M}`,

.. math::  \boldsymbol{R} = \boldsymbol{M} \times \boldsymbol{S}.

A forecast :math:`\boldsymbol{\Lambda}` can be fully specified as the
expected number of events (or rate) in each space-magnitude bin
(:math:`m_i, s_j`) covering the region :math:`\boldsymbol{R}` and
therefore can be written as

.. math::  \boldsymbol{\Lambda} = \{ \lambda_{m_i, s_j}| m_i \in \boldsymbol{M}, s_j \in \boldsymbol{S} \},

where :math:`\lambda_{m_i, s_j}` is the expected rate of events in
magnitude bin :math:`m_i` and spatial bin :math:`s_j`. The observed
catalogue of events :math:`\boldsymbol{\Omega}` we use to evaluate the
forecast is similarly discretised into the same space-magnitude bins,
and can be described as

.. math::  \boldsymbol{\Omega} = \{ \omega_{m_i, s_j}| m_i \in \boldsymbol{M}, s_j \in \boldsymbol{S} \},

where :math:`\omega_{m_i, s_j}` is the observed number of
events in spatial cell :math:`s_j` and magnitude bin :math:`m_i`. The
magnitude bins are specified in the forecast: typically these are in 0.1
increments and this is the case in the examples we use here. These
examples use the Helmstetter et al (2007) smoothed seismicity forecast
(including aftershocks), testing over a 5 year period between 2010 and
2015.

.. code:: ipython3

    # Set up experiment parameters
    start_date = time_utils.strptime_to_utc_datetime('2010-01-01 00:00:00.0')
    end_date = time_utils.strptime_to_utc_datetime('2015-01-01 00:00:00.0')

    # Loads from the PyCSEP package
    helmstetter = csep.load_gridded_forecast(
        datasets.helmstetter_aftershock_fname,
        start_date=start_date,
        end_date=end_date,
        name='helmstetter_aftershock'
    )

    # Set up evaluation catalog
    catalog = csep.query_comcat(helmstetter.start_time, helmstetter.end_time,
                                min_magnitude=helmstetter.min_magnitude)

    # Filter evaluation catalog
    catalog = catalog.filter_spatial(helmstetter.region)

    # Add seed for reproducibility in simulations
    seed = 123456

    # Number of simulations for Poisson consistency tests
    nsim = 100000


.. parsed-literal::

    Fetched ComCat catalog in 5.9399449825286865 seconds.

    Downloaded catalog from ComCat with following parameters
    Start Date: 2010-01-10 00:27:39.320000+00:00
    End Date: 2014-08-24 10:20:44.070000+00:00
    Min Latitude: 31.9788333 and Max Latitude: 41.1431667
    Min Longitude: -125.3308333 and Max Longitude: -115.0481667
    Min Magnitude: 4.96
    Found 24 events in the ComCat catalog.


Consistency tests
~~~~~~~~~~~~~~~~~

The consistency tests evaluate the consistency of a forecast against
observed earthquakes. These tests were developed across a range of
experiments and publications (Schorlemmer et al, 2007; Zechar et al
2010; Werner et al, 2011a). The consistency tests are based on the
likelihood of observing the catalogue (actual recorded events) given the
forecast. Since the space-magnitude bins are assumed to be independent,
the joint-likelihood of observing the events in each individual bin
given the specified forecast can be written as

.. math::  Pr(\omega_1 | \lambda_1) Pr(\omega_2 | \lambda_2)...Pr(\omega_n | \lambda_n) = \prod_{m_i , s_j \in \boldsymbol{R}} f_{m_i, s_j}(\omega(m_i, s_j)),

where :math:`f_{m_i, s_j}` specifies the probability distribution in
each space-magnitude bin. We prefer to use the joint log-likelihood in
order to sum log-likelihoods rather than multiply the likelihoods. The
joint log-likelihood can be written as:

.. math::  L(\boldsymbol{\Omega} | \boldsymbol{\Lambda}) = \sum_{m_i , s_j \in \boldsymbol{R}} log(f_{m_i, s_j}(\omega(m_i, s_j)).

The likelihood of the observations, :math:`\boldsymbol{\Omega}`, given
the forecast :math:`\boldsymbol{\Lambda}` is the sum over all
space-magnitude bins of the log probabilities in individual cells of the
forecast. Grid-based forecasts are specified by the expected number of
events in a discrete space-magnitude bin. From the maximum entropy
principle, we assign a Poisson distribution in each bin. In this case,
the probability of an event occurring is independent of the time since
the last event, and events occur at a rate :math:`\lambda`. The
Poissonian joint log-likelihood can be written as

.. math::  L(\boldsymbol{\Omega} | \boldsymbol{\Lambda}) = \sum_{m_i , s_j \in \boldsymbol{R}} -\lambda(m_i, s_j) + \omega(m_i, s_j)\log(\lambda(m_i, s_j)) - log(\omega(m_i, s_j)!),

where :math:`\lambda(m_i, s_j)` and :math:`\omega(m_i, s_j)` are the
expected counts from the forecast and observed counts in cell
:math:`m_i, s_j` respectively. We can calculate the likelihood directly
given the forecast and discretised observations.

Forecast uncertainty

A simulation based approach is used to account for uncertainty in the
forecast. We simulate realizations of catalogs that are consistent with
the forecast to obtain distributions of scores. In the pyCSEP package,
as in the original CSEP tests, simulation is carried out using the
cumulative probability density of the forecast obtained by ordering the
rates in each bin. We shall call :math:`F_{m_is_j}` the cumulative
probability density in cell :math:`(m_i, s_j)`. The simulation approach
then works as follows:

-  For each forecast bin, draw a random number :math:`z` from a uniform
   distribution between 0 and 1
-  Assign this event to a space-magnitude bin through the inverse
   cumulative density distribution at this point
   :math:`F^{-1}_{m_i, s_j}(z)`
-  Iterate over all simulated events to generate a catalog containing
   :math:`N_{sim}` events consistent with the forecast

For each of these tests, we can plot the distribution of likelihoods
computed from theses simulated catalogs relative to the observations
using the ``plots.plot_poisson_consistency_test`` function. We also
calculate a quantile score to diagnose a particular forecast with
repsect. The number of simulations can be supplied to the Poisson
consistency test functions using the ``num_simulations`` argument: for
best results we suggest 100,000 simulations to ensure convergence.

Scoring the tests

Through simulation (as described above), we obtain a set of simulated
catalogs :math:`\{\hat{\boldsymbol{\Omega}}\}`. Each catalogue can be
written as

.. math:: \hat{\boldsymbol{\Omega}}_x =\{ \hat{\lambda}_x(m_i, s_j)|(m_i, s_j) \in \boldsymbol{R}\},

where :math:`\hat{\lambda}_x(m_i, s_j)` is the number of
simulated earthquakes in cell :math:`(m_i, s_j)` of (simulated) catalog
:math:`x` that is consistent with the forecast :math:`\Lambda`. We then
compute the joint log-likelihood for each simulated catalogue
:math:`\hat{L}_x = L(\hat{\Omega}_x|\Lambda)`. The joint log-likelihood
for each simulated catalogue given the forecast gives us a set of
log-likelihoods :math:`\{\hat{\boldsymbol{L}}\}` that represents the
range of log-likelihoods consistent with the forecast. We then compare
our simulated log-likelihoods with the observed log-likelihood
:math:`L_{obs} = L(\boldsymbol{\Omega}|\boldsymbol{\Lambda})` using a
quantile score.

The quantile score is defined by the fraction of simulated joint
log-likelihoods less than or equal to the observed likelihood.

.. math:: \gamma = \frac{ |\{ \hat{L}_x | \hat{L}_x \le L_{obs}\} |}{|\{ \hat{\boldsymbol{L}} \}|}

Whether a forecast can be said to pass an evaluation depends on the
significance level chosen for the testing process. The quantile score
explicitly tells us something about the significance of the result: the
observation is consistent with the forecast with :math:`100(1-\gamma)\%`
confidence (Zechar, 2011). Low :math:`\gamma` values demonstrate that
the observed likelihood score is less than most of the simulated
catalogs. The consistency tests, excluding the N-test, are considered to
be one-sided tests: values which are too small are ruled inconsistent
with the forecast, but very large values may not necessarily be
inconsistent with the forecast and additional testing should be used to
further clarify this (Schorlemmer et al, 2007).

Different CSEP experiments have used different sensitivity values.
Schorlemmer et al (2010b) consider :math:`\gamma \lt 0.05` while the
implementation in the Italian CSEP testing experiment uses
:math:`\gamma` < 0.01 (Taroni et al, 2018). However, the consistency
tests are most useful as diagnostic tools where the quantile score
assesses the level of consistency between observations and data.
Temporal variations in seismicity make it difficult to formally reject a
model from a consistency test over a single evaluation period.

Likelihood-test (L-test)
^^^^^^^^^^^^^^^^^^^^^^^^

Aim: Evaluate the likelihood of observed events given the provided
forecast - this includes the rate, spatial distribution and magnitude
components of the forecast.

Method: The L-test is one of the original forecast tests described in
Schorlemmer et al, 2007. The likelihood of the observation given the
model is described by a Poisson likelihood function in each cell and the
total joint likelihood described by the product over all bins, or the
sum of the log-likelihoods (see above, or Zechar 2011 for more details).

Note: The likelihood scores are dominated by the rate-component of the
forecast. This causes issues in scoring forecasts where the expected
number of events are different from the observed number of events. We
suggest to use the N-test (below) and CL-test (below) independently to
score the rate component, and spatial-magnitude components of the
forecast. This behavior can be observed by comparing the CL-test and
N-test results with the L-test results in this notebook. Since the
forecast overpredicts the rate of events during this testing period, the
L-test provides a passing score even though the space-magnitude and rate
components perform poorly during this evaluation period.

pyCSEP implementation

pyCSEP uses the forecast and catalog and returns the test distribution,
observed statistic and quantile score, which can be accessed from the
``likelihood_test_result`` object. We can pass this directly to the
plotting function, specifying that the test should be one-sided.

.. code:: ipython3

    likelihood_test_result = poisson.likelihood_test(
        helmstetter,
        catalog,
        seed=seed,
        num_simulations=nsim
    )
    ax = plots.plot_poisson_consistency_test(
        likelihood_test_result,
        one_sided_lower=True,
        plot_args={'title': r'$\mathcal{L}-\mathrm{test}$', 'xlabel': 'Log-likelihood'}
    )



.. image:: static/output_6_0.png



pyCSEP plots the resulting :math:`95\%` range of likelihoods returned by
the simulation with the black bar by default. The observed likelihood
score is shown by a green square where the forecast passes the test and
a red circle where the observed likelihood is outside the likelihood
distribution.

CL-test
^^^^^^^

Aim: The original likelihood test described above gives a result that
combines the spatial, magnitude and number components of a forecast. The
conditional likelihood or CL-Test was developed to test the spatial and
magnitude performance of a forecast without the influence of the number
of events (Werner et al. 2011a, 2011b). By conditioning the test
distribution on the observed number of events we elimiate the dependency
with the forecasted number of events as described above.

| Method
| The CL-test is computed in the same way as the L-test, but with the
  number of events normalised to the observed catalog :math:`N_{obs}`
  during the simulation stage. The quantile score is then calculated
  similarly such that

.. math:: \gamma_{CL} = \frac{ |\{ \hat{CL}_x | \hat{CL}_x \le CL_{obs}\} |}{|\{ \hat{\boldsymbol{CL}} \}|}.

Implementation in pyCSEP

.. code:: ipython3

    cond_likelihood_test_result = poisson.conditional_likelihood_test(
        helmstetter,
        catalog,
        seed=seed,
        num_simulations=nsim
    )
    ax = plots.plot_poisson_consistency_test(
        cond_likelihood_test_result,
        one_sided_lower=True,
        plot_args = {'title': r'$CL-\mathrm{test}$', 'xlabel': 'conditional log-likelihood'}
    )



.. image:: static/output_9_0.png


Again, the :math:`95\%` confidence range of likelihoods is shown by the
black bar, and the symbol reflects the observed conditional-likelihood
score. In this case, the observed conditional-likelihood is shown with
the red circle, which falls outside the range of likelihoods simulated
from the forecast. To understand why the L- and CL-tests give different
results, consider the results of the N-test and S-test in the following
sections.

N-test
^^^^^^

Aim: The number or N-test is the most conceptually simple test of a
forecast: To test whether the number of observed events is consistent
with that of the forecast.

Method: The originial N-test was introduced by Schorlemmer et al (2007)
and modified by Zechar et al (2010). The observed number of events is
given by,

.. math:: N_{obs} = \sum_{m_i, s_j \in R} \omega(m_i, s_j).

Using the simulations described above, the expected number of events is
calculated by summing the simulated number of events over all grid cells

.. math:: \hat{N_x} = \sum_{m_i, s_j \in R} \hat{\omega}_x(m_i, s_j),

where :math:`\hat{\omega}_x(m_i, s_j)` is the simulated number of events
in catalog :math:`x` in spatial cell :math:`s_j` and magnitude cell
:math:`m_i`, generating a set of simulated rates :math:`\{ \hat{N} \}`.
We can then calculate the probability of i) observing at most
:math:`N_{obs}` events and ii) of observing at least :math:`N_{obs}`
events. These probabilities can be written as:

.. math:: \delta_1 =  \frac{ |\{ \hat{N_x} | \hat{N_x} \le N_{obs}\} |}{|\{ \hat{N} \}|}

and

.. math:: \delta_2 =  \frac{ |\{ \hat{N_x} | \hat{N_x} \ge N_{obs}\} |}{|\{ \hat{N} \}|}

If a forecast is Poisson, the expected number of events in the forecast
follows a Poisson distribution with expectation
:math:`N_{fore} = \sum_{m_i, s_j \in R} \lambda(m_i, s_j)`. The
cumulative distribution is then a Poisson cumulative distribution:

.. math:: F(x|N_{fore}) = \exp(-N_{fore}) \sum^{x}_{i=0} \frac{(N_{fore})^i}{i!}

which can be used directly without the need for simulations. The N-test
quantile score is then

.. math:: \delta_1 =  1 - F((N_{obs}-1)|N_{fore}),

and

.. math:: \delta_2 = F(N_{obs}|N_{fore}).

The original N-test considered only :math:`\delta_2` and it’s complement
:math:`1-\delta_2`, which effectively tested the probability of at most
:math:`N_{obs}` events and more than :math:`N_{obs}` events. Very small
or very large values (<0.025 or > 0.975 respectively) were considered to
be inconsistent with the forecast in Schorlemmer et al (2010). However
the approach above aims to test something subtely different, that is at
least :math:`N_{obs}` events and at most :math:`N_{obs}` events. Zechar
et al (2010a) recommends testing both :math:`\delta_1` and
:math:`\delta_2` with an effective significance of have the required
significance level, so for a required significance level of 0.05, a
forecast is consistent if both :math:`\delta_1` and :math:`\delta_2` are
greater than 0.025. A very small :math:`\delta_1` suggest the rate is
too low while a very low :math:`\delta_2` suggests a rate which is too
high to be consistent with observations.

Implementation in pyCSEP

pyCSEP uses the Zechar et al (2010) version of the N-test and the
cumulative Poisson approach to estimate the range of expected events
from the forecasts, so does not implement a simulation in this case. The
upper and lower bounds for the test are determined from the cumulative
Poisson distribution. ``number_test_result.quantile`` will return both
:math:`\delta_1` and :math:`\delta_2` values.

.. code:: ipython3

    number_test_result = poisson.number_test(helmstetter, catalog)
    ax = plots.plot_poisson_consistency_test(
        number_test_result,
        plot_args={'xlabel':'Number of events'}
    )



.. image:: static/output_13_0.png


In this case, the black bar shows the :math:`95\%` interval for the
number of events in the forecast. The actual observed number of events
is shown by the green box, which just passes the N-test in this case:
the forecast generallly expects more events than are observed in
practice, but the observed number falls just within the lower limits of
what is expected so the forecast (just!) passes the N-test.

M-test
^^^^^^

Aim: Establish consistency (or lack thereof) of observed event
magnitudes with forecast magnitudes.

Method: The M-test is first described in Zechar et al. (2010) and aims
to isolate the magnitude component of a forecast. To do this, we sum
over the spatial bins and normalise so that the sum of events matches
the observations.

.. math:: \hat{\boldsymbol{\Omega}}^m = \big{\{}\omega^{m}(m_i)| m_i \in \boldsymbol{M}\big{\}},

where

.. math::  \omega^m(m_i) = \sum_{s_j \in \boldsymbol{S}} \omega(m_i, s_j),

and

.. math:: \boldsymbol{\Lambda}^m = \big{\{} \lambda^m(m_i)| m_i \in \boldsymbol{M} \big{\}},

where

.. math::  \lambda^m(m_i) = \frac{N_{obs}}{N_{fore}}\sum_{s_j \in \boldsymbol{S}} \lambda\big{(}m_i, s_j\big{)}.

Then we compute the joint log-likelihood as we did for the L-test:

.. math::  M = L(\boldsymbol{\Omega}^m | \boldsymbol{\Lambda}^m)

We then wish to compare this with the distribution of simulated
log-likelihoods, this time keep the number of events fixed to

:math:`N_{obs}`. Then for each simulated catalogue,
:math:`\hat{M}_x = L(\hat{\boldsymbol{\Omega}}^m | \boldsymbol{\Lambda}^m)`

Quantile score: The final test statistic is again the fraction of
observed log likelihoods within the range of the simulated log
likelihood values:

.. math:: \kappa =  \frac{ |\{ \hat{M_x} | \hat{M_x} \le M\} |}{|\{ \hat{M} \}|}

and the observed magnitudes are inconsistent with the forecast if
:math:`\kappa` is less than the significance level.

pyCSEP implementation

.. code:: ipython3

    mag_test_result = poisson.magnitude_test(
        helmstetter,
        catalog,
        seed=seed,
        num_simulations=nsim
    )
    ax = plots.plot_poisson_consistency_test(
        mag_test_result,
        one_sided_lower=True,
        plot_args={'xlabel':'Normalized likelihood'}
    )



.. image:: static/output_16_0.png


In this example, the forecast passes the M-test, demonstrating that the
magnitude distribution in the forecast is consistent with observed
events. This is shown by the green square marking the joint
log-likelihood for the observed events.

S-test
^^^^^^

Aim: The spatial or S-test aims to establish consistency (or lack
thereof) of observed event locations with a forecast. It is originally
defined in Zechar et al (2010).

Method: Similar to the M-test, but in this case we sum over all
magnitude bins.

.. math:: \hat{\boldsymbol{\Omega}^s} = \{\omega^s(s_j)| s_j \in \boldsymbol{S}\},

where

.. math::  \omega^s(s_j) = \sum_{m_i \in \boldsymbol{M}} \omega(m_i, s_j),

and

.. math:: \boldsymbol{\Lambda}^s = \{ \lambda^s(s_j)| s_j \in \boldsymbol{S} \},

where

.. math::  \lambda^s(s_j) = \frac{N_{obs}}{N_{fore}}\sum_{m_i \in M} \lambda(m_i, s_j).

Then we compute the joint log-likelihood as we did for the L-test or the
M-test:

.. math::  S = L(\boldsymbol{\Omega}^s | \boldsymbol{\Lambda}^s)

We then wish to compare this with the distribution of simulated
log-likelihoods, this time keeping the number of events fixed to
:math:`N_{obs}`. Then for each simulated catalogue,
:math:`\hat{S}_x = L(\hat{\boldsymbol{\Omega}}^s | \boldsymbol{\Lambda}^s)`

The final test statistic is again the fraction of observed log
likelihoods within the range of the simulated log likelihood values:

.. math:: \zeta =  \frac{ |\{ \hat{S_x} | \hat{S_x} \le S\} |}{|\{ \hat{S} \}|}

and again the distinction between a forecast passing or failing the test
depends on our significance level.

pyCSEP implementation

The S-test is again a one-sided test, so we specify this when plotting
the result.

.. code:: ipython3

    spatial_test_result = poisson.spatial_test(
        helmstetter,
        catalog,
        seed=seed,
        num_simulations=nsim
    )
    ax = plots.plot_poisson_consistency_test(
        spatial_test_result,
        one_sided_lower=True,
        plot_args = {'xlabel':'normalized spatial likelihood'}
    )


.. image:: static/output_19_0.png


The Helmstetter model fails the S-test as the observed spatial
likelihood falls in the tail of the simulated likelihood distribution.
Again this is shown by a coloured symbol which highlights whether the
forecast model passes or fails the test.

Forecast comparison tests
~~~~~~~~~~~~~~~~~~~~~~~~~

The consistency tests above check whether a forecast is consistent with
observations, but do not provide a straightforward way to compare two
different forecasts. A few suggestions for this focus on the information
gain of one forecast relative to another (Harte and Vere-Jones 2005,
Imoto and Hurukawa, 2006, Imoto and Rhoades, 2010, Rhoades et al 2011).
The T-test and W-test implementations for earthquake forecast comparison
are first described in Rhoades et al. (2011).

The information gain per earthquake (IGPE) of model A compared to model
B is defined by :math:`I_{N}(A, B) = R/N` where R is the rate-corrected
log-likelihood ratio of models A and B gven by

.. math::  R = \sum_{k=1}^{N}\big{(}\log\lambda_A(i_k) - \log \lambda_B(i_k)\big{)} - \big{(}\hat{N}_A - \hat{N}_B\big{)}

If we set :math:`X_i=\log\lambda_A(k_i)` and
:math:`Y_i=\log\lambda_B(k_i)` then we can define the information gain
per earthquake (IGPE) as

.. math:: I_N(A, B) = \frac{1}{N}\sum^N_{i=1}\big{(}X_i - Y_i\big{)} - \frac{\hat{N}_A - \hat{N}_B}{N}

If :math:`I(A, B)` differs significantly from 0, the model with the
lower likelihood can be rejected in favour of the other.

t-test

If :math:`X_i - Y_i` are independent and come from the same normal
population with mean :math:`\mu` then we can use the classic paired
t-test to evaluate the null hypothesis that
:math:`\mu = (\hat{N}_A - \hat{N}_B)/N` against the alternative
hypothesis :math:`\mu \ne (\hat{N}_A - \hat{N}_B)/N`. To implement this,
we let :math:`s` denote the sample variance of :math:`(X_i - Y_i)` such
that

.. math::  s^2 = \frac{1}{N-1}\sum^N_{i=1}\big{(}X_i - Y_i\big{)}^2 - \frac{1}{N^2 - N}\bigg{(}\sum^N_{i=1}\big{(}X_i - Y_i\big{)}\bigg{)}^2

Under the null hypothesis
:math:`T = I_N(A, B)\big{/}\big{(}s/\sqrt{N}\big{)}` has a
t-distribution with :math:`N-1` degrees of freedom and the null
hypothesis can be rejected if :math:`|T|` exceeds a critical value of
the :math:`t_{N-1}` distribution. The confidence intervals for
:math:`\mu - (\hat{N}_A - \hat{N}_B)/N` can then be constructed with the
form :math:`I_N(A,B) \pm ts/\sqrt{N}` where t is the appropriate
quantile of the :math:`t_{N-1}` distribution.

W-test

An alternative to the t-test is the Wilcoxan signed-rank test or W-test.
This is a non-parameteric alternative to the t-test which can be used if
we do not feel the assumption of normally distributed differences in
:math:`X_i - Y_i` is valid. This assumption might b particularly poor
when we have small sample sizes. The W-test instead depends on the
(weaker) assumption that :math:`X_i - Y_i` is symmetric and tests
whether the meadian of :math:`X_i - Y_i` is equal to
:math:`(\hat{N}_A - \hat{N}_B)/N`. The W-test is less powerful than the
T-test for normally distributed differences and cannot reject the null
hypothesis (with :math:`95\%` confidence) for very small sample sizes
(:math:`N \leq 5`).

The t-test becomes more accurate as :math:`N \rightarrow \infty` due to
the central limit theorem and therefore the t-test is considered
dependable for large :math:`N`. Where :math:`N` is small, a model might
only be considered more informative if both the t- and W-test results
agree.

Implementation in pyCSEP

The t-test and W-tests are implemented in pyCSEP as below.

.. code:: ipython3

    helmstetter_ms = csep.load_gridded_forecast(
        datasets.helmstetter_mainshock_fname,
        name = "Helmstetter Mainshock"
    )

    t_test = poisson.paired_t_test(helmstetter, helmstetter_ms, catalog)
    w_test = poisson.w_test(helmstetter, helmstetter_ms, catalog)
    comp_args = {'title': 'Paired T-test Result',
                 'ylabel': 'Information gain',
                 'xlabel': '',
                 'xticklabels_rotation': 0,
                 'figsize': (6,4)}

    ax = plots.plot_comparison_test([t_test], [w_test], plot_args=comp_args)



.. image:: static/output_22_0.png


The first argument to the ``paired_t_test`` function is taken as model A
and the second as our basline model, or model B. When plotting the
result, the horizontal dashed line indicates the performance of model B
and the vertical bar shows the confidence bars for the information gain
:math:`I_N(A, B)` associated with model A relative to model B. In this
case, the model with aftershocks performs statistically worse than the
benchmark model. We note that this comparison is used for demonstation
purposes only.

Catalog-based forecast tests
----------------------------

Catalog-based forecast tests evaluate forecasts using simulated outputs
in the form of synthetic earthquake catalogs. Thus, removing the need
for the Poisson approximation and simulation procedure used with
grid-based forecasts. We know that observed seismicity is overdispersed
with respect to a Poissonian model due to spatio-temporal clustering.
Overdispersed models are more likely to be rejected by the original
Poisson-based CSEP tests (Werner et al, 2011a). This modification of the
testing framework allows for a broader range of forecast models. The
distribution of realizations is then compared with observations, similar
to in the grid-based case. These tests were developed by Savran et al
2020, who applied them to test forecasts following the 2019 Ridgecrest
earthquake in Southern California.

In the following text, we show how catalog-based forecasts are defined.
Again we begin by defining a region :math:`\boldsymbol{R}` as a function
of some magnitude range :math:`\boldsymbol{M}`, spatial domain
:math:`\boldsymbol{S}` and time period :math:`\boldsymbol{T}`

.. math::  \boldsymbol{R} = \boldsymbol{M} \times \boldsymbol{S} \times \boldsymbol{T}.

An earthquake :math:`e` can be described by a magnitude :math:`m_i` at
some location :math:`s_j` and time :math:`t_k`. A catalog is simply a
collection of earthquakes, thus the observed catalog can be written as

.. math:: \Omega = \big{\{}e_n \big{|} n= 1...N_{obs}; e_n \in \boldsymbol{R} \big{\}},

and a forecast is then specified as a collection of synthetic catalogs
containing events :math:`\hat{e}_{nj}` in domain :math:`\boldsymbol{R}`,
as

.. math::  \boldsymbol{\Lambda} \equiv \Lambda_j = \{\hat{e}_{nj} | n = 1... N_j, j= 1....J ;\hat{e}_{nj} \in \boldsymbol{R} \}.

That is, a forecast consists of :math:`J` simulated catalogs each
containing :math:`N_j` events, described in time, space and magnitude
such that :math:`\hat{e}_{nj}` describes the :math:`n`\ th synthetic
event in the :math:`j`\ th synthetic catalog :math:`\Lambda_j`

When using simulated forecasts in pyCSEP, we must first explicitly
specify the forecast region by specifying the spatial domain and
magnitude regions as below. In effect, these are filters applied to the
forecast and observations to retain only the events in
:math:`\boldsymbol{R}`. The examples in this section are catalog-based
forecast simulations for the Landers earthquake and aftershock sequence
generated using UCERF3-ETAS (Field et al, 2017).

.. code:: ipython3

    # Define the start and end times of the forecasts
    start_time = time_utils.strptime_to_utc_datetime("1992-06-28 11:57:35.0")
    end_time = time_utils.strptime_to_utc_datetime("1992-07-28 11:57:35.0")

    # Magnitude bins properties
    min_mw = 4.95
    max_mw = 8.95
    dmw = 0.1

    # Create space and magnitude regions.
    magnitudes = regions.magnitude_bins(min_mw, max_mw, dmw)
    region = regions.california_relm_region()
    space_magnitude_region = regions.create_space_magnitude_region(
        region,
        magnitudes
    )

    # Load forecast
    forecast = csep.load_catalog_forecast(
        datasets.ucerf3_ascii_format_landers_fname,
        start_time = start_time,
        end_time = end_time,
        region = space_magnitude_region,
        apply_filters = True
    )

    # Compute expected rates
    forecast.filters = [
        f'origin_time >= {forecast.start_epoch}',
        f'origin_time < {forecast.end_epoch}'
    ]
    _ = forecast.get_expected_rates(verbose=False)

    # Obtain Comcat catalog and filter to region
    comcat_catalog = csep.query_comcat(
        start_time,
        end_time,
        min_magnitude=forecast.min_magnitude
    )

    # Filter observed catalog using the same region as the forecast
    comcat_catalog = comcat_catalog.filter_spatial(forecast.region)


.. parsed-literal::

    Fetched ComCat catalog in 0.31937098503112793 seconds.

    Downloaded catalog from ComCat with following parameters
    Start Date: 1992-06-28 12:00:45+00:00
    End Date: 1992-07-24 18:14:36.250000+00:00
    Min Latitude: 33.901 and Max Latitude: 36.705
    Min Longitude: -118.067 and Max Longitude: -116.285
    Min Magnitude: 4.95
    Found 19 events in the ComCat catalog.


Number Test
~~~~~~~~~~~

Aim: As above, the number test aims to evaluate if the number of
observed events is consistent with the forecast.

Method: The observed statistic in this case is given by
:math:`N_{obs} = |\Omega|`, which is simply the number of events in the
observed catalog. To build the test distribution from the forecast, we
simply count the number of events in each simulated catalog.

.. math::  N_{j} = |\Lambda_c|; j = 1...J

As in the gridded test above, we can then evaluate the probabilities of
at least and at most N events, in this case using the empirical
cumlative distribution function of :math:`F_N`:

.. math:: \delta_1 = P(N_j \geq N_{obs}) = 1 - F_N(N_{obs}-1)

and

.. math:: \delta_2 = P(N_j \leq N_{obs}) = F_N(N_{obs})

Implementation in pyCSEP

.. code:: ipython3

    number_test_result = catalog_evaluations.number_test(
        forecast,
        comcat_catalog,
        verbose=False
    )
    ax = number_test_result.plot()



.. image:: static/output_27_0.png


Plotting the number test result of a simulated catalog forecast displays
a histogram of the numbers of events :math:`\hat{N}_j` in each simulated
catalog :math:`j`, which makes up the test distribution. The test
statistic is shown by the dashed line - in this case it is the number of
observed events in the catalog :math:`N_{obs}`.

Magnitude Test
~~~~~~~~~~~~~~

Aim: The magnitude test aims to test the consistency of the observed
frequency-magnitude distribution with that in the simulated catalogs
that make up the forecast.

Method: The catalog-based magnitude test is implemented quite
differently to the grid-based equivalent. We first define the union
catalog :math:`\Lambda_U` as the union of all simulated catalogs in the
forecast. Formally:

.. math::  \Lambda_U = \{ \lambda_1 \cup \lambda_2 \cup ... \cup \lambda_j \}

| so that the union catalog contains all events across all simulated
  catalogs for a total of
  :math:`N_U = \sum_{j=1}^{J} \big{|}\lambda_j\big{|}` events.
| We then compute the following histograms discretised to the magnitude
  range and magnitude step size (specified earlier for pyCSEP): 1. the
  histogram of the union catalog magnitudes :math:`\Lambda_U^{(m)}` 2.
  Histograms of magnitudes in each of the individual simulated catalogs
  :math:`\lambda_j^{(m)}` 3. the histogram of the observed catalog
  magnitudes :math:`\Omega^{(m)}`

The histograms are normalized so that the total number of events across
all bins is equal to the observed number. The observed statistic is then
calculated as the sum of squared logarithmic residuals between the
normalised observed magnitudes and the union histograms. This statistic
is related to the Kramer von-Mises statistic.

.. math:: d_{obs}= \sum_{k}\Bigg(\log\Bigg[\frac{N_{obs}}{N_U} \Lambda_U^{(m)}(k) + 1\Bigg]- \log\Big[\Omega^{(m)}(k) + 1\Big]\Bigg)^2

where :math:`\Lambda_U^{(m)}(k)` and :math:`\Omega^{(m)}(k)`
represent the count in the :math:`k`\ th bin of the magnitude-frequency
distribution in the union and observed catalogs respectively. We add
unity to each bin to avoid :math:`\log(0)`. We then build the test
distribution from the catalogs in :math:`\boldsymbol{\Lambda}`:

.. math::  D_j =  \sum_{k}\Bigg(\log\Bigg[\frac{N_{obs}}{N_U} \Lambda_U^{(m)}(k) + 1\Bigg]- \log\Bigg[\frac{N_{obs}}{N_j}\Lambda_j^{(m)}(k) + 1\Bigg]\Bigg)^2; j= 1...J

where :math:`\lambda_j^{(m)}(k)` represents the count in the
:math:`k`\ th bin of the magnitude-frequency distribution of the
:math:`j`\ th catalog.

The quantile score can then be calculated using the empirical CDF such
that

.. math::  \gamma_m = F_D(d_{obs})= P(D_j \leq d_{obs})

|  Implementation in pyCSEP
| Hopefully you now see why it was necessary to specify our magnitude
  range explicitly when we set up the catalog-type testing - we need to
  makes sure the magnitudes are properly discretised for the model we
  want to test.

.. code:: ipython3

    magnitude_test_result = catalog_evaluations.magnitude_test(
        forecast,
        comcat_catalog,verbose=False
    )
    ax = magnitude_test_result.plot(plot_args={'xy': (0.6,0.7)})



.. image:: static/output_30_0.png


The histogram shows the resulting test distribution with :math:`D^*`
calculated for each simulated catalog as described in the method above.
The test statistic :math:`\omega = d_{obs}` is shown with the dashed
horizontal line. The quantile score for this forecast is
:math:`\gamma = 0.66`.

Pseudo-likelihood test
~~~~~~~~~~~~~~~~~~~~~~

Aim : The pseudo-likelihood test aims to evaluate the likelihood of a
forecast given an observed catalog.

Method : The pseudo-likelihood test has similar aims to the grid-based
likelihood test above, but its implementation differs in a few
significant ways. Firstly, it does not compute an actual likelihood
(hence the name pseudo-likelihood), and instead of aggregating over
cells as in the grid-based case, the pseudo-likelihood test aggregates
likelihood over target event likelihood scores (so likelihood score per
target event, rather than likelihood score per grid cell). The most
important difference, however, is that the pseudo-likelihood tests do
not use a Poisson likelihood.

The pseudo-likelihood approach is based on the continuous point process
likelihood function. A continuous marked space-time point process can be
specified by a conditional intensity function
:math:`\lambda(\boldsymbol{e}|H_t)`, in which :math:`H_t` describes the
history of the process in time. The log-likelihood function for any
point process in :math:`\boldsymbol{R}` is given by

.. math::  L = \sum_{i=1}^{N} \log \lambda(e_i|H_t) - \int_{\boldsymbol{R}}\lambda(\boldsymbol{e}|H_t)d\boldsymbol{R}

Not all models will have an explicit likelihood function, so instead we
approximate the expectation of :math:`\lambda(e|H_t)` using the forecast
catalogs. The approximate rate density is defined as the conditional
expectation given a discretised region :math:`R_d` of the continuous
rate

.. math:: \hat{\lambda}(\boldsymbol{e}|H_t) = E\big[\lambda(\boldsymbol{e}|H_t)|R_d\big]

We still regard the model as continuous, but the rate density is
approximated within a single cell. This is analogous to the gridded
approach where we count the number of events in discrete cells. The
pseudo-loglikelihood is then

.. math:: \hat{L} = \sum_{i=1}^N \log \hat{\lambda}(e_i|H_t) - \int_R \hat{\lambda}(\boldsymbol{e}|H_t) dR

and we can write the approximate rate density as

.. math:: \hat{\lambda}(\boldsymbol{e}|H_t) = \sum_M \hat{\lambda}(\boldsymbol{e}|H_t),

where we take the sum over all magnitude bins :math:`M`. We can
calculate observed pseudolikelihood as

.. math::  \hat{L}_{obs} = \sum_{i=1}^{N_{obs}} \log \hat{\lambda}_s(k_i) - \bar{N},

where :math:`\hat{\lambda}_s(k_i)` is the approximate rate density in
the :math:`k`\ th spatial cell and :math:`k_i` denotes the spatil cell
in which the :math:`i`\ th event occurs. :math:`\bar{N}` is the expected
number of events in :math:`R_d`. Similarly, we calculate the test
distribution as

.. math:: \hat{L}_{j} = \Bigg[\sum_{i=1}^{N_{j}} \log\hat{\lambda}_s(k_{ij}) - \bar{N}\Bigg]; j = 1....J,

where :math:`\hat{\lambda}_s(k_{ij})` describes the approximate rate
density of the :math:`i`\ th event in the :math:`j`\ th catalog. We can
then calculate the quantile score as

.. math::  \gamma_L = F_L(\hat{L}_{obs})= P(\hat{L}_j \leq \hat{L}_{obs}).

Implementation in pyCSEP

.. code:: ipython3

    pseudolikelihood_test_result = catalog_evaluations.pseudolikelihood_test(
        forecast,
        comcat_catalog,
        verbose=False
    )
    ax = pseudolikelihood_test_result.plot()



.. image:: static/output_33_0.png


The histogram shows the test distribution of pseudolikelihood as
calculated above for each catalog :math:`j`. The dashed vertical line
shows the observed statistic :math:`\hat{L}_{obs} = \omega`. It is clear
that the observed statistic falls within the critical region of test
distribution, as reflected in the quantile score of
:math:`\gamma_L = 0.02`.

Spatial test
~~~~~~~~~~~~

Aim: The spatial test again aims to isolate the spatial component of the
forecast and test the consistency of spatial rates with observed events.

Method We perform the spatial test in the catalog-based approach in a
similar way to the grid-based spatial test approach: by normalising the
approximate rate density. In this case, we use the normalisation
:math:`\hat{\lambda}_s = \hat{\lambda}_s \big/ \sum_{R} \hat{\lambda}_s`.
Then the observed spatial test statistic is calculated as

.. math::  S_{obs} = \Bigg[\sum_{i=1}^{N_{obs}} \log \hat{\lambda}_s^*(k_i)\Bigg]N_{obs}^{-1}

in which :math:`\hat{\lambda}_s^*(k_i)` is the normalised approximate
rate density in the :math:`k`\ th cell corresponding to the
:math:`i`\ th event in the observed catalog :math:`\Omega`. Similarly,
we define the test distribution using

.. math::  S_{c} = \bigg[\sum_{i=1}^{N_{j}} \log \hat{\lambda}_s^*(k_{ij})\bigg]N_{j}^{-1}; j= 1...J

for each catalog j. Finally, the quantile score for the spatial test is
determined by once again comparing the observed and test distribution
statistics:

.. math:: \gamma_s = F_s(\hat{S}_{obs}) = P (\hat{S}_j \leq \hat{S}_{obs})

Implementation in pyCSEP

.. code:: ipython3

    spatial_test_result = catalog_evaluations.spatial_test(
        forecast,
        comcat_catalog,
        verbose=False
    )
    ax = spatial_test_result.plot()



.. image:: static/output_36_0.png


The histogram shows the test distribution of normalised
pseduo-likelihood computed for each simulated catalog :math:`j`. The
dashed vertical line shows the observed test statistic
:math:`s_{obs} = \omega = -5.88`, which is clearly within the test
distribution. The quantile score :math:`\gamma_s = 0.36` is also printed
on the figure by default.

References
----------

Field, E. H., K. R. Milner, J. L. Hardebeck, M. T. Page, N. J. van der
Elst, T. H. Jordan, A. J. Michael, B. E. Shaw, and M. J. Werner (2017).
A spatiotemporal clustering model for the third Uniform California
Earthquake Rupture Forecast (UCERF3-ETAS): Toward an operational
earthquake forecast, Bull. Seismol. Soc. Am. 107, 1049–1081.

Harte, D., and D. Vere-Jones (2005), The entropy score and its uses in
earthquake forecasting, Pure Appl. Geophys. 162 , 6-7, 1229-1253, DOI:
10.1007/ s00024-004-2667-2.

Helmstetter, A., Y. Y. Kagan, and D. D. Jackson (2006). Comparison of
short-term and time-independent earthquake forecast models for southern
California, Bulletin of the Seismological Society of America 96 90-106.

Imoto, M., and N. Hurukawa (2006), Assessing potential seismic activity
in Vrancea, Romania, using a stress-release model, Earth Planets Space
58 , 1511-1514.

Imoto, M., and D.A. Rhoades (2010), Seismicity models of moderate
earthquakes in Kanto, Japan utilizing multiple predictive parameters,
Pure Appl. Geophys. 167, 6-7, 831-843, DOI: 10.1007/s00024-010-0066-4.

Rhoades, D.A, D., Schorlemmer, M.C.Gerstenberger, A. Christophersen, J.
D. Zechar & M. Imoto (2011) Efficient testing of earthquake forecasting
models, Acta Geophysica 59

Savran, W., M. J. Werner, W. Marzocchi, D. Rhoades, D. D. Jackson, K. R.
Milner, E. H. Field, and A. J. Michael (2020). Pseudoprospective
evaluation of UCERF3-ETAS forecasts during the 2019 Ridgecrest Sequence,
Bulletin of the Seismological Society of America.

Schorlemmer, D., and M.C. Gerstenberger (2007), RELM testing center,
Seismol. Res. Lett. 78, 30–36.

Schorlemmer, D., M.C. Gerstenberger, S. Wiemer, D.D. Jackson, and D.A.
Rhoades (2007), Earthquake likelihood model testing, Seismol. Res. Lett.
78, 17–29.

Schorlemmer, D., A. Christophersen, A. Rovida, F. Mele, M. Stucci and W.
Marzocchi (2010a). Setting up an earthquake forecast experiment in
Italy, Annals of Geophysics, 53, no.3

Schorlemmer, D., J.D. Zechar, M.J. Werner, E.H. Field, D.D. Jackson, and
T.H. Jordan (2010b), First results of the Regional Earthquake Likelihood
Models experiment, Pure Appl. Geophys., 167, 8/9,
doi:10.1007/s00024-010-0081-5.

M. Taroni, W. Marzocchi, D. Schorlemmer, M. J. Werner, S. Wiemer, J. D.
Zechar, L. Heiniger, F. Euchner; Prospective CSEP Evaluation of 1‐Day,
3‐Month, and 5‐Yr Earthquake Forecasts for Italy. Seismological Research
Letters 2018;; 89 (4): 1251–1261. doi:
https://doi.org/10.1785/0220180031

Werner, M. J., A. Helmstetter, D. D. Jackson, and Y. Y. Kagan (2011a).
High-Resolution Long-Term and Short-Term Earthquake Forecasts for
California, Bulletin of the Seismological Society of America 101
1630-1648

Werner, M.J. J.D. Zechar, W. Marzocchi, and S. Wiemer (2011b),
Retrospective evaluation of the five-year and ten-year CSEP-Italy
earthquake forecasts, Annals of Geophysics 53, no. 3, 11–30,
doi:10.4401/ag-4840.

Zechar, 2011: Evaluating earthquake predictions and earthquake
forecasts: a guide for students and new researchers, CORSSA
(http://www.corssa.org/en/articles/theme_6/)

Zechar, J.D., M.C. Gerstenberger, and D.A. Rhoades (2010a),
Likelihood-based tests for evaluating space-rate-magnitude forecasts,
Bull. Seis. Soc. Am., 100(3), 1184—1195, doi:10.1785/0120090192.

Zechar, J.D., D. Schorlemmer, M. Liukis, J. Yu, F. Euchner, P.J.
Maechling, and T.H. Jordan (2010b), The Collaboratory for the Study of
Earthquake Predictability perspective on computational earthquake
science, Concurr. Comp-Pract. E., doi:10.1002/cpe.1519.
===========================
Core Concepts for Beginners
===========================

If you are reading this documentation, there is a good chance that you are developing/evaluating an earthquake forecast or
implementing an experiment at a CSEP testing center. This section will help you understand how we conceptualize forecasts,
evaluations, and earthquake catalogs. These components make up the majority of the PyCSEP package. We also include some
prewritten visualizations along with some utilities that might be useful in your work.

Catalogs
========
Earthquake catalogs are fundamental to both forecasts and evaluations and make up a core component of the PyCSEP package.
At some point you will be working with catalogs if you are evaluating earthquake forecasts.

One major difference between PyCSEP and a project like `ObsPy <https://docs.obspy.org/>`_ is that typical 'CSEP' calculations
operate on an entire catalog at once to perform methods like filtering and binning that are required to evaluate an earthquake
forecast. We provide earthquake catalog classes that follow the interface defined by
:class:`AbstractBaseCatalog <csep.core.catalogs.AbstractBaseCatalog>`.

The catalog data are stored internally as a `structured Numpy array <https://numpy.org/doc/stable/user/basics.rec.html>`_
which effectively treats events contiguously in memory like a c-style struct. This allows us to accelerate calculations
using the vectorized operations provided by Numpy. The necessary attributes for an event to be used
in an evaluation are the spatial location (lat, lon), magnitude, and origin time. Additionally, depth and other identifying
characteristics can be used. The `default storage format <https://scec.usc.edu/scecpedia/CSEP_2_CATALOG_FORMAT>`_ for
an earthquake catalog is an ASCII/utf-8 text file with events stored in CSV format.

The :class:`AbstractBaseCatalog <csep.core.catalogs.AbstractBaseCatalog>` can be extended to accommodate different catalog formats
or input and output routines. For example :class:`UCERF3Catalog <csep.core.catalogs.UCERF3Catalog>` extends this class to deal
with the big-endian storage routine from the `UCERF3-ETAS <https://github.com/opensha/opensha-ucerf3>`_ forecasting model. More
information will be included in the :ref:`catalogs-reference` section of the documentation.

Forecasts
=========

PyCSEP provides objects for interacting with :ref:`earthquake forecasts <earthquake_forecast>`. PyCSEP supports two types
of earthquake forecasts, and provides separate objects for interacting with both. The forecasts share similar
characteristics, but, conceptually, they should be treated differently because they require different types of evaluations.

Both time-independent and time-dependent forecasts are represented using the same PyCSEP forecast objects. Typically, for
time-dependent forecasts, one would create separate forecast objects for each time period. As the name suggests,
time-independent forecasts do not change with time.


Grid-based forecast
-------------------

Grid-based earthquake forecasts are specified by the expected rate of earthquakes within discrete, independent
space-time-magnitude bins. Within each bin, the expected rate represents the parameter of a Poisson distribution. For, details
about the forecast objects visit the :ref:`forecast-reference` section of the documentation.

The forecast object contains three main components: (1) the expected earthquake rates, (2) the
:class:`spatial region <csep.utils.spatial.CartesianGrid2D>` associated with the rates, and (3) the magnitude range
associated with the expected rates. The spatial bins are usually discretized according to the geographical coordinates
latitude and longitude with most previous CSEP spatial regions defining a spatial size of 0.1° x 0.1°. Magnitude bins are
also discretized similarly, with 0.1 magnitude units being a standard choice. PyCSEP does not enforce constraints on the
bin-sizes for both space and magnitude, but the discretion must be regular.


Catalog-based forecast
----------------------

Catalog-based forecasts are specified by families of synthetic earthquake catalogs that are generated through simulation
by probabilistic models. Each catalog represents a stochastic representation of seismicity consistent with the forecasting
model. Probabilistic statements are made by computing statistics (usually by counting) within the family of synthetic catalogs,
which can be as simple as counted the number of events in each catalog. These statistics represent the full-distribution of outcomes as
specified by the forecasting models, thereby allowing for more direct assessments of the models that produce them.

Within PyCSEP catalog forecasts are effectively lists of earthquake catalogs, no different than those obtained from
authoritative sources. Thus, any operation that can be performed on an observed earthquake catalog can be performed on a
synthetic catalog from a catalog-based forecast.

It can be useful to count the numbers of forecasted earthquakes within discrete space-time bins (like those used for
grid-based forecasts). Therefore, it's common to have a :class:`spatial region <csep.utils.spatial.CartesianGrid2D>` and
set of magnitude bins associated with a forecast. Again, the only rules that PyCSEP enforces are that the space-magnitude
regions are regularly discretized.

Evaluations
===========

PyCSEP provides implementations of statistical tests used to evaluate both grid-based and catalog-based earthquake forecasts.
The former use parametric evaluations based on Poisson likelihood functions, while the latter use so-called 'likelihood-free'
evaluations that are computed from empirical distributions provided by the forecasts. Details on the specific implementation
of the evaluations will be provided in the :ref:`evaluation-reference` section.

Every evaluation can be different, but in general, the evaluations need the following information:

1. Earthquake forecast(s)

    * Spatial region
    * Magnitude range

2. Authoritative earthquake catalog

PyCSEP does not produce earthquake forecasts, but provides the ability to represent them using internal data models to
facilitate their evaluation. General advice on how to administer the statistical tests will be provided in the
:ref:`evaluation-reference` section.