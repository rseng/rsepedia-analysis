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
