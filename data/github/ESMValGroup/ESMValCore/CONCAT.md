# ESMValCore package

[![Documentation Status](https://readthedocs.org/projects/esmvaltool/badge/?version=latest)](https://esmvaltool.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3387139.svg)](https://doi.org/10.5281/zenodo.3387139)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ESMValGroup?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![CircleCI](https://circleci.com/gh/ESMValGroup/ESMValCore/tree/main.svg?style=svg)](https://circleci.com/gh/ESMValGroup/ESMValCore/tree/main)
[![codecov](https://codecov.io/gh/ESMValGroup/ESMValCore/branch/main/graph/badge.svg?token=wQnDzguwq6)](https://codecov.io/gh/ESMValGroup/ESMValCore)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/5d496dea9ef64ec68e448a6df5a65783)](https://www.codacy.com/gh/ESMValGroup/ESMValCore?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ESMValGroup/ESMValCore&amp;utm_campaign=Badge_Grade)
[![Docker Build Status](https://img.shields.io/docker/cloud/build/esmvalgroup/esmvalcore)](https://hub.docker.com/r/esmvalgroup/esmvalcore/)
[![Anaconda-Server Badge](https://anaconda.org/esmvalgroup/esmvalcore/badges/installer/conda.svg)](https://conda.anaconda.org/esmvalgroup)

![esmvaltoollogo](https://github.com/ESMValGroup/ESMValCore/blob/main/doc/figures/ESMValTool-logo-2.png)

ESMValCore: core functionalities for the ESMValTool, a community diagnostic
and performance metrics tool for routine evaluation of Earth System Models
in the Climate Model Intercomparison Project (CMIP).

# Getting started

Please have a look at the
[documentation](https://docs.esmvaltool.org/projects/esmvalcore/en/latest/quickstart/install.html)
to get started.

## Using the ESMValCore package to run recipes

The ESMValCore package provides the `esmvaltool` command, which can be used to run
[recipes](https://docs.esmvaltool.org/projects/esmvalcore/en/latest/recipe/overview.html)
for working with CMIP-like data.
A large collection of ready to use
[recipes and diagnostics](https://docs.esmvaltool.org/en/latest/recipes/index.html)
is provided by the
[ESMValTool](https://github.com/ESMValGroup/ESMValTool)
package.

## Using ESMValCore as a Python library

The ESMValCore package provides various functions for:

-   Finding data in a directory structure typically used for CMIP data.

-   Reading CMIP/CMOR tables and using those to check model and observational data.

-   ESMValTool preprocessor functions based on
    [iris](https://scitools-iris.readthedocs.io) for e.g. regridding,
    vertical interpolation, statistics, correcting (meta)data errors, extracting
    a time range, etcetera.

read all about it in the
[API documentation](https://docs.esmvaltool.org/projects/esmvalcore/en/latest/api/esmvalcore.html).

## Getting help

The easiest way to get help if you cannot find the answer in the documentation
on [readthedocs](https://docs.esmvaltool.org), is to open an
[issue on GitHub](https://github.com/ESMValGroup/ESMValCore/issues).

## Contributing

Contributions are very welcome, please read our
[contribution guidelines](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html)
to get started.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
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

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community. Examples of
representing a project or community include using an official project e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event. Representation of a project may be
further defined and clarified by project maintainers.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting the project team at birgit.hassler@dlr.de or
alistair.sellar@metoffice.gov.uk. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see
https://www.contributor-covenant.org/faq
# Contributions are very welcome

Please read our [contribution guidelines](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html).
<!--
    Thank you for contributing to our project!

    Please do not delete this text completely, but read the text below and keep
    items that seem relevant. If in doubt, just keep everything and add your
    own text at the top, a reviewer will update the checklist for you.

-->

## Description

<!--
    Please describe your changes here, especially focusing on why this pull
    request makes ESMValCore better and what problem it solves.

    Before you start, please read our contribution guidelines: https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html

    Please fill in the GitHub issue that is closed by this pull request,
    e.g. Closes #1903
-->

Closes #issue_number

Link to documentation:

***

## [Before you get started](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#getting-started)

- [ ] [☝ Create an issue](https://github.com/ESMValGroup/ESMValCore/issues) to discuss what you are going to do

## [Checklist](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#checklist-for-pull-requests)

It is the responsibility of the author to make sure the pull request is ready to review. The icons indicate whether the item will be subject to the [🛠 Technical][1] or [🧪 Scientific][2] review.

<!-- The next two lines turn the 🛠 and 🧪 below into hyperlinks -->
[1]: https://docs.esmvaltool.org/en/latest/community/review.html#technical-review
[2]: https://docs.esmvaltool.org/en/latest/community/review.html#scientific-review

- [ ] [🧪][2] The new functionality is [relevant and scientifically sound](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#scientific-relevance)
- [ ] [🛠][1] This pull request has a [descriptive title and labels](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#pull-request-title-and-label)
- [ ] [🛠][1] Code is written according to the [code quality guidelines](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#code-quality)
- [ ] [🧪][2] and [🛠][1] [Documentation](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#documentation) is available
- [ ] [🛠][1] [Unit tests](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#tests) have been added
- [ ] [🛠][1] Changes are [backward compatible](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#backward-compatibility)
- [ ] [🛠][1] Any changed [dependencies have been added or removed](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#dependencies) correctly
- [ ] [🛠][1] The [list of authors](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#list-of-authors) is up to date
- [ ] [🛠][1] All [checks below this pull request](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/contributing.html#pull-request-checks) were successful

***

To help with the number pull requests:

- 🙏 We kindly ask you to [review](https://docs.esmvaltool.org/en/latest/community/review.html#review-of-pull-requests) two other [open pull requests](https://github.com/ESMValGroup/ESMValCore/pulls) in this repository
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is. If you are developing a new diagnostic script, please provide a link to the code/branch on GitHub that you are working in.

**Please attach**
  - The recipe that you are trying to run, you can find a copy in the `run` directory in the output directory
  - The `main_log_debug.txt` file, this can also be found in the `run` directory in the output directory
---
name: Dataset issue report
about: Create a report about a problematic dataset
title: 'Dataset problem: '
labels: ''
assignees: valeriupredoi, zklaus

---

**Describe the dataset issue**
A clear and concise description of what the problem with the data is is. Here are some guidelines that you can use
right out the box by filling in the blanks (note that, if needed we will assist you with writing a fix; also note that
if the problem is with a CMIP6 dataset, we will alert the ESGF and/or the model developers for them to fix the problem
in the long run; fixes for CMIP5 by the model developers are no longer possible):

- Data file has been changed by you in any way (if you answer yes the issue will be void and closed, we are not
supporting custom-made datasets that present problems, it is your resposability to fix those problems):
- Project (CMIP5/CMIP6/OBS/obs4MIPs/etc):
- Full dataset description (fill out as many as you know, please):
  - dataset:
  - experiment:
  - mip:
  - grid:
  - type:
  - version:
  - tier:
  - variable used:
  - data version:
- Problems encountered (please describe what the actual problems are: e.g. wrong standard name, issue with dimensions):
- Pointer to existing copy of data on ESGF node (it would be very useful if you could provide a physical
fullpath to the file(s) that are causing the problem, e.g. on CEDA Jasmin or DKRZ):
- Other notes and mentions:

**Assign the right people**
If you are already a member of the ESMValTool GitHub project, please assign Valeriu Predoi (valeriupredoi) and
Klaus Zimmermann (zklaus) to the issue. They will then check the issue raised and propagate it further to the
data model developers.

**Please attach**
  - The recipe that you are trying to run, you can find a copy in the `run` directory in the output directory
  - The `main_log_debug.txt` file, this can also be found in the `run` directory in the output directory
  - Attach a text file with the output of `ncdump -h /path/to/file.nc`

Cheers :beer:
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: enhancement
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is.

**Would you be able to help out?**
Would you have the time and skills to implement the solution yourself?
# cmip3-cmor-tables
Text Tables for CMOR1 to create a CMIP3 dataset
# cmip6-cmor-tables

Note:
----
Only `CMOR` and `PrePARE` should use the `JSON` file of concatenated `CMIP6` Control Vocabulary [CMIP6_CV.json](https://github.com/PCMDI/cmip6-cmor-tables/blob/master/Tables/CMIP6_CV.json) found in this repository.  

All other software should get the CV's directly from the [WCRP](https://github.com/WCRP-CMIP/CMIP6_CVs) repository. 

CMIP6 Data Request can be found in the following [mip tables](http://clipc-services.ceda.ac.uk/dreq/index/miptable.html).
