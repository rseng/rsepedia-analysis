[![Maintenance](https://img.shields.io/badge/Maintained%3F-yes-green.svg)](https://GitHub.com/Naereen/StrapDown.js/graphs/commit-activity)
[![made-with-python](https://img.shields.io/badge/Made%20with-Python-1f425f.svg)](https://www.python.org/)
[![Documentation Status](https://readthedocs.org/projects/esmvaltool/badge/?version=latest)](https://esmvaltool.readthedocs.io/en/latest/?badge=latest)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3401363.svg)](https://doi.org/10.5281/zenodo.3401363)
[![Gitter](https://badges.gitter.im/Join%20Chat.svg)](https://gitter.im/ESMValGroup?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![CircleCI](https://circleci.com/gh/ESMValGroup/ESMValTool/tree/main.svg?style=svg)](https://circleci.com/gh/ESMValGroup/ESMValTool/tree/main)
[![example branch parameter](https://github.com/github/docs/actions/workflows/main.yml/badge.svg?branch=main)](https://github.com/ESMValGroup/ESMValTool/actions)
[![Codacy Badge](https://app.codacy.com/project/badge/Coverage/79bf6932c2e844eea15d0fb1ed7e415c)](https://www.codacy.com/gh/ESMValGroup/ESMValTool?utm_source=github.com&utm_medium=referral&utm_content=ESMValGroup/ESMValTool&utm_campaign=Badge_Coverage)
[![Codacy Badge](https://app.codacy.com/project/badge/Grade/79bf6932c2e844eea15d0fb1ed7e415c)](https://www.codacy.com/gh/ESMValGroup/ESMValTool?utm_source=github.com&amp;utm_medium=referral&amp;utm_content=ESMValGroup/ESMValTool&amp;utm_campaign=Badge_Grade)
[![Docker Build Status](https://img.shields.io/docker/cloud/build/esmvalgroup/esmvaltool.svg)](https://hub.docker.com/r/esmvalgroup/esmvaltool/)
[![Anaconda-Server Badge](https://anaconda.org/esmvalgroup/esmvaltool/badges/installer/conda.svg)](https://conda.anaconda.org/esmvalgroup)

![esmvaltoollogo](https://github.com/ESMValGroup/ESMValTool/blob/main/doc/sphinx/source/figures/ESMValTool-logo-2.png)

- [**Documentation**](https://docs.esmvaltool.org/en/latest/)
- [**ESMValTool Website**](https://www.esmvaltool.org/)
- [**ESMValTool Tutorial**](https://esmvalgroup.github.io/ESMValTool_Tutorial/index.html)
- [**ESMValGroup Project on GitHub**](https://github.com/ESMValGroup)
- [**Gallery**](https://docs.esmvaltool.org/en/latest/gallery.html)

# Introduction

ESMValTool is a community-developed climate model diagnostics and evaluation software package, driven
both by computational performance and scientific accuracy and reproducibility. ESMValTool is open to both
users and developers, encouraging open exchange of diagnostic source code and evaluation results from the
Coupled Model Intercomparison Project [CMIP](https://www.wcrp-climate.org/wgcm-cmip) ensemble. For a
comprehensive introduction to ESMValTool please visit our
[documentation](https://docs.esmvaltool.org/en/latest/introduction.html) page.

# Running esmvaltool

Diagnostics from ESMValTool are run using [recipe](https://docs.esmvaltool.org/en/latest/recipes/index.html)
files that contain pointers to the requested data types, directives for the preprocessing steps that data
will be subject to, and directives for the actual diagnostics that will be run with the now preprocessed data.
Data preprocessing is done via the [ESMValCore](https://docs.esmvaltool.org/projects/ESMValCore/en/latest/quickstart/index.html) package, a pure Python, highly-optimized scientific library, developed by the ESMValTool core developers,
and that performs a number of common analysis tasks
such as regridding, masking, levels extraction etc. [Diagnostics](https://docs.esmvaltool.org/en/latest/develop/diagnostic.html) are written in a variety of programming languages (Python, NCL, R, Julia) and are developed by the wider
scientific community, and included after a scientific and technical review process.

# Input data

ESMValTool can run with the following types of data as input:

- CMIP5
- CMIP6
- OBS, OBS6
- obs4MIPs
- ana4mips
- CORDEX

# Getting started

Please see [getting started](https://docs.esmvaltool.org/en/latest/quickstart/index.html) on readthedocs as well as [ESMValTool tutorial](https://esmvalgroup.github.io/ESMValTool_Tutorial/index.html). The tutorial is a set of lessons that together teach skills needed to work with ESMValTool in climate-related domains.

## Getting help

The easiest way to get help if you cannot find the answer in the documentation on [readthedocs](https://docs.esmvaltool.org), is to open an [issue on GitHub](https://github.com/ESMValGroup/ESMValTool/issues).

## Contributing

If you would like to contribute a new diagnostic or feature, please have a look at our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/index.html).
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

Please read our [contribution guidelines](https://docs.esmvaltool.org/en/latest/community/index.html).
<!--
    Thank you for contributing to our project!

    Please do not delete this text completely, but read the text below and keep
    items that seem relevant. If in doubt, just keep everything and add your
    own text at the top, a reviewer will update the checklist for you.

    While the checklist is intended to be filled in by the technical and scientific
    reviewers, it is the responsibility of the author of the pull request to make
    sure all items on it are properly implemented.
-->

## Description

<!--
    Please describe your changes here, especially focusing on why this pull request makes
    ESMValTool better and what problem it solves.

    Before you start, please read our contribution guidelines: https://docs.esmvaltool.org/en/latest/community/

    Please fill in the GitHub issue that is closed by this pull request, e.g. Closes #1903
-->
- Closes #issue_number
- Link to documentation:

* * *

## Before you get started

<!--
    Please discuss your idea with the development team before getting started,
    to avoid disappointment or unnecessary work later. The way to do this is
    to open a new issue on GitHub.
-->

- [ ] [☝ Create an issue](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#contributing-code-and-documentation) to discuss what you are going to do

## Checklist

It is the responsibility of the author to make sure the pull request is ready to review. The icons indicate whether the item will be subject to the [🛠 Technical][1] or [🧪 Scientific][2] review.

<!-- The next two lines turn the 🛠 and 🧪 below into hyperlinks -->
[1]: https://docs.esmvaltool.org/en/latest/community/review.html#technical-review
[2]: https://docs.esmvaltool.org/en/latest/community/review.html#scientific-review

- [ ] [🛠][1] This pull request has a [descriptive title](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#pull-request-title)
- [ ] [🛠][1] Code is written according to the [code quality guidelines](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#code-quality)
- [ ] [🛠][1] [Documentation](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#documentation) is available
- [ ] [🛠][1] [Tests](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#tests) run successfully
- [ ] [🛠][1] The [list of authors](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#list-of-authors) is up to date
- [ ] [🛠][1] Any changed dependencies have been [added or removed correctly](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#dependencies)
- [ ] [🛠][1] All [checks below this pull request](https://docs.esmvaltool.org/en/latest/community/code_documentation.html#pull-request-checks) were successful

### [New or updated recipe/diagnostic](https://docs.esmvaltool.org/en/latest/community/diagnostic.html)

- [ ] [🧪][2] [Recipe runs successfully](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#testing-recipes)
- [ ] [🧪][2] [Recipe is well documented](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recipe-and-diagnostic-documentation)
- [ ] [🧪][2] [Figure(s) and data](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#diagnostic-output) look as expected from literature
- [ ] [🛠][1] [Provenance information](https://docs.esmvaltool.org/en/latest/community/diagnostic.html#recording-provenance) has been added

### [New or updated data reformatting script](https://docs.esmvaltool.org/en/latest/develop/dataset.html)

- [ ] [🛠][1] [Documentation](https://docs.esmvaltool.org/en/latest/community/dataset.html#dataset-documentation) is available
- [ ] [🛠][1] The dataset has been [added to the CMOR check recipe](https://docs.esmvaltool.org/en/latest/community/dataset.html#testing)
- [ ] [🧪][2] Numbers and units of the data look [physically meaningful](https://docs.esmvaltool.org/en/latest/community/dataset.html#scientific-sanity-check)

***

To help with the number of pull requests:

-  🙏 We kindly ask you to [review](https://docs.esmvaltool.org/en/latest/community/review.html#review-of-pull-requests) two other [open pull requests](https://github.com/ESMValGroup/ESMValTool/pulls) in this repository

<!--
If you need help with any of the items on the checklists above, please do not hesitate to ask by commenting in the issue or pull request.
-->
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
name: New diagnostic
about: Develop a new diagnostic.
title: ''
labels: diagnostic
assignees: ''

---

**Short description of the diagnostic**
Add a short description of the diagnostic that you would like to add.

**Branch and pull request**
Once you've started working, add the branch (and pull request)
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
ESMVal python library
=====================

General structure
-----------------

This directory contains usefull scripts and tools for ESMVal which are to support the writing of diagnostics in a clean and short way.

In general you have two options to include code into this library.

A) you have some small code snippet that is needed somewhere else (e.g. a special function).
In this case you should put your code into either a new file in the 'python' directory or an already existing file.
Let's say that you have developed some statistics test functionality, then it would be good to put it in a file statistics.py

However, more importantly it is important to recognize that python is in general organized in modules (http://docs.python.org/2/tutorial/modules.html#more-on-modules).
Modules shall ever been used if you are working with bigger projects or code which encapsulated some functionality in different code parts (e.g. a bigger statistical library, a plotting library ...)

B) Modules are simply organized in directories which and contain an additional __init__.py file which tells the python interpreter which files are part of the module.
Using modules allows you to load entire packages as a whole or parts of it into your code using e.g.

from scipy import stats

Thus the ESMVal python library is supposed to look like e.g.

/lib+
    |_python+
            |_my_cool_collection_of_snippets.py
            |_plots.py
            |_statlib+
                     |___init__.py
                     |_linear_regressions.py
                     |_pdf.py
                     |_stattests.py


Further reading
---------------
For further reading reagarding modules in python look at the following links:

* http://docs.python.org/2/tutorial/modules.html
