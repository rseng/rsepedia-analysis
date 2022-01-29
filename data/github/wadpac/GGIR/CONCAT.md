![](vignettes/GGIR-MASTERLOGO-RGB.png)

![GitHub Actions R-CMD-check](https://github.com/wadpac/GGIR/workflows/R-CMD-check-full/badge.svg)
[![codecov](https://codecov.io/gh/wadpac/GGIR/branch/master/graph/badge.svg)](https://app.codecov.io/gh/wadpac/GGIR)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1051064.svg)](https://doi.org/10.5281/zenodo.1051064)
[![](https://cranlogs.r-pkg.org/badges/last-month/GGIR)](https://cran.r-project.org/package=GGIR)

## Getting started:
The package [vignette](https://CRAN.R-project.org/package=GGIR/vignettes/GGIR.html) and [this](https://youtu.be/S8YPTrYNWdU) short tutorial video provide an introduction to GGIR, including: How it can be installed, Key software features, and where to get help.

## Contribution guidelines:
We always welcome contributions to the package.
If you want to contribute to the development of GGIR, have a look at the [contribution guidelines](https://github.com/wadpac/GGIR/blob/master/CONTRIBUTING.md).

### Images usaged
The copyright of the GGIR logo as contained in the file vignettes/GGIR-MASTERLOGO-RGB.png lies with Accelting (Almere, The Netherlands), please contact v.vanhees@accelting.com to ask for permission to use this logo.

All other images in this repository are released under the Creative Commons Attribution 4.0 International (CC BY 4.0) license.

### Research notice

Please note that this repository is participating in a study into sustainability of open source projects. Data will be gathered about this repository for approximately the next 12 months, starting from June 2021.

Data collected will include number of contributors, number of PRs, time taken to close/merge these PRs, and issues closed.

For more information, please visit [the informational page](https://sustainable-open-science-and-software.github.io/) or download the [participant information sheet](https://sustainable-open-science-and-software.github.io/assets/PIS_sustainable_software.pdf).

# Version numbering

We use version encoding **A.B-C**:

- A increases with major changes that affect backward compatibility with previous releases like changes in function names, function arguments or file format.
- B increases with every CRAN release. We aim to avoid more than four CRAN releases per year.
- C increases with every GitHub release. We aim to avoid more than one GitHub release per month.

# GitHub releases

Before releasing, please make sure to check the following:

1. Create GitHub issue at least 1 weeks before the intended release to announce the release and indicate what will be in the release.
2. Make sure the change log `inst/NEWS.Rd` is up to date and that it says "GitHub-only-release date" rather than "release date"
3. Make sure the third (last) digit in the version number is incremented by one relative to the master branch and the date is the present date. This applies to the files `DESCRIPTION`, `CITATION.cff` (not the cff-version, but the version on line 56 of the .cff-file), `GGIR-package.Rd` and `NEWS.Rd` file. Use function `prepareNewRelease.R` in the root of GGIR to double check that version number and date are consistent between these files.
4. Update package contributor list if new people have contributed.
5. Run `R CMD check --as-cran` to make sure all tests and checks pass.

Note that GitHub releases require a release name. We typically choose a random name of a city or town in South America. Whatever you choose this should be an easy to read and remember word.

# CRAN releases

To do a CRAN release, follow the following steps:

1. Create GitHub issue at least 4 weeks before the intended CRAN release announcing the release and indicating what will be in the release and a to do list.
2. A CRAN release should not come with major changes that have not been covered by any of the GitHub-only releases.
3. When everything looks ready for the release, repeat the same process as for the GitHub release with a few differences:
    - In the change log it should now say "release date" rather than "GitHub-only-release date".
    - Second digit in the version number is incremented by 1 relative to the current CRAN version.
    - Check whether a new R version has been released or is coming up and make sure GGIR is also tested with that version.
    - Run in RStudio `devtools::check( manual = TRUE, remote = TRUE, incoming = TRUE)` which will help to check urls
4. Ask Vincent (GitHub tag: vincentvanhees) to submit the release to CRAN as it needs to come from my e-mail address.
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
education, socio-economic status, nationality, personal appearance, race,
religion, or sexual identity and orientation.

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
reported by contacting the main project maintainer v.vanhees@accelting.com. All
complaints will be reviewed and investigated and will result in a response that
is deemed necessary and appropriate to the circumstances. The project team is
obligated to maintain confidentiality with regard to the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant](https://www.contributor-covenant.org), version 1.4,
available at https://www.contributor-covenant.org/version/1/4/code-of-conduct.html
# Contributing guidelines

We welcome any kind of contribution to our software, from simple comment or question to a full fledged [pull request](https://help.github.com/articles/about-pull-requests/). Please read and follow our [Code of Conduct](CODE_OF_CONDUCT.md).

A contribution can be one of the following cases:

1. [you have a question](#questions);
2. [you think you may have found a bug](#bugs) (including unexpected behavior);
3. [you want to make some kind of change to the code base](#changes-or-additions) (e.g. to fix a bug, to add a new feature, to update documentation);
4. [you want to make a new release of the code base](#new-release).

The sections below outline the steps in each case.

## Questions

1. use the search functionality [here](https://groups.google.com/g/RpackageGGIR) to see if someone already experienced the same issue;
2. if your search did not yield any relevant results, start a new conversation.

## Bugs

1. use the search functionality [here](https://github.com/wadpac/GGIR/issues) to see if someone already filed the same issue;
2. if your issue search did not yield any relevant results, make a new issue, and choose the Bug report type. This includes a checklist to make sure you provide enough information to the rest of the community to understand the cause and context of the problem.

## Changes or additions

1. (**important**) announce your plan to the rest of the community *before you start working*. This announcement should be in the form of a (new) issue. Choose the Feature request type, which includes a checklist of things to consider to get the discussion going;
2. (**important**) wait until some kind of consensus is reached about your idea being a good idea;
3. if needed, fork the repository to your own Github profile and create your own feature branch off of the latest master commit. While working on your feature branch, make sure to stay up to date with the master branch by pulling in changes, possibly from the 'upstream' repository (follow the instructions [here](https://help.github.com/articles/configuring-a-remote-for-a-fork/) and [here](https://help.github.com/articles/syncing-a-fork/));
4. make sure the existing tests still work by running the test suite from RStudio;
5. add your own tests (if necessary);
6. update or expand the documentation;
7. make sure the release notes in `inst/NEWS.Rd` are updated;
8. add your name to the contributors lists in the `DESCRIPTION` and `CITATION.cff` files;
9. push your feature branch to (your fork of) the GGIR repository on GitHub;
10. create the pull request, e.g. following the instructions [here](https://help.github.com/articles/creating-a-pull-request/). The pull request template includes a checklist with the above items.

In case you feel like you've made a valuable contribution, but you don't know how to write or run tests for it, or how to generate the documentation: don't let this discourage you from making the pull request; we can help you! Just go ahead and submit the pull request, but keep in mind that you might be asked to append additional commits to your pull request.

### Coding style

We loosely follow the [tidyverse style guide](https://style.tidyverse.org/), but do not enforce every rule strictly.
For instance, we prefer `=` instead of `<-` as the default assignment operator.
When in doubt about what style to use, don't hesitate to get in touch.

Some general guidelines that we try to adhere to:

- Use standard R as much as possible, keep dependencies to a minimum.
- Keep loops to a minimum.
- Don't make lines too long.

If you are a first time contributor, don't worry about coding style too much.
We will help you get things in shape.

## New release

GGIR follows the [release cycle process described in this document](RELEASE_CYCLE.md).<!-- Describe your PR here -->

<!-- Please, make sure the following items are checked -->
Checklist before merging:

- [ ] Existing tests still work (check by running the test suite, e.g. from RStudio).
- [ ] Added tests (if you added functionality) or fixed existing test (if you fixed a bug).
- [ ] Updated or expanded the documentation.
- [ ] Updated release notes in `inst/NEWS.Rd` with a user-readable summary. Please, include references to relevant issues or PR discussions.
- [ ] Added your name to the contributors lists in the `DESCRIPTION` and `CITATION.cff` files.
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A short, clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior.

1. Sensor brand: '...'
2. Data format: '...'
3. Approximate recording duration '...' days
4. Are you using a sleep diary to guide the sleep detection: YES / NO
5. Copy of R command used: '....'
6. Have you tried processing your data based on GGIR's default argument values? Does the issue you report still appear? YES / NO
7.

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem. Note that usually we are not only interested in see the error message in red, but all GGIR output to the console.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - GGIR Version [e.g. 2.2-0]

**Additional context**
Add any other context about the problem here.

**Before submitting**
- [ ] Have you tried the steps to reproduce? Do they include all relevant data and configuration? Does the issue you report still appear there?
- [ ] Have you tried this on the latest `master` branch from GitHub?
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
