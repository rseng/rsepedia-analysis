# mauricer

[![Peer Review Status](https://badges.ropensci.org/209_status.svg)](https://github.com/ropensci/onboarding/issues/209)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/mauricer)](https://cran.r-project.org/package=mauricer)
[![](http://cranlogs.r-pkg.org/badges/grand-total/mauricer)]( https://CRAN.R-project.org/package=mauricer)
[![](http://cranlogs.r-pkg.org/badges/mauricer)](https://CRAN.R-project.org/package=mauricer)

Branch   |[![GitHub Actions logo](man/figures/GitHubActions.png)](https://github.com/ropensci/mauricer/actions)|[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
---------|-----------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------
`master` |![R-CMD-check](https://github.com/ropensci/mauricer/workflows/R-CMD-check/badge.svg?branch=master)   |[![codecov.io](https://codecov.io/github/ropensci/mauricer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/mauricer/branch/master)
`develop`|![R-CMD-check](https://github.com/ropensci/mauricer/workflows/R-CMD-check/badge.svg?branch=develop)  |[![codecov.io](https://codecov.io/github/ropensci/mauricer/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/mauricer/branch/develop)

Work with BEAST2 packages from R.

Related packages:

 * [babette](https://github.com/ropensci/babette) do a full BEAST2 workflow.
 * [beautier](https://github.com/ropensci/beautier) creates BEAST2 input (`.xml`) files.
 * [beastier](https://github.com/ropensci/beastier) runs BEAST2
 * [lumier](https://github.com/ropensci/lumier) helps to create the `babette` function call needed
 * [tracerer](https://github.com/ropensci/tracerer) parses BEAST2 output (`.log`, `.trees`, etc) files.

Non-CRAN extensions:

 * [beastierinstall](https://github.com/richelbilderbeek/beastierinstall) Install and uninstall BEAST2
 * [mauricerinstall](https://github.com/richelbilderbeek/mauricerinstall) Install and uninstall BEAST2 packages

## Examples

To install the BEAST2 NS package:

```
remotes::install_github("richelbilderbeek/mauricerinstall")
mauricerinstall::install_beast2_pkg("NS")
```

An introduction video:

 * [YouTube video about mauricer](https://youtu.be/Yk737gorcrw)


## Package dependencies

Package                                                                     |[![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.com)                                                 |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
----------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------
[beastier](https://github.com/ropensci/beastier)                            |[![Build Status](https://travis-ci.com/ropensci/beastier.svg?branch=master)](https://travis-ci.com/ropensci/beastier)|[![codecov.io](https://codecov.io/github/ropensci/beastier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/beastier/branch/master)

## Related packages

Package                                         |[![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.com)                                                 |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------
[babette](https://github.com/ropensci/babette)  |[![Build Status](https://travis-ci.com/ropensci/babette.svg?branch=master)](https://travis-ci.com/ropensci/babette)  |[![codecov.io](https://codecov.io/github/ropensci/babette/coverage.svg?branch=master)](https://codecov.io/github/ropensci/babette/branch/master)
[beautier](https://github.com/ropensci/beautier)|[![Build Status](https://travis-ci.com/ropensci/beautier.svg?branch=master)](https://travis-ci.com/ropensci/beautier)|[![codecov.io](https://codecov.io/github/ropensci/beautier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/beautier/branch/master)
[lumier](https://github.com/ropensci/lumier)    |[![Build Status](https://travis-ci.com/ropensci/lumier.svg?branch=master)](https://travis-ci.com/ropensci/lumier)    |[![codecov.io](https://codecov.io/github/ropensci/lumier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/lumier/branch/master)
[tracerer](https://github.com/ropensci/tracerer)|[![Build Status](https://travis-ci.com/ropensci/tracerer.svg?branch=master)](https://travis-ci.com/ropensci/tracerer)|[![codecov.io](https://codecov.io/github/ropensci/tracerer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/tracerer/branch/master)

Package                                                                        | Status
-------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
[mauricer_on_windows](https://github.com/richelbilderbeek/mauricer_on_windows) |[![Build status](https://ci.appveyor.com/api/projects/status/bc43iwp68xo2dduh/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/mauricer-on-windows/branch/master)

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

# News

Newest versions at top.

## `mauricer` 2.5.1 (2021-05-29)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Creates no temporary files

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.5 (2021-05-22)

### NEW FEATURES

  * Use the GitHub Actions continuous integration service

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * `install_beast2_pkg` and  `uninstall_beast2_pkg` are deprecated,
    as these violated CRAN policy.
    The deprecation message will point users to the non-official
    `mauricerinstall` package at 
    `https://github.com/richelbilderbeek/mauricerinstall`

## `mauricer` 2.4 (2021-05-14)

### NEW FEATURES

  * Use the GitHub Actions continuous integration service

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.3 (2020-10-31)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * `install_beast2_pkg` and `uninstall_beast2_pkg` can be verbose

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.2 (2020-10-21)

### NEW FEATURES

  * Can install BEAST2 packages when BEAST2 is installed at a custom location

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.1 (2020-08-05)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Use https of BEAST2 website
  * Builds on MacOS again
  * `install_beast2_pkg` gives better error message when a 
    non-existing/misspelled package is requested to be installed

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.0.6 (2020-01-06)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Documentation at [rOpenSci](https://docs.ropensci.org/mauricer) to
    shows pictures

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.0.5 (2019-12-02)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Add vignette

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.0.4 (2019-10-27)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Removed all deprecated functions; those starting with `mrc_`
  * Added `is_beast2_ns_pkg_installed` that specifically checks if the BEAST
    NS package is installed, without needing an internet connection

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.0.3 (2019-06-03)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Depends on `beastier` being on CRAN

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.0.2 (2019-03-27)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Builds with and without BEAST2 installed

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.0.1 (2019-01-04)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Processed CRAN feedback on `beastier` on this repo 

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 2.0 (2019-01-04)

### NEW FEATURES

  * Transferred ownership to `ropensci`

### MINOR IMPROVEMENTS

  * Use more readable function names

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 1.2 (2018-10-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tested to work on macOS
  * Improved DESCRIPTION
  * Improved function names

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 1.1.0 (2018-10-20)

### NEW FEATURES

  * Works on Windows

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## `mauricer` 1.0.0 (2018-09-11)

### NEW FEATURES

  * First version

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None
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
reported by contacting the project team at richel@richelbilderbeek.nl. All
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
# Contributing

Awesome that you are reading this.

This GitHub follows the [Contributor Covenant Code of Conduct](code_of_conduct.md).

 * For questions, you can create an Issue
 * Code changes go via Pull Requests

## Which package to contribute to?

`maurcer` is part of the `babette` package suite,
which consists out of five packages.
Here is how to determine which package is best suited for your contribution:

If you want to contribute to the creation of BEAST2 XML input file,
go to [beautier](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md).

If you want to contribute to how BEAST2 is run,
go to [beautier](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md).

If you want to contribute to how BEAST2 output is parsed,
go to [tracerer](https://github.com/ropensci/tracerer/blob/master/CONTRIBUTING.md)

If you want to contribute how these packages work together,
go to [babette](https://github.com/ropensci/babette/blob/master/CONTRIBUTING.md)

If you want to contribute regarding the BEAST2 package management, 
you are at the right spot :-) 

## Submitting code

Submitted code should follow these quality guidelines:

 * All tests pass cleanly/silently
 * Code coverage must be 100%
 * Coding style should follow the default style by `lintr`

These are all checked by Travis CI when submitting
a Pull Request. 

Emails with code will not be accepted.

## Submitting bugs

Awesome. These are your options:

 * Add an Issue, with the test that fails
 * Submit a Pull Request, where the test is added to the `tests/testthat` folder
 * Send @richelbilderbeek an email (@richelbilderbeek will make an Issue of it)

Pull Requests should follow the same guidelines as 'Submitting code'.

## Branching policy

 * The `master` branch should always build successfully
 * The `development` branch is for developers

## git usage

To get started working on `babette` do:

```
git clone https://github.com/ropensci/babette
```

Development is done on the `develop` branch. 
To download and checkout the `develop` branch, 
first go into the `beautier` folder (`cd babette`), then do:

```
git checkout develop
```

Then the workflow is the common `git` workflow:

```
git pull
git add --all :/
git commit -m "Did something awesome"
git push
```

Hi @richelbilderbeek,

With this Pull Request I'd would like to [add reason].

Sure, I've read [CONTRIBUTING.md](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md) :+1:

Cheers, [your name]

---
name: Custom issue
about: Anything else
title: ''
labels: ''
assignees: ''

---


---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Script to reproduce the behavior:

```r
# Your R script here, without this comment :-)
```

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Environment:**
Show the results of running the following script:

```r
library(beastier)
beastier::beastier_report()
```

**Additional context**
Add any other context about the problem here.

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
