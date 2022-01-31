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
---
title: "mauricer demo"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{mauricer demo}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Introduction

`mauricer` is an R package to handle BEAST2 package management from R,
similar to the BEAST2 package manager.

![The BEAST2 package manager](beast2_pkg_manager.png)

In this vignette, we will:

 * List all BEAST2 packages
 * Install the first non-installed BEAST2 package
 * Verify the BEAST2 package is indeed installed
 * Uninstall the BEAST2 package

To do so, we'll need to load `mauricer`. We'll also load `testthat` for testing:

```{r}
library(mauricer)
library(testthat)
```

To install BEAST2 packages, one needs

 * an internet connection
 * BEAST2 must be installed

If there is no internet connection, or if BEAST2 is not installed,
this vignette will be close to empty.

```{r}
if (!curl::has_internet()) {
  print("No internet connection")
}
if (!beastier::is_beast2_installed()) {
  print("No BEAST2 installed")
}
```

## List all BEAST2 packages

Use `get_beast2_pkg_names` to get a data frame with the name,
version and install status of all BEAST2 packages:

```{r}
if (curl::has_internet() && beastier::is_beast2_installed()) {
  beast2_packages <- get_beast2_pkg_names()
  knitr::kable(head(beast2_packages))
}
```

## Install the first non-installed BEAST2 package

Find a package that is not installed:

```{r}
if (curl::has_internet() && beastier::is_beast2_installed()) {
  package_name <- beast2_packages[
    beast2_packages$installed_version == "NA",
  ]$name[1]
  print(package_name)
}
```

Install that package:

```
if (curl::has_internet() && beastier::is_beast2_installed()) {
  mauricerinstall::install_beast2_pkg(package_name)
}
```

Note that this installation uses a non-CRAN package,
called `mauricerinstall`, as installing software
is against CRAN policy.

## Uninstall the BEAST2 package

Use `uninstall_beast2_pkg`:

```
if (curl::has_internet() && beastier::is_beast2_installed()) {
  mauricerinstall::uninstall_beast2_pkg(package_name)
}
```

Also this code uses a non-CRAN package,
called `mauricerinstall`, as installing software
is against CRAN policy.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_path.R
\name{get_mrc_path}
\alias{get_mrc_path}
\title{Get the full path of a \code{mauricer} file}
\usage{
get_mrc_path(filename)
}
\arguments{
\item{filename}{the file's name, without the path}
}
\value{
the full path of the filename, if and only if
  the file is present. Will stop otherwise.
}
\description{
Get the full path of a file in the \code{inst/extdata} folder.
If there is no \code{mauricer} file, \link{get_mrc_path} will \link{stop}.
}
\examples{
get_mrc_path("anthus_aco_sub.fas")
}
\seealso{
for more files, use \code{\link{get_mrc_paths}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_paths.R
\name{get_mrc_paths}
\alias{get_mrc_paths}
\title{Get the full path of one or more \code{mauricer} files}
\usage{
get_mrc_paths(filenames)
}
\arguments{
\item{filenames}{the files' names, without the path}
}
\value{
the filenames' full paths
}
\description{
Get the full paths of files in the \code{inst/extdata} folder
If there is a \code{mauricer} file absent,
\link{get_mrc_paths} will \link{stop}.
}
\examples{
get_mrc_paths(c("anthus_aco_sub.fas", "anthus_nd2_sub.fas"))
}
\seealso{
for one file, use \code{\link{get_mrc_path}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/install_beast2_pkg.R
\name{install_beast2_pkg}
\alias{install_beast2_pkg}
\title{Install a BEAST2 package}
\usage{
install_beast2_pkg(
  name,
  beast2_folder = beastier::get_default_beast2_folder(),
  verbose = FALSE,
  has_internet = curl::has_internet()
)
}
\arguments{
\item{name}{the package's name}

\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}

\item{verbose}{set to TRUE for extra output, as can be used in debugging}

\item{has_internet}{boolean to indicate if the user has access to the
internet. By default, this value equals the result
of \code{curl::has_internet}}
}
\value{
nothing. It does install the BEAST2 package
}
\description{
Install a BEAST2 package. If the package is already installed,
(see \link{is_beast2_pkg_installed}), this function \link{stop}s.
}
\note{
Installing or uninstalling a BEAST2 package for a
(singular) BEAST2 installation, does so for all BEAST2
installations
}
\examples{
\dontrun{
  install_beast2_pkg("NS")
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_params_doc.R
\name{default_params_doc}
\alias{default_params_doc}
\title{This function does nothing. It is intended to inherit is parameters'
documentation.}
\usage{
default_params_doc(beast2_folder, has_internet, name, show_warnings, verbose)
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}

\item{has_internet}{boolean to indicate if the user has access to the
internet. By default, this value equals the result
of \code{curl::has_internet}}

\item{name}{the package's name}

\item{show_warnings}{set to TRUE to show warnings}

\item{verbose}{set to TRUE for extra output, as can be used in debugging}
}
\description{
This function does nothing. It is intended to inherit is parameters'
documentation.
}
\note{
This is an internal function, so it should be marked with
  \code{@noRd}. This is not done, as this will disallow all
  functions to find the documentation parameters
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beast2_pkg_names.R
\name{get_beast2_pkg_names}
\alias{get_beast2_pkg_names}
\title{Get all BEAST2 package names}
\usage{
get_beast2_pkg_names(
  beast2_folder = beastier::get_default_beast2_folder(),
  has_internet = curl::has_internet(),
  verbose = FALSE
)
}
\arguments{
\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}

\item{has_internet}{boolean to indicate if the user has access to the
internet. By default, this value equals the result
of \code{curl::has_internet}}

\item{verbose}{set to TRUE for extra output, as can be used in debugging}
}
\value{
a data frame with columns\cr
\enumerate{
  \item name package name, for example, code{bdmm}
  \item installed_version the installed version, for example, \code{2.6.2}.
    \code{installed_version} will be NA if the package is not installed
  \item latest_version version number of the latest version, for example,
    \code{2.6.3}
  \item dependencies packages the package depends on, for example
    \code{BEASTLabs, GEO_SPHERE}. \code{dependencies} will be empty if there
    are no dependencies
  \item description description of the package, for example
    \code{Nested sampling for model selection and posterior inference}
}
}
\description{
List all BEAST2 packages that are available and installed.
Will \link{stop} if there is no internet connection
}
\examples{
if (is_beast2_installed() && curl::has_internet()) {
  get_beast2_pkg_names()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uninstall_beast2_pkg.R
\name{uninstall_beast2_pkg}
\alias{uninstall_beast2_pkg}
\title{Uninstall a BEAST2 package}
\usage{
uninstall_beast2_pkg(
  name,
  beast2_folder = beastier::get_default_beast2_folder(),
  verbose = FALSE,
  has_internet = curl::has_internet()
)
}
\arguments{
\item{name}{the package's name}

\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}

\item{verbose}{set to TRUE for extra output, as can be used in debugging}

\item{has_internet}{boolean to indicate if the user has access to the
internet. By default, this value equals the result
of \code{curl::has_internet}}
}
\value{
nothing. It does install the BEAST2 package
}
\description{
Uninstall a BEAST2 package
}
\note{
Installing or uninstalling a BEAST2 package for a
(singular) BEAST2 installation, does so for all BEAST2
installations
}
\examples{
\dontrun{
  uninstall_beast2_pkg("Beasy")
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mauricer.R
\docType{package}
\name{mauricer}
\alias{mauricer}
\title{mauricer: Install BEAST2 packages.}
\description{
'BEAST2' (<http://www.beast2.org>) is a widely used
Bayesian phylogenetic tool, that uses DNA/RNA/protein data
and many model priors to create a posterior of jointly estimated
phylogenies and parameters.
}
\details{
'BEAST2' is commonly accompanied by 'BEAUti 2' (<http://www.beast2.org>),
which, among others, allows one to install 'BEAST2' package.
This package allows to install 'BEAST2' packages from 'R'.

* \link{get_beast2_pkg_names}: get the names of all BEAST2 packages
* \link{install_beast2_pkg}: install a BEAST2 package
* \link{is_beast2_pkg_installed}: determine if a BEAST2 package is installed
* \link{uninstall_beast2_pkg}: uninstall a BEAST2 package
}
\seealso{
\link{mauricer} is part of the \code{babette} package suite:
\itemize{
  \item{
    \code{babette}: work with BEAST2
  }
  \item{
    \link[beautier]{beautier}: create BEAST2 input files
  }
  \item{
    \link[beastier]{beastier}: run BEAST2
  }
  \item{
    \link{mauricer}: install BEAST2 packages
  }
  \item{
    \link[tracerer]{tracerer}: parse and analyse BEAST2 output
  }
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_beast2_pkg_installed.R
\name{is_beast2_pkg_installed}
\alias{is_beast2_pkg_installed}
\title{Is a BEAST2 package installed?}
\usage{
is_beast2_pkg_installed(
  name,
  beast2_folder = beastier::get_default_beast2_folder(),
  has_internet = curl::has_internet()
)
}
\arguments{
\item{name}{the package's name}

\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}

\item{has_internet}{boolean to indicate if the user has access to the
internet. By default, this value equals the result
of \code{curl::has_internet}}
}
\value{
\itemize{
  \item \code{TRUE} if the BEAST2 package is installed
  \item \code{FALSE} if the BEAST2 package is not installed
  \item \code{NULL} if there is no internet connection
}
}
\description{
Checks if a BEAST2 package is installed.
}
\details{
To be able to check this, an internet connection is needed.
Without an internet connection, \code{NULL} is returned.
}
\examples{
\dontrun{
  is_beast2_pkg_installed("Beasy")
}
}
\seealso{
Use \link{is_beast2_ns_pkg_installed}
  to see if the NS package is installed without an internet connection
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_beast2_ns_pkg_installed.R
\name{is_beast2_ns_pkg_installed}
\alias{is_beast2_ns_pkg_installed}
\title{Is the BEAST2 NS package installed?}
\usage{
is_beast2_ns_pkg_installed(
  show_warnings = FALSE,
  verbose = FALSE,
  beast2_folder = beastier::get_default_beast2_folder()
)
}
\arguments{
\item{show_warnings}{set to TRUE to show warnings}

\item{verbose}{set to TRUE for extra output, as can be used in debugging}

\item{beast2_folder}{the folder where the BEAST2 is installed.
Note that this is not the folder where the BEAST2 executable is installed:
the BEAST2 executable is in a subfolder.
Use \link[beastier]{get_default_beast2_folder}
  to get the default BEAST2 folder.
Use \link[beastier]{get_default_beast2_bin_path}
  to get the full path to the default BEAST2 executable.
Use \link[beastier]{get_default_beast2_jar_path}
  to get the full path to the default BEAST2 jar file.}
}
\value{
TRUE if the BEAST2 NS package is installed, FALSE otherwise
}
\description{
Determine if the BEAST2 NS package is installed.
}
\details{
Unlike \link{is_beast2_pkg_installed},
this function does not need an internet connection.
Instead, the function calls BEAST2 to read a BEAST2 XML file that
uses NS.
}
\examples{
\dontrun{
  is_beast2_ns_pkg_installed()
}
}
\author{
Richèl J.C. Bilderbeek
}
