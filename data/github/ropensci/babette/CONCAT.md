# babette

[![Peer Review Status](https://badges.ropensci.org/209_status.svg)](https://github.com/ropensci/onboarding/issues/209)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/babette)](https://cran.r-project.org/package=babette)
[![](http://cranlogs.r-pkg.org/badges/grand-total/babette)]( https://CRAN.R-project.org/package=babette)
[![](http://cranlogs.r-pkg.org/badges/babette)](https://CRAN.R-project.org/package=babette)

Branch   |[![GitHub Actions logo](man/figures/GitHubActions.png)](https://github.com/ropensci/babette/actions)|[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
---------|----------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------
`master` |![R-CMD-check](https://github.com/ropensci/babette/workflows/R-CMD-check/badge.svg?branch=master)   |[![codecov.io](https://codecov.io/github/ropensci/babette/coverage.svg?branch=master)](https://codecov.io/github/ropensci/babette/branch/master)
`develop`|![R-CMD-check](https://github.com/ropensci/babette/workflows/R-CMD-check/badge.svg?branch=develop)  |[![codecov.io](https://codecov.io/github/ropensci/babette/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/babette/branch/develop)

[![DOI](https://zenodo.org/badge/118616108.svg)](https://zenodo.org/badge/latestdoi/118616108)

`babette` is an R package that combines:

 * [beautier](https://github.com/ropensci/beautier) creates a BEAST2 input (`.xml`) file from an inference model
 * [tiebeaur](https://github.com/richelbilderbeek/tiebeaur) creates an inference model from a BEAST2 input (`.xml`) file :warning: experimental :warning:
 * [beastier](https://github.com/ropensci/beastier) runs BEAST2
 * [mauricer](https://github.com/ropensci/mauricer) install BEAST2 packages
 * [tracerer](https://github.com/ropensci/tracerer) parses BEAST2 output (`.log`, `.trees`, etc) files.

![babette logo](man/figures/babette_logo.png)

Non-CRAN extensions:

 * [beastierinstall](https://github.com/richelbilderbeek/beastierinstall) Install and uninstall BEAST2
 * [mauricerinstall](https://github.com/richelbilderbeek/mauricerinstall) Install and uninstall BEAST2 packages

Related R packages:

 * [lumier](https://github.com/ropensci/lumier): Shiny app to help create the function call needed

Related R functions:

Related function                                                      |`babette` function
----------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------------------
`[ips::rbeauti](https://github.com/heibl/ips/blob/master/R/rbeauti.R)`|`[beautier::create_beast2_input_from_model](https://github.com/ropensci/beautier/blob/master/R/create_beast2_input_from_model.R)`

Related software:

 * [BEASTifier](https://github.com/josephwb/BEASTifier): command-line tool to generate BEAST2 XML input files

## Examples

See:

 * [lumier](https://github.com/ropensci/lumier): R Shiny app to help create the R function call needed
 * [examples](https://github.com/richelbilderbeek/babette_examples): examples tested by Travis CI and AppVeyor
 * [article](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.13032), in 'Methods in Ecology and Evolution'
 * [Methods.blog post: The babette R Package: How to Sooth the Phylogenetic BEAST2](https://methodsblog.wordpress.com/2018/06/25/babette-beast2/)
 * [rOpenSci blog post: Call BEAST2 for Bayesian evolutionary analysis from R](https://ropensci.org/blog/2020/01/28/babette/)
 * [pre-print article](https://doi.org/10.1101/271866), in bioRxiv
 * ['babette' YouTube channel](https://www.youtube.com/watch?v=nA-0-Fc95xY&list=PLu8_ZyzXyRDFIRx-kdDI5Q6xVr-HnY7TB)

## Installation

See [doc/install.md](doc/install.md) (or just click [here](doc/install.md))

## FAQ

See [FAQ](doc/faq.md).

## Missing features/unsupported

`babette` cannot do everything `BEAUti` and `BEAST2` and `Tracer` can.

See [beautier](https://github.com/ropensci/beautier) 
for missing features in creating a BEAST2 input (`.xml`) file.

See [beastier](https://github.com/ropensci/beastier) for missing
features in running BEAST2

See [mauricer](https://github.com/ropensci/mauricer) for missing
features in installing BEAST2 packages.

See [tracerer](https://github.com/ropensci/tracerer) 
for missing features in parsing BEAST2 output (`.log`, `.trees`, etc) files.

## There is a feature I miss

See [CONTRIBUTING](CONTRIBUTING.md), at `Submitting use cases`

## I want to collaborate

See [CONTRIBUTING](CONTRIBUTING.md), at `Submitting code`

## I think I have found a bug

See [CONTRIBUTING](CONTRIBUTING.md), at `Submitting bugs` 

## There's something else I want to say

Sure, just add an Issue. Or send an email.

## Dependencies

Branch                                          |[![GitHub Actions logo](man/figures/GitHubActions.png)](https://github.com/ropensci/babette/actions) `master`|[![GitHub Actions logo](man/figures/GitHubActions.png)](https://github.com/ropensci/babette/actions) `develop`|[![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.com) `master`                                        |[![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.com) `develop`                                        |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io) `master`                                                                       |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io) `develop`                                                                        
------------------------------------------------|-------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------|---------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|--------------------------------------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------
[beautier](https://github.com/ropensci/beautier)|![R-CMD-check](https://github.com/ropensci/beautier/workflows/R-CMD-check/badge.svg?branch=master)           |![R-CMD-check](https://github.com/ropensci/beautier/workflows/R-CMD-check/badge.svg?branch=develop)           |[![Build Status](https://travis-ci.com/ropensci/beautier.svg?branch=master)](https://travis-ci.com/ropensci/beautier)|[![Build Status](https://travis-ci.com/ropensci/beautier.svg?branch=develop)](https://travis-ci.com/ropensci/beautier)|[![codecov.io](https://codecov.io/github/ropensci/beautier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/beautier/branch/master)|[![codecov.io](https://codecov.io/github/ropensci/beautier/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/beautier/branch/develop)
[beastier](https://github.com/ropensci/beastier)|![R-CMD-check](https://github.com/ropensci/beastier/workflows/R-CMD-check/badge.svg?branch=master)           |![R-CMD-check](https://github.com/ropensci/beastier/workflows/R-CMD-check/badge.svg?branch=develop)           |[![Build Status](https://travis-ci.com/ropensci/beastier.svg?branch=master)](https://travis-ci.com/ropensci/beastier)|[![Build Status](https://travis-ci.com/ropensci/beastier.svg?branch=develop)](https://travis-ci.com/ropensci/beastier)|[![codecov.io](https://codecov.io/github/ropensci/beastier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/beastier/branch/master)|[![codecov.io](https://codecov.io/github/ropensci/beastier/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/beastier/branch/develop)
[mauricer](https://github.com/ropensci/mauricer)|![R-CMD-check](https://github.com/ropensci/mauricer/workflows/R-CMD-check/badge.svg?branch=master)           |![R-CMD-check](https://github.com/ropensci/mauricer/workflows/R-CMD-check/badge.svg?branch=develop)           |[![Build Status](https://travis-ci.com/ropensci/mauricer.svg?branch=master)](https://travis-ci.com/ropensci/mauricer)|[![Build Status](https://travis-ci.com/ropensci/mauricer.svg?branch=develop)](https://travis-ci.com/ropensci/mauricer)|[![codecov.io](https://codecov.io/github/ropensci/mauricer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/mauricer/branch/master)|[![codecov.io](https://codecov.io/github/ropensci/mauricer/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/mauricer/branch/develop)
[tracerer](https://github.com/ropensci/tracerer)|![R-CMD-check](https://github.com/ropensci/tracerer/workflows/R-CMD-check/badge.svg?branch=master)           |![R-CMD-check](https://github.com/ropensci/tracerer/workflows/R-CMD-check/badge.svg?branch=develop)           |[![Build Status](https://travis-ci.com/ropensci/tracerer.svg?branch=master)](https://travis-ci.com/ropensci/tracerer)|[![Build Status](https://travis-ci.com/ropensci/tracerer.svg?branch=develop)](https://travis-ci.com/ropensci/tracerer)|[![codecov.io](https://codecov.io/github/ropensci/tracerer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/tracerer/branch/master)|[![codecov.io](https://codecov.io/github/ropensci/tracerer/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/tracerer/branch/develop)

## Windows

Package                                                                       | Status
------------------------------------------------------------------------------|------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------
[babette_on_windows](https://github.com/richelbilderbeek/babette_on_windows)  |[![Build status](https://ci.appveyor.com/api/projects/status/jv76errjocm5d5yq/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/babette-on-windows/branch/master)
[beastier_on_windows](https://github.com/richelbilderbeek/beastier_on_windows)|[![Build status](https://ci.appveyor.com/api/projects/status/ralex9sdnnxlwbgx/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/beastier-on-windows/branch/master)
[beautier_on_windows](https://github.com/richelbilderbeek/beautier_on_windows)|[![Build status](https://ci.appveyor.com/api/projects/status/blvjo5pulbkqxrhb/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/beautier-on-windows/branch/master)
[mauricer_on_windows](https://github.com/richelbilderbeek/mauricer_on_windows)|[![Build status](https://ci.appveyor.com/api/projects/status/bc43iwp68xo2dduh/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/mauricer-on-windows/branch/master)
[tracerer_on_windows](https://github.com/richelbilderbeek/tracerer_on_windows)|[![Build status](https://ci.appveyor.com/api/projects/status/jyhck66d6yrbr12h/branch/master?svg=true)](https://ci.appveyor.com/project/richelbilderbeek/tracerer-on-windows/branch/master)

## Related packages

Package                                     |[![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.com)                                             |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
--------------------------------------------|-----------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------
[lumier](https://github.com/ropensci/lumier)|[![Build Status](https://travis-ci.com/ropensci/lumier.svg?branch=master)](https://travis-ci.com/ropensci/lumier)|[![codecov.io](https://codecov.io/github/ropensci/lumier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/lumier/branch/master)

## External links

 * [BEAST2 GitHub](https://github.com/CompEvol/beast2)

## References

Article about `babette`:

 * Bilderbeek, Richèl JC, and Rampal S. Etienne. "babette: BEAUti 2, BEAST 2 and Tracer for R." Methods in Ecology and Evolution (2018). https://doi.org/10.1111/2041-210X.13032

```
@article{bilderbeek2018babette,
  title={babette: BEAUti 2, BEAST 2 and Tracer for R},
  author={Bilderbeek, Richèl JC and Etienne, Rampal S},
  journal={Methods in Ecology and Evolution},
  year={2018},
  publisher={Wiley Online Library}
}
```

FASTA files `anthus_aco.fas` and `anthus_nd2.fas` from:
 
 * Van Els, Paul, and Heraldo V. Norambuena. "A revision of species limits in Neotropical pipits Anthus based on multilocus genetic and vocal data." Ibis.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)


# News

Newest versions at top.

## babette 2.3 (2021-05-14)

### NEW FEATURES

  * Use GitHub Actions as a continuous integration service
  
### MINOR IMPROVEMENTS

  * Add URL to DESCRIPTION

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 2.2 (2020-11-06)

### NEW FEATURES

  * None
  
### MINOR IMPROVEMENTS

  * Do not run example code of internal function 'prepare_file_creation'
    on CRAN flavor `r-oldrel-macos-x86_64`

### BUG FIXES

  * Can run BEAST2 and BEAST2 packages when these are installed at
    a custom location

### DEPRECATED AND DEFUNCT

  * None

## babette 2.1.4 (2020-08-11)

### NEW FEATURES

  * CRAN version
  
### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 2.1.3 (2020-08-05)

### NEW FEATURES

  * None
  
### MINOR IMPROVEMENTS

  * No testthat usage in documentation examples

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 2.1.2 (2020-02-01)

### NEW FEATURES

  * None
  
### MINOR IMPROVEMENTS

  * Use https of BEAST2 website
  * Building the package does not create temporary files in `~/.cache`, fix #85

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 2.1.1 (2019-12-30)

### NEW FEATURES

  * None
  
### MINOR IMPROVEMENTS

  * Process feedback CRAN (Issue 81)
  * Renamed internal function names to be more descriptive

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * Remove arguments that are deprecated, as there is no CRAN
    version yet

## babette 2.1 (2019-10-27)

### NEW FEATURES

  * Follow the `beastier` v2.3 new interface, in which the log filenames
    are part of the MCMC (as created by `create_mcmc`)
  
### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * In `create_beast2_options`, 
    `output_log_filename` and `output_trees_filenames` are removed
    as arguments. These values are set by the BEAST2 XML

## babette 2.0.3 (2019-08-26)

### NEW FEATURES

  * babette can run when the BEAST2 options request temporary files
    in sub-sub-subfolders

### MINOR IMPROVEMENTS

  * Builds on Travis CI's Bionic distribution

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 2.0.2 (2019-03-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Fix spelling error

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 2.0.1 (2019-03-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Travis CI uses `phangorn`s binary R package  

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 2.0 (2019-01-06)

### NEW FEATURES

  * Transferred GitHub repository ownership to `ropensci`
  * CRAN and rOpenSci release candidate

### MINOR IMPROVEMENTS

  * Tagged for [the academic article about babette](https://github.com/ropensci/babette_article)

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None


## babette 1.2.5 (2018-10-25)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * BEAST2 package management using `mauricer` works under Windows as well
  * `babette` indicates that BEAST2.exe cannot be run under Windows
  * Simplified interface for parameters:

```
# Old
distr <- create_distr_poisson(id = 1, lambda = create_lambda_param(value = 1.2))
# Added
distr <- create_distr_poisson(id = 1, lambda = 1.2)
```

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * Cannot use hyper parameters

```
# Deprecated: create a hyper parameter
create_lambda_param(value = 1.2, estimate = TRUE)
```

## babette 1.2.4 (2018-05-17)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tagged for [the academic article about babette](https://github.com/ropensci/babette_article)

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 1.2.3 (2018-05-04)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Renamed `run` to `bbt_run` to reduce conflicts with other packages

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## babette 1.2.2 (2018-04-05)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Follow all [rOpenSci packaging guidelines](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md)

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

This GitHub follows the [Contributor Covenant Code of Conduct](doc/code_of_conduct.md).

 * For questions, you can create an Issue
 * Code changes go via Pull Requests

## Which package to contribute to?

`babette` is part of the `babette` package suite,
which consists out of five packages.
Here is how to determine which package is best suited for your contribution:

If you want to contribute to the creation of BEAST2 XML input file,
go to [beautier](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md).

If you want to contribute to how BEAST2 is run,
go to [beautier](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md).

If you want to contribute to how BEAST2 output is parsed,
go to [tracerer](https://github.com/ropensci/tracerer/blob/master/CONTRIBUTING.md)

If you want to contribute regarding the BEAST2 package management,
go to [mauricer](https://github.com/ropensci/mauricer/blob/master/CONTRIBUTING.md)

If you want to contribute with an overarching idea, you are at the right spot :-) 

## Submitting code

Submitted code should follow these quality guidelines:

 * All tests pass cleanly/silently
 * Code coverage must be 100%
 * Coding style should follow the default style by `lintr`

These are all checked by Travis CI when submitting
a Pull Request. 

Emails with code will not be accepted.

## Submitting bugs

Awesome that your are reading this!

First, know that `babette` is thouroughly tested. 
Most reported bugs are due to users' local settings,
caused by `babette`'s creation of temporary files.
This may cause a problem if:

 * (Linux) the home folder is encrypted
 * (Linux) [CentOS is used](https://github.com/ropensci/babette/issues/31).
 * (MacOS) if `babette` is run in a folder monitored by DropBox

If this is not the case, these are your options:

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
# `pics`

Pictures for `babette`.

## Attribution

### Swan

![Swan](swan.png)

The swan is drawn by Jose Scholte, who kindly allowed her work to be used for free, by attribution.

### `babette` logo

The logo is created by Richèl J.C. Bilderbeek, by modifying the swan picture and merging it with the R logo.

![babette logo](babette_logo.png)

### Orange-fronted barbet 

![Orange-fronted barbet](Capitosquamatus.jpg)

 * By John Gerrard Keulemans - The Ibis, ser. 3, vol. 6, Public Domain, https://commons.wikimedia.org/w/index.php?curid=14887149

## Convert the fuzzy white background of swan picture to one single color

```
convert swan.png -fuzz 15% -fill white -opaque white swan_mono_background.png
convert swan_mono_background.png -background white -alpha remove swan_mono_background_2.png
```
# Install

This package describes the `babette` installation and installation problems.

## Install `rJava`

If you have problems installing rJava, [Duck](http://www.duckduckgo.com) 
or [view my rJava notes](rjava.md).

## `babette` installation

There are two types of `babette` installations:

 * Stable (recommended)
 * Bleeding edge

`babette` assumes that BEAST2 is installed. 
See below at `Install BEAST2` how to install BEAST2.

### Stable

`babette` is on CRAN. 

```
install.packages("babette")
```

### Bleeding endge

 * Video how to install `babette`: [download (.mkv)](http://richelbilderbeek.nl/babette_install_windows.mkv) [YouTube](https://youtu.be/SiJlssZeeaM)

`babette` is installed most easily using the `remotes` package:

```
remotes::install_github("ropensci/beautier")
remotes::install_github("ropensci/tracerer")
remotes::install_github("ropensci/beastier")
remotes::install_github("ropensci/mauricer")
remotes::install_github("ropensci/babette")

```

## Install BEAST2

`babette` assumes that BEAST2 is installed. 

If not, one can install BEAST2 from R using `beastierinstall`:

To install BEAST2:

```
remotes::install_github("richelbilderbeek/beastierinstall")
beastierinstall::install_beast2()
```

To install a BEAST2 package, for example, the `NS` package:

```
remotes::install_github("richelbilderbeek/mauricerinstall")
mauricerinstall::install_beast2_pkg("NS")
```

# `doc`

The `babette` documentation.

## Installation

See [install.md](install.md).

# Examples

See [the babette examples](https://github.com/richelbilderbeek/babette_examples).

## Create the dependency graph from the `.dot` file

```
dot -Tps dependencies.dot -o dependencies.ps
convert dependencies.ps dependencies.png
```


# `rJava` questions

How to install `rJava` under different operating systems

 * Installation
 * Troubleshooting

## Installation

### Mac OS

On Travis, I do:

```
R CMD javareconf
R --quiet -e 'install.packages("rJava", type="source", repos="http://cran.us.r-project.org")'
```

May work for you as well.

### Ubuntu 14.5 (Trusty Tahr)

The `.travis.yml` file shows a Trusty install:

```
 - sudo apt-get install -qq oracle-java8-installer # Java 8
 - sudo apt-get install oracle-java8-set-default
```



So I assume the same can be achieved with:

```
sudo add-apt-repository -y ppa:webupd8team/java 
sudo apt-get update -qq
sudo apt-get install oracle-java8-installer
sudo apt-get install oracle-java8-set-default
```

### Ubuntu 17.10 (Artful Aardvark)

```
sudo apt-get install r-cran-rjava openjdk-8-jdk
R CMD javareconf
```

Do not use `openjdk-9-jdk`.

### Ubuntu 18.4 (Bionic Beaver)

The `.travis.yml` file shows a Trusty install:

```
# - sudo apt install -qq oracle-java8-installer # Java 8
# - sudo apt install oracle-java8-set-default
```

On Bionic, I achieved the same with approx:

```
sudo add-apt-repository ppa:marutter/c2d4u3.5
sudo apt update
sudo apt grade
sudo apt install r-cran-rjava
sudo apt-get install openjdk-11-jdk
sudo R CMD javareconf
```

```
#sudo add-apt-repository -y ppa:webupd8team/java 
#sudo apt-get update -qq
#sudo apt-get install oracle-java8-installer
#sudo apt-get install oracle-java8-set-default
```


## Troubleshooting

### Error: `libjvm.so: cannot open shared object file: No such file or directory`

#### Random solution 1

Sometimes works:

For me, [this Stack Overflow post](https://stackoverflow.com/a/25932828) helped me out:

```
sudo mousepad /etc/ld.so.conf.d/java.conf
```

In that file put:

```
/usr/lib/jvm/java-8-oracle/jre/lib/amd64
/usr/lib/jvm/java-8-oracle/jre/lib/amd64/server
```

Save, close, restart R studio, fixed!

#### Random solution 2

Random notes:

Else [this Stack Overflow post may be helpful](https://stackoverflow.com/a/43466434):

Ruthlessly install all JDK stuff:

```
sudo apt-get install jdk-*
```

```
sudo R CMD javareconf
```

```
sudo R CMD javareconf -e
export LD_LIBRARY_PATH=$JAVA_LD_LIBRARY_PATH
sudo apt-get install r-cran-rjava
```

### BEAST2 cannot find Java

![BEAST2 cannot find Java](beast_cannot_find_java.png)


Download the Oracle Java SDK:

![](download_oracle_java_sdk.png)

Open the Oracle Java SDK with the package installer:

![](open_oracle_java_sdk.png)

Install the Oracle Java SDK with the package installer:

![](install_oracle_java_sdk.png)

Pick the right `java`:


```
sudo update-alternatives --config java
```

I picked:

```
There are 5 choices for the alternative java (providing /usr/bin/java).

  Selection    Path                                            Priority   Status
------------------------------------------------------------
  0            /usr/lib/jvm/java-9-openjdk-amd64/bin/java       1091      auto mode
  1            /usr/bin/gij-4.8                                 1048      manual mode
  2            /usr/bin/gij-5                                   1050      manual mode
* 3            /usr/bin/gij-6                                   1060      manual mode
  4            /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java   1081      manual mode
  5            /usr/lib/jvm/java-9-openjdk-amd64/bin/java       1091      manual mode

Press <enter> to keep the current choice[*], or type selection number: 
```

Reconfig:

```
sudo R CMD javareconf
```
# FAQ

## BEAST2

### Which version of BEAUti do you use as a guideline?

Run this code:

```
beastier::get_beast2_version()
```

Currently, this returns `2.6.4` which means that the
BEAUti version `2.6` (the patch version is irrelevant)
is ignored.

### How to install BEAST2?

See [install.md](install.md)

## `babette` development 

### What's the road map?

Currently, `babette` does not have a road map itself, but `beautier` does:

 * [beautier road map](https://github.com/richelbilderbeek/beautier/blob/master/road_map.md)

### How can I indicate a feature that I miss?

Submit an Issue.

### How can I submit code?

See [CONTRIBUTING](../CONTRIBUTING.md), at 'Submitting code'

### How can I submit a bug?

See [CONTRIBUTING](../CONTRIBUTING.md), at 'Submitting bugs' 

### How can I indicate something else?

Submit an Issue. Or send an email to Richèl Bilderbeek.

### What are the `babette` dependencies?

![babette dependencies](dependencies.png)

## `babette` technical questions

### How can I inspect a generated BEAST2 `.xml` file?

The path to the BEAST2 `.xml` file, 
which is the file that contains the inference model,
can be found in the BEAST2 options.
A typical user is unaware these options exists,
as these are generated by default.

Here is an example that calls `bbt_run_from_model` with
an explicit BEAST2 options:

```
# Default BEAST2 options
beast2_options <- create_beast2_options()

# The path where the BEAST2 XML file will be created:
beast2_options$input_filename

bbt_run_from_model(
  fasta_filename = get_babette_path("anthus_aco.fas"),
  beast2_options = beast2_options
)

# Show BEAST2 XML file
readLines(beast2_options$input_filename)
```

### How can I inspect a generated BEAST2 `.xml.state` file?

The path to the BEAST2 `.xml.state` file, 
which is the file that contains the current/final state of the MCMC machinery,
can be found in the BEAST2 options.
A typical user is unaware these options exists,
as these are generated by default.

A typical user is unaware these options exists,
as these are generated by default.

Here is an example that calls `bbt_run_from_model` with
an explicit BEAST2 options:

```
# Default BEAST2 options
beast2_options <- create_beast2_options()

# The path where the final BEAST2 state file will be created:
beast2_options$output_state_filename

bbt_run_from_model(
  fasta_filename = get_babette_path("anthus_aco.fas"),
  beast2_options = beast2_options
)

# Show BEAST2 state file
readLines(beast2_options$output_state_filename)
```

### How can I inspect a generated BEAST2 `.trees` file?

This `.trees` file, which holds the estimated trees,
is part of the `inference_model`:

```
inference_model <- create_test_inference_model()

# Path to .trees file
inference_model$mcmc$treelog$filename

bbt_run_from_model(
  fasta_filename = get_babette_path("anthus_aco.fas"),
  inference_model = inference_model
)

# Show the generated .trees file
readLines(inference_model$mcmc$tracelog$filename)
```

### How can I inspect a generated BEAST2 `.log` file?

This `.log` file holds the parameter estimations
is part of the `inference_model`:

```
inference_model <- create_test_inference_model()

# Path to .log file
inference_model$mcmc$tracelog$filename

bbt_run_from_model(
  fasta_filename = get_babette_path("anthus_aco.fas"),
  inference_model = inference_model
)

# Show the generated .log file
readLines(inference_model$mcmc$tracelog$filename)
```

### Can `babette` handle missing data in the FASTA files?

Yes. BEAST2 can missing data in the FASTA files (which
contain the DNA/RNA or protein sequences) and so does `babette`.

In your FASTA file, use a dash `-` (not NA)
for the sequences that are missing.
This is similar to what the `phangorn` package does.

### If I set a fixed crown age with multiple alignments, only the first alignment has so

Correct. This is a feature of BEAST2. Using the `create_mrca` prior 
gives prettier results.

## `babette` in academia

### How do I reference to this work?

```
Bilderbeek, Richèl JC, and Rampal S. Etienne. "babette: BEAUti 2, BEAST 2 and Tracer for R." Methods in Ecology and Evolution (2018).
```

or

```
@article{bilderbeek2018babette,
  title={babette: BEAUti 2, BEAST 2 and Tracer for R},
  author={Bilderbeek, Richel JC and Etienne, Rampal S},
  journal={Methods in Ecology and Evolution},
  year={2018},
  publisher={Wiley Online Library}
}
```

### Are there any related packages?

 * [lumier](https://github.com/richelbilderbeek/lumier): Shiny app to help create the function call needed
 * [BEASTmasteR](https://github.com/nmatzke/BEASTmasteR): tip-dating analyses using fossils as dated terminal taxa
 * [RBeast](https://github.com/beast-dev/RBeast): misc other things


### What are the FASTA files?

FASTA files `anthus_aco.fas` and `anthus_nd2.fas` from:
 
 * Van Els, Paul, and Heraldo V. Norambuena. "A revision of species limits in Neotropical pipits Anthus based on multilocus genetic and vocal data." Ibis.

Thanks to Paul van Els.

## `babette` misc

### Why the name?

 * `ba`: BeAutier
 * `bet`: BEasTier
 * `te`: TracErer

Later on, `mauricer` got added.

### Why the logo?

Initially, the logo was a low-tech remake of Babette, a maid in Beauty and the Beast. 
To prevent problems with Disney, a different logo was picked.

The current logo shows a swan, an animal considered to be graceful.
The swan is drawn by Jose Scholte, who kindly allowed her work to
be used for free, by attribution.

## Error: `libjvm.so: cannot open shared object file: No such file or directory`

For me, [this Stack Overflow post](https://stackoverflow.com/a/25932828) helped me out:

```
sudo mousepad /etc/ld.so.conf.d/java.conf
```

In that file put:

```
/usr/lib/jvm/java-8-oracle/jre/lib/amd64
/usr/lib/jvm/java-8-oracle/jre/lib/amd64/server
```

Save, close, restart R studio, fixed!

Else [this Stack Overflow post may be helpful](https://stackoverflow.com/a/43466434):

```
sudo apt-get install jdk-*
sudo R CMD javareconf
sudo R CMD javareconf -e
export LD_LIBRARY_PATH=$JAVA_LD_LIBRARY_PATH
sudo apt-get install r-cran-rjava
```
# Contributor Covenant Code of Conduct

## Our Pledge

In the interest of fostering an open and welcoming environment, we as
contributors and maintainers pledge to making participation in our project and
our community a harassment-free experience for everyone, regardless of age, body
size, disability, ethnicity, gender identity and expression, level of experience,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

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
obligated to maintain confidentiality regarding the reporter of an incident.
Further details of specific enforcement policies may be posted separately.

Project maintainers who do not follow or enforce the Code of Conduct in good
faith may face temporary or permanent repercussions as determined by other
members of the project's leadership.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage], version 1.4,
available at [http://contributor-covenant.org/version/1/4][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/4/
# Examples

See [the babette examples](https://github.com/richelbilderbeek/babette_examples).
---
title: "babette demo"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{babette demo}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r check_empty_cache_at_start, include = FALSE}
# This is only needed to pass the CRAN Windows build.
#
# This vignette is tested to clean up nicely on GitHub Actions
# and r-hub on Windows
#
unlink(
  dirname(beastier::get_beastier_tempfilename()),
  recursive = TRUE
)

beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
# beastierinstall::clear_beautier_cache() ; beastierinstall::clear_beastier_cache() # nolint
```


## Introduction

![](babette_logo.png)

This vignette briefly demonstrates multiple features of `babette`,
without going into much detail.

First, load the library:

```{r load_babette, results='hide', warning=FALSE, error=FALSE, message=FALSE}
library(babette)
```

This vignette shows how to:

 * Let `babette` run 'BEAST2'
 * Plot the posterior estimates
 * Show the effective sample sizes (ESS)
 * Show the summary statistics
 * Plot the posterior phylogenies

In all cases, this is done for a short MCMC chain length of 10K:

```{r}
inference_model <- create_test_inference_model()
```
Also, in all cases, we use the same BEAST2 options:

```{r}
beast2_options <- create_beast2_options(verbose = TRUE)
```


## Let `babette` run 'BEAST2'

For an alignment, we'll use a `babette` example alignment.

```{r}
if (is_beast2_installed()) {
  out <- bbt_run_from_model(
    fasta_filename = get_babette_path("anthus_aco_sub.fas"),
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

## Plot the posterior estimates

```{r}
if (is_beast2_installed()) {
  library(ggplot2)
  p <- ggplot(
    data = out$estimates,
    aes(x = Sample)
  )
  p + geom_line(aes(y = TreeHeight), color = "green")
  p + geom_line(aes(y = YuleModel), color = "red")
  p + geom_line(aes(y = birthRate), color = "blue")
}
```

## Show the effective sample sizes (ESS)

Effective sample sizes, with 20% burn-in removed:

```{r}
if (is_beast2_installed()) {
  traces <- remove_burn_ins(
    traces = out$estimates,
    burn_in_fraction = 0.2
  )
  esses <- t(
    calc_esses(
      traces,
      sample_interval = inference_model$mcmc$tracelog$log_every
    )
  )
  colnames(esses) <- "ESS"
  knitr::kable(esses)
}
```

For a reliable inference, use an ESS of at least 200.

## Show the summary statistics

```{r}
if (is_beast2_installed()) {
  sum_stats <- t(
    calc_summary_stats(
      traces$posterior,
      sample_interval = inference_model$mcmc$tracelog$log_every
    )
  )
  colnames(sum_stats) <- "Statistic"
  knitr::kable(sum_stats)
}
```

## Plot the posterior phylogenies

```{r fig.width=7, fig.height=7}
if (is_beast2_installed()) {
  plot_densitree(out$anthus_aco_sub_trees, width = 2)
}
```

```{r check_empty_cache_at_end, include = FALSE}
# This is only needed to pass the CRAN Windows build.
#
# This vignette is tested to clean up nicely on GitHub Actions
# and r-hub on Windows
#
unlink(
  dirname(beastier::get_beastier_tempfilename()),
  recursive = TRUE
)
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
# beastierinstall::clear_beautier_cache() ; beastierinstall::clear_beastier_cache() # nolint
```
---
title: "babette Tutorial"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{babette Tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r check_empty_cache_at_start, include = FALSE}
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
```

## Introduction

![](babette_logo.png)

This vignette is a tutorial how to use `babette` and its most
important `bbt_run_from_model` function.

First, load `babette`:

```{r load_babette}
library(babette)
```

The main function of `babette` is `bbt_run_from_model`. Here is part of its help:

```
Do a full run: create a 'BEAST2' configuration file (like BEAUti 2), 
run 'BEAST2', parse results (like Tracer)

Usage

bbt_run_from_model(
  fasta_filename,
  inference_model,
  beast2_options
)
```

Simplifying this to all arguments that do not have a default:

```
bbt_run_from_model(
  fasta_filename
)
```

## `fasta_filename`

`fasta_filename` is the argument to specify which
FASTA file to work on. `babette` is bundled with some FASTA files,
so obtaining a path to a FASTA file is easy:

```{r}
fasta_filename <- get_babette_path("anthus_aco_sub.fas")
library(testthat)
expect_true(file.exists(fasta_filename))
```

With `fasta_filename` available, we have the minimal
requirements to call `bbt_run_from_model` like this:

```
out <- bbt_run_from_model(fasta_filename)
```

Note that this code is not ran, as it would take too long.
The reason this would take too long, is that
the MCMC run that will be executed is set to one million states by default.
To specify the MCMC options and shorten this run,
the `mcmc` argument is used.

## `inference_model` and `mcmc`

The inference run's MCMC is part of the inference model.
To get an inference model with a short MCMC, create
a test inference model like this:

```{r}
inference_model <- create_test_inference_model()
names(inference_model)
```

`mcmc` is the `inference_model` argument to specify the MCMC run options:

```{r}
print(inference_model$mcmc$chain_length)
```

With these MCMC options, we can now call `bbt_run_from_model` in way that
it will finish fast:

```{r cache=TRUE}
if (is_beast2_installed()) {
  beast2_options <- create_beast2_options()
  out <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

The return value, `out` contains the results of the MCMC run.
For this tutorial, visualizing `out` is ignored, as the 'Demo' vignette
discusses this. 
Instead, we will work through the other `bbt_run_from_model` parameters.

## `site_model`

`site_model` is the `inference_model` parameter for a site model.
As this tutorial works on a DNA alignment, such a site model can also
be called a nucleotide substitution model.

Picking a site model is easy: just type:

```
create_site_model_
```

This will trigger auto-complete to show all site models.

The simplest site model is the Jukes-Cantor DNA substitution model.
To use this model in `babette`, do:

```{r}
inference_model <- create_test_inference_model(
  site_model = create_jc69_site_model()
)
```

Using this site model:

```{r cache=TRUE}
if (is_beast2_installed()) {
  beast2_options <- create_beast2_options()
  out <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

## `clock_model`

`clock_models` is the `inference_model` parameter for a clock model.

Picking a clock model is easy: just type:

```
create_clock_model_
```

This will trigger auto-complete to show all clock models.

The simplest site model is the strict clock model.
To use this model in `babette`, do:

```{r}
inference_model <- create_test_inference_model(
  clock_model = create_strict_clock_model()
)
```

Using this clock model:

```{r cache=TRUE}
if (is_beast2_installed()) {
  beast2_options <- create_beast2_options()
  out <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

## `tree_prior`

`tree_prior` is the `inference_model` parameter to select a tree prior.

Picking a tree prior is easy: just type:

```
create_tree_prior_
```

This will trigger auto-complete to show all tree priors.

The simplest tree prior is the Yule (pure-birth) tree prior.
To use this model in `babette`, do:

```{r}
inference_model <- create_test_inference_model(
  tree_prior = create_yule_tree_prior()
)
```

Using this tree prior:

```{r cache=TRUE}
if (is_beast2_installed()) {
  beast2_options <- create_beast2_options()
  out <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

## `mrca_prior`

`mrca_priors` is the `inference_model` parameter to use 
a Most Recent Common Ancestor (hence, MRCA) prior.
With such a prior, it can be specified which taxa have a shared
common ancestor and when it existed.

Here is how to specify that the first two taxa in a FASTA file
are sister species:

```{r}
mrca_prior <- create_mrca_prior(
  alignment_id = get_alignment_id(fasta_filename = fasta_filename),
  taxa_names = get_taxa_names(filename = fasta_filename)[1:2],
  is_monophyletic = TRUE
)
```

To specify when the MRCA of all taxa was present, we'll first
create a prior distribution of the crown age, after which we can
use that distribution.

To assume the crown age to follow a normal distribution,
with a mean of 15.0 (time units), with a standard deviation of 1.0,
use `create_normal_distr`:

```{r}
mrca_distr <- create_normal_distr(
  mean = 15.0,
  sigma = 1.0
)
```

To use that distribution in our MRCA prior:

```{r}
mrca_prior <- create_mrca_prior(
  alignment_id = get_alignment_id(fasta_filename = fasta_filename),
  taxa_names = get_taxa_names(filename = fasta_filename),
  mrca_distr = mrca_distr
)
```

Using such an MRCA prior:

```{r}
inference_model <- create_test_inference_model(
  mrca_prior = mrca_prior
)
```

```{r cache=TRUE}
if (is_beast2_installed()) {
  beast2_options <- create_beast2_options()
  out <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

```{r check_empty_cache_at_end, include = FALSE}
unlink(
  dirname(beastier::get_beastier_tempfilename()),
  recursive = TRUE
)
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
# beastierinstall::clear_beautier_cache() ; beastierinstall::clear_beastier_cache() # nolint
```
---
title: "Examples"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Examples}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r check_empty_cache_at_start, include = FALSE}
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
```


## Introduction

![](babette_logo.png)

This vignette shows some examples how to set up different inference models
in `babette`.

For all examples, do load `babette`:

```{r load_babette, results='hide', warning=FALSE, error=FALSE, message=FALSE}
library(babette)
```

All these examples check that `BEAST2` is installed at the default
location at `r get_default_beast2_path()`. If this is not the case,
we will use some fabricated output:

```{r}
posterior <- create_test_bbt_run_output()
posterior$anthus_aco_sub_trees <- posterior$anthus_aco_trees
names(posterior)
```

All examples read the alignment from a FASTA file (usually `my_fasta.fas`).

```{r}
fasta_filename <- get_babette_path("anthus_aco_sub.fas")
```

Instead of a full run, the MCMC chain length is shortened to 10K states,
with a measurement every 1K states:

```{r}
mcmc <- create_test_mcmc(chain_length = 10000)
```

We will re-create this MCMC setup, 
as doing so initializes it with new filenames for
temporary files. 
These temporary files should not exist before a run and should exist
after a run. Sure, there is the option to overwrite...

## Example #1: all default

Using all default settings, only specify a DNA alignment.

![Example #1: all default](all_default.png)

```{r example_1, cache=TRUE}
if (is_beast2_installed()) {
  inference_model <- create_inference_model(
    mcmc = mcmc
  )
  beast2_options <- create_beast2_options()
  posterior <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

All other parameters are set to their defaults, as in BEAUti.

```{r fig.width=7, fig.height=7}
plot_densitree(posterior$anthus_aco_sub_trees, width = 2)
```

## Example #2: using an MRCA prior to specify a crown age

![Example #2: using an MRCA prior to specify a crown age](mrca_crown_age.png)

An alternative is to date the node of the most recent common ancestor
of all taxa.

Create the MCMC:

```{r}
mcmc <- create_test_mcmc(chain_length = 10000)
```



```{r example_2_mrca, cache=FALSE}
if (is_beast2_installed()) {
  inference_model <- create_inference_model(
    mcmc = mcmc,
    mrca_prior = create_mrca_prior(
      taxa_names = sample(get_taxa_names(fasta_filename), size = 3),
      alignment_id = get_alignment_id(fasta_filename),
      is_monophyletic = TRUE,
      mrca_distr = create_normal_distr(
        mean = 15.0,
        sigma = 0.025
      )
    )
  )
  beast2_options <- create_beast2_options()
  posterior <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

Here we use an MRCA prior with fixed (non-estimated) values of the mean
and standard deviation for the common ancestor node's time.

```{r fig.width=7, fig.height=7}
plot_densitree(posterior$anthus_aco_sub_trees, width = 2)
```

## Example #3: JC69 site model

![Example #3: JC69 site model](jc69_2_4.png)



```{r example_3, cache=TRUE}
if (is_beast2_installed()) {
  inference_model <- create_inference_model(
    site_model = create_jc69_site_model(),
    mcmc = mcmc
  )
  beast2_options <- create_beast2_options()
  posterior <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

```{r fig.width=7, fig.height=7}
plot_densitree(posterior$anthus_aco_sub_trees, width = 2)
```

## Example #4: Relaxed clock log normal

![Example #4: Relaxed clock log normal](rln_2_4.png)



```{r example_4, cache=TRUE}
if (is_beast2_installed()) {
  inference_model <- create_inference_model(
    clock_model = create_rln_clock_model(),
    mcmc = mcmc
  )
  beast2_options <- create_beast2_options()
  posterior <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

```{r fig.width=7, fig.height=7}
plot_densitree(posterior$anthus_aco_sub_trees, width = 2)
```

## Example #5: Birth-Death tree prior

![Example #5: Birth-Death tree prior](bd_2_4.png)



```{r example_5, cache=TRUE}
if (is_beast2_installed()) {
  inference_model <- create_inference_model(
    tree_prior = create_bd_tree_prior(),
    mcmc = mcmc
  )
  beast2_options <- create_beast2_options()
  posterior <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

```{r fig.width=7, fig.height=7}
plot_densitree(posterior$anthus_aco_sub_trees, width = 2)
```

## Example #6: Yule tree prior with a normally distributed birth rate

![Example #6: Yule tree prior with a normally distributed birth rate](birth_rate_normal_2_4.png)

```{r example_6, cache=TRUE}
if (is_beast2_installed()) {
  inference_model <- create_inference_model(
    tree_prior = create_yule_tree_prior(
      birth_rate_distr = create_normal_distr(
        mean = 1.0,
        sigma = 0.1
      )
    ),
    mcmc = mcmc
  )
  beast2_options <- create_beast2_options()
  posterior <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

```{r fig.width=7, fig.height=7}
plot_densitree(posterior$anthus_aco_sub_trees, width = 2)
```

Thanks to Yacine Ben Chehida for this use case

## Example #7: HKY site model with a non-zero proportion of invariants

![Example #7: HKY site model with a non-zero proportion of invariants](hky_prop_invariant_0_5_2_4.png)

```{r example_7, cache=TRUE}
if (is_beast2_installed()) {
  inference_model <- create_inference_model(
    site_model = create_hky_site_model(
      gamma_site_model = create_gamma_site_model(prop_invariant = 0.5)
    ),
    mcmc = mcmc
  )
  beast2_options <- create_beast2_options()
  posterior <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

```{r fig.width=7, fig.height=7}
plot_densitree(posterior$anthus_aco_sub_trees, width = 2)
```

Thanks to Yacine Ben Chehida for this use case

## Example #8: Strict clock with a known clock rate

![Example #8: Strict clock with a known clock rate](strict_clock_rate_0_5_2_4.png)



```{r example_8, cache=TRUE}
if (is_beast2_installed()) {
  inference_model <- create_inference_model(
    clock_model = create_strict_clock_model(
      clock_rate_param = 0.5
    ),
    mcmc = mcmc
  )
  beast2_options <- create_beast2_options()
  posterior <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

```{r fig.width=7, fig.height=7}
plot_densitree(posterior$anthus_aco_sub_trees, width = 2)
```

Thanks to Paul van Els and Yacine Ben Chehida for this use case.

```{r check_empty_cache_at_end, include = FALSE}
unlink(
  dirname(beastier::get_beastier_tempfilename()),
  recursive = TRUE
)
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
# beastierinstall::clear_beautier_cache() ; beastierinstall::clear_beastier_cache() # nolint
```
---
title: "babette: Step by Step"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{babette: Step by Step}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r check_empty_cache_at_start, include = FALSE}
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
# beastierinstall::clear_beautier_cache() ; beastierinstall::clear_beastier_cache() # nolint
```


## Introduction

![](babette_logo.png)

This step-by-step demo shows how to run the `babette` pipeline in detail.

First, load `babette`:

```{r load_babette, results='hide', warning=FALSE, error=FALSE, message=FALSE}
library(babette)
```

In all cases, this is done for a short MCMC chain length of 10K:

```{r}
inference_model <- create_inference_model()
inference_model$mcmc$chain_length <- 10000
inference_model$mcmc$tracelog$filename <- normalizePath(
  get_beautier_tempfilename(
    pattern = "tracelog_", fileext = ".log"
  ),
  mustWork = FALSE
)
inference_model$mcmc$treelog$filename <- normalizePath(
  get_beautier_tempfilename(
    pattern = "treelog_", fileext = ".trees"  
  ),
  mustWork = FALSE
)
```

## Create a 'BEAST2' input file

This step is commonly done using BEAUti.
With `babette`, this can be done as follows:

```{r}
beast2_input_file <- tempfile(pattern = "beast2_", fileext = ".xml")
create_beast2_input_file_from_model(
  input_filename = get_babette_path("anthus_aco.fas"),
  inference_model = inference_model,
  output_filename = beast2_input_file
)
```

## Display (part of) the 'BEAST2' input file

```{r}
print(head(readLines(beast2_input_file)))
```

This file can both be loaded by BEAUti and be used by 'BEAST2'.

The file can be checked if it is indeed a valid input file:

```{r}
if (is_beast2_installed()) {
  is_beast2_input_file(beast2_input_file)
}
```

## Run MCMC

This step is commonly done using 'BEAST2' from the command-line or using its GUI.
With `babette`, this can be done as follows:

```{r}
if (is_beast2_installed()) {
  beast2_options <- create_beast2_options(
    input_filename = beast2_input_file
  )
  beastier::check_can_create_file(beast2_options$output_state_filename)
  beastier::check_can_create_treelog_file(beast2_options)
  run_beast2_from_options(
    beast2_options = beast2_options
  )
  testthat::expect_true(file.exists(beast2_options$output_state_filename))
}
```

## Display (part of) the 'BEAST2' output files

The `.log` file contains the model parameters and parameter estimates:

```{r}
if (is_beast2_installed()) {
  print(head(readLines(inference_model$mcmc$tracelog$filename)))
  print(tail(readLines(inference_model$mcmc$tracelog$filename)))
}
```

The `.trees` file contains the alignment, taxa and posterior trees:

```{r}
if (is_beast2_installed()) {
  print(head(readLines(inference_model$mcmc$treelog$filename)))
  print(tail(readLines(inference_model$mcmc$treelog$filename)))
}
```

The `.xml.state` file contains the final state of the MCMC run and the
MCMC operator acceptances thus far:

```{r}
if (is_beast2_installed()) {
  print(head(readLines(beast2_options$output_state_filename)))
  print(tail(readLines(beast2_options$output_state_filename)))
}
```

## Parse output

This step is commonly done using Tracer.
With `babette`, this can be done as follows.

Parsing `.log` file to obtain the parameter estimates:

```{r}
if (is_beast2_installed()) {
  knitr::kable(head(parse_beast_tracelog_file(inference_model$mcmc$tracelog$filename)))
}
```

Parsing `.trees` file to obtain the posterior phylogenies:

```{r fig.width = 7, fig.height = 7}
if (is_beast2_installed()) {
  plot_densitree(parse_beast_trees(inference_model$mcmc$treelog$filename))
}
```

Parsing `.xml.state` file to obtain the MCMC operator acceptances:

```{r}
if (is_beast2_installed()) {
  knitr::kable(head(parse_beast_state_operators(beast2_options$output_state_filename)))
}
```

```{r cleunup, include = FALSE}
if (is_beast2_installed()) {
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```


```{r check_empty_cache_at_end, include = FALSE}
unlink(
  dirname(beastier::get_beastier_tempfilename()),
  recursive = TRUE
)
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
# beastierinstall::clear_beautier_cache() ; beastierinstall::clear_beastier_cache() # nolint
```
---
title: "Nested Sampling"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Nested Sampling}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r check_empty_cache_at_start, include = FALSE}
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
```

## Introduction

![](babette_logo.png)

This vignette demonstrates how to use the Nested Sampling approach
to obtain the marginal likelihood and a Bayes factor,
as described in [1]

## Setup

```{r}
library(babette)
library(testthat)
```

`babette` needs to have 'BEAST2' installed to work.

In the case 'BEAST2' is not installed, we'll use
this fabricated data:

```{r}
out_jc69 <- create_test_bbt_run_output()
out_jc69$ns$marg_log_lik <- c(-1.1)
out_jc69$ns$marg_log_lik_sd <- c(0.1)
out_gtr <- out_jc69
```

Here we setup how to interpret the Bayes factor:

```{r}
interpret_bayes_factor <- function(bayes_factor) {
  if (bayes_factor < 10^-2.0) {
    "decisive for GTR"
  } else if (bayes_factor < 10^-1.5) {
    "very strong for GTR"
  } else if (bayes_factor < 10^-1.0) {
    "strong for GTR"
  } else if (bayes_factor < 10^-0.5) {
    "substantial for GTR"
  } else if (bayes_factor < 10^0.0) {
    "barely worth mentioning for GTR"
  } else if (bayes_factor < 10^0.5) {
    "barely worth mentioning for JC69"
  } else if (bayes_factor < 10^1.0) {
    "substantial for JC69"
  } else if (bayes_factor < 10^1.5) {
    "strong for JC69"
  } else if (bayes_factor < 10^2.0) {
    "very strong for JC69"
  } else {
    "decisive for JC69"
  }
}
expect_equal(interpret_bayes_factor(1 / 123.0), "decisive for GTR")
expect_equal(interpret_bayes_factor(1 / 85.0), "very strong for GTR")
expect_equal(interpret_bayes_factor(1 / 12.5), "strong for GTR")
expect_equal(interpret_bayes_factor(1 / 8.5), "substantial for GTR")
expect_equal(interpret_bayes_factor(1 / 1.5), "barely worth mentioning for GTR")
expect_equal(interpret_bayes_factor(0.99), "barely worth mentioning for GTR")
expect_equal(interpret_bayes_factor(1.01), "barely worth mentioning for JC69")
expect_equal(interpret_bayes_factor(1.5), "barely worth mentioning for JC69")
expect_equal(interpret_bayes_factor(8.5), "substantial for JC69")
expect_equal(interpret_bayes_factor(12.5), "strong for JC69")
expect_equal(interpret_bayes_factor(85.0), "very strong for JC69")
expect_equal(interpret_bayes_factor(123.0), "decisive for JC69")
```

## Experiment

In this experiment, we will use the same DNA alignment to see which
DNA nucleotide substitution model is the better fit.

Load the DNA alignment, a subset of taxa from [2]:


```{r fig.width=7}
fasta_filename <- get_babette_path("anthus_aco_sub.fas")
image(ape::read.FASTA(fasta_filename))
```

## Do the run

In this vignette, the MCMC run is set up to be short:

```{r}
mcmc <- beautier::create_test_ns_mcmc()
```

  For academic research, better use a longer MCMC chain (with an effective
sample size above 200).

Here we do two 'BEAST2' runs, with both site models:

```{r}
if (is_beast2_installed() && is_beast2_pkg_installed("NS")) {
  inference_model <- create_inference_model(
    site_model = beautier::create_jc69_site_model(),
    mcmc = mcmc
  )  
  beast2_options <- create_beast2_options(
    beast2_path = beastier::get_default_beast2_bin_path()
  )
  out_jc69 <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  inference_model <- create_inference_model(
    site_model = beautier::create_gtr_site_model(),
    mcmc = mcmc
  )
  beast2_options <- create_beast2_options(
    beast2_path = beastier::get_default_beast2_bin_path()
  )
  out_gtr <- bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
```

I will display the marginal likelihoods in a nice table:

```{r}
if (is_beast2_installed() && is_beast2_pkg_installed("NS")) {
  df <- data.frame(
    model = c("JC69", "GTR"),
    mar_log_lik = c(out_jc69$ns$marg_log_lik, out_gtr$ns$marg_log_lik),
    mar_log_lik_sd = c(out_jc69$ns$marg_log_lik_sd, out_gtr$ns$marg_log_lik_sd)
  )
  knitr::kable(df)
}
```

The Bayes factor is ratio between the marginal (non-log) likelihoods.
In this case, we use the simpler JC69 model as a focus:

```{r}
if (is_beast2_installed() && is_beast2_pkg_installed("NS")) {
  bayes_factor <- exp(out_jc69$ns$marg_log_lik) / exp(out_gtr$ns$marg_log_lik)
  print(interpret_bayes_factor(bayes_factor))
}
```

Whatever the support is, be sure to take the error of the marginal
likelihood estimation into account.

## References

 * [1] Maturana, P., Brewer, B. J., Klaere, S., & Bouckaert, R. (2017).
   Model selection and parameter inference in phylogenetics
   using Nested Sampling. arXiv preprint arXiv:1703.05471.
 * [2] Van Els, Paul, and Heraldo V. Norambuena. "A revision of species
   limits in Neotropical pipits Anthus based on multilocus genetic and
   vocal data." Ibis.

```{r check_empty_cache_at_end, include = FALSE}
unlink(
  dirname(beastier::get_beastier_tempfilename()),
  recursive = TRUE
)
beautier::check_empty_beautier_folder()
beastier::check_empty_beastier_folder()
# beastierinstall::clear_beautier_cache() ; beastierinstall::clear_beastier_cache() # nolint
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beast2_pkgs.R
\name{check_beast2_pkgs}
\alias{check_beast2_pkgs}
\title{Checks if \code{\link{bbt_run}} has the 'BEAST2' packages needed to process
its arguments. Will \link{stop} if not.}
\usage{
check_beast2_pkgs(mcmc, beast2_path = get_default_beast2_bin_path())
}
\arguments{
\item{mcmc}{the MCMC options,
see \link[beautier]{create_mcmc}}

\item{beast2_path}{name of either a 'BEAST2'  binary file
(usually simply \code{beast})
or a 'BEAST2'  jar file
(usually has a \code{.jar} extension).
Use \code{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \code{get_default_beast2_jar_path} to get
the default BEAST jar file's path}
}
\description{
For example, to use a Nested Sampling MCMC, the 'BEAST2' 'NS' package
needs to be installed.
}
\examples{
if (is_beast2_installed()) {
  # Minimal BEAST2 setup
  check_beast2_pkgs(mcmc = create_mcmc())

  # BEAST2 with NS package installed
  if (is_beast2_ns_pkg_installed()) {
    check_beast2_pkgs(mcmc = create_ns_mcmc())
  }
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_babette_path.R
\name{get_babette_path}
\alias{get_babette_path}
\title{Get the full path of a file in the \code{inst/extdata} folder}
\usage{
get_babette_path(filename)
}
\arguments{
\item{filename}{the file's name, without the path}
}
\value{
the full path of the filename, if and only if
  the file is present. Will stop otherwise.
}
\description{
Get the full path of a file in the \code{inst/extdata} folder
}
\examples{
get_babette_path("anthus_aco.fas")
}
\seealso{
for more files, use \code{\link{get_babette_paths}}
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
default_params_doc(
  beast2_input_filename,
  beast2_options,
  beast2_output_log_filename,
  beast2_output_state_filename,
  beast2_output_trees_filenames,
  beast2_path,
  beast2_working_dir,
  cleanup,
  clock_model,
  clock_models,
  fasta_filename,
  fasta_filenames,
  inference_model,
  mcmc,
  mrca_prior,
  mrca_priors,
  overwrite,
  rng_seed,
  site_model,
  site_models,
  tipdates_filename,
  tree_prior,
  tree_priors,
  verbose
)
}
\arguments{
\item{beast2_input_filename}{path of the 'BEAST2'  configuration file.
By default, this file is put in a temporary folder with a random
filename, as the user needs not read it: it is used as input of 'BEAST2'.
Specifying a \code{beast2_input_filename} allows to store that file
in a more permanently stored location.}

\item{beast2_options}{'BEAST2'  options,
as can be created by \link[beastier]{create_beast2_options}}

\item{beast2_output_log_filename}{name of the log file created by 'BEAST2',
containing the parameter estimates in time. By default, this
file is put a temporary folder with a random
filename, as the user needs not read it: its content
is parsed and returned by this function.
Specifying a \code{beast2_output_log_filename} allows to store that file
in a more permanently stored location.}

\item{beast2_output_state_filename}{name of the final state file created
by 'BEAST2', containing the operator acceptances. By default, this
file is put a temporary folder with a random
filename, as the user needs not read it: its content
is parsed and returned by this function.
Specifying a \code{beast2_output_state_filename} allows to store that file
in a more permanently stored location.}

\item{beast2_output_trees_filenames}{name of the one or more trees
files created by 'BEAST2', one per alignment. By default, these
files are put a temporary folder with a random
filename, as the user needs not read it: their content
is parsed and returned by this function.
Specifying \code{beast2_output_trees_filenames} allows to store these
one or more files in a more permanently stored location.}

\item{beast2_path}{name of either a 'BEAST2'  binary file
(usually simply \code{beast})
or a 'BEAST2'  jar file
(usually has a \code{.jar} extension).
Use \code{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \code{get_default_beast2_jar_path} to get
the default BEAST jar file's path}

\item{beast2_working_dir}{the folder 'BEAST2'  will work in. This is
an (empty) temporary folder by default. This allows to call
'BEAST2'  in multiple parallel processes, as each process can have
its own working directory}

\item{cleanup}{set to FALSE to keep all temporary files}

\item{clock_model}{one clock model,
see \link[beautier]{create_clock_model}}

\item{clock_models}{one or more clock models,
see \link[beautier]{create_clock_models}}

\item{fasta_filename}{a FASTA filename}

\item{fasta_filenames}{one or more FASTA filename, each with one alignment}

\item{inference_model}{a Bayesian phylogenetic inference model,
as returned by \link[beautier]{create_inference_model}}

\item{mcmc}{the MCMC options,
see \link[beautier]{create_mcmc}}

\item{mrca_prior}{one Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}

\item{overwrite}{will 'BEAST2'  overwrite files? Like 'BEAST2',
this is set to \link{TRUE} by default.
If \link{TRUE}, 'BEAST2'  will overwrite the
\code{beast2_options$output_state_filename} if its present.
If \link{FALSE}, 'BEAST2'  will not overwrite the
\code{beast2_options$output_state_filename} if its present
and \link{babette} will give an error message.
Note that if \code{overwrite} is set to \link{FALSE} when
a \code{tracelog} (see \link{create_tracelog}),
\code{screenlog} (see \link{create_screenlog})
or \code{treelog} (see \link{create_treelog})
file already exists,
'BEAST2'  (and thus \link{babette}) will freeze.}

\item{rng_seed}{the random number generator seed. Must be either
\code{NA} or a positive non-zero value. An RNG seed of \code{NA}
results in 'BEAST2'  picking a random seed.}

\item{site_model}{one site model,
see \link[beautier]{create_site_models}}

\item{site_models}{one or more site models,
see \link[beautier]{create_site_models}}

\item{tipdates_filename}{name of the file containing tip dates}

\item{tree_prior}{one tree priors,
as created by \link[beautier]{create_tree_prior}}

\item{tree_priors}{one or more tree priors,
see \link[beautier]{create_tree_priors}}

\item{verbose}{set to TRUE for more output}
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
% Please edit documentation in R/create_test_bbt_run_output.R
\name{create_test_bbt_run_output}
\alias{create_test_bbt_run_output}
\title{Get an example output of \code{\link{bbt_run}}
or \code{\link{bbt_run_from_model}}.}
\usage{
create_test_bbt_run_output()
}
\value{
the same results as \code{\link{bbt_run}}
  or \code{\link{bbt_run_from_model}}
}
\description{
This output is used in testing.
}
\examples{
create_test_bbt_run_output()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbt_run.R
\name{bbt_run}
\alias{bbt_run}
\title{Run BEAST2}
\usage{
bbt_run(
  fasta_filename,
  tipdates_filename = NA,
  site_model = beautier::create_jc69_site_model(),
  clock_model = beautier::create_strict_clock_model(),
  tree_prior = beautier::create_yule_tree_prior(),
  mrca_prior = NA,
  mcmc = beautier::create_mcmc(),
  beast2_input_filename = beastier::create_temp_input_filename(),
  rng_seed = 1,
  beast2_output_state_filename = beastier::create_temp_state_filename(),
  beast2_path = beastier::get_default_beast2_path(),
  overwrite = TRUE,
  verbose = FALSE
)
}
\arguments{
\item{fasta_filename}{a FASTA filename}

\item{tipdates_filename}{name of the file containing tip dates}

\item{site_model}{one site model,
see \link[beautier]{create_site_models}}

\item{clock_model}{one clock model,
see \link[beautier]{create_clock_model}}

\item{tree_prior}{one tree priors,
as created by \link[beautier]{create_tree_prior}}

\item{mrca_prior}{one Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{mcmc}{the MCMC options,
see \link[beautier]{create_mcmc}}

\item{beast2_input_filename}{path of the 'BEAST2'  configuration file.
By default, this file is put in a temporary folder with a random
filename, as the user needs not read it: it is used as input of 'BEAST2'.
Specifying a \code{beast2_input_filename} allows to store that file
in a more permanently stored location.}

\item{rng_seed}{the random number generator seed. Must be either
\code{NA} or a positive non-zero value. An RNG seed of \code{NA}
results in 'BEAST2'  picking a random seed.}

\item{beast2_output_state_filename}{name of the final state file created
by 'BEAST2', containing the operator acceptances. By default, this
file is put a temporary folder with a random
filename, as the user needs not read it: its content
is parsed and returned by this function.
Specifying a \code{beast2_output_state_filename} allows to store that file
in a more permanently stored location.}

\item{beast2_path}{name of either a 'BEAST2'  binary file
(usually simply \code{beast})
or a 'BEAST2'  jar file
(usually has a \code{.jar} extension).
Use \code{get_default_beast2_bin_path} to get
the default BEAST binary file's path
Use \code{get_default_beast2_jar_path} to get
the default BEAST jar file's path}

\item{overwrite}{will 'BEAST2'  overwrite files? Like 'BEAST2',
this is set to \link{TRUE} by default.
If \link{TRUE}, 'BEAST2'  will overwrite the
\code{beast2_options$output_state_filename} if its present.
If \link{FALSE}, 'BEAST2'  will not overwrite the
\code{beast2_options$output_state_filename} if its present
and \link{babette} will give an error message.
Note that if \code{overwrite} is set to \link{FALSE} when
a \code{tracelog} (see \link{create_tracelog}),
\code{screenlog} (see \link{create_screenlog})
or \code{treelog} (see \link{create_treelog})
file already exists,
'BEAST2'  (and thus \link{babette}) will freeze.}

\item{verbose}{set to TRUE for more output}
}
\value{
a list with the following elements:\cr
\itemize{
  \item{
    \code{estimates}: a data frame with 'BEAST2'
    parameter estimates
  }
  \item{
    \code{[alignment_id]_trees}: a \code{multiPhylo}
    containing the phylogenies
    in the 'BEAST2' posterior. \code{[alignment_id]}
    is the ID of the alignment. For example,
    when running \code{\link{bbt_run}} with
    \code{anthus_aco.fas}, this element will have
    name \code{anthus_aco_trees}
  }
  \item{
    \code{operators}: a data frame with the
    'BEAST2' MCMC operator acceptances
  }
  \item{
    \code{output}: a numeric vector with the output
    sent to standard output and error streams
  }
  \item{
    \code{ns}: (optional) the results of a marginal likelihood estimation,
    will exist only when \code{\link[beautier]{create_ns_mcmc}} was
    used for the MCMC.
    This structure will contain the following elements:
    \itemize{
      \item \code{marg_log_lik} the marginal log likelihood estimate
      \item \code{marg_log_lik_sd} the standard deviation around the estimate
      \item \code{estimates} the parameter estimates
        created during the marginal likelihood estimation
      \item \code{trees} the trees
        created during the marginal likelihood estimation
    }
  }
}
}
\description{
Do a full BEAST2 run: create a 'BEAST2' configuration file (like 'BEAUti 2'),
run 'BEAST2', parse results (like 'Tracer')
}
\details{
Prefer using \code{\link{bbt_run_from_model}}, as it has a cleaner interface.
}
\examples{
if (is_beast2_installed()) {

  # Setup for a short run
  mcmc <- create_test_mcmc()

  # Store filenames for cleanup.
  # Note that 'bbt_run_from_model allows for easier cleanup
  mcmc$tracelog$filename <- tempfile()
  mcmc$treelog$filename <- tempfile()
  mcmc$screenlog$filename <- tempfile()
  beast2_input_filename <- tempfile()
  beast2_output_state_filename <- tempfile()

  bbt_run(
    fasta_filename = get_babette_path("anthus_aco.fas"),
    beast2_input_filename = beast2_input_filename,
    beast2_output_state_filename = beast2_output_state_filename,
    mcmc = mcmc
  )

  # Cleanup
  # Again, note that 'bbt_run_from_model allows for easier cleanup
  file.remove(mcmc$tracelog$filename)
  file.remove(mcmc$treelog$filename)
  file.remove(mcmc$screenlog$filename)
  file.remove(beast2_input_filename)
  file.remove(beast2_output_state_filename)
}
}
\seealso{
Use \code{\link[tracerer]{remove_burn_ins}}
  to remove the burn-ins from
  the posterior's estimates (\code{posterior$estimates})
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbt_self_test.R
\name{bbt_self_test}
\alias{bbt_self_test}
\title{Do a self test to verify \link{babette} that works correctly.}
\usage{
bbt_self_test(beast2_options = beastier::create_beast2_options())
}
\arguments{
\item{beast2_options}{'BEAST2'  options,
as can be created by \link[beastier]{create_beast2_options}}
}
\description{
Do a self test to verify \link{babette} that works correctly.
}
\examples{
# Will stop if BEAST2 is not installed
if (is_beast2_installed()) {
  bbt_self_test()
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_beast2_output.R
\name{parse_beast2_output}
\alias{parse_beast2_output}
\title{Process the 'BEAST2' output dependent on 'BEAST2' package specifics}
\usage{
parse_beast2_output(out, inference_model)
}
\arguments{
\item{out}{a list with the complete babette output, with elements:
\itemize{
  \item \code{output} textual output of a 'BEAST2' run
}}

\item{inference_model}{a Bayesian phylogenetic inference model,
as returned by \link[beautier]{create_inference_model}}
}
\value{
complete babette output with added attributes,
  which depends on the 'BEAST2' package.
  \itemize{
    \item \code{marg_log_lik} the marginal log likelihood estimate
    \item \code{marg_log_lik_sd} the standard deviation around the estimate
    \item \code{estimates} the parameter estimates
      created during the marginal likelihood estimation
    \item \code{trees} the trees
      created during the marginal likelihood estimation
  }
}
\description{
Process the 'BEAST2' output dependent on 'BEAST2' package specifics
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbt_delete_temp_files.R
\name{bbt_delete_temp_files}
\alias{bbt_delete_temp_files}
\title{Delete all the temporary files created by \link{bbt_run_from_model}}
\usage{
bbt_delete_temp_files(inference_model, beast2_options)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model,
as returned by \link[beautier]{create_inference_model}}

\item{beast2_options}{'BEAST2'  options,
as can be created by \link[beastier]{create_beast2_options}}
}
\description{
Delete all the temporary files created by \link{bbt_run_from_model}
}
\examples{
if (is_beast2_installed()) {
  # Do a minimal run
  inference_model <- create_test_inference_model()
  beast2_options <- create_beast2_options()
  bbt_run_from_model(
    fasta_filename = get_fasta_filename(),
    inference_model = inference_model,
    beast2_options = beast2_options
  )

  # Cleanup
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_babette_paths.R
\name{get_babette_paths}
\alias{get_babette_paths}
\title{Get the full paths of files in the \code{inst/extdata} folder}
\usage{
get_babette_paths(filenames)
}
\arguments{
\item{filenames}{the files' names, without the path}
}
\value{
the filenames' full paths, if and only if
  all files are present. Will stop otherwise.
}
\description{
Get the full paths of files in the \code{inst/extdata} folder
}
\examples{
get_babette_paths(c("anthus_aco.fas", "anthus_nd2.fas"))
}
\seealso{
for one file, use \code{\link{get_babette_path}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_beast2_output_to_ns.R
\name{parse_beast2_output_to_ns}
\alias{parse_beast2_output_to_ns}
\title{Parse BEAST2 NS output}
\usage{
parse_beast2_output_to_ns(output)
}
\arguments{
\item{output}{screen output}
}
\value{
a list with the following elements:
  \itemize{
    \item \code{marg_log_lik} the marginal log likelihood estimate
    \item \code{marg_log_lik_sd} the standard deviation around the estimate
  }
}
\description{
Parse the BEAST2 output when run with the BEAST2 NS ('Nested Sampling')
package.
}
\examples{
parse_beast2_output_to_ns(
  output = create_test_ns_output()
)
}
\seealso{
use \code{\link{create_test_ns_output}} to obtain
a test screen output.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbt_continue.R
\name{bbt_continue}
\alias{bbt_continue}
\title{Continue a BEAST2 run}
\usage{
bbt_continue(fasta_filename, inference_model, beast2_options)
}
\arguments{
\item{fasta_filename}{a FASTA filename}

\item{inference_model}{a Bayesian phylogenetic inference model,
as returned by \link[beautier]{create_inference_model}}

\item{beast2_options}{'BEAST2'  options,
as can be created by \link[beastier]{create_beast2_options}}
}
\value{
a list with the following elements:\cr
\itemize{
  \item{
    \code{estimates}: a data frame with 'BEAST2'
    parameter estimates
  }
  \item{
    \code{[alignment_id]_trees}: a \code{multiPhylo}
    containing the phylogenies
    in the 'BEAST2' posterior. \code{[alignment_id]}
    is the ID of the alignment. For example,
    when running \code{bbt_run_from_model} with
    \code{anthus_aco.fas}, this element will have
    name \code{anthus_aco_trees}
  }
  \item{
    \code{operators}: a data frame with the
    'BEAST2' MCMC operator acceptances
  }
  \item{
    \code{output}: a numeric vector with the output
    sent to standard output and error streams
  }
  \item{
    \code{ns}: (optional) the results of a marginal likelihood estimation,
    will exist only when \code{create_ns_mcmc} was
    used for \code{mcmc}.
    This structure will contain the following elements:
    \itemize{
      \item \code{marg_log_lik} the marginal log likelihood estimate
      \item \code{marg_log_lik_sd} the standard deviation around the estimate
      \item \code{estimates} the parameter estimates
        created during the marginal likelihood estimation
      \item \code{trees} the trees
        created during the marginal likelihood estimation
    }
  }
}
}
\description{
Do a full run: create a 'BEAST2' configuration file (like 'BEAUti 2'),
run 'BEAST2', parse results (like 'Tracer')
}
\examples{
if (is_beast2_installed()) {

  # A simple FASTA file
  fasta_filename <- beautier::get_beautier_path("test_output_0.fas")

  # Simple short inference
  inference_model <- create_test_inference_model()

  # Default BEAST2 options
  beast2_options <- create_beast2_options()

  bbt_run_from_model(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )

  bbt_continue(
    fasta_filename = fasta_filename,
    inference_model = inference_model,
    beast2_options = beast2_options
  )

  # Cleanup
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
}
\seealso{
Use \code{\link[tracerer]{remove_burn_ins}}
  to remove the burn-ins from
  the posterior's estimates (\code{posterior$estimates})
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prepare_file_creation.R
\name{prepare_file_creation}
\alias{prepare_file_creation}
\title{Internal function to prepare for 'BEAST2' creating files}
\usage{
prepare_file_creation(inference_model, beast2_options)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model,
as returned by \link[beautier]{create_inference_model}}

\item{beast2_options}{'BEAST2'  options,
as can be created by \link[beastier]{create_beast2_options}}
}
\description{
The inference model and 'BEAST2' options contain paths that may point
to sub-sub-sub folders. Create those folders and test
if these folders can be written to
}
\examples{
# This example will fail on the CRAN
# r-oldrel-macos-x86_64 platform
if (rappdirs::app_dir()$os != "mac") {
  # For a test inference model, the files can be prepared
  inference_model <- create_test_inference_model()
  beast2_options <- create_beast2_options()
  prepare_file_creation(inference_model, beast2_options)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alignment_ids_from_xml.R
\name{get_alignment_ids_from_xml}
\alias{get_alignment_ids_from_xml}
\title{Get the alignment IDs from one or more 'BEAST2' XML input files.}
\usage{
get_alignment_ids_from_xml(xml_filename)
}
\arguments{
\item{xml_filename}{name of a 'BEAST2' XML input file.}
}
\value{
a character vector with one or more alignment IDs.
}
\description{
Get the alignment IDs from one or more 'BEAST2' XML input files.
}
\examples{
alignment_ids <- get_alignment_ids_from_xml(
  get_babette_path("anthus_2_4.xml")
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bbt_run_from_model.R
\name{bbt_run_from_model}
\alias{bbt_run_from_model}
\title{Run BEAST2}
\usage{
bbt_run_from_model(
  fasta_filename,
  inference_model = beautier::create_inference_model(),
  beast2_options = beastier::create_beast2_options()
)
}
\arguments{
\item{fasta_filename}{a FASTA filename}

\item{inference_model}{a Bayesian phylogenetic inference model,
as returned by \link[beautier]{create_inference_model}}

\item{beast2_options}{'BEAST2'  options,
as can be created by \link[beastier]{create_beast2_options}}
}
\value{
a list with the following elements:\cr
\itemize{
  \item{
    \code{estimates}: a data frame with 'BEAST2'
    parameter estimates
  }
  \item{
    \code{[alignment_id]_trees}: a \code{multiPhylo}
    containing the phylogenies
    in the 'BEAST2' posterior. \code{[alignment_id]}
    is the ID of the alignment. For example,
    when running \code{bbt_run_from_model} with
    \code{anthus_aco.fas}, this element will have
    name \code{anthus_aco_trees}
  }
  \item{
    \code{operators}: a data frame with the
    'BEAST2' MCMC operator acceptances
  }
  \item{
    \code{output}: a numeric vector with the output
    sent to standard output and error streams
  }
  \item{
    \code{ns}: (optional) the results of a marginal likelihood estimation,
    will exist only when \code{create_ns_mcmc} was
    used for \code{mcmc}.
    This structure will contain the following elements:
    \itemize{
      \item \code{marg_log_lik} the marginal log likelihood estimate
      \item \code{marg_log_lik_sd} the standard deviation around the estimate
      \item \code{estimates} the parameter estimates
        created during the marginal likelihood estimation
      \item \code{trees} the trees
        created during the marginal likelihood estimation
    }
  }
}
}
\description{
Do a full run: create a 'BEAST2' configuration file (like 'BEAUti 2'),
run 'BEAST2', parse results (like 'Tracer')
}
\examples{
if (is_beast2_installed()) {

  # Simple short inference
  inference_model <- create_test_inference_model()

  # Default BEAST2 options
  beast2_options <- create_beast2_options()

  bbt_run_from_model(
    fasta_filename = get_babette_path("anthus_aco.fas"),
    inference_model = inference_model,
    beast2_options = beast2_options
  )

  # Cleanup
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
}
\seealso{
Use \code{\link[tracerer]{remove_burn_ins}}
  to remove the burn-ins from
  the posterior's estimates (\code{posterior$estimates})
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/babette.R
\docType{package}
\name{babette}
\alias{babette}
\title{babette: A package for Bayesian phylogenetics.}
\description{
'babette' provides for an alternative workflow of using
the popular phylogenetics tool 'BEAST2', including
it peripheral tools. From an alignment and
inference model, a posterior of jointly estimated
phylogenies and parameter estimates is generated.
}
\examples{
if (is_beast2_installed()) {

  inference_model <- create_test_inference_model()
  beast2_options <- create_beast2_options()

  bbt_run_from_model(
    fasta_filename = get_babette_path("anthus_aco.fas"),
    inference_model = inference_model,
    beast2_options = beast2_options
  )

  # Clean up temporary files created by babette
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
}
\seealso{
Use \link{bbt_self_test} to do verify \link{babette} is installed
correctly.\cr

These are packages associated with 'babette':
\itemize{
  \item{
    '\link[beautier]{beautier}' creates 'BEAST2' input files.
  }
  \item{
    '\link[beastier]{beastier}' runs 'BEAST2'.
  }
  \item{
    '\link[mauricer]{mauricer}' does 'BEAST2' package management.
  }
  \item{
    '\link[tracerer]{tracerer}' parses 'BEAST2' output files.
  }
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_densitree.R
\name{plot_densitree}
\alias{plot_densitree}
\title{Draw multiple trees on top of one another.}
\usage{
plot_densitree(phylos, ...)
}
\arguments{
\item{phylos}{one or more phylogenies, must be of class \code{multiPhylo}}

\item{...}{options to be passed to \code{phangorn}'s
\link[phangorn]{densiTree} function}
}
\value{
nothing. Will produce a plot.
}
\description{
Draw multiple trees on top of one another.
}
\examples{
if (is_beast2_installed()) {
  inference_model <- create_test_inference_model()
  beast2_options <- create_beast2_options()

   out <- bbt_run_from_model(
    get_babette_path("anthus_aco.fas"),
    inference_model = inference_model,
    beast2_options = beast2_options
  )
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )

  plot_densitree(out$anthus_aco_trees)

  # Clean up temporary files created by babette
  bbt_delete_temp_files(
    inference_model = inference_model,
    beast2_options = beast2_options
  )
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_test_ns_output.R
\name{create_test_ns_output}
\alias{create_test_ns_output}
\title{Create NS testing output}
\usage{
create_test_ns_output()
}
\description{
Create testing output similar to when running a 'BEAST2' run
with nested sampling
}
\examples{
create_test_ns_output()
}
\seealso{
Use \link{parse_beast2_output_to_ns} to parse
this output to a Nested Sampling result.
See \link[beautier]{create_ns_mcmc} to see how to do a marginal
likelihood estimation using Nested Sampling.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/update_babette.R
\name{update_babette}
\alias{update_babette}
\title{Update all babette dependencies, by installing their
latest versions}
\usage{
update_babette(upgrade = "default")
}
\arguments{
\item{upgrade}{Should package dependencies be upgraded? One of "default", "ask", "always", or "never". "default"
respects the value of the \code{R_REMOTES_UPGRADE} environment variable if set,
and falls back to "ask" if unset. "ask" prompts the user for which out of
date packages to upgrade. For non-interactive sessions "ask" is equivalent
to "always". \code{TRUE} and \code{FALSE} are also accepted and correspond to
"always" and "never" respectively.}
}
\description{
Update all babette dependencies, by installing their
latest versions
}
\examples{
\dontrun{
  # Updates the babette dependencies without asking
}
}
\author{
Giovanni Laudanno, Richèl J.C. Bilderbeek
}
