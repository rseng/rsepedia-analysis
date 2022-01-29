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
