# beautier

[![Peer Review Status](https://badges.ropensci.org/209_status.svg)](https://github.com/ropensci/onboarding/issues/209)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/beautier)](https://cran.r-project.org/package=beautier)
[![](http://cranlogs.r-pkg.org/badges/grand-total/beautier)]( https://CRAN.R-project.org/package=beautier)
[![](http://cranlogs.r-pkg.org/badges/beautier)](https://CRAN.R-project.org/package=beautier)
[![DOI](https://zenodo.org/badge/53443354.svg)](https://zenodo.org/badge/latestdoi/53443354)

Branch   |[![GitHub Actions logo](man/figures/GitHubActions.png)](https://github.com/ropensci/beautier/actions)|[![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.com)                                                  |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
---------|-----------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------
`master` |![R-CMD-check](https://github.com/ropensci/beautier/workflows/R-CMD-check/badge.svg?branch=master)   |[![Build Status](https://travis-ci.com/ropensci/beautier.svg?branch=master)](https://travis-ci.com/ropensci/beautier) |[![codecov.io](https://codecov.io/github/ropensci/beautier/coverage.svg?branch=master)](https://codecov.io/github/ropensci/beautier/branch/master)
`develop`|![R-CMD-check](https://github.com/ropensci/beautier/workflows/R-CMD-check/badge.svg?branch=develop)  |[![Build Status](https://travis-ci.com/ropensci/beautier.svg?branch=develop)](https://travis-ci.com/ropensci/beautier)|[![codecov.io](https://codecov.io/github/ropensci/beautier/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/beautier/branch/develop)

`beautier` is `BEAUti` for R.

![beautier logo](man/figures/beautier_logo.png)

The purpose of `beautier` is to create 
[a valid BEAST2 XML input file](inst/extdata/2_4.xml)
from a n inference model. In this way, a scientific pipeline using 
`BEAST2` can be fully scripted, instead of using `BEAUti`'s GUI.

`beautier` is part of the [`babette`](https://github.com/ropensci/babette) package suite:

 * [`beautier`](https://github.com/ropensci/beautier) create a BEAST2 input (`.xml`) file from an inference model.
 * [`tiebeaur`](https://github.com/richelbilderbeek/tiebeaur) creates an inference model from a BEAST2 input (`.xml`) file :warning: experimental :warning:
 * [`beastier`](https://github.com/ropensci/beastier) runs BEAST2
 * [`tracerer`](https://github.com/ropensci/tracerer) pastes BEAST2 output (`.log`, `.trees`, etc) files.
 * [`mauricer`](https://github.com/ropensci/mauricer) install BEAST2 packages

Related R packages:

 * [`beautier_on_windows`](https://github.com/richelbilderbeek/beautier_on_windows): verifies
   `beautier` builds on Windows
 * [`lumier`](https://github.com/ropensci/lumier): Shiny app to help create the function call needed

## Examples

See [examples](doc/examples.md).

## Installation

`beautier` can be installed:

 * Latest CRAN version: CRAN
 * Latest stable version: GitHub, `master` branch
 * Bleeding-edge version: GitHub, `develop` branch

### CRAN

For the latest CRAN version:

```r
install.packages("beautier")
```

### GitHub, `master` branch

For the latest stable version: 

```r
remotes::install_github("ropensci/beautier")
```

### GitHub, `develop` branch

For the bleeding-edge version: 

```r
remotes::install_github("ropensci/beautier", ref = "develop")
```

## [FAQ](doc/faq.md)

See [FAQ](doc/faq.md).

## Supported

This works, and the interface is unlikely to change.

 * 1 DNA alignment
 * Site models:
    * JC69
    * HKY
    * TN93
    * GTR
 * Clock models:
    * Strickt
    * Relaxed log-normal
 * Tree models:
    * Yule
    * Birth-Death
    * Coalescent Bayesian Skyline 
    * Coalescent Constant Population
    * Coalescent Exponential Population
 * Handle missing data: simply use a dash (´-´) as a sequence
   in a FASTA file

## Experimental

This works partially, and the interface may change as well.

### Tip dating

The tip dates file is a file 
that needs to not have column, nor row names.
The columns need to be tab separated.

See [here](https://github.com/ropensci/beautier/blob/master/inst/extdata/G_VII_pre2003_dates_4.txt)
for an example, of which the first rows are shown here:

```
KF767106_Indonesia_1976_VII	1976
KF767104_Indonesia_1988_VII	1988
KF767105_Indonesia_1988_VII	1988
AY288998_Indonesia_1990_VII	1990
```

In the future, there probably will be a ´to_tipdates_file´ function,
to create a temporary tipdates file from a table.

## Missing features/unsupported

`beautier` cannot do everything `BEAUti` can. 

Here are some missing or (yet) unsupported features:

 * Two or more DNA alignments
 * Two or more site, clock or tree models
 * Two or more MRCA priors
 * Shared site, clock and/or tree models
 * Using an amino acid alignment
 * Support for hyper parameters
 * Clock models
   * Relaxed exponential
   * Random local
 * Tree priors
   * Calibrated Yule model
   * Coalescent Extended Bayesian Skyline
 * Initialization (this is a tab that is hidden by default in `BEAUti`)

## There is a feature I miss

See [CONTRIBUTING](CONTRIBUTING.md), at `Submitting use cases`

## I want to collaborate

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting code'

## I think I have found a bug

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting bugs' 

## There's something else I want to say

Sure, just add an Issue. Or send an email.

## External links

 * [BEAST2 GitHub](https://github.com/CompEvol/beast2)

## References

Article about `babette`:

 * Bilderbeek, Richèl JC, and Rampal S. Etienne. "`babette`: BEAUti 2, BEAST 2 and Tracer for R." Methods in Ecology and Evolution (2018). https://doi.org/10.1111/2041-210X.13032

FASTA files `anthus_aco.fas` and `anthus_nd2.fas` from:
 
 * Van Els, Paul, and Heraldo V. Norambuena. "A revision of species limits in Neotropical pipits Anthus based on multilocus genetic and vocal data." Ibis.

FASTA file `G_VII_pre2003_msa.fas` from:

 * Durr, PA; Wibowo, MH; Tabbu, CR; Asmara, W; Selleck, P; Wang, J; Broz, I; Graham, K.; Dimitrov, K and Afonso, C. (in preparation). Phylodynamics of Genotype VII Newcastle disease virus in Indonesia. 

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)


# News

Newest versions at top.

## beautier 2.6.3 (unreleased)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * None

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.6.2 (2021-07-24)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * `check_empty_beautier_folder` works correctly on Windows

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.6.1 (2021-07-23)

### NEW FEATURES

 * Added `check_empty_beautier_folder` to make sure
   no temporary files are created

### MINOR IMPROVEMENTS

 * None

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * Removed all deprecated functions and function arguments
   that warn the user in the current CRAN version (v2.6)
 * Internal functions `clock_model_to_xml_prior_distr`,
   `tipdate_taxa_to_xml_tree`, `tree_priors_to_xml_operators`
   work with `infererence_model`,
   deprecated other arguments
 * Deprecated internal function `clock_models_to_xml_prior_distr`,
   redirect user to `clock_model_to_xml_prior_distr`
 * Deprecated internal function `tree_priors_to_xml_operators`,
   redirect user to `tree_prior_to_xml_operators`
 * Deprecated `clock_models_to_xml_state`,
   redirect user to `clock_model_to_xml_state`
 * Deprecated `clock_models_to_xml_operators`,
   redirect user to `clock_model_to_xml_operators`
 * Deprecated `clock_models_to_xml_tracelog`,
   redirect user to `clock_model_to_xml_tracelog`
 * Deprecated `mrca_priors_to_xml_tracelog`,
   redirect user to `mrca_prior_to_xml_tracelog`
 * Deprecated `tree_models_to_xml_tracelog`,
   redirect user to `tree_model_to_tracelog_xml`

## beautier 2.6 (2021-05-22)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Make sure that no files are created in `~/.cache`, fix #129

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.5.2 (2021-05-22)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Removed `LazyData: true` from DESCRIPTION

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.5.1 (2021-05-22)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * GitHub Actions scheduler tests package once per week

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.5 (2021-05-14)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Use GitHub Actions to check build
 * Check more setups, thanks @GaryNapier
 * Better synergy with `system2`

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * `phylo_to_xml_state` removed 
 * `site_models_to_xml_state` gives a deprecation message
 * `tree_priors_to_xml_state` gives a deprecation message
 * `mrca_priors_to_xml_state` has the arguments `mrca_prior`
   and `has_non_strict_clock_model` replaced by `inference_model`

## beautier 2.4 (2020-10-15)

### NEW FEATURES

 * Add `inference_models` vignette

### MINOR IMPROVEMENTS

 * Added `create_ns_inference_model`
 * No `testthat` tests in examples

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.3.7 (2020-08-05)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * Internal functions have better names
 * More internal functions use an `inference_model` (instead of 
   multiple elements of it)

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.3.6 (2020-04-21)

### NEW FEATURES

 * None

### MINOR IMPROVEMENTS

 * None

### BUG FIXES

 * Fix bug when using tipdating, thanks Katherine S. Walter, @ksw9

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.3.5 (2020-02-25)

### NEW FEATURES

 * Add `get_mcmc_filenames` and `get_inference_model_filenames`
 * Add `rename_mcmc_filenames` and `rename_inference_model_filenames`
 * Add rename functions `get_remove_dir_fun`, `get_remove_hex_fun` and `get_replace_dir_fun`

### MINOR IMPROVEMENTS

 * Link to https of BEAST2

### BUG FIXES

 * None

### DEPRECATED AND DEFUNCT

 * None

## beautier 2.3.4 (2020-01-15)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Lost dependency on orphaned `geiger` package

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * Moved `fasta_to_phylo` and `fastas_to_phylos` to repository at
    `https://github.com/richelbilderbeek/ribir` to lose dependency on
    orphaned `geiger` package

## beautier 2.3.3 (2019-12-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Implemented all @lintr-bot's suggestions

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 2.3.2 (2019-11-30)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * Remove files for CRAN submission

### DEPRECATED AND DEFUNCT

  * None

## beautier 2.3.1 (2019-11-29)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Use `testthat` in documentation, instead of the less known `testit`

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 2.3 (2019-10-27)

### NEW FEATURES

  * Added `create_tracelog`, `create_screenlog`, `create_treelog`, that
    indicate where BEAST2 will store its output file
  * `create_mcmc` has all elements that BEAUti has
  * Added testing functions, that use a short MCMC chain length and/or
    a simple inference model: `create_test_mcmc`, `create_test_tracelog`,
    `create_test_screenlog`, `create_test_treelog`
    and `create_test_inference_model`

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None


## beautier 2.2.4 (2019-10-15)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Added `create_beast2_input_from_model` function

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None



## beautier 2.2.3 (2019-09-10)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Added `check_mcmc_nested_sampling` function

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 2.2.2 (2019-09-10)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Exported and documented more `is_` functions
  * Added `check_file_exists` function
  * Added `is_one_int` and `is_one_double` functions

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 2.2.1 (2019-03-08)

### NEW FEATURES

  * `beautier` has passed rOpenSci peer review
  * `beautier` is on CRAN

### MINOR IMPROVEMENTS

  * Exported multiple `is_x` and `check_x` functions

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 1.15 (2018-11-10)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * Correct BEAST2 files are created when using a RLN clock
    model, with an MRCA prior with a distribution, thanks @rscherrer
    and Jana Riederer

### DEPRECATED AND DEFUNCT

  * Support for multiple alignments, site models, clock models,
    tree models has decreased 

## beautier 1.14.2 (2018-10-30)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tested to build under macOS

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 1.14.1 (2018-10-29)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Documented the simplified interface for parameters better
  * Use the simplified interface for parameters as defaults more often

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 1.14.0 (2018-10-26)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Simplified interface for parameters:

```
# Old
distr <- create_distr_poisson(id = 1, lambda = create_lambda_param(value = 1.2))
# Added
distr <- create_distr_poisson(id = 1, lambda = 1.2)
```

### BUG FIXES

  * When using MRCA priors with/without monophyly with/without
    a distribution, resulted in incorrect BEAST2 `.xml` files

### DEPRECATED AND DEFUNCT

  * Cannot set parameters to be estimated, when the resulting
    BEAST2 `.xml` would be incomplete

## beautier 1.13.4 (2018-09-11)

### NEW FEATURES

  * Support for [a Nested-Sampling MCMC run](https://github.com/BEAST2-Dev/nested-sampling)

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 1.13.3 (2018-05-17)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tagged for [the academic article](https://github.com/ropensci/babette_article) about `babette`

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## beautier 1.13.2 (2018-04-05)

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

This GitHub follows the [Contributor Covenant Code of Conduct](CODE_OF_CONDUCT.md).

 * For questions, you can create an Issue
 * Code changes go via Pull Requests

## Which package to contribute to?

`beautier` is part of the `babette` package suite,
which consists out of five packages.
Here is how to determine which package is best suited for your contribution:

If you want to contribute to how BEAST2 is run,
go to [beautier](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md).

If you want to contribute to how BEAST2 output is parsed,
go to [tracerer](https://github.com/ropensci/tracerer/blob/master/CONTRIBUTING.md)

If you want to contribute regarding the BEAST2 package management,
go to [mauricer](https://github.com/ropensci/mauricer/blob/master/CONTRIBUTING.md)

If you want to contribute with an overarching idea,
go to [babette](https://github.com/ropensci/babette/blob/master/CONTRIBUTING.md).

If you want to contribute to the creation of BEAST2 XML input files, 
you are at the right spot :-) 

## Submitting use cases

Use cases within the default BEAUti environment are welcomed.

Please send all that is needed to reproduce the use case:

 * the alignment file
 * screenshots of BEAUti settings you've changed
 * the resulting XML file
 * (optional) the desired call to `beautier`

BEAUti plugins are not supported (for now). See 'Submitting code'
if you'd like `beautier` to do so.

You can do so by:

 * Add an Issue
 * Send @richelbilderbeek an email (@richelbilderbeek will make an Issue of it)

## Submitting code

Submitted code should follow these quality guidelines:

 * All tests pass cleanly/silently
 * Code coverage above 95%
 * Coding style should follow the default style by `lintr`

These are all checked by Travis CI when submitting
a Pull Request. 

For BEAUti plugins, there needs to be added:

 * Enough documentation so an example use case can be reproduced
 * Stating that maintenance will be done by the Collaborator

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

To get started working on `beautier` do:

```
git clone https://github.com/ropensci/beautier
```

Development is done on the `develop` branch. 
To download and checkout the `develop` branch, 
first go into the `beautier` folder (`cd beautier`), then do:

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
utils::sessionInfo()
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

`beautier` pictures.

## How did you convert the fuzzy white background to one single color?

```
convert butterfly.png -fuzz 15% -fill white -opaque white butterfly_mono_background.png
convert butterfly_mono_background.png -background white -alpha remove butterfly_mono_background_2.png
```
# Demo

See [the 'Demo' vignette](https://github.com/ropensci/beautier/blob/master/vignettes/demo.Rmd).


# Road map

## Future milestones

Version  | Goal
---------|----------------------------------------------------------------------
?        | Use S3 classes instead of a lists
?        | Use S4 classes instead of S3 classes

There are some minor things that can be added, see 'How can I indicate a feature that I miss?'
to request so.

## Past milestones

Version  | Goal
---------|----------------------------------------------------------------------
`v2.0`   | allow tip dating
`v1.12`  | ~~support two alignments~~
`v1.12.1`| simpler interface for specifying a fixed crown age
`v1.12.2`| remove cyclic dependencies on other related packages
`v1.13`  | allow one or two MRCA nodes (also called 'calibration nodes')
# FAQ

## Shouldn't the slogan be 'beautier: BEAUti 2 for R'?

That slogan would indeed be more precise. That
extra precision would come at the cost of 
readability (the extra '2'). As there is no `BEAUti 1`,
there is no possible confusion and that extra number 
does not add extra information.

## Which version of BEAUti do you use as a guideline?

The first BEAST2 XML files created by `beautier`
followed BEAST2 v2.4. `beautier` follows the
BEAST2 versions, which is now at v2.6.0.

The BEAST2 version actually used by `babette`
can be found in the [beastier::install_beast2](https://github.com/ropensci/beastier/blob/master/R/install_beast2.R) function.

## Why does AppVeyor only check the `master` branch?

Because `ropensci` does not have AppVeyor.

To do check for Windows, 
[the beautier_on_windows repo](https://github.com/richelbilderbeek/beautier_on_windows)
is created. That repo only checks the `master` branch of `beautier`.

## What's the [road map](road_map.md)?

See [road map](road_map.md).

## How can I indicate a feature that I miss?

Submit an Issue.

## How can I submit code?

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting code'

## How can I submit a bug?

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting bugs' 

## How can I indicate something else?

Submit an Issue. Or send an email to Richèl Bilderbeek.

## How do I reference to this work?

Cite:

```
Bilderbeek, Richèl JC, and Rampal S. Etienne. "babette: BEAUti 2, BEAST 2 and Tracer for R." Methods in Ecology and Evolution (2018).
```

or

```
@article{bilderbeek2018babette,
  title={babette: BEAUti 2, BEAST 2 and Tracer for R},
  author={Bilderbeek, Richèl JC and Etienne, Rampal S},
  journal={Methods in Ecology and Evolution},
  year={2018},
  publisher={Wiley Online Library}
}
```

### Are there any related packages?

 * [lumier](https://github.com/ropensci/lumier): Shiny app to help create the function call needed
 * [BEASTmasteR](https://github.com/nmatzke/BEASTmasteR): tip-dating analyses using fossils as dated terminal taxa
 * [BEASTifier](https://github.com/josephwb/BEASTifier): generate BEAST input files from a NEXUS file, similar to [beautier](https://github.com/ropensci/beautier)
 * [RBeast](https://github.com/beast-dev/RBeast): misc other things

## What is the idea behind the logo?

The butterfly symbolizes beauty.
Then it was combined with an R logo. 

## What are the FASTA files?

Filename               |Reference
-----------------------|------------
`anthus_aco.fas`       |[1]
`anthus_nd2.fas`       |[1]
`G_VII_pre2003_msa.fas`|[2]
Others                 |Artificial

 * [1] Van Els, Paul, and Heraldo V. Norambuena. "A revision of species limits in Neotropical pipits Anthus based on multilocus genetic and vocal data." Ibis.
 * [2] Durr, PA; Wibowo, MH; Tabbu, CR; Asmara, W; Selleck, P; Wang, J; Broz, I; Graham, K.; Dimitrov, K and Afonso, C. (in preparation). Phylodynamics of Genotype VII Newcastle disease virus in Indonesia.
 
Thanks to Peter A. Durr and Paul van Els for supplying the FASTA files.

## Why are the functions prefixed with `create_`?

Or, why is this chosen:

```{r}
out <- create_beast2_input(
  "alignment.fas",
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = create_exp_distr()    
  )
)
```

over this:

```{r}
out <- create_beast2_input(
  "alignment.fas",
  tree_prior = yule_tree_prior(
    birth_rate_distr = exp_distr()    
  )
)
```

Answer: because function names should start with a 
verb (see e.g. [https://style.tidyverse.org/functions.html#naming](The Tidyverse Style Guide))

## Why the name?

`beautier` is 'BEAUti for R'. 

Additionally, it is a joke that suggests `beautier` would have more beauty than `BEAUti`.
This suggestion benefits the image of author of `beautier`, who, however, thinks that
both tools are equally valuable and beautiful.

## Why the logo?

Initially, the logo was a low-tech remake of Belle, for Beauty and the Beast. 
To prevent problems with Disney, a different logo was picked.

The current logo shows a butterfly, an animal considered to be beautiful.
The butterfly is drawn by Jose Scholte, who kindly allowed her work to
be used for free, by attribution.

## `BEAUti` problems

### Not enough memory

```
./beauti
```

```
Can't start up: not enough memory
```

On Artful Aardvark, remove `-Xms256m -Xmx4g` from the `bin/beauti` file's last line. Change:

```
"$JAVA" -Dlauncher.wait.for.exit=true -Xms256m -Xmx4g -Djava.library.path="$BEAST_LIB" -Duser.language=en -cp "$BEAST_LIB/launcher.jar" beast.app.beauti.BeautiLauncher -capture $*
```

to

```
"$JAVA" -Dlauncher.wait.for.exit=true -Djava.library.path="$BEAST_LIB" -Duser.language=en -cp "$BEAST_LIB/launcher.jar" beast.app.beauti.BeautiLauncher -capture $*
```

### BEAUti requires Java version at least 8

![](beauti_requires_java_8_or_more.png)

Do:

```
sudo update-alternatives --config java
```

Pick `/usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java`:

```
There are 5 choices for the alternative java (providing /usr/bin/java).

  Selection    Path                                            Priority   Status
------------------------------------------------------------
  0            /usr/lib/jvm/java-9-openjdk-amd64/bin/java       1091      auto mode
  1            /usr/bin/gij-4.8                                 1048      manual mode
  2            /usr/bin/gij-5                                   1050      manual mode
  3            /usr/bin/gij-6                                   1060      manual mode
* 4            /usr/lib/jvm/java-8-openjdk-amd64/jre/bin/java   1081      manual mode
  5            /usr/lib/jvm/java-9-openjdk-amd64/bin/java       1091      manual mode
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

See [the 'Examples' vignette](https://github.com/ropensci/beautier/blob/master/vignettes/examples.Rmd).

---
title: "beautier demo"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{beautier demo}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)
```

# Introduction

![](beautier_logo.png)

The purpose of `beautier` is to create a valid `BEAST2` XML input file
from its function argument. In this way, a scientific pipeline using
`BEAST2` can be fully scripted, instead of using `BEAUti`'s GUI.

`beautier` is part of the `babette` package suite (website at [https://github.com/ropensci/babette](https://github.com/ropensci/babette)).
`babette` allows to use BEAST2 (and its tools) from R.

# Demonstration

First, `beautier` is loaded:

```{r load_beautier}
library(beautier)
```

A BEAST2 XML input file needs an alignment (as BEAST2 infers
phylogenies and parameters on DNA sequences). This demonstration uses
a testing FASTA file used by `beautier`:

```{r get_fasta_filename}
fasta_filename <- get_beautier_path("test_output_0.fas")
```

We can display the alignment in the file:


```{r show_alignment}
image(ape::read.FASTA(fasta_filename))
```

Specify the filename for our XML file, here I use a temporary filename,
so there is no need to clean up afterwards:

```{r create_output_filename}
output_filename <- get_beautier_tempfilename()
output_filename
```

Now we can create our XML file. We do not specify any inference model,
and just use the BEAUti default settings:

![](all_default.png)

```{r create_beast2_input_file}
create_beast2_input_file(
  fasta_filename,
  output_filename
)
```

The file indeed is a BEAST2 input file:

```{r show_beast2_input_file}
readLines(output_filename)
```

This XML input file can be read by BEAST2.

You can use `beastier` to run BEAST2 from R, see [https://github.com/ropensci/beastier](https://github.com/ropensci/beastier).
You can use `babette` to do a BEAST2 inference directly,
see [https://github.com/ropensci/babette](https://github.com/ropensci/babette).


## Cleanup

```{r cleanup}
file.remove(output_filename)
```
---
title: "Inference Models"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Inference Models}
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

![](beautier_logo.png)

`beautier` allows to specify an inference model, 
from which, together with a DNA alignment, a posterior (holding
a distribution of phylogenies and jointly-estimated parameter estimates) 
can be inferred.

The inference model entails the evolutionary model used in inference,
as well as some settings to do the inference.

This vignette shows the options how to specify an inference model,
and which options are possible.

## Getting started

First, load `beautier`:

```{r}
library(beautier)
```

Now, we'll look at the default inference model to get an idea of what 
an inference model entails.

## Default inference model

Creating a default inference model is easy:

```{r}
inference_model <- create_inference_model()
```

An inference model is a list, that uses the BEAST2 defaults. Here
are the elements of that list:

```{r}
names(inference_model)
```

As we can see, an inference model entails these elements, which
we will cover below in more detail:

 * a site model: how the alignment changes
 * a clock model: how the mutation rates differ over the branches
 * a tree prior: the speciation model
 * (optional) an MRCA prior: a 'Most Recent Common Ancestor'
 * an MCMC: the Markov Chain Monte Carlo setup
 * (optional) BEAUti options: version-specific options
 * (experimental) a tipdates filename: a filename for tip-dating

### Site models

The site model entails how the alignment changes over time.
Currently, `beautier` supplies a gamma site model. 
One element of a site model, for DNA/RNA, is the nucleotide 
substitution model ('NSM'). 

Due to historical reasons, `beautier` confuses the site model and
NSM: `beautier` has no functions with
the word 'nucleotide substitution' (nor `nsm`) `in it. Instead, it is as if
these are specific site models.

To see the available site models, use `?create_site_model` to
see a list, or use:

```{r}
get_site_model_names()
```

The simplest NSM is the JC69 NSM, which assumes all nucleotides
are substituted by one another at the same rate.
As an example, to use a gamma site model with the JC69 NSM model
in inference:

```{r}
inference_model <- create_inference_model(
  site_model = create_jc69_site_model()
)
```

### Clock models

The clock model entails how the mutation rates differ over the branches.

To see the available site models, use `?create_clock_model` to
see a list, or use:

```{r}
get_clock_model_names()
```

The simplest clock model is the strict clock model, 
which assumes all branches have the same mutation rate.
As an example, to use a strict clock model in inference:

```{r}
inference_model <- create_inference_model(
  clock_model = create_strict_clock_model()
)
```

## Tree prior

The tree prior is the tree model used in inference.
It is called 'tree prior' instead of 'tree model', as
this follow the BEAUti naming.
The tree model specifies the branching process of a tree. 

To see the available tree models, use `?create_tree_prior` to
see a list, or use:

```{r}
get_tree_prior_names()
```

The simplest tree model is the Yule (aka pure-birth) tree model, 
which assumes that branching events occur at a constant rate,
and there are no extinctions.
As an example, to use a Yule tree model in inference:

```{r}
inference_model <- create_inference_model(
  tree_prior = create_yule_tree_prior()
)
```

## (optional) MRCA prior: a 'Most Recent Common Ancestor'

With the MRCA ('Most Recent Common Ancestor') prior, one can
specify which tips share a common ancestor.

```{r}
# The alignmet
fasta_filename <- get_beautier_path("anthus_aco.fas")

# The alignment's ID
alignment_id <- get_alignment_id(fasta_filename)

# Get the first two taxa's names
taxa_names <- get_taxa_names(fasta_filename)[1:2]

# Specify that the first two taxa share a common ancestor
mrca_prior <- create_mrca_prior(
  alignment_id = alignment_id,
  taxa_names = taxa_names
)

# Use the MRCA prior in inference
inference_model <- create_inference_model(
  mrca_prior = mrca_prior
)
```

## MCMC: the Markov Chain Monte Carlo setup

The MCMC ('Markov Chain Monte Carlo') specifies how the
inference algorithm does its work.

The available MCMC's can be found using `?create_mcmc`
and are:

 * `create_mcmc`: regular MCMC
 * `create_test_mcmc`: shorter regular MCMC, to be used in testing
 * `create_ns_mcmc`: MCMC to estimate a marginal likelihood 
   using nested sampling

## (optional) BEAUti options: version-specific options

The BEAUti options entail version-specific options
to store an inference model as a BEAST2 XML input file.

The available BEAUti options can be found using `?create_beauti_options`
and are:

 * `create_beauti_options_v2_4`: BEAUti v2.4
 * `create_beauti_options_v2_6`: BEAUti v2.6

Using a specific version for an inference:

```{r}
inference_model <- create_inference_model(
  beauti_options = create_beauti_options_v2_4()
)
```

## (experimental) a tipdates filename: a filename for tip-dating

A tipdates filename is an experimental feature for tip-dating:

```{r}
inference_model <- create_inference_model(
  tipdates_filename = get_beautier_path("G_VII_pre2003_dates_4.txt")
)
```

The tipdates filename and the alignment must be compatible.
Here is an example:

```{r}
output_filename <- get_beautier_tempfilename()

create_beast2_input_file_from_model(
  input_filename = get_beautier_path("G_VII_pre2003_msa.fas"),
  inference_model = inference_model,
  output_filename = output_filename
)
# Cleanup
file.remove(output_filename)
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
  comment = ""
)
```

# Introduction

![](beautier_logo.png)

The purpose of `beautier` is to create a valid `BEAST2` XML input file
from its function argument. In this way, a scientific pipeline using
`BEAST2` can be fully scripted, instead of using `BEAUti`'s GUI.

`beautier` is part of the `babette` package suite (website at [https://github.com/ropensci/babette](https://github.com/ropensci/babette)).
`babette` allows to use BEAST2 (and its tools) from R.

# Getting started

For all examples, do load `beautier`:

```{r load_beautier}
library(beautier)
```

Each example shows a picture of a BEAUti dialog to achieve the same.
BEAUti is part of the BEAST2 tool suite and it's a GUI to create BEAST2
input files. `beautier` is an R package to supplement BEAUti, by providing
to do the same from an R script.

Each example reads the alignment from a FASTA file called `anthus_aco_sub.fas`,
which is part of the files supplied with `beautier`:

```{r show_input_file}
input_filename <- get_beautier_path("anthus_aco_sub.fas")
```

In this vignette, the generated BEAST2 XML is shown.
Use `create_beast2_input_file` to save the resulting
XML directly to file instead.

## Example #1: all default

![](all_default.png)

Using all default settings, only specify a DNA alignment.

```{r example_1}
create_beast2_input(
  input_filename
)
```

All other parameters are set to their defaults, as in BEAUti.

## Example #2: JC69 site model

![](jc69_2_4.png)

```{r example_2}
create_beast2_input(
  input_filename,
  site_model = create_jc69_site_model()
)
```

## Example #3: Relaxed clock log normal

![](rln_2_4.png)

```{r example_3}
create_beast2_input(
  input_filename,
  clock_model = create_rln_clock_model()
)
```

## Example #4: Birth-Death tree prior

![](bd_2_4.png)

```{r example_4}
create_beast2_input(
  input_filename,
  tree_prior = create_bd_tree_prior()
)
```

## Example #5: Yule tree prior with a normally distributed birth rate

![](birth_rate_normal_2_4.png)

```{r example_5}
create_beast2_input(
  input_filename,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = create_normal_distr()
  )
)
```

Thanks to Yacine Ben Chehida for this use case

## Example #6: HKY site model with a non-zero proportion of invariants

![](hky_prop_invariant_0_5_2_4.png)

```{r example_6}
create_beast2_input(
  input_filename,
  site_model = create_hky_site_model(
    gamma_site_model = create_gamma_site_model(prop_invariant = 0.5)
  )
)
```

Thanks to Yacine Ben Chehida for this use case

## Example #7: Strict clock with a known clock rate

![](strict_clock_rate_0_5_2_4.png)

```{r example_7}
create_beast2_input(
  input_filename,
  clock_model = create_strict_clock_model(
    clock_rate_param = 0.5
  )
)
```

Thanks to Paul van Els and Yacine Ben Chehida for this use case.

## Example #8: Use MRCA prior

![](mrca_prior_all.png)

```{r example_8}
create_beast2_input(
  input_filename,
  mrca_prior = create_mrca_prior()
)
```

## Example #9: Use MRCA prior to specify a close-to-fixed crown age

![](mrca_prior_crown_age.png)

With an MRCA  prior, it is possible
to specify a close-to-fixed crown age:

```{r example_9}
crown_age <- 15

create_beast2_input(
  input_filename,
  mrca_prior = create_mrca_prior(
    is_monophyletic = TRUE,
    mrca_distr = create_normal_distr(
      mean = crown_age,
      sigma = 0.001
    )
  )
)
```
---
title: "beautier demo"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{beautier demo}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = ""
)
```

# Introduction

![](beautier_logo.png)

The purpose of `beautier` is to create a valid `BEAST2` XML input file
from its function argument. In this way, a scientific pipeline using
`BEAST2` can be fully scripted, instead of using `BEAUti`'s GUI.

`beautier` is part of the `babette` package suite (website at [https://github.com/ropensci/babette](https://github.com/ropensci/babette)).
`babette` allows to use BEAST2 (and its tools) from R.

# Demonstration

First, `beautier` is loaded:

```{r load_beautier}
library(beautier)
```

A BEAST2 XML input file needs an alignment (as BEAST2 infers
phylogenies and parameters on DNA sequences). This demonstration uses
a testing FASTA file used by `beautier`:

```{r get_fasta_filename}
fasta_filename <- get_beautier_path("test_output_0.fas")
```

We can display the alignment in the file:


```{r show_alignment}
image(ape::read.FASTA(fasta_filename))
```

Specify the filename for our XML file, here I use a temporary filename,
so there is no need to clean up afterwards:

```{r create_output_filename}
output_filename <- get_beautier_tempfilename()
output_filename
```

Now we can create our XML file. We do not specify any inference model,
and just use the BEAUti default settings:

![](all_default.png)

```{r create_beast2_input_file}
create_beast2_input_file(
  fasta_filename,
  output_filename
)
```

The file indeed is a BEAST2 input file:

```{r show_beast2_input_file}
readLines(output_filename)
```

This XML input file can be read by BEAST2.

You can use `beastier` to run BEAST2 from R, see [https://github.com/ropensci/beastier](https://github.com/ropensci/beastier).
You can use `babette` to do a BEAST2 inference directly,
see [https://github.com/ropensci/babette](https://github.com/ropensci/babette).

## Clean up

```{r show_beast2_input_file}
file.remove(output_filename)
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
  comment = ""
)
```

# Introduction

![](beautier_logo.png)

The purpose of `beautier` is to create a valid `BEAST2` XML input file
from its function argument. In this way, a scientific pipeline using
`BEAST2` can be fully scripted, instead of using `BEAUti`'s GUI.

`beautier` is part of the `babette` package suite (website at [https://github.com/ropensci/babette](https://github.com/ropensci/babette)).
`babette` allows to use BEAST2 (and its tools) from R.

# Getting started

For all examples, do load `beautier`:

```{r load_beautier}
library(beautier)
```

Each example shows a picture of a BEAUti dialog to achieve the same.
BEAUti is part of the BEAST2 tool suite and it's a GUI to create BEAST2
input files. `beautier` is an R package to supplement BEAUti, by providing
to do the same from an R script.

Each example reads the alignment from a FASTA file called `anthus_aco_sub.fas`,
which is part of the files supplied with `beautier`:

```{r show_input_file}
input_filename <- get_beautier_path("anthus_aco_sub.fas")
```

In this vignette, the generated BEAST2 XML is shown.
Use `create_beast2_input_file` to save the resulting
XML directly to file instead.

## Example #1: all default

![](all_default.png)

Using all default settings, only specify a DNA alignment.

```{r example_1}
create_beast2_input(
  input_filename
)
```

All other parameters are set to their defaults, as in BEAUti.

## Example #2: JC69 site model

![](jc69_2_4.png)

```{r example_2}
create_beast2_input(
  input_filename,
  site_model = create_jc69_site_model()
)
```

## Example #3: Relaxed clock log normal

![](rln_2_4.png)

```{r example_3}
create_beast2_input(
  input_filename,
  clock_model = create_rln_clock_model()
)
```

## Example #4: Birth-Death tree prior

![](bd_2_4.png)

```{r example_4}
create_beast2_input(
  input_filename,
  tree_prior = create_bd_tree_prior()
)
```

## Example #5: Yule tree prior with a normally distributed birth rate

![](birth_rate_normal_2_4.png)

```{r example_5}
create_beast2_input(
  input_filename,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = create_normal_distr()
  )
)
```

Thanks to Yacine Ben Chehida for this use case

## Example #6: HKY site model with a non-zero proportion of invariants

![](hky_prop_invariant_0_5_2_4.png)

```{r example_6}
create_beast2_input(
  input_filename,
  site_model = create_hky_site_model(
    gamma_site_model = create_gamma_site_model(prop_invariant = 0.5)
  )
)
```

Thanks to Yacine Ben Chehida for this use case

## Example #7: Strict clock with a known clock rate

![](strict_clock_rate_0_5_2_4.png)

```{r example_7}
create_beast2_input(
  input_filename,
  clock_model = create_strict_clock_model(
    clock_rate_param = 0.5
  )
)
```

Thanks to Paul van Els and Yacine Ben Chehida for this use case.

## Example #8: Use MRCA prior

![](mrca_prior_all.png)

```{r example_8}
create_beast2_input(
  input_filename,
  mrca_prior = create_mrca_prior()
)
```

## Example #9: Use MRCA prior to specify a close-to-fixed crown age

![](mrca_prior_crown_age.png)

With an MRCA  prior, it is possible
to specify a close-to-fixed crown age:

```{r example_9}
crown_age <- 15

create_beast2_input(
  input_filename,
  mrca_prior = create_mrca_prior(
    is_monophyletic = TRUE,
    mrca_distr = create_normal_distr(
      mean = crown_age,
      sigma = 0.001
    )
  )
)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mcmc.R
\name{check_mcmc_values}
\alias{check_mcmc_values}
\title{Check if the MCMC has the list elements with valid values
for being a valid MCMC object.}
\usage{
check_mcmc_values(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
nothing
}
\description{
Calls \code{stop} if a value is invalid
}
\seealso{
Use \code{\link{create_mcmc}} to create a valid MCMC
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_strict_clock_model.R
\name{has_strict_clock_model}
\alias{has_strict_clock_model}
\title{Determine if the \code{inference_model} uses a strict clock model.}
\usage{
has_strict_clock_model(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
TRUE if the \code{inference_model} uses a strict clock model,
FALSE otherwise
}
\description{
Determine if the \code{inference_model} uses a strict clock model
}
\examples{
# Yes, has a strict clock model
has_strict_clock_model(
  create_inference_model(clock_model = create_strict_clock_model())
)

# No strict clock model
has_strict_clock_model(
  create_inference_model(clock_model = create_rln_clock_model())
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_init_site_models.R
\name{are_init_site_models}
\alias{are_init_site_models}
\title{Determine if x consists out of initialized site_models objects}
\usage{
are_init_site_models(x)
}
\arguments{
\item{x}{the object to check if it consists out of
initialized site_models objects}
}
\value{
TRUE if x, or all elements of x, are initialized site_model objects
}
\description{
Determine if x consists out of initialized site_models objects
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tree_prior.R
\name{create_ccp_tree_prior}
\alias{create_ccp_tree_prior}
\alias{create_tree_prior_ccp}
\title{Create a Coalescent Constant Population tree prior}
\usage{
create_ccp_tree_prior(
  id = NA,
  pop_size_distr = beautier::create_one_div_x_distr(value = 0.3)
)
}
\arguments{
\item{id}{the ID of the alignment}

\item{pop_size_distr}{the population distribution,
as created by a \code{\link{create_distr}} function}
}
\value{
a Coalescent Constant Population tree_prior
}
\description{
Create a Coalescent Constant Population tree prior
}
\examples{
ccp_tree_prior <- create_ccp_tree_prior()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = ccp_tree_prior
)
file.remove(beast2_input_file)
}
\seealso{
An alignment ID can be extracted from
  its FASTA filename using \code{\link{get_alignment_id}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_last_xml_closing_tag_line.R
\name{find_last_xml_closing_tag_line}
\alias{find_last_xml_closing_tag_line}
\title{Find the highest line number of a section's closing tag}
\usage{
find_last_xml_closing_tag_line(lines, section)
}
\arguments{
\item{lines}{the lines of an XML text}

\item{section}{the name of the XML section}
}
\value{
the line number's index (which is 1 for the first line) if the
  opening tag is found, else NA
}
\description{
Find the highest line number of a section's closing tag
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_model.R
\name{create_jc69_site_model}
\alias{create_jc69_site_model}
\alias{create_site_model_jc69}
\title{Create a JC69 site model}
\usage{
create_jc69_site_model(id = NA, gamma_site_model = create_gamma_site_model())
}
\arguments{
\item{id}{the IDs of the alignment (can be extracted from
the FASTA filename using \code{\link{get_alignment_id}})}

\item{gamma_site_model}{a gamma site model, as created
by \code{\link{create_gamma_site_model}}}
}
\value{
a JC69 site_model
}
\description{
Create a JC69 site model
}
\examples{
 jc69_site_model <- create_jc69_site_model()

 output_filename <- get_beautier_tempfilename()
 create_beast2_input_file(
   input_filename = get_fasta_filename(),
   output_filename = output_filename,
   site_model = jc69_site_model
 )
file.remove(output_filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_rate_ac}
\alias{parameter_to_xml_rate_ac}
\title{Internal function}
\usage{
parameter_to_xml_rate_ac(
  parameter,
  beauti_options = create_beauti_options(),
  which_name = "state_node"
)
}
\arguments{
\item{parameter}{a 'rate AC' parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_rate_ac_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{which_name}{the name, can be \code{state_node} or \code{rate_name}}
}
\value{
the parameter as XML text
}
\description{
Converts a 'rate AC' parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_rate_cg}
\alias{parameter_to_xml_rate_cg}
\title{Internal function}
\usage{
parameter_to_xml_rate_cg(
  parameter,
  beauti_options = create_beauti_options(),
  which_name = "state_node"
)
}
\arguments{
\item{parameter}{a 'rate CG' parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_rate_cg_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{which_name}{the name, can be \code{state_node} or \code{rate_name}}
}
\value{
the parameter as XML text
}
\description{
Converts a 'rate CG' parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml}
\alias{parameter_to_xml}
\title{Internal function}
\usage{
parameter_to_xml(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a parameter,
as created by \code{\link{create_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a parameter to XML
}
\examples{
parameter_to_xml(create_alpha_param(id = 1))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_models_from_names.R
\name{create_site_models_from_names}
\alias{create_site_models_from_names}
\title{Create site models from their names}
\usage{
create_site_models_from_names(site_model_names)
}
\arguments{
\item{site_model_names}{one or more names of a site model,
must be name among those returned by \code{\link{get_site_model_names}}}
}
\value{
one or more site models
}
\description{
Create site models from their names
}
\examples{
create_site_models_from_names(get_site_model_names())
}
\seealso{
Use \link{create_site_model} to create a site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_clock_models.R
\name{init_strict_clock_model}
\alias{init_strict_clock_model}
\title{Initializes a strict clock model}
\usage{
init_strict_clock_model(strict_clock_model, distr_id = 0, param_id = 0)
}
\arguments{
\item{strict_clock_model}{a strict clock model,
as returned by \code{\link{create_strict_clock_model}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized strict clock model
}
\description{
Initializes a strict clock model
}
\examples{

strict_clock_model <- create_strict_clock_model()
# FALSE: not yet initialized
is_init_strict_clock_model(strict_clock_model)
strict_clock_model <- init_strict_clock_model(strict_clock_model)
# TRUE: initialized
is_init_strict_clock_model(strict_clock_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tree_prior.R
\name{create_bd_tree_prior}
\alias{create_bd_tree_prior}
\alias{create_tree_prior_bd}
\title{Create a Birth-Death tree prior}
\usage{
create_bd_tree_prior(
  id = NA,
  birth_rate_distr = create_uniform_distr(),
  death_rate_distr = create_uniform_distr()
)
}
\arguments{
\item{id}{the ID of the alignment}

\item{birth_rate_distr}{the birth rate distribution,
as created by a \code{\link{create_distr}} function}

\item{death_rate_distr}{the death rate distribution,
as created by a \code{\link{create_distr}} function}
}
\value{
a Birth-Death tree_prior
}
\description{
Create a Birth-Death tree prior
}
\examples{
bd_tree_prior <- create_bd_tree_prior()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = bd_tree_prior
)
file.remove(beast2_input_file)

bd_tree_prior_exp <- create_bd_tree_prior(
  birth_rate_distr = create_exp_distr()
)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = bd_tree_prior_exp
)
file.remove(beast2_input_file)
}
\seealso{
An alignment ID can be extracted from
  its FASTA filename using \code{\link{get_alignment_id}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_operators.R
\name{create_beast2_input_operators}
\alias{create_beast2_input_operators}
\title{Creates the operators section of a BEAST2 XML parameter file}
\usage{
create_beast2_input_operators(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Creates the operators section of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site_models_n_params.R
\name{get_site_models_n_params}
\alias{get_site_models_n_params}
\title{Get the number of distributions one or more site models have}
\usage{
get_site_models_n_params(site_models)
}
\arguments{
\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}
}
\value{
the number of parameters the site models have
}
\description{
Get the number of distributions one or more site models have
}
\examples{
  testit::assert(
    get_site_models_n_params(list(create_gtr_site_model())) == 10
  )
  testit::assert(
    get_site_models_n_params(list(create_hky_site_model())) == 2
  )
  testit::assert(
    get_site_models_n_params(list(create_jc69_site_model())) == 0
  )
  testit::assert(
    get_site_models_n_params(list(create_tn93_site_model())) == 4
  )
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_models_to_xml_tracelog.R
\name{clock_models_to_xml_tracelog}
\alias{clock_models_to_xml_tracelog}
\title{Deprecated internal function}
\usage{
clock_models_to_xml_tracelog(clock_models, mrca_priors = NA)
}
\arguments{
\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}

\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}
}
\value{
a character vector of XML strings
}
\description{
Creates the clock models' XML for the tracelog section
}
\examples{
# <logger id="tracelog" ...>
#'   # Here
# </logger>
}
\seealso{
the complete tracelog section is created
  by \code{\link{create_tracelog_xml}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_tree_priors.R
\name{are_tree_priors}
\alias{are_tree_priors}
\title{Determine if x consists out of tree_priors objects}
\usage{
are_tree_priors(x)
}
\arguments{
\item{x}{the object to check if it consists out of tree_priors objects}
}
\value{
TRUE if x, or all elements of x, are tree_prior objects
}
\description{
Determine if x consists out of tree_priors objects
}
\examples{

yule_tree_prior <- create_yule_tree_prior()
bd_tree_prior <- create_bd_tree_prior()
both_tree_priors <- list(yule_tree_prior, bd_tree_prior)
# TRUE
are_tree_priors(yule_tree_prior)
# TRUE
are_tree_priors(bd_tree_prior)
# TRUE
are_tree_priors(both_tree_priors)
}
\seealso{
Use \link{create_yule_tree_prior} to create a Yule tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_clock_models.R
\name{are_rln_clock_models}
\alias{are_rln_clock_models}
\title{Are the clock models Relaxed Log-Normal clock models?}
\usage{
are_rln_clock_models(clock_models)
}
\arguments{
\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}
}
\value{
vector of booleans with the same length
  as the number of clock models in \code{clock_models}.
  Each nth element is TRUE if the nth element
  in \code{clock_models} is a relaxed log-normal
  clock model, FALSE otherwise
}
\description{
Are the clock models Relaxed Log-Normal clock models?
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_s_param}
\alias{is_s_param}
\title{Determine if the object is a valid
s parameter}
\usage{
is_s_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
s parameter}
}
\value{
TRUE if x is a valid s parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
s parameter
}
\examples{

is_s_param(create_alpha_param())
is_s_param(create_beta_param())
is_s_param(create_clock_rate_param())
is_s_param(create_kappa_1_param())
is_s_param(create_kappa_2_param())
is_s_param(create_lambda_param())
is_s_param(create_m_param())
is_s_param(create_mean_param())
is_s_param(create_mu_param())
is_s_param(create_rate_ac_param())
is_s_param(create_rate_ag_param())
is_s_param(create_rate_at_param())
is_s_param(create_rate_cg_param())
is_s_param(create_rate_ct_param())
is_s_param(create_rate_gt_param())
is_s_param(create_s_param())
is_s_param(create_scale_param())
is_s_param(create_sigma_param())

is_s_param(NA)
is_s_param(NULL)
is_s_param("nonsense")
is_s_param(create_jc69_site_model())
is_s_param(create_strict_clock_model())
is_s_param(create_yule_tree_prior())
is_s_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_rate_ag_param}
\alias{create_rate_ag_param}
\alias{create_param_rate_ag}
\title{Create a parameter called 'rate AG'}
\usage{
create_rate_ag_param(id = NA, estimate = TRUE, value = "1.0", lower = "0.0")
}
\arguments{
\item{id}{the parameter's ID}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}

\item{value}{value of the parameter}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}
}
\value{
a parameter called 'rate AG'
}
\description{
Create a parameter called 'rate AG'
}
\examples{
# Create parameter
rate_ag_param <- create_rate_ag_param(value = 1, estimate = FALSE)

# Use the parameter to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  site_model = create_gtr_site_model(
    rate_ag_param = rate_ag_param
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_run.R
\name{create_beast2_input_run}
\alias{create_beast2_input_run}
\title{Creates the '\code{run}' section of a BEAST2 XML parameter file}
\usage{
create_beast2_input_run(
  input_filename,
  inference_model = beautier::create_inference_model()
)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Creates the '\code{run}' section of a BEAST2 XML parameter file,
without being indented.
}
\details{
The \code{run} tag has these elements:
\preformatted{
   <run[...]>
       <state[...]>
       [...]
       </state>
       <init[...]>
       [...]
       </init>
       <distribution[...]>
       [...]
       </distribution>
       [operator ids]
       [loggers]
    </run>
}
}
\seealso{
Use \link{create_beast2_input_state}
to create the XML text of the \code{state} tag.
Use \link{create_beast2_input_init}
to create the XML text of the \code{init} tag.
Use \link{create_beast2_input_distr}
to create the XML text of the \code{distribution} tag.
Use \link{create_beast2_input_operators}
to create the XML text of the \code{[operator ids]} section.
Use \link{create_loggers_xml}
to create the XML text of the \code{[loggers]} part.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_models_to_xml_operators.R
\name{clock_models_to_xml_operators}
\alias{clock_models_to_xml_operators}
\title{Deprecated}
\usage{
clock_models_to_xml_operators()
}
\description{
Deprecated
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_param.R
\name{check_param_names}
\alias{check_param_names}
\title{Check if the \code{param} has the list elements
of a valid \code{param} object.}
\usage{
check_param_names(param)
}
\arguments{
\item{param}{a parameter, as can be created by \code{\link{create_param}}.}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element is missing
}
\seealso{
Use \link{create_param} to create a valid \code{param}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_mean_param}
\alias{create_mean_param}
\alias{create_param_mean}
\title{Create a parameter called mean}
\usage{
create_mean_param(id = NA, value = 0)
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}
}
\value{
a parameter called mean
}
\description{
Create a parameter called mean
}
\note{
this parameter is used in an exponential distribution
  (as returned by \code{\link{create_exp_distr}})
  and normal distribution
  (as returned by \code{\link{create_normal_distr}}).
  It cannot be estimated (as a hyper parameter) yet.
}
\examples{
# Create the parameter
mean_param <- create_mean_param(value = 1.0)

# Use the parameter in a distribution
exp_distr <- create_exp_distr(
  mean = mean_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = exp_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_model_to_xml_operators.R
\name{site_model_to_xml_operators}
\alias{site_model_to_xml_operators}
\title{Converts a site model to XML,
  used in the \code{operators} section}
\usage{
site_model_to_xml_operators(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the site model as XML text
}
\description{
Converts a site model to XML,
  used in the \code{operators} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_xml.R
\name{is_xml}
\alias{is_xml}
\title{Checks if the text is a valid XML node, that is,
it has a opening and matching closing tag}
\usage{
is_xml(text)
}
\arguments{
\item{text}{text to be determined to be valid}
}
\value{
TRUE if the text is valid XML, FALSE otherwise
}
\description{
Checks if the text is a valid XML node, that is,
it has a opening and matching closing tag
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_models_to_xml_tracelog.R
\name{site_models_to_xml_tracelog}
\alias{site_models_to_xml_tracelog}
\title{Creates the site models' XML for the tracelog section}
\usage{
site_models_to_xml_tracelog(site_models)
}
\arguments{
\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}
}
\value{
lines of XML text
}
\description{
Creates the site models' XML for the tracelog section
}
\examples{
# <logger id="tracelog" ...>
#'   # Here
# </logger>
}
\seealso{
the complete tracelog section is created
  by \code{\link{create_tracelog_xml}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_to_xml_run.R
\name{mcmc_to_xml_run_nested_sampling}
\alias{mcmc_to_xml_run_nested_sampling}
\title{Converts an MCMC object to the run section's XML for a Nested-Sampling MCMC}
\usage{
mcmc_to_xml_run_nested_sampling(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
the XML as text
}
\description{
Converts an MCMC object to the run section's XML for a Nested-Sampling MCMC
}
\examples{
#  "<run id=\"mcmc\" spec=\"beast.gss.NS\" chainLength=\"1e+07\" "
#  "particleCount=\"1\" subChainLength=\"5000\" epsilon=\"1e-12\">"
mcmc_to_xml_run_nested_sampling(create_ns_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_site_model.R
\name{is_init_jc69_site_model}
\alias{is_init_jc69_site_model}
\title{Determine if x is an initialized JC69 site model
as created by \code{\link{create_jc69_site_model}}}
\usage{
is_init_jc69_site_model(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized JC69 site model}
}
\value{
TRUE if x is an initialized JC69 site model
}
\description{
Determine if x is an initialized JC69 site model
as created by \code{\link{create_jc69_site_model}}
}
\examples{

jc69_site_model <- create_jc69_site_model(
  gamma_site_model = create_gamma_site_model(
    gamma_cat_count = 2,
    gamma_shape_prior_distr = create_normal_distr()
  )
)
# FALSE: not yet initialized
is_init_jc69_site_model(jc69_site_model)
jc69_site_model <- init_jc69_site_model(jc69_site_model)
# TRUE: now it is initialized
is_init_jc69_site_model(jc69_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tipdate_taxa_to_xml_trait.R
\name{tipdate_taxa_to_xml_trait}
\alias{tipdate_taxa_to_xml_trait}
\title{Internal function}
\usage{
tipdate_taxa_to_xml_trait(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Internal function to creates the '\code{trait}' section
of a BEAST2 XML parameter file,
which is part of a '\code{tree}' section,
without being indented.
}
\details{
The \code{tree} tag has these elements:
\preformatted{
<run[...]>
  <state[...]>
    <tree[...]>
      <trait[...]>
      This part
      </trait>
    </tree>
  </run>
</state>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_inference_models.R
\name{check_inference_models}
\alias{check_inference_models}
\title{Check if the \code{inference_model} is a valid BEAUti inference model.}
\usage{
check_inference_models(inference_models)
}
\arguments{
\item{inference_models}{a list of one or more inference models,
as can be created by \link{create_inference_model}}
}
\value{
nothing
}
\description{
Calls \code{stop} if not.
}
\examples{
check_inference_models(list(create_inference_model()))
}
\seealso{
Use \link{create_inference_model} to create a valid
  BEAST2 options object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_tree_prior.R
\name{is_init_yule_tree_prior}
\alias{is_init_yule_tree_prior}
\title{Determine if x is an initialized Yule tree_prior object}
\usage{
is_init_yule_tree_prior(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized Yule tree prior object}
}
\value{
TRUE if x is an initialized Yule tree_prior object
}
\description{
Determine if x is an initialized Yule tree_prior object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_rate_at_param}
\alias{is_rate_at_param}
\title{Determine if the object is a valid
'rate AT' parameter}
\usage{
is_rate_at_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
'rate AT' parameter}
}
\value{
TRUE if x is a valid 'rate AT' parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
'rate AT' parameter
}
\examples{

is_rate_at_param(create_alpha_param())
is_rate_at_param(create_beta_param())
is_rate_at_param(create_clock_rate_param())
is_rate_at_param(create_kappa_1_param())
is_rate_at_param(create_kappa_2_param())
is_rate_at_param(create_lambda_param())
is_rate_at_param(create_m_param())
is_rate_at_param(create_mean_param())
is_rate_at_param(create_mu_param())
is_rate_at_param(create_rate_ac_param())
is_rate_at_param(create_rate_ag_param())
is_rate_at_param(create_rate_at_param())
is_rate_at_param(create_rate_cg_param())
is_rate_at_param(create_rate_ct_param())
is_rate_at_param(create_rate_gt_param())
is_rate_at_param(create_s_param())
is_rate_at_param(create_scale_param())
is_rate_at_param(create_sigma_param())

is_rate_at_param(NA)
is_rate_at_param(NULL)
is_rate_at_param("nonsense")
is_rate_at_param(create_jc69_site_model())
is_rate_at_param(create_strict_clock_model())
is_rate_at_param(create_yule_tree_prior())
is_rate_at_param(create_mcmc())
}
\seealso{
\code{\link{create_rate_at_param}} creates a 'rate AT' parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beauti_options_v2_6.R
\name{create_beauti_options_v2_6}
\alias{create_beauti_options_v2_6}
\title{Function to create the BEAUti options for version 2.6.}
\usage{
create_beauti_options_v2_6()
}
\value{
a BEAUti options structure
}
\description{
Function to create the BEAUti options for version 2.6, by
calling \link{create_beauti_options}.
}
\examples{
beauti_options <- create_beauti_options_v2_6()
xml <- create_beast2_input(
  get_fasta_filename(),
  beauti_options = beauti_options
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_ids.R
\name{are_ids}
\alias{are_ids}
\title{Determine if x consists out of IDs}
\usage{
are_ids(x)
}
\arguments{
\item{x}{the object to check if it consists out of IDs}
}
\value{
TRUE if x, or all elements of x, are IDs
}
\description{
Determine if x consists out of IDs
}
\seealso{
to check one ID, use \link{is_id}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_distr.R
\name{ccp_tree_prior_to_xml_prior_distr}
\alias{ccp_tree_prior_to_xml_prior_distr}
\title{Creates the tree prior section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a
Coalescent Constant Population tree prior}
\usage{
ccp_tree_prior_to_xml_prior_distr(ccp_tree_prior)
}
\arguments{
\item{ccp_tree_prior}{a Coalescent Constant Population tree prior,
as returned by \code{\link{create_ccp_tree_prior}}}
}
\description{
Creates the tree prior section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a
Coalescent Constant Population tree prior
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tree_prior.R
\name{create_cep_tree_prior}
\alias{create_cep_tree_prior}
\alias{create_tree_prior_cep}
\title{Create a Coalescent Exponential Population tree prior}
\usage{
create_cep_tree_prior(
  id = NA,
  pop_size_distr = create_one_div_x_distr(),
  growth_rate_distr = create_laplace_distr()
)
}
\arguments{
\item{id}{the ID of the alignment}

\item{pop_size_distr}{the population distribution,
as created by a \code{\link{create_distr}} function}

\item{growth_rate_distr}{the growth rate distribution,
as created by a \code{\link{create_distr}} function}
}
\value{
a Coalescent Exponential Population tree_prior
}
\description{
Create a Coalescent Exponential Population tree prior
}
\examples{
cep_tree_prior <- create_cep_tree_prior()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = cep_tree_prior
)
file.remove(beast2_input_file)
}
\seealso{
An alignment ID can be extracted from
  its FASTA filename using \code{\link{get_alignment_id}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tree_priors_n_params.R
\name{get_tree_priors_n_params}
\alias{get_tree_priors_n_params}
\title{Get the number of parameters a list of tree priors has}
\usage{
get_tree_priors_n_params(tree_priors)
}
\arguments{
\item{tree_priors}{one or more tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
the number of parameters the tree priors have
}
\description{
Get the number of parameters a list of tree priors has
}
\examples{
# Two
get_tree_priors_n_params(
  list(
    create_bd_tree_prior(), # zero
    create_cep_tree_prior() # two
  )
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_operator_id_pre.R
\name{get_operator_id_pre}
\alias{get_operator_id_pre}
\title{Get the prefix of operator IDs}
\usage{
get_operator_id_pre(tree_prior)
}
\arguments{
\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
the prefix of operator IDs, similar to the name of a tree prior
}
\description{
Get the prefix of operator IDs
}
\examples{
# BirthDeath
get_operator_id_pre(
  tree_prior = create_bd_tree_prior()
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gamma_site_model_n_distrs.R
\name{get_gamma_site_model_n_distrs}
\alias{get_gamma_site_model_n_distrs}
\title{Get the number of distributions in a gamma site model}
\usage{
get_gamma_site_model_n_distrs(gamma_site_model)
}
\arguments{
\item{gamma_site_model}{a site model's gamma site model,
as returned by \code{\link{create_gamma_site_model}}}
}
\value{
the number of distributions a gamma site model has
}
\description{
Get the number of distributions in a gamma site model
}
\examples{
  gamma_site_model <- create_gamma_site_model()
  n_distrs <- get_gamma_site_model_n_distrs(
    gamma_site_model
   )
  testit::assert(n_distrs == 0)

  gamma_site_model <- create_gamma_site_model(
    gamma_cat_count = 2,
    gamma_shape_prior_distr = create_exp_distr()
  )
  n_distrs <- get_gamma_site_model_n_distrs(gamma_site_model)
  testit::assert(n_distrs == 1)
}
\seealso{
Use \link{create_gamma_site_model} to create a gamma site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_distr.R
\name{create_beast2_input_distr_lh}
\alias{create_beast2_input_distr_lh}
\title{Creates the XML text for the \code{distribution} tag
with the \code{likelihood} ID,
of a BEAST2 parameter file.}
\usage{
create_beast2_input_distr_lh(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
Creates the XML text for the \code{distribution} tag
with the \code{likelihood} ID,
of a BEAST2 parameter file,
in an unindented form
}
\details{
The \code{distribution} tag (with ID equals \code{likelihood})
has these elements:

\preformatted{
  <distribution id="likelihood"[...]>
     <distribution id="treeLikelihood"[...]>
       [...]
     </distribution>
  </distribution>
}

The \code{distribution} section with ID \code{treeLikelihood}
is created by \link{create_tree_likelihood_distr_xml}.

Zooming out:

\preformatted{
  <beast[...]>
    <run[...]>
      <distribution id="posterior"[...]>
        <distribution id="likelihood"[...]>
          [this section]
        </distribution>
      </distribution>
    </run>
  </beast>
}
}
\note{
this function is not intended for regular use, thus its
  long name length is accepted
}
\seealso{
this function is called by \code{create_beast2_input_distr},
  together with \code{create_beast2_input_distr_prior}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_alignment_id.R
\name{check_alignment_id}
\alias{check_alignment_id}
\title{Check if the \code{alignment_id} is valid.}
\usage{
check_alignment_id(alignment_id)
}
\arguments{
\item{alignment_id}{ID of the alignment,
as returned by \link{get_alignment_id}.
Keep at \code{NA} to have it initialized automatically}
}
\description{
Will \link{stop} if not.
}
\examples{
# anthus_aco_sub
created <- get_alignment_id("/home/homer/anthus_aco_sub.fas")
check_alignment_id(created)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/yule_tree_prior_to_xml_operators.R
\name{yule_tree_prior_to_xml_operators}
\alias{yule_tree_prior_to_xml_operators}
\title{Internal function}
\usage{
yule_tree_prior_to_xml_operators(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the tree prior as XML text
}
\description{
Creates the XML of a Yule tree prior,
  as used in the \code{operators} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_poisson_distr}
\alias{is_init_poisson_distr}
\title{Determine if x is an initialized Poisson distribution object
  as created by \code{\link{create_poisson_distr}}}
\usage{
is_init_poisson_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized Poisson distribution object}
}
\value{
TRUE if x is an initialized Poisson distribution object
}
\description{
Determine if x is an initialized Poisson distribution object
  as created by \code{\link{create_poisson_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_tip_dating.R
\name{has_tip_dating}
\alias{has_tip_dating}
\title{Determine if the \code{inference_model} uses tip dating.}
\usage{
has_tip_dating(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
TRUE if the \code{inference_model} uses tip dating,
FALSE otherwise
}
\description{
Determine if the \code{inference_model} uses tip dating
}
\examples{
# Yes, has tip dating
has_strict_clock_model(
  create_inference_model(
    tipdates_filename = get_beautier_path("test_output_0_tipdates.tsv")
  )
)

# No tip dating
has_strict_clock_model(
  create_inference_model()
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tree_prior.R
\name{create_cbs_tree_prior}
\alias{create_cbs_tree_prior}
\alias{create_tree_prior_cbs}
\title{Create a Coalescent Bayesian Skyline tree prior}
\usage{
create_cbs_tree_prior(id = NA, group_sizes_dimension = 5)
}
\arguments{
\item{id}{an alignment's IDs.
An ID can be extracted from its FASTA filename
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{group_sizes_dimension}{the group sizes' dimension,
as used by the CBS tree prior (see \code{\link{create_cbs_tree_prior}})}
}
\value{
a Coalescent Bayesian Skyline tree_prior
}
\description{
Create a Coalescent Bayesian Skyline tree prior
}
\examples{
cbs_tree_prior <- create_cbs_tree_prior()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_beautier_path("test_output_6.fas"),
  beast2_input_file,
  tree_prior = cbs_tree_prior
)
file.remove(beast2_input_file)
}
\seealso{
An alignment ID can be extracted from
  its FASTA filename using \code{\link{get_alignment_id}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa_to_xml_tree.R
\name{taxa_to_xml_tree}
\alias{taxa_to_xml_tree}
\title{Internal function}
\usage{
taxa_to_xml_tree(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Internal function to creates the '\code{tree}' section
of a BEAST2 XML parameter file,
which is part of a '\code{state}' section,
without being indented.
}
\details{
The \code{tree} tag has these elements:
\preformatted{
   <tree[...]>
       <taxonset[...]>
       [...]
       </taxonset>
    </run>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa_to_xml_tree.R
\name{no_taxa_to_xml_tree}
\alias{no_taxa_to_xml_tree}
\title{Internal function}
\usage{
no_taxa_to_xml_tree(inference_model, id = "deprecated")
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}

\item{id}{an alignment's IDs.
An ID can be extracted from its FASTA filename
with \code{\link{get_alignment_ids_from_fasta_filenames}})}
}
\value{
the random phylogeny as XML text
}
\description{
Creates the '\code{tree}' section of a BEAST2 XML parameter file,
which is part of a '\code{state}' section,
without being indented,
when there is no tip-dating
}
\details{
The \code{tree} tag has these elements:
\preformatted{
   <tree[...]>
       <taxonset[...]>
       [...]
       </taxonset>
    </run>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_log_normal_distr}
\alias{is_log_normal_distr}
\title{Determine if the object is a valid
log-normal distribution,
as created by \code{\link{create_log_normal_distr}}}
\usage{
is_log_normal_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
log-normal distribution}
}
\value{
TRUE if x is a valid log-normal distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
log-normal distribution,
as created by \code{\link{create_log_normal_distr}}
}
\examples{
# TRUE
is_log_normal_distr(create_log_normal_distr())
# FALSE
is_log_normal_distr(create_normal_distr())
is_distr(NA)
is_distr(NULL)
is_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unindent.R
\name{unindent}
\alias{unindent}
\title{Unindents text}
\usage{
unindent(text)
}
\arguments{
\item{text}{one or more lines of text}
}
\value{
unindented lines of text
}
\description{
Unindents text
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_file_from_model.R
\name{create_beast2_input_file_from_model}
\alias{create_beast2_input_file_from_model}
\title{Create a BEAST2 input file from an inference model}
\usage{
create_beast2_input_file_from_model(
  input_filename,
  output_filename,
  inference_model = beautier::create_inference_model()
)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{output_filename}{Name of the XML parameter file created by this
function. BEAST2 uses this file as input.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
nothing
}
\description{
Create a BEAST2 input file from an inference model
}
\examples{
output_filename <- get_beautier_tempfilename()
create_beast2_input_file_from_model(
  input_filename = get_fasta_filename(),
  output_filename = output_filename,
  inference_model = create_inference_model()
)
file.remove(output_filename)
}
\seealso{
use \link{create_beast2_input_from_model} to
get the BEAST2 input file as text

See \code{\link{create_site_model}} for examples with
different site models.
See \code{\link{create_clock_model}} for examples
with clock models.
See \code{\link{create_tree_prior}} for examples with
different tree priors.
See \code{\link{create_mcmc}} for examples with
a different MCMC setup.
Use \link{create_beast2_input_file} to do the same with the elements
of an inference model.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_mrca_align_id_in_fasta.R
\name{is_mrca_align_id_in_fasta}
\alias{is_mrca_align_id_in_fasta}
\title{Determine if an MRCA prior's alignment IDs is present in the FASTA file}
\usage{
is_mrca_align_id_in_fasta(mrca_prior, fasta_filename)
}
\arguments{
\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{fasta_filename}{a FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.
Note that BEAST2 also supports missing data,
by using a dash (\code{-}) or question mark (\code{?})
as a sequence.}
}
\value{
TRUE if the MRCA prior's alignment IDs
  is present in the FASTA file.
  Returns FALSE otherwise
}
\description{
Determine if an MRCA prior's alignment IDs is present in the FASTA file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equal_treelogs.R
\name{are_equal_treelogs}
\alias{are_equal_treelogs}
\title{Determine if two treelogs are equal.}
\usage{
are_equal_treelogs(treelog_1, treelog_2)
}
\arguments{
\item{treelog_1}{an treelog, as created by \link{create_treelog}}

\item{treelog_2}{an treelog, as created by \link{create_treelog}}
}
\value{
TRUE if the two treelogs are equal
}
\description{
Will \link{stop} if the arguments are not treelogs.
}
\examples{

treelog_1 <- create_treelog(log_every = 1000)
treelog_2 <- create_treelog(log_every = 314)
# TRUE
are_equal_treelogs(treelog_1, treelog_1)
# FALSE
are_equal_treelogs(treelog_1, treelog_2)
}
\seealso{
Use \link{create_treelog} to create an treelog
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_clock_model_from_name.R
\name{create_clock_model_from_name}
\alias{create_clock_model_from_name}
\title{Create a clock model from name}
\usage{
create_clock_model_from_name(clock_model_name)
}
\arguments{
\item{clock_model_name}{name of a clock model,
must be a name as returned by \code{\link{get_clock_model_names}}}
}
\value{
a clock model,
  as can be created by \link{create_clock_model}
}
\description{
Create a clock model from name
}
\examples{
create_clock_model_from_name(get_clock_model_names()[1])
}
\seealso{
Use \code{\link{create_clock_model}} to create a clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_poisson}
\alias{distr_to_xml_poisson}
\title{Internal function}
\usage{
distr_to_xml_poisson(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{a Poisson distribution,
as created by \code{\link{create_poisson_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a Poisson distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_subst_model_xml.R
\name{create_tn93_subst_model_xml}
\alias{create_tn93_subst_model_xml}
\title{Converts a TN93 site model to XML,
  used in the \code{substModel} section}
\usage{
create_tn93_subst_model_xml(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the site model as XML text
}
\description{
Converts a TN93 site model to XML,
  used in the \code{substModel} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_one_div_x_distr}
\alias{is_init_one_div_x_distr}
\title{Determine if x is an initialized one_div_x distribution object
  as created by \code{\link{create_one_div_x_distr}}}
\usage{
is_init_one_div_x_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized one_div_x distribution object}
}
\value{
TRUE if x is an initialized one_div_x distribution object
}
\description{
Determine if x is an initialized one_div_x distribution object
  as created by \code{\link{create_one_div_x_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_clock_models.R
\name{init_rln_clock_model}
\alias{init_rln_clock_model}
\title{Initializes a Relaxed Log-Normal clock model}
\usage{
init_rln_clock_model(rln_clock_model, distr_id = 0, param_id = 0)
}
\arguments{
\item{rln_clock_model}{a Relaxed Log-Normal clock model,
as returned by \code{\link{create_rln_clock_model}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized Relaxed Log-Normal clock model
}
\description{
Initializes a Relaxed Log-Normal clock model
}
\examples{

rln_clock_model <- create_rln_clock_model()
# FALSE: not yet initialized
is_init_rln_clock_model(rln_clock_model)
rln_clock_model <- init_rln_clock_model(rln_clock_model)
# Dimension is set to NA by default, for unknown reasons.
# Because 'init_rln_clock_model' does not initialize it (for
# unknown reasons), set it manually
rln_clock_model$dimension <- 42
# TRUE: now it is initialized
is_init_rln_clock_model(rln_clock_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tree_prior_names.R
\name{get_tree_prior_names}
\alias{get_tree_prior_names}
\title{Get the tree prior names}
\usage{
get_tree_prior_names()
}
\value{
the tree prior names
}
\description{
Get the tree prior names
}
\examples{
get_tree_prior_names()
}
\seealso{
Use \link{create_tree_priors} to get all tree priors
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_clock_models_from_names.R
\name{create_clock_models_from_names}
\alias{create_clock_models_from_names}
\title{Create clock models from their names}
\usage{
create_clock_models_from_names(clock_model_names)
}
\arguments{
\item{clock_model_names}{one or more names of a clock model,
must be name among those returned by \code{\link{get_clock_model_names}}}
}
\value{
a list of one or more clock models
}
\description{
Create clock models from their names
}
\examples{
clock_models <- create_clock_models_from_names(get_clock_model_names())
}
\seealso{
Use \link{create_clock_models} to get all clock models
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_subst_model_xml.R
\name{create_gtr_subst_model_xml}
\alias{create_gtr_subst_model_xml}
\title{Converts a GTR site model to XML,
  used in the \code{substModel} section}
\usage{
create_gtr_subst_model_xml(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the site model as XML text
}
\description{
Converts a GTR site model to XML,
  used in the \code{substModel} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_mu}
\alias{parameter_to_xml_mu}
\title{Internal function}
\usage{
parameter_to_xml_mu(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a mu parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_mu_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a mu parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equal_xml_files.R
\name{are_equal_xml_files}
\alias{are_equal_xml_files}
\title{Determine if XML files result in equal trees}
\usage{
are_equal_xml_files(filename_1, filename_2, section)
}
\arguments{
\item{filename_1}{name of a first XML file}

\item{filename_2}{name of a second XML file}

\item{section}{name of an XML section.
Assumes that there is one line that starts with \code{<section}
(excluding whitespace)
and one line that is \code{</section>} (also excluding whitespace)}
}
\value{
TRUE if the two sections of the XML files are equal,
  FALSE otherwise
}
\description{
Determine if XML files result in equal trees
}
\seealso{
to check for equivalence, use \code{\link{are_equivalent_xml_files}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_lambda_param}
\alias{create_lambda_param}
\alias{create_param_lambda}
\title{Create a parameter called lambda}
\usage{
create_lambda_param(id = NA, value = 0)
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}
}
\value{
a parameter called lambda
}
\description{
Create a parameter called lambda
}
\note{
this parameter is used in a Poisson distribution
  (as returned by \code{\link{create_poisson_distr}})
}
\examples{
# Create the parameter
lambda_param <- create_lambda_param()

# Use the parameter in a distribution
poisson_distr <- create_poisson_distr(
  lambda = lambda_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = poisson_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_model_to_xml_tracelog.R
\name{site_model_to_xml_tracelog}
\alias{site_model_to_xml_tracelog}
\title{Creates the site model's XML for the tracelog section}
\usage{
site_model_to_xml_tracelog(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
lines of XML text
}
\description{
Creates the site model's XML for the tracelog section
}
\examples{
# <logger id="tracelog" ...>
#'   # Here
# </logger>
}
\seealso{
all site models' tracelog section is created
  by \code{\link{site_models_to_xml_tracelog}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_site_model.R
\name{is_init_site_model}
\alias{is_init_site_model}
\title{Determine if x is an initialized site model,
as created by \code{\link{create_site_model}}}
\usage{
is_init_site_model(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized site_models object}
}
\value{
TRUE if x is an initialized site model
}
\description{
Determine if x is an initialized site model,
as created by \code{\link{create_site_model}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_xml_opening_tag.R
\name{get_xml_opening_tag}
\alias{get_xml_opening_tag}
\title{Get the XML opening tag}
\usage{
get_xml_opening_tag(text)
}
\arguments{
\item{text}{text to be determined to be valid}
}
\value{
the opening tag if found, else NA
}
\description{
Get the XML opening tag
}
\examples{
# my_tag
get_xml_opening_tag("<my_tag text=something/>")

# NA when there is no opening tag
get_xml_opening_tag("no_xml")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_beta_param}
\alias{create_beta_param}
\alias{create_param_beta}
\title{Create a parameter called beta}
\usage{
create_beta_param(id = NA, value = 1)
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}
}
\value{
a parameter called beta
}
\description{
Create a parameter called beta
}
\note{
this parameter is used in a beta distribution
  (as returned by \code{\link{create_beta_distr}})
and gamma distribution
  (as returned by \code{\link{create_gamma_distr}})
and inverse-gamma distribution
  (as returned by \code{\link{create_inv_gamma_distr}}).
It cannot be estimated (as a hyper parameter) yet.
}
\examples{
# Create the parameter
beta_param <- create_beta_param()

# Use the parameter in a distribution
gamma_distr <- create_gamma_distr(
  beta = beta_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = gamma_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrca_priors_to_xml_tracelog.R
\name{mrca_priors_to_xml_tracelog}
\alias{mrca_priors_to_xml_tracelog}
\title{Creates the MRCA priors' XML for the tracelog section}
\usage{
mrca_priors_to_xml_tracelog(clock_models, mrca_priors, tipdates_filename = NA)
}
\arguments{
\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}

\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
lines of XML text
}
\description{
Creates the MRCA priors' XML for the tracelog section.
}
\details{
\code{
  <logger id="tracelog" ...>
   # Here
  </logger>
}
}
\seealso{
the complete tracelog section is created
  by \code{\link{create_tracelog_xml}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tracelog_xml.R
\name{create_tracelog_xml}
\alias{create_tracelog_xml}
\title{Internal function}
\usage{
create_tracelog_xml(input_filename, inference_model)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
Creates the \code{tracelog} section of the \code{logger} section
of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equal_screenlogs.R
\name{are_equal_screenlogs}
\alias{are_equal_screenlogs}
\title{Determine if two screenlogs are equal.}
\usage{
are_equal_screenlogs(screenlog_1, screenlog_2)
}
\arguments{
\item{screenlog_1}{an screenlog, as created by \link{create_screenlog}}

\item{screenlog_2}{an screenlog, as created by \link{create_screenlog}}
}
\value{
TRUE if the two screenlogs are equal
}
\description{
Will \link{stop} if the arguments are not screenlogs.
}
\examples{
screenlog_1 <- create_screenlog(log_every = 1000)
screenlog_2 <- create_screenlog(log_every = 314)
# TRUE
are_equal_screenlogs(screenlog_1, screenlog_1)
# FALSE
are_equal_screenlogs(screenlog_1, screenlog_2)
}
\seealso{
Use \link{create_screenlog} to create an screenlog
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_one_int.R
\name{is_one_int}
\alias{is_one_int}
\title{Determines if the argument is a whole number}
\usage{
is_one_int(x, tolerance = .Machine$double.eps^0.5)
}
\arguments{
\item{x}{the object to be determined of if it is one integer}

\item{tolerance}{tolerance to rounding errors}
}
\description{
Determines if the argument is a whole number
}
\examples{
# TRUE
is_one_int(314)
is_one_int(0)
is_one_int(-314)
# FALSE
is_one_int(3.14)
is_one_int(NULL)
is_one_int(NA)
is_one_int(Inf)
is_one_int("nonsense")
is_one_int(c())
is_one_int(c(1, 2))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_prior_to_xml_state.R
\name{tree_prior_to_xml_state}
\alias{tree_prior_to_xml_state}
\title{Creates the XML of a tree prior,
  as used in the \code{state} section}
\usage{
tree_prior_to_xml_state(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the tree prior as XML text
}
\description{
Creates the XML of a tree prior,
  as used in the \code{state} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_site_models.R
\name{init_tn93_site_model}
\alias{init_tn93_site_model}
\title{Initializes a TN93 site model}
\usage{
init_tn93_site_model(tn93_site_model, distr_id = 0, param_id = 0)
}
\arguments{
\item{tn93_site_model}{a TN93 site model,
as returned by \code{\link{create_tn93_site_model}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized TN93 site model
}
\description{
Initializes a TN93 site model
}
\examples{

tn93_site_model <- create_tn93_site_model()
is_init_tn93_site_model(tn93_site_model)
tn93_site_model <- init_tn93_site_model(tn93_site_model)
is_init_tn93_site_model(tn93_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_subst_model_xml.R
\name{create_jc69_subst_model_xml}
\alias{create_jc69_subst_model_xml}
\title{Converts a JC69 site model to XML,
used in the \code{substModel} section}
\usage{
create_jc69_subst_model_xml(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the site model as XML text
}
\description{
Converts a JC69 site model to XML,
used in the \code{substModel} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_store_every.R
\name{check_store_every}
\alias{check_store_every}
\title{Check if \code{store_every} holds a valid value}
\usage{
check_store_every(store_every)
}
\arguments{
\item{store_every}{number of states the MCMC will process
before the posterior's state will be saved to file.
Use -1 or \code{NA} to use the default frequency.}
}
\description{
Will \link{stop} if not
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_tree_prior_name.R
\name{is_tree_prior_name}
\alias{is_tree_prior_name}
\title{Determines if the name is a valid tree prior name}
\usage{
is_tree_prior_name(name)
}
\arguments{
\item{name}{the name to be tested}
}
\value{
TRUE if the name is a valid tree_prior name, FALSE otherwise
}
\description{
Determines if the name is a valid tree prior name
}
\examples{
# TRUE
is_tree_prior_name("birth_death")
is_tree_prior_name("coalescent_bayesian_skyline")
is_tree_prior_name("coalescent_constant_population")
is_tree_prior_name("coalescent_exp_population")
is_tree_prior_name("yule")
# FALSE
is_tree_prior_name("nonsense")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_beauti_options.R
\name{is_beauti_options}
\alias{is_beauti_options}
\title{Determine if the object is a valid \code{beauti_options}}
\usage{
is_beauti_options(x)
}
\arguments{
\item{x}{an object, to be determined if it is a \code{beauti_options}}
}
\value{
\link{TRUE} if the object is a valid \code{beauti_options},
  \link{FALSE} otherwise
}
\description{
Determine if the object is a valid \code{beauti_options}
}
\examples{
# TRUE
is_beauti_options(create_beauti_options())

# FALSE
is_beauti_options("nonsense")
is_beauti_options(NA)
is_beauti_options(NULL)
is_beauti_options("")
is_beauti_options(c())
}
\seealso{
use \link{create_beauti_options} to create a valid
\code{beauti_options} object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_gamma_site_model.R
\name{create_gamma_site_model}
\alias{create_gamma_site_model}
\title{Create a gamma site model, part of a site model}
\usage{
create_gamma_site_model(
  gamma_cat_count = "0",
  gamma_shape = "1.0",
  prop_invariant = "0.0",
  gamma_shape_prior_distr = NA,
  freq_equilibrium = "estimated"
)
}
\arguments{
\item{gamma_cat_count}{the number of gamma categories, must
be an integer with value zero or more}

\item{gamma_shape}{gamma curve shape parameter}

\item{prop_invariant}{the proportion invariant, must be a value
from 0.0 to 1.0}

\item{gamma_shape_prior_distr}{the distribution of the gamma shape prior.
\code{gamma_shape_prior_distr} must be \code{NA} for
a \code{gamma_cat_count} of zero or one.
For a \code{gamma_cat_count} of two or more,
leaving \code{gamma_shape_prior_distr} equal to its default
value of \code{NA}, a default distribution is used.
Else \code{gamma_shape_prior_distr} must be a
distribution, as can be created by \code{\link{create_distr}}}

\item{freq_equilibrium}{the frequency in which the rates are at equilibrium
are either \code{estimated}, \code{empirical} or \code{all_equal}.
\code{get_freq_equilibrium_names} returns the possible values
for \code{freq_equilibrium}}
}
\value{
a gamma site model
}
\description{
Create a gamma site model, part of a site model
}
\examples{
gamma_site_model <- create_gamma_site_model(prop_invariant = 0.5)

site_model <- create_hky_site_model(gamma_site_model = gamma_site_model)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  site_model = site_model
)
file.remove(beast2_input_file)
}
\seealso{
Use \code{\link{create_gamma_site_model}}
  to create a gamma site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_inv_gamma_distr}
\alias{is_init_inv_gamma_distr}
\title{Determine if x is an initialized inverse-gamma distribution
  as created by \code{\link{create_inv_gamma_distr}}}
\usage{
is_init_inv_gamma_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized inverse-gamma distribution}
}
\value{
TRUE if x is an initialized inverse-gamma distribution
}
\description{
Determine if x is an initialized inverse-gamma distribution
  as created by \code{\link{create_inv_gamma_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_model.R
\name{create_hky_site_model}
\alias{create_hky_site_model}
\alias{create_site_model_hky}
\title{Create an HKY site model}
\usage{
create_hky_site_model(
  id = NA,
  kappa = "2.0",
  gamma_site_model = create_gamma_site_model(),
  kappa_prior_distr = create_log_normal_distr(m = create_m_param(value = "1.0"), s =
    1.25),
  freq_equilibrium = "estimated"
)
}
\arguments{
\item{id}{the IDs of the alignment (can be extracted from
the FASTA filename using \code{\link{get_alignment_id}})}

\item{kappa}{the kappa}

\item{gamma_site_model}{a gamma site model, as created
by \code{\link{create_gamma_site_model}}}

\item{kappa_prior_distr}{the distribution of the kappa prior,
which is a log-normal distribution
(as created by \code{\link{create_log_normal_distr}})
by default}

\item{freq_equilibrium}{the frequency in which the rates are at equilibrium
are either \code{estimated}, \code{empirical} or \code{all_equal}.
\code{get_freq_equilibrium_names} returns the possible values
for \code{freq_equilibrium}}
}
\value{
an HKY site_model
}
\description{
Create an HKY site model
}
\examples{
 hky_site_model <- create_hky_site_model()
 output_filename <- get_beautier_tempfilename()
 create_beast2_input_file(
   input_filename = get_fasta_filename(),
   output_filename = output_filename,
   site_model = hky_site_model
 )
file.remove(output_filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_rln_clock_model.R
\name{has_rln_clock_model}
\alias{has_rln_clock_model}
\title{Determine if the \code{inference_model} uses
a relaxed log-normal clock model.}
\usage{
has_rln_clock_model(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
TRUE if the \code{inference_model} uses
a relaxed log-normal clock model,
FALSE otherwise
}
\description{
Determine if the \code{inference_model} uses
a relaxed log-normal clock model.
}
\examples{
# Yes, has a RLN clock model
has_rln_clock_model(
  create_inference_model(clock_model = create_rln_clock_model())
)

# No RLN clock model
has_rln_clock_model(
  create_inference_model(clock_model = create_strict_clock_model())
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_test_tracelog.R
\name{create_test_tracelog}
\alias{create_test_tracelog}
\title{Create a \code{tracelog} object}
\usage{
create_test_tracelog(
  filename = create_temp_tracelog_filename(),
  log_every = 1000,
  mode = "autodetect",
  sanitise_headers = TRUE,
  sort = "smart"
)
}
\arguments{
\item{filename}{name of the file to store the posterior traces.
Use \link{NA} to use the filename \code{[alignment_id].log},
where \code{alignment_id} is obtained using \link{get_alignment_id}}

\item{log_every}{number of MCMC states between writing to file}

\item{mode}{mode how to log.
Valid values are the ones returned by \link{get_log_modes}}

\item{sanitise_headers}{set to \link{TRUE} to sanitise the headers of the
log file}

\item{sort}{how to sort the log.
Valid values are the ones returned by \link{get_log_sorts}}
}
\description{
Create a \code{tracelog} object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_tree_prior.R
\name{is_init_ccp_tree_prior}
\alias{is_init_ccp_tree_prior}
\title{Determine if x is an initialized Coalescent Constant Population
  tree_prior object}
\usage{
is_init_ccp_tree_prior(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized Coalescent Constant Population tree prior object}
}
\value{
TRUE if x is an initialized Coalescent Constant Population
  tree prior object
}
\description{
Determine if x is an initialized Coalescent Constant Population
  tree_prior object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/freq_equilibrium_to_xml.R
\name{freq_equilibrium_to_xml}
\alias{freq_equilibrium_to_xml}
\title{Creates the \code{freq_equilibrium} as XML}
\usage{
freq_equilibrium_to_xml(freq_equilibrium, id)
}
\arguments{
\item{freq_equilibrium}{a \code{freq_equilibrium} name}

\item{id}{a site model's name}
}
\value{
the \code{freq_equilibrium} as XML
}
\description{
Creates the \code{freq_equilibrium} as XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_rate_param_to_xml.R
\name{clock_rate_param_to_xml}
\alias{clock_rate_param_to_xml}
\title{Internal function}
\usage{
clock_rate_param_to_xml(
  clock_rate_param,
  beauti_options = create_beauti_options()
)
}
\arguments{
\item{clock_rate_param}{a \code{clockRate} parameter,
a numeric value,
as created by \link{create_clock_rate_param}}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a \code{clockRate} parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equivalent_xmls.R
\name{are_equivalent_xml_lines}
\alias{are_equivalent_xml_lines}
\title{Determine if XML lines result in equivalent trees}
\usage{
are_equivalent_xml_lines(lines_1, lines_2, section = NA, verbose = FALSE)
}
\arguments{
\item{lines_1}{lines of a first XML file}

\item{lines_2}{lines of a second XML file}

\item{section}{the name of the XML section}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
TRUE if the two XML lines result in equivalent trees,
  FALSE otherwise
}
\description{
Determine if XML lines result in equivalent trees
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_mrca_priors.R
\name{init_mrca_priors}
\alias{init_mrca_priors}
\title{Initializes all MRCA priors}
\usage{
init_mrca_priors(mrca_priors, distr_id = 0, param_id = 0)
}
\arguments{
\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}

\item{distr_id}{the first distributions' ID}

\item{param_id}{the first parameter's ID}
}
\value{
a list of initialized MRCA priors
}
\description{
Initializes all MRCA priors
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_mcmc.R
\name{is_mcmc}
\alias{is_mcmc}
\title{Determine if the object is a valid MCMC}
\usage{
is_mcmc(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid MCMC}
}
\value{
TRUE if x is a valid MCMC, FALSE otherwise
}
\description{
Determine if the object is a valid MCMC
}
\examples{
# Returns TRUE
is_mcmc(create_mcmc())
is_mcmc(create_ns_mcmc())

# Returns FALSE
is_mcmc("nonsense")
is_mcmc(NULL)
is_mcmc(NA)
is_mcmc("")
is_mcmc(c())
}
\seealso{
Use \code{\link{create_mcmc}} to create an MCMC
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alignment_id.R
\name{get_alignment_id}
\alias{get_alignment_id}
\title{Conclude the ID from a FASTA filename.}
\usage{
get_alignment_id(fasta_filename, capitalize_first_char_id = FALSE)
}
\arguments{
\item{fasta_filename}{a FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.
Note that BEAST2 also supports missing data,
by using a dash (\code{-}) or question mark (\code{?})
as a sequence.}

\item{capitalize_first_char_id}{if TRUE, the first character will
be capitalized}
}
\value{
an alignment's ID
}
\description{
This is done in the same way as BEAST2 will do so.
}
\examples{
# Path need not exist, use UNIX path as example
# anthus_aco_sub
created <- get_alignment_id("/home/homer/anthus_aco_sub.fas")
check_alignment_id(created)
}
\seealso{
Use \link{check_alignment_id} to check if an alignment
ID is valid.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_scale_param}
\alias{create_scale_param}
\alias{create_param_scale}
\title{Create a parameter called scale}
\usage{
create_scale_param(id = NA, value = 0)
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}
}
\value{
a parameter called scale
}
\description{
Create a parameter called scale
}
\note{
this parameter is used in a Laplace distribution
  (as returned by \code{\link{create_laplace_distr}})
}
\examples{
# Create the parameter
scale_param <- create_scale_param()

# Use the parameter in a distribution
laplace_distr <- create_laplace_distr(
  scale = scale_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = laplace_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rnd_phylo_xml_init.R
\name{rnd_phylo_to_xml_init}
\alias{rnd_phylo_to_xml_init}
\title{Creates the XML of a random phylogeny,
  as used in the \code{init} section}
\usage{
rnd_phylo_to_xml_init(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the phylogeny as XML text
}
\description{
Creates the XML text for the \code{beast} tag of a BEAST2 parameter file,
which is directly after the XML
declaration (created by \link{create_xml_declaration}.
}
\details{
The \code{init} tag has these elements:

\preformatted{
  <init id=\"RandomTree.t:[...]>
      <populationModel[...]>
      [...]
      </populationModel>
  </init>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrca_prior_to_xml_state.R
\name{mrca_prior_to_xml_state}
\alias{mrca_prior_to_xml_state}
\title{Internal function to create the XML of an MRCA prior,
  as used in the \code{state} section}
\usage{
mrca_prior_to_xml_state(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the tree prior as XML text
}
\description{
Internal function to create the XML of an MRCA prior,
  as used in the \code{state} section
}
\examples{
mrca_prior_to_xml_state(
  inference_model = create_inference_model(
    mrca_prior = create_mrca_prior(
      alignment_id = "test_output_0",
      mrca_distr = create_normal_distr(id = 42)
    ),
    clock_model = create_strict_clock_model()
  )
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_beast_xml.R
\name{create_beast2_beast_xml}
\alias{create_beast2_beast_xml}
\title{Create the \code{<beast ...>} XML}
\usage{
create_beast2_beast_xml(beauti_options)
}
\arguments{
\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the XML
}
\description{
The \code{<beast ...>} XML is the XML at the start of a BEAST2
XML input file, directly after the general XML declaration (as
created by \link{create_xml_declaration}).
}
\examples{
create_beast2_beast_xml(
  beauti_options = create_beauti_options_v2_6()
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_rename_fun.R
\name{check_rename_fun}
\alias{check_rename_fun}
\title{Check if the rename function is a valid filename rename function}
\usage{
check_rename_fun(rename_fun)
}
\arguments{
\item{rename_fun}{a function to rename a filename,
as can be checked by \link{check_rename_fun}. This function should
have one argument, which will be a filename or \link{NA}. The
function should \link{return} one filename (when passed one filename) or
one \link{NA} (when passed one \link{NA}).
Example rename functions are:
\itemize{
  \item \link{get_remove_dir_fun} get a function that removes the directory
    paths from the filenames, in effect turning these into local files
  \item \link{get_replace_dir_fun} get a function that replaces the directory
    paths from the filenames
  \item \link{get_remove_hex_fun} get a function that removes the
    hex string from filenames.
    For example, \code{tracelog_82c1a522040.log} becomes \code{tracelog.log}
}}
}
\description{
Will \link{stop} if not
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_clock_model.R
\name{is_init_strict_clock_model}
\alias{is_init_strict_clock_model}
\title{Determine if x is an initialized strict clock_model object}
\usage{
is_init_strict_clock_model(strict_clock_model)
}
\arguments{
\item{strict_clock_model}{a strict clock model,
as returned by \code{\link{create_strict_clock_model}}}
}
\value{
TRUE if x is an initialized strict clock_model object
}
\description{
Determine if x is an initialized strict clock_model object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_uniform_distr}
\alias{is_init_uniform_distr}
\title{Determine if x is an initialized uniform distribution object
  as created by \code{\link{create_uniform_distr}}}
\usage{
is_init_uniform_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized uniform distribution object}
}
\value{
TRUE if x is an initialized uniform distribution object
}
\description{
Determine if x is an initialized uniform distribution object
  as created by \code{\link{create_uniform_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_inference_model_filenames.R
\name{get_inference_model_filenames}
\alias{get_inference_model_filenames}
\title{Get the filenames stored in an inference model.}
\usage{
get_inference_model_filenames(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
If there is no name for a \code{tipdates} file specified (as done by
setting \code{inference_model$tipdates_filename} to \link{NA},
there will be one filename less returned
}
\examples{
inference_model <- create_inference_model()
filenames <- get_inference_model_filenames(inference_model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_tn93_site_model.R
\name{check_tn93_site_model}
\alias{check_tn93_site_model}
\title{Check if the \code{tn93_site_model} is a valid
TN93 nucleotide substitution model.}
\usage{
check_tn93_site_model(tn93_site_model)
}
\arguments{
\item{tn93_site_model}{a TN93 site model,
as returned by \code{\link{create_tn93_site_model}}}
}
\description{
Use \link{create_tn93_site_model} to create a valid
TN93 nucleotide substitution model.
}
\examples{
check_tn93_site_model(create_tn93_site_model())
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_gamma_site_model.R
\name{check_gamma_site_model_names}
\alias{check_gamma_site_model_names}
\title{Checks if the gamma site model has the right list elements' names}
\usage{
check_gamma_site_model_names(gamma_site_model)
}
\arguments{
\item{gamma_site_model}{a site model's gamma site model,
as returned by \code{\link{create_gamma_site_model}}}
}
\value{
nothing. Will call \code{stop} if the argument is not a valid
  gamma site model
}
\description{
Checks if the gamma site model has the right list elements' names
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tree_prior_n_params.R
\name{get_tree_prior_n_params}
\alias{get_tree_prior_n_params}
\title{Get the number of parameters a tree prior has}
\usage{
get_tree_prior_n_params(tree_prior)
}
\arguments{
\item{tree_prior}{a tree_prior,
as created by \code{\link{create_tree_prior}}}
}
\value{
the number of parameters a tree prior has
}
\description{
Get the number of parameters a tree prior has
}
\examples{
 # birth_rate_distr is uniform, which has zero parameters
 # death_rate_distr is uniform, which has zero parameters
 testit::assert(
   get_tree_prior_n_params(create_bd_tree_prior()) == 0
 )

 # no distributions, no parameters
 testit::assert(
   get_tree_prior_n_params(create_cbs_tree_prior()) == 0
 )

 # pop_size_distr is 1/x, which has zero parameters
 testit::assert(
   get_tree_prior_n_params(create_ccp_tree_prior()) == 0
 )

 # pop_size_distr is 1/x, which has zero parameters
 # growth_rate_distr is Laplace, which has two parameters
 testit::assert(
   get_tree_prior_n_params(create_cep_tree_prior()) == 2
 )

 # birth_rate_distr is uniform, which has zero parameters
 testit::assert(
   get_tree_prior_n_params(create_yule_tree_prior()) == 0
 )
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_clock_model.R
\name{is_init_clock_model}
\alias{is_init_clock_model}
\title{Determine if x is an initialized clock_model object,
as created by \code{\link{create_clock_model}}}
\usage{
is_init_clock_model(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized clock_models object}
}
\value{
TRUE if x is an initialized clock_model object
}
\description{
Determine if x is an initialized clock_model object,
as created by \code{\link{create_clock_model}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_laplace}
\alias{distr_to_xml_laplace}
\title{Internal function}
\usage{
distr_to_xml_laplace(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{a Laplace distribution
as created by \code{\link{create_laplace_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a Laplace distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_clock_rate_param}
\alias{is_clock_rate_param}
\title{Determine if the object is a valid
clock_rate parameter}
\usage{
is_clock_rate_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
clock_rate parameter}
}
\value{
TRUE if x is a valid clock_rate parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
clock_rate parameter
}
\examples{

is_clock_rate_param(create_alpha_param())
is_clock_rate_param(create_beta_param())
is_clock_rate_param(create_clock_rate_param())
is_clock_rate_param(create_kappa_1_param())
is_clock_rate_param(create_kappa_2_param())
is_clock_rate_param(create_lambda_param())
is_clock_rate_param(create_m_param())
is_clock_rate_param(create_mean_param())
is_clock_rate_param(create_mu_param())
is_clock_rate_param(create_rate_ac_param())
is_clock_rate_param(create_rate_ag_param())
is_clock_rate_param(create_rate_at_param())
is_clock_rate_param(create_rate_cg_param())
is_clock_rate_param(create_rate_ct_param())
is_clock_rate_param(create_rate_gt_param())
is_clock_rate_param(create_s_param())
is_clock_rate_param(create_scale_param())
is_clock_rate_param(create_sigma_param())

is_clock_rate_param(NA)
is_clock_rate_param(NULL)
is_clock_rate_param("nonsense")
is_clock_rate_param(create_jc69_site_model())
is_clock_rate_param(create_strict_clock_model())
is_clock_rate_param(create_yule_tree_prior())
is_clock_rate_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_treelog.R
\name{check_treelog_names}
\alias{check_treelog_names}
\title{Check if the \code{treelog} has the list elements
of a valid \code{treelog} object.}
\usage{
check_treelog_names(treelog)
}
\arguments{
\item{treelog}{a \code{treelog},
as created by \link{create_treelog}}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element is missing
}
\seealso{
Use \link{create_treelog} to create a valid \code{treelog}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_rate_ac_param}
\alias{create_rate_ac_param}
\alias{create_param_rate_ac}
\title{Create a parameter called 'rate AC'}
\usage{
create_rate_ac_param(id = NA, estimate = TRUE, value = "1.0", lower = "0.0")
}
\arguments{
\item{id}{the parameter's ID}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}

\item{value}{value of the parameter}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}
}
\value{
a parameter called 'rate AC'
}
\description{
Create a parameter called 'rate AC'
}
\examples{
# Create parameter
rate_ac_param <- create_rate_ac_param(value = 1, estimate = FALSE)

# Use the parameter to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  site_model = create_gtr_site_model(
    rate_ac_param = rate_ac_param
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxa_names.R
\name{get_taxa_names}
\alias{get_taxa_names}
\title{Extract the names of taxa from a file}
\usage{
get_taxa_names(filename)
}
\arguments{
\item{filename}{name of a FASTA file}
}
\value{
the taxa names
}
\description{
Extract the names of taxa from a file
}
\examples{
  created <- get_taxa_names(get_beautier_path("anthus_aco_sub.fas"))
  expected <- c(
    "61430_aco", "626029_aco", "630116_aco", "630210_aco", "B25702_aco"
   )
  testit::assert(created == expected)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_param.R
\name{check_param_types}
\alias{check_param_types}
\title{Check if the \code{param} has the list elements
of the right type for a valid \code{param} object.}
\usage{
check_param_types(param)
}
\arguments{
\item{param}{a parameter, as can be created by \code{\link{create_param}}.}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element has the incorrect type
}
\seealso{
Use \link{create_param} to create a valid \code{param}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_init_clock_models.R
\name{are_init_clock_models}
\alias{are_init_clock_models}
\title{Determine if x consists out of initialized clock_models objects}
\usage{
are_init_clock_models(x)
}
\arguments{
\item{x}{the object to check if it consists out of
initialized clock_models objects}
}
\value{
TRUE if x, or all elements of x, are initialized clock_model objects
}
\description{
Determine if x consists out of initialized clock_models objects
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_site_models.R
\name{init_jc69_site_model}
\alias{init_jc69_site_model}
\title{Initializes a JC69 site model}
\usage{
init_jc69_site_model(jc69_site_model, distr_id = 0, param_id = 0)
}
\arguments{
\item{jc69_site_model}{a JC69 site model,
as returned by \code{\link{create_jc69_site_model}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized HKY site model
}
\description{
Initializes a JC69 site model
}
\examples{

hky_site_model <- create_hky_site_model()
is_init_hky_site_model(hky_site_model)
hky_site_model <- init_hky_site_model(hky_site_model)
is_init_hky_site_model(hky_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_m_param}
\alias{create_m_param}
\alias{create_param_m}
\title{Create a parameter called m}
\usage{
create_m_param(id = NA, estimate = FALSE, lower = NA, upper = NA, value = 0)
}
\arguments{
\item{id}{the parameter's ID}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}

\item{upper}{upper value of the parameter}

\item{value}{value of the parameter}
}
\value{
a parameter called m
}
\description{
Create a parameter called m
}
\note{
this parameter is used in a log-normal distribution
  (as returned by \code{\link{create_log_normal_distr}})
  It cannot be estimated (as a hyper parameter) yet.
}
\examples{
# Create the parameter
m_param <- create_m_param()

# Use the parameter in a distribution
log_normal_distr <- create_log_normal_distr(
  m = m_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = log_normal_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_init_tree_priors.R
\name{are_init_tree_priors}
\alias{are_init_tree_priors}
\title{Determine if x consists out of initialized tree_priors objects}
\usage{
are_init_tree_priors(x)
}
\arguments{
\item{x}{the object to check if it consists out of
initialized tree_priors objects}
}
\value{
TRUE if x, or all elements of x, are initialized tree_prior objects
}
\description{
Determine if x consists out of initialized tree_priors objects
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_distr.R
\name{cep_tree_prior_to_xml_prior_distr}
\alias{cep_tree_prior_to_xml_prior_distr}
\title{Creates the tree prior section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a
Coalescent Exponential Population tree prior}
\usage{
cep_tree_prior_to_xml_prior_distr(cep_tree_prior)
}
\arguments{
\item{cep_tree_prior}{a Coalescent Exponential Population tree prior,
as returned by \code{\link{create_cep_tree_prior}}}
}
\description{
Creates the tree prior section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a
Coalescent Exponential Population tree prior
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mcmc.R
\name{check_mcmc}
\alias{check_mcmc}
\title{Check if the MCMC is a valid MCMC object.}
\usage{
check_mcmc(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
nothing
}
\description{
Calls \code{stop} if the MCMC is invalid
}
\examples{
check_mcmc(create_mcmc())
}
\seealso{
Use \code{\link{create_mcmc}} to create a valid MCMC
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_s}
\alias{parameter_to_xml_s}
\title{Internal function}
\usage{
parameter_to_xml_s(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a s parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_s_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a s parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equivalent_xmls.R
\name{are_equivalent_xml_files}
\alias{are_equivalent_xml_files}
\title{Internal function}
\usage{
are_equivalent_xml_files(filename_1, filename_2, section = NA)
}
\arguments{
\item{filename_1}{name of a first XML file}

\item{filename_2}{name of a second XML file}

\item{section}{the name of the XML section, use NA to check the whole file}
}
\value{
TRUE if the two XML files result in equivalent trees,
  FALSE otherwise
}
\description{
Internal function used for debugging to
determine if XML files result in equivalent trees
}
\seealso{
to check for equality, use \code{are_equal_xml_files}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_exp}
\alias{distr_to_xml_exp}
\title{Internal function}
\usage{
distr_to_xml_exp(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{an exponential distribution,
as created by \code{\link{create_exp_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts an exponential distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_inference_model.R
\name{init_inference_model}
\alias{init_inference_model}
\title{Initialize an inference model}
\usage{
init_inference_model(input_filename, inference_model)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
Initialize an inference model
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_alpha_param}
\alias{is_alpha_param}
\title{Determine if the object is a valid
alpha parameter}
\usage{
is_alpha_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
alpha parameter}
}
\value{
TRUE if x is a valid alpha parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
alpha parameter
}
\examples{

is_alpha_param(create_alpha_param())
is_alpha_param(create_beta_param())
is_alpha_param(create_clock_rate_param())
is_alpha_param(create_kappa_1_param())
is_alpha_param(create_kappa_2_param())
is_alpha_param(create_lambda_param())
is_alpha_param(create_m_param())
is_alpha_param(create_mean_param())
is_alpha_param(create_mu_param())
is_alpha_param(create_rate_ac_param())
is_alpha_param(create_rate_ag_param())
is_alpha_param(create_rate_at_param())
is_alpha_param(create_rate_cg_param())
is_alpha_param(create_rate_ct_param())
is_alpha_param(create_rate_gt_param())
is_alpha_param(create_s_param())
is_alpha_param(create_scale_param())
is_alpha_param(create_sigma_param())

is_alpha_param(NA)
is_alpha_param(NULL)
is_alpha_param("nonsense")
is_alpha_param(create_jc69_site_model())
is_alpha_param(create_strict_clock_model())
is_alpha_param(create_yule_tree_prior())
is_alpha_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_clock_model.R
\name{check_rln_clock_model}
\alias{check_rln_clock_model}
\title{Check if the clock model is a valid clock model.}
\usage{
check_rln_clock_model(clock_model)
}
\arguments{
\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}
}
\value{
TRUE if \code{clock_model} is a valid clock model
}
\description{
Calls \code{stop} if the clock model is invalid
}
\examples{
check_rln_clock_model(create_rln_clock_model())
}
\seealso{
Use \link{create_clock_model} to create a valid clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_xml_closing_tag.R
\name{get_xml_closing_tag}
\alias{get_xml_closing_tag}
\title{Get the XML closing tag}
\usage{
get_xml_closing_tag(text)
}
\arguments{
\item{text}{lines of XML to extract the XML closing tag from}
}
\value{
the closing tag if found, else NA
}
\description{
Get the XML closing tag
}
\examples{

# my_tag
get_xml_closing_tag("<my_tag text=something></my_tag>")

# Will return NA
get_xml_closing_tag("<my_tag text=something/>")
get_xml_closing_tag("no_xml")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_inv_gamma}
\alias{distr_to_xml_inv_gamma}
\title{Internal function}
\usage{
distr_to_xml_inv_gamma(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{an inverse-gamma distribution,
as created by \code{\link{create_inv_gamma_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts an inverse-gamma distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alignment_ids.R
\name{get_alignment_ids}
\alias{get_alignment_ids}
\title{Get the alignment IDs from one or more files.}
\usage{
get_alignment_ids(filenames)
}
\arguments{
\item{filenames}{names of the files to be checked}
}
\value{
the IDs extracted from the one or more files
}
\description{
This is done in the same way as BEAST2 does by default
The file extension will be used to determine which
type of file is worked on.
}
\examples{
  created <- get_alignment_ids(
    get_beautier_paths(c("anthus_aco.fas", "anthus_nd2.fas"))
  )
  expected <- c(
    get_alignment_id(get_beautier_path("anthus_aco.fas")),
    get_alignment_id(get_beautier_path("anthus_nd2.fas"))
  )
  testit::assert(created == expected)
}
\seealso{
Use \link{get_alignment_ids_from_fasta_filenames} to
get the alignment IDs from files known to be FASTA files
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_tree_prior.R
\name{is_cbs_tree_prior}
\alias{is_cbs_tree_prior}
\title{Determine if the object is a valid constant coalescent Bayesian skyline prior}
\usage{
is_cbs_tree_prior(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid constant coalescent
Bayesian skyline prior}
}
\value{
TRUE if x is a valid constant coalescent Bayesian skyline prior,
  FALSE otherwise
}
\description{
Determine if the object is a valid constant coalescent Bayesian skyline prior
}
\examples{
  testit::assert(!is_cbs_tree_prior(create_bd_tree_prior()))
  testit::assert( is_cbs_tree_prior(create_cbs_tree_prior()))
  testit::assert(!is_cbs_tree_prior(create_ccp_tree_prior()))
  testit::assert(!is_cbs_tree_prior(create_cep_tree_prior()))
  testit::assert(!is_cbs_tree_prior(create_yule_tree_prior()))
}
\seealso{
Use \code{\link{create_cbs_tree_prior}} to create a valid
  coalescent Bayes skyline tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_inference_model.R
\name{check_inference_model}
\alias{check_inference_model}
\title{Check if the supplied object is a valid
Bayesian phylogenetic inference model.}
\usage{
check_inference_model(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
nothing
}
\description{
Calls \code{stop} if the supplied object is not a valid
  Bayesian phylogenetic inference model.
}
\examples{
check_inference_model(create_inference_model())
}
\seealso{
Use \link{create_inference_model} to create a valid Bayesian
  phylogenetic inference model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beautier_folder.R
\name{get_beautier_folder}
\alias{get_beautier_folder}
\title{Get the path to the \link{beautier} temporary files folder}
\usage{
get_beautier_folder()
}
\value{
the path to the \link{beautier} temporary files folder
}
\description{
Get the path to the \link{beautier} temporary files folder
}
\examples{
get_beautier_folder()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_mrca_align_ids_in_fasta.R
\name{are_mrca_align_ids_in_fasta}
\alias{are_mrca_align_ids_in_fasta}
\title{Determine if the MRCA priors' alignment IDs are present in the FASTA files}
\usage{
are_mrca_align_ids_in_fasta(mrca_prior, fasta_filename)
}
\arguments{
\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{fasta_filename}{a FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.
Note that BEAST2 also supports missing data,
by using a dash (\code{-}) or question mark (\code{?})
as a sequence.}
}
\value{
TRUE if all the MRCA priors' alignment IDs
  are present in the FASTA files.
  Returns FALSE otherwise
}
\description{
Determine if the MRCA priors' alignment IDs are present in the FASTA files
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_clock_models_ids.R
\name{get_clock_models_ids}
\alias{get_clock_models_ids}
\title{Collect the IDs of the list of clock models}
\usage{
get_clock_models_ids(clock_models)
}
\arguments{
\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}
}
\value{
IDs of the clock models
}
\description{
Collect the IDs of the list of clock models
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_gamma_distr}
\alias{is_gamma_distr}
\title{Determine if the object is a valid
gamma distribution,
as created by \code{\link{create_gamma_distr}}}
\usage{
is_gamma_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
gamma distribution}
}
\value{
TRUE if x is a valid gamma distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
gamma distribution,
as created by \code{\link{create_gamma_distr}}
}
\examples{
# TRUE
is_gamma_distr(create_gamma_distr())
# FALSE
is_gamma_distr(create_inv_gamma_distr())
is_gamma_distr(NA)
is_gamma_distr(NULL)
is_gamma_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_params_doc.R
\name{default_params_doc}
\alias{default_params_doc}
\title{Documentation of general function arguments.
This function does nothing.
It is intended to inherit function argument documentation.}
\usage{
default_params_doc(
  alignment_id,
  alpha_parameter,
  bd_tree_prior,
  beautier_folder,
  cbs_tree_prior,
  beast2_version,
  beauti_options,
  beta_parameter,
  ccp_tree_prior,
  cep_tree_prior,
  chain_length,
  clock_model,
  clock_model_name,
  clock_model_names,
  clock_models,
  clock_prior_distr_id,
  clock_rate_param,
  crown_age,
  crown_ages,
  distr_id,
  fasta_filename,
  fasta_filenames,
  fixed_crown_age,
  fixed_crown_ages,
  gamma_distr,
  gamma_site_model,
  group_sizes_dimension,
  gtr_site_model,
  has_non_strict_clock_model,
  has_tip_dating,
  hky_site_model,
  id,
  ids,
  inference_model,
  inference_models,
  initial_phylogenies,
  input_filename,
  input_filenames,
  is_monophyletic,
  jc69_site_model,
  log_every,
  m_param,
  mcmc,
  mode,
  mrca_prior,
  mrca_priors,
  mrca_prior_name,
  n_init_attempts,
  output_filename,
  param,
  param_id,
  phylogeny,
  pre_burnin,
  rename_fun,
  rln_clock_model,
  sample_from_prior,
  sanitise_headers,
  screenlog,
  sequence_length,
  site_model,
  site_model_name,
  site_model_names,
  site_models,
  sort,
  store_every,
  strict_clock_model,
  taxa_names,
  tipdates_filename,
  tn93_site_model,
  tracelog,
  treelog,
  tree_prior,
  tree_prior_name,
  tree_prior_names,
  tree_priors,
  verbose,
  yule_tree_prior
)
}
\arguments{
\item{alignment_id}{ID of the alignment,
as returned by \link{get_alignment_id}.
Keep at \code{NA} to have it initialized automatically}

\item{alpha_parameter}{an alpha parameter,
as created by \link{create_alpha_param}}

\item{bd_tree_prior}{a Birth-Death tree prior, as created
by \code{\link{create_bd_tree_prior}}}

\item{beautier_folder}{the path to
the \link{beautier} temporary files folder}

\item{cbs_tree_prior}{a Coalescent Bayesian Skyline tree prior,
as returned by \code{\link{create_cbs_tree_prior}}}

\item{beast2_version}{BEAST2 version, for example, code{"2.5"}}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{beta_parameter}{a beta parameter,
as created by \link{create_beta_param}}

\item{ccp_tree_prior}{a Coalescent Constant Population tree prior,
as returned by \code{\link{create_ccp_tree_prior}}}

\item{cep_tree_prior}{a Coalescent Exponential Population tree prior,
as returned by \code{\link{create_cep_tree_prior}}}

\item{chain_length}{length of the MCMC chain}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{clock_model_name}{name of a clock model,
must be a name as returned by \code{\link{get_clock_model_names}}}

\item{clock_model_names}{one or more names of a clock model,
must be name among those returned by \code{\link{get_clock_model_names}}}

\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}

\item{clock_prior_distr_id}{ID of an MRCA clock model's distribution.
Keep at \code{NA} to have it initialized automatically}

\item{clock_rate_param}{a \code{clockRate} parameter,
a numeric value,
as created by \link{create_clock_rate_param}}

\item{crown_age}{the crown age of the phylogeny}

\item{crown_ages}{the crown ages of the phylogenies. Set to NA
if the crown age needs to be estimated}

\item{distr_id}{a distributions' ID}

\item{fasta_filename}{a FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.
Note that BEAST2 also supports missing data,
by using a dash (\code{-}) or question mark (\code{?})
as a sequence.}

\item{fasta_filenames}{One or more FASTA filenames.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{fixed_crown_age}{determines if the phylogeny's crown age is
fixed. If FALSE, crown age is estimated by BEAST2. If TRUE,
the crown age is fixed to the crown age
of the initial phylogeny.}

\item{fixed_crown_ages}{one or more booleans to determine if the
phylogenies' crown ages are fixed.
If FALSE, crown age is estimated by BEAST2. If TRUE,
the crown age is fixed to the crown age
of the initial phylogeny.}

\item{gamma_distr}{a gamma distribution,
as created by \code{\link{create_gamma_distr}})}

\item{gamma_site_model}{a site model's gamma site model,
as returned by \code{\link{create_gamma_site_model}}}

\item{group_sizes_dimension}{the group sizes' dimension,
as used by the CBS tree prior (see \code{\link{create_cbs_tree_prior}})}

\item{gtr_site_model}{a GTR site model,
as returned by \code{\link{create_gtr_site_model}}}

\item{has_non_strict_clock_model}{boolean to indicate that the is
already at least one non-strict (i.e. relaxed log-normal) clock model}

\item{has_tip_dating}{TRUE if the user has supplied tip dates,
FALSE otherwise}

\item{hky_site_model}{an HKY site model,
as returned by \code{\link{create_hky_site_model}}}

\item{id}{an alignment's IDs.
An ID can be extracted from its FASTA filename
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{ids}{one or more alignments' IDs.
IDs can be extracted from their FASTA filenames
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}

\item{inference_models}{a list of one or more inference models,
as can be created by \link{create_inference_model}}

\item{initial_phylogenies}{one or more MCMC chain's initial phylogenies.
Each one set to \code{NA} will result in BEAST2 using a random phylogeny.
Else the phylogeny is assumed to be of class \code{phylo} from the
\code{ape} package}

\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{input_filenames}{One or more FASTA filenames.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{is_monophyletic}{boolean to indicate monophyly is assumed in
a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{jc69_site_model}{a JC69 site model,
as returned by \code{\link{create_jc69_site_model}}}

\item{log_every}{number of MCMC states between writing to file}

\item{m_param}{an m parameter,
as created by \link{create_m_param}}

\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}

\item{mode}{mode how to log.
Valid values are the ones returned by \link{get_log_modes}}

\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}

\item{mrca_prior_name}{the unique name of the MRCA prior,
for example a genus, family,
order or even class name.
Leave at \link{NA} to have it named automatically.}

\item{n_init_attempts}{number of initialization attempts before failing}

\item{output_filename}{Name of the XML parameter file created by this
function. BEAST2 uses this file as input.}

\item{param}{a parameter, as can be created by \code{\link{create_param}}.}

\item{param_id}{a parameter's ID}

\item{phylogeny}{a phylogeny of type \code{phylo} from the \code{ape}
package}

\item{pre_burnin}{number of burn in samples taken before entering
the main loop}

\item{rename_fun}{a function to rename a filename,
as can be checked by \link{check_rename_fun}. This function should
have one argument, which will be a filename or \link{NA}. The
function should \link{return} one filename (when passed one filename) or
one \link{NA} (when passed one \link{NA}).
Example rename functions are:
\itemize{
  \item \link{get_remove_dir_fun} get a function that removes the directory
    paths from the filenames, in effect turning these into local files
  \item \link{get_replace_dir_fun} get a function that replaces the directory
    paths from the filenames
  \item \link{get_remove_hex_fun} get a function that removes the
    hex string from filenames.
    For example, \code{tracelog_82c1a522040.log} becomes \code{tracelog.log}
}}

\item{rln_clock_model}{a Relaxed Log-Normal clock model,
as returned by \code{\link{create_rln_clock_model}}}

\item{sample_from_prior}{set to \link{TRUE} to sample from the prior}

\item{sanitise_headers}{set to \link{TRUE} to sanitise the headers of the
log file}

\item{screenlog}{a \code{screenlog},
as created by \link{create_screenlog}}

\item{sequence_length}{a DNA sequence length, in base pairs}

\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}

\item{site_model_name}{name of a site model,
must be a name as returned by \code{\link{get_site_model_names}}}

\item{site_model_names}{one or more names of a site model,
must be name among those returned by \code{\link{get_site_model_names}}}

\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}

\item{sort}{how to sort the log.
Valid values are the ones returned by \link{get_log_sorts}}

\item{store_every}{number of states the MCMC will process
before the posterior's state will be saved to file.
Use -1 or \code{NA} to use the default frequency.}

\item{strict_clock_model}{a strict clock model,
as returned by \code{\link{create_strict_clock_model}}}

\item{taxa_names}{names of the taxa,
as returned by \code{\link{get_taxa_names}}.
Keep at \code{NA} to have it initialized automatically,
using all taxa in the alignment}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}

\item{tn93_site_model}{a TN93 site model,
as returned by \code{\link{create_tn93_site_model}}}

\item{tracelog}{a \code{tracelog},
as created by \link{create_tracelog}}

\item{treelog}{a \code{treelog},
as created by \link{create_treelog}}

\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{tree_prior_name}{name of a tree prior,
must be a name as returned by \code{\link{get_tree_prior_names}}}

\item{tree_prior_names}{one or more names of a tree prior,
must be a name among those returned by \code{\link{get_tree_prior_names}}}

\item{tree_priors}{one or more tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}

\item{yule_tree_prior}{a Yule tree_prior,
as created by \code{\link{create_yule_tree_prior}}}
}
\description{
Documentation of general function arguments.
This function does nothing.
It is intended to inherit function argument documentation.
}
\note{
This is an internal function, so it should be marked with
  \code{@export}. This is not done, as this will disallow all
  functions to find the documentation parameters
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_site_model.R
\name{is_init_gtr_site_model}
\alias{is_init_gtr_site_model}
\title{Determine if x is an initialized GTR site model
as created by \code{\link{create_gtr_site_model}}}
\usage{
is_init_gtr_site_model(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized GTR site model}
}
\value{
TRUE if x is an initialized GTR site model
}
\description{
Determine if x is an initialized GTR site model
as created by \code{\link{create_gtr_site_model}}
}
\examples{

gtr_site_model <- create_gtr_site_model()
# FALSE: not yet initialized
is_init_gtr_site_model(gtr_site_model)
gtr_site_model <- init_gtr_site_model(gtr_site_model)
# TRUE: now it is initialized
is_init_gtr_site_model(gtr_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_remove_hex_fun.R
\name{get_remove_hex_fun}
\alias{get_remove_hex_fun}
\title{Get a function that removes the hex string from filenames.}
\usage{
get_remove_hex_fun()
}
\description{
The default filenames created by \link{beautier} are temporary files,
such as \code{/home/john/.cache/tracelog_82c5888db98.log} (on Linux),
where \code{/home/john/.cache} is the location to a temporary folder
(on Linux) and \code{tracelog_82c5888db98.log} the filename.
The filename ends with a hex string (as is common for temporary files,
as \link{tempfile} does so).
Because \link{beautier} puts an underscore
between the filename description (\code{tracelog}) and the hex
string, this function removes both.
}
\examples{

f <- get_remove_hex_fun()
# /home/john/beast2.xml.state
f("/home/john/beast2_186c7404208c.xml.state")

# beast2.xml.state
f("beast2_186c7404208c.xml.state")

# NA
f(NA)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_mrca_prior.R
\name{is_mrca_prior}
\alias{is_mrca_prior}
\title{Determine of the object is an empty (\code{NA}) or valid MRCA prior.}
\usage{
is_mrca_prior(mrca_prior)
}
\arguments{
\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}
}
\value{
TRUE if \code{x} is an MRCA prior, FALSE otherwise
}
\description{
Determine of the object is an empty (\code{NA}) or valid MRCA prior.
}
\examples{
# TRUE
is_mrca_prior(create_mrca_prior())
# Also 'NA' is a valid MRCA prior,
# denoting that there no MRCA priors
is_mrca_prior(NA)

# FALSE
is_mrca_prior(NULL)
is_mrca_prior("nonsense")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_param}
\alias{create_param}
\title{General function to create a parameter.}
\usage{
create_param(name, id, value, ...)
}
\arguments{
\item{name}{the parameters' name. Valid
names can be found in \code{get_param_names}}

\item{id}{the parameter's ID}

\item{value}{value of the parameter}

\item{...}{specific parameter parameters}
}
\value{
a parameter
}
\description{
General function to create a parameter.
}
\note{
Prefer using the
  named functions
  \code{\link{create_alpha_param}},
  \code{\link{create_beta_param}},
  \code{\link{create_clock_rate_param}},
  \code{\link{create_kappa_1_param}},
  \code{\link{create_kappa_2_param}},
  \code{\link{create_lambda_param}},
  \code{\link{create_m_param}},
  \code{\link{create_mean_param}},
  \code{\link{create_mu_param}},
  \code{\link{create_rate_ac_param}},
  \code{\link{create_rate_ag_param}},
  \code{\link{create_rate_at_param}},
  \code{\link{create_rate_cg_param}},
  \code{\link{create_rate_ct_param}},
  \code{\link{create_rate_gt_param}},
  \code{\link{create_s_param}},
  \code{\link{create_scale_param}},
  and \code{\link{create_sigma_param}}
}
\examples{
# Create an alpha parameter
alpha_param <- create_alpha_param()

# Use the parameter in a distribution
beta_distr <- create_beta_distr(
  alpha = alpha_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = beta_distr
  )
)
file.remove(beast2_input_file)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_tracelog.R
\name{check_tracelog}
\alias{check_tracelog}
\title{Check if a \code{tracelog} is valid.}
\usage{
check_tracelog(tracelog)
}
\arguments{
\item{tracelog}{a \code{tracelog},
as created by \link{create_tracelog}}
}
\description{
Will call \link{stop} if not.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_from_model.R
\name{create_beast2_input_from_model}
\alias{create_beast2_input_from_model}
\title{Create a BEAST2 XML input text from an inference model}
\usage{
create_beast2_input_from_model(input_filename, inference_model)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
a character vector of XML strings
}
\description{
The main two XML tags are these:
\preformatted{
  <?xml[...]?><beast[...]>
  [...]
  </beast>
}
}
\examples{
text <- create_beast2_input_from_model(
  input_filename = get_fasta_filename(),
  inference_model = create_inference_model()
)
}
\seealso{
Use \link{create_beast2_input_file_from_model} to also save it to file.
Use \link{create_xml_declaration}
to create the XML text of the XML declaration.
Use \link{create_beast2_input_beast} to create
to create the XML text of the \code{beast} tag.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_models_to_xml_prior_distr.R
\name{site_models_to_xml_prior_distr}
\alias{site_models_to_xml_prior_distr}
\title{Represent the site models as XML}
\usage{
site_models_to_xml_prior_distr(site_models)
}
\arguments{
\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}
}
\value{
lines of XML text
}
\description{
Represent the site models as XML
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_uniform_distr}
\alias{create_uniform_distr}
\alias{create_distr_uniform}
\title{Create a uniform distribution}
\usage{
create_uniform_distr(id = NA, value = NA, lower = NA, upper = Inf)
}
\arguments{
\item{id}{the distribution's ID}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
a uniform distribution
}
\description{
Create a uniform distribution
}
\examples{
uniform_distr <- create_uniform_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = uniform_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_xml_short_closing_tag.R
\name{has_xml_short_closing_tag}
\alias{has_xml_short_closing_tag}
\title{Is an XML closing tag with short closing text in
one of the lines of the text?}
\usage{
has_xml_short_closing_tag(lines)
}
\arguments{
\item{lines}{lines of an XML text}
}
\value{
TRUE if there is an XML tag that also closes present in the lines
  of text, FALSE otherwise
}
\description{
Is an XML closing tag with short closing text in
one of the lines of the text?
}
\examples{
# TRUE
has_xml_short_closing_tag("<my_tag id=1/>")
# FALSE
has_xml_short_closing_tag("<my_tag id=1>text</my_tag>")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_replace_dir_fun.R
\name{get_replace_dir_fun}
\alias{get_replace_dir_fun}
\title{Get a function to replace the directory of a filename}
\usage{
get_replace_dir_fun(new_dir_name = "")
}
\arguments{
\item{new_dir_name}{the new directory name}
}
\description{
Get a function to replace the directory of a filename
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_ucld_stdev_state_node_param_xml.R
\name{create_ucld_stdev_state_node_param_xml}
\alias{create_ucld_stdev_state_node_param_xml}
\title{Internal function}
\usage{
create_ucld_stdev_state_node_param_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the following XML:
\code{
  <parameter id="ucldStdev.c:[id]" lower="0.0" name="stateNode">
    0.1
  </parameter>
}
}
\description{
Creates the \code{ucldStdev} parameter with the name \code{stateNode},
such as:
\code{
  <parameter id="ucldStdev.c:[id]" [...] name="stateNode">0.1</parameter>
}
}
\examples{
create_ucld_stdev_state_node_param_xml(
  create_inference_model(
    clock_model = create_rln_clock_model(id = 314)
  )
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_clock_model.R
\name{create_rln_clock_model}
\alias{create_rln_clock_model}
\alias{create_clock_model_rln}
\title{Create a relaxed log-normal clock model}
\usage{
create_rln_clock_model(
  id = NA,
  mean_rate_prior_distr = create_uniform_distr(),
  ucldstdev_distr = create_gamma_distr(),
  mparam_id = NA,
  mean_clock_rate = "1.0",
  n_rate_categories = -1,
  normalize_mean_clock_rate = FALSE,
  dimension = NA
)
}
\arguments{
\item{id}{an alignment's IDs.
An ID can be extracted from its FASTA filename
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{mean_rate_prior_distr}{the mean clock rate prior distribution,
as created by a \code{\link{create_distr}} function}

\item{ucldstdev_distr}{the standard deviation of the uncorrelated
log-normal distribution,
as created by a \code{\link{create_distr}} function}

\item{mparam_id}{the ID of the M parameter in the \code{branchRateModel},
set to NA to have it initialized}

\item{mean_clock_rate}{the mean clock rate, 1.0 by default
(is called \code{ucld_stdev} in XML, where \code{ucld_stdev} is always 0.1)}

\item{n_rate_categories}{the number of rate categories.
-1 is default,
0 denotes as much rates as branches}

\item{normalize_mean_clock_rate}{normalize the mean clock rate}

\item{dimension}{the dimensionality of the relaxed clock model.
Leave NA to let beautier calculate it.
Else, the dimensionality of the clock
equals twice the number of taxa minus two.}
}
\value{
a relaxed log-normal clock_model
}
\description{
Create a relaxed log-normal clock model
}
\examples{
rln_clock_model <- create_rln_clock_model()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  clock_model = rln_clock_model
)
file.remove(beast2_input_file)

rln_clock_model_exp <- create_rln_clock_model(
  mean_rate_prior_distr = create_exp_distr()
)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  clock_model = rln_clock_model_exp
)
file.remove(beast2_input_file)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_clock_model.R
\name{is_init_rln_clock_model}
\alias{is_init_rln_clock_model}
\title{Determine if x is an initialized relaxed log-normal clock_model object}
\usage{
is_init_rln_clock_model(rln_clock_model)
}
\arguments{
\item{rln_clock_model}{a Relaxed Log-Normal clock model,
as returned by \code{\link{create_rln_clock_model}}}
}
\value{
TRUE if x is an initialized relaxed log-normal clock_model object,
  FALSE otherwise
}
\description{
Determine if x is an initialized relaxed log-normal clock_model object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equivalent_xmls.R
\name{are_equivalent_xml_lines_all}
\alias{are_equivalent_xml_lines_all}
\title{Determine if XML lines result in equivalent trees}
\usage{
are_equivalent_xml_lines_all(lines_1, lines_2, verbose = FALSE)
}
\arguments{
\item{lines_1}{lines of a first XML file}

\item{lines_2}{lines of a second XML file}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
TRUE if the two XML lines result in equivalent trees,
  FALSE otherwise
}
\description{
Determine if XML lines result in equivalent trees
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_model_to_xml_tracelog.R
\name{clock_model_to_xml_tracelog}
\alias{clock_model_to_xml_tracelog}
\title{Internal function}
\usage{
clock_model_to_xml_tracelog(
  inference_model,
  clock_model = "deprecated",
  mrca_priors = "deprecated"
)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}
}
\value{
a character vector of XML strings
}
\description{
Creates the clock model's XML for the tracelog section
}
\examples{
# <logger id="tracelog" ...>
#'   # Here
# </logger>
}
\seealso{
all clock models' tracelog section is created
  by \code{\link{clock_models_to_xml_tracelog}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site_model_names.R
\name{get_site_model_names}
\alias{get_site_model_names}
\title{Get the site models' names}
\usage{
get_site_model_names()
}
\value{
the site model names
}
\description{
Get the site models' names
}
\examples{
  # Check all names
  names <- get_site_model_names()
  testit::assert("JC69" \%in\% names)
  testit::assert("HKY" \%in\% names)
  testit::assert("TN93" \%in\% names)
  testit::assert("GTR" \%in\% names)
}
\seealso{
Use \link{create_site_models} to get all site models
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_param.R
\name{is_init_param}
\alias{is_init_param}
\title{Determine if x is an initialized parameter,
  as created by \link{create_param}}
\usage{
is_init_param(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized parameter}
}
\value{
\link{TRUE} if x is an initialized parameter,
\link{FALSE} otherwise
}
\description{
Determine if x is an initialized parameter,
  as created by \link{create_param}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_beta_param}
\alias{is_beta_param}
\title{Determine if the object is a valid
beta parameter}
\usage{
is_beta_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
beta parameter}
}
\value{
TRUE if x is a valid beta parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
beta parameter
}
\examples{

is_beta_param(create_alpha_param())
is_beta_param(create_beta_param())
is_beta_param(create_clock_rate_param())
is_beta_param(create_kappa_1_param())
is_beta_param(create_kappa_2_param())
is_beta_param(create_lambda_param())
is_beta_param(create_m_param())
is_beta_param(create_mean_param())
is_beta_param(create_mu_param())
is_beta_param(create_rate_ac_param())
is_beta_param(create_rate_ag_param())
is_beta_param(create_rate_at_param())
is_beta_param(create_rate_cg_param())
is_beta_param(create_rate_ct_param())
is_beta_param(create_rate_gt_param())
is_beta_param(create_s_param())
is_beta_param(create_scale_param())
is_beta_param(create_sigma_param())

is_beta_param(NA)
is_beta_param(NULL)
is_beta_param("nonsense")
is_beta_param(create_jc69_site_model())
is_beta_param(create_strict_clock_model())
is_beta_param(create_yule_tree_prior())
is_beta_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equal_mcmcs.R
\name{are_equal_mcmcs}
\alias{are_equal_mcmcs}
\title{Determine if two MCMCs are equal.}
\usage{
are_equal_mcmcs(mcmc_1, mcmc_2)
}
\arguments{
\item{mcmc_1}{an MCMC, as created by \code{\link{create_mcmc}}}

\item{mcmc_2}{an MCMC, as created by \code{\link{create_mcmc}}}
}
\value{
TRUE if the two MCMCs are equal
}
\description{
Will \link{stop} if the arguments are not MCMCs.
}
\examples{
mcmc_1 <- create_mcmc(chain_length = 1000)
mcmc_2 <- create_mcmc(chain_length = 314)
# TRUE
are_equal_mcmcs(mcmc_1, mcmc_1)
# FALSE
are_equal_mcmcs(mcmc_1, mcmc_2)
}
\seealso{
Use \code{\link{create_mcmc}} to create an MCMC
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_lambda}
\alias{parameter_to_xml_lambda}
\title{Internal function}
\usage{
parameter_to_xml_lambda(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a lambda parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_lambda_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a lambda parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_mrca_prior.R
\name{has_mrca_prior}
\alias{has_mrca_prior}
\title{Determines if the inference model has an MRCA prior.}
\usage{
has_mrca_prior(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
TRUE if the inference model has an MRCA prior,
  FALSE otherwise
}
\description{
Will \link{stop} if the inference model is invalid
}
\note{
MRCA: 'Most Recent Common Ancestor'
}
\examples{
# No MRCA prior
inference_model <- create_inference_model(
  mrca_prior = NA
)
has_mrca_prior(inference_model) # Returns FALSE

# A default MRCA prior
inference_model <- create_inference_model(
  mrca_prior = create_mrca_prior()
)
has_mrca_prior(inference_model) # Returns TRUE
}
\seealso{
\itemize{
    \item \code{\link{create_inference_model}}: create an inference model
    \item \code{\link{create_mrca_prior}}: create an MRCA prior
  }
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_mu_param}
\alias{create_mu_param}
\alias{create_param_mu}
\title{Create a parameter called mu}
\usage{
create_mu_param(id = NA, value = 0)
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}
}
\value{
a parameter called mu
}
\description{
Create a parameter called mu
}
\note{
this parameter is used in a Laplace distribution
  (as returned by \code{\link{create_laplace_distr}}).
  It cannot be estimated (as a hyper parameter) yet.
}
\examples{
# Create the parameter
mu_param <- create_mu_param()

# Use the parameter in a distribution
laplace_distr <- create_laplace_distr(
  mu = mu_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = laplace_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_to_xml_run.R
\name{mcmc_to_xml_run}
\alias{mcmc_to_xml_run}
\title{Converts an MCMC object to the run section's XML}
\usage{
mcmc_to_xml_run(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
the XML as text
}
\description{
Converts an MCMC object to the run section's XML
}
\examples{
# <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"1e+07\">
mcmc_to_xml_run(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_param.R
\name{check_param}
\alias{check_param}
\title{Check if the parameter is a valid parameter}
\usage{
check_param(param)
}
\arguments{
\item{param}{a parameter, as can be created by \code{\link{create_param}}.}
}
\value{
nothing
}
\description{
Calls \code{stop} if the parameter is invalid
}
\examples{
check_param(create_alpha_param())
check_param(create_beta_param())
}
\seealso{
Use \link{create_param} to create a valid parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_distr}
\alias{init_distr}
\title{Initializes a distribution}
\usage{
init_distr(distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{distr}{a distribution,
using \code{\link{create_distr}}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized distribution
}
\description{
Initializes a distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mcmc_nested_sampling.R
\name{check_ns_mcmc}
\alias{check_ns_mcmc}
\alias{check_mcmc_nested_sampling}
\alias{check_nested_sampling_mcmc}
\title{Check if this an MCMC that uses Nested Sampling
to estimate a marginal likelihood.}
\usage{
check_ns_mcmc(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\description{
Will \link{stop} if not, else will do nothing
}
\seealso{
use \code{\link{create_ns_mcmc}}
to create an MCMC that uses Nested Sampling
to estimate a marginal likelihood
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_rate_categories_state_node_xml.R
\name{create_rate_categories_state_node_xml}
\alias{create_rate_categories_state_node_xml}
\title{Internal function}
\usage{
create_rate_categories_state_node_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the following XML:
\code{
  "<stateNode id=\"rateCategories.c:[id]\"
  spec=\"parameter.IntegerParameter\" dimension=\"[dimension]\">
  1
  </stateNode>"
}
}
\description{
Creates the \code{rateCategories} state node,
such as:
\code{
  "<stateNode id=\"rateCategories.c:[id]\"
    spec=\"parameter.IntegerParameter\"
    dimension=\"[dimension]\">
  1
  </stateNode>"
}
}
\examples{
create_rate_categories_state_node_xml(
  create_inference_model(
    clock_model = create_rln_clock_model(
      id = 314,
      dimension = 1
    )
  )
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_last_regex_line.R
\name{find_last_regex_line}
\alias{find_last_regex_line}
\title{Find the index of the last line that matches a regex}
\usage{
find_last_regex_line(lines, regex)
}
\arguments{
\item{lines}{lines of text}

\item{regex}{regex string}
}
\value{
index of the line
}
\description{
Find the index of the last line that matches a regex
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_laplace_distr}
\alias{init_laplace_distr}
\title{Initializes an Laplace distribution}
\usage{
init_laplace_distr(laplace_distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{laplace_distr}{a Laplace distribution,
using \link{create_laplace_distr}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized Laplace distribution
}
\description{
Initializes an Laplace distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_clock_model.R
\name{find_clock_model}
\alias{find_clock_model}
\title{Finds a clock model with a certain ID}
\usage{
find_clock_model(clock_models, id)
}
\arguments{
\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}

\item{id}{the ID of the clock model}
}
\value{
the clock models with the desired ID, NULL if such a clock model is
  absent
}
\description{
Finds a clock model with a certain ID
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_mcmc.R
\name{is_mcmc_nested_sampling}
\alias{is_mcmc_nested_sampling}
\alias{is_nested_sampling_mcmc}
\title{Determine if the object is a valid Nested-Sampling MCMC,
  as used in [1]}
\usage{
is_mcmc_nested_sampling(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid MCMC}
}
\value{
TRUE if x is a valid Nested-Sampling MCMC, FALSE otherwise
}
\description{
Determine if the object is a valid Nested-Sampling MCMC,
  as used in [1]
}
\examples{
# TRUE
is_nested_sampling_mcmc(create_ns_mcmc())
# FALSE
is_nested_sampling_mcmc(create_mcmc())
is_nested_sampling_mcmc("nonsense")
}
\references{
* [1] Patricio Maturana Russel, Brendon J Brewer, Steffen Klaere,
    Remco R Bouckaert; Model Selection and Parameter Inference in
    Phylogenetics Using Nested Sampling, Systematic Biology, 2018,
    syy050, https://doi.org/10.1093/sysbio/syy050
}
\seealso{
Use \link{create_ns_mcmc} to create an NS MCMC
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equivalent_xmls.R
\name{are_equivalent_xml_lines_operators}
\alias{are_equivalent_xml_lines_operators}
\title{Determine if XML operator lines result in equivalent trees}
\usage{
are_equivalent_xml_lines_operators(lines_1, lines_2, verbose = FALSE)
}
\arguments{
\item{lines_1}{lines of a first XML file}

\item{lines_2}{lines of a second XML file}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
TRUE if the two XML lines result in equivalent trees,
  FALSE otherwise
}
\description{
Determine if XML operator lines result in equivalent trees
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_clock_model.R
\name{is_clock_model}
\alias{is_clock_model}
\title{Determine if the object is a valid clock_model}
\usage{
is_clock_model(x)
}
\arguments{
\item{x}{an object, to be determined if it is a clock_model}
}
\value{
TRUE if the clock_model is a valid clock_model, FALSE otherwise
}
\description{
Determine if the object is a valid clock_model
}
\examples{
# TRUE
is_clock_model(create_strict_clock_model())
is_clock_model(create_rln_clock_model())

# FALSE
is_clock_model(NA)
is_clock_model(NULL)
is_clock_model("nonsense")
is_clock_model(create_jc69_site_model())
is_clock_model(create_mcmc())
}
\seealso{
see \code{\link{create_clock_model}} for an overview of functions
  to create valid clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_tracelog.R
\name{check_tracelog_values}
\alias{check_tracelog_values}
\title{Check if the tracelog has the list elements with valid values
for being a valid tracelog object.}
\usage{
check_tracelog_values(tracelog)
}
\arguments{
\item{tracelog}{a \code{tracelog},
as created by \link{create_tracelog}}
}
\value{
nothing
}
\description{
Calls \code{stop} if a value is invalid
}
\seealso{
Use \link{create_tracelog} to create a valid tracelog
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tree_prior.R
\name{create_tree_prior}
\alias{create_tree_prior}
\title{Internal function to create a tree prior}
\usage{
create_tree_prior(name, id, ...)
}
\arguments{
\item{name}{the tree prior name. Can be any name
in \code{get_tree_prior_names}}

\item{id}{the ID of the alignment}

\item{...}{specific tree prior parameters}
}
\value{
a tree_prior
}
\description{
Internal function to create a tree prior
}
\note{
Prefer the use the named functions
  \code{\link{create_bd_tree_prior}},
  \code{\link{create_cbs_tree_prior}},
  \code{\link{create_ccp_tree_prior}}
  \code{\link{create_cep_tree_prior}}
  and \code{\link{create_yule_tree_prior}}
  instead
}
\examples{
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_bd_tree_prior()
)
file.remove(beast2_input_file)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_beautier_path("test_output_6.fas"),
  beast2_input_file,
  tree_prior = create_cbs_tree_prior()
)
file.remove(beast2_input_file)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_ccp_tree_prior()
)
file.remove(beast2_input_file)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_cep_tree_prior()
)
file.remove(beast2_input_file)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior()
)
file.remove(beast2_input_file)
}
\seealso{
See
  \code{\link{create_bd_tree_prior}},
  \code{\link{create_cbs_tree_prior}},
  \code{\link{create_ccp_tree_prior}}
  \code{\link{create_cep_tree_prior}}
  and \code{\link{create_yule_tree_prior}}
  for more examples using those functions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_kappa_2_param}
\alias{create_kappa_2_param}
\alias{create_param_kappa_2}
\title{Create a parameter called kappa 2}
\usage{
create_kappa_2_param(id = NA, lower = "0.0", value = "2.0", estimate = TRUE)
}
\arguments{
\item{id}{the parameter's ID}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}

\item{value}{value of the parameter}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}
}
\value{
a parameter called kappa 2
}
\description{
Create a parameter called kappa 2
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_clock_model.R
\name{is_strict_clock_model}
\alias{is_strict_clock_model}
\title{Determine if the object is a valid strict clock model,
  as returned by \code{\link{create_strict_clock_model}}}
\usage{
is_strict_clock_model(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid strict clock model}
}
\value{
TRUE if x is a valid strict clock model, FALSE otherwise
}
\description{
Determine if the object is a valid strict clock model,
  as returned by \code{\link{create_strict_clock_model}}
}
\examples{

is_strict_clock_model(create_strict_clock_model())
is_strict_clock_model(create_rln_clock_model())

is_strict_clock_model(NA)
is_strict_clock_model(NULL)
is_strict_clock_model("nonsense")
is_strict_clock_model(create_jc69_site_model())
is_strict_clock_model(create_mcmc())
}
\seealso{
\code{\link{create_clock_model}} shows an overview of
  functions to create a clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_poisson_distr}
\alias{is_poisson_distr}
\title{Determine if the object is a valid
Poisson distribution
as created by \code{\link{create_poisson_distr}}}
\usage{
is_poisson_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
Poisson distribution}
}
\value{
TRUE if x is a valid Poisson distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
Poisson distribution
as created by \code{\link{create_poisson_distr}}
}
\examples{
# TRUE
is_poisson_distr(create_poisson_distr())
# FALSE
is_poisson_distr(create_uniform_distr())
is_distr(NA)
is_distr(NULL)
is_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_priors_to_xml_operators.R
\name{tree_priors_to_xml_operators}
\alias{tree_priors_to_xml_operators}
\title{Deprecated}
\usage{
tree_priors_to_xml_operators(
  tree_priors = "deprecated",
  fixed_crown_ages = "deprecated"
)
}
\arguments{
\item{tree_priors}{one or more tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{fixed_crown_ages}{one or more booleans to determine if the
phylogenies' crown ages are fixed.
If FALSE, crown age is estimated by BEAST2. If TRUE,
the crown age is fixed to the crown age
of the initial phylogeny.}
}
\value{
the tree priors as XML text
}
\description{
Creates the XML of a list of one or more tree priors,
  as used in the \code{operators} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_poisson_distr}
\alias{init_poisson_distr}
\title{Initializes an Poisson distribution}
\usage{
init_poisson_distr(poisson_distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{poisson_distr}{a Poisson distribution,
using \link{create_poisson_distr}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized Poisson distribution
}
\description{
Initializes an Poisson distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_subst_model_xml.R
\name{create_subst_model_xml}
\alias{create_subst_model_xml}
\title{Internal function to create the \code{substModel} section}
\usage{
create_subst_model_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the \code{substModel} section as XML text
}
\description{
Internal function to create the \code{substModel} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_mrca_taxon_names_in_fasta.R
\name{are_mrca_taxon_names_in_fasta}
\alias{are_mrca_taxon_names_in_fasta}
\title{Determine if the MRCA priors' taxa names are present in the FASTA files}
\usage{
are_mrca_taxon_names_in_fasta(mrca_prior, fasta_filename)
}
\arguments{
\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{fasta_filename}{a FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.
Note that BEAST2 also supports missing data,
by using a dash (\code{-}) or question mark (\code{?})
as a sequence.}
}
\value{
TRUE if the MRCA priors' taxa names are
  present in the FASTA files. FALSE otherwise.
}
\description{
Determine if the MRCA priors' taxa names are present in the FASTA files
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_tree_prior.R
\name{is_init_tree_prior}
\alias{is_init_tree_prior}
\title{Determine if x is an initialized tree_prior objects}
\usage{
is_init_tree_prior(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized tree_priors object}
}
\value{
TRUE if x is an initialized tree_prior object
}
\description{
Determine if x is an initialized tree_prior objects
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_site_model.R
\name{is_jc69_site_model}
\alias{is_jc69_site_model}
\title{Determine if the object is a valid JC69 site model}
\usage{
is_jc69_site_model(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid JC69 site model}
}
\value{
TRUE if x is a valid JC69 site model, FALSE otherwise
}
\description{
Determine if the object is a valid JC69 site model
}
\examples{

# site models
is_jc69_site_model(create_gtr_site_model())
is_jc69_site_model(create_hky_site_model())
is_jc69_site_model(create_jc69_site_model())
is_jc69_site_model(create_tn93_site_model())

# other models
is_jc69_site_model(NA)
is_jc69_site_model(NULL)
is_jc69_site_model("nonsense")
is_jc69_site_model(create_strict_clock_model())
is_jc69_site_model(create_bd_tree_prior())
is_jc69_site_model(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_mu_param}
\alias{is_mu_param}
\title{Determine if the object is a valid
mu parameter}
\usage{
is_mu_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
mu parameter}
}
\value{
TRUE if x is a valid mu parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
mu parameter
}
\examples{

is_mu_param(create_alpha_param())
is_mu_param(create_beta_param())
is_mu_param(create_clock_rate_param())
is_mu_param(create_kappa_1_param())
is_mu_param(create_kappa_2_param())
is_mu_param(create_lambda_param())
is_mu_param(create_m_param())
is_mu_param(create_mean_param())
is_mu_param(create_mu_param())
is_mu_param(create_rate_ac_param())
is_mu_param(create_rate_ag_param())
is_mu_param(create_rate_at_param())
is_mu_param(create_rate_cg_param())
is_mu_param(create_rate_ct_param())
is_mu_param(create_rate_gt_param())
is_mu_param(create_s_param())
is_mu_param(create_scale_param())
is_mu_param(create_sigma_param())

is_mu_param(NA)
is_mu_param(NULL)
is_mu_param("nonsense")
is_mu_param(create_jc69_site_model())
is_mu_param(create_strict_clock_model())
is_mu_param(create_yule_tree_prior())
is_mu_param(create_mcmc())
}
\seealso{
\code{\link{create_mu_param}} creates a mu parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_kappa_2}
\alias{parameter_to_xml_kappa_2}
\title{Internal function}
\usage{
parameter_to_xml_kappa_2(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a kappa 2 parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_kappa_2_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a kappa 2 parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_branch_rate_model_xml.R
\name{create_branch_rate_model_sc_xml}
\alias{create_branch_rate_model_sc_xml}
\title{Internal function to call \link{create_branch_rate_model_xml}
for a strict clock.}
\usage{
create_branch_rate_model_sc_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
a character vector of XML strings
}
\description{
Internal function to call \link{create_branch_rate_model_xml}
for a strict clock.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_clock_model.R
\name{check_strict_clock_model}
\alias{check_strict_clock_model}
\title{Check if the clock model is a valid clock model.}
\usage{
check_strict_clock_model(clock_model)
}
\arguments{
\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}
}
\value{
TRUE if \code{clock_model} is a valid clock model
}
\description{
Calls \code{stop} if the clock model is invalid
}
\examples{
 check_strict_clock_model(create_strict_clock_model())
}
\seealso{
Use \link{create_clock_model} to create a valid clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_rate_ct_param}
\alias{create_rate_ct_param}
\alias{create_param_rate_ct}
\title{Create a parameter called 'rate CT'}
\usage{
create_rate_ct_param(id = NA, value = "1.0", lower = "0.0")
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}
}
\value{
a parameter called 'rate CT'
}
\description{
Create a parameter called 'rate CT'
}
\examples{
# Create parameter
rate_ct_param <- create_rate_ct_param(value = 1)

# Use the parameter to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  site_model = create_gtr_site_model(
    rate_ct_param = rate_ct_param
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_tree_priors.R
\name{init_bd_tree_prior}
\alias{init_bd_tree_prior}
\title{Initializes a Birth-Death tree prior}
\usage{
init_bd_tree_prior(bd_tree_prior, distr_id, param_id)
}
\arguments{
\item{bd_tree_prior}{a Birth-Death tree prior, as created
by \code{\link{create_bd_tree_prior}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized Birth-Death tree prior
}
\description{
Initializes a Birth-Death tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_empty_lines.R
\name{remove_empty_lines}
\alias{remove_empty_lines}
\title{Remove all lines that are only whitespace}
\usage{
remove_empty_lines(lines, trim = FALSE)
}
\arguments{
\item{lines}{character vector with text}

\item{trim}{FALSE if indentation must be preserved,
TRUE will remove all surrounding whitespace}
}
\value{
the lines with text
}
\description{
Remove all lines that are only whitespace
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_distr.R
\name{yule_tree_prior_to_xml_prior_distr}
\alias{yule_tree_prior_to_xml_prior_distr}
\title{Creates the \code{prior} section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a Yule tree prior}
\usage{
yule_tree_prior_to_xml_prior_distr(yule_tree_prior)
}
\arguments{
\item{yule_tree_prior}{a Yule tree_prior,
as created by \code{\link{create_yule_tree_prior}}}
}
\description{
Creates the \code{prior} section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a Yule tree prior
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_distr}
\alias{is_init_distr}
\title{Determine if x is an initialized distribution object
  as created by \code{\link{create_distr}}}
\usage{
is_init_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized distribution object}
}
\value{
TRUE if x is an initialized distribution object
}
\description{
Determine if x is an initialized distribution object
  as created by \code{\link{create_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename_inference_model_filenames.R
\name{rename_inference_model_filenames}
\alias{rename_inference_model_filenames}
\title{Rename the filenames in an inference model}
\usage{
rename_inference_model_filenames(inference_model, rename_fun)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}

\item{rename_fun}{a function to rename a filename,
as can be checked by \link{check_rename_fun}. This function should
have one argument, which will be a filename or \link{NA}. The
function should \link{return} one filename (when passed one filename) or
one \link{NA} (when passed one \link{NA}).
Example rename functions are:
\itemize{
  \item \link{get_remove_dir_fun} get a function that removes the directory
    paths from the filenames, in effect turning these into local files
  \item \link{get_replace_dir_fun} get a function that replaces the directory
    paths from the filenames
  \item \link{get_remove_hex_fun} get a function that removes the
    hex string from filenames.
    For example, \code{tracelog_82c1a522040.log} becomes \code{tracelog.log}
}}
}
\value{
an inference model with the renamed filenames
}
\description{
Rename the filenames in an inference model
}
\examples{

inference_model <- create_inference_model()
inference_model$mcmc$tracelog$filename <- "trace.log"
inference_model$mcmc$screenlog$filename <- "screen.log"
inference_model$mcmc$treelog$filename <- "tree.log"
inference_model$tipdates_filename <- "tipdates.csv"

# Nah, put the files in a folder
inference_model <- rename_inference_model_filenames(
  inference_model = inference_model,
  rename_fun = get_replace_dir_fun("/home/john")
)

# Nah, put the files in anoth folder
inference_model <- rename_inference_model_filenames(
  inference_model = inference_model,
  rename_fun = get_replace_dir_fun("/home/doe")
)

# Nah, store the files locally
rename_inference_model_filenames(
  inference_model = inference_model,
  rename_fun = get_remove_dir_fun()
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_uniform_distr}
\alias{init_uniform_distr}
\title{Initializes a uniform distribution}
\usage{
init_uniform_distr(uniform_distr, distr_id = 0)
}
\arguments{
\item{uniform_distr}{a uniform distribution,
using \link{create_uniform_distr}}

\item{distr_id}{the first distribution's ID}
}
\value{
an initialized uniform distribution
}
\description{
Initializes a uniform distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mcmc_to_xml_run.R
\name{mcmc_to_xml_run_default}
\alias{mcmc_to_xml_run_default}
\title{Converts an MCMC object to the run section's XML for a default MCMC}
\usage{
mcmc_to_xml_run_default(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
the XML as text
}
\description{
Converts an MCMC object to the run section's XML for a default MCMC
}
\examples{
# <run id=\"mcmc\" spec=\"MCMC\" chainLength=\"1e+07\">
xml <- mcmc_to_xml_run_default(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_xml_opening_tag.R
\name{has_xml_opening_tag}
\alias{has_xml_opening_tag}
\title{Is an XML opening tag with value 'section' present among the lines of
  the text?}
\usage{
has_xml_opening_tag(lines, section = NA)
}
\arguments{
\item{lines}{lines of an XML text}

\item{section}{if NA, this function returns TRUE if there is any
XML opening tag. If \code{section} is set to a certain word,
this function returns TRUE if that tag matches \code{section}}
}
\value{
lines of XML text
}
\description{
Is an XML opening tag with value 'section' present among the lines of
  the text?
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_log_normal}
\alias{distr_to_xml_log_normal}
\title{Internal function}
\usage{
distr_to_xml_log_normal(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{a log-normal distribution,
as created by \code{\link{create_log_normal_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a log-normal distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_mrca_align_ids_in_fastas.R
\name{is_mrca_align_ids_in_fastas}
\alias{is_mrca_align_ids_in_fastas}
\title{Determine if an MRCA prior's alignment IDs are present in the FASTA files}
\usage{
is_mrca_align_ids_in_fastas(mrca_prior, fasta_filenames)
}
\arguments{
\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{fasta_filenames}{One or more FASTA filenames.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}
}
\value{
TRUE if the MRCA prior's alignment IDs
  is present in the FASTA files.
  Returns FALSE otherwise
}
\description{
Determine if an MRCA prior's alignment IDs are present in the FASTA files
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_mrca_prior_with_distr.R
\name{has_mrca_prior_with_distr}
\alias{has_mrca_prior_with_distr}
\title{See if the inference model has one MRCA prior with a distribution}
\usage{
has_mrca_prior_with_distr(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
TRUE if the inference model has one MRCA prior with a distribution,
  FALSE otherwise
}
\description{
See if the inference model has one MRCA prior with a distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_crown_age.R
\name{get_crown_age}
\alias{get_crown_age}
\title{Obtain the crown age of a phylogeny.}
\usage{
get_crown_age(phylogeny)
}
\arguments{
\item{phylogeny}{The phylogeny to obtain the crown age of}
}
\value{
the crown age of the phylogeny
}
\description{
The crown age of a phylogeny is the time between
the present and the moment of at which the first
diversification (resulting in two lineages) happened.
}
\examples{
  phylogeny <- ape::read.tree(text = "(a:15,b:15):1;")
  created <- get_crown_age(phylogeny = phylogeny)
  testit::assert(created == 15)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_screenlog.R
\name{check_screenlog_names}
\alias{check_screenlog_names}
\title{Check if the \code{screenlog} has the list elements
of a valid \code{screenlog} object.}
\usage{
check_screenlog_names(screenlog)
}
\arguments{
\item{screenlog}{a \code{screenlog},
as created by \link{create_screenlog}}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element is missing
}
\seealso{
Use \link{create_screenlog} to create a valid \code{screenlog}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_m_param.R
\name{is_m_param}
\alias{is_m_param}
\title{Determine if the object is a valid
m parameter}
\usage{
is_m_param(m_param)
}
\arguments{
\item{m_param}{an m parameter,
as created by \link{create_m_param}}
}
\value{
TRUE if x is a valid m parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
m parameter
}
\examples{

is_m_param(create_alpha_param())
is_m_param(create_beta_param())
is_m_param(create_clock_rate_param())
is_m_param(create_kappa_1_param())
is_m_param(create_kappa_2_param())
is_m_param(create_lambda_param())
is_m_param(create_m_param())
is_m_param(create_mean_param())
is_m_param(create_mu_param())
is_m_param(create_rate_ac_param())
is_m_param(create_rate_ag_param())
is_m_param(create_rate_at_param())
is_m_param(create_rate_cg_param())
is_m_param(create_rate_ct_param())
is_m_param(create_rate_gt_param())
is_m_param(create_s_param())
is_m_param(create_scale_param())
is_m_param(create_sigma_param())

is_m_param(NA)
is_m_param(NULL)
is_m_param("nonsense")
is_m_param(create_jc69_site_model())
is_m_param(create_strict_clock_model())
is_m_param(create_yule_tree_prior())
is_m_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mcmc.R
\name{create_mcmc}
\alias{create_mcmc}
\title{Create an MCMC configuration.}
\usage{
create_mcmc(
  chain_length = 1e+07,
  store_every = -1,
  pre_burnin = 0,
  n_init_attempts = 10,
  sample_from_prior = FALSE,
  tracelog = beautier::create_tracelog(),
  screenlog = beautier::create_screenlog(),
  treelog = beautier::create_treelog()
)
}
\arguments{
\item{chain_length}{length of the MCMC chain}

\item{store_every}{number of states the MCMC will process
before the posterior's state will be saved to file.
Use -1 or \code{NA} to use the default frequency.}

\item{pre_burnin}{number of burn in samples taken before entering
the main loop}

\item{n_init_attempts}{number of initialization attempts before failing}

\item{sample_from_prior}{set to \link{TRUE} to sample from the prior}

\item{tracelog}{a \code{tracelog},
as created by \link{create_tracelog}}

\item{screenlog}{a \code{screenlog},
as created by \link{create_screenlog}}

\item{treelog}{a \code{treelog},
as created by \link{create_treelog}}
}
\value{
an MCMC configuration
}
\description{
Create an MCMC configuration, as in the BEAUti MCMC tab.
}
\details{
There are four things that can be saved:
 * \code{store_every}: saves the state of the MCMC to file,
   as a \code{.state.xml} file
 * \code{tracelog}: stores the trace of the state of the MCMC
   to file. See \code{create_tracelog}
   how to specify the filename
 * \code{screenlog}: stores the screen output
   to file. See \code{create_screenlog}
   how to specify the filename
 * \code{treelog}: stores the estimated phylogenies
   to file. See \code{create_treelog}
   how to specify the filename
}
\examples{
# Create an MCMC chain with 50 states
mcmc <- create_mcmc(chain_length = 50000, store_every = 1000)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  mcmc = mcmc
)
file.remove(beast2_input_file)
}
\seealso{
Use \code{\link{create_test_mcmc}} to create a short regular MCMC,
that can be used for testing runs.
Use \code{\link{create_ns_mcmc}} to create an MCMC for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/m_param_to_xml.R
\name{m_param_to_xml}
\alias{m_param_to_xml}
\title{Internal function}
\usage{
m_param_to_xml(m_param, beauti_options = create_beauti_options())
}
\arguments{
\item{m_param}{an m parameter,
as created by \link{create_m_param}}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts an m parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_test_mcmc.R
\name{create_test_mcmc}
\alias{create_test_mcmc}
\title{Create an MCMC configuration for testing.}
\usage{
create_test_mcmc(
  chain_length = 3000,
  store_every = 1000,
  pre_burnin = 0,
  n_init_attempts = 10,
  sample_from_prior = FALSE,
  tracelog = beautier::create_test_tracelog(),
  screenlog = beautier::create_test_screenlog(),
  treelog = beautier::create_test_treelog()
)
}
\arguments{
\item{chain_length}{length of the MCMC chain}

\item{store_every}{number of states the MCMC will process
before the posterior's state will be saved to file.
Use -1 or \code{NA} to use the default frequency.}

\item{pre_burnin}{number of burn in samples taken before entering
the main loop}

\item{n_init_attempts}{number of initialization attempts before failing}

\item{sample_from_prior}{set to \link{TRUE} to sample from the prior}

\item{tracelog}{a \code{tracelog},
as created by \link{create_tracelog}}

\item{screenlog}{a \code{screenlog},
as created by \link{create_screenlog}}

\item{treelog}{a \code{treelog},
as created by \link{create_treelog}}
}
\value{
an MCMC configuration
}
\description{
Create an MCMC configuration for testing.
}
\examples{
# Create an MCMC chain with 50 states
mcmc <- create_test_mcmc()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  mcmc = mcmc
)
file.remove(beast2_input_file)
}
\seealso{
Use \code{\link{create_mcmc}} to create a default BEAST2 MCMC
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_tree_prior.R
\name{check_tree_prior}
\alias{check_tree_prior}
\title{Check if the tree prior is a valid tree prior}
\usage{
check_tree_prior(tree_prior)
}
\arguments{
\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if the tree priors are invalid
}
\examples{
check_tree_prior(create_yule_tree_prior())
check_tree_prior(create_bd_tree_prior())
check_tree_prior(create_cbs_tree_prior())
check_tree_prior(create_ccp_tree_prior())
check_tree_prior(create_cep_tree_prior())

# Can use list of one tree prior
check_tree_prior(list(create_yule_tree_prior()))
}
\seealso{
Use \link{create_tree_prior} to create a valid tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tree_prior.R
\name{create_yule_tree_prior}
\alias{create_yule_tree_prior}
\alias{create_tree_prior_yule}
\title{Create a Yule tree prior}
\usage{
create_yule_tree_prior(
  id = NA,
  birth_rate_distr = create_uniform_distr()
)
}
\arguments{
\item{id}{the ID of the alignment}

\item{birth_rate_distr}{the birth rate distribution,
as created by a \code{\link{create_distr}} function}
}
\value{
a Yule tree_prior
}
\description{
Create a Yule tree prior
}
\examples{
 yule_tree_prior <- create_yule_tree_prior()

 beast2_input_file <- get_beautier_tempfilename()
 create_beast2_input_file(
   input_filename = get_fasta_filename(),
   beast2_input_file,
   tree_prior = yule_tree_prior
)
file.remove(beast2_input_file)
}
\seealso{
An alignment ID can be extracted from
  its FASTA filename using \code{\link{get_alignment_id}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_sigma_param}
\alias{create_sigma_param}
\alias{create_param_sigma}
\title{Create a parameter called sigma}
\usage{
create_sigma_param(id = NA, value = 1)
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}
}
\value{
a parameter called sigma
}
\description{
Create a parameter called sigma
}
\note{
this parameter is used in a normal distribution
  (as returned by \code{\link{create_normal_distr}})
}
\examples{
# Create the parameter
sigma_param <- create_sigma_param()

# Use the parameter in a distribution
normal_distr <- create_normal_distr(
  sigma = sigma_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = normal_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_models_to_xml_operators.R
\name{site_models_to_xml_operators}
\alias{site_models_to_xml_operators}
\title{Write the XML \code{operators} section from the site models.}
\usage{
site_models_to_xml_operators(site_models)
}
\arguments{
\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}
}
\value{
lines of XML text
}
\description{
Write the XML \code{operators} section from the site models.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_branch_rate_model_rln_xml.R
\name{create_branch_rate_model_rln_xml}
\alias{create_branch_rate_model_rln_xml}
\title{Internal function}
\usage{
create_branch_rate_model_rln_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
a character vector of XML strings
}
\description{
Internal function to call \link{create_branch_rate_model_xml}
for a relaxed log-normal clock.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_site_model.R
\name{is_site_model}
\alias{is_site_model}
\title{Determine if the object is a valid site_model}
\usage{
is_site_model(x)
}
\arguments{
\item{x}{an object, to be determined if it is a site_model}
}
\value{
TRUE if the site_model is a valid site_model, FALSE otherwise
}
\description{
Determine if the object is a valid site_model
}
\examples{
# TRUE
is_site_model(create_gtr_site_model())
is_site_model(create_hky_site_model())
is_site_model(create_jc69_site_model())
is_site_model(create_tn93_site_model())

# FALSE
is_site_model(NA)
is_site_model(NULL)
is_site_model("nonsense")
is_site_model(create_strict_clock_model())
is_site_model(create_bd_tree_prior())
is_site_model(create_mcmc())
}
\seealso{
A site model can be created using \code{\link{create_site_model}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_inv_gamma_distr}
\alias{create_inv_gamma_distr}
\alias{create_distr_inv_gamma}
\title{Create an inverse-gamma distribution}
\usage{
create_inv_gamma_distr(
  id = NA,
  alpha = 0,
  beta = 1,
  value = NA,
  lower = NA,
  upper = NA
)
}
\arguments{
\item{id}{the distribution's ID}

\item{alpha}{the alpha shape parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_alpha_param}}}

\item{beta}{the beta shape parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_beta_param}}}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
an inverse-gamma distribution
}
\description{
Create an inverse-gamma distribution
}
\examples{
inv_gamma_distr <- create_inv_gamma_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = inv_gamma_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_xml_operators_from_lines.R
\name{extract_xml_operators_from_lines}
\alias{extract_xml_operators_from_lines}
\title{Extract everything between first operators and last operators line}
\usage{
extract_xml_operators_from_lines(lines)
}
\arguments{
\item{lines}{lines of text}
}
\value{
lines of text from the first to and including the last operators line
}
\description{
Extract everything between first operators and last operators line
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_log_sorts.R
\name{get_log_sorts}
\alias{get_log_sorts}
\title{Get the possible log sorts}
\usage{
get_log_sorts()
}
\description{
Get the possible log sorts
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gamma_site_model_n_params.R
\name{get_gamma_site_model_n_params}
\alias{get_gamma_site_model_n_params}
\title{Get the number of distributions a site model has}
\usage{
get_gamma_site_model_n_params(gamma_site_model)
}
\arguments{
\item{gamma_site_model}{a site model's gamma site model,
as returned by \code{\link{create_gamma_site_model}}}
}
\value{
the number of parameters a site model has
}
\description{
Get the number of distributions a site model has
}
\examples{
  testit::assert(
    get_gamma_site_model_n_params(
      create_gamma_site_model(gamma_cat_count = 0)
    ) == 0
  )
  testit::assert(
    get_gamma_site_model_n_params(
      create_gamma_site_model(gamma_cat_count = 1)
    ) == 0
  )
  testit::assert(
    get_gamma_site_model_n_params(
      create_gamma_site_model(
        gamma_cat_count = 2,
        gamma_shape_prior_distr = create_exp_distr()
      )
    ) == 1
  )
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_trailing_spaces.R
\name{count_trailing_spaces}
\alias{count_trailing_spaces}
\title{Count the number of spaces before the first character}
\usage{
count_trailing_spaces(line)
}
\arguments{
\item{line}{line of text}
}
\value{
the number of spaces before the first character
}
\description{
Count the number of spaces before the first character
}
\examples{
# 0
count_trailing_spaces("x")
# 1
count_trailing_spaces(" y")
# 2
count_trailing_spaces("  <")
# 0
count_trailing_spaces("")
# 1
count_trailing_spaces(" ")
# 2
count_trailing_spaces("  ")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_normal_distr}
\alias{create_normal_distr}
\alias{create_distr_normal}
\title{Create an normal distribution}
\usage{
create_normal_distr(
  id = NA,
  mean = 0,
  sigma = 1,
  value = NA,
  lower = NA,
  upper = NA
)
}
\arguments{
\item{id}{the distribution's ID}

\item{mean}{the mean parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_mean_param}}}

\item{sigma}{the sigma parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_sigma_param}}}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
a normal distribution
}
\description{
Create an normal distribution
}
\examples{
normal_distr <- create_normal_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = normal_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beta_parameter_to_xml.R
\name{beta_parameter_to_xml}
\alias{beta_parameter_to_xml}
\title{Internal function}
\usage{
beta_parameter_to_xml(beta_parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{beta_parameter}{a beta parameter,
as created by \link{create_beta_param}}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a beta parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_rate_at_param}
\alias{create_rate_at_param}
\alias{create_param_rate_at}
\title{Create a parameter called 'rate AT'}
\usage{
create_rate_at_param(id = NA, estimate = TRUE, value = "1.0", lower = "0.0")
}
\arguments{
\item{id}{the parameter's ID}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}

\item{value}{value of the parameter}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}
}
\value{
a parameter called 'rate AT'
}
\description{
Create a parameter called 'rate AT'
}
\examples{
# Create parameter
rate_at_param <- create_rate_at_param(value = 1, estimate = FALSE)

# Use the parameter to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  site_model = create_gtr_site_model(
    rate_at_param = rate_at_param
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_rate_ac_param}
\alias{is_rate_ac_param}
\title{Determine if the object is a valid
'rate AC' parameter}
\usage{
is_rate_ac_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
'rate AC' parameter}
}
\value{
TRUE if x is a valid 'rate AC' parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
'rate AC' parameter
}
\examples{

is_rate_ac_param(create_alpha_param())
is_rate_ac_param(create_beta_param())
is_rate_ac_param(create_clock_rate_param())
is_rate_ac_param(create_kappa_1_param())
is_rate_ac_param(create_kappa_2_param())
is_rate_ac_param(create_lambda_param())
is_rate_ac_param(create_m_param())
is_rate_ac_param(create_mean_param())
is_rate_ac_param(create_mu_param())
is_rate_ac_param(create_rate_ac_param())
is_rate_ac_param(create_rate_ag_param())
is_rate_ac_param(create_rate_at_param())
is_rate_ac_param(create_rate_cg_param())
is_rate_ac_param(create_rate_ct_param())
is_rate_ac_param(create_rate_gt_param())
is_rate_ac_param(create_s_param())
is_rate_ac_param(create_scale_param())
is_rate_ac_param(create_sigma_param())

is_rate_ac_param(NA)
is_rate_ac_param(NULL)
is_rate_ac_param("nonsense")
is_rate_ac_param(create_jc69_site_model())
is_rate_ac_param(create_strict_clock_model())
is_rate_ac_param(create_yule_tree_prior())
is_rate_ac_param(create_mcmc())
}
\seealso{
\code{\link{create_rate_ac_param}} creates a 'rate AC' parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_mrca_priors.R
\name{are_mrca_priors}
\alias{are_mrca_priors}
\title{Determine if x consists out of MRCA priors}
\usage{
are_mrca_priors(mrca_priors)
}
\arguments{
\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}
}
\value{
TRUE if x, or all elements of x, are MRCA priors.
  Returns FALSE otherwise
}
\description{
Determine if x consists out of MRCA priors
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_clock_model.R
\name{is_rln_clock_model}
\alias{is_rln_clock_model}
\title{Determine if the object is a valid relaxed log normal clock model}
\usage{
is_rln_clock_model(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
relaxed log normal clock model,
as created by \code{\link{create_rln_clock_model}})}
}
\value{
TRUE if x is a valid relaxed log normal clock model, FALSE otherwise
}
\description{
Determine if the object is a valid relaxed log normal clock model
}
\examples{

is_rln_clock_model(create_strict_clock_model())
is_rln_clock_model(create_rln_clock_model())

is_rln_clock_model(NA)
is_rln_clock_model(NULL)
is_rln_clock_model("nonsense")
is_rln_clock_model(create_jc69_site_model())
is_rln_clock_model(create_mcmc())
}
\seealso{
\code{\link{create_clock_model}} shows an overview of
  functions to create a clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_beta_distr}
\alias{create_beta_distr}
\alias{create_distr_beta}
\title{Create a beta distribution}
\usage{
create_beta_distr(
  id = NA,
  alpha = 0,
  beta = 1,
  value = NA,
  lower = NA,
  upper = NA
)
}
\arguments{
\item{id}{the distribution's ID}

\item{alpha}{the alpha shape parameter,
a numeric value.
The value
of alpha must be at least 0.0.
For advanced usage, use the structure
as returned by \code{\link{create_alpha_param}}.}

\item{beta}{the beta shape parameter,
a numeric value.
The value
of beta must be at least 1.0.
For advanced usage, use the structure
as returned by \code{\link{create_beta_param}}.}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
a beta distribution
}
\description{
Create a beta distribution
}
\examples{
beta_distr <- create_beta_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = beta_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alpha_parameter_to_xml.R
\name{alpha_parameter_to_xml}
\alias{alpha_parameter_to_xml}
\title{Internal function}
\usage{
alpha_parameter_to_xml(
  alpha_parameter,
  beauti_options = create_beauti_options()
)
}
\arguments{
\item{alpha_parameter}{an alpha parameter,
as created by \link{create_alpha_param}}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts an alpha parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_tree_priors.R
\name{init_yule_tree_prior}
\alias{init_yule_tree_prior}
\title{Initializes a Yule tree prior}
\usage{
init_yule_tree_prior(yule_tree_prior, distr_id, param_id)
}
\arguments{
\item{yule_tree_prior}{a Yule tree_prior,
as created by \code{\link{create_yule_tree_prior}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized Yule tree prior
}
\description{
Initializes a Yule tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_site_model.R
\name{check_site_model}
\alias{check_site_model}
\title{Check if the site model is a valid site model}
\usage{
check_site_model(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if the site models are invalid
}
\examples{

check_site_model(create_jc69_site_model())
check_site_model(create_hky_site_model())
check_site_model(create_tn93_site_model())
check_site_model(create_gtr_site_model())

# Can use list of one site model
check_site_model(list(create_jc69_site_model()))
}
\seealso{
Use \link{create_site_model} to create a valid site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_data_xml.R
\name{create_data_xml}
\alias{create_data_xml}
\title{Create the \code{<data ..>} XML}
\usage{
create_data_xml(id, beast2_version)
}
\arguments{
\item{id}{an alignment's IDs.
An ID can be extracted from its FASTA filename
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{beast2_version}{BEAST2 version, for example, code{"2.5"}}
}
\value{
lines of XML text
}
\description{
Create the \code{<data ..>} XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_rate_cg_param}
\alias{is_rate_cg_param}
\title{Determine if the object is a valid
'rate CG' parameter}
\usage{
is_rate_cg_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
'rate CG' parameter}
}
\value{
TRUE if x is a valid 'rate CG' parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
'rate CG' parameter
}
\examples{

is_rate_cg_param(create_alpha_param())
is_rate_cg_param(create_beta_param())
is_rate_cg_param(create_clock_rate_param())
is_rate_cg_param(create_kappa_1_param())
is_rate_cg_param(create_kappa_2_param())
is_rate_cg_param(create_lambda_param())
is_rate_cg_param(create_m_param())
is_rate_cg_param(create_mean_param())
is_rate_cg_param(create_mu_param())
is_rate_cg_param(create_rate_ac_param())
is_rate_cg_param(create_rate_ag_param())
is_rate_cg_param(create_rate_at_param())
is_rate_cg_param(create_rate_cg_param())
is_rate_cg_param(create_rate_ct_param())
is_rate_cg_param(create_rate_gt_param())
is_rate_cg_param(create_s_param())
is_rate_cg_param(create_scale_param())
is_rate_cg_param(create_sigma_param())

is_rate_cg_param(NA)
is_rate_cg_param(NULL)
is_rate_cg_param("nonsense")
is_rate_cg_param(create_jc69_site_model())
is_rate_cg_param(create_strict_clock_model())
is_rate_cg_param(create_yule_tree_prior())
is_rate_cg_param(create_mcmc())
}
\seealso{
\code{\link{create_rate_cg_param}} creates a 'rate CG' parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_param_names.R
\name{get_param_names}
\alias{get_param_names}
\title{Get the parameter names}
\usage{
get_param_names()
}
\value{
the parameter names
}
\description{
Get the parameter names
}
\examples{
names <- get_param_names()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma_site_model_to_xml_prior_distr.R
\name{gamma_site_model_to_xml_prior_distr}
\alias{gamma_site_model_to_xml_prior_distr}
\title{Creates the gamma site models section in the distribution section
of a BEAST2 XML parameter file}
\usage{
gamma_site_model_to_xml_prior_distr(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
lines of XML text
}
\description{
Creates the gamma site models section in the distribution section
of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rln_clock_model_to_xml_mean_rate_prior.R
\name{rln_clock_model_to_xml_mean_rate_prior}
\alias{rln_clock_model_to_xml_mean_rate_prior}
\title{Used by \code{\link{clock_models_to_xml_prior_distr}}}
\usage{
rln_clock_model_to_xml_mean_rate_prior(rln_clock_model)
}
\arguments{
\item{rln_clock_model}{a Relaxed Log-Normal clock model,
as returned by \code{\link{create_rln_clock_model}}}
}
\value{
lines of XML text
}
\description{
Used by \code{\link{clock_models_to_xml_prior_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_kappa_1_param}
\alias{create_kappa_1_param}
\alias{create_param_kappa_1}
\title{Create a parameter called kappa 1}
\usage{
create_kappa_1_param(id = NA, lower = "0.0", value = "2.0", estimate = TRUE)
}
\arguments{
\item{id}{the parameter's ID}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}

\item{value}{value of the parameter}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}
}
\value{
a parameter called kappa 1
}
\description{
Create a parameter called kappa 1
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_distr.R
\name{create_beast2_input_distr_prior}
\alias{create_beast2_input_distr_prior}
\title{Creates the prior section in the distribution section
of a BEAST2 XML parameter file}
\usage{
create_beast2_input_distr_prior(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
Creates the prior section in the distribution section
of a BEAST2 XML parameter file
}
\note{
this function is not intended for regular use, thus its
  long name length is accepted
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\seealso{
this function is called by \code{create_beast2_input_distr},
  together with \code{create_beast2_input_distr_lh}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_rate_at}
\alias{parameter_to_xml_rate_at}
\title{Internal function}
\usage{
parameter_to_xml_rate_at(
  parameter,
  beauti_options = create_beauti_options(),
  which_name = "state_node"
)
}
\arguments{
\item{parameter}{a 'rate AT' parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_rate_at_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{which_name}{the name, can be \code{state_node} or \code{rate_name}}
}
\value{
the parameter as XML text
}
\description{
Converts a 'rate AT' parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_model_to_tracelog_xml.R
\name{tree_model_to_tracelog_xml}
\alias{tree_model_to_tracelog_xml}
\title{Internal function}
\usage{
tree_model_to_tracelog_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Creates the tree models' XML for the tracelog section.
That is, all XML tags that have the word 'tree' in them.
}
\note{
use site_models just because it contains all IDs
}
\examples{
# <logger id="tracelog" ...>
#'   # Here
# </logger>
}
\seealso{
the complete tracelog section is created
  by \code{\link{create_tracelog_xml}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_distr.R
\name{bd_tree_prior_to_xml_prior_distr}
\alias{bd_tree_prior_to_xml_prior_distr}
\title{Creates the tree prior section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a Birth-Death tree prior}
\usage{
bd_tree_prior_to_xml_prior_distr(bd_tree_prior)
}
\arguments{
\item{bd_tree_prior}{a Birth-Death tree prior, as created
by \code{\link{create_bd_tree_prior}}}
}
\description{
Creates the tree prior section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a Birth-Death tree prior
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_data.R
\name{create_beast2_input_data}
\alias{create_beast2_input_data}
\title{Creates the \code{data} section of a BEAST2 XML parameter file}
\usage{
create_beast2_input_data(
  input_filename,
  beauti_options = beautier::create_beauti_options()
)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
lines of XML text
}
\description{
Creates the \code{data} section of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_inv_gamma_distr}
\alias{init_inv_gamma_distr}
\title{Initializes an inverse gamma distribution}
\usage{
init_inv_gamma_distr(inv_gamma_distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{inv_gamma_distr}{an inverse gamma distribution,
using \link{create_inv_gamma_distr}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized inverse gamma distribution
}
\description{
Initializes an inverse gamma distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_s_param}
\alias{create_s_param}
\alias{create_param_s}
\title{Create a parameter called s}
\usage{
create_s_param(id = NA, value = 0, lower = 0, upper = Inf)
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}

\item{upper}{upper value of the parameter}
}
\value{
a parameter called s
}
\description{
Create a parameter called s
}
\note{
this parameter is used in a log-normal distribution
  (as returned by \code{\link{create_log_normal_distr}})
}
\examples{
# Create the parameter
s_param <- create_s_param()

# Use the parameter in a distribution
log_normal_distr <- create_log_normal_distr(
  s = s_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = log_normal_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_prior_to_xml_prior_distr.R
\name{tree_prior_to_xml_prior_distr}
\alias{tree_prior_to_xml_prior_distr}
\title{Creates the distribution section in the prior section of the
distribution section of a BEAST2 XML parameter file.}
\usage{
tree_prior_to_xml_prior_distr(tree_prior)
}
\arguments{
\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
lines of XML text
}
\description{
These lines start with '<distribution id='
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_exp_distr}
\alias{is_exp_distr}
\title{Determine if the object is a valid
exponential distribution
as created by \code{\link{create_exp_distr}}}
\usage{
is_exp_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
exponential distribution}
}
\value{
TRUE if x is a valid exponential distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
exponential distribution
as created by \code{\link{create_exp_distr}}
}
\examples{
# TRUE
is_exp_distr(create_exp_distr())
# FALSE
is_exp_distr(create_gamma_distr())
is_exp_distr(NA)
is_exp_distr(NULL)
is_exp_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_phylo.R
\name{is_phylo}
\alias{is_phylo}
\title{Checks if the input is a phylogeny}
\usage{
is_phylo(x)
}
\arguments{
\item{x}{input to be checked}
}
\value{
TRUE or FALSE
}
\description{
Checks if the input is a phylogeny
}
\examples{
  phylogeny <- ape::read.tree(text = "(a:15,b:15):1;")
  testit::assert(is_phylo(phylogeny))

  testit::assert(!is_phylo("nonsense"))
  testit::assert(!is_phylo(NA))
  testit::assert(!is_phylo(NULL))
}
\seealso{
Use \link{check_phylogeny} to check for a phylogeny
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_temp_treelog_filename.R
\name{create_temp_treelog_filename}
\alias{create_temp_treelog_filename}
\title{Create a filename for a temporary treelog file}
\usage{
create_temp_treelog_filename()
}
\description{
Create a filename for a temporary treelog file
}
\seealso{
use \link{create_treelog} to create a treelog.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_strict_clock_rate_scaler_operator_xml.R
\name{create_strict_clock_rate_scaler_operator_xml}
\alias{create_strict_clock_rate_scaler_operator_xml}
\title{Internal function}
\usage{
create_strict_clock_rate_scaler_operator_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the following XML:
\code{
  ...
}
}
\description{
Creates the \code{StrictClockRateScaler} operator
such as:
\code{
  ...
}
}
\examples{
create_strict_clock_rate_scaler_operator_xml(
  create_inference_model(
    clock_model = create_strict_clock_model(id = 314)
  )
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mrca_prior_taxa_names.R
\name{check_mrca_prior_taxa_names}
\alias{check_mrca_prior_taxa_names}
\title{Check the MRCA prior's taxon names are valid.}
\usage{
check_mrca_prior_taxa_names(taxa_names)
}
\arguments{
\item{taxa_names}{names of the taxa,
as returned by \code{\link{get_taxa_names}}.
Keep at \code{NA} to have it initialized automatically,
using all taxa in the alignment}
}
\description{
Will \link{stop} if not.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_ucld_mean_state_node_param_xml.R
\name{create_ucld_mean_state_node_param_xml}
\alias{create_ucld_mean_state_node_param_xml}
\title{Internal function}
\usage{
create_ucld_mean_state_node_param_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the XML
\code{
  <parameter id=\"ucldMean.c:[id]\" spec=\"parameter.RealParameter\"
    name=\"stateNode\">1.0</parameter>
}
}
\description{
Creates the \code{ucldMean.c} parameter with the name \code{stateNode},
such as:
\code{
  <parameter id=\"ucldMean.c:[id]\" spec=\"parameter.RealParameter\"
    name=\"stateNode\">1.0</parameter>
}
}
\examples{
create_ucld_mean_state_node_param_xml(
  create_inference_model(
    clock_model = create_rln_clock_model(id = 314),
    beauti_options = create_beauti_options_v2_6()
  )
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_model_to_xml_state.R
\name{clock_model_to_xml_state}
\alias{clock_model_to_xml_state}
\title{Internal function}
\usage{
clock_model_to_xml_state(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text, without indentation nor \code{state}
  tags
}
\description{
Converts a clock model to the \code{state} section of the
XML as text
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_tree_prior.R
\name{is_init_cep_tree_prior}
\alias{is_init_cep_tree_prior}
\title{Determine if x is an initialized Coalescent Exponential Population
  tree_prior object}
\usage{
is_init_cep_tree_prior(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized Coalescent Exponential Population tree prior object}
}
\value{
TRUE if x is an initialized Coalescent Exponential Population
  tree prior object
}
\description{
Determine if x is an initialized Coalescent Exponential Population
  tree_prior object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_gamma_site_model.R
\name{check_gamma_site_model}
\alias{check_gamma_site_model}
\title{Checks if the parameter is a valid gamma site model}
\usage{
check_gamma_site_model(gamma_site_model)
}
\arguments{
\item{gamma_site_model}{a site model's gamma site model,
as returned by \code{\link{create_gamma_site_model}}}
}
\value{
nothing. Will call \code{stop} if the argument is not a valid
  gamma site model
}
\description{
Checks if the parameter is a valid gamma site model
}
\examples{
check_gamma_site_model(create_gamma_site_model())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_beauti_options.R
\name{check_beauti_options}
\alias{check_beauti_options}
\title{Check if the \code{beauti_options} is a valid \code{beauti_options} object.}
\usage{
check_beauti_options(beauti_options)
}
\arguments{
\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if the \code{beauti_options} object is invalid
}
\examples{
check_beauti_options(create_beauti_options())
}
\seealso{
Use \link{create_beauti_options} to create a valid
  BEAUti options setup
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_clock_models.R
\name{are_clock_models}
\alias{are_clock_models}
\title{Determine if x consists out of clock_models objects}
\usage{
are_clock_models(x)
}
\arguments{
\item{x}{the object to check if it consists out of clock_models objects}
}
\value{
TRUE if x, or all elements of x, are clock_model objects
}
\description{
Determine if x consists out of clock_models objects
}
\examples{

rln_clock_model <- create_rln_clock_model()
strict_clock_model <- create_strict_clock_model()
both_clock_models <- list(rln_clock_model, strict_clock_model)
# TRUE
are_clock_models(rln_clock_model)
are_clock_models(strict_clock_model)
are_clock_models(both_clock_models)

# FALSE
are_clock_models(NA)
are_clock_models(NULL)
are_clock_models("nonsense")
are_clock_models(create_jc69_site_model())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_inv_gamma_distr}
\alias{is_inv_gamma_distr}
\title{Determine if the object is a valid
inverse-gamma distribution
as created by \code{\link{create_inv_gamma_distr}}}
\usage{
is_inv_gamma_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
inverse-gamma distribution}
}
\value{
TRUE if x is a valid inverse-gamma distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
inverse-gamma distribution
as created by \code{\link{create_inv_gamma_distr}}
}
\examples{
# TRUE
is_inv_gamma_distr(create_inv_gamma_distr())
# FALSE
is_inv_gamma_distr(create_laplace_distr())
is_inv_gamma_distr(NA)
is_inv_gamma_distr(NULL)
is_inv_gamma_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_model.R
\name{create_gtr_site_model}
\alias{create_gtr_site_model}
\alias{create_site_model_gtr}
\title{Create a GTR site model}
\usage{
create_gtr_site_model(
  id = NA,
  gamma_site_model = create_gamma_site_model(),
  rate_ac_prior_distr = create_gamma_distr(alpha = 0.05, beta = create_beta_param(value
    = "10.0")),
  rate_ag_prior_distr = create_gamma_distr(alpha = 0.05, beta = create_beta_param(value
    = "20.0")),
  rate_at_prior_distr = create_gamma_distr(alpha = 0.05, beta = create_beta_param(value
    = "10.0")),
  rate_cg_prior_distr = create_gamma_distr(alpha = 0.05, beta = create_beta_param(value
    = "10.0")),
  rate_gt_prior_distr = create_gamma_distr(alpha = 0.05, beta = create_beta_param(value
    = "10.0")),
  rate_ac_param = create_rate_ac_param(),
  rate_ag_param = create_rate_ag_param(),
  rate_at_param = create_rate_at_param(),
  rate_cg_param = create_rate_cg_param(),
  rate_ct_param = create_rate_ct_param(),
  rate_gt_param = create_rate_gt_param(),
  freq_equilibrium = "estimated"
)
}
\arguments{
\item{id}{the IDs of the alignment (can be extracted from
the FASTA filename using \code{\link{get_alignment_id}})}

\item{gamma_site_model}{a gamma site model, as created
by \code{\link{create_gamma_site_model}}}

\item{rate_ac_prior_distr}{the AC rate prior distribution,
as returned by \code{\link{create_distr}}}

\item{rate_ag_prior_distr}{the AG rate prior distribution,
as returned by \code{\link{create_distr}}}

\item{rate_at_prior_distr}{the AT rate prior distribution,
as returned by \code{\link{create_distr}}}

\item{rate_cg_prior_distr}{the CG rate prior distribution,
as returned by \code{\link{create_distr}}}

\item{rate_gt_prior_distr}{the GT rate prior distribution,
as returned by \code{\link{create_distr}}}

\item{rate_ac_param}{the 'rate AC' parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_rate_ac_param}}}

\item{rate_ag_param}{the 'rate AG' parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_rate_ag_param}}}

\item{rate_at_param}{the 'rate AT' parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_rate_at_param}}}

\item{rate_cg_param}{the 'rate CG' parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_rate_cg_param}}}

\item{rate_ct_param}{the 'rate CT' parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_rate_ct_param}}}

\item{rate_gt_param}{the 'rate GT' parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_rate_gt_param}}}

\item{freq_equilibrium}{the frequency in which the rates are at equilibrium
are either \code{estimated}, \code{empirical} or \code{all_equal}.
\code{get_freq_equilibrium_names} returns the possible values
for \code{freq_equilibrium}}
}
\value{
a GTR site_model
}
\description{
Create a GTR site model
}
\examples{
gtr_site_model <- create_gtr_site_model(
  rate_ac_param = 1.2,
  rate_ag_param = 2.3,
  rate_at_param = 3.4,
  rate_cg_param = 4.5,
  rate_ct_param = 5.6,
  rate_gt_param = 6.7
)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  site_model = gtr_site_model
)
file.remove(beast2_input_file)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_first_regex_line.R
\name{find_first_regex_line}
\alias{find_first_regex_line}
\title{Find the first line that satisfies a regex}
\usage{
find_first_regex_line(lines, regex)
}
\arguments{
\item{lines}{lines of text}

\item{regex}{the regex as text}
}
\value{
index of the line
}
\description{
Find the first line that satisfies a regex
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_gamma_distr}
\alias{create_gamma_distr}
\alias{create_distr_gamma}
\title{Create a gamma distribution}
\usage{
create_gamma_distr(
  id = NA,
  alpha = 0.5396,
  beta = 0.3819,
  value = NA,
  lower = NA,
  upper = NA
)
}
\arguments{
\item{id}{the distribution's ID}

\item{alpha}{the alpha shape parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_alpha_param}}}

\item{beta}{the beta shape parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_beta_param}}}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
a gamma distribution
}
\description{
Create a gamma distribution
}
\examples{
gamma_distr <- create_gamma_distr(
   alpha = 0.05,
   beta = 10.0
)

gtr_site_model <- create_gtr_site_model(
  rate_ac_prior_distr = gamma_distr
)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  site_model = gtr_site_model
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_model_to_xml_treelogger.R
\name{clock_model_to_xml_treelogger}
\alias{clock_model_to_xml_treelogger}
\title{Convert a clock model to the XML of the \code{TreeLogger}}
\usage{
clock_model_to_xml_treelogger(clock_model)
}
\arguments{
\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}
}
\value{
a character vector of XML strings
}
\description{
Convert a clock model to the XML of the \code{TreeLogger}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beauti_options.R
\name{create_beauti_options}
\alias{create_beauti_options}
\title{Function to create a set of BEAUti options.}
\usage{
create_beauti_options(
  capitalize_first_char_id = FALSE,
  nucleotides_uppercase = FALSE,
  beast2_version = "2.4",
  required = "",
  sequence_indent = 20
)
}
\arguments{
\item{capitalize_first_char_id}{must the ID of alignment start with a
capital? TRUE if yes, FALSE if it can be left lower case (if it is
lowercase)}

\item{nucleotides_uppercase}{must the nucleotides of the DNA sequence be
in uppercase?}

\item{beast2_version}{the BEAST2 version}

\item{required}{things that may be required,
for example \code{BEAST v2.5.0}}

\item{sequence_indent}{the number of spaces the XML \code{sequence}
lines are indented}
}
\value{
a BEAUti options structure
}
\description{
BEAUti options are settings that differ between BEAUti
version. The use of these options is mostly for testing
older versions
Whatever option chosen here, the created XML file will be valid.
}
\details{
Available BEAUti options are:\cr
\itemize{
  \item \link{create_beauti_options_v2_4}
  \item \link{create_beauti_options_v2_6}
}
}
\examples{
beauti_options <- create_beauti_options_v2_4()
xml <- create_beast2_input(
  get_fasta_filename(),
  beauti_options = beauti_options
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_clock_models.R
\name{create_clock_models}
\alias{create_clock_models}
\title{Creates all supported clock models,
  which is a list of the types returned by
  \code{\link{create_rln_clock_model}},
  and \code{\link{create_strict_clock_model}}}
\usage{
create_clock_models()
}
\value{
a list of site_models
}
\description{
Creates all supported clock models,
  which is a list of the types returned by
  \code{\link{create_rln_clock_model}},
  and \code{\link{create_strict_clock_model}}
}
\examples{
clock_models <- create_clock_models()
is_rln_clock_model(clock_models[[1]])
is_strict_clock_model(clock_models[[2]])
}
\seealso{
Use \link{create_clock_model} to create a clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_model.R
\name{create_tn93_site_model}
\alias{create_tn93_site_model}
\alias{create_site_model_tn93}
\title{Create a TN93 site model}
\usage{
create_tn93_site_model(
  id = NA,
  gamma_site_model = create_gamma_site_model(),
  kappa_1_param = create_kappa_1_param(),
  kappa_2_param = create_kappa_2_param(),
  kappa_1_prior_distr = create_log_normal_distr(m = 1, s = 1.25),
  kappa_2_prior_distr = create_log_normal_distr(m = 1, s = 1.25),
  freq_equilibrium = "estimated"
)
}
\arguments{
\item{id}{the IDs of the alignment (can be extracted from
the FASTA filename using \code{\link{get_alignment_id}})}

\item{gamma_site_model}{a gamma site model, as created
by \code{\link{create_gamma_site_model}}}

\item{kappa_1_param}{the 'kappa 1' parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_kappa_1_param}}}

\item{kappa_2_param}{the 'kappa 2' parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_kappa_2_param}}}

\item{kappa_1_prior_distr}{the distribution of the kappa 1 prior,
which is a log-normal distribution
(as created by \code{\link{create_log_normal_distr}})
by default}

\item{kappa_2_prior_distr}{the distribution of the kappa 2 prior,
which is a log-normal distribution
(as created by \code{\link{create_log_normal_distr}})
by default}

\item{freq_equilibrium}{the frequency in which the rates are at equilibrium
are either \code{estimated}, \code{empirical} or \code{all_equal}.
\code{get_freq_equilibrium_names} returns the possible values
for \code{freq_equilibrium}}
}
\value{
a TN93 site_model
}
\description{
Create a TN93 site model
}
\examples{
 tn93_site_model <- create_tn93_site_model(
   kappa_1_param = 2.0,
   kappa_2_param = 2.0
 )

 output_filename <- get_beautier_tempfilename()
 create_beast2_input_file(
   input_filename = get_fasta_filename(),
   output_filename = output_filename,
   site_model = tn93_site_model
 )
file.remove(output_filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_kappa_2_param}
\alias{is_kappa_2_param}
\title{Determine if the object is a valid kappa 2 parameter}
\usage{
is_kappa_2_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
kappa 2 parameter}
}
\value{
TRUE if x is a valid kappa_2 parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid kappa 2 parameter
}
\examples{

is_kappa_2_param(create_alpha_param())
is_kappa_2_param(create_beta_param())
is_kappa_2_param(create_clock_rate_param())
is_kappa_2_param(create_kappa_1_param())
is_kappa_2_param(create_kappa_2_param())
is_kappa_2_param(create_lambda_param())
is_kappa_2_param(create_m_param())
is_kappa_2_param(create_mean_param())
is_kappa_2_param(create_mu_param())
is_kappa_2_param(create_rate_ac_param())
is_kappa_2_param(create_rate_ag_param())
is_kappa_2_param(create_rate_at_param())
is_kappa_2_param(create_rate_cg_param())
is_kappa_2_param(create_rate_ct_param())
is_kappa_2_param(create_rate_gt_param())
is_kappa_2_param(create_s_param())
is_kappa_2_param(create_scale_param())
is_kappa_2_param(create_sigma_param())

is_kappa_2_param(NA)
is_kappa_2_param(NULL)
is_kappa_2_param("nonsense")
is_kappa_2_param(create_jc69_site_model())
is_kappa_2_param(create_strict_clock_model())
is_kappa_2_param(create_yule_tree_prior())
is_kappa_2_param(create_mcmc())
}
\seealso{
kappa 2 parameters are returned by
  \code{\link{create_kappa_2_param}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_file_exists.R
\name{check_file_exists}
\alias{check_file_exists}
\title{Function to check if a file exists.
Calls \code{stop} if the file is absent}
\usage{
check_file_exists(filename, filename_description = NA)
}
\arguments{
\item{filename}{name of the file}

\item{filename_description}{description of the filename}
}
\value{
nothing. Will \code{stop} if the file is absent,
  with a proper error message
}
\description{
Function to check if a file exists.
Calls \code{stop} if the file is absent
}
\examples{
check_file_exists(get_beautier_path("anthus_aco_sub.fas"))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_distr}
\alias{is_distr}
\title{Determine if the object is a valid distribution}
\usage{
is_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
distribution}
}
\value{
TRUE if x is a valid distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid distribution
}
\examples{
# TRUE
is_distr(create_beta_distr())
is_distr(create_exp_distr())
is_distr(create_gamma_distr())
is_distr(create_inv_gamma_distr())
is_distr(create_laplace_distr())
is_distr(create_log_normal_distr())
is_distr(create_normal_distr())
is_distr(create_one_div_x_distr())
is_distr(create_poisson_distr())
is_distr(create_uniform_distr())

# FALSE
is_distr(NA)
is_distr(NULL)
is_distr("nonsense")
}
\seealso{
use
 \code{\link{is_beta_distr}},
 \code{\link{is_exp_distr}},
 \code{\link{is_gamma_distr}},
 \code{\link{is_inv_gamma_distr}},
 \code{\link{is_laplace_distr}},
 \code{\link{is_log_normal_distr}},
 \code{\link{is_normal_distr}},
 \code{\link{is_one_div_x_distr}},
 \code{\link{is_poisson_distr}},
 or \code{\link{is_uniform_distr}},
 to check for more specific distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_treelog.R
\name{check_treelog_values}
\alias{check_treelog_values}
\title{Check if the treelog has the list elements with valid values
for being a valid treelog object.}
\usage{
check_treelog_values(treelog)
}
\arguments{
\item{treelog}{a \code{treelog},
as created by \link{create_treelog}}
}
\value{
nothing
}
\description{
Calls \code{stop} if a value is invalid
}
\seealso{
Use \link{create_treelog} to create a valid treelog
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_alignment_ids.R
\name{get_alignment_ids_from_fasta_filenames}
\alias{get_alignment_ids_from_fasta_filenames}
\title{Get the alignment ID from one or more FASTA filenames.}
\usage{
get_alignment_ids_from_fasta_filenames(fasta_filenames)
}
\arguments{
\item{fasta_filenames}{One or more FASTA filenames.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}
}
\value{
the IDs from one or more FASTA files
}
\description{
This is done in the same way as BEAST2 does by default.
The files are assumed to be FASTA. If this is not the case, there
may be any kind of error message when calling this function.
}
\examples{
  created <- get_alignment_ids_from_fasta_filenames(
    get_beautier_paths(c("anthus_aco.fas", "anthus_nd2.fas"))
  )
  expected <- c(
    get_alignment_id(get_beautier_path("anthus_aco.fas")),
    get_alignment_id(get_beautier_path("anthus_nd2.fas"))
  )
  testit::assert(created == expected)
}
\seealso{
Use \link{get_alignment_ids} to get the alignment IDs from multiple
kids of files.
Use \link{are_fasta_filenames} to
see if the filenames all have a common FASTA filename extension.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_clock_model_names.R
\name{get_clock_model_names}
\alias{get_clock_model_names}
\title{Get the clock model names}
\usage{
get_clock_model_names()
}
\value{
the clock model names
}
\description{
Get the clock model names
}
\examples{
names <- get_clock_model_names()
}
\seealso{
Use \link{create_clock_models} to create a list
  with all clock models
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_site_model.R
\name{check_site_model_names}
\alias{check_site_model_names}
\title{Check if the \code{site_model} has the list elements
of a valid \code{site_model} object.}
\usage{
check_site_model_names(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element is missing
}
\seealso{
Use \link{create_site_model} to create a valid \code{site_model}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_prior_to_xml_tracelog.R
\name{tree_prior_to_xml_tracelog}
\alias{tree_prior_to_xml_tracelog}
\title{Creates the tree prior's XML for the tracelog section}
\usage{
tree_prior_to_xml_tracelog(tree_prior)
}
\arguments{
\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
lines of XML text
}
\description{
Creates the tree prior's XML for the tracelog section
}
\examples{
# <logger id="tracelog" ...>
#'   # Here
# </logger>
}
\seealso{
all tree priors' tracelog section is created
  by \code{\link{tree_priors_to_xml_tracelog}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_test_treelog.R
\name{create_test_treelog}
\alias{create_test_treelog}
\title{Create a \code{treelog} object}
\usage{
create_test_treelog(
  filename = create_temp_treelog_filename(),
  log_every = 1000,
  mode = "tree",
  sanitise_headers = FALSE,
  sort = "none"
)
}
\arguments{
\item{filename}{name of the file to store the posterior trees}

\item{log_every}{number of MCMC states between writing to file}

\item{mode}{mode how to log.
Valid values are the ones returned by \link{get_log_modes}}

\item{sanitise_headers}{set to \link{TRUE} to sanitise the headers of the
log file}

\item{sort}{how to sort the log.
Valid values are the ones returned by \link{get_log_sorts}}
}
\description{
Create a \code{treelog} object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_mean_param}
\alias{is_mean_param}
\title{Determine if the object is a valid mean parameter}
\usage{
is_mean_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid mean parameter,
as created by \code{\link{create_mean_param}})}
}
\value{
TRUE if x is a valid mean parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid mean parameter
}
\examples{

is_mean_param(create_alpha_param())
is_mean_param(create_beta_param())
is_mean_param(create_clock_rate_param())
is_mean_param(create_kappa_1_param())
is_mean_param(create_kappa_2_param())
is_mean_param(create_lambda_param())
is_mean_param(create_m_param())
is_mean_param(create_mean_param())
is_mean_param(create_mu_param())
is_mean_param(create_rate_ac_param())
is_mean_param(create_rate_ag_param())
is_mean_param(create_rate_at_param())
is_mean_param(create_rate_cg_param())
is_mean_param(create_rate_ct_param())
is_mean_param(create_rate_gt_param())
is_mean_param(create_s_param())
is_mean_param(create_scale_param())
is_mean_param(create_sigma_param())

is_mean_param(NA)
is_mean_param(NULL)
is_mean_param("nonsense")
is_mean_param(create_jc69_site_model())
is_mean_param(create_strict_clock_model())
is_mean_param(create_yule_tree_prior())
is_mean_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_file_base_sans_ext.R
\name{get_file_base_sans_ext}
\alias{get_file_base_sans_ext}
\title{Get the base of the filename base without extension}
\usage{
get_file_base_sans_ext(filename)
}
\arguments{
\item{filename}{A filename}
}
\value{
That filename without its full path and extension
}
\description{
The path need not to actually exist
}
\examples{
# Path need not exist, use UNIX path as example
# test
get_file_base_sans_ext("/home/homer/test.txt")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_site_models.R
\name{are_site_models}
\alias{are_site_models}
\title{Determine if x consists out of site_models objects}
\usage{
are_site_models(x)
}
\arguments{
\item{x}{the object to check if it consists out of site_models objects}
}
\value{
TRUE if x, or all elements of x, are site_model objects
}
\description{
Determine if x consists out of site_models objects
}
\examples{
  jc69_site_model <- create_jc69_site_model()
  gtr_site_model <- create_gtr_site_model()
  both_site_models <- list(jc69_site_model, gtr_site_model)
  testit::assert(are_site_models(jc69_site_model))
  testit::assert(are_site_models(gtr_site_model))
  testit::assert(are_site_models(both_site_models))
}
\seealso{
Use \link{create_site_model} to create a site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_exp_distr}
\alias{is_init_exp_distr}
\title{Determine if x is an initialized exponential distribution object
  as created by \code{\link{create_exp_distr}}}
\usage{
is_init_exp_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized exponential distribution object}
}
\value{
TRUE if x is an initialized exponential distribution object
}
\description{
Determine if x is an initialized exponential distribution object
  as created by \code{\link{create_exp_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_one_double.R
\name{is_one_double}
\alias{is_one_double}
\title{Determines if the argument is a double}
\usage{
is_one_double(x)
}
\arguments{
\item{x}{the object to be determined of if it is one double}
}
\description{
Determines if the argument is a double
}
\examples{
# TRUE
is_one_double(314)
is_one_double(0)
is_one_double(-314)
is_one_double(3.14)

# FALSE
is_one_double(NULL)
is_one_double(NA)
is_one_double(Inf)
is_one_double("nonsense")
is_one_double(is_one_double)
is_one_double(c())
is_one_double(c(1, 2))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_param}
\alias{is_param}
\title{Determine if the object is a valid parameter}
\usage{
is_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid parameter,
as created by \code{\link{create_param}})}
}
\value{
TRUE if x is a valid parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid parameter
}
\examples{
# TRUE
is_param(create_alpha_param())
is_param(create_beta_param())
is_param(create_clock_rate_param())
is_param(create_kappa_1_param())
is_param(create_kappa_2_param())
is_param(create_lambda_param())
is_param(create_m_param())
is_param(create_mean_param())
is_param(create_mu_param())
is_param(create_rate_ac_param())
is_param(create_rate_ag_param())
is_param(create_rate_at_param())
is_param(create_rate_cg_param())
is_param(create_rate_ct_param())
is_param(create_rate_gt_param())
is_param(create_s_param())
is_param(create_scale_param())
is_param(create_sigma_param())

# FALSE
is_param(NA)
is_param(NULL)
is_param("nonsense")
is_param(create_jc69_site_model())
is_param(create_strict_clock_model())
is_param(create_yule_tree_prior())
is_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_state.R
\name{create_beast2_input_state}
\alias{create_beast2_input_state}
\title{Creates the '\code{state}' section of a BEAST2 XML parameter file}
\usage{
create_beast2_input_state(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Creates the '\code{state}' section of a BEAST2 XML parameter file,
without being indented.
}
\details{
The \code{state} tag has these elements:
\preformatted{
   <state[...]>
       <tree[...]>
       [...]
       </tree>
       [parameters]
    </run>
}
}
\seealso{
Use \link{create_beast2_input_state}
to create the XML text of the \code{tree} tag.
to create the XML text of the \code{[parameters]} section.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_tree_prior.R
\name{is_init_bd_tree_prior}
\alias{is_init_bd_tree_prior}
\title{Determine if x is an initialized Birth-Death tree_prior object}
\usage{
is_init_bd_tree_prior(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized Birth-Death tree prior object}
}
\value{
TRUE if x is an initialized Birth-Death tree_prior object
}
\description{
Determine if x is an initialized Birth-Death tree_prior object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_tree_prior.R
\name{is_ccp_tree_prior}
\alias{is_ccp_tree_prior}
\title{Determine if the object is a valid
  constant coalescence population tree prior}
\usage{
is_ccp_tree_prior(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
constant coalescence population tree prior}
}
\value{
TRUE if x is a valid constant coalescence population tree prior,
  FALSE otherwise
}
\description{
Determine if the object is a valid
  constant coalescence population tree prior
}
\examples{
  testit::assert(!is_ccp_tree_prior(create_bd_tree_prior()))
  testit::assert(!is_ccp_tree_prior(create_cbs_tree_prior()))
  testit::assert( is_ccp_tree_prior(create_ccp_tree_prior()))
  testit::assert(!is_ccp_tree_prior(create_cep_tree_prior()))
  testit::assert(!is_ccp_tree_prior(create_yule_tree_prior()))
}
\seealso{
Use \code{\link{create_ccp_tree_prior}} to create a valid
  constant coalescence population tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_mrca_prior.R
\name{is_init_mrca_prior}
\alias{is_init_mrca_prior}
\title{Determine if x is an initialized MRCA prior}
\usage{
is_init_mrca_prior(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized MRCA prior}
}
\value{
TRUE if x is an initialized MRCA prior
}
\description{
Determine if x is an initialized MRCA prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_priors_to_xml_distr.R
\name{tree_priors_to_xml_prior_distr}
\alias{tree_priors_to_xml_prior_distr}
\title{Creates the distribution section in the prior section of the
distribution section of a BEAST2 XML parameter file.}
\usage{
tree_priors_to_xml_prior_distr(tree_priors)
}
\arguments{
\item{tree_priors}{one or more tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
lines of XML text
}
\description{
These lines start with '<distribution id='
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_distr.R
\name{cbs_tree_prior_to_xml_prior_distr}
\alias{cbs_tree_prior_to_xml_prior_distr}
\title{Creates the tree prior section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a Birth-Death tree prior}
\usage{
cbs_tree_prior_to_xml_prior_distr(cbs_tree_prior)
}
\arguments{
\item{cbs_tree_prior}{a Coalescent Bayesian Skyline tree prior,
as returned by \code{\link{create_cbs_tree_prior}}}
}
\description{
Creates the tree prior section in the prior section of
the prior section of the distribution section
of a BEAST2 XML parameter file for a Birth-Death tree prior
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_site_models.R
\name{check_site_models}
\alias{check_site_models}
\title{Check if the object is a list of one or more site models.}
\usage{
check_site_models(site_models)
}
\arguments{
\item{site_models}{the object to be checked if it is a list of one
or more valid site models}
}
\value{
nothing.
  Will \link{stop} if the object is not a list of one or more site models.
}
\description{
Will \link{stop} if the object is not a list of one or more site models.
}
\examples{
check_site_models(create_jc69_site_model())
check_site_models(list(create_jc69_site_model()))
check_site_models(
  list(create_jc69_site_model(), create_gtr_site_model())
)
}
\seealso{
Use \link{create_site_model} to create a valid site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beauti_options_v2_4.R
\name{create_beauti_options_v2_4}
\alias{create_beauti_options_v2_4}
\title{Function to create the BEAUti options for version 2.4.}
\usage{
create_beauti_options_v2_4()
}
\value{
a BEAUti options structure
}
\description{
Function to create the BEAUti options for version 2.4, by
calling \link{create_beauti_options}.
}
\examples{
beauti_options <- create_beauti_options_v2_4()
xml <- create_beast2_input(
  get_fasta_filename(),
  beauti_options = beauti_options
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_rate_gt}
\alias{parameter_to_xml_rate_gt}
\title{Internal function}
\usage{
parameter_to_xml_rate_gt(
  parameter,
  beauti_options = create_beauti_options(),
  which_name = "state_node"
)
}
\arguments{
\item{parameter}{a 'rate GT' parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_rate_gt_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{which_name}{the name, can be \code{state_node} or \code{rate_name}}
}
\value{
the parameter as XML text
}
\description{
Converts a 'rate GT' parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_scale_param}
\alias{is_scale_param}
\title{Determine if the object is a valid
scale parameter}
\usage{
is_scale_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
scale parameter}
}
\value{
TRUE if x is a valid scale parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
scale parameter
}
\examples{

is_scale_param(create_alpha_param())
is_scale_param(create_beta_param())
is_scale_param(create_clock_rate_param())
is_scale_param(create_kappa_1_param())
is_scale_param(create_kappa_2_param())
is_scale_param(create_lambda_param())
is_scale_param(create_m_param())
is_scale_param(create_mean_param())
is_scale_param(create_mu_param())
is_scale_param(create_rate_ac_param())
is_scale_param(create_rate_ag_param())
is_scale_param(create_rate_at_param())
is_scale_param(create_rate_cg_param())
is_scale_param(create_rate_ct_param())
is_scale_param(create_rate_gt_param())
is_scale_param(create_s_param())
is_scale_param(create_scale_param())
is_scale_param(create_sigma_param())

is_scale_param(NA)
is_scale_param(NULL)
is_scale_param("nonsense")
is_scale_param(create_jc69_site_model())
is_scale_param(create_strict_clock_model())
is_scale_param(create_yule_tree_prior())
is_scale_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_clock_model.R
\name{check_clock_model}
\alias{check_clock_model}
\title{Check if the clock model is a valid clock model.}
\usage{
check_clock_model(clock_model)
}
\arguments{
\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}
}
\value{
TRUE if \code{clock_model} is a valid clock model
}
\description{
Calls \code{stop} if the clock model is invalid
}
\examples{
check_clock_model(create_strict_clock_model())
check_clock_model(create_rln_clock_model())
}
\seealso{
Use \link{create_clock_model} to create a valid clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_laplace_distr}
\alias{is_laplace_distr}
\title{Determine if the object is a valid
Laplace distribution,
as created by \code{\link{create_laplace_distr}}}
\usage{
is_laplace_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
Laplace distribution}
}
\value{
TRUE if x is a valid Laplace distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
Laplace distribution,
as created by \code{\link{create_laplace_distr}}
}
\examples{
# TRUE
is_laplace_distr(create_laplace_distr())
# FALSE
is_laplace_distr(create_log_normal_distr())
is_laplace_distr(NA)
is_laplace_distr(NULL)
is_laplace_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site_model_n_params.R
\name{get_site_model_n_params}
\alias{get_site_model_n_params}
\title{Get the number of distributions a site model has}
\usage{
get_site_model_n_params(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the number of parameters a site model has
}
\description{
Get the number of distributions a site model has
}
\examples{
  testit::assert(
    get_site_model_n_params(create_gtr_site_model()) == 10
  )
  testit::assert(
    get_site_model_n_params(create_hky_site_model()) == 2
  )
  testit::assert(
    get_site_model_n_params(create_jc69_site_model()) == 0
  )
  testit::assert(
    get_site_model_n_params(create_tn93_site_model()) == 4
  )
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_data_sequences.R
\name{create_beast2_input_data_sequences}
\alias{create_beast2_input_data_sequences}
\title{Creates the data section of a BEAST2 XML parameter file}
\usage{
create_beast2_input_data_sequences(
  input_fasta_filename,
  beauti_options = beautier::create_beauti_options()
)
}
\arguments{
\item{input_fasta_filename}{one FASTA filename}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
lines of XML text
}
\description{
Creates the data section of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_distr_names.R
\name{get_distr_names}
\alias{get_distr_names}
\title{Get the distribution names}
\usage{
get_distr_names()
}
\value{
the distribution names
}
\description{
Get the distribution names
}
\examples{
get_distr_names()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_model_to_xml_prior.R
\name{site_model_to_xml_prior_distr}
\alias{site_model_to_xml_prior_distr}
\title{Converts a site model to XML,
  used in the \code{prior} section}
\usage{
site_model_to_xml_prior_distr(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the site model as XML text
}
\description{
Converts a site model to XML,
  used in the \code{prior} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mcmc.R
\name{check_mcmc_list_element_names}
\alias{check_mcmc_list_element_names}
\title{Check if the MCMC has the list elements of a valid MCMC object.}
\usage{
check_mcmc_list_element_names(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element is missing
}
\seealso{
Use \code{\link{create_mcmc}} to create a valid MCMC
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_log_mode.R
\name{check_log_mode}
\alias{check_log_mode}
\title{Check if the supplied \code{mode} is a valid logging mode.}
\usage{
check_log_mode(mode)
}
\arguments{
\item{mode}{mode how to log.
Valid are \code{tree}, \code{autodetect} and \code{compound}}
}
\description{
Check if the supplied \code{mode} is a valid logging mode.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_gamma_distr}
\alias{init_gamma_distr}
\title{Initializes a gamma distribution}
\usage{
init_gamma_distr(gamma_distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{gamma_distr}{a gamma distribution,
using \link{create_gamma_distr}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized gamma distribution
}
\description{
Initializes a gamma distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_sigma}
\alias{parameter_to_xml_sigma}
\title{Internal function}
\usage{
parameter_to_xml_sigma(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a sigma parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_sigma_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a sigma parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_ns_inference_model.R
\name{create_ns_inference_model}
\alias{create_ns_inference_model}
\title{Create an inference model to measure the evidence of.}
\usage{
create_ns_inference_model(
  site_model = beautier::create_jc69_site_model(),
  clock_model = beautier::create_strict_clock_model(),
  tree_prior = beautier::create_yule_tree_prior(),
  mcmc = beautier::create_ns_mcmc()
)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
an inference model
}
\description{
Create an inference model to measure the evidence of.
To do so, the inference model is created as usual (see
\link{create_inference_model}), except
for using a Nested Sampling MCMC (see \link{create_ns_mcmc})
}
\examples{
inference_model <- create_ns_inference_model()
}
\seealso{
Use \link{create_inference_model} to create a
regular  inference model.
Use \link{create_test_ns_inference_model} to create an inference model
to estimate the marginal likelihood with a short MCMC,
to be used in testing.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_beta_distr}
\alias{is_beta_distr}
\title{Determine if the object is a valid
beta distribution,
as created by \code{\link{create_beta_distr}}}
\usage{
is_beta_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
beta distribution,}
}
\value{
TRUE if x is a valid beta distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
beta distribution,
as created by \code{\link{create_beta_distr}}
}
\examples{
# TRUE
is_beta_distr(create_beta_distr())
# FALSE
is_beta_distr(create_exp_distr())
is_beta_distr(NA)
is_beta_distr(NULL)
is_beta_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_model_from_name.R
\name{create_site_model_from_name}
\alias{create_site_model_from_name}
\title{Create a site model from name}
\usage{
create_site_model_from_name(site_model_name)
}
\arguments{
\item{site_model_name}{name of a site model,
must be a name as returned by \code{\link{get_site_model_names}}}
}
\value{
a site model
}
\description{
Create a site model from name
}
\examples{
site_model <- create_site_model_from_name(get_site_model_names()[1])
}
\seealso{
Use \link{create_site_model} to create a site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_rate_cg_param}
\alias{create_rate_cg_param}
\alias{create_param_rate_cg}
\title{Create a parameter called 'rate CG'}
\usage{
create_rate_cg_param(id = NA, estimate = TRUE, value = "1.0", lower = "0.0")
}
\arguments{
\item{id}{the parameter's ID}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}

\item{value}{value of the parameter}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}
}
\value{
a parameter called 'rate CG'
}
\description{
Create a parameter called 'rate CG'
}
\examples{
# Create parameter
rate_cg_param <- create_rate_cg_param(value = 1, estimate = FALSE)

# Use the parameter to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  site_model = create_gtr_site_model(
    rate_cg_param = rate_cg_param
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_clock_model_name.R
\name{get_clock_model_name}
\alias{get_clock_model_name}
\title{Get the BEAUti name for a clock model}
\usage{
get_clock_model_name(clock_model)
}
\arguments{
\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}
}
\value{
name of the clock model
}
\description{
Will \link{stop} if the clock model is an invalid clock model
}
\examples{
# StrictClock
get_clock_model_name(create_strict_clock_model())

# RelaxedClock
get_clock_model_name(create_rln_clock_model())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_rate_ag_param}
\alias{is_rate_ag_param}
\title{Determine if the object is a valid
'rate AG' parameter}
\usage{
is_rate_ag_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
'rate AG' parameter}
}
\value{
TRUE if x is a valid 'rate AG' parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
'rate AG' parameter
}
\examples{

is_rate_ag_param(create_alpha_param())
is_rate_ag_param(create_beta_param())
is_rate_ag_param(create_clock_rate_param())
is_rate_ag_param(create_kappa_1_param())
is_rate_ag_param(create_kappa_2_param())
is_rate_ag_param(create_lambda_param())
is_rate_ag_param(create_m_param())
is_rate_ag_param(create_mean_param())
is_rate_ag_param(create_mu_param())
is_rate_ag_param(create_rate_ac_param())
is_rate_ag_param(create_rate_ag_param())
is_rate_ag_param(create_rate_at_param())
is_rate_ag_param(create_rate_cg_param())
is_rate_ag_param(create_rate_ct_param())
is_rate_ag_param(create_rate_gt_param())
is_rate_ag_param(create_s_param())
is_rate_ag_param(create_scale_param())
is_rate_ag_param(create_sigma_param())

is_rate_ag_param(NA)
is_rate_ag_param(NULL)
is_rate_ag_param("nonsense")
is_rate_ag_param(create_jc69_site_model())
is_rate_ag_param(create_strict_clock_model())
is_rate_ag_param(create_yule_tree_prior())
is_rate_ag_param(create_mcmc())
}
\seealso{
\code{\link{create_rate_ag_param}} creates a 'rate AG' parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_log_modes.R
\name{get_log_modes}
\alias{get_log_modes}
\title{Get the possible log modes}
\usage{
get_log_modes()
}
\description{
Get the possible log modes
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_site_models.R
\name{init_gtr_site_model}
\alias{init_gtr_site_model}
\title{Initializes a GTR site model}
\usage{
init_gtr_site_model(gtr_site_model, distr_id = 0, param_id = 0)
}
\arguments{
\item{gtr_site_model}{a GTR site model,
as returned by \code{\link{create_gtr_site_model}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized GTR site model
}
\description{
Initializes a GTR site model
}
\examples{
gtr_site_model <- create_gtr_site_model()
# FALSE
is_init_gtr_site_model(gtr_site_model)
gtr_site_model <- init_gtr_site_model(gtr_site_model)
# TRUE
is_init_gtr_site_model(gtr_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_sigma_param}
\alias{is_sigma_param}
\title{Determine if the object is a valid
sigma parameter}
\usage{
is_sigma_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
sigma parameter}
}
\value{
TRUE if x is a valid sigma parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
sigma parameter
}
\examples{

is_sigma_param(create_alpha_param())
is_sigma_param(create_beta_param())
is_sigma_param(create_clock_rate_param())
is_sigma_param(create_kappa_1_param())
is_sigma_param(create_kappa_2_param())
is_sigma_param(create_lambda_param())
is_sigma_param(create_m_param())
is_sigma_param(create_mean_param())
is_sigma_param(create_mu_param())
is_sigma_param(create_rate_ac_param())
is_sigma_param(create_rate_ag_param())
is_sigma_param(create_rate_at_param())
is_sigma_param(create_rate_cg_param())
is_sigma_param(create_rate_ct_param())
is_sigma_param(create_rate_gt_param())
is_sigma_param(create_s_param())
is_sigma_param(create_scale_param())
is_sigma_param(create_sigma_param())

is_sigma_param(NA)
is_sigma_param(NULL)
is_sigma_param("nonsense")
is_sigma_param(create_jc69_site_model())
is_sigma_param(create_strict_clock_model())
is_sigma_param(create_yule_tree_prior())
is_sigma_param(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/site_model_to_xml_state.R
\name{site_model_to_xml_state}
\alias{site_model_to_xml_state}
\title{Converts a site model to XML,
  used in the \code{state} section}
\usage{
site_model_to_xml_state(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the site model as XML text
}
\description{
Converts a site model to XML,
  used in the \code{state} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tree_prior_n_distrs.R
\name{get_tree_prior_n_distrs}
\alias{get_tree_prior_n_distrs}
\title{Get the number of distributions a tree prior has}
\usage{
get_tree_prior_n_distrs(tree_prior)
}
\arguments{
\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
the number of distributions a tree prior has
}
\description{
Get the number of distributions a tree prior has
}
\examples{
# 2: birth_rate_distr and death_rate_distr
get_tree_prior_n_distrs(create_bd_tree_prior())

# 0:none
get_tree_prior_n_distrs(create_cbs_tree_prior())

# 1: pop_size_distr
get_tree_prior_n_distrs(create_ccp_tree_prior())

 # 2:pop_size_distr and growth_rate_distr
get_tree_prior_n_distrs(create_cep_tree_prior())

# 1: birth_rate_distr
get_tree_prior_n_distrs(create_yule_tree_prior())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_normal_distr}
\alias{is_normal_distr}
\title{Determine if the object is a valid
normal distribution
as created by \code{\link{create_normal_distr}}}
\usage{
is_normal_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
normal distribution}
}
\value{
TRUE if x is a valid normal distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
normal distribution
as created by \code{\link{create_normal_distr}}
}
\examples{
# TRUE
is_normal_distr(create_normal_distr())
# FALSE
is_normal_distr(create_one_div_x_distr())
is_normal_distr(NA)
is_normal_distr(NULL)
is_normal_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_id.R
\name{is_id}
\alias{is_id}
\title{Determine if the object is a valid ID}
\usage{
is_id(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid ID}
}
\value{
TRUE if x is a valid ID, FALSE otherwise
}
\description{
Determine if the object is a valid ID
}
\examples{
# TRUE
is_id("anthus_aco")
is_id(3)
# FALSE
is_id(ape::rcoal(3))
is_id(NULL)
is_id(NA)
}
\seealso{
to check multiple IDs, use \link{are_ids}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_one_div_x}
\alias{distr_to_xml_one_div_x}
\title{Internal function}
\usage{
distr_to_xml_one_div_x(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{a 1/x distribution,
as created by \code{\link{create_one_div_x_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a 1/x distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mrca_prior_name.R
\name{check_mrca_prior_name}
\alias{check_mrca_prior_name}
\title{Check if \code{mrca_prior_name} is a valid MRCA prior name.}
\usage{
check_mrca_prior_name(mrca_prior_name)
}
\arguments{
\item{mrca_prior_name}{the unique name of the MRCA prior,
for example a genus, family,
order or even class name.
Leave at \link{NA} to have it named automatically.}
}
\description{
A valid MRCA prior name is either \link{NA} or one character string.
Will \link{stop} if not.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/compare_lines.R
\name{compare_lines}
\alias{compare_lines}
\title{Internal function}
\usage{
compare_lines(
  lines,
  expected,
  section = NA,
  created_lines_filename = "created.xml",
  expected_lines_filename = "expected.xml"
)
}
\arguments{
\item{lines}{the created lines}

\item{expected}{the expected/goal/target lines}

\item{section}{the XML section. Leave at NA to compare all lines}

\item{created_lines_filename}{name of the file where the (section of
the) created lines are stored}

\item{expected_lines_filename}{name of the file where the (section of
the) expected lines are stored}
}
\value{
nothing. Instead, two files are created, with the
  names \code{created_lines_filename}
  and \code{expected_lines_filename} that contain the
  section under investigation, so that a diff tool
  can compare these
}
\description{
Internal debug function to compare the actually created
lines to expected lines using any diff tool
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_site_models.R
\name{init_site_models}
\alias{init_site_models}
\title{Initializes all site models}
\usage{
init_site_models(site_models, ids, distr_id = 0, param_id = 0)
}
\arguments{
\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}

\item{ids}{one or more alignments' IDs.
IDs can be extracted from their FASTA filenames
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{distr_id}{the first distributions' ID}

\item{param_id}{the first parameter's ID}
}
\value{
a list of initialized site models
}
\description{
Initializes all site models
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_freq_equilibrium_names.R
\name{get_freq_equilibrium_names}
\alias{get_freq_equilibrium_names}
\title{Returns valid values for the \code{freq_equilibrium} argument}
\usage{
get_freq_equilibrium_names()
}
\value{
the valid values for the \code{freq_equilibrium} argument
}
\description{
Returns valid values for the \code{freq_equilibrium} argument
}
\examples{
get_freq_equilibrium_names()
}
\seealso{
the \code{freq_equilibrium} argument is used in
  \code{\link{create_gtr_site_model}},
  \code{\link{create_hky_site_model}},
  and \code{\link{create_tn93_site_model}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_one_div_x_distr}
\alias{create_one_div_x_distr}
\alias{create_distr_one_div_x}
\title{Create a 1/x distribution}
\usage{
create_one_div_x_distr(id = NA, value = NA, lower = NA, upper = NA)
}
\arguments{
\item{id}{the distribution's ID}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
a 1/x distribution
}
\description{
Create a 1/x distribution
}
\examples{
one_div_x_distr <- create_one_div_x_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = one_div_x_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrca_priors_to_xml_distr.R
\name{mrca_priors_to_xml_prior_distr}
\alias{mrca_priors_to_xml_prior_distr}
\title{Creates the the \code{distribution}'s prior section (which is part of
a posterior distribution section) of a BEAST2 XML parameter file.}
\usage{
mrca_priors_to_xml_prior_distr(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
These lines start with '\code{<distribution id="prior"}'
}
\details{
\code{
   <distribution id="posterior" spec="util.CompoundDistribution">
       <distribution id="prior" spec="util.CompoundDistribution">
         HERE, where the ID of the distribution is 'prior'
       </distribution>
       <distribution id="likelihood" ...>
       </distribution>
  </distribution>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_phylogeny.R
\name{check_phylogeny}
\alias{check_phylogeny}
\title{Check if the phylogeny is a valid phylogeny object.}
\usage{
check_phylogeny(phylogeny)
}
\arguments{
\item{phylogeny}{a phylogeny of type \code{phylo} from the \code{ape}
package}
}
\value{
nothing
}
\description{
Calls \code{stop} if the phylogeny is invalid
}
\examples{

# Must do nothing on phylogenies
phylogeny <- ape::read.tree(text = "(A:1, B:1):1;")
check_phylogeny(phylogeny)
}
\seealso{
Use \code{ape::read.tree} to create a phylogeny
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_clock_model.R
\name{create_clock_model}
\alias{create_clock_model}
\title{General function to create a clock model}
\usage{
create_clock_model(name, id, ...)
}
\arguments{
\item{name}{the clock model name. Valid
names can be found in \code{get_clock_model_names}}

\item{id}{a clock model's ID}

\item{...}{specific clock model parameters}
}
\value{
a valid clock model
}
\description{
General function to create a clock model
}
\note{
Prefer using the named function
  \code{\link{create_rln_clock_model}}
  and \code{\link{create_strict_clock_model}}
}
\examples{
rln_clock_model <- create_rln_clock_model()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  clock_model = rln_clock_model
)
file.remove(beast2_input_file)

strict_clock_model <- create_strict_clock_model()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  clock_model = strict_clock_model
)
file.remove(beast2_input_file)
}
\seealso{
An alignment ID can be extracted from
  its FASTA filename using \code{\link{get_alignment_id}}.
  For more examples about creating a relaxed log-normal clock
  model, see \code{\link{create_rln_clock_model}}.
  For more examples about creating a strict clock
  model, see \code{\link{create_strict_clock_model}}.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equivalent_xmls.R
\name{are_equivalent_xml_lines_section}
\alias{are_equivalent_xml_lines_section}
\title{Determine if XML lines result in equivalent trees}
\usage{
are_equivalent_xml_lines_section(lines_1, lines_2, section, verbose = FALSE)
}
\arguments{
\item{lines_1}{lines of a first XML file}

\item{lines_2}{lines of a second XML file}

\item{section}{the name of the XML section}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
TRUE if the two XML lines result in equivalent trees,
  FALSE otherwise
}
\description{
Determine if XML lines result in equivalent trees
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_mrca_prior.R
\name{init_mrca_prior}
\alias{init_mrca_prior}
\title{Initialize the MRCA prior.}
\usage{
init_mrca_prior(input_filename, inference_model)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
Initialized by
\itemize{
  \item if no alignment ID is set,
    it is set by reading it from the alignment file
  \item if no taxa names are set,
    these are set by reading these from the alignment file
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_branch_rate_model_stuff_xml.R
\name{create_branch_rate_model_stuff_xml}
\alias{create_branch_rate_model_stuff_xml}
\title{Internal function called by \link{create_branch_rate_model_xml}}
\usage{
create_branch_rate_model_stuff_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
It generates the desired XML for some circumstances.
Yes, that is a vague description.
Would be nice if someone would untangle this :-)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_in_patterns.R
\name{is_in_patterns}
\alias{is_in_patterns}
\title{Is there at least one regular expression having a match with the line?}
\usage{
is_in_patterns(line, patterns)
}
\arguments{
\item{line}{a line of text}

\item{patterns}{one or more regular expression patterns}
}
\value{
TRUE if there is at least one match found
}
\description{
Is there at least one regular expression having a match with the line?
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_normal}
\alias{distr_to_xml_normal}
\title{Internal function}
\usage{
distr_to_xml_normal(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{a normal distribution,
as created by \code{\link{create_normal_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a normal distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_rate_gt_param}
\alias{create_rate_gt_param}
\alias{create_param_rate_gt}
\title{Create a parameter called 'rate GT'}
\usage{
create_rate_gt_param(id = NA, estimate = TRUE, value = "1.0", lower = "0.0")
}
\arguments{
\item{id}{the parameter's ID}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}

\item{value}{value of the parameter}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}
}
\value{
a parameter called 'rate GT'
}
\description{
Create a parameter called 'rate GT'
}
\examples{
# Create parameter
rate_gt_param <- create_rate_gt_param(value = 1, estimate = FALSE)

# Use the parameter to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  site_model = create_gtr_site_model(
    rate_gt_param = rate_gt_param
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_site_model_name.R
\name{is_site_model_name}
\alias{is_site_model_name}
\title{Determines if the name is a valid site_model name}
\usage{
is_site_model_name(name)
}
\arguments{
\item{name}{the name to be tested}
}
\value{
TRUE if the name is a valid site_model name, FALSE otherwise
}
\description{
Determines if the name is a valid site_model name
}
\examples{
# TRUE
is_site_model_name("JC69")
is_site_model_name("HKY")
is_site_model_name("TN93")
is_site_model_name("GTR")
# FALSE
is_site_model_name("nonsense")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_tree_prior.R
\name{is_bd_tree_prior}
\alias{is_bd_tree_prior}
\title{Determine if the object is a valid Birth Death tree prior}
\usage{
is_bd_tree_prior(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid birth death tree prior}
}
\value{
TRUE if x is a valid birth death tree prior, FALSE otherwise
}
\description{
Determine if the object is a valid Birth Death tree prior
}
\examples{
  testit::assert(is_bd_tree_prior(create_bd_tree_prior()))
  testit::assert(!is_bd_tree_prior(create_cbs_tree_prior()))
  testit::assert(!is_bd_tree_prior(create_ccp_tree_prior()))
  testit::assert(!is_bd_tree_prior(create_cep_tree_prior()))
  testit::assert(!is_bd_tree_prior(create_yule_tree_prior()))
}
\seealso{
Use \code{\link{create_bd_tree_prior}} to create a valid
  Birth-Death tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_screenlog.R
\name{check_screenlog}
\alias{check_screenlog}
\title{Check if a \code{screenlog} is valid.}
\usage{
check_screenlog(screenlog)
}
\arguments{
\item{screenlog}{a \code{screenlog},
as created by \link{create_screenlog}}
}
\description{
Will call \link{stop} if not.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_tree_priors.R
\name{init_tree_priors}
\alias{init_tree_priors}
\title{Initializes all tree priors}
\usage{
init_tree_priors(tree_priors, ids, distr_id = 0, param_id = 0)
}
\arguments{
\item{tree_priors}{one or more tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{ids}{one or more alignments' IDs.
IDs can be extracted from their FASTA filenames
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{distr_id}{the first distributions' ID}

\item{param_id}{the first parameter's ID}
}
\value{
a list of initialized tree priors
}
\description{
Initializes all tree priors
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beautier_tempfilename.R
\name{get_beautier_tempfilename}
\alias{get_beautier_tempfilename}
\title{Get a temporary filename}
\usage{
get_beautier_tempfilename(pattern = "file", fileext = "")
}
\arguments{
\item{pattern}{a non-empty character vector
giving the initial part of the name.}

\item{fileext}{a non-empty character vector
giving the file extension}
}
\value{
name for a temporary file
}
\description{
Get a temporary filename, similar to \link{tempfile},
except that it always writes to a temporary folder
named \link{beautier}.
}
\note{
this function is added to make sure no temporary
cache files are left undeleted
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/default_params_doc.R
\name{default_parameters_doc}
\alias{default_parameters_doc}
\title{Documentation of parameters (for example, \code{create_param}.
This function does nothing. It is intended to inherit documentation from.}
\usage{
default_parameters_doc(estimate, id, lower, name, upper, value, ...)
}
\arguments{
\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}

\item{id}{the parameter's ID}

\item{lower}{lowest possible value of the parameter. If the parameter
is estimated, \code{lower} must be less than \code{value}}

\item{name}{the parameters' name. Valid
names can be found in \code{get_param_names}}

\item{upper}{upper value of the parameter}

\item{value}{value of the parameter}

\item{...}{specific parameter parameters}
}
\description{
Documentation of parameters (for example, \code{create_param}.
This function does nothing. It is intended to inherit documentation from.
}
\note{
This is an internal function, so it should be marked with
  \code{@export}. This is not done, as this will disallow all
  functions to find the documentation parameters
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_gamma_site_model.R
\name{init_gamma_site_model}
\alias{init_gamma_site_model}
\title{Initializes a gamma site model}
\usage{
init_gamma_site_model(gamma_site_model, distr_id = 0, param_id = 0)
}
\arguments{
\item{gamma_site_model}{a site model's gamma site model,
as returned by \code{\link{create_gamma_site_model}}}

\item{distr_id}{the first distributions' ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized gamma site model
}
\description{
Initializes a gamma site model
}
\examples{
gamma_site_model <- create_gamma_site_model(
  gamma_cat_count = 2,
  gamma_shape_prior_distr = create_one_div_x_distr(id = NA)
)
# FALSE: not yet initialized
is_init_gamma_site_model(gamma_site_model)
gamma_site_model <- init_gamma_site_model(gamma_site_model)
# TRUE: now it is initialized
is_init_gamma_site_model(gamma_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beautier_paths.R
\name{get_beautier_paths}
\alias{get_beautier_paths}
\title{Get the full paths of files in the \code{inst/extdata} folder}
\usage{
get_beautier_paths(filenames)
}
\arguments{
\item{filenames}{the files' names, without the path}
}
\value{
the filenames' full paths
}
\description{
Get the full paths of files in the \code{inst/extdata} folder
}
\examples{
  testit::assert(
    length(
      get_beautier_paths(
        c("test_output_0.fas", "anthus_aco.fas", "anthus_nd2.fas")
      )
     ) == 3
   )
}
\seealso{
Use \link{get_beautier_path} to get the path of one file

for one file, use \code{\link{get_beautier_path}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_distr_n_params.R
\name{get_distr_n_params}
\alias{get_distr_n_params}
\title{Get the number of parameters a distribution uses}
\usage{
get_distr_n_params(distr)
}
\arguments{
\item{distr}{a distribution,
as created by \code{\link{create_distr}} or (preferable)
its named functions}
}
\value{
the number of parameters that distribution uses
}
\description{
Get the number of parameters a distribution uses
}
\examples{
get_distr_n_params(create_beta_distr())
get_distr_n_params(create_exp_distr())
get_distr_n_params(create_gamma_distr())
get_distr_n_params(create_inv_gamma_distr())
get_distr_n_params(create_laplace_distr())
get_distr_n_params(create_log_normal_distr())
get_distr_n_params(create_normal_distr())
get_distr_n_params(create_one_div_x_distr())
get_distr_n_params(create_poisson_distr())
get_distr_n_params(create_uniform_distr())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_rate_ct}
\alias{parameter_to_xml_rate_ct}
\title{Internal function}
\usage{
parameter_to_xml_rate_ct(
  parameter,
  beauti_options = create_beauti_options(),
  which_name = "state_node"
)
}
\arguments{
\item{parameter}{a 'rate CT' parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_rate_ct_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{which_name}{the name, can be \code{state_node} or \code{rate_name}}
}
\value{
the parameter as XML text
}
\description{
Converts a 'rate CT' parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_one_div_x_distr}
\alias{init_one_div_x_distr}
\title{Initializes an one-divided-by-x distribution}
\usage{
init_one_div_x_distr(one_div_x_distr, distr_id = 0)
}
\arguments{
\item{one_div_x_distr}{a one-divided-by-x distribution,
using \link{create_one_div_x_distr}}

\item{distr_id}{the first distribution's ID}
}
\value{
an initialized one-divided-by-x distribution
}
\description{
Initializes an one-divided-by-x distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_beta_distr}
\alias{is_init_beta_distr}
\title{Determine if x is an initialized beta distribution object
  as created by \code{\link{create_beta_distr}}}
\usage{
is_init_beta_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized beta distribution object}
}
\value{
TRUE if x is an initialized beta distribution object
}
\description{
Determine if x is an initialized beta distribution object
  as created by \code{\link{create_beta_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_model_to_xml_prior_distr.R
\name{clock_model_to_xml_prior_distr}
\alias{clock_model_to_xml_prior_distr}
\title{Internal function}
\usage{
clock_model_to_xml_prior_distr(
  inference_model,
  clock_model = "deprecated",
  mrca_priors = "deprecated",
  tipdates_filename = "deprecated"
)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
a character vector of XML strings
}
\description{
Internal function to converts a clock model
to the \code{prior} section of the XML as text
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rename_mcmc_filenames.R
\name{rename_mcmc_filenames}
\alias{rename_mcmc_filenames}
\title{Rename the filenames within an MCMC}
\usage{
rename_mcmc_filenames(mcmc, rename_fun)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}

\item{rename_fun}{a function to rename a filename,
as can be checked by \link{check_rename_fun}. This function should
have one argument, which will be a filename or \link{NA}. The
function should \link{return} one filename (when passed one filename) or
one \link{NA} (when passed one \link{NA}).
Example rename functions are:
\itemize{
  \item \link{get_remove_dir_fun} get a function that removes the directory
    paths from the filenames, in effect turning these into local files
  \item \link{get_replace_dir_fun} get a function that replaces the directory
    paths from the filenames
  \item \link{get_remove_hex_fun} get a function that removes the
    hex string from filenames.
    For example, \code{tracelog_82c1a522040.log} becomes \code{tracelog.log}
}}
}
\description{
Rename the filenames within an MCMC
}
\examples{

# Create an MCMC with local filenames
mcmc <- create_mcmc()
mcmc$tracelog$filename <- "trace.log"
mcmc$screenlog$filename <- "screen.log"
mcmc$treelog$filename <- "tree.log"

# Nah, files should be put in '/home/john' folder
mcmc <- rename_mcmc_filenames(
  mcmc = mcmc,
  rename_fun = get_replace_dir_fun("/home/john")
)

# Nah, files should be put in '/home/doe' folder instead
mcmc <- rename_mcmc_filenames(
  mcmc = mcmc,
  rename_fun = get_replace_dir_fun("/home/doe")
)

# Nah, files should be put in local folder instead
mcmc <- rename_mcmc_filenames(
  mcmc = mcmc,
  rename_fun = get_remove_dir_fun()
)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interspace.R
\name{interspace}
\alias{interspace}
\title{Puts spaces in between the lines}
\usage{
interspace(lines)
}
\arguments{
\item{lines}{lines of text}
}
\value{
interspaced lines of text
}
\description{
Puts spaces in between the lines
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tracelog.R
\name{create_tracelog}
\alias{create_tracelog}
\title{Create a \code{tracelog} object}
\usage{
create_tracelog(
  filename = NA,
  log_every = 1000,
  mode = "autodetect",
  sanitise_headers = TRUE,
  sort = "smart"
)
}
\arguments{
\item{filename}{name of the file to store the posterior traces.
Use \link{NA} to use the filename \code{[alignment_id].log},
where \code{alignment_id} is obtained using \link{get_alignment_id}}

\item{log_every}{number of MCMC states between writing to file}

\item{mode}{mode how to log.
Valid values are the ones returned by \link{get_log_modes}}

\item{sanitise_headers}{set to \link{TRUE} to sanitise the headers of the
log file}

\item{sort}{how to sort the log.
Valid values are the ones returned by \link{get_log_sorts}}
}
\description{
Create a \code{tracelog} object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_rate_gt_param}
\alias{is_rate_gt_param}
\title{Determine if the object is a valid
'rate GT' parameter}
\usage{
is_rate_gt_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
'rate GT' parameter}
}
\value{
TRUE if x is a valid 'rate GT' parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
'rate GT' parameter
}
\examples{

is_rate_gt_param(create_alpha_param())
is_rate_gt_param(create_beta_param())
is_rate_gt_param(create_clock_rate_param())
is_rate_gt_param(create_kappa_1_param())
is_rate_gt_param(create_kappa_2_param())
is_rate_gt_param(create_lambda_param())
is_rate_gt_param(create_m_param())
is_rate_gt_param(create_mean_param())
is_rate_gt_param(create_mu_param())
is_rate_gt_param(create_rate_ac_param())
is_rate_gt_param(create_rate_ag_param())
is_rate_gt_param(create_rate_at_param())
is_rate_gt_param(create_rate_cg_param())
is_rate_gt_param(create_rate_ct_param())
is_rate_gt_param(create_rate_gt_param())
is_rate_gt_param(create_s_param())
is_rate_gt_param(create_scale_param())
is_rate_gt_param(create_sigma_param())

is_rate_gt_param(NA)
is_rate_gt_param(NULL)
is_rate_gt_param("nonsense")
is_rate_gt_param(create_jc69_site_model())
is_rate_gt_param(create_strict_clock_model())
is_rate_gt_param(create_yule_tree_prior())
is_rate_gt_param(create_mcmc())
}
\seealso{
\code{\link{create_rate_gt_param}} creates a 'rate GT' parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrca_prior_to_xml_prior_distr.R
\name{mrca_prior_to_xml_prior_distr}
\alias{mrca_prior_to_xml_prior_distr}
\title{Creates the distribution section in the prior section of the
distribution section of a BEAST2 XML parameter file.}
\usage{
mrca_prior_to_xml_prior_distr(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
These lines start with '<distribution id='
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxa_to_xml_tree.R
\name{tipdate_taxa_to_xml_tree}
\alias{tipdate_taxa_to_xml_tree}
\title{Internal function}
\usage{
tipdate_taxa_to_xml_tree(
  inference_model,
  id = "deprecated",
  tipdates_filename = "deprecated"
)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}

\item{id}{an alignment's IDs.
An ID can be extracted from its FASTA filename
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
the random phylogeny as XML text
}
\description{
Creates the \code{tree} section
(part of the \code{state} section)
when there is tip-dating
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_log_sort.R
\name{check_log_sort}
\alias{check_log_sort}
\title{Check if the supplied \code{sort} is a valid logging sorting option.}
\usage{
check_log_sort(sort)
}
\arguments{
\item{sort}{how to sort the entries in a log.
Valid are \code{smart}, \code{none} and \code{alphabetic}}
}
\description{
Check if the supplied \code{sort} is a valid logging sorting option.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_lambda_param}
\alias{is_lambda_param}
\title{Determine if the object is a valid lambda parameter}
\usage{
is_lambda_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
lambda parameter}
}
\value{
TRUE if x is a valid lambda parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid lambda parameter
}
\examples{

is_lambda_param(create_alpha_param())
is_lambda_param(create_beta_param())
is_lambda_param(create_clock_rate_param())
is_lambda_param(create_kappa_1_param())
is_lambda_param(create_kappa_2_param())
is_lambda_param(create_lambda_param())
is_lambda_param(create_m_param())
is_lambda_param(create_mean_param())
is_lambda_param(create_mu_param())
is_lambda_param(create_rate_ac_param())
is_lambda_param(create_rate_ag_param())
is_lambda_param(create_rate_at_param())
is_lambda_param(create_rate_cg_param())
is_lambda_param(create_rate_ct_param())
is_lambda_param(create_rate_gt_param())
is_lambda_param(create_s_param())
is_lambda_param(create_scale_param())
is_lambda_param(create_sigma_param())

is_lambda_param(NA)
is_lambda_param(NULL)
is_lambda_param("nonsense")
is_lambda_param(create_jc69_site_model())
is_lambda_param(create_strict_clock_model())
is_lambda_param(create_yule_tree_prior())
is_lambda_param(create_mcmc())
}
\seealso{
lambda parameters are returned by \code{\link{create_lambda_param}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equivalent_xmls.R
\name{are_equivalent_xml_lines_loggers}
\alias{are_equivalent_xml_lines_loggers}
\title{Determine if XML operator lines result in equivalent trees}
\usage{
are_equivalent_xml_lines_loggers(lines_1, lines_2, verbose = FALSE)
}
\arguments{
\item{lines_1}{lines of a first XML file}

\item{lines_2}{lines of a second XML file}

\item{verbose}{if TRUE, additional information is displayed, that
is potentially useful in debugging}
}
\value{
TRUE if the two XML lines result in equivalent trees,
  FALSE otherwise
}
\description{
Determine if XML operator lines result in equivalent trees
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_fasta_filenames.R
\name{are_fasta_filenames}
\alias{are_fasta_filenames}
\title{Checks if all filenames have a FASTA filename extension}
\usage{
are_fasta_filenames(filenames)
}
\arguments{
\item{filenames}{filenames}
}
\value{
TRUE if all filenames have a FASTA filename extension
}
\description{
Checks if all filenames have a FASTA filename extension
}
\examples{
# TRUE
are_fasta_filenames("1.fas")
are_fasta_filenames("1.fasta")
are_fasta_filenames("1.FAS")
are_fasta_filenames("1.FASTA")
are_fasta_filenames(c("1.fas", "2.fas"))

# FALSE
are_fasta_filenames("")
are_fasta_filenames(NA)
are_fasta_filenames(NULL)
are_fasta_filenames(Inf)
are_fasta_filenames("1.fasX")
are_fasta_filenames(c("1.fas", "2.exe"))
are_fasta_filenames(c("1.bat", "2.exe"))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_rate_ct_param}
\alias{is_rate_ct_param}
\title{Determine if the object is a valid
'rate CT' parameter}
\usage{
is_rate_ct_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
'rate CT' parameter}
}
\value{
TRUE if x is a valid 'rate CG' parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid
'rate CT' parameter
}
\examples{

is_rate_ct_param(create_alpha_param())
is_rate_ct_param(create_beta_param())
is_rate_ct_param(create_clock_rate_param())
is_rate_ct_param(create_kappa_1_param())
is_rate_ct_param(create_kappa_2_param())
is_rate_ct_param(create_lambda_param())
is_rate_ct_param(create_m_param())
is_rate_ct_param(create_mean_param())
is_rate_ct_param(create_mu_param())
is_rate_ct_param(create_rate_ac_param())
is_rate_ct_param(create_rate_ag_param())
is_rate_ct_param(create_rate_at_param())
is_rate_ct_param(create_rate_cg_param())
is_rate_ct_param(create_rate_ct_param())
is_rate_ct_param(create_rate_gt_param())
is_rate_ct_param(create_s_param())
is_rate_ct_param(create_scale_param())
is_rate_ct_param(create_sigma_param())

is_rate_ct_param(NA)
is_rate_ct_param(NULL)
is_rate_ct_param("nonsense")
is_rate_ct_param(create_jc69_site_model())
is_rate_ct_param(create_strict_clock_model())
is_rate_ct_param(create_yule_tree_prior())
is_rate_ct_param(create_mcmc())
}
\seealso{
\code{\link{create_rate_ct_param}} creates a 'rate CT' parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tree_priors_n_distrs.R
\name{get_tree_priors_n_distrs}
\alias{get_tree_priors_n_distrs}
\title{Get the number of distributions a tree prior has}
\usage{
get_tree_priors_n_distrs(tree_priors)
}
\arguments{
\item{tree_priors}{one or more tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
the number of distributions a tree prior has
}
\description{
Get the number of distributions a tree prior has
}
\examples{
# Three
get_tree_priors_n_distrs(
  list(
    create_bd_tree_prior(), # has two distributions
    create_ccp_tree_prior() # has one distribution
  )
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_model_to_xml_operators.R
\name{clock_model_to_xml_operators}
\alias{clock_model_to_xml_operators}
\title{Converts a clock model to the \code{operators} section of the
XML as text}
\usage{
clock_model_to_xml_operators(
  inference_model,
  clock_model = "deprecated",
  mrca_priors = "deprecated",
  tipdates_filename = "deprecated"
)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
a character vector of XML strings
}
\description{
Converts a clock model to the \code{operators} section of the
XML as text
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/has_xml_closing_tag.R
\name{has_xml_closing_tag}
\alias{has_xml_closing_tag}
\title{Is an XML closing tag with the value of \code{section}
  present among the lines of
  the text?}
\usage{
has_xml_closing_tag(lines, section)
}
\arguments{
\item{lines}{lines of the XML text}

\item{section}{the XML section}
}
\value{
TRUE if there is an XML closing tag with the value of
  \code{section} present. FALSE otherwise
}
\description{
Is an XML closing tag with the value of \code{section}
  present among the lines of
  the text?
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_xml_section_from_lines.R
\name{extract_xml_section_from_lines}
\alias{extract_xml_section_from_lines}
\title{Get the lines of an XML section, including the section tags}
\usage{
extract_xml_section_from_lines(lines, section)
}
\arguments{
\item{lines}{lines of the XML text}

\item{section}{the XML section name}
}
\value{
the section's lines of XML text, including the tags
}
\description{
Get the lines of an XML section, including the section tags
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equal_tracelogs.R
\name{are_equal_tracelogs}
\alias{are_equal_tracelogs}
\title{Determine if two tracelogs are equal.}
\usage{
are_equal_tracelogs(tracelog_1, tracelog_2)
}
\arguments{
\item{tracelog_1}{an tracelog, as created by \link{create_tracelog}}

\item{tracelog_2}{an tracelog, as created by \link{create_tracelog}}
}
\value{
TRUE if the two tracelogs are equal
}
\description{
Will \link{stop} if the arguments are not tracelogs.
}
\examples{

tracelog_1 <- create_tracelog(log_every = 1000)
tracelog_2 <- create_tracelog(log_every = 314)
# TRUE
are_equal_tracelogs(tracelog_1, tracelog_1)
# FALSE
are_equal_tracelogs(tracelog_1, tracelog_2)
}
\seealso{
Use \link{create_tracelog} to create an tracelog
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_clock_rate_param}
\alias{create_clock_rate_param}
\alias{create_param_clock_rate}
\title{Create a parameter called \code{clock_rate},
  as needed by \code{\link{create_strict_clock_model}}}
\usage{
create_clock_rate_param(value = "1.0", estimate = FALSE, id = NA)
}
\arguments{
\item{value}{value of the parameter}

\item{estimate}{TRUE if this parameter is to be estimated by BEAST2,
FALSE otherwise}

\item{id}{the parameter's ID}
}
\value{
a parameter called rate
}
\description{
Create a parameter called \code{clock_rate},
  as needed by \code{\link{create_strict_clock_model}}
}
\note{
It cannot be estimated (as a hyper parameter) yet.
}
\examples{
clock_rate_param <- create_clock_rate_param(
  id = "anthus_aco", value = 1.0
)

# Use the parameter in a clock model
strict_clock_model <- create_strict_clock_model(
  clock_rate_param = clock_rate_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  clock_model = strict_clock_model
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_clock_models.R
\name{check_clock_models}
\alias{check_clock_models}
\title{Check if the object is a list of one or more clock models.}
\usage{
check_clock_models(clock_models)
}
\arguments{
\item{clock_models}{the object to be checked if it is a list of one
or more valid clock models}
}
\value{
nothing.
  Will \link{stop} if the object is not a list of one or more clock models.
}
\description{
Will \link{stop} if the object is not a list of one or more clock models.
}
\examples{
check_clock_models(create_strict_clock_model())
check_clock_models(list(create_strict_clock_model()))
check_clock_models(
  list(create_strict_clock_model(), create_rln_clock_model())
)
}
\seealso{
Use \link{create_clock_model} to create a valid clock model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_tree_prior.R
\name{is_tree_prior}
\alias{is_tree_prior}
\title{Determine if an object is a valid tree prior}
\usage{
is_tree_prior(x)
}
\arguments{
\item{x}{an object}
}
\value{
TRUE if x is a valid tree_prior, FALSE otherwise
}
\description{
Determine if an object is a valid tree prior
}
\examples{
  testit::assert(is_tree_prior(create_bd_tree_prior()))
  testit::assert(is_tree_prior(create_yule_tree_prior()))
  testit::assert(!is_tree_prior("nonsense"))
}
\seealso{
tree priors can be created by \code{\link{create_tree_prior}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/indent.R
\name{indent}
\alias{indent}
\title{Indent text for a certain number of spaces.
If the text is only whitespace, leave it as such}
\usage{
indent(text, n_spaces = 4)
}
\arguments{
\item{text}{the text to indent}

\item{n_spaces}{the number of spaces to add before the text. BEAUti uses
four spaces by default}
}
\value{
the indented text
}
\description{
Indent text for a certain number of spaces.
If the text is only whitespace, leave it as such
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mcmc_nested_sampling.R
\name{create_ns_mcmc}
\alias{create_ns_mcmc}
\alias{create_mcmc_nested_sampling}
\title{Create an MCMC object to estimate the marginal likelihood
using Nested Sampling.}
\usage{
create_ns_mcmc(
  chain_length = 1e+07,
  store_every = -1,
  pre_burnin = 0,
  n_init_attempts = 3,
  particle_count = 1,
  sub_chain_length = 5000,
  epsilon = "1e-12",
  tracelog = beautier::create_tracelog(),
  screenlog = beautier::create_screenlog(),
  treelog = beautier::create_treelog()
)
}
\arguments{
\item{chain_length}{upper bound to the length of the MCMC chain}

\item{store_every}{number of states the MCMC will process
before the posterior's state will be saved to file.
Use -1 or \code{NA} to use the default frequency.}

\item{pre_burnin}{number of burn in samples taken before entering
the main loop}

\item{n_init_attempts}{number of initialization attempts before failing}

\item{particle_count}{number of particles}

\item{sub_chain_length}{sub-chain length}

\item{epsilon}{epsilon}

\item{tracelog}{a \code{tracelog},
as created by \link{create_tracelog}}

\item{screenlog}{a \code{screenlog},
as created by \link{create_screenlog}}

\item{treelog}{a \code{treelog},
as created by \link{create_treelog}}
}
\value{
an MCMC object
}
\description{
This will result in a BEAST run that estimates the marginal
likelihood until convergence is achieved.
In this context, \code{chain_length} is only an upper bound
to the length of that run.
}
\examples{
mcmc <- create_ns_mcmc(
  chain_length = 1e7,
  store_every = 1000,
  particle_count = 1,
  sub_chain_length = 1000,
  epsilon = 1e-12
)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  mcmc = mcmc
)
file.remove(beast2_input_file)
}
\references{
* [1] Patricio Maturana Russel, Brendon J Brewer, Steffen Klaere,
    Remco R Bouckaert; Model Selection and Parameter Inference in
    Phylogenetics Using Nested Sampling, Systematic Biology, 2018,
    syy050, https://doi.org/10.1093/sysbio/syy050
}
\seealso{
Use \code{\link{create_mcmc}} to create a regular MCMC.
Use \code{\link{create_test_ns_mcmc}} to create an NS MCMC for testing,
  with, among others, a short MCMC chain length.
Use \code{\link{check_ns_mcmc}} to check that an NS MCMC object is valid.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_xml.R
\name{create_xml_declaration}
\alias{create_xml_declaration}
\title{Create the XML declaration of the BEAST2 XML input file}
\usage{
create_xml_declaration()
}
\value{
one line of XML text
}
\description{
Create the XML declaration of the BEAST2 XML input file
}
\examples{
create_xml_declaration()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr_name.R
\name{is_distr_name}
\alias{is_distr_name}
\title{Determines if the name is a valid distribution name}
\usage{
is_distr_name(name)
}
\arguments{
\item{name}{the name to be tested}
}
\value{
TRUE if the name is a valid distribution name, FALSE otherwise
}
\description{
Determines if the name is a valid distribution name
}
\examples{
# TRUE
is_distr_name("uniform")
is_distr_name("normal")
is_distr_name("one_div_x")
is_distr_name("log_normal")
is_distr_name("exponential")
is_distr_name("gamma")
is_distr_name("beta")
is_distr_name("laplace")
is_distr_name("inv_gamma")
is_distr_name("poisson")
# FALSE
is_distr_name("nonsense")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_freq_equilibrium_name.R
\name{is_freq_equilibrium_name}
\alias{is_freq_equilibrium_name}
\title{Checks if \code{name} is a valid \code{freq_equilibrium} argument value}
\usage{
is_freq_equilibrium_name(name)
}
\arguments{
\item{name}{the name to check if it is a valid \code{freq_equilibrium}
argument value}
}
\value{
TRUE if the name is a valid \code{freq_equilibrium} value
}
\description{
Checks if \code{name} is a valid \code{freq_equilibrium} argument value
}
\examples{
# TRUE
is_freq_equilibrium_name("estimated")
is_freq_equilibrium_name("empirical")
is_freq_equilibrium_name("all_equal")
# FALSE
is_freq_equilibrium_name("nonsense")
}
\seealso{
the \code{freq_equilibrium} argument is used by
  \code{\link{create_gtr_site_model}},
  \code{\link{create_hky_site_model}},
  and \code{\link{create_tn93_site_model}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/beautier.R
\docType{package}
\name{beautier}
\alias{beautier}
\title{\code{beautier}: A package to create a \code{BEAST2} input file.}
\description{
\code{beautier} allows to create a \code{BEAST2} input file, using
an R interface. \code{beautier} closely follows the interface
of \code{BEAUti 2}, a GUI tool bundled with \code{BEAST2}, including
its default settings.
}
\details{
See the documentation of \code{create_inference_model} to see the
features of BEAST2 that \code{beautier} supports.
}
\examples{
# Get an example FASTA file
input_filename <- get_fasta_filename()

# The file created by beautier, a BEAST2 input file
output_filename <- get_beautier_tempfilename()

# Use the default BEAUti settings to create a BEAST2 input file
create_beast2_input_file_from_model(
  input_filename,
  output_filename,
  inference_model = create_inference_model()
)
file.remove(output_filename)
}
\seealso{
These are packages associated with \code{beautier}:
\itemize{
  \item{
    The package \code{beastier} can run
    BEAST2 from R
  }
  \item{
    The package \code{tracerer} can parse
    BEAST2 output files from R
  }
  \item{
    The package \code{mauricer} manages
    BEAST2 packages from R
  }
  \item{
    The package \code{babette} combines the
    functionality of \code{beautier},
    \code{beastier}, \code{mauricer} and \code{tracerer}
    into a single workflow
  }
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_site_model.R
\name{is_hky_site_model}
\alias{is_hky_site_model}
\title{Determine if the object is a valid HKY site model,
as created by \code{\link{create_hky_site_model}}}
\usage{
is_hky_site_model(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid HKY site model}
}
\value{
TRUE if x is a valid HKY site model, FALSE otherwise
}
\description{
Determine if the object is a valid HKY site model,
as created by \code{\link{create_hky_site_model}}
}
\examples{

# site models
is_hky_site_model(create_hky_site_model())
is_hky_site_model(create_gtr_site_model())
is_hky_site_model(create_jc69_site_model())
is_hky_site_model(create_tn93_site_model())

# other models
is_hky_site_model(NA)
is_hky_site_model(NULL)
is_hky_site_model("nonsense")
is_hky_site_model(create_strict_clock_model())
is_hky_site_model(create_bd_tree_prior())
is_hky_site_model(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_inference_model.R
\name{is_inference_model}
\alias{is_inference_model}
\title{Determine if the input is an inference model}
\usage{
is_inference_model(x)
}
\arguments{
\item{x}{object to be determined of if it is an inference model}
}
\value{
TRUE if the object is an inference model
}
\description{
Determine if the input is an inference model
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_n_taxa.R
\name{get_n_taxa}
\alias{get_n_taxa}
\title{Extract the number of taxa from a file}
\usage{
get_n_taxa(filename)
}
\arguments{
\item{filename}{name of a FASTA file}
}
\value{
the number of taxa
}
\description{
Extract the number of taxa from a file
}
\examples{
fasta_filename <- get_beautier_path("test_output_5.fas")
# 5
get_n_taxa(fasta_filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_site_model.R
\name{is_init_hky_site_model}
\alias{is_init_hky_site_model}
\title{Determine if x is an initialized HKY site model
as created by \code{\link{create_hky_site_model}}}
\usage{
is_init_hky_site_model(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized HKY site model}
}
\value{
TRUE if x is an initialized HKY site model
}
\description{
Determine if x is an initialized HKY site model
as created by \code{\link{create_hky_site_model}}
}
\examples{

hky_site_model <- create_hky_site_model()
# FALSE: not yet initialized
is_init_hky_site_model(hky_site_model)
hky_site_model <- init_hky_site_model(hky_site_model)
# TRUE: now it is initialized
is_init_hky_site_model(hky_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_beast.R
\name{create_beast2_input_beast}
\alias{create_beast2_input_beast}
\title{Creates the XML text for the \code{beast} tag of a BEAST2 parameter file.}
\usage{
create_beast2_input_beast(
  input_filename,
  inference_model = beautier::create_inference_model()
)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Creates the XML text for the \code{beast} tag of a BEAST2 parameter file,
which is directly after the XML
declaration (created by \link{create_xml_declaration}.
}
\details{
The \code{beast} tag has these elements:
\preformatted{
  <beast[...]>
      <data
      [...]
      </data>
      [map names]
      <run[...]>
      [...]
      </run>
  </beast>
}
}
\seealso{
Use \link{create_beast2_input_from_model} to create the complete XML text.
Use \link{create_beast2_input_data} to create the XML text for
  the \code{data} tag only.
Use \link{create_beast2_input_map} to create the XML text for
  the \code{[map names]} part.
Use \link{create_beast2_input_run} to create the XML text for
  the \code{run} tag only.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_treelog.R
\name{check_treelog}
\alias{check_treelog}
\title{Check if a \code{treelog} is valid.}
\usage{
check_treelog(treelog)
}
\arguments{
\item{treelog}{a \code{treelog},
as created by \link{create_treelog}}
}
\description{
Will call \link{stop} if not.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_temp_screenlog_filename.R
\name{create_temp_screenlog_filename}
\alias{create_temp_screenlog_filename}
\title{Create a filename for a temporary screenlog file}
\usage{
create_temp_screenlog_filename()
}
\description{
Create a filename for a temporary screenlog file
}
\seealso{
use \link{create_screenlog} to create a screenlog.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_test_inference_model.R
\name{create_test_inference_model}
\alias{create_test_inference_model}
\title{Create a testing inference model.}
\usage{
create_test_inference_model(
  site_model = beautier::create_jc69_site_model(),
  clock_model = beautier::create_strict_clock_model(),
  tree_prior = beautier::create_yule_tree_prior(),
  mrca_prior = NA,
  mcmc = beautier::create_test_mcmc(),
  beauti_options = beautier::create_beauti_options(),
  tipdates_filename = NA
)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
an inference model
}
\description{
Creates a simple inference model with a short MCMC chain,
to be used in testing.
}
\examples{
inference_model <- create_test_inference_model()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file_from_model(
  get_fasta_filename(),
  beast2_input_file,
  inference_model = inference_model
)
file.remove(beast2_input_file)
}
\seealso{
Use \link{create_inference_model} to create a
regular inference model.
Use \link{create_test_ns_inference_model} to create an inference model
to estimate the marginal likelihood with a short MCMC, to be
used in testing
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fasta_to_sequences.R
\name{fasta_file_to_sequences}
\alias{fasta_file_to_sequences}
\title{Convert a FASTA file to a table of sequences}
\usage{
fasta_file_to_sequences(fasta_filename)
}
\arguments{
\item{fasta_filename}{One existing FASTA filenames}
}
\value{
a table of sequences
}
\description{
Convert a FASTA file to a table of sequences
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_one_na.R
\name{is_one_na}
\alias{is_one_na}
\title{Determines if x is one NA}
\usage{
is_one_na(x)
}
\arguments{
\item{x}{the object to be determined if it is one NA}
}
\value{
TRUE if x is one NA, FALSE otherwise
}
\description{
Determines if x is one NA
}
\examples{
  testit::assert(is_one_na(NA))
  testit::assert(!is_one_na(NULL))
  testit::assert(!is_one_na(42))
  testit::assert(!is_one_na("Hello"))
  testit::assert(!is_one_na(3.14))
  testit::assert(!is_one_na(c(NA, NA)))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input.R
\name{create_beast2_input}
\alias{create_beast2_input}
\title{Create a BEAST2 XML input text}
\usage{
create_beast2_input(
  input_filename,
  tipdates_filename = NA,
  site_model = beautier::create_jc69_site_model(),
  clock_model = beautier::create_strict_clock_model(),
  tree_prior = beautier::create_yule_tree_prior(),
  mrca_prior = NA,
  mcmc = beautier::create_mcmc(),
  beauti_options = beautier::create_beauti_options()
)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}

\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
a character vector of XML strings
}
\description{
Create a BEAST2 XML input text
}
\examples{
  text <- create_beast2_input(
    input_filename = get_fasta_filename()
  )
  testit::assert(substr(text[1], 1, 5) == "<?xml")
  text[1]
  testit::assert(tail(text, n = 1) == "</beast>")
}
\seealso{
Use \link{create_beast2_input_from_model} to create the BEAST2 XML
  input text from an inference model
  Use \link{create_beast2_input_file} to also save it to file.

\code{\link{create_beast2_input_file}} shows more examples
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_site_model.R
\name{is_tn93_site_model}
\alias{is_tn93_site_model}
\title{Determine if the object is a valid TN93 site model,}
\usage{
is_tn93_site_model(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid TN93 site model,
as created by \code{\link{create_tn93_site_model}}}
}
\value{
TRUE if x is a valid TN93 site model, FALSE otherwise
}
\description{
Determine if the object is a valid TN93 site model,
}
\examples{

# site models
is_tn93_site_model(create_gtr_site_model())
is_tn93_site_model(create_hky_site_model())
is_tn93_site_model(create_jc69_site_model())
is_tn93_site_model(create_tn93_site_model())

# other models
is_tn93_site_model(NA)
is_tn93_site_model(NULL)
is_tn93_site_model("nonsense")
is_tn93_site_model("")
is_tn93_site_model(c())
is_tn93_site_model(create_strict_clock_model())
is_tn93_site_model(create_bd_tree_prior())
is_tn93_site_model(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_model_xml.R
\name{create_site_model_xml}
\alias{create_site_model_xml}
\title{Internal function to creates the XML text for the \code{siteModel} tag
of a BEAST2 parameter file.}
\usage{
create_site_model_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the site model as XML text
}
\description{
Creates the XML text for the \code{siteModel} tag of
a BEAST2 parameter file,
which is part of the \code{distribution} node for the
\code{treeLikelihood} ID.
}
\details{
The \code{siteModel} tag has these elements:

\preformatted{
  <siteModel[...]>

      [parameters]

      <substModel[...]>
        [...]
      </substModel>
  </siteModel>
}

The \code{parameter} section is created by
\link{create_site_model_parameters_xml}
The \code{substModel} section is created by \link{create_subst_model_xml}
}
\examples{
 # <distribution id="posterior"[...]">
 #     <distribution id="likelihood" [...]>
 #       <siteModel...>
 #         [parameters]
 #       </siteModel>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_treelog.R
\name{create_treelog}
\alias{create_treelog}
\title{Create a \code{treelog} object}
\usage{
create_treelog(
  filename = "$(tree).trees",
  log_every = 1000,
  mode = "tree",
  sanitise_headers = FALSE,
  sort = "none"
)
}
\arguments{
\item{filename}{name of the file to store the posterior trees}

\item{log_every}{number of MCMC states between writing to file}

\item{mode}{mode how to log.
Valid values are the ones returned by \link{get_log_modes}}

\item{sanitise_headers}{set to \link{TRUE} to sanitise the headers of the
log file}

\item{sort}{how to sort the log.
Valid values are the ones returned by \link{get_log_sorts}}
}
\description{
Create a \code{treelog} object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_remove_dir_fun.R
\name{get_remove_dir_fun}
\alias{get_remove_dir_fun}
\title{Get a function that, from a filename, returns the part
without the directory.}
\usage{
get_remove_dir_fun()
}
\description{
Or: get a function that returns the local version of a filename.
Also, the function will return \link{NA} if the filename is \link{NA}
}
\seealso{
see \link{check_rename_fun}
  for an overview of file renaming functions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_mrca_prior.R
\name{create_mrca_prior}
\alias{create_mrca_prior}
\title{Create a Most Recent Common Ancestor prior}
\usage{
create_mrca_prior(
  alignment_id = NA,
  taxa_names = NA,
  is_monophyletic = FALSE,
  mrca_distr = NA,
  name = NA,
  clock_prior_distr_id = NA
)
}
\arguments{
\item{alignment_id}{ID of the alignment,
as returned by \link{get_alignment_id}.
Keep at \code{NA} to have it initialized automatically}

\item{taxa_names}{names of the taxa,
as returned by \code{\link{get_taxa_names}}.
Keep at \code{NA} to have it initialized automatically,
using all taxa in the alignment}

\item{is_monophyletic}{boolean to indicate monophyly is assumed in
a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{mrca_distr}{the distribution used by the MRCA prior.
Can be NA (the default) or any distribution
returned by \code{\link{create_distr}}}

\item{name}{the unique name of the MRCA prior, for example a genus, family,
order or even class name.
Leave at \link{NA} to have it named automatically.}

\item{clock_prior_distr_id}{ID of an MRCA clock model's distribution.
Keep at \code{NA} to have it initialized automatically}
}
\value{
an MRCA prior
}
\description{
Create a Most Recent Common Ancestor prior
}
\examples{
 fasta_filename <- get_beautier_path("anthus_aco.fas")

 # The first two taxa are sister species
 mrca_prior <- create_mrca_prior(
   taxa_names = get_taxa_names(filename = fasta_filename)[1:2]
 )

 # The taxa are monophyletic
 mrca_prior <- create_mrca_prior(
   taxa_names = get_taxa_names(filename = fasta_filename),
   is_monophyletic = TRUE
 )

 # Set the crown age to 10
 mrca_prior <- create_mrca_prior(
   taxa_names = get_taxa_names(fasta_filename),
   mrca_distr = create_normal_distr(mean = 10, sigma = 0.1)
 )
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_beta_distr}
\alias{init_beta_distr}
\title{Initializes a beta distribution}
\usage{
init_beta_distr(beta_distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{beta_distr}{a beta distribution,
using \link{create_beta_distr}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized beta distribution
}
\description{
Initializes a beta distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_file_and_model_agree.R
\name{check_file_and_model_agree}
\alias{check_file_and_model_agree}
\title{Checks if the input FASTA file and the inference model agree.}
\usage{
check_file_and_model_agree(input_filename, inference_model)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
Will \link{stop} if not
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_site_model.R
\name{is_gtr_site_model}
\alias{is_gtr_site_model}
\title{Determine if the object is a valid GTR site model,
as created by \code{\link{create_gtr_site_model}}}
\usage{
is_gtr_site_model(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid GTR site model}
}
\value{
TRUE if x is a valid GTR site model, FALSE otherwise
}
\description{
Determine if the object is a valid GTR site model,
as created by \code{\link{create_gtr_site_model}}
}
\examples{

# site models
is_gtr_site_model(create_gtr_site_model())
is_gtr_site_model(create_hky_site_model())
is_gtr_site_model(create_jc69_site_model())
is_gtr_site_model(create_tn93_site_model())

# other models
is_gtr_site_model(NA)
is_gtr_site_model(NULL)
is_gtr_site_model("nonsense")
is_gtr_site_model(create_strict_clock_model())
is_gtr_site_model(create_bd_tree_prior())
is_gtr_site_model(create_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_clock_models.R
\name{init_clock_models}
\alias{init_clock_models}
\title{Initializes all clock models}
\usage{
init_clock_models(fasta_filenames, clock_models, distr_id = 0, param_id = 0)
}
\arguments{
\item{fasta_filenames}{One or more FASTA filenames.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}

\item{distr_id}{the first distributions' ID}

\item{param_id}{the first parameter's ID}
}
\value{
a list of initialized clock models
}
\description{
Initializes all clock models
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_multiline.R
\name{remove_multiline}
\alias{remove_multiline}
\title{Remove consecutive lines}
\usage{
remove_multiline(text, lines_to_remove)
}
\arguments{
\item{text}{lines of characters}

\item{lines_to_remove}{lines of character that need to be removed from text}
}
\value{
lines of text
}
\description{
Remove consecutive lines
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_gamma_site_model.R
\name{is_init_gamma_site_model}
\alias{is_init_gamma_site_model}
\title{Determine if x is an initialized gamma site model,
as created by \code{\link{create_gamma_site_model}}}
\usage{
is_init_gamma_site_model(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized gamma site_models object}
}
\value{
TRUE if x is an initialized gamma site model
}
\description{
Determine if x is an initialized gamma site model,
as created by \code{\link{create_gamma_site_model}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_loggers_xml.R
\name{create_loggers_xml}
\alias{create_loggers_xml}
\title{Creates the three logger sections of a BEAST2 XML parameter file}
\usage{
create_loggers_xml(input_filename, inference_model)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
The logger section has these elements:
\preformatted{
 <logger id="tracelog" [...]>
     [...]
 </logger>
 <logger id="screenlog" [...]>
     [...]
 </logger>
 <logger id="treelog.t:[alignment ID]"  [...]>
     [...]
 </logger>
}
}
\seealso{
Use \link{create_tracelog_xml} to create the XML text
of the logger with the \code{tracelog} ID.
Use \link{create_screenlog_xml} to create the XML text
of the logger with the \code{screenlog} ID.
Use \link{create_treelog_xml} to create the XML text
of the loggers with the \code{treelog} ID.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_xml_loggers_from_lines.R
\name{extract_xml_loggers_from_lines}
\alias{extract_xml_loggers_from_lines}
\title{Extract everything between first loggers and last loggers line}
\usage{
extract_xml_loggers_from_lines(lines)
}
\arguments{
\item{lines}{lines of text}
}
\value{
lines of text from the first to and including the last operators line
}
\description{
Extract everything between first loggers and last loggers line
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_is_monophyletic.R
\name{check_is_monophyletic}
\alias{check_is_monophyletic}
\title{Check if \code{is_monophyletic} has a valid value.}
\usage{
check_is_monophyletic(is_monophyletic)
}
\arguments{
\item{is_monophyletic}{boolean to indicate monophyly is assumed in
a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}
}
\description{
Will \link{stop} if not.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrca_prior_to_xml_tracelog.R
\name{mrca_prior_to_xml_tracelog}
\alias{mrca_prior_to_xml_tracelog}
\title{Internal function}
\usage{
mrca_prior_to_xml_tracelog(
  inference_model,
  clock_models = "deprecated",
  mrca_prior = "deprecated",
  tipdates_filename = "deprecated"
)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}

\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}

\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
lines of XML text
}
\description{
Internal function to creates the MRCA prior's XML for the tracelog section.
}
\details{
\code{
  <logger id="tracelog" ...>
    # Here
  </logger>
}
}
\seealso{
all MRCA priors' tracelog section is created
  by \code{\link{mrca_priors_to_xml_tracelog}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_scale}
\alias{parameter_to_xml_scale}
\title{Internal function}
\usage{
parameter_to_xml_scale(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a scale parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_scale_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a scale parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_gtr_site_model.R
\name{check_gtr_site_model_names}
\alias{check_gtr_site_model_names}
\title{Check if the \code{gtr_site_model} has the list elements
of a valid \code{gtr_site_model} object.}
\usage{
check_gtr_site_model_names(gtr_site_model)
}
\arguments{
\item{gtr_site_model}{a GTR site model,
as returned by \code{\link{create_gtr_site_model}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element is missing
}
\seealso{
Use \link{create_gtr_site_model}
to create a valid \code{gtr_site_model}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_priors_to_xml_tracelog.R
\name{tree_priors_to_xml_tracelog}
\alias{tree_priors_to_xml_tracelog}
\title{Creates the tree priors' XML for the tracelog section}
\usage{
tree_priors_to_xml_tracelog(tree_priors)
}
\arguments{
\item{tree_priors}{one or more tree priors,
as returned by \code{\link{create_tree_prior}}}
}
\value{
lines of XML text
}
\description{
Creates the tree priors' XML for the tracelog section
}
\examples{
# <logger id="tracelog" ...>
#'   # Here
# </logger>
}
\seealso{
the complete tracelog section is created
  by \code{\link{create_tracelog_xml}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_screenlog_xml.R
\name{create_screenlog_xml}
\alias{create_screenlog_xml}
\title{Creates the \code{screenlog} section of the \code{logger} section
of a BEAST2 XML parameter file}
\usage{
create_screenlog_xml(inference_model = beautier::create_inference_model())
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the XML text
}
\description{
Creates the \code{screenlog} section of the \code{logger} section
of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_exp_distr}
\alias{create_exp_distr}
\alias{create_distr_exp}
\title{Create an exponential distribution}
\usage{
create_exp_distr(id = NA, mean = 1, value = NA, lower = NA, upper = NA)
}
\arguments{
\item{id}{the distribution's ID}

\item{mean}{the mean parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_mean_param}}}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
an exponential distribution
}
\description{
Create an exponential distribution
}
\examples{
exp_distr <- create_exp_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = exp_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tree_likelihood_distr_xml.R
\name{create_tree_likelihood_distr_xml}
\alias{create_tree_likelihood_distr_xml}
\title{Creates the XML text for the \code{distribution} tag
with the \code{treeLikelihood} ID,
of a BEAST2 parameter file.}
\usage{
create_tree_likelihood_distr_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
Creates the XML text for the \code{distribution} tag
with the \code{treeLikelihood} ID,
of a BEAST2 parameter file,
in an unindented form
}
\details{
The \code{distribution} tag (with ID equals \code{treeLikelihood})
has these elements:

\preformatted{
   <distribution id="treeLikelihood"[...]>
      <siteModel[...]>
        [...]
      </siteModel>
      <branchRateModel[...]>
        [...]
      </branchRateModel>
   </distribution>
}

The \code{siteModel} section
is created by \link{create_site_model_xml}.
The \code{branchRateModel} section
is created by \link{create_branch_rate_model_xml}.

Zooming out:

\preformatted{
  <beast[...]>
    <run[...]>
      <distribution id="posterior"[...]>
        <distribution id="likelihood"[...]>
          [this section]
        </distribution>
      </distribution>
    </run>
  </beast>
}
}
\note{
this function is not intended for regular use, thus its
  long name length is accepted
}
\seealso{
this function is called by \code{create_beast2_input_distr},
  together with \code{create_beast2_input_distr_prior}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_site_model.R
\name{is_init_tn93_site_model}
\alias{is_init_tn93_site_model}
\title{Determine if x is an initialized tn93 site model
as created by \code{\link{create_tn93_site_model}}}
\usage{
is_init_tn93_site_model(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized TN93 site model}
}
\value{
TRUE if x is an initialized TN93 site model
}
\description{
Determine if x is an initialized tn93 site model
as created by \code{\link{create_tn93_site_model}}
}
\examples{

tn93_site_model <- create_tn93_site_model()
# FALSE: not yet initialized
is_init_tn93_site_model(tn93_site_model)
tn93_site_model <- init_tn93_site_model(tn93_site_model)
# TRUE: now it is initialized
is_init_tn93_site_model(tn93_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_normal_distr}
\alias{is_init_normal_distr}
\title{Determine if x is an initialized normal distribution object
  as created by \code{\link{create_normal_distr}}}
\usage{
is_init_normal_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized normal distribution object}
}
\value{
TRUE if x is an initialized normal distribution object
}
\description{
Determine if x is an initialized normal distribution object
  as created by \code{\link{create_normal_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_one_bool.R
\name{is_one_bool}
\alias{is_one_bool}
\title{Check if the argument is one boolean}
\usage{
is_one_bool(x)
}
\arguments{
\item{x}{the argument to be tested to be boolean}
}
\description{
Check if the argument is one boolean
}
\examples{
# TRUE
is_one_bool(TRUE)
is_one_bool(FALSE)

# FALSE
is_one_bool(NULL)
is_one_bool(NA)
is_one_bool(c())
is_one_bool("nonsense")
is_one_bool(is_one_bool)
is_one_bool(c(TRUE, FALSE))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_param.R
\name{init_param}
\alias{init_param}
\title{Initializes a parameter}
\usage{
init_param(param, id)
}
\arguments{
\item{param}{a parameter,
using \code{\link{create_param}}}

\item{id}{the parameter's ID. Will be ignored if the parameter already
has an ID}
}
\value{
an initialized parameter
}
\description{
Initializes a parameter
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site_models_n_distrs.R
\name{get_site_models_n_distrs}
\alias{get_site_models_n_distrs}
\title{Get the number of distributions a site model has}
\usage{
get_site_models_n_distrs(site_models)
}
\arguments{
\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}
}
\value{
the number of distributions the site models have
}
\description{
Get the number of distributions a site model has
}
\examples{
# 5
get_site_models_n_distrs(list(create_gtr_site_model()))
# 1
get_site_models_n_distrs(list(create_hky_site_model()))
# 0
get_site_models_n_distrs(list(create_jc69_site_model()))
# 2
get_site_models_n_distrs(list(create_tn93_site_model()))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_beautier_path.R
\name{get_beautier_path}
\alias{get_beautier_path}
\title{Get the full path of a file in the \code{inst/extdata} folder}
\usage{
get_beautier_path(filename)
}
\arguments{
\item{filename}{the file's name, without the path}
}
\value{
the full path of the filename
}
\description{
Get the full path of a file in the \code{inst/extdata} folder
}
\examples{
get_beautier_path("test_output_0.fas")
get_beautier_path("anthus_aco.fas")
get_beautier_path("anthus_nd2.fas")
}
\seealso{
for more files, use \code{\link{get_beautier_paths}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_tree_priors.R
\name{check_tree_priors}
\alias{check_tree_priors}
\title{Check if the object is a list of one or more tree priors.}
\usage{
check_tree_priors(tree_priors)
}
\arguments{
\item{tree_priors}{the object to be checked if it is a list of one
or more valid tree priors}
}
\value{
nothing.
  Will \link{stop} if the object is not a list of one or more tree priors.
}
\description{
Will \link{stop} if the object is not a list of one or more tree priors.
}
\examples{
check_tree_priors(create_yule_tree_prior())
check_tree_priors(list(create_yule_tree_prior()))
check_tree_priors(list(create_yule_tree_prior(), create_bd_tree_prior()))
}
\seealso{
Use \link{create_tree_prior} to create a valid tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mrca_prior.R
\name{check_mrca_prior_names}
\alias{check_mrca_prior_names}
\title{Check if the MRCA prior,
which is a list, has all the named elements.}
\usage{
check_mrca_prior_names(mrca_prior)
}
\arguments{
\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if not.
}
\seealso{
Use \link{check_mrca_prior} to check the entire MRCA prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_tree_prior.R
\name{is_init_cbs_tree_prior}
\alias{is_init_cbs_tree_prior}
\title{Determine if x is an initialized Coalescent Bayesian Skyline
  tree_prior object}
\usage{
is_init_cbs_tree_prior(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized Coalescent Bayesian Skyline tree prior object}
}
\value{
TRUE if x is an initialized Coalescent Bayesian Skyline
  tree prior object
}
\description{
Determine if x is an initialized Coalescent Bayesian Skyline
  tree_prior object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_clock_model_name.R
\name{is_clock_model_name}
\alias{is_clock_model_name}
\title{Determines if the name is a valid clock model name}
\usage{
is_clock_model_name(name)
}
\arguments{
\item{name}{the name to be tested}
}
\value{
TRUE if the name is a valid clock_model name, FALSE otherwise
}
\description{
Determines if the name is a valid clock model name
}
\examples{
# TRUE
is_clock_model_name("relaxed_log_normal")
is_clock_model_name("strict")
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_tree_priors.R
\name{init_ccp_tree_prior}
\alias{init_ccp_tree_prior}
\title{Initializes a Coalescent Constant Population tree prior}
\usage{
init_ccp_tree_prior(ccp_tree_prior, distr_id, param_id)
}
\arguments{
\item{ccp_tree_prior}{a Coalescent Constant Population tree prior,
as returned by \code{\link{create_ccp_tree_prior}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized Coalescent Constant Population tree prior
}
\description{
Initializes a Coalescent Constant Population tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_test_screenlog.R
\name{create_test_screenlog}
\alias{create_test_screenlog}
\title{Create a \code{screenlog} object}
\usage{
create_test_screenlog(
  filename = beautier::create_temp_screenlog_filename(),
  log_every = 1000,
  mode = "autodetect",
  sanitise_headers = FALSE,
  sort = "none"
)
}
\arguments{
\item{filename}{name of the file to store the posterior screens
phylogenies to. By default, this is \code{$(screen).screens}}

\item{log_every}{number of MCMC states between writing to file}

\item{mode}{mode how to log.
Valid values are the ones returned by \link{get_log_modes}}

\item{sanitise_headers}{set to \link{TRUE} to sanitise the headers of the
log file}

\item{sort}{how to sort the log.
Valid values are the ones returned by \link{get_log_sorts}}
}
\description{
Create a \code{screenlog} object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_rate_ag}
\alias{parameter_to_xml_rate_ag}
\title{Internal function}
\usage{
parameter_to_xml_rate_ag(
  parameter,
  beauti_options = create_beauti_options(),
  which_name = "state_node"
)
}
\arguments{
\item{parameter}{a 'rate AG' parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_rate_ag_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{which_name}{the name, can be \code{state_node} or \code{rate_name}}
}
\value{
the parameter as XML text
}
\description{
Converts a 'rate AG' parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_uniform}
\alias{distr_to_xml_uniform}
\title{Internal function}
\usage{
distr_to_xml_uniform(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{a uniform distribution,
as created by \code{\link{create_uniform_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a uniform distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_test_ns_inference_model.R
\name{create_test_ns_inference_model}
\alias{create_test_ns_inference_model}
\title{Create an inference model to be tested by Nested Sampling}
\usage{
create_test_ns_inference_model(
  site_model = beautier::create_jc69_site_model(),
  clock_model = beautier::create_strict_clock_model(),
  tree_prior = beautier::create_yule_tree_prior(),
  mcmc = beautier::create_test_ns_mcmc()
)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
an inference model
}
\description{
Create an inference model to be tested by Nested Sampling
}
\examples{
inference_model <- create_test_ns_inference_model()
}
\seealso{
Use \link{create_test_inference_model} to create a
regular inference model with a short MCMC, to be used in testing.
Use \link{create_ns_inference_model} to create an inference model
to estimate the marginal likelihood.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_test_filenames.R
\name{get_fasta_filename}
\alias{get_fasta_filename}
\title{Get the path of a FASTA file used in testing}
\usage{
get_fasta_filename()
}
\value{
the path of a FASTA file used in testing
}
\description{
Get the path of a FASTA file used in testing
}
\examples{
input_filename <- beautier::get_fasta_filename()
output_filename <- get_beautier_tempfilename()

create_beast2_input_file(
  input_filename = input_filename,
  output_filename = output_filename
)

file.remove(output_filename)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param_name.R
\name{is_param_name}
\alias{is_param_name}
\title{Determines if the name is a valid parameter name}
\usage{
is_param_name(name)
}
\arguments{
\item{name}{the name to be tested}
}
\value{
TRUE if the name is a valid parameter name, FALSE otherwise
}
\description{
Determines if the name is a valid parameter name
}
\examples{
# TRUE
is_param_name("alpha")
is_param_name("beta")
is_param_name("clock_rate")
is_param_name("kappa_1")
is_param_name("kappa_2")
is_param_name("lambda")
is_param_name("m")
is_param_name("mean")
is_param_name("mu")
is_param_name("rate_ac")
is_param_name("rate_ag")
is_param_name("rate_at")
is_param_name("rate_cg")
is_param_name("rate_ct")
is_param_name("rate_gt")
is_param_name("s")
is_param_name("scale")
is_param_name("sigma")

# FALSE
is_param_name("nonsense")
is_param_name(NA)
is_param_name(NULL)
is_param_name("")
is_param_name(c())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_models_to_xml_state.R
\name{clock_models_to_xml_state}
\alias{clock_models_to_xml_state}
\title{Deprecated internal function}
\usage{
clock_models_to_xml_state(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text, without indentation nor \code{state}
  tags
}
\description{
Converts one or more clock models to the \code{state} section of the
XML as text
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_tree_priors.R
\name{init_cep_tree_prior}
\alias{init_cep_tree_prior}
\title{Initializes a Coalescent Exponential Population tree prior}
\usage{
init_cep_tree_prior(cep_tree_prior, distr_id, param_id)
}
\arguments{
\item{cep_tree_prior}{a Coalescent Exponential Population tree prior,
as returned by \code{\link{create_cep_tree_prior}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized Coalescent Exponential Population tree prior
}
\description{
Initializes a Coalescent Exponential Population tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_empty_beautier_folder.R
\name{check_empty_beautier_folder}
\alias{check_empty_beautier_folder}
\title{Check there are no files in the default \link{beautier} folder}
\usage{
check_empty_beautier_folder(beautier_folder = get_beautier_folder())
}
\arguments{
\item{beautier_folder}{the path to
the \link{beautier} temporary files folder}
}
\value{
Nothing.
}
\description{
Check there are no files in the default \link{beautier} folder.
The goal is to make sure no temporary files are left undeleted.
Will \link{stop} if there are files in the \link{beautier} folder
}
\examples{
check_empty_beautier_folder()
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_prior_to_xml_operators.R
\name{tree_prior_to_xml_operators}
\alias{tree_prior_to_xml_operators}
\title{Internal function}
\usage{
tree_prior_to_xml_operators(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the tree prior as XML text
}
\description{
Creates the XML of a tree prior,
  as used in the \code{operators} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_site_model.R
\name{check_site_model_types}
\alias{check_site_model_types}
\title{Check if the \code{site_model} has the list elements
of the right type for a valid \code{site_model} object.}
\usage{
check_site_model_types(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element has the incorrect type
}
\seealso{
Use \link{create_site_model} to create a valid \code{site_model}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_test_ns_mcmc.R
\name{create_test_ns_mcmc}
\alias{create_test_ns_mcmc}
\title{Create an NS MCMC object for testing}
\usage{
create_test_ns_mcmc(
  chain_length = 2000,
  store_every = 1000,
  pre_burnin = 0,
  n_init_attempts = 3,
  particle_count = 1,
  sub_chain_length = 500,
  epsilon = 1e-12,
  tracelog = create_test_tracelog(),
  screenlog = create_test_screenlog(),
  treelog = create_test_treelog()
)
}
\arguments{
\item{chain_length}{upper bound to the length of the MCMC chain}

\item{store_every}{number of states the MCMC will process
before the posterior's state will be saved to file.
Use -1 or \code{NA} to use the default frequency.}

\item{pre_burnin}{number of burn in samples taken before entering
the main loop}

\item{n_init_attempts}{number of initialization attempts before failing}

\item{particle_count}{number of particles}

\item{sub_chain_length}{sub-chain length}

\item{epsilon}{epsilon}

\item{tracelog}{a \code{tracelog},
as created by \link{create_tracelog}}

\item{screenlog}{a \code{screenlog},
as created by \link{create_screenlog}}

\item{treelog}{a \code{treelog},
as created by \link{create_treelog}}
}
\value{
an MCMC object
}
\description{
Create an NS MCMC object for testing
}
\examples{
mcmc <- create_test_ns_mcmc()
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  mcmc = mcmc
)
file.remove(beast2_input_file)
}
\seealso{
Use \code{\link{create_ns_mcmc}} to create a default
nested sampling MCMC
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_subst_model_xml.R
\name{create_hky_subst_model_xml}
\alias{create_hky_subst_model_xml}
\title{Converts a site model to XML,
  used in the \code{substModel} section}
\usage{
create_hky_subst_model_xml(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the site model as XML text
}
\description{
Converts a site model to XML,
  used in the \code{substModel} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_treelog_xml.R
\name{create_treelog_xml}
\alias{create_treelog_xml}
\title{Creates the XML text for the \code{logger} tag with ID \code{treelog}.
This section has these elements:
\preformatted{
<logger id="treelog.t:test_output_0" spec="Logger" fileName="my_treelog.trees" logEvery="345000" mode="tree" sanitiseHeaders="true" sort="smart"> # nolint indeed long
    <log id="TreeWithMetaDataLogger.t:test_output_0" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:test_output_0"/> # nolint indeed long
</logger>
}}
\usage{
create_treelog_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\description{
Creates the XML text for the \code{logger} tag with ID \code{treelog}.
This section has these elements:
\preformatted{
<logger id="treelog.t:test_output_0" spec="Logger" fileName="my_treelog.trees" logEvery="345000" mode="tree" sanitiseHeaders="true" sort="smart"> # nolint indeed long
    <log id="TreeWithMetaDataLogger.t:test_output_0" spec="beast.evolution.tree.TreeWithMetaDataLogger" tree="@Tree.t:test_output_0"/> # nolint indeed long
</logger>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ccp_tree_prior_to_xml_state.R
\name{ccp_tree_prior_to_xml_state}
\alias{ccp_tree_prior_to_xml_state}
\title{Convert a CCP tree prior
to the XML as part of the \code{state} section}
\usage{
ccp_tree_prior_to_xml_state(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
XML as text
}
\description{
Convert a CCP tree prior
to the XML as part of the \code{state} section
}
\examples{
# Need an ID and inital value
inference_model <- create_inference_model(
  tree_prior = create_ccp_tree_prior(
    id = "anthus_nd2_sub",
    pop_size_distr = create_normal_distr(
      id = 123,
      value = 3.14
    )
  )
)

ccp_tree_prior_to_xml_state(inference_model)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_laplace_distr}
\alias{is_init_laplace_distr}
\title{Determine if x is an initialized Laplace distribution
  as created by \code{\link{create_laplace_distr}}}
\usage{
is_init_laplace_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized Laplace distribution}
}
\value{
TRUE if x is an initialized Laplace distribution
}
\description{
Determine if x is an initialized Laplace distribution
  as created by \code{\link{create_laplace_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_map.R
\name{create_beast2_input_map}
\alias{create_beast2_input_map}
\title{Creates the map section of a BEAST2 XML parameter file}
\usage{
create_beast2_input_map(beauti_options)
}
\arguments{
\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
lines of XML text
}
\description{
Creates the map section of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rln_clock_model_to_xml_prior_distr.R
\name{rln_clock_model_to_xml_prior_distr}
\alias{rln_clock_model_to_xml_prior_distr}
\title{Internal function}
\usage{
rln_clock_model_to_xml_prior_distr(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
a character vector of XML strings
}
\description{
Internal function to converts a relaxed log-normal clock model
to the \code{prior} section of the XML as text
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_has_non_strict_clock_model.R
\name{get_has_non_strict_clock_model}
\alias{get_has_non_strict_clock_model}
\title{Determines if there is at least one non-strict clock model
in the list of one or more clock models}
\usage{
get_has_non_strict_clock_model(clock_models)
}
\arguments{
\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}
}
\value{
TRUE if there is at least one non-strict clock model
}
\description{
Determines if there is at least one non-strict clock model
in the list of one or more clock models
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_branch_rate_model_xml.R
\name{create_branch_rate_model_xml}
\alias{create_branch_rate_model_xml}
\title{Internal function to create the \code{branchRateModel} section
of the XML as text.}
\usage{
create_branch_rate_model_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
a character vector of XML strings
}
\description{
Creates the \code{branchRateModel} section
of the XML as text.
}
\details{
This function will be called only if there are no MRCA priors.

The \code{distribution} tag (with ID equals \code{treeLikelihood})
has these elements:

\preformatted{
  <branchRateModel[...]>
    [...]
  </branchRateModel>
}

When there is a strict clock,
  \link{create_branch_rate_model_sc_xml} is called.
When there is an RLN clock,
  \link{create_branch_rate_model_rln_xml} is called.

Zooming out:

\preformatted{
  <beast[...]>
    <run[...]>
      <distribution id="posterior"[...]>
        <distribution id="likelihood"[...]>
          <distribution id="treeLikelihood"[...]>
             [...]

             [this section]
          </distribution>
        </distribution>
      </distribution>
    </run>
  </beast>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_screenlog.R
\name{create_screenlog}
\alias{create_screenlog}
\title{Create a \code{screenlog} object}
\usage{
create_screenlog(
  filename = "",
  log_every = 1000,
  mode = "autodetect",
  sanitise_headers = FALSE,
  sort = "none"
)
}
\arguments{
\item{filename}{name of the file to store the posterior screens
phylogenies to. By default, this is \code{$(screen).screens}}

\item{log_every}{number of MCMC states between writing to file}

\item{mode}{mode how to log.
Valid values are the ones returned by \link{get_log_modes}}

\item{sanitise_headers}{set to \link{TRUE} to sanitise the headers of the
log file}

\item{sort}{how to sort the log.
Valid values are the ones returned by \link{get_log_sorts}}
}
\description{
Create a \code{screenlog} object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma_site_model_to_xml_state.R
\name{gamma_site_model_to_xml_state}
\alias{gamma_site_model_to_xml_state}
\title{Converts a gamma site model to XML,
  used in the \code{state} section}
\usage{
gamma_site_model_to_xml_state(gamma_site_model, id)
}
\arguments{
\item{gamma_site_model}{a gamma site model,
as created by \code{\link{create_gamma_site_model}})}

\item{id}{the site model's ID}
}
\value{
the gamma_site model as XML text
}
\description{
Converts a gamma site model to XML,
  used in the \code{state} section
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_clock_rate_state_node_parameter_xml.R
\name{create_clock_rate_state_node_parameter_xml}
\alias{create_clock_rate_state_node_parameter_xml}
\title{Internal function}
\usage{
create_clock_rate_state_node_parameter_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the following XML:
\code{
  <parameter id="ucldStdev.c:[id]" lower="0.0" name="stateNode">
    0.1
  </parameter>
}
}
\description{
Creates the \code{clockRate} parameter with the name \code{stateNode},
such as:
\code{
  <parameter id="ucldStdev.c:[id]" [...] name="stateNode">0.1</parameter>
}
}
\examples{
create_ucld_stdev_state_node_param_xml(
  create_inference_model(
    clock_model = create_rln_clock_model(id = 314)
  )
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_clock_model.R
\name{create_strict_clock_model}
\alias{create_strict_clock_model}
\alias{create_clock_model_strict}
\title{Create a strict clock model}
\usage{
create_strict_clock_model(
  id = NA,
  clock_rate_param = create_clock_rate_param(),
  clock_rate_distr = create_uniform_distr()
)
}
\arguments{
\item{id}{an alignment's IDs.
An ID can be extracted from its FASTA filename
with \code{\link{get_alignment_ids_from_fasta_filenames}})}

\item{clock_rate_param}{the clock rate's parameter,
a numeric value.
For advanced usage, use the structure
as created by the \code{\link{create_clock_rate_param}} function}

\item{clock_rate_distr}{the clock rate's distribution,
as created by a \code{\link{create_distr}} function}
}
\value{
a strict clock_model
}
\description{
Create a strict clock model
}
\examples{
strict_clock_model <- create_strict_clock_model(
  clock_rate_param = 1.0,
  clock_rate_distr = create_uniform_distr()
)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  clock_model = strict_clock_model
)
file.remove(beast2_input_file)

strict_clock_model_gamma <- create_strict_clock_model(
  clock_rate_distr = create_gamma_distr()
)

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  get_fasta_filename(),
  beast2_input_file,
  clock_model = strict_clock_model_gamma
)
file.remove(beast2_input_file)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_tree_prior.R
\name{is_yule_tree_prior}
\alias{is_yule_tree_prior}
\title{Determine if the object is a valid Yule tree prior,}
\usage{
is_yule_tree_prior(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid Yule tree prior}
}
\value{
TRUE if x is a valid Yule tree prior, FALSE otherwise
}
\description{
Determine if the object is a valid Yule tree prior,
}
\examples{
  testit::assert(!is_yule_tree_prior(create_bd_tree_prior()))
  testit::assert(!is_yule_tree_prior(create_cbs_tree_prior()))
  testit::assert(!is_yule_tree_prior(create_ccp_tree_prior()))
  testit::assert(!is_yule_tree_prior(create_cep_tree_prior()))
  testit::assert( is_yule_tree_prior(create_yule_tree_prior()))
}
\seealso{
Use \code{\link{create_yule_tree_prior}} to create a valid
  Yule tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_file.R
\name{create_beast2_input_file}
\alias{create_beast2_input_file}
\title{Create a BEAST2 input file}
\usage{
create_beast2_input_file(
  input_filename,
  output_filename,
  site_model = beautier::create_jc69_site_model(),
  clock_model = beautier::create_strict_clock_model(),
  tree_prior = beautier::create_yule_tree_prior(),
  mrca_prior = NA,
  mcmc = beautier::create_mcmc(),
  beauti_options = beautier::create_beauti_options(),
  tipdates_filename = NA
)
}
\arguments{
\item{input_filename}{A FASTA filename.
Use \code{\link{get_fasta_filename}} to obtain a testing FASTA filename.}

\item{output_filename}{Name of the XML parameter file created by this
function. BEAST2 uses this file as input.}

\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
nothing
}
\description{
Create a BEAST2 input file
}
\examples{
# Get an example FASTA file
input_filename <- get_fasta_filename()

# The file created by beautier, a BEAST2 input file
output_filename <- get_beautier_tempfilename()

create_beast2_input_file(
  input_filename,
  output_filename
)
file.remove(output_filename)
}
\seealso{
Use \link{create_beast2_input_file_from_model} to do the same with an
  inference model.
  See \code{\link{create_site_model}} for examples with
  different site models. See \code{\link{create_clock_model}} for examples
  with clock models. See \code{\link{create_tree_prior}} for examples with
  different tree priors. See \code{\link{create_mcmc}} for examples with
  a different MCMC setup.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_param.R
\name{create_alpha_param}
\alias{create_alpha_param}
\alias{create_param_alpha}
\title{Create a parameter called alpha}
\usage{
create_alpha_param(id = NA, value = 0)
}
\arguments{
\item{id}{the parameter's ID}

\item{value}{value of the parameter}
}
\value{
a parameter called alpha
}
\description{
Create a parameter called alpha
}
\note{
this parameter is used in a beta distribution
  (as returned by \code{\link{create_beta_distr}})
and gamma distribution
  (as returned by \code{\link{create_gamma_distr}})
and inverse-gamma distribution
  (as returned by \code{\link{create_inv_gamma_distr}}).
It cannot be estimated (as a hyper parameter) yet.
}
\examples{
# Create the parameter
alpha_param <- create_alpha_param()

# Use the parameter in a distribution
beta_distr <- create_beta_distr(
  alpha = alpha_param
)

# Use the distribution to create a BEAST2 input file
beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = beta_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_param}} contains a list
  of all parameters that can be created
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_gtr_site_model.R
\name{check_gtr_site_model}
\alias{check_gtr_site_model}
\title{Check if the \code{gtr_site_model} is a valid
GTR nucleotide substitution model.}
\usage{
check_gtr_site_model(gtr_site_model)
}
\arguments{
\item{gtr_site_model}{a GTR site model,
as returned by \code{\link{create_gtr_site_model}}}
}
\description{
Use \link{create_gtr_site_model} to create a valid
GTR nucleotide substitution model.
}
\examples{
check_gtr_site_model(create_gtr_site_model())
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_kappa_1}
\alias{parameter_to_xml_kappa_1}
\title{Internal function}
\usage{
parameter_to_xml_kappa_1(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a kappa 1 parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_kappa_1_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a kappa 1 parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_default_mcmc.R
\name{is_default_mcmc}
\alias{is_default_mcmc}
\title{Determine if the MCMC is a default MCMC}
\usage{
is_default_mcmc(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\value{
TRUE if the MCMC is a default MCMC
}
\description{
Determine if the MCMC is a default MCMC
}
\examples{

# TRUE: An MCMC created by 'create_mcmc' is default.
is_default_mcmc(create_mcmc())

# FALSE: An MCMC created by 'create_ns_mcmc' is not
is_default_mcmc(create_ns_mcmc())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml}
\alias{distr_to_xml}
\title{Internal function}
\usage{
distr_to_xml(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{a distribution,
as created by \code{\link{create_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a distribution to XML
}
\examples{
distr_to_xml(create_uniform_distr(id = 1))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_tree_prior.R
\name{is_cep_tree_prior}
\alias{is_cep_tree_prior}
\title{Determine if the object is a valid
coalescent exponential population tree prior}
\usage{
is_cep_tree_prior(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
constant coalescent exponential population tree prior}
}
\value{
TRUE if x is a valid coalescent exponential population tree prior,
  FALSE otherwise
}
\description{
Determine if the object is a valid
coalescent exponential population tree prior
}
\examples{
  testit::assert(!is_cep_tree_prior(create_bd_tree_prior()))
  testit::assert(!is_cep_tree_prior(create_cbs_tree_prior()))
  testit::assert(!is_cep_tree_prior(create_ccp_tree_prior()))
  testit::assert( is_cep_tree_prior(create_cep_tree_prior()))
  testit::assert(!is_cep_tree_prior(create_yule_tree_prior()))
}
\seealso{
Use \code{\link{create_cep_tree_prior}} to create a valid
  coalescent exponential population tree prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_screenlog.R
\name{check_screenlog_values}
\alias{check_screenlog_values}
\title{Check if the screenlog has the list elements with valid values
for being a valid screenlog object.}
\usage{
check_screenlog_values(screenlog)
}
\arguments{
\item{screenlog}{a \code{screenlog},
as created by \link{create_screenlog}}
}
\value{
nothing
}
\description{
Calls \code{stop} if a value is invalid
}
\seealso{
Use \link{create_screenlog} to create a valid screenlog
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_site_model_n_distrs.R
\name{get_site_model_n_distrs}
\alias{get_site_model_n_distrs}
\title{Get the number of distributions a site model has}
\usage{
get_site_model_n_distrs(site_model)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}
}
\value{
the number of distributions a site model has
}
\description{
Get the number of distributions a site model has
}
\examples{
# 5: rates AC, AG, AT, CG and GT
get_site_model_n_distrs(create_gtr_site_model())

# 1: kappa
get_site_model_n_distrs(create_hky_site_model())

# 0: npne
get_site_model_n_distrs(create_jc69_site_model())

# 2: kappa 1 and kappa 2
get_site_model_n_distrs(create_tn93_site_model())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_log_normal_distr}
\alias{create_log_normal_distr}
\alias{create_distr_log_normal}
\title{Create a log-normal distribution}
\usage{
create_log_normal_distr(
  id = NA,
  m = 0,
  s = 0,
  value = NA,
  lower = NA,
  upper = NA
)
}
\arguments{
\item{id}{the distribution's ID}

\item{m}{the m parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_m_param}}}

\item{s}{the s parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_s_param}}}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
a log-normal distribution
}
\description{
Create a log-normal distribution
}
\examples{
log_normal_distr <- create_log_normal_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = log_normal_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mrca_prior_to_xml_taxonset.R
\name{mrca_prior_to_xml_taxonset}
\alias{mrca_prior_to_xml_taxonset}
\title{Creates the \code{taxonset} section in the prior section of the
distribution section of a BEAST2 XML parameter file.}
\usage{
mrca_prior_to_xml_taxonset(mrca_prior, taxa_names_with_ids = NULL)
}
\arguments{
\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{taxa_names_with_ids}{taxa names that already have received
an ID. Causes the XML to \code{idref} these}
}
\value{
lines of XML text
}
\description{
Creates the \code{taxonset} section in the prior section of the
distribution section of a BEAST2 XML parameter file.
}
\details{
\code{
  <taxonset id="all" spec="TaxonSet">
      <taxon id="626029_aco" spec="Taxon"/>
      <taxon id="630116_aco" spec="Taxon"/>
      <taxon id="630210_aco" spec="Taxon"/>
      <taxon id="B25702_aco" spec="Taxon"/>
      <taxon id="61430_aco" spec="Taxon"/>
  </taxonset>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_temp_tracelog_filename.R
\name{create_temp_tracelog_filename}
\alias{create_temp_tracelog_filename}
\title{Create a filename for a temporary tracelog file}
\usage{
create_temp_tracelog_filename()
}
\description{
Create a filename for a temporary tracelog file
}
\seealso{
use \link{create_tracelog} to create a tracelog.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_trait_set_string.R
\name{create_trait_set_string}
\alias{create_trait_set_string}
\title{Create a trait set string.}
\usage{
create_trait_set_string(df)
}
\arguments{
\item{df}{a data frame with two columns}
}
\value{
the trait set string
}
\description{
For example, a data frame with row \code{A 1}
and another row \code{B 2}, the trait set string will be
\code{A=1,B=2}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_one_div_x_distr}
\alias{is_one_div_x_distr}
\title{Determine if the object is a valid
1/x distribution,
as created by \code{\link{create_one_div_x_distr}}}
\usage{
is_one_div_x_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
1/x distribution}
}
\value{
TRUE if x is a valid 1/x distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
1/x distribution,
as created by \code{\link{create_one_div_x_distr}}
}
\examples{
# TRUE
is_one_div_x_distr(create_one_div_x_distr())
# FALSE
is_one_div_x_distr(create_poisson_distr())
is_one_div_x_distr(NA)
is_one_div_x_distr(NULL)
is_one_div_x_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tree_models_to_xml_tracelog.R
\name{tree_models_to_xml_tracelog}
\alias{tree_models_to_xml_tracelog}
\title{Creates the tree models' XML for the tracelog section}
\usage{
tree_models_to_xml_tracelog(site_models)
}
\arguments{
\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}
}
\value{
lines of XML text
}
\description{
Creates the tree models' XML for the tracelog section
}
\note{
use site_models just because it contains all IDs
}
\examples{
# <logger id="tracelog" ...>
#'   # Here
# </logger>
}
\seealso{
the complete tracelog section is created
  by \code{\link{create_tracelog_xml}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_poisson_distr}
\alias{create_poisson_distr}
\alias{create_distr_poisson}
\title{Create a Poisson distribution}
\usage{
create_poisson_distr(id = NA, lambda = 0, value = NA, lower = NA, upper = NA)
}
\arguments{
\item{id}{the distribution's ID}

\item{lambda}{the lambda parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_lambda_param}}}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
a Poisson distribution
}
\description{
Create a Poisson distribution
}
\examples{
poisson_distr <- create_poisson_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = poisson_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parameter_to_xml.R
\name{parameter_to_xml_mean}
\alias{parameter_to_xml_mean}
\title{Internal function}
\usage{
parameter_to_xml_mean(parameter, beauti_options = create_beauti_options())
}
\arguments{
\item{parameter}{a mean parameter,
a numeric value.
For advanced usage, use the structure
as created by \code{\link{create_mean_param}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the parameter as XML text
}
\description{
Converts a mean parameter to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_first_xml_opening_tag_line.R
\name{find_first_xml_opening_tag_line}
\alias{find_first_xml_opening_tag_line}
\title{Find the line number of the first section's opening tag}
\usage{
find_first_xml_opening_tag_line(lines, section)
}
\arguments{
\item{lines}{the lines of an XML text}

\item{section}{the name of the XML section}
}
\value{
the line number's index (which is 1 for the first line) if the
  opening tag is found, else NA
}
\description{
Find the line number of the first section's opening tag
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_model.R
\name{create_site_model}
\alias{create_site_model}
\title{General function to create a site model.}
\usage{
create_site_model(name, id, gamma_site_model = create_gamma_site_model(), ...)
}
\arguments{
\item{name}{the site model name. Valid
names can be found in \code{get_site_model_names}}

\item{id}{the IDs of the alignment (can be extracted from
the FASTA filename using \code{\link{get_alignment_id}})}

\item{gamma_site_model}{a gamma site model, as created
by \code{\link{create_gamma_site_model}}}

\item{...}{specific site model parameters}
}
\value{
a site_model
}
\description{
General function to create a site model.
}
\note{
Prefer using the
  named functions
  \code{\link{create_gtr_site_model}},
  \code{\link{create_hky_site_model}},,
  \code{\link{create_jc69_site_model}},
  and \code{\link{create_tn93_site_model}}
}
\examples{
input_filename <- get_fasta_filename()

# GTR
output_filename <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = input_filename,
  output_filename = output_filename,
  site_model = create_gtr_site_model()
)
file.remove(output_filename)

# HKY
output_filename <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = input_filename,
  output_filename = output_filename,
  site_model = create_hky_site_model()
)
file.remove(output_filename)

# JC69
output_filename <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = input_filename,
  output_filename = output_filename,
  site_model = create_jc69_site_model()
)
file.remove(output_filename)

# TN93
output_filename <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = input_filename,
  output_filename = output_filename,
  site_model = create_tn93_site_model()
)
file.remove(output_filename)
}
\seealso{
See \code{\link{create_gtr_site_model}} for more examples
  with a GTR site model. See \code{\link{create_hky_site_model}}
  for more examples with an HKY site model. See
  \code{\link{create_jc69_site_model}} for more examples with a JC69
  site model. See \code{\link{create_tn93_site_model}} for more
  examples with a TN93 site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_tn93_site_model.R
\name{check_tn93_site_model_names}
\alias{check_tn93_site_model_names}
\title{Check if the \code{tn93_site_model} has the list elements
of a valid \code{tn93_site_model} object.}
\usage{
check_tn93_site_model_names(tn93_site_model)
}
\arguments{
\item{tn93_site_model}{a TN93 site model,
as returned by \code{\link{create_tn93_site_model}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element is missing
}
\seealso{
Use \link{create_tn93_site_model}
to create a valid \code{tn93_site_model}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_mcmc_filenames.R
\name{get_mcmc_filenames}
\alias{get_mcmc_filenames}
\title{Get the filenames stored in an MCMC.}
\usage{
get_mcmc_filenames(mcmc)
}
\arguments{
\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}
}
\description{
If a filename is set to an empty string, to indicate a certain log file
need not be created, this (non-)filename will not be returned.
}
\examples{

mcmc <- create_mcmc()
mcmc$tracelog$filename <- "/home/john/trace.log"
mcmc$screenlog$filename <- "/home/john/screen.log"
mcmc$treelog$filename <- "/home/john/tree.log"

# 3 filenames
filenames <- get_mcmc_filenames(mcmc)

# If there is no need to write to the screenlog file ...
mcmc$screenlog$filename <- ""

# 2 filenames
# ... one file less will be created
filenames <- get_mcmc_filenames(mcmc)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_gamma_site_model.R
\name{is_gamma_site_model}
\alias{is_gamma_site_model}
\title{Is object x a gamma site model?}
\usage{
is_gamma_site_model(x)
}
\arguments{
\item{x}{the object to be determined if it is a valid gamma site object}
}
\value{
TRUE if x is a valid gamma site object, FALSE otherwise
}
\description{
Is object x a gamma site model?
}
\examples{
# TRUE
is_gamma_site_model(create_gamma_site_model())

# FALSE
is_gamma_site_model("nonsense")
is_gamma_site_model(NA)
is_gamma_site_model(NULL)
is_gamma_site_model("")
is_gamma_site_model(c())
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_normal_distr}
\alias{init_normal_distr}
\title{Initializes an normal distribution}
\usage{
init_normal_distr(normal_distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{normal_distr}{a normal distribution,
using \link{create_normal_distr}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized normal distribution
}
\description{
Initializes an normal distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_model_parameters_xml.R
\name{create_site_model_parameters_xml}
\alias{create_site_model_parameters_xml}
\title{Internal function to creates the XML text for the
\code{parameter}s within the \code{siteModel} section
of a BEAST2 parameter file.}
\usage{
create_site_model_parameters_xml(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
the site model as XML text
}
\description{
Internal function to creates the XML text for the
\code{parameter}s within the \code{siteModel} section,
which is part of the \code{siteModel} section
of a BEAST2 parameter file.
}
\details{
The \code{parameter}s sections has these elements:

\preformatted{
   [parameters]
}

\code{[parameters]} can be a combination of these:

\preformatted{
  <parameter id="mutationRate.s[...]>
  <parameter id="gammaShape.s[...]>
  <parameter id="proportionInvariant.s[...]>
}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_distr.R
\name{is_uniform_distr}
\alias{is_uniform_distr}
\title{Determine if the object is a valid
uniform distribution
as created by \code{\link{create_uniform_distr}}}
\usage{
is_uniform_distr(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
uniform distribution}
}
\value{
TRUE if x is a valid uniform distribution,
  FALSE otherwise
}
\description{
Determine if the object is a valid
uniform distribution
as created by \code{\link{create_uniform_distr}}
}
\examples{
# TRUE
is_uniform_distr(create_uniform_distr())
# FALSE
is_uniform_distr(create_beta_distr())
is_uniform_distr(NA)
is_uniform_distr(NULL)
is_uniform_distr("nonsense")
}
\seealso{
use \code{\link{is_distr}} to see if x is any
  distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma_site_models_to_xml_prior_distr.R
\name{gamma_site_models_to_xml_prior_distr}
\alias{gamma_site_models_to_xml_prior_distr}
\title{Creates the gamma site models section in the distribution section
of a BEAST2 XML parameter file}
\usage{
gamma_site_models_to_xml_prior_distr(site_models)
}
\arguments{
\item{site_models}{one or more site models,
as returned by \code{\link{create_site_model}}}
}
\value{
lines of XML text
}
\description{
Creates the gamma site models section in the distribution section
of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_equal_xml_lines.R
\name{are_equal_xml_lines}
\alias{are_equal_xml_lines}
\title{Determine if XML lines result in equal trees}
\usage{
are_equal_xml_lines(lines_1, lines_2, section)
}
\arguments{
\item{lines_1}{lines of a first XML file}

\item{lines_2}{lines of a second XML file}

\item{section}{name of an XML section.
Assumes that there is one line that starts with \code{<section}
(excluding whitespace)
and one line that is \code{</section>} (also excluding whitespace)}
}
\value{
TRUE if the two sections of the XML files are equal,
  FALSE otherwise
}
\description{
Determine if XML lines result in equal trees
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_exp_distr}
\alias{init_exp_distr}
\title{Initializes an exponential distribution}
\usage{
init_exp_distr(exp_distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{exp_distr}{a exponential distribution,
using \link{create_exp_distr}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized exponential distribution
}
\description{
Initializes an exponential distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_log_normal_distr}
\alias{is_init_log_normal_distr}
\title{Determine if x is an initialized log_normal distribution object
  as created by \code{\link{create_log_normal_distr}}}
\usage{
is_init_log_normal_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized log_normal distribution object}
}
\value{
TRUE if x is an initialized log_normal distribution object
}
\description{
Determine if x is an initialized log_normal distribution object
  as created by \code{\link{create_log_normal_distr}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_distr.R
\name{init_log_normal_distr}
\alias{init_log_normal_distr}
\title{Initializes an log-normal distribution}
\usage{
init_log_normal_distr(log_normal_distr, distr_id = 0, param_id = 0)
}
\arguments{
\item{log_normal_distr}{a log-normal distribution,
using \link{create_log_normal_distr}}

\item{distr_id}{the first distribution's ID}

\item{param_id}{the first parameter's ID}
}
\value{
an initialized log-normal distribution
}
\description{
Initializes an log-normal distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_site_models.R
\name{create_site_models}
\alias{create_site_models}
\title{Creates all supported site models
  which is a list of the types returned by
  \code{\link{create_gtr_site_model}},
  \code{\link{create_hky_site_model}},
  \code{\link{create_jc69_site_model}}
  and \code{\link{create_tn93_site_model}}}
\usage{
create_site_models()
}
\value{
a list of site_models
}
\description{
Creates all supported site models
  which is a list of the types returned by
  \code{\link{create_gtr_site_model}},
  \code{\link{create_hky_site_model}},
  \code{\link{create_jc69_site_model}}
  and \code{\link{create_tn93_site_model}}
}
\examples{
# All created site models are a kind of site model
site_models <- create_site_models()

# TRUE
is_gtr_site_model(site_models[[1]])
is_hky_site_model(site_models[[2]])
is_jc69_site_model(site_models[[3]])
is_tn93_site_model(site_models[[4]])
}
\seealso{
Use \link{create_site_model} to create a site model
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clock_models_to_xml_prior_distr.R
\name{clock_models_to_xml_prior_distr}
\alias{clock_models_to_xml_prior_distr}
\title{Deprecated function}
\usage{
clock_models_to_xml_prior_distr(
  clock_models = "deprecated",
  mrca_priors = "deprecated",
  tipdates_filename = "deprecated"
)
}
\arguments{
\item{clock_models}{a list of one or more clock models,
as returned by \code{\link{create_clock_model}}}

\item{mrca_priors}{a list of one or more Most Recent Common Ancestor priors,
as returned by \code{\link{create_mrca_prior}}}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
a character vector of XML strings
}
\description{
Internal function to represent the clock models as XML
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_distr}
\alias{create_distr}
\title{General function to create a distribution.}
\usage{
create_distr(name, id, value = NA, lower = NA, upper = NA, ...)
}
\arguments{
\item{name}{the distribution name. Valid
names can be found in \code{get_distr_names}}

\item{id}{the distribution's ID}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}

\item{...}{specific distribution parameters}
}
\value{
a distribution
}
\description{
General function to create a distribution.
}
\note{
Prefer using the
  named functions
  \code{\link{create_beta_distr}},
  \code{\link{create_exp_distr}},
  \code{\link{create_gamma_distr}},
  \code{\link{create_inv_gamma_distr}},
  \code{\link{create_laplace_distr}},
  \code{\link{create_log_normal_distr}},
  \code{\link{create_normal_distr}},
  \code{\link{create_one_div_x_distr}},
  \code{\link{create_poisson_distr}}
  and \code{\link{create_uniform_distr}}

See
  \code{\link{create_beta_distr}},
  \code{\link{create_exp_distr}},
  \code{\link{create_gamma_distr}},
  \code{\link{create_inv_gamma_distr}},
  \code{\link{create_laplace_distr}},
  \code{\link{create_log_normal_distr}},
  \code{\link{create_normal_distr}},
  \code{\link{create_one_div_x_distr}},
  \code{\link{create_poisson_distr}}
  and \code{\link{create_uniform_distr}}
  for examples how to use those distributions
}
\examples{
# Use any distribution
distr <- create_beta_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = distr
  )
)
file.remove(beast2_input_file)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_mrca_prior_with_distr.R
\name{is_mrca_prior_with_distr}
\alias{is_mrca_prior_with_distr}
\title{See if x is one MRCA prior with a distribution}
\usage{
is_mrca_prior_with_distr(x)
}
\arguments{
\item{x}{the object to be tested}
}
\value{
TRUE if x is one MRCA prior with a distribution,
  FALSE otherwise
}
\description{
See if x is one MRCA prior with a distribution
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/are_init_mrca_priors.R
\name{are_init_mrca_priors}
\alias{are_init_mrca_priors}
\title{Determine if x consists out of initialized MRCA priors}
\usage{
are_init_mrca_priors(x)
}
\arguments{
\item{x}{the object to check if it consists out of
initialized MRCA priors}
}
\value{
TRUE if x, or all elements of x, are initialized MRCA priors
}
\description{
Determine if x consists out of initialized MRCA priors
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_tree_priors.R
\name{create_tree_priors}
\alias{create_tree_priors}
\title{Creates all supported tree priors,
  which is a list of the types returned by
  \code{\link{create_bd_tree_prior}},
  \code{\link{create_cbs_tree_prior}},
  \code{\link{create_ccp_tree_prior}},
  \code{\link{create_cep_tree_prior}}
  and \code{\link{create_yule_tree_prior}}}
\usage{
create_tree_priors()
}
\value{
a list of tree_priors
}
\description{
Creates all supported tree priors,
  which is a list of the types returned by
  \code{\link{create_bd_tree_prior}},
  \code{\link{create_cbs_tree_prior}},
  \code{\link{create_ccp_tree_prior}},
  \code{\link{create_cep_tree_prior}}
  and \code{\link{create_yule_tree_prior}}
}
\examples{
tree_priors <- create_tree_priors()
# TRUE
is_bd_tree_prior(tree_priors[[1]])
is_cbs_tree_prior(tree_priors[[2]])
is_ccp_tree_prior(tree_priors[[3]])
is_cep_tree_prior(tree_priors[[4]])
is_yule_tree_prior(tree_priors[[5]])
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/gamma_distr_to_xml.R
\name{gamma_distr_to_xml}
\alias{gamma_distr_to_xml}
\title{Internal function}
\usage{
gamma_distr_to_xml(gamma_distr, beauti_options = create_beauti_options())
}
\arguments{
\item{gamma_distr}{a gamma distribution,
as created by \code{\link{create_gamma_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a gamma distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_mrca_prior.R
\name{check_mrca_prior}
\alias{check_mrca_prior}
\title{Check if the MRCA prior is a valid MRCA prior.}
\usage{
check_mrca_prior(mrca_prior)
}
\arguments{
\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}
}
\value{
nothing
}
\description{
Calls \code{stop} if the MRCA prior is invalid.
}
\examples{
fasta_filename <- get_beautier_path("anthus_aco.fas")
mrca_prior <- create_mrca_prior(
  alignment_id = get_alignment_id(fasta_filename = fasta_filename),
  taxa_names = get_taxa_names(filename = fasta_filename)
)
mrca_prior <- create_mrca_prior(
 alignment_id = get_alignment_id(fasta_filename = fasta_filename),
 taxa_names = get_taxa_names(filename = fasta_filename)
)
check_mrca_prior(mrca_prior)
}
\seealso{
Use \link{create_mrca_prior} to create a valid MRCA prior
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_inference_model.R
\name{create_inference_model}
\alias{create_inference_model}
\title{Create a Bayesian phylogenetic inference model.}
\usage{
create_inference_model(
  site_model = beautier::create_jc69_site_model(),
  clock_model = beautier::create_strict_clock_model(),
  tree_prior = beautier::create_yule_tree_prior(),
  mrca_prior = NA,
  mcmc = beautier::create_mcmc(),
  beauti_options = beautier::create_beauti_options(),
  tipdates_filename = NA
)
}
\arguments{
\item{site_model}{a site model,
as returned by \code{\link{create_site_model}}}

\item{clock_model}{a clock model,
as returned by \code{\link{create_clock_model}}}

\item{tree_prior}{a tree priors,
as returned by \code{\link{create_tree_prior}}}

\item{mrca_prior}{a Most Recent Common Ancestor prior,
as returned by \code{\link{create_mrca_prior}}}

\item{mcmc}{one MCMC.
Use \code{\link{create_mcmc}} to create an MCMC.
Use \code{\link{create_ns_mcmc}} to create an MCMC
  for a Nested Sampling run.
Use \code{\link{check_mcmc}} to check if an MCMC is valid.
Use \code{\link{rename_mcmc_filenames}} to rename the filenames in an MCMC.}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}

\item{tipdates_filename}{name of the file containing the tip dates.
This file is assumed to have two columns, separated by a tab.
The first column contains the taxa names, the second column contains
the date.}
}
\value{
an inference model
}
\description{
Create a Bayesian phylogenetic inference model,
as can be done by BEAUti.
}
\examples{
# Create an MCMC chain with 50 states
inference_model <- create_inference_model(
  mcmc = create_mcmc(chain_length = 50000, store_every = 1000)
)

output_filename <- get_beautier_tempfilename()
create_beast2_input_file_from_model(
  input_filename = get_fasta_filename(),
  output_filename = output_filename,
  inference_model = inference_model
)
file.remove(output_filename)
}
\seealso{
Use \link{create_test_inference_model} to create an inference model
with a short MCMC, to be used in testing.
Use \link{create_ns_inference_model} to create an inference model
to estimate the marginal likelihood
(aka evidence) using a nested sampling approach.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_init.R
\name{create_beast2_input_init}
\alias{create_beast2_input_init}
\title{Creates the \code{init} section of a BEAST2 XML parameter file}
\usage{
create_beast2_input_init(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Creates the \code{init} section of a BEAST2 XML parameter file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_init_distr.R
\name{is_init_gamma_distr}
\alias{is_init_gamma_distr}
\title{Determine if x is an initialized gamma distribution object}
\usage{
is_init_gamma_distr(x)
}
\arguments{
\item{x}{the object to check if it is an
initialized gamma distribution object}
}
\value{
TRUE if x is an initialized gamma distribution object
}
\description{
Determine if x is an initialized gamma distribution object
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/init_site_models.R
\name{init_hky_site_model}
\alias{init_hky_site_model}
\title{Initializes an HKY site model}
\usage{
init_hky_site_model(hky_site_model, distr_id = 0, param_id = 0)
}
\arguments{
\item{hky_site_model}{an HKY site model,
as returned by \code{\link{create_hky_site_model}}}

\item{distr_id}{a distributions' ID}

\item{param_id}{a parameter's ID}
}
\value{
an initialized HKY site model
}
\description{
Initializes an HKY site model
}
\examples{

hky_site_model <- create_hky_site_model()
is_init_hky_site_model(hky_site_model)
hky_site_model <- init_hky_site_model(hky_site_model)
is_init_hky_site_model(hky_site_model)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_distr.R
\name{create_laplace_distr}
\alias{create_laplace_distr}
\alias{create_distr_laplace}
\title{Create a Laplace distribution}
\usage{
create_laplace_distr(
  id = NA,
  mu = 0,
  scale = 1,
  value = NA,
  lower = NA,
  upper = NA
)
}
\arguments{
\item{id}{the distribution's ID}

\item{mu}{the mu parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_mu_param}}}

\item{scale}{the scale parameter,
a numeric value.
For advanced usage, use the structure
as returned by \code{\link{create_scale_param}}}

\item{value}{the initial value for the MCMC}

\item{lower}{the lower bound, the lowest possible value}

\item{upper}{an upper limit of the uniform distribution.
If the upper limits needs to be infinity, set \code{upper} to \code{Inf}.}
}
\value{
a Laplace distribution
}
\description{
Create a Laplace distribution
}
\examples{
laplace_distr <- create_laplace_distr()

beast2_input_file <- get_beautier_tempfilename()
create_beast2_input_file(
  input_filename = get_fasta_filename(),
  beast2_input_file,
  tree_prior = create_yule_tree_prior(
    birth_rate_distr = laplace_distr
  )
)
file.remove(beast2_input_file)
}
\seealso{
the function \code{\link{create_distr}} shows an overview
  of all supported distributions
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/create_beast2_input_distr.R
\name{create_beast2_input_distr}
\alias{create_beast2_input_distr}
\title{Creates the distribution section of a BEAST2 XML parameter file.}
\usage{
create_beast2_input_distr(inference_model)
}
\arguments{
\item{inference_model}{a Bayesian phylogenetic inference model.
An inference model is the complete model setup in which a site model,
clock model, tree prior and more are specified.
Use \link{create_inference_model} to create an inference model.
Use \link{check_inference_model} to check if  an inference model is valid.
Use \link{rename_inference_model_filenames} to rename the files in an
inference model.}
}
\value{
lines of XML text
}
\description{
Creates the distribution section of a BEAST2 XML parameter file.
}
\note{
this function is not intended for regular use, thus its
  long name length is accepted
}
\examples{
 # <distribution id="posterior" spec="util.CompoundDistribution">
 #     <distribution id="prior" spec="util.CompoundDistribution">
 #       HERE, where the ID of the distribution is 'prior'
 #     </distribution>
 #     <distribution id="likelihood" ...>
 #     </distribution>
 # </distribution>
}
\seealso{
\code{\link{create_beast2_input}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_param.R
\name{is_kappa_1_param}
\alias{is_kappa_1_param}
\title{Determine if the object is a valid kappa 1 parameter}
\usage{
is_kappa_1_param(x)
}
\arguments{
\item{x}{an object, to be determined if it is a valid
kappa 1 parameter}
}
\value{
TRUE if x is a valid kappa 1 parameter,
  FALSE otherwise
}
\description{
Determine if the object is a valid kappa 1 parameter
}
\examples{

is_kappa_1_param(create_alpha_param())
is_kappa_1_param(create_beta_param())
is_kappa_1_param(create_clock_rate_param())
is_kappa_1_param(create_kappa_1_param())
is_kappa_1_param(create_kappa_2_param())
is_kappa_1_param(create_lambda_param())
is_kappa_1_param(create_m_param())
is_kappa_1_param(create_mean_param())
is_kappa_1_param(create_mu_param())
is_kappa_1_param(create_rate_ac_param())
is_kappa_1_param(create_rate_ag_param())
is_kappa_1_param(create_rate_at_param())
is_kappa_1_param(create_rate_cg_param())
is_kappa_1_param(create_rate_ct_param())
is_kappa_1_param(create_rate_gt_param())
is_kappa_1_param(create_s_param())
is_kappa_1_param(create_scale_param())
is_kappa_1_param(create_sigma_param())

is_kappa_1_param(NA)
is_kappa_1_param(NULL)
is_kappa_1_param("nonsense")
is_kappa_1_param(create_jc69_site_model())
is_kappa_1_param(create_strict_clock_model())
is_kappa_1_param(create_yule_tree_prior())
is_kappa_1_param(create_mcmc())
}
\seealso{
kappa 1 parameters are returned by
  \code{\link{create_kappa_1_param}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_tracelog.R
\name{check_tracelog_names}
\alias{check_tracelog_names}
\title{Check if the \code{tracelog} has the list elements
of a valid \code{tracelog} object.}
\usage{
check_tracelog_names(tracelog)
}
\arguments{
\item{tracelog}{a \code{tracelog},
as created by \link{create_tracelog}}
}
\value{
nothing
}
\description{
Calls \code{stop} if an element is missing
}
\seealso{
Use \link{create_tracelog} to create a valid \code{tracelog}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distr_to_xml.R
\name{distr_to_xml_beta}
\alias{distr_to_xml_beta}
\title{Internal function}
\usage{
distr_to_xml_beta(distr, beauti_options = create_beauti_options())
}
\arguments{
\item{distr}{a beta distribution,
as created by \code{\link{create_beta_distr}})}

\item{beauti_options}{one BEAUti options object,
as returned by \code{\link{create_beauti_options}}}
}
\value{
the distribution as XML text
}
\description{
Converts a beta distribution to XML
}
\author{
Richèl J.C. Bilderbeek
}
