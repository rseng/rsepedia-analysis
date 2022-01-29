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

