# tracerer

[![Peer Review Status](https://badges.ropensci.org/209_status.svg)](https://github.com/ropensci/onboarding/issues/209)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/tracerer)](https://cran.r-project.org/package=tracerer)
[![](http://cranlogs.r-pkg.org/badges/grand-total/tracerer)]( https://CRAN.R-project.org/package=tracerer)
[![](http://cranlogs.r-pkg.org/badges/tracerer)](https://CRAN.R-project.org/package=tracerer)
[![DOI](https://zenodo.org/badge/114987588.svg)](https://zenodo.org/badge/latestdoi/114987588)

Branch   |[![GitHub Actions logo](man/figures/GitHubActions.png)](https://github.com/ropensci/tracerer/actions)|[![Travis CI logo](man/figures/TravisCI.png)](https://travis-ci.com)                                                  |[![Codecov logo](man/figures/Codecov.png)](https://www.codecov.io)
---------|-----------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------|----------------------------------------------------------------------------------------------------------------------------------------------------
`master` |![R-CMD-check](https://github.com/ropensci/tracerer/workflows/R-CMD-check/badge.svg?branch=master)   |[![Build Status](https://travis-ci.com/ropensci/tracerer.svg?branch=master)](https://travis-ci.com/ropensci/tracerer) |[![codecov.io](https://codecov.io/github/ropensci/tracerer/coverage.svg?branch=master)](https://codecov.io/github/ropensci/tracerer/branch/master)
`develop`|![R-CMD-check](https://github.com/ropensci/tracerer/workflows/R-CMD-check/badge.svg?branch=develop)  |[![Build Status](https://travis-ci.com/ropensci/tracerer.svg?branch=develop)](https://travis-ci.com/ropensci/tracerer)|[![codecov.io](https://codecov.io/github/ropensci/tracerer/coverage.svg?branch=develop)](https://codecov.io/github/ropensci/tracerer/branch/develop)

`tracerer`, 'Tracer for R' is an R package to work with BEAST2 output files. 

![tracerer logo](man/figures/tracerer_logo.png)

`tracerer` is part of the [`babette`](https://github.com/ropensci/babette) package suite:

 * [`beautier`](https://github.com/ropensci/beautier) creates BEAST2 input (`.xml`) files.
 * [`beastier`](https://github.com/ropensci/beastier) runs BEAST2
 * [`mauricer`](https://github.com/ropensci/mauricer): install BEAST2 packages
 * [`tracerer`](https://github.com/ropensci/tracerer) pastes BEAST2 output (`.log`, `.trees`, etc) files.

Related R packages:

 * [`lumier`](https://github.com/ropensci/lumier): Shiny app to help create the function call needed
 * [`BEASTmasteR`](https://github.com/nmatzke/BEASTmasteR): tip-dating using fossils as dated terminal taxa
 * [`RBeast`](https://github.com/beast-dev/RBeast): misc other things
 * [`tracerer_on_windows`](https://github.com/richelbilderbeek/tracerer_on_windows): verifies `tracerer` builds on Windows

## Example

```{r}
library(tracerer)

# Obtain an example log file its name
filename <- get_tracerer_path("beast2_example_output.log")

# Parse that log file
beast_log_full <- parse_beast_tracelog_file(filename)

# Remove the burn-in
beast_log <- remove_burn_ins(
  beast_log_full,
  burn_in_fraction = 0.1
)

# Calculates the effective sample sizes of all parameter estimates
esses <- calc_esses(beast_log, sample_interval = 1000)
```

## Installation

You can install:

 * (recommended) The CRAN version
 * The stable development version
 * The bleeding edge development version

### CRAN

`tracerer` is on CRAN:

```{r}
install.packages("tracerer")
```

### Stable development version

Install the `tracerer` `master` branch using `remotes`:

```{r}
remotes::install_github("ropensci/tracerer")
```

### Bleeding edge development version

Install the `tracerer` `develop` branch using `remotes`:

```{r}
remotes::install_github("ropensci/tracerer", ref = "develop")
```

## FAQ

See [FAQ](faq.md)

## There is a feature I miss

See [CONTRIBUTING](CONTRIBUTING.md), at `Submitting use cases`

## I want to collaborate

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting code'

## I think I have found a bug

See [CONTRIBUTING](CONTRIBUTING.md), at 'Submitting bugs' 

## There's something else I want to say

Sure, just add an Issue. Or send an email.

## External links

 * [BEAST2 GitHub repository](https://github.com/CompEvol/beast2)
 * [Tracer GitHub repository](https://github.com/beast-dev/tracer)

## References

Article about `babette`:

 * Bilderbeek, Richèl JC, and Rampal S. Etienne. "`babette`: BEAUti 2, BEAST 2 and Tracer for R." Methods in Ecology and Evolution (2018). https://doi.org/10.1111/2041-210X.13032

FASTA files `anthus_aco.fas` and `anthus_nd2.fas` from:
 
 * Van Els, Paul, and Heraldo V. Norambuena. "A revision of species limits in Neotropical pipits Anthus based on multilocus genetic and vocal data." Ibis.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# News

Newest versions at top.

## tracerer 2.2.2 (2021-05-30)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * No `LazyData` in DESCRIPTION

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## tracerer 2.2.1 (2021-05-29)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * `save_beast_estimates` can save in sub-sub-folder
  * All temporary files are cleaned up

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## tracerer 2.2 (2021-05-14)

### NEW FEATURES

  * Use the GitHub Actions continuous integration service

### MINOR IMPROVEMENTS

  * Update documentation

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * `parse_beast_log` will be deprecated and gives a warning. 
    Use `parse_beast_tracelog_file` instead"

## tracerer 2.1 (2020-04-27)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Use `message` instead of `print`
  * Use https of BEAST2 website

### BUG FIXES

  * `calc_summary_stats` shows median. Thanks to Hongkai Zhang, @HKyleZhang

### DEPRECATED AND DEFUNCT

  * None

## tracerer 2.0.4 (2020-01-23)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Removed needless dependency on `geiger`

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## tracerer 2.0.3 (2020-01-06)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Documentation at [rOpenSci](https://docs.ropensci.org/tracerer) to
    shows pictures

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## tracerer 2.0.2 (2019-12-02)

### NEW FEATURES

  * Increased error checking in saving and loading `.trees` files
  * Add vignettes

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## tracerer 2.0.1 (2019-03-08)

### NEW FEATURES

  * `tracerer` has passed rOpenSci peer review
  * `tracerer` is on CRAN

### MINOR IMPROVEMENTS

  * None

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## tracerer 1.5.2 (2018-11-01)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tested to work on macOS
  * Tests to not create local temporary files

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## tracerer 1.4.3 (2018-05-17)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Tagged for [the academic article about `babette`](https://github.com/ropensci/babette_article)

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None

## tracerer 1.4.2 (2018-04-05)

### NEW FEATURES

  * None

### MINOR IMPROVEMENTS

  * Follow all [rOpenSci packaging guidelines](https://github.com/ropensci/onboarding/blob/master/packaging_guide.md)

### BUG FIXES

  * None

### DEPRECATED AND DEFUNCT

  * None
# FAQ

## Which version of BEAUti do you use as a guideline?

Version 2.5.0, as can be found in the [install_beast2](https://github.com/ropensci/beastier/blob/master/R/install_beast2.R) function.

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

 * Bilderbeek, Richèl J.C., Etienne, Rampal S., "babette: BEAUti 2, BEAST2 and Tracer for R". bioRxiv 271866; doi: https://doi.org/10.1101/271866

## What is the idea behind the logo?

The logo consists of a rough redraw of Tracer, from the Blizzard
game 'Overwatch', and the R logo. 

## What are the FASTA files?

FASTA files `anthus_aco.fas` and `anthus_nd2.fas` from:
 
 * Van Els, Paul, and Heraldo V. Norambuena. "A revision of species limits in Neotropical pipits Anthus based on multilocus genetic and vocal data." Ibis.

Thanks to Paul van Els.

## Why the logo?

Initially, the logo was a low-tech remake of Tracer, from the game Overwatch by Blizzard. 
To prevent problems with Blizzard, a different logo was picked.

The current logo shows an ant, an animal that leaves a trace of pheromones.
The any is drawn by Jose Scholte, who kindly allowed her work to
be used for free, by attribution.
# Contributing

Awesome that you are reading this.

This GitHub follows the [Contributor Covenant Code of Conduct](code_of_conduct.md).

 * For questions, you can create an Issue
 * Code changes go via Pull Requests

## Submitting use cases

New use cases are welcome.

You can do so by:

 * Add an Issue
 * Send @richelbilderbeek an email (@richelbilderbeek will make an Issue of it)

## Which package to contribute to?

`tracerer` is part of the `babette` package suite,
which consists out of five packages.
Here is how to determine which package is best suited for your contribution:

If you want to contribute to the creation of BEAST2 XML input file,
go to [beautier](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md).

If you want to contribute to how BEAST2 is run,
go to [beautier](https://github.com/ropensci/beautier/blob/master/CONTRIBUTING.md).

If you want to contribute regarding the BEAST2 package management,
go to [mauricer](https://github.com/ropensci/mauricer/blob/master/CONTRIBUTING.md)

If you want to contribute with an overarching idea,
go to [babette](https://github.com/ropensci/babette/blob/master/CONTRIBUTING.md)

If you want to contribute to how BEAST2 output is parsed, you are at the right spot :-) 

## Submitting code

Submitted code should follow these quality guidelines:

 * All tests pass cleanly/silently
 * Code coverage above 95%
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

To get started working on `tracerer` do:

```
git clone https://github.com/ropensci/tracerer
```

Development is done on the `develop` branch. 
To download and checkout the `develop` branch, 
first go into the `tracerer` folder (`cd tracerer`), then do:

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

`tracerer` pictures.

## How did you convert the fuzzy white background to one single color?

```
convert ant.png -fuzz 15% -fill white -opaque white ant_mono_background.png
convert ant_mono_background.png -background white -alpha remove ant_mono_background_2.png
```