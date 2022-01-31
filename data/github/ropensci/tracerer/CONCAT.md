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
```---
title: "tracerer versus Tracer demo"
author: "Richèl J.C. Bilderbeek"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{tracerer versus Tracer demo}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

![](tracerer_logo.png)

`tracerer`: 'Tracer for R' is an R package that does the same as
Tracer does, from within R. 

To use `tracerer`, it needs to be loaded:

```{r}
library(tracerer)
```

When loading `beast2_example_output.log` in Tracer, the following is displayed:

![Tracer output](tracer_example_output.png)

Most prominently, at the left, the effective sample sizes (ESSes) are shown.

The show the ESSes using `tracerer`:

```{r}
estimates <- parse_beast_tracelog_file(
  get_tracerer_path("beast2_example_output.log")
)
estimates <- remove_burn_ins(estimates, burn_in_fraction = 0.1)
esses <- calc_esses(estimates, sample_interval = 1000)
table <- t(esses)
colnames(table) <- c("ESS")
knitr::kable(table)
```

At the top-right, some measures of the variable `posterior` is shown.
To reproduce these measures in `tracerer`:

```{r}
sum_stats <- calc_summary_stats(
  estimates$posterior,
  sample_interval = 1000
)
table <- t(sum_stats)
colnames(table) <- c("sum_stat")
knitr::kable(table)
```

Unlike Tracer, in `tracerer` all summary statistics can be obtained at once:

```{r}
sum_stats <- calc_summary_stats(
  estimates,
  sample_interval = 1000
)
knitr::kable(sum_stats)
```

At the bottom-right, a histogram of the posterior estimates is shown.
To reproduce these measures in `tracerer`:

```{r}
ggplot2::ggplot(
  data = remove_burn_ins(estimates, burn_in_fraction = 0.1),
  ggplot2::aes(posterior)
) + ggplot2::geom_histogram(binwidth = 0.21) +
  ggplot2::scale_x_continuous(breaks = seq(-75, -68))

```

Tracer can also show the trace of each estimated variable:

![Tracer shows the trace of the posterior likelihood](tracer_trace_posterior.png)

Same can be done with `tracerer`:

```{r}
ggplot2::ggplot(
  data = remove_burn_ins(estimates, burn_in_fraction = 0.1),
  ggplot2::aes(x = Sample)
) + ggplot2::geom_line(ggplot2::aes(y = posterior))

```

`tracerer` can also use part of `DensiTree`'s functionality.
Here is `beast2_example_output.trees` displayed by `DensiTree`:

![DensiTree output](densitree_example_output.png)

The same is achieved in `tracerer` with:

```{r fig.width=7, fig.height=7}
trees <- parse_beast_trees(
  get_tracerer_path("beast2_example_output.trees")
)
phangorn::densiTree(trees, width = 2)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_posterior.R
\name{is_posterior}
\alias{is_posterior}
\title{Determines if the input is a BEAST2 posterior}
\usage{
is_posterior(x)
}
\arguments{
\item{x}{the input}
}
\value{
TRUE if the input contains all information of
  a BEAST2 posterior. Returns FALSE otherwise.
}
\description{
Determines if the input is a BEAST2 posterior
}
\examples{
trees_filename <- get_tracerer_path("beast2_example_output.trees")
tracelog_filename <- get_tracerer_path("beast2_example_output.log")
posterior <- parse_beast_posterior(
  trees_filename = trees_filename,
  tracelog_filename = tracelog_filename
)
is_posterior(posterior)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cs_std_dev.R
\name{cs_std_dev}
\alias{cs_std_dev}
\title{Calculate the corrected sample standard deviation.}
\usage{
cs_std_dev(values)
}
\arguments{
\item{values}{numeric values}
}
\value{
the corrected sample standard deviation
}
\description{
Calculate the corrected sample standard deviation.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_operators_lines.R
\name{extract_operators_lines}
\alias{extract_operators_lines}
\title{Extract the JSON lines out of a \code{.xml.state} with
  the unparsed BEAST2 MCMC operator acceptances
file with the operators}
\usage{
extract_operators_lines(filename)
}
\arguments{
\item{filename}{name of the BEAST2 .xml.state output file}
}
\value{
the JSON lines of a \code{.xml.state} file with
  the unparsed BEAST2 MCMC operator acceptances
}
\description{
Extract the JSON lines out of a \code{.xml.state} with
  the unparsed BEAST2 MCMC operator acceptances
file with the operators
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_beast_estimates.R
\name{save_beast_estimates}
\alias{save_beast_estimates}
\title{Save the BEAST2 estimates as a BEAST2 \code{.log} file.
There will be some differences: a BEAST2 \code{.log} file also saves
the model as comments and formats the numbers in a way non-standard to R}
\usage{
save_beast_estimates(estimates, filename)
}
\arguments{
\item{estimates}{a data frame of BEAST2 parameter estimates}

\item{filename}{name of the \code{.log} file to save to}
}
\value{
nothing
}
\description{
Save the BEAST2 estimates as a BEAST2 \code{.log} file.
There will be some differences: a BEAST2 \code{.log} file also saves
the model as comments and formats the numbers in a way non-standard to R
}
\seealso{
Use \code{\link{parse_beast_log}} to read a BEAST2 \code{.log} file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_beast_posterior.R
\name{parse_beast_posterior}
\alias{parse_beast_posterior}
\title{Parses BEAST2 output files to a posterior}
\usage{
parse_beast_posterior(
  trees_filenames,
  tracelog_filename,
  log_filename = "deprecated"
)
}
\arguments{
\item{trees_filenames}{the names of one or more a BEAST2
posterior \code{.trees} file.
Each \code{.trees} file can be read using \link{parse_beast_trees}}

\item{tracelog_filename}{name of the BEAST2 tracelog \code{.log} output file,
as can be read using \link{parse_beast_tracelog_file}}

\item{log_filename}{deprecated name
of the BEAST2 tracelog \code{.log} output file.
Use \code{tracelog_filename} instead}
}
\value{
a list with the following elements:\cr
  \itemize{
    item{\code{estimates}: parameter estimates}
    item{
      \code{[alignment_id]_trees}: the phylogenies in the
      BEAST2 posterior. \code{[alignment_id]} is the ID
      of the alignment.
    }
  }
}
\description{
Parses BEAST2 output files to a posterior
}
\examples{
trees_filenames <- get_tracerer_path("beast2_example_output.trees")
tracelog_filename <- get_tracerer_path("beast2_example_output.log")
posterior <- parse_beast_posterior(
  trees_filenames = trees_filenames,
  tracelog_filename = tracelog_filename
)
}
\seealso{
Use \code{\link{remove_burn_ins}} to remove the burn-ins from
  the posterior's estimates (\code{posterior$estimates})
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_stderr_mean.R
\name{calc_stderr_mean}
\alias{calc_stderr_mean}
\title{Calculate the standard error of the mean}
\usage{
calc_stderr_mean(trace)
}
\arguments{
\item{trace}{the values}
}
\value{
the standard error of the mean
}
\description{
Calculate the standard error of the mean
}
\examples{
trace <- sin(seq(from = 0.0, to = 2.0 * pi, length.out = 100))
calc_stderr_mean(trace) # 0.4347425
}
\seealso{
Java code can be found here: \url{https://github.com/beast-dev/beast-mcmc/blob/800817772033c13061f026226e41128d21fd14f3/src/dr/inference/trace/TraceCorrelation.java#L159} # nolint URLs can be long
}
\author{
The original Java version of the algorithm was from Remco Bouckaert,
  ported to R and adapted by Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_mode.R
\name{calc_mode}
\alias{calc_mode}
\title{Calculate the mode of values
If the distribution is bi or multimodal or uniform, NA is returned}
\usage{
calc_mode(values)
}
\arguments{
\item{values}{numeric vector to calculate the mode of}
}
\value{
the mode of the trace
}
\description{
Calculate the mode of values
If the distribution is bi or multimodal or uniform, NA is returned
}
\examples{
# In a unimodal distribution, find the value that occurs most
calc_mode(c(1, 2, 2))
calc_mode(c(1, 1, 2))

# For a uniform distribution, NA is returned
tracerer:::calc_mode(c(1, 2))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_act_cpp}
\alias{calc_act_cpp}
\title{Calculate the auto correlation time
from \url{https://github.com/beast-dev/beast-mcmc/blob/800817772033c13061f026226e41128d21fd14f3/src/dr/inference/trace/TraceCorrelation.java#L159} # nolint}
\usage{
calc_act_cpp(sample, sample_interval)
}
\arguments{
\item{sample}{sample}

\item{sample_interval}{sample interval}
}
\value{
the auto correlation time
}
\description{
Calculate the auto correlation time
from \url{https://github.com/beast-dev/beast-mcmc/blob/800817772033c13061f026226e41128d21fd14f3/src/dr/inference/trace/TraceCorrelation.java#L159} # nolint
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
  log_filename,
  sample_interval,
  state_filename,
  trace,
  tracelog_filename,
  trees_filename,
  trees_filenames,
  verbose
)
}
\arguments{
\item{log_filename}{deprecated name
of the BEAST2 tracelog \code{.log} output file.
Use \code{tracelog_filename} instead}

\item{sample_interval}{the interval in timesteps between samples}

\item{state_filename}{name of the BEAST2 state \code{.xml.state} output file}

\item{trace}{the values}

\item{tracelog_filename}{name of the BEAST2 tracelog \code{.log} output file,
as can be read using \link{parse_beast_tracelog_file}}

\item{trees_filename}{name of a BEAST2 posterior \code{.trees} file,
as can be read using \link{parse_beast_trees}}

\item{trees_filenames}{the names of one or more a BEAST2
posterior \code{.trees} file.
Each \code{.trees} file can be read using \link{parse_beast_trees}}

\item{verbose}{set to TRUE for more output}
}
\description{
Documentation of general function arguments.
This function does nothing.
It is intended to inherit function argument documentation.
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
% Please edit documentation in R/get_path.R
\name{get_tracerer_path}
\alias{get_tracerer_path}
\title{Get the full path of a file in the \code{inst/extdata} folder}
\usage{
get_tracerer_path(filename)
}
\arguments{
\item{filename}{the file's name, without the path}
}
\value{
the full path to the filename
}
\description{
Get the full path of a file in the \code{inst/extdata} folder
}
\examples{
get_tracerer_path("beast2_example_output.log")
get_tracerer_path("beast2_example_output.trees")
get_tracerer_path("beast2_example_output.xml")
get_tracerer_path("beast2_example_output.xml.state")
}
\seealso{
for more files, use \code{\link{get_tracerer_paths}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/save_beast_trees.R
\name{save_beast_trees}
\alias{save_beast_trees}
\title{Save the BEAST2 trees as a BEAST2 \code{.log} file.
There will be some differences: a BEAST2 \code{.log} file also saves
the model as comments and formats the numbers in a way non-standard to R}
\usage{
save_beast_trees(trees, filename)
}
\arguments{
\item{trees}{BEAST2 posterior trees, of type \code{ape::multiPhylo}}

\item{filename}{name of the \code{.trees} file to save to}
}
\value{
nothing
}
\description{
Save the BEAST2 trees as a BEAST2 \code{.log} file.
There will be some differences: a BEAST2 \code{.log} file also saves
the model as comments and formats the numbers in a way non-standard to R
}
\seealso{
Use \code{\link{parse_beast_log}} to read a BEAST2 \code{.log} file
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count_trees_in_file.R
\name{count_trees_in_file}
\alias{count_trees_in_file}
\title{Count the number of trees in a \code{.trees} file}
\usage{
count_trees_in_file(trees_filename)
}
\arguments{
\item{trees_filename}{name of a BEAST2 posterior \code{.trees} file,
as can be read using \link{parse_beast_trees}}
}
\value{
the number of trees
}
\description{
Count the number of trees in a \code{.trees} file
}
\seealso{
if the \code{.trees} file is invalid,
  use \link{is_trees_file} with \code{verbose = TRUE}
  for the reason
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_geom_mean.R
\name{calc_geom_mean}
\alias{calc_geom_mean}
\title{Calculate the geometric mean}
\usage{
calc_geom_mean(values)
}
\arguments{
\item{values}{a numeric vector of values}
}
\value{
returns the geometric mean if all values are at least zero,
  else returns NA
}
\description{
Calculate the geometric mean
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_esses.R
\name{calc_esses}
\alias{calc_esses}
\title{Calculates the Effective Sample Sizes from a parsed BEAST2 log file}
\usage{
calc_esses(traces, sample_interval)
}
\arguments{
\item{traces}{a dataframe with traces with removed burn-in}

\item{sample_interval}{the interval in timesteps between samples}
}
\value{
the effective sample sizes
}
\description{
Calculates the Effective Sample Sizes from a parsed BEAST2 log file
}
\examples{
# Parse an example log file
estimates <- parse_beast_tracelog_file(
  get_tracerer_path("beast2_example_output.log")
)

# Calculate the effective sample sizes of all parameter estimates
calc_esses(estimates, sample_interval = 1000)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_summary_stats.R
\name{calc_summary_stats_traces}
\alias{calc_summary_stats_traces}
\title{Calculates the Effective Sample Sizes of the traces of multiple
  estimated variables.}
\usage{
calc_summary_stats_traces(traces, sample_interval)
}
\arguments{
\item{traces}{a data frame with traces of estimated parameters.
Assumes the burn-ins are removed.}

\item{sample_interval}{the interval in timesteps between samples}
}
\value{
the effective sample sizes
}
\description{
Calculates the Effective Sample Sizes of the traces of multiple
  estimated variables.
}
\examples{
estimates_all <- parse_beast_tracelog_file(
  get_tracerer_path("beast2_example_output.log")
)
estimates <- remove_burn_ins(estimates_all, burn_in_fraction = 0.1)

calc_summary_stats_traces(
  estimates,
  sample_interval = 1000
)
}
\seealso{
Use \code{\link{remove_burn_ins}} to remove the burn-ins
  of all traces
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_paths.R
\name{get_tracerer_paths}
\alias{get_tracerer_paths}
\title{Get the full paths of files in the \code{inst/extdata} folder}
\usage{
get_tracerer_paths(filenames)
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
get_tracerer_paths(
  c(
    "beast2_example_output.log",
    "beast2_example_output.trees",
    "beast2_example_output.xml",
    "beast2_example_output.xml.state"
  )
)
}
\seealso{
for one file, use \code{\link{get_tracerer_path}}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_ess.R
\name{calc_ess}
\alias{calc_ess}
\title{Calculates the Effective Sample Size}
\usage{
calc_ess(trace, sample_interval)
}
\arguments{
\item{trace}{the values without burn-in}

\item{sample_interval}{the interval in timesteps between samples}
}
\value{
the effective sample size
}
\description{
Calculates the Effective Sample Size
}
\examples{
filename <- get_tracerer_path("beast2_example_output.log")
estimates <- parse_beast_tracelog_file(filename)
calc_ess(estimates$posterior, sample_interval = 1000)
}
\seealso{
Java code can be found here: \url{https://github.com/CompEvol/beast2/blob/9f040ed0357c4b946ea276a481a4c654ad4fff36/src/beast/core/util/ESS.java#L161} # nolint URLs can be long
}
\author{
The original Java version of the algorithm was from Remco Bouckaert,
  ported to R and adapted by Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_burn_ins.R
\name{remove_burn_ins}
\alias{remove_burn_ins}
\title{Removed the burn-ins from a data frame}
\usage{
remove_burn_ins(traces, burn_in_fraction = 0.1)
}
\arguments{
\item{traces}{a data frame with traces}

\item{burn_in_fraction}{the fraction that needs to be removed,
must be \code{[0,1>}. Its default value of 10% is the same
as of Tracer}
}
\value{
the data frame with the burn-in removed
}
\description{
Removed the burn-ins from a data frame
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_act.R
\name{calc_act_r}
\alias{calc_act_r}
\title{Calculate the auto-correlation time using only R. Consider using
\link{calc_act} instead, as it is orders of magnitude faster}
\usage{
calc_act_r(trace, sample_interval)
}
\arguments{
\item{trace}{the values}

\item{sample_interval}{the interval in timesteps between samples}
}
\value{
the auto correlation time
}
\description{
Calculate the auto-correlation time using only R. Consider using
\link{calc_act} instead, as it is orders of magnitude faster
}
\examples{
trace <- sin(seq(from = 0.0, to = 2.0 * pi, length.out = 100))
calc_act_r(trace = trace, sample_interval = 1) # 38.18202
}
\seealso{
Java code can be found here: \url{https://github.com/CompEvol/beast2/blob/9f040ed0357c4b946ea276a481a4c654ad4fff36/src/beast/core/util/ESS.java#L161} # nolint URLs can be long
}
\author{
The original Java version of the algorithm was from Remco Bouckaert,
  ported to R and adapted by Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_beast_tracelog_file.R
\name{parse_beast_tracelog_file}
\alias{parse_beast_tracelog_file}
\title{Parses a BEAST2 tracelog \code{.log} output file}
\usage{
parse_beast_tracelog_file(tracelog_filename)
}
\arguments{
\item{tracelog_filename}{name of the BEAST2 tracelog \code{.log} output file,
as can be read using \link{parse_beast_tracelog_file}}
}
\value{
data frame with the parameter estimates
}
\description{
Parses a BEAST2 tracelog \code{.log} output file
}
\examples{
parse_beast_tracelog_file(
  tracelog_filename = get_tracerer_path("beast2_example_output.log")
)
}
\seealso{
Use \code{\link{remove_burn_ins}} to remove the burn-in from
  the returned parameter estimates.
  Use \code{\link{save_beast_estimates}} to save the estimates
  to a \code{.log} file.
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_beast_trees.R
\name{parse_beast_trees}
\alias{parse_beast_trees}
\title{Parses a BEAST2 .trees output file}
\usage{
parse_beast_trees(filename)
}
\arguments{
\item{filename}{name of the BEAST2 .trees output file}
}
\value{
the phylogenies in the posterior
}
\description{
Parses a BEAST2 .trees output file
}
\examples{
trees_filename <- get_tracerer_path("beast2_example_output.trees")
parse_beast_trees(trees_filename)
}
\seealso{
Use \code{\link{save_beast_trees}} to save the phylogenies
  to a \code{.trees} file.
  Use \link{is_trees_file} with \code{verbose = TRUE} to find out
  why a file is invalid
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_trees_posterior.R
\name{is_trees_posterior}
\alias{is_trees_posterior}
\title{Determines if the input is a BEAST2 posterior,
as parsed by parse_beast_trees}
\usage{
is_trees_posterior(x)
}
\arguments{
\item{x}{the input}
}
\value{
TRUE or FALSE
}
\description{
Determines if the input is a BEAST2 posterior,
as parsed by parse_beast_trees
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_hpd_interval.R
\name{calc_hpd_interval}
\alias{calc_hpd_interval}
\title{Calculate the Highest Probability Density of an MCMC trace that
  has its burn-in removed}
\usage{
calc_hpd_interval(trace, proportion = 0.95)
}
\arguments{
\item{trace}{a numeric vector of parameter estimates obtained from
an MCMC run. Must have its burn-in removed}

\item{proportion}{the proportion of numbers within the interval.
For example, use 0.95 for a 95 percentage interval}
}
\value{
a numeric vector, with at index 1 the lower boundary of the
  interval, and at index 2 the upper boundary of the interval
}
\description{
Calculate the Highest Probability Density of an MCMC trace that
  has its burn-in removed
}
\examples{
estimates <- parse_beast_tracelog_file(
  get_tracerer_path("beast2_example_output.log")
)
tree_height_trace <- remove_burn_in(
  estimates$TreeHeight,
  burn_in_fraction = 0.1
)

# Values will be 0.453 and 1.816
calc_hpd_interval(tree_height_trace, proportion = 0.95)
}
\seealso{
The function \code{\link{remove_burn_in}} removes
  a burn-in.
  The Java code that inspired this function can be found here:
  \url{https://github.com/beast-dev/beast-mcmc/blob/98705c59db65e4f406a420bbade949aeecfe05d0/src/dr/stats/DiscreteStatistics.java#L317} # nolint URLs can be long
}
\author{
The original Java version of the algorithm was from J. Heled,
  ported to R and adapted by Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tracerer_tempfilename.R
\name{get_tracerer_tempfilename}
\alias{get_tracerer_tempfilename}
\title{Get a temporary filename}
\usage{
get_tracerer_tempfilename(pattern = "file", fileext = "")
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
named \link{tracerer}.
}
\note{
this function is added to make sure no temporary
cache files are left undeleted
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_summary_stats.R
\name{calc_summary_stats}
\alias{calc_summary_stats}
\title{Calculates the Effective Sample Sizes of one estimated variable's trace.}
\usage{
calc_summary_stats(traces, sample_interval)
}
\arguments{
\item{traces}{one or more traces, supplies as either, (1) a numeric
vector or, (2) a data frame of numeric values.}

\item{sample_interval}{the interval (the number of state
transitions between samples) of the MCMC run that produced the trace.
Using a different \code{sample_interval} than the actually used
sampling interval will result in bogus return values.}
}
\value{
the summary statistics of the traces. If one numeric
vector is supplied, a list is returned with the elements
listed below. If the traces are supplied as a data frame,
a data frame is returned with the elements listed
below as column names.\cr
The elements are:\cr
\itemize{
  \item{\code{mean}: mean}
  \item{\code{stderr_mean}: standard error of the mean}
  \item{\code{stdev}: standard deviation}
  \item{\code{variance}: variance}
  \item{\code{mode}: mode}
  \item{\code{geom_mean}: geometric mean}
  \item{\code{hpd_interval_low}:
    lower bound of 95\% highest posterior density}
  \item{\code{hpd_interval_high}:
    upper bound of 95\% highest posterior density}
  \item{\code{act}: auto correlation time}
  \item{\code{ess}: effective sample size}
}
}
\description{
Calculates the Effective Sample Sizes of one estimated variable's trace.
}
\note{
This function assumes the burn-in is removed.
  Use \code{\link{remove_burn_in}} (on a vector) or
  \code{\link{remove_burn_ins}} (on a data frame) to remove
  the burn-in.
}
\examples{
estimates_all <- parse_beast_tracelog_file(
  get_tracerer_path("beast2_example_output.log")
)
estimates <- remove_burn_ins(estimates_all, burn_in_fraction = 0.1)

# From a single variable's trace
calc_summary_stats(
  estimates$posterior,
  sample_interval = 1000
)

# From all variables' traces
calc_summary_stats(
  estimates,
  sample_interval = 1000
)
}
\seealso{
Use \code{\link{calc_summary_stats_trace}} to calculate the
  summary statistics of one trace (stored as a numeric vector). Use
  \code{\link{calc_summary_stats_traces}} to calculate the
  summary statistics of more traces (stored as a data frame).
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_beast_log.R
\name{parse_beast_log}
\alias{parse_beast_log}
\title{Deprecated function to parse a BEAST2 \code{.log} output file.
Use \link{parse_beast_tracelog_file} instead}
\usage{
parse_beast_log(tracelog_filename, filename = "deprecated")
}
\arguments{
\item{tracelog_filename}{name of the BEAST2 tracelog \code{.log} output file,
as can be read using \link{parse_beast_tracelog_file}}

\item{filename}{deprecated name of the BEAST2 .log output file}
}
\value{
data frame with the parameter estimates
}
\description{
Deprecated function to parse a BEAST2 \code{.log} output file.
Use \link{parse_beast_tracelog_file} instead
}
\examples{
# Deprecated
parse_beast_log(
  tracelog_filename = get_tracerer_path("beast2_example_output.log")
)
# Use the function 'parse_beast_tracelog_file' instead
parse_beast_tracelog_file(
  tracelog_filename = get_tracerer_path("beast2_example_output.log")
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{calc_std_error_of_mean_cpp}
\alias{calc_std_error_of_mean_cpp}
\title{Calculates the standard error of the mean}
\usage{
calc_std_error_of_mean_cpp(sample)
}
\arguments{
\item{sample}{numeric vector of values}
}
\value{
the standard error of the mean
}
\description{
Calculates the standard error of the mean
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_summary_stats.R
\name{calc_summary_stats_trace}
\alias{calc_summary_stats_trace}
\title{Calculates the Effective Sample Sizes of one estimated variable's trace.}
\usage{
calc_summary_stats_trace(trace, sample_interval)
}
\arguments{
\item{trace}{a numeric vector of values. Assumes the burn-in
is removed.}

\item{sample_interval}{the interval in timesteps between samples}
}
\value{
the effective sample sizes
}
\description{
Calculates the Effective Sample Sizes of one estimated variable's trace.
}
\examples{
estimates_all <- parse_beast_tracelog_file(
  get_tracerer_path("beast2_example_output.log")
)
estimates <- remove_burn_ins(estimates_all, burn_in_fraction = 0.1)

calc_summary_stats_trace(
  estimates$posterior,
  sample_interval = 1000
)
}
\seealso{
Use \code{\link{remove_burn_in}} to remove the burn-in
  of a trace
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_beast_state.R
\name{parse_beast_state_operators}
\alias{parse_beast_state_operators}
\title{Parses a BEAST2 state \code{.xml.state} output file to get only the operators
  acceptances}
\usage{
parse_beast_state_operators(
  state_filename = get_tracerer_path("beast2_example_output.xml.state"),
  filename = "deprecated"
)
}
\arguments{
\item{state_filename}{name of the BEAST2 state \code{.xml.state} output file}

\item{filename}{deprecated name of the BEAST2 .xml.state output file,
use \code{state_filename} instead}
}
\value{
data frame with all the operators' success rates
}
\description{
Parses a BEAST2 state \code{.xml.state} output file to get only the operators
  acceptances
}
\examples{
parse_beast_state_operators(
  state_filename = get_tracerer_path("beast2_example_output.xml.state")
)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/check_trace.R
\name{check_trace}
\alias{check_trace}
\title{Check if the trace is a valid. Will \link{stop} if not}
\usage{
check_trace(trace)
}
\arguments{
\item{trace}{the values}
}
\description{
Check if the trace is a valid. Will \link{stop} if not
}
\examples{
check_trace(seq(1, 2))
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/remove_burn_in.R
\name{remove_burn_in}
\alias{remove_burn_in}
\title{Removed the burn-in from a trace}
\usage{
remove_burn_in(trace, burn_in_fraction)
}
\arguments{
\item{trace}{the values}

\item{burn_in_fraction}{the fraction that needs to be removed, must be [0,1>}
}
\value{
the values with the burn-in removed
}
\description{
Removed the burn-in from a trace
}
\examples{
# Create a trace from one to and including ten
v <- seq(1, 10)

# Remove the first ten percent of its values,
# in this case removes the first value, which is one
w <- remove_burn_in(trace = v, burn_in_fraction = 0.1)
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calc_act.R
\name{calc_act}
\alias{calc_act}
\title{Calculate the auto-correlation time, alternative implementation}
\usage{
calc_act(trace, sample_interval)
}
\arguments{
\item{trace}{the values}

\item{sample_interval}{the interval in timesteps between samples}
}
\value{
the auto_correlation time
}
\description{
Calculate the auto-correlation time, alternative implementation
}
\examples{
trace <- sin(seq(from = 0.0, to = 2.0 * pi, length.out = 100))
# 38.18202
calc_act(trace = trace, sample_interval = 1)
}
\seealso{
Java code can be found here: \url{https://github.com/CompEvol/beast2/blob/9f040ed0357c4b946ea276a481a4c654ad4fff36/src/beast/core/util/ESS.java#L161} # nolint URLs can be long
}
\author{
The original Java version of the algorithm was from Remco Bouckaert,
  ported to R and adapted by Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_trees_file.R
\name{is_trees_file}
\alias{is_trees_file}
\title{Measure if a file a valid BEAST2 \code{.trees} file}
\usage{
is_trees_file(trees_filename, verbose = FALSE)
}
\arguments{
\item{trees_filename}{name of a BEAST2 posterior \code{.trees} file,
as can be read using \link{parse_beast_trees}}

\item{verbose}{set to TRUE for more output}
}
\value{
TRUE if \code{trees_filename} is a valid \code{.trees} file
}
\description{
Measure if a file a valid BEAST2 \code{.trees} file
}
\examples{
# TRUE
is_trees_file(get_tracerer_path("beast2_example_output.trees"))
is_trees_file(get_tracerer_path("unplottable_anthus_aco.trees"))
is_trees_file(get_tracerer_path("anthus_2_4_a.trees"))
is_trees_file(get_tracerer_path("anthus_2_4_b.trees"))
# FALSE
is_trees_file(get_tracerer_path("mcbette_issue_8.trees"))
}
\seealso{
Most of the work is done by \link[ape]{read.nexus}
}
\author{
Richèl J.C. Bilderbeek
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tracerer.R
\docType{package}
\name{tracerer}
\alias{tracerer}
\title{\code{tracerer}: A package to parse BEAST2 output files.}
\description{
\code{tracerer} allows to parse BEAST2 input files, using
an R interface. 'tracerer' closely follows the functionality
of Tracer, a GUI tool bundled with BEAST and BEAST2,
including its default settings.
}
\seealso{
These are packages associated with \code{tracerer}:
\itemize{
  \item{
    The package \code{beautier} can create
    BEAST2 input files from R
  }
  \item{
    The package \code{beastier} can run
    BEAST2 from R
  }
  \item{
    The package \code{mauricer} manages
    BEAST2 packages from R
  }
  \item{
    The package \code{babette} combines the
    functionality of \code{beautier},
    \code{beastier}, \code{tracerer} and
    \code{mauricer} and
    into a single workflow
  }
}
If something is (still) missing from \code{tracerer},
the \code{coda} package may have the functionality you need.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_beast_output_files.R
\name{parse_beast_output_files}
\alias{parse_beast_output_files}
\title{Parse all BEAST2 output files}
\usage{
parse_beast_output_files(log_filename, trees_filenames, state_filename)
}
\arguments{
\item{log_filename}{deprecated name
of the BEAST2 tracelog \code{.log} output file.
Use \code{tracelog_filename} instead}

\item{trees_filenames}{the names of one or more a BEAST2
posterior \code{.trees} file.
Each \code{.trees} file can be read using \link{parse_beast_trees}}

\item{state_filename}{name of the BEAST2 state \code{.xml.state} output file}
}
\value{
a list with the following elements:\cr
  \itemize{
    item{\code{estimates}: parameter estimates}
    item{
      \code{[alignment_id]_trees}: the phylogenies in the
      BEAST2 posterior. \code{[alignment_id]} is the ID
      of the alignment.
    }
    item{\code{operators}: the BEAST2 MCMC operator
      acceptances
    }
  }
}
\description{
Parse all BEAST2 output files
}
\examples{
trees_filenames <- get_tracerer_path("beast2_example_output.trees")
log_filename <- get_tracerer_path("beast2_example_output.log")
state_filename <- get_tracerer_path("beast2_example_output.xml.state")
parse_beast_output_files(
  log_filename = log_filename,
  trees_filenames = trees_filenames,
  state_filename = state_filename
)
}
\seealso{
Use \code{\link{remove_burn_ins}} to remove the burn-in from
  \code{out$estimates}
}
\author{
Richèl J.C. Bilderbeek
}
