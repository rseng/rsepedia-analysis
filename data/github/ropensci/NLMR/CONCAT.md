# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(http:contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/

[![Build Status](https://travis-ci.org/ropensci/NLMR.svg?branch=master)](https://travis-ci.org/ropensci/NLMR)[![Build status](https://ci.appveyor.com/api/projects/status/djw840fitcvolbxg?svg=true)](https://ci.appveyor.com/project/ropensci/NLMR) [![codecov](https://codecov.io/gh/ropensci/NLMR/branch/develop/graph/badge.svg?token=MKCm2fVrDa)](https://codecov.io/gh/ropensci/NLMR) [![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/NLMR)](https://cran.r-project.org/package=NLMR) [![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing) [![](http://cranlogs.r-pkg.org/badges/grand-total/NLMR)](http://cran.rstudio.com/web/packages/NLMR/index.html) [![](https://badges.ropensci.org/188_status.svg)](https://github.com/ropensci/onboarding/issues/188) [![DOI:10.1111/2041-210X.13076](https://zenodo.org/badge/DOI/10.1111/2041-210X.13076.svg)](https://doi.org/10.1111/2041-210X.13076)

NLMR <img src="man/figures/logo.png" align="right" width="150" />
=================================================================

**NLMR** is an `R` package for simulating **n**eutral **l**andscape **m**odels (NLM). Designed to be a generic framework like [NLMpy](https://pypi.python.org/pypi/nlmpy), it leverages the ability to simulate the most common NLM that are described in the ecological literature. **NLMR** builds on the advantages of the **raster** package and returns all simulation as `RasterLayer` objects, thus ensuring a direct compatibility to common GIS tasks and a flexible and simple usage. Furthermore, it simulates NLMs within a self-contained, reproducible framework.

Installation
------------

Install the release version from CRAN:

``` r
install.packages("NLMR")
```

To install the developmental version of **NLMR**, use the following R code:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/NLMR")
```

Example
-------

Each neutral landscape models is simulated with a single function (all starting with `nlm_`) in `NLMR`, e.g.:

``` r
random_cluster <- NLMR::nlm_randomcluster(nrow = 100,
                                      ncol = 100,
                                      p    = 0.5,
                                      ai   = c(0.3, 0.6, 0.1),
                                      rescale = FALSE)

random_curdling <- NLMR::nlm_curds(curds = c(0.5, 0.3, 0.6),
                              recursion_steps = c(32, 6, 2))


midpoint_displacememt <- NLMR::nlm_mpd(ncol = 100,
                                 nrow = 100,
                                 roughness = 0.61)
```

Overview
--------

**NLMR** supplies 15 NLM algorithms, with several options to simulate derivatives of them. The algorithms differ from each other in spatial auto-correlation, from no auto-correlation (random NLM) to a constant gradient (planar gradients):

<table class="table table-striped table-hover table-condensed" style="margin-left: auto; margin-right: auto;">
<thead>
<tr>
<th style="text-align:left;">
Function
</th>
<th style="text-align:left;">
Description
</th>
<th style="text-align:left;">
Crossreference
</th>
<th style="text-align:left;">
Reference
</th>
</tr>
</thead>
<tbody>
<tr>
<td style="text-align:left;">
nlm\_curds
</td>
<td style="text-align:left;">
Simulates a randomly curdled or wheyed neutral landscape model. Random curdling recursively subdivides the landscape into blocks. At each level of the recursion, a fraction of these blocks is declared as habitat while the remaining stays matrix. When option q is set, it simulates a wheyed curdling model, where previously selected cells that were declared matrix during recursion, can now contain a proportion of habitat cells
</td>
<td style="text-align:left;">
Figure 1a,p
</td>
<td style="text-align:left;">
O’Neill, Gardner, and Turner (1992); Keitt (2000)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_distancegradient
</td>
<td style="text-align:left;">
Simulates a distance gradient neutral landscape model. The gradient is always measured from a rectangle that one has to specify in the function (parameter origin)
</td>
<td style="text-align:left;">
Figure 1b
</td>
<td style="text-align:left;">
Etherington, Holland, and O’Sullivan (2015)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_edgegradient
</td>
<td style="text-align:left;">
Simulates a linear gradient orientated neutral model. The gradient has a specified or random direction that has a central peak, which runs perpendicular to the gradient direction
</td>
<td style="text-align:left;">
Figure 1c
</td>
<td style="text-align:left;">
Travis and Dytham (2004); Schlather et al. (2015)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_fbm
</td>
<td style="text-align:left;">
Simulates neutral landscapes using fractional Brownian motion (fBm). fBm is an extension of Brownian motion in which the amount of spatial autocorrelation between steps is controlled by the Hurst coefficient H
</td>
<td style="text-align:left;">
Figure 1d
</td>
<td style="text-align:left;">
Schlather et al. (2015)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_gaussianfield
</td>
<td style="text-align:left;">
Simulates a spatially correlated random fields (Gaussian random fields) model, where one can control the distance and magnitude of spatial autocorrelation
</td>
<td style="text-align:left;">
Figure 1e
</td>
<td style="text-align:left;">
Schlather et al. (2015)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_mosaicfield
</td>
<td style="text-align:left;">
Simulates a mosaic random field neutral landscape model. The algorithm imitates fault lines by repeatedly bisecting the landscape and lowering the values of cells in one half and increasing the values in the other half. If one sets the parameter infinite to TRUE, the algorithm approaches a fractal pattern
</td>
<td style="text-align:left;">
Figure 1f
</td>
<td style="text-align:left;">
Schlather et al. (2015)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_neigh
</td>
<td style="text-align:left;">
Simulates a neutral landscape model with land cover classes and clustering based on neighbourhood characteristics. The cluster are based on the surrounding cells. If there is a neighbouring cell of the current value/type, the target cell will more likely turned into a cell of that type/value
</td>
<td style="text-align:left;">
Figure 1g
</td>
<td style="text-align:left;">
Scherer et al. (2016)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_percolation
</td>
<td style="text-align:left;">
Simulates a binary neutral landscape model based on percolation theory. The probability for a cell to be assigned habitat is drawn from a uniform distribution
</td>
<td style="text-align:left;">
Figure 1h
</td>
<td style="text-align:left;">
Gardner et al. (1989)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_planargradient
</td>
<td style="text-align:left;">
Simulates a planar gradient neutral landscape model. The gradient is sloping in a specified or (by default) random direction between 0 and 360 degree
</td>
<td style="text-align:left;">
Figure 1i
</td>
<td style="text-align:left;">
Palmer (1992)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_mosaictess
</td>
<td style="text-align:left;">
Simulates a patchy mosaic neutral landscape model based on the tessellation of a random point process. The algorithm randomly places points (parameter germs) in the landscape, which are used as the centroid points for a voronoi tessellation. A higher number of points therefore leads to a more fragmented landscape
</td>
<td style="text-align:left;">
Figure 1k
</td>
<td style="text-align:left;">
Gaucherel (2008), Method 1
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_mosaicgibbs
</td>
<td style="text-align:left;">
Simulates a patchy mosaic neutral landscape model based on the tessellation of an inhibition point process. This inhibition point process starts with a given number of points and uses a minimisation approach to fit a point pattern with a given interaction parameter (0 - hardcore process; 1 - Poisson process) and interaction radius (distance of points/germs being apart)
</td>
<td style="text-align:left;">
Figure 1l
</td>
<td style="text-align:left;">
Gaucherel (2008), Method 2
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_random
</td>
<td style="text-align:left;">
Simulates a spatially random neutral landscape model with values drawn a uniform distribution
</td>
<td style="text-align:left;">
Figure 1m
</td>
<td style="text-align:left;">
With and Crist (1995)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_randomcluster
</td>
<td style="text-align:left;">
Simulates a random cluster nearest-neighbour neutral landscape. The parameter ai controls for the number and abundance of land cover classes and p controls for proportion of elements randomly selected to form clusters
</td>
<td style="text-align:left;">
Figure 1n
</td>
<td style="text-align:left;">
Saura and Martínez-Millán (2000)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_mpd
</td>
<td style="text-align:left;">
Simulates a midpoint displacement neutral landscape model where the parameter roughness controls the level of spatial autocorrelation
</td>
<td style="text-align:left;">
Figure 1n
</td>
<td style="text-align:left;">
Peitgen and Saupe (1988)
</td>
</tr>
<tr>
<td style="text-align:left;">
nlm\_randomrectangularcluster
</td>
<td style="text-align:left;">
Simulates a random rectangular cluster neutral landscape model. The algorithm randomly distributes overlapping rectangles until the landscape is filled
</td>
<td style="text-align:left;">
Figure 1o
</td>
<td style="text-align:left;">
Gustafson and Parker (1992)
</td>
</tr>
</tbody>
</table>
<img src="https://wol-prod-cdn.literatumonline.com/cms/attachment/b963a726-ed88-4ede-863c-a65451f91d0f/mee313076-fig-0001-m.jpg"  width="100%" />

See also
--------

**NLMR** was split during its development process - to have a minimal dependency version for simulating neutral landscape models and an utility toolbox to facilitate workflows with raster data. If you are interested in merging, visualizing or further handling neutral landscape models have a look at the [landscapetools](https://github.com/ropensci/landscapetools/) package.

Meta
----

-   Please [report any issues or bugs](https://github.com/ropensci/NLMR/issues/new/).
-   License: GPL3
-   Get citation information for `NLMR` in R doing `citation(package = 'NLMR')`
    -   Additionally, we keep a [record of publications](https://ropensci.github.io/NLMR/articles/publication_record.html/) that use **NLMR**. Hence, if you used **NLMR** please [file an issue on GitHub](https://github.com/ropensci/NLMR/issues/new/) so we can add it to the list.
-   We are very open to contributions - if you are interested check out our [Contributor Guidelines](CONTRIBUTING.md).
    -   Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
____________________________________________________________________________________
## NLMR 0.4.2 Release Notes

- Bugfix in nlm_mosaicfield to rely on new version of RandomFields

## NLMR 0.4.1 Release Notes

- Bugfix in nlm_mpd to not rely on landscapetools

## NLMR 0.4 Release Notes

- nlm_neigh, nlm_mpd and nlm_randomrectangularcluster are now implemented in Rcpp
- all of the Rcpp also take the R random seed
- Minor bug fixes
- Improvements to documentation
- More examples on the package website

## NLMR 0.3.2 Release Notes

- Update citation 

## NLMR 0.3.1 Release Notes

- Minor bug fixes
- Updated documentation
- removed purrr as dependency

## NLMR 0.3.0 Release Notes

- successful review through rOpenSci
- split package into two packages:
  - `NLMR` 
    - contains now only the neutral landscape models, minimal dependencies
  - [`landscapetools`](https://github.com/ropensci/landscapetools)
    - contains now only utility functions
- small bug fixes
- `nlm_fBm` is now `nlm_fbm`

## NLMR 0.2.1 Release Notes

- Skip one test on CRAN to keep the Roboto font available
- Function `show_landscape` to plot a list of rasters as ggplot2 facet
- Small updates to the webpage

## NLMR 0.2 Release Notes

- Small bug fixes
- New neutral landscape models
    - `nlm_wheys`: Simulates a wheyed neutral landscape model
- Parameter `p` in `nlm_curds` now controls the proportion of habitat instead of 
  the amount of matrix
- Implemented new theme `theme_nlm`
- Functions to coerce raster to tibbles and vice versa (for facetting with `ggplot2`)
- We now have unit tests covering the main functionality of the package
- Removed several packages as dependencies 

## NLMR 0.1.0 Release Notes

v0.1.0 was released on 30/11/2017

- First stable release of NLMR
# CONTRIBUTING #

### Please contribute!

We love collaboration.

### Bugs?

* Submit an issue on the Issues page [here](https://github.com/ropensci/nlmr/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/nlmr.git`
* Make sure to track progress upstream (i.e., on our version of `scrubr` at `ropensci/scrubr`) by doing `git remote add upstream https://github.com/ropensci/nlmr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new branch)
* If you alter package functionality at all (e.g., the code itself, not just documentation)
please do write some tests to cover the new functionality.
* Push up to your account
* Submit a pull request to home base at `ropensci/nlmr`

### Questions? Get in touch: [sciaini.marco@gmail.com](mailto:sciaini.marco@gmail.com)

### Thanks for contributing!
## Version update

Update spatstat dep 

## Test environments

* local Ubuntu Linux 18.10 LTS install
* Ubuntu 14.04 (on travis-ci)
* Windows Server 2012 R2 x64 (build 9600) (on appveyor)
* rhub
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 note

## Reverse dependencies

There are currently no reverse dependencies.
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/{owner}/{repo}/issues)

### Issues and Pull Requests

If you are considering a pull request, you may want to open an issue first to discuss with the maintainer(s).

### Code contributions

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/{repo}.git`
* Make sure to track progress upstream (i.e., on our version of `{repo}` at `{owner}/{repo}`) by doing `git remote add upstream https://github.com/{owner}/{repo}.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch - see <https://guides.github.com/introduction/flow/> for how to contribute by branching, making changes, then submitting a pull request)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `{owner}/{repo}`

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.

### Prefer to Email?

Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Thanks for contributing!
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
output:
  github_document:
    html_preview: false
editor_options: 
  chunk_output_type: console
always_allow_html: yes
---

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "vignettes/README-"
)
```

 <!-- badges: start -->
 
 [![R-CMD-check](https://github.com/ropensci/NLMR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/NLMR/actions)
[![codecov](https://codecov.io/gh/ropensci/NLMR/branch/develop/graph/badge.svg?token=MKCm2fVrDa)](https://codecov.io/gh/ropensci/NLMR)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/NLMR)](https://cran.r-project.org/package=NLMR) 
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](http://cranlogs.r-pkg.org/badges/grand-total/NLMR)](http://cran.rstudio.com/web/packages/NLMR/index.html) 
[![](https://badges.ropensci.org/188_status.svg)](https://github.com/ropensci/onboarding/issues/188)
[![DOI:10.1111/2041-210X.13076](https://zenodo.org/badge/DOI/10.1111/2041-210X.13076.svg)](https://doi.org/10.1111/2041-210X.13076)

<!-- badges: end -->

# NLMR <img src="man/figures/logo.png" align="right" width="150" />

**NLMR** is an ``R`` package for simulating **n**eutral **l**andscape **m**odels (NLM). Designed to be a generic framework like [NLMpy](https://pypi.python.org/pypi/nlmpy), it leverages the ability to simulate the most common NLM that are described in the ecological literature. 
**NLMR** builds on the advantages of the **raster** package and returns all simulation as ``RasterLayer`` objects, thus ensuring a direct compatibility to common GIS tasks and a flexible and simple usage.
Furthermore, it simulates NLMs within a self-contained, reproducible framework.

## Installation

Install the release version from CRAN:

```{r eval = FALSE}
install.packages("NLMR")
```

To install the developmental version of **NLMR**, use the following R code:

```{r eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/NLMR")
```

## Example

Each neutral landscape models is simulated with a single function (all starting with `nlm_`) in `NLMR`, e.g.:

```{r eval=FALSE}
random_cluster <- NLMR::nlm_randomcluster(nrow = 100,
                                      ncol = 100,
                                      p    = 0.5,
                                      ai   = c(0.3, 0.6, 0.1),
                                      rescale = FALSE)

random_curdling <- NLMR::nlm_curds(curds = c(0.5, 0.3, 0.6),
                              recursion_steps = c(32, 6, 2))


midpoint_displacememt <- NLMR::nlm_mpd(ncol = 100,
                                 nrow = 100,
                                 roughness = 0.61)
```

## Overview

**NLMR** supplies 15 NLM algorithms, with several options to simulate derivatives of
them. The algorithms differ from each other in spatial auto-correlation, from 
no auto-correlation (random NLM) to a constant gradient (planar gradients):  

```{r warning=FALSE, message = FALSE, results='asis', echo=FALSE, cache=FALSE}
library(tibble)
library(magrittr)
library(knitr)
library(kableExtra)

function_tibble <- tibble(Function = character(), Description = character(), 	Crossreference = character(), Reference = character())

# nlm_curds
function_tibble[1,1] <- "nlm_curds"
function_tibble[1,2] <- "Simulates a randomly curdled or wheyed neutral landscape model. Random curdling recursively subdivides the landscape into blocks. At each level of the recursion, a fraction of these blocks is declared as habitat while the remaining stays matrix. When option q is set, it simulates a wheyed curdling model, where previously selected cells that were declared matrix during recursion, can now contain a proportion of habitat cells"
function_tibble[1,3] <- "Figure 1a,p"
function_tibble[1,4] <- "O’Neill, Gardner, and Turner (1992); Keitt (2000)"

# nlm_distancegradient
function_tibble[2,1] <- "nlm_distancegradient"
function_tibble[2,2] <- "Simulates a distance gradient neutral landscape model. The gradient is always measured from a rectangle that one has to specify in the function (parameter origin)"
function_tibble[2,3] <- "Figure 1b"
function_tibble[2,4] <- "Etherington, Holland, and O’Sullivan (2015)"

# nlm_edgegradient
function_tibble[3,1] <- "nlm_edgegradient"
function_tibble[3,2] <- "Simulates a linear gradient orientated neutral model. The gradient has a specified or random direction that has a central peak, which runs perpendicular to the gradient direction"
function_tibble[3,3] <- "Figure 1c"
function_tibble[3,4] <- "Travis and Dytham (2004); Schlather et al. (2015)"

# nlm_edgegradient
function_tibble[4,1] <- "nlm_fbm"
function_tibble[4,2] <- "Simulates neutral landscapes using fractional Brownian motion (fBm). fBm is an extension of Brownian motion in which the amount of spatial autocorrelation between steps is controlled by the Hurst coefficient H"
function_tibble[4,3] <- "Figure 1d"
function_tibble[4,4] <- "Schlather et al. (2015)"

# nlm_gaussianfield
function_tibble[5,1] <- "nlm_gaussianfield"
function_tibble[5,2] <- "Simulates a spatially correlated random fields (Gaussian random fields) model, where one can control the distance and magnitude of spatial autocorrelation	"
function_tibble[5,3] <- "Figure 1e"
function_tibble[5,4] <- "Schlather et al. (2015)"

# nlm_mosaicfield
function_tibble[6,1] <- "nlm_mosaicfield"
function_tibble[6,2] <- "Simulates a mosaic random field neutral landscape model. The algorithm imitates fault lines by repeatedly bisecting the landscape and lowering the values of cells in one half and increasing the values in the other half. If one sets the parameter infinite to TRUE, the algorithm approaches a fractal pattern"
function_tibble[6,3] <- "Figure 1f"
function_tibble[6,4] <- "Schlather et al. (2015)"

# nlm_neigh
function_tibble[7,1] <- "nlm_neigh"
function_tibble[7,2] <- "Simulates a neutral landscape model with land cover classes and clustering based on neighbourhood characteristics. The cluster are based on the surrounding cells. If there is a neighbouring cell of the current value/type, the target cell will more likely turned into a cell of that type/value"
function_tibble[7,3] <- "Figure 1g"
function_tibble[7,4] <- "Scherer et al. (2016)"

# nlm_percolation
function_tibble[8,1] <- "nlm_percolation"
function_tibble[8,2] <- "Simulates a binary neutral landscape model based on percolation theory. The probability for a cell to be assigned habitat is drawn from a uniform distribution"
function_tibble[8,3] <- "Figure 1h"
function_tibble[8,4] <- "Gardner et al. (1989)"

# nlm_planargradient
function_tibble[9,1] <- "nlm_planargradient"
function_tibble[9,2] <- "Simulates a planar gradient neutral landscape model. The gradient is sloping in a specified or (by default) random direction between 0 and 360 degree"
function_tibble[9,3] <- "Figure 1i"
function_tibble[9,4] <- "Palmer (1992)"

# nlm_mosaictess
function_tibble[10,1] <- "nlm_mosaictess"
function_tibble[10,2] <- "Simulates a patchy mosaic neutral landscape model based on the tessellation of a random point process. The algorithm randomly places points (parameter germs) in the landscape, which are used as the centroid points for a voronoi tessellation. A higher number of points therefore leads to a more fragmented landscape"
function_tibble[10,3] <- "Figure 1k"
function_tibble[10,4] <- "Gaucherel (2008), Method 1"

# nlm_mosaicgibbs	
function_tibble[11,1] <- "nlm_mosaicgibbs	"
function_tibble[11,2] <- "Simulates a patchy mosaic neutral landscape model based on the tessellation of an inhibition point process. This inhibition point process starts with a given number of points and uses a minimisation approach to fit a point pattern with a given interaction parameter (0 ‐ hardcore process; 1 ‐ Poisson process) and interaction radius (distance of points/germs being apart)"
function_tibble[11,3] <- "Figure 1l"
function_tibble[11,4] <- "Gaucherel (2008), Method 2"

# nlm_random
function_tibble[12,1] <- "nlm_random"
function_tibble[12,2] <- "Simulates a spatially random neutral landscape model with values drawn a uniform distribution"
function_tibble[12,3] <- "Figure 1m"
function_tibble[12,4] <- "With and Crist (1995)"
# nlm_randomcluster
function_tibble[13,1] <- "nlm_randomcluster"
function_tibble[13,2] <- "Simulates a random cluster nearest‐neighbour neutral landscape. The parameter ai controls for the number and abundance of land cover classes and p controls for proportion of elements randomly selected to form clusters"
function_tibble[13,3] <- "Figure 1n"
function_tibble[13,4] <- "Saura and Martínez-Millán (2000)"

# nlm_mpd
function_tibble[14,1] <- "nlm_mpd"
function_tibble[14,2] <- "Simulates a midpoint displacement neutral landscape model where the parameter roughness controls the level of spatial autocorrelation"
function_tibble[14,3] <- "Figure 1n"
function_tibble[14,4] <- "Peitgen and Saupe (1988)"

# nlm_randomrectangularcluster
function_tibble[15,1] <- "nlm_randomrectangularcluster"
function_tibble[15,2] <- "Simulates a random rectangular cluster neutral landscape model. The algorithm randomly distributes overlapping rectangles until the landscape is filled"
function_tibble[15,3] <- "Figure 1o"
function_tibble[15,4] <- "Gustafson and Parker (1992)"

kable(function_tibble) %>%
  kable_styling(bootstrap_options = c("striped", "hover", "condensed"))
```

<img src="https://wol-prod-cdn.literatumonline.com/cms/attachment/b963a726-ed88-4ede-863c-a65451f91d0f/mee313076-fig-0001-m.jpg"  width="100%" />

## See also

**NLMR** was split during its development process - to have a minimal dependency version
for simulating neutral landscape models and an utility toolbox to facilitate workflows
with raster data.
If you are interested in merging, visualizing or further handling neutral landscape models
have a look at the [landscapetools](https://github.com/ropensci/landscapetools/) package.

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/NLMR/issues/new/).
* License: GPL3
* Get citation information for `NLMR` in R doing `citation(package = 'NLMR')`
    * Additionally, we keep a [record of publications](https://ropensci.github.io/NLMR/articles/publication_record.html/) that use **NLMR**. Hence, if you used **NLMR** please [file an issue on GitHub](https://github.com/ropensci/NLMR/issues/new/) so we can add it to the list.
* We are very open to contributions - if you are interested check out our [Contributor Guidelines](CONTRIBUTING.md).
    * Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
title: "Basic Usage of NLMR"
author: "Marco Sciaini & Craig E.Simpkins"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
bibliography: Citations.bib
vignette: >
  %\VignetteIndexEntry{Basic Usage of NLMR}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r global_options, include=FALSE}
library(raster)
library(NLMR)
```

`NLMR` is a `R` package designed to generate neutral landscape models (NLMs), 
simulated landscapes used to explore landscape scale ecological patterns and 
processes. The `NLMR` package was designed with a similar philosophy to the 
Python package `NLMpy` [see @Etheringtonnlmr2015], offering a general numeric 
framework allowing for a high degree of flexibility. Most of the common NLMs, 
as described by the relevant literature, can be produced using NLMR. 
Additionally, NLMR allows users to merge multiple landscapes, classify landscape
elements categorically and measure basic landscape level metrics. All NLMs 
produced take the form of two-dimensional raster arrays with specified row and 
column dimensions and cell values ranging between 0 and 1. By returning raster 
arrays, NLMs are easily integrated into the workflow of many useful spatial
analysis packages, notably the `raster` package.

For further information on neutral landscape models, the authors goals for
this package, and additional use case examples please see the associated
publication Sciani, Fritsch, Scherer and Simpkins [-@Sciani2018]

## Basic landscape generation

`NLMR` supplies 16 NLM algorithms. The algorithms differ from each other in 
spatial auto-correlation, from no auto-correlation (random NLM) to a constant 
gradient (planar gradients) [see @Palmer1992].  

The 16 NLM algorithms are:

1. distance gradient
1. edge gradient  
1. hierarchical curdling
1. wheyed hierarchical curdling
1. midpoint displacement 
1. neighbourhood clustering
1. planar gradient  
1. random  
1. random cluster nearest-neighbour 
1. random element 
1. random mosaic fields
1. random polygonal landscapes
1. random percolation  
1. random rectangular cluster 
1. spatially correlated random fields (Gaussian random fields)  
1. two-dimensional fractional Brownian motion

The basic syntax used to produce a NLM landscape is:
```
nlm_modeltype(ncol, nrow, resolution, ...)
```
    
For example, to produce a simple random neutral landscape one could use the 
following code:

```{r, fig.height=7, fig.width=7, fig.align='center'}
x <- NLMR::nlm_random(20,20)
plot(x)
```

## Merging landscapes

Multiple NLM rasters can be merged or merged together to create new landscape 
patterns. A single primary or base raster can be merged with any number of 
additional secondary rasters, with optional scaling factors used to control the 
influence of the secondary rasters.  

The `util_merge` function is used to merge the rasters as in the example below:

```{r, fig.height=7, fig.width=7, fig.align='center'}
  #Create primary landscape raster
  pL <- NLMR::nlm_edgegradient(ncol = 100,
                               nrow = 100)

  plot(pL)

  #Create secondary landscape rasters
  sL1 <- NLMR::nlm_distancegradient(ncol = 100,
                                    nrow = 100,
                                    origin = c(10, 10, 10, 10))
  sL2 <- NLMR::nlm_random(ncol = 100,
                          nrow = 100)

  mL1 <- pL + (sL1 + sL2)
  
  plot(mL1)
```

## Classifying categories

Landscape rasters generated by `NLMR` contain continuous values between 0 and 1,
though these can be converted into categorical values using `util_classify` from 
**landscapetools**. By default classes are numerical starting from 1. 
If non-numerical levels are required, `level_names` can be 
specified. These classes can be plotted by selecting `discrete = TRUE` in 
`show_landscape`.

```{r fig.height=7, fig.width=7, fig.align='center'}
nr <- NLMR::nlm_fbm(50, 100, fract_dim = 1.2)
                              
nr_classified <- landscapetools::util_classify(nr, weighting = c(0.3, 0.3, 0.3))

plot(nr_classified)
```

## References
---
title: "Visualize Neutral Landscape Models"
author: "Marco Sciaini"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Visualize Neutral Landscape Models}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

As the ever growing R package environment can be a rough terrain to navigate and find the appropriate tools to achieve one's goals, this vignette is meant to point out some ways to overcome initial problems with visualizing neutral landscape models or more general raster data. This is probably a heavily biased view on packages and functions and I am sure there are other good R packages out there to achieve the same (if so - feel free to point that out to me and I will include it!). However, I am also sure this collection can at least be a kickstart for quickly visualizing your results and help you to communicate them.

## Static plots

### landscapetools

**landscapetools** function `show_landscape` was developed to help users to adhere to
some standards concerning color scales and typography. This means for example 
that by default the [viridis color scale](https://cran.r-project.org/web/packages/viridis/vignettes/intro-to-viridis.html) 
is applied which makes your plots easier to read by those with colorblindness.

```{r , fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library("NLMR")
library("landscapetools")

landscape <- nlm_mosaictess(200, 200, germs = 444)

# default theme
show_landscape(landscape)

# ... chose another color scale from viridis ("E" = cividis)
show_landscape(landscape, viridis_scale = "E")

# ... chose any other scale:
# show_landscape returns a ggplot2 object, so you can follow your usual ggplot2
# workflow and change the color, axis labels, ...
library(ggplot2)
library(pals)
show_landscape(landscape)  + 
  scale_fill_gradientn(colours=pals::parula(100)) + # parula color scale
  theme_void() +  # minimal theme
  guides(fill = FALSE) # remove legend
```

### rasterVis

**rasterVis** also offers some convenience functions to plot raster, for example:

```{r , fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library("NLMR")
library("rasterVis")

landscape <- nlm_mosaictess(200, 200, germs = 444)

levelplot(landscape, , margin = FALSE)
```

Another nice function from **rasterVis** is `gplot()`, a wrapper to use ggplot2
with raster data without reshaping your data as long data.frame:

```{r , fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library("NLMR")
library("rasterVis")

landscape <- nlm_mosaictess(200, 200, germs = 444)

gplot(landscape) + 
  geom_tile(aes(fill = value)) + 
  coord_equal()
```


### ggplot2

If you want to start from scratch with ggplot2:

```{r , fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library("NLMR")
library("raster")
library("ggplot2")

landscape <- nlm_mosaictess(200, 200, germs = 444)

# transform to long format for ggplot2
landscape_long <- as.data.frame(landscape, xy = TRUE)

# plot with ggplot2
ggplot(landscape_long, aes(x,y)) + 
  geom_tile(aes(fill = layer)) + 
  coord_equal()
```

### raster + plot()

... if you are in a lot of hurry, raster itself also has a plot method for raster:

```{r , fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library("NLMR")
library("raster")

landscape <- nlm_mosaictess(200, 200, germs = 444)

plot(landscape)
```

### Perspective plot

```{r fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library("raster")
library("NLMR")
landscape <- nlm_fbm(ncol = 50, nrow = 50, fract_dim = 1.3)

persp(landscape,
      exp=0.5,
      maxpixels = 5000,
      theta = 125,
      phi=45,
      xlab="Longitude",
      ylab="Latitude",
      zlab="Z",
      shade = 0.45)
```

### Contour plots

```{r fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library("NLMR")
library("rasterVis")

landscape <- nlm_mpd(ncol = 50, nrow = 50, roughness = 0.6)


contourplot(landscape,
            pretty = TRUE) 

levelplot(landscape,
          contour = TRUE,
          pretty = TRUE)
```


## Interactive plots

### rgl + rasterVis

```{r fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library("rgl")
library("rasterVis")
library("viridis")
library("NLMR")
landscape <- nlm_mpd(ncol = 100, nrow = 100, roughness = 0.6)

plot3D(landscape,
       zfac=2,
       lit=FALSE,
       col=colorRampPalette(magma(11)))

rglwidget()
```


### highcharter + plotly

```{r message=FALSE, warning=FALSE}
library("highcharter")
library("magrittr")
library("plotly")
library("NLMR")

# create a NLM to work with
landscape <- nlm_mosaicfield(ncol = 100, nrow = 100, n = 20)

# coerce to matrix
landscape_matrix <- raster::as.matrix(landscape)

# plot interactive graph
hchart(landscape_matrix) %>%
  # changing default color
  hc_colorAxis(stops = color_stops(colors = viridis::inferno(10))) %>%
  hc_exporting(
    enabled = TRUE
  )

# With plotly we can combine the interactive approach with the 3D Visualization
plot_ly(z = as.matrix(landscape_matrix), type = "surface", colors = viridis::magma(8))
```



### rayshader
```{r fig.height=7, fig.width=7, eval=FALSE}
library(rayshader)
library(NLMR)
library(raster)
library(rgl)

set.seed(123)

landscape <- nlm_mpd(1000, 1000, roughness = 0.6, rescale = FALSE) * 500
landscape <- raster::focal(landscape, w=matrix(1, 31, 31), mean, pad = TRUE, padValue=0)
landscape <- raster::as.matrix(landscape)

shadow = ray_shade(landscape,
                   zscale=1,
                   lambert=FALSE)
amb = ambient_shade(landscape,
                    zscale=1,
                    sunbreaks = 15, 
                    maxsearch = 100)

landscape %>%
  sphere_shade(zscale=5,texture = "imhof1") %>% 
  add_water(detect_water(landscape, min_area = 4000)) %>%
  add_shadow(shadow,0.7) %>%
  add_shadow(amb) %>%
  add_shadow(lamb_shade(landscape)) %>%
  plot_3d(landscape,
          zscale=5,
          fov=0,
          theta=-45,
          phi=45,
          windowsize=c(1200,1200),
          zoom=1.2,
          water=TRUE, 
          wateralpha = 0.8,
          watercolor = "lightblue",
          waterlinecolor = "white",
          waterlinealpha = 0.3,
          solid = FALSE) 
```
![](rayshader.gif)

```{r fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center',results='hide',fig.keep='all'}
library(rayshader)
library(NLMR)
library(raster)
library(rgl)

set.seed(123)

landscape <- nlm_mpd(1000, 1000, roughness = 0.6, rescale = FALSE) * 500
landscape <- raster::focal(landscape, w=matrix(1, 31, 31), mean, pad = TRUE, padValue=0)
landscape <- raster::as.matrix(landscape)

shadow = ray_shade(landscape,
                   zscale=1,
                   lambert=FALSE)
amb = ambient_shade(landscape,
                    zscale=1,
                    sunbreaks = 15, 
                    maxsearch = 100)

landscape %>%
  sphere_shade(zscale=5,texture = "imhof1") %>% 
  add_water(detect_water(landscape, min_area = 4000)) %>%
  add_shadow(shadow,0.7) %>%
  add_shadow(amb) %>%
  add_shadow(lamb_shade(landscape)) %>%
  plot_3d(landscape,
          zscale = 5,
          fov = 0,
          theta = -45,
          phi = 45,
          windowsize = c(1200, 1200),
          zoom = 1.2,
          water = TRUE,
          wateralpha = 0.8,
          watercolor = "lightblue",
          waterlinecolor = "white",
          waterlinealpha = 0.3,
          solid = TRUE,
          solidcolor = "grey75")

render_depth(
  focallength = 30,
  fstop = 2,
  bokehshape = "hex",
  bokehintensity = 5,
  progbar = FALSE
)
```
---
output: github_document
---

We try to maintain a list here of publications that have used the `NLMR` package. If you have used `NLMR` in your own work, we would love to hear from you. Please [file an issue on GitHub](https://github.com/marcosci/nlmr/issues/new/) so we can add your work to this list.

* Fletcher R., Fortin MJ. (2018) Land-Cover Pattern and Change. In: Spatial Ecology and Conservation Modeling. Springer, Cham. https://doi.org/10.1007/978-3-030-01989-1_3
* Langhammer, Maria, Jule Thober, Martin Lange, Karin Frank, and Volker Grimm. “Agricultural Landscape Generators for Simulation Models: A Review of Existing Solutions and an Outline of Future Directions.” Ecological Modelling 393 (February 2019): 135–51. https://doi.org/10.1016/j.ecolmodel.2018.12.010.
* Poggi S, Papaïx J, Lavigne C et al. Issues and challenges in landscape models for agriculture: from the representation of agroecosystems to the design of management strategies. Landscape Ecol. 2018;00:1-12. https://doi.org/10.1007/s10980-018-0699-8
* Nowosad, J., & Stepinski, T. (2018). Information-theoretical approach to measuring landscape complexity. https://doi.org/10.1101/383281
* Sciaini M, Fritsch M, Scherer C, Simpkins CE. NLMR and landscapetools: An integrated environment for simulating and modifying neutral landscape models in R. Methods Ecol Evol. 2018;00:1–9. https://doi.org/10.1111/2041-210X.13076
---
title: "NLMR Overview and Tips"
author: "Marco Sciaini"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{NLMR Overview and Tips}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r, echo = FALSE}
knitr::opts_chunk$set(
  cache = TRUE
)
```

## Overview

### Core models

<img src="https://wol-prod-cdn.literatumonline.com/cms/attachment/b963a726-ed88-4ede-863c-a65451f91d0f/mee313076-fig-0001-m.jpg"  width="100%" />

### Selection of possible merges

```{r, warning=FALSE, dpi = 200, fig.align='center'}

# 1.
edge_nlm <- nlm_edgegradient(100, 100)
distance_nlm <- nlm_distancegradient(100, 100, origin = c(20, 20,10, 10))
random_nlm <- nlm_random(100, 100)

# 2.
gauss_nlm <- nlm_gaussianfield(100, 100)
rectan_nlm <- nlm_randomrectangularcluster(100, 100, maxl = 30, minl = 10)

# 3.
mosaic_nlm <- nlm_mosaicfield(100, 100)

# 4.
planar_nlm <- nlm_planargradient(100, 100)
tess_nlm <- nlm_mosaictess(100, 100, germs = 200)

# plot it
landscapetools::show_landscape(list("a" = landscapetools::util_merge(edge_nlm, list(distance_nlm, random_nlm)),
                                    "b" = landscapetools::util_merge(gauss_nlm, rectan_nlm),
                                    "c" = landscapetools::util_merge(mosaic_nlm, list(random_nlm)),
                                    "d" = landscapetools::util_merge(planar_nlm, list(distance_nlm, tess_nlm))))
```

## Tips

If you are new to the raster package, I hope to collect here some useful tips
how to handle raster data in general. Furthermore, this section also serves
as a place to collect workflows on how to use **NLMR** and other R packages to
simulate specific patterns one can find in the literature.

### Basics

#### Counting cells

Ecologists are for example often interested in how much habitat one actually finds
in the study area you are looking at and there are a couple of nice ways to do
that with rasters. One of them is:

```{r}
library(NLMR)
library(raster)
library(dplyr)

landscape <- nlm_curds(curds = c(0.5, 0.3, 0.6),
                       recursion_steps = c(32, 6, 2),
                       wheyes = c(0.1, 0.05, 0.2))

# count cells for each category (0 = Matrix, 1 = Habitat)
landscape %>% 
  freq()
```

### More specific tips

#### Use matrix of parameter to simulate landscapes

```{r, fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library(NLMR)
library(landscapetools)
library(raster)
library(dplyr)
library(purrr)
library(tibble)

# simulation function that has the parameters we want to alter as input
simulate_landscape = function(roughness, weighting){
    nlm_mpd(ncol = 33,
            nrow = 33,
            roughness = roughness,
            rescale = TRUE) %>%
        util_classify(weighting = weighting)
}

# paramter combinations we are interested in
param_df = expand.grid(roughness = c(0.2, 0.9),
                       weighting = list(c(0.2, 0.8), c(0.2, 0.3, 0.5))) %>%
    as.tibble()

# map over the nested tibble and use each row as input for our simulation function
nlm_list = param_df %>% pmap(simulate_landscape)

# look at the results
show_landscape(nlm_list)
```


#### Simulate ecotones

Merging different types of NLMs, such as a planar gradient with a less autocorrelated landscape, provide a means of generating more complex landscapes and realistic-looking ecotones (Travis 2004):

```{r, fig.height=7, fig.width=7, message=FALSE, warning=FALSE, fig.align='center'}
library(NLMR)
library(landscapetools)

# landscape with higher autocorrelation 
high_autocorrelation <- nlm_edgegradient(ncol = 100, nrow = 100, direction = 80)

# landscape with lower autocorrelation 
low_autocorrelation <- nlm_fbm(ncol = 100, nrow = 100, fract_dim = 0.5)

# merge to derive ecotone
ecotones <- util_merge(low_autocorrelation, high_autocorrelation)

# look at the results
show_landscape(list("Low autocorrelation" = low_autocorrelation,
                    "High autocorrelation" = high_autocorrelation,
                    "Ecotones" = ecotones
                    ))
```

---
title: "NLMR Software Heritage"
author: "Marco Sciaini"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{NLMR Software Heritage}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

One of the compelling reasons to actually use NLMR is to substitute any of the
other software out there to simulate neutral landscape models. 
This vignette therefore aims to be a collection of code snippets that can be used
to imitate the functionality of previous tools (e.g. QRULE).

The list is very much work in progress, if this is of concerning to you it makes
probably sense to revisit here in the future.

## Replicate maps

One of the reasons to actually use neutral landscape models is the statistical
comparison and testing and development of landscape metrics. Therefore, one needs
replicates 

### Built in parameters

The neutral landscape models from **NLMR** can be simulated with a varying 
number of parameters to control for spatial autocorrelation. Simulating 
replicates that have the same properties as the parameter of interest can 
therefore be achieved like this:

```{r}
library(NLMR)
library(landscapetools)
library(purrr)
library(raster)

# write a simulation function with the desired parameters
simulate_landscape = function(x) {
  nlm_randomcluster(
  ncol = 30,
  nrow = 30,
  p = 0.4,
  ai = c(0.25, 0.25, 0.5),
  rescale = FALSE
  )
}

# rerun it for example 5 times
landscape_list <- rerun(5, simulate_landscape()) 

# look at the result
show_landscape(stack(landscape_list))
```

### Landscape metrics

If you are interested in landscapes that share a metric which is not a built-in 
parameter, the most clever way I can came up with is to simulate models as long as 
it takes to have the desired number of landscapes with a certain metric.

An exemplary workflow for this could look like this:

```{r warning=FALSE}
library(NLMR)
library(landscapetools)
library(landscapemetrics)
library(dplyr)

# simulation helpers
n <- 1 # counter and index variable
sim_results <- list() # list to store simulation results

# loop until we have 5 landscapes with the metric we are interested in
while (n < 6) {
  # In this use case we are interested in a categorical metric,
  # which is why we reclassify the continous result in three categories
  landscape <- nlm_mosaictess(100, 100, germs = 50) %>%
  util_classify(n = 3)
  
  # We are interested in the Euclidean Nearest Neighbor Distance Distribution
  enn_value <- lsm_l_enn_mn(landscape) %>%
  pull(value)
  
  # ... and we want to keep simulation results that have a mean ENN
  # between 7.5 and 8
  if (enn_value > 7.5 | enn_value < 8) {
    sim_results[n] <- landscape
    n = n + 1
  }
}

# look at the result
show_landscape(stack(sim_results))
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_planargradient.R
\name{nlm_planargradient}
\alias{nlm_planargradient}
\title{nlm_planargradient}
\usage{
nlm_planargradient(ncol, nrow, resolution = 1, direction = NA, rescale = TRUE)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{direction}{[\code{numerical(1)}]\cr
Direction of the gradient in degrees, if unspecified the direction is randomly
determined.}

\item{rescale}{[\code{logical(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\value{
RasterLayer
}
\description{
Simulates a planar gradient neutral landscape model.
}
\details{
Simulates a linear gradient sloping in a specified or random direction.
}
\examples{
# simulate planar gradient
planar_gradient <- nlm_planargradient(ncol = 200, nrow = 200)

\dontrun{
# visualize the NLM
landscapetools::show_landscape(planar_gradient)
}

}
\references{
Palmer, M.W. (1992) The coexistence of species in fractal landscapes.
\emph{The American Naturalist}, 139, 375 - 397.
}
\seealso{
\code{\link{nlm_distancegradient}},
\code{\link{nlm_edgegradient}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_randomrectangularcluster.R
\name{nlm_randomrectangularcluster}
\alias{nlm_randomrectangularcluster}
\title{nlm_randomrectangularcluster}
\usage{
nlm_randomrectangularcluster(
  ncol,
  nrow,
  resolution = 1,
  minl,
  maxl,
  rescale = TRUE
)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr Resolution of the raster.}

\item{minl}{[\code{numerical(1)}]\cr The minimum possible width and height for each random rectangular cluster.}

\item{maxl}{[\code{numerical(1)}]\cr The maximum possible width and height for each random rectangular cluster.}

\item{rescale}{[\code{logical(1)}]\cr If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\value{
RasterLayer
}
\description{
Simulates a random rectangular clusters neutral landscape model with values ranging 0-1.
}
\details{
The random rectangular cluster algorithm starts to fill a raster randomly
with rectangles defined by \code{minl} and \code{maxl} until the surface
of the landscape is completely covered.
This is one type of realisation of a "falling/dead leaves" algorithm,
for more details see Galerne & Gousseau (2012).
}
\examples{
# simulate random rectangular cluster
randomrectangular_cluster <- nlm_randomrectangularcluster(ncol = 50,
                                                          nrow = 30,
                                                          minl = 5,
                                                          maxl = 10)
\dontrun{
# visualize the NLM
landscapetools::show_landscape(randomrectangular_cluster)
}

}
\references{
Gustafson, E.J. & Parker, G.R. (1992). Relationships between landcover
proportion and indices of landscape spatial pattern. \emph{Landscape ecology},
7, 101–110.
Galerne B. & Gousseau Y. (2012). The Transparent Dead Leaves Model. Advances in
Applied Probability, \emph{Applied Probability Trust}, 44, 1–20.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/NLMR.R
\docType{package}
\name{NLMR-package}
\alias{NLMR}
\alias{NLMR-package}
\title{Simulating Neutral Landscape Models}
\description{
\emph{NLMR} is an R package for simulating neutral landscape models (NLMs).
}
\details{
This package contains vignettes that introduce NLM and basic usage
of the \emph{NLMR} package. The vignettes in this package are listed below.

\describe{
\item{\href{https://ropensci.github.io/NLMR/articles/getstarted.html}{
Quickstart Guide}}{Short walk-through of the \emph{NLMR} package and how to
handle the simulations.}
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://ropensci.github.io/NLMR/}
  \item Report bugs at \url{https://github.com/ropensci/NLMR/issues/}
}

}
\author{
\strong{Maintainer}: Marco Sciaini \email{marco.sciaini@posteo.net} (\href{https://orcid.org/0000-0002-3042-5435}{ORCID})

Authors:
\itemize{
  \item Matthias Fritsch \email{matthias.fritsch@forst.uni-goettingen.de}
  \item Maximilian Hesselbarth \email{mhk.hesselbarth@gmail.com}
  \item Craig Simpkins \email{simpkinscraig063@gmail.com} (\href{https://orcid.org/0000-0003-3212-1379}{ORCID})
  \item Cédric Scherer \email{cedricphilippscherer@gmail.com} (\href{https://orcid.org/0000-0003-0465-2543}{ORCID})
  \item Sebastian Hanß (\href{https://orcid.org/0000-0002-3990-4897}{ORCID})
}

Other contributors:
\itemize{
  \item Laura Graham (Laura reviewed the package for rOpenSci, see https://github.com/ropensci/onboarding/issues/188) [reviewer]
  \item Jeffrey Hollister (Jeffrey reviewed the package for rOpenSci, see https://github.com/ropensci/onboarding/issues/188) [reviewer]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_mosaicgibbs.R
\name{nlm_mosaicgibbs}
\alias{nlm_mosaicgibbs}
\title{nlm_mosaicgibbs}
\usage{
nlm_mosaicgibbs(
  ncol,
  nrow,
  resolution = 1,
  germs,
  R,
  patch_classes,
  rescale = TRUE
)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{germs}{[\code{numerical(1)}]\cr
Intensity parameter (non-negative integer).}

\item{R}{[\code{numerical(1)}]\cr
Interaction radius (non-negative integer) for the fitting of the spatial point
pattern process - the min. distance between germs in map units.}

\item{patch_classes}{[\code{numerical(1)}]\cr
Number of classes for germs.}

\item{rescale}{[\code{logical(1)}]\cr If \code{TRUE} (default), the values
are rescaled between 0-1.}
}
\value{
RasterLayer
}
\description{
Simulate a neutral landscape model using the Gibbs algorithm introduced in Gaucherel (2008).
}
\details{
\code{nlm_mosaicgibbs} offers the second option of simulating a neutral landscape model
described in Gaucherel (2008).
The method works in principal like the tessellation method (\code{nlm_mosaictess}),
but instead of a random point pattern the algorithm fits a simulated realization of the Strauss
process. The Strauss process starts with a given number of points and
uses a minimization approach to fit a point pattern with a given interaction
parameter (0 - hardcore process; 1 - Poisson process) and interaction radius
(distance of points/germs being apart).
}
\examples{
# simulate polygonal landscapes
mosaicgibbs <- nlm_mosaicgibbs(ncol = 40,
                              nrow = 30,
                              germs = 20,
                              R = 0.02,
                              patch_classes = 12)

\dontrun{
# visualize the NLM
landscapetools::show_landscape(mosaicgibbs)
}

}
\references{
Gaucherel, C. (2008) Neutral models for polygonal landscapes with linear
networks. \emph{Ecological Modelling}, 219, 39 - 48.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_mosaicfield.R
\name{nlm_mosaicfield}
\alias{nlm_mosaicfield}
\title{nlm_mosaicfield}
\usage{
nlm_mosaicfield(
  ncol,
  nrow,
  resolution = 1,
  n = 20,
  mosaic_mean = 0.5,
  mosaic_sd = 0.5,
  collect = FALSE,
  infinit = FALSE,
  rescale = TRUE
)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{n}{[\code{numerical(1)}]\cr
Number of steps over which the mosaic random field algorithm is run}

\item{mosaic_mean}{[\code{numerical(1)}]\cr
Mean value of the mosaic displacement distribution}

\item{mosaic_sd}{[\code{numerical(1)}]\cr
Standard deviation of the mosaic displacement distribution}

\item{collect}{[\code{logical(1)}]\cr
return \code{RasterBrick} of all steps 1:\code{n}}

\item{infinit}{[\code{logical(1)}]\cr
return raster of the random mosaic field algorithm with infinite steps}

\item{rescale}{[\code{logical(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\value{
RasterLayer or List with RasterLayer/s and/or RasterBrick
}
\description{
Simulates a mosaic random field neutral landscape model.
}
\examples{

# simulate mosaic random field
mosaic_field <- nlm_mosaicfield(ncol = 100,
                                nrow = 200,
                                n = NA,
                                infinit = TRUE,
                                collect = FALSE)
\dontrun{
# visualize the NLM
landscapetools::show_landscape(mosaic_field)
}

}
\references{
Schwab, Dimitri, Martin Schlather, and Jürgen Potthoff. "A general class of
mosaic random fields." arXiv preprint arXiv:1709.01441 (2017). \cr
Baddeley, Adrian, Ege Rubak, and Rolf Turner. Spatial point patterns:
methodology and applications with R. CRC Press, 2015.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_fbm.R
\name{nlm_fbm}
\alias{nlm_fbm}
\title{nlm_fbm}
\usage{
nlm_fbm(
  ncol,
  nrow,
  resolution = 1,
  fract_dim = 1,
  user_seed = NULL,
  rescale = TRUE,
  ...
)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{fract_dim}{[\code{numerical(1)}]\cr
The fractal dimension of the process (0,2)}

\item{user_seed}{[\code{numerical(1)}]\cr
Set random seed for the simulation}

\item{rescale}{[\code{numeric(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.}

\item{...}{Other options to RandomFields::RFoptions, especially if using
a fractal dimension between ~ 1.6 and 1.9 one must set the option
\code{modus_operandi = "sloppy"}.}
}
\value{
RasterLayer
}
\description{
Creates a two-dimensional fractional Brownian motion neutral landscape model.
}
\details{
Neutral landscapes are generated using fractional Brownian motion,
 an extension of Brownian motion in which the amount of correlation between
  steps is controlled by \code{frac_dim}. A high value of \code{frac_dim} produces a
   relatively smooth, correlated surface while a low value produces a rough, uncorrelated one.
}
\examples{
# simulate fractional brownian motion
fbm_raster  <- nlm_fbm(ncol = 20, nrow = 30, fract_dim = 0.8)

\dontrun{

# visualize the NLM
landscapetools::show_landscape(fbm_raster)

}

}
\references{
Travis, J.M.J. & Dytham, C. (2004). A method for simulating patterns of
habitat availability at static and dynamic range margins. \emph{Oikos} , 104, 410–416.

Martin Schlather, Alexander Malinowski, Peter J. Menck, Marco Oesting,
Kirstin Strokorb (2015). nlm_fBm. \emph{Journal of Statistical
Software}, 63(8), 1-25. URL http://www.jstatsoft.org/v63/i08/.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_neigh.R
\name{nlm_neigh}
\alias{nlm_neigh}
\title{nlm_neigh}
\usage{
nlm_neigh(
  ncol,
  nrow,
  resolution = 1,
  p_neigh,
  p_empty,
  categories = 3,
  neighbourhood = 4,
  proportions = NA,
  rescale = TRUE
)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{p_neigh}{[\code{numerical(1)}]\cr
Probability of a cell will turning into a value if there is any neighbor with the same or a
 higher value.}

\item{p_empty}{[\code{numerical(1)}]\cr
Probability a cell receives a value if all neighbors have no value (i.e.
zero).}

\item{categories}{[\code{numerical(1)}]\cr
Number of categories used.}

\item{neighbourhood}{[\code{numerical(1)}]\cr
The neighbourhood used to determined adjacent cells: `8 ("Moore")` takes the eight
surrounding cells, while `4 ("Von-Neumann")` takes the four adjacent cells
(i.e. left, right, upper and lower cells).}

\item{proportions}{[\code{vector(1)}]\cr
The algorithm uses uniform proportions for each category by default. A vector
with as many proportions as categories and that sums up to 1 can be used for
other distributions.}

\item{rescale}{[\code{logical(1)}]\cr If \code{TRUE} (default), the values
are rescaled between 0-1.}
}
\value{
RasterLayer
}
\description{
Create a neutral landscape model with categories and clustering
 based on neighborhood characteristics.
}
\details{
The algorithm draws a random cell and turns it into a given category based on
 the probabilities \code{p_neigh} and \code{p_empty}, respectively. The decision is
 based on the probability \code{p_neigh}, if there is any cell in the Moore- (8 cells) or
 Von-Neumann-neighborhood (4 cells), otherwise it is based on \code{p_empty}. To create
 clustered neutral landscape models, \code{p_empty} should be (significantly) smaller than
 \code{p_neigh}. By default, the Von-Neumann-neighborhood is used to check adjacent
 cells. The algorithm starts with the highest categorical value. If the
 proportion of cells with this value is reached, the categorical value is
 reduced by 1. By default, a uniform distribution of the categories is
 applied.
}
\examples{
# simulate neighborhood model
neigh_raster <- nlm_neigh(ncol = 50, nrow = 50, p_neigh = 0.7, p_empty = 0.1,
                    categories = 5, neighbourhood = 4)

\dontrun{
# visualize the NLM
landscapetools::show_landscape(neigh_raster)
}

}
\references{
Scherer, Cédric, et al. "Merging trait-based and individual-based modelling:
An animal functional type approach to explore the responses of birds to
climatic and land use changes in semi-arid African savannas."
\emph{Ecological Modelling} 326 (2016): 75-89.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_mpd.R
\name{nlm_mpd}
\alias{nlm_mpd}
\title{nlm_mpd}
\usage{
nlm_mpd(
  ncol,
  nrow,
  resolution = 1,
  roughness = 0.5,
  rand_dev = 1,
  torus = FALSE,
  rescale = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{roughness}{[\code{numerical(1)}]\cr
Controls the level of spatial autocorrelation (!= Hurst exponent)}

\item{rand_dev}{[\code{numerical(1)}]\cr
Initial standard deviation for the displacement step (default == 1), sets the
scale of the overall variance in the resulting landscape.}

\item{torus}{[\code{logical(1)}]\cr  Logical value indicating wether the algorithm should be simulated on a torus (default FALSE)}

\item{rescale}{[\code{logical(1)}]\cr If \code{TRUE} (default), the values
are rescaled between 0-1.}

\item{verbose}{[\code{logical(1)}]\cr If \code{TRUE} (default), the user gets
a warning that the functions changes the dimensions to an appropriate one for
the algorithm.}
}
\value{
RasterLayer
}
\description{
Simulates a midpoint displacement neutral landscape model.
}
\details{
The algorithm is a direct implementation of the midpoint displacement
algorithm.
It performs the following steps:

\itemize{
 \item{Initialization: }{ Determine the smallest fit of
 \code{max(ncol, nrow)} in \emph{n^2 + 1} and assign value to n.
 Setup matrix of size (n^2 + 1)*(n^2 + 1).
 Afterwards, assign a random value to the four corners of the matrix.}
 \item{Diamond Step: }{ For each square in the matrix, assign the average of
 the four corner points plus a random value to the midpoint of that square.}
 \item{Diamond Step: }{ For each diamond in the matrix, assign the average
  of the four corner points plus a random value to the midpoint of that
  diamond.}
}

At each iteration the roughness, an approximation to common Hurst exponent,
is reduced.
}
\examples{

# simulate midpoint displacement
midpoint_displacememt <- nlm_mpd(ncol = 100,
                                 nrow = 100,
                                 roughness = 0.3)
\dontrun{
# visualize the NLM
landscapetools::show_landscape(midpoint_displacememt)
}
}
\references{
\url{https://en.wikipedia.org/wiki/Diamond-square_algorithm}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_distancegradient.R
\name{nlm_distancegradient}
\alias{nlm_distancegradient}
\title{nlm_distancegradient}
\usage{
nlm_distancegradient(ncol, nrow, resolution = 1, origin, rescale = TRUE)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{origin}{[\code{numerical(4)}]\cr
Edge coordinates of the origin (raster::extent with xmin, xmax, ymin, ymax)
of the distance measurement.}

\item{rescale}{[\code{logical(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.
Otherwise, the distance in raster units is calculated.}
}
\value{
RasterLayer
}
\description{
Simulates a distance-gradient neutral landscape model.
}
\details{
The function takes the number of columns and rows as input and creates a
\code{RasterLayer} with the same extent. \code{Origin} is a numeric vector of
xmin, xmax, ymin, ymax for a rectangle inside the raster from which the
distance is measured.
}
\examples{

# simulate a distance gradient
distance_gradient <- nlm_distancegradient(ncol = 100, nrow = 100,
                                           origin = c(20, 30, 10, 15))
\dontrun{
# visualize the NLM
landscapetools::show_landscape(distance_gradient)
}
}
\seealso{
\code{\link{nlm_edgegradient}},
\code{\link{nlm_planargradient}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_edgegradient.R
\name{nlm_edgegradient}
\alias{nlm_edgegradient}
\title{nlm_edgegradient}
\usage{
nlm_edgegradient(ncol, nrow, resolution = 1, direction = NA, rescale = TRUE)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{direction}{[\code{numerical(1)}]\cr
Direction of the gradient (between 0 and 360 degrees), if unspecified the
direction is randomly determined.}

\item{rescale}{[\code{logical(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\value{
RasterLayer
}
\description{
Simulates an edge-gradient neutral landscape model.
}
\details{
Simulates a linear gradient orientated on a specified or random direction
that has a central peak running perpendicular to the gradient direction.
}
\examples{

# simulate random curdling
edge_gradient <- nlm_edgegradient(ncol = 100, nrow = 100, direction = 80)

\dontrun{
# visualize the NLM
landscapetools::show_landscape(edge_gradient)
}

}
\references{
Travis, J.M.J. & Dytham, C. (2004) A method for simulating patterns of
habitat availability at static and dynamic range margins. \emph{Oikos}, 104,
410–416.
}
\seealso{
\code{\link{nlm_distancegradient}},
\code{\link{nlm_planargradient}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_mosaictess.R
\name{nlm_mosaictess}
\alias{nlm_mosaictess}
\alias{nlm_polylands}
\title{nlm_mosaictess}
\usage{
nlm_mosaictess(ncol, nrow, resolution = 1, germs, rescale = TRUE)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{germs}{[\code{numerical(1)}]\cr
Intensity parameter (non-negative integer).}

\item{rescale}{[\code{logical(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\value{
RasterLayer
}
\description{
Simulate a neutral landscape model using the tesselation approach introduced in Gaucherel (2008).
}
\details{
\code{nlm_mosaictess} offers the first option of simulating a neutral landscape model
described in Gaucherel (2008). It generates a random point pattern (germs)
with an independent distribution and uses the Voronoi tessellation to simulate mosaic landscapes.
}
\examples{
# simulate polygonal landscapes
mosaictess <- nlm_mosaictess(ncol = 30, nrow = 60, germs = 200)

\dontrun{
# visualize the NLM
landscapetools::show_landscape(mosaictess)
}

}
\references{
Gaucherel, C. (2008) Neutral models for polygonal landscapes with linear
networks. \emph{Ecological Modelling}, 219, 39 - 48.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/util_rescale.R
\name{util_rescale}
\alias{util_rescale}
\title{util_rescale}
\usage{
util_rescale(x)
}
\arguments{
\item{x}{[\code{Raster* object}]}
}
\value{
Raster* object with values ranging from 0-1
}
\description{
Linearly rescale element values in a raster to a range between 0 and 1
}
\details{
Rasters generated by \code{nlm_*} functions are scaled between 0 and 1 as default, this option can be set to \code{FALSE} if needed.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_percolation.R
\name{nlm_percolation}
\alias{nlm_percolation}
\title{nlm_percolation}
\usage{
nlm_percolation(ncol, nrow, resolution = 1, prob = 0.5)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{prob}{[\code{numerical(1)}]\cr
Probability value for setting a cell to 1.}
}
\value{
RasterLayer
}
\description{
Generates a random percolation neutral landscape model.
}
\details{
The simulation of a random percolation map is accomplished in two steps:

\itemize{
 \item{Initialization: }{ Setup matrix of size (\code{ncol}*\code{nrow})}
 \item{Map generation: }{ For each cell in the matrix a single uniformly
 distributed random number is generated and tested against a probability
 \code{prob}. If the random number is smaller than \code{prob}, the cell is set to
 TRUE - if it is higher the cell is set to FALSE.}
}
}
\examples{
# simulate percolation model
percolation <- nlm_percolation(ncol = 100, nrow = 100, prob = 0.5)
\dontrun{
# visualize the NLM
landscapetools::show_landscape(percolation)
}
}
\references{
1. Gardner RH, O'Neill R V, Turner MG, Dale VH. 1989. Quantifying
scale-dependent effects of animal movement with simple percolation models.
\emph{ Landscape Ecology} 3:217 - 227.

2. Gustafson, E.J. & Parker, G.R. (1992) Relationships between landcover
proportion and indices of landscape spatial pattern. \emph{Landscape Ecology}
, 7, 101 - 110.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_randomcluster.R
\name{nlm_randomcluster}
\alias{nlm_randomcluster}
\title{nlm_randomcluster}
\usage{
nlm_randomcluster(
  ncol,
  nrow,
  resolution = 1,
  p,
  ai = c(0.5, 0.5),
  neighbourhood = 4,
  rescale = TRUE
)
}
\arguments{
\item{ncol}{[\code{integer(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{integer(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{p}{[\code{numerical(1)}]\cr
Defines the proportion of elements randomly selected to form
clusters.}

\item{ai}{Vector with the cluster type distribution (percentages of occupancy).
This directly controls the number of types via the given length.}

\item{neighbourhood}{[\code{numerical(1)}]\cr
Clusters are defined using a set of neighbourhood structures,
 4 (Rook's or von Neumann neighbourhood) (default), 8 (Queen's or Moore neighbourhood).}

\item{rescale}{[\code{logical(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\value{
Raster with random values ranging from 0-1.
}
\description{
Simulates a random cluster nearest-neighbour neutral landscape.
}
\details{
This is a direct implementation of steps A - D of the modified random clusters algorithm
by Saura & Martínez-Millán (2000), which creates naturalistic patchy patterns.
}
\examples{
# simulate random clustering
random_cluster <- nlm_randomcluster(ncol = 30, nrow = 30,
                                     p = 0.4,
                                     ai = c(0.25, 0.25, 0.5))
\dontrun{
# visualize the NLM
landscapetools::show_landscape(random_cluster)
}

}
\references{
Saura, S. & Martínez-Millán, J. (2000) Landscape patterns simulation with a
modified random clusters method. \emph{Landscape Ecology}, 15, 661 – 678.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_random.R
\name{nlm_random}
\alias{nlm_random}
\title{nlm_random}
\usage{
nlm_random(ncol, nrow, resolution = 1, rescale = TRUE)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{rescale}{[\code{logical(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\value{
RasterLayer
}
\description{
Simulates a spatially random neutral landscape model with values
drawn a uniform distribution.
}
\details{
The function takes the number of columns and rows as input and creates a
RasterLayer with the same extent. Each raster cell is randomly assigned a
value between 0 and 1 drawn from an uniform distribution (\code{runif(1,0,1)}).
}
\examples{
# simulate spatially random model
random <- nlm_random(ncol = 200, nrow = 100)

\dontrun{
# visualize the NLM
landscapetools::show_landscape(random)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_curds.R
\name{nlm_curds}
\alias{nlm_curds}
\title{nlm_curds}
\usage{
nlm_curds(curds, recursion_steps, wheyes = NULL, resolution = 1)
}
\arguments{
\item{curds}{[\code{numerical(x)}]\cr
Vector with percentage/s to fill with curds (fill with habitat (value ==
TRUE)).}

\item{recursion_steps}{[\code{numerical(x)}]\cr
Vector of successive cutting steps for the blocks (split 1 block into \code{x}
blocks).}

\item{wheyes}{[\code{numerical(x)}]\cr
Vector with percentage/s to fill with wheys, which fill matrix in an
additional step with habitat.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the resulting raster.}
}
\value{
raster
}
\description{
Simulates a random curd neutral landscape model with optional wheys.
}
\details{
Random curdling recursively subdivides the plane into blocks.
At each level of the recursion, a fraction of the blocks are declared as
habitat (value == TRUE) while the remaining blocks continue to be defined as matrix (value == FALSE) and enter the next recursive cycle.

The optional argument (\code{wheyes}) allows wheys to be added, in which a set proportion of cells that were
declared matrix (value == FALSE) during recursion, are now set as habitat cells (value == TRUE).

If \deqn{curds_{1} = curds_{2} = recursion_steps_{2} = ... = curds_{n} =
recursion_steps_{n}} the models resembles a binary random map.

Note that you can not set ncol and nrow with this landscape algorithm.
The amount of cells and hence dimension of the raster is given by the vector product of the recursive steps.
}
\examples{

# simulate random curdling
(random_curdling <- nlm_curds(curds = c(0.5, 0.3, 0.6),
                              recursion_steps = c(32, 6, 2)))

# simulate wheyed curdling
(wheyed_curdling <- nlm_curds(curds = c(0.5, 0.3, 0.6),
                              recursion_steps = c(32, 6, 2),
                              wheyes = c(0.1, 0.05, 0.2)))
\dontrun{
# Visualize the NLMs
landscapetools::show_landscape(random_curdling)
landscapetools::show_landscape(wheyed_curdling)
}

}
\references{
Keitt TH. 2000. Spectral representation of neutral landscapes.
\emph{Landscape Ecology} 15:479-493.

Szaro, Robert C., and David W. Johnston, eds. Biodiversity in managed
landscapes: theory and practice. \emph{Oxford University Press}, USA, 1996.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nlm_gaussianfield.R
\name{nlm_gaussianfield}
\alias{nlm_gaussianfield}
\title{nlm_gaussianfield}
\usage{
nlm_gaussianfield(
  ncol,
  nrow,
  resolution = 1,
  autocorr_range = 10,
  mag_var = 5,
  nug = 0.2,
  mean = 0.5,
  user_seed = NULL,
  rescale = TRUE
)
}
\arguments{
\item{ncol}{[\code{numerical(1)}]\cr
Number of columns forming the raster.}

\item{nrow}{[\code{numerical(1)}]\cr
Number of rows forming the raster.}

\item{resolution}{[\code{numerical(1)}]\cr
Resolution of the raster.}

\item{autocorr_range}{[\code{numerical(1)}]\cr
Maximum range (raster units) of spatial autocorrelation.}

\item{mag_var}{[\code{numerical(1)}]\cr
Magnitude of variation over the entire landscape.}

\item{nug}{[\code{numerical(1)}]\cr
Magnitude of variation in the scale of \code{autocorr_range},
smaller values lead to more homogeneous landscapes.}

\item{mean}{[\code{numerical(1)}]\cr
Mean value over the field.}

\item{user_seed}{[\code{numerical(1)}]\cr
Set random seed for the simulation}

\item{rescale}{[\code{numeric(1)}]\cr
If \code{TRUE} (default), the values are rescaled between 0-1.}
}
\description{
Simulates a spatially correlated random fields (Gaussian random
fields) neutral landscape model.
}
\details{
Gaussian random fields are a collection of random numbers on a spatially
discrete set of coordinates (landscape raster). Natural sciences often apply
them with spatial autocorrelation, meaning that objects which distant are more
distinct from one another than they are to closer objects.
}
\examples{
# simulate random gaussian field
gaussian_field <- nlm_gaussianfield(ncol = 90, nrow = 90,
                                    autocorr_range = 60,
                                    mag_var = 8,
                                    nug = 5)

\dontrun{
# visualize the NLM
landscapetools::show_landscape(gaussian_field)
}

}
\references{
Kéry & Royle (2016) \emph{Applied Hierarchical Modeling in Ecology}
Chapter 20
}
