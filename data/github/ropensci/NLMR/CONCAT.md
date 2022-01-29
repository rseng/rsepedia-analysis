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
