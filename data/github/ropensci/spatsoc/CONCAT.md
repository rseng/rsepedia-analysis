
<!-- badges: start -->

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![peer-review](https://badges.ropensci.org/237_status.svg)](https://github.com/ropensci/software-review/issues/237)
[![CRAN](https://www.r-pkg.org/badges/version/spatsoc)](https://cran.r-project.org/package=spatsoc)
[![](https://img.shields.io/badge/devel%20version-0.1.16-blue.svg)](https://github.com/robitalec/spatsoc)
[![cran
checks](https://cranchecks.info/badges/summary/spatsoc)](https://cran.r-project.org/web/checks/check_results_spatsoc.html)
[![codecov](https://codecov.io/gh/ropensci/spatsoc/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/spatsoc)
[![R-CMD-check](https://github.com/ropensci/spatsoc/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/spatsoc/actions)
<!-- badges: end -->

# spatsoc

### [News](#news) - [Installation](#installation) - [Usage](#usage) - [Contributing](#contributing)

`spatsoc` is an R package for detecting spatial and temporal groups in
GPS relocations. It can be used to convert GPS relocations to
gambit-of-the-group format to build proximity-based social networks with
grouping and edge-list generating functions. In addition, the
`randomizations` function provides data-stream randomization methods
suitable for GPS data and the `get_gbi` function generates group by
individual matrices useful for building networks with
`asnipe::get_network`.

See below for [installation](#installation) and basic [usage](#usage).

For more details, see the [blog
post](https://ropensci.org/blog/2018/12/04/spatsoc/) and vignettes:

-   [Introduction to
    spatsoc](https://docs.ropensci.org/spatsoc/articles/intro-spatsoc.html)
-   [Frequently asked
    questions](https://docs.ropensci.org/spatsoc/articles/faq.html)
-   [Using spatsoc in social network
    analysis](https://docs.ropensci.org/spatsoc/articles/using-in-sna.html)
-   [Using edge list and dyad id
    functions](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html)

## News

We wrote a [`targets`](https://github.com/ropensci/targets) workflow,
available at
[github.com/robitalec/targets-spatsoc-networks](https://github.com/ropensci/targets).
`targets` is an incredible package for designing workflows in R and,
with it, we can reproducibly run all steps from raw telemetry data to
output networks and metrics. Check it out and let us know how it works
for you!

Edge-list generating functions added:

-   `edge_nn`
-   `edge_dist`

and dyad id function:

-   `dyad_id`

(feedback welcome as always!)

Both documented further in a vignette: [Using edge list and dyad id
functions](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html).

Also, our article describing `spatsoc` is published at Methods in
Ecology and Evolution. [Link
here](https://doi.org/10.1111/2041-210X.13215). Thanks to reviewers and
editors at
[rOpenSci](https://github.com/ropensci/software-review/issues/237) and
at [MEE](https://besjournals.onlinelibrary.wiley.com/journal/2041210x).

More detailed news
[here](https://docs.ropensci.org/spatsoc/news/index.html).

## Installation

``` r
# Stable release
install.packages('spatsoc')

# Development version
remotes::install_github('ropensci/spatsoc')
```

`spatsoc` depends on `rgeos` and requires
[GEOS](https://trac.osgeo.org/geos/) installed on the system.

-   Debian/Ubuntu: `apt-get install libgeos-dev`
-   Arch: `pacman -S geos`
-   Fedora: `dnf install geos geos-devel`
-   Mac: `brew install geos`
-   Windows: see [here](https://trac.osgeo.org/osgeo4w/)

## Usage

### Load package, import data

`spatsoc` expects a `data.table` for all of its functions. If you have a
`data.frame`, you can use `data.table::setDT()` to convert it by
reference. If your data is a text file (e.g.: CSV), you can use
`data.table::fread()` to import it as a `data.table`.

``` r
library(spatsoc)
library(data.table)
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]
```

### Temporal grouping

`group_times` groups rows temporally using a threshold defined in units
of minutes (B), hours (C) or days (D).

<img src="man/figures/fig1.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

### Spatial grouping

`group_pts` groups points spatially using a distance matrix (B) and a
spatial threshold defined by the user (50m in this case). Combined with
`group_times`, the returned ‘group’ column represents spatiotemporal,
point based groups (D).

<img src="man/figures/fig2.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

`group_lines` groups sequences of points (forming a line) spatially by
buffering each line (A) by the user defined spatial threshold. Combined
with `group_times`, the returned ‘group’ column represents
spatiotemporal, line overlap based groups (B).

<img src="man/figures/fig3.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

`group_polys` groups home ranges by spatial and proportional overlap.
Combined with `group_times`, the returned ‘group’ column represents
spatiotemporal, polygon overlap based groups.

<img src="man/figures/fig4.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

### Edge-list generating functions

`edge_dist` and `edge_nn` generate edge-lists. `edge_dist` measures the
spatial distance between individuals (A) and returns all pairs within
the user specified distance threshold (B). `edge_nn` measures the
distance between individuals (C) and returns the nearest neighbour to
each individual (D).

<img src="man/figures/fig5.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

### Social network analysis functions

`randomizations` for data-stream randomization and `get_gbi` for
generating group by individual matrices.

# Contributing

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree
to abide by its terms.

Development of `spatsoc` welcomes contribution of feature requests, bug
reports and suggested improvements through the [issue
board](https://github.com/ropensci/spatsoc/issues).

See details in [CONTRIBUTING.md](CONTRIBUTING.md).

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# v 0.1.17 (unreleased)

* added a link to our `spatsoc` + `targets` workflow example
* changed the error and underlying check for `group_polys` from alphanumeric to
spaces in input DT's id column

# v 0.1.16 (2021-03-23)
* added an option for `edge_dist` to handle threshold = NULL. If NULL, `edge_dist` will return all neighbours observed (eg. useful if one wanted to calculated mean nearest neighbour distance at each timegroup). 
* updated EPSG argument according to newest recommendations in tests, man and vignettes ([PR 38](https://github.com/ropensci/spatsoc/pull/38)
* removed expect_silent tests ([PR 37](https://github.com/ropensci/spatsoc/pull/37))
* switched CI for tests and code coverage to GitHub Actions ([PR 36](https://github.com/ropensci/spatsoc/pull/36))


# v 0.1.15 (2020-10-21)
* fix TZ=UTC data.table tests ([Issue 32](https://github.com/ropensci/spatsoc/issues/32))

# v 0.1.14 (2020-07-03)
* updated tests, man and vignettes following new handling of projections in sp ([PR 31](https://github.com/ropensci/spatsoc/pull/31), [R spatial information](https://www.r-spatial.org/r/2020/03/17/wkt.html))
* clarified explicit drop of NAs in dyadID in edge list vignette

# v 0.1.13 (2020-03-25)
* added `dyad_id` function for generating dyad IDs with edge functions ([PR 27](https://github.com/ropensci/spatsoc/pull/25))
* added a vignette describing `edge_dist`, `edge_nn` and `dyad_id` functions [here](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html) ([PR 14](https://github.com/ropensci/spatsoc/pull/14))

# v 0.1.12 (2020-03-02)
* fixed `data.table` error in `edge_dist` and `edge_nn` ([PR 25](https://github.com/ropensci/spatsoc/pull/25))


# v 0.1.11 (2020-02-20)
* removed default NULL from 'timegroup' arguments in `group_pts`, `edge_dist` and `edge_nn` ([PR 24](https://github.com/ropensci/spatsoc/pull/24))


# v 0.1.10 (2019-06-06)
* added optional return of distance between individuals with `edge_dist` ([PR 19](https://github.com/ropensci/spatsoc/pull/19)) and `edge_nn` ([PR 21](https://github.com/ropensci/spatsoc/pull/21))


# v 0.1.9 (2019-05-14)
* fixed bug for randomizations type 'step' and 'daily' ([PR 13](https://github.com/ropensci/spatsoc/pull/13)). 
* clarified `SIMPLIFY=FALSE` in SNA vignette. 


# v 0.1.8 (2019-04-05)
* update [FAQ](https://docs.ropensci.org/spatsoc/articles/faq.html) and [Introduction to spatsoc](https://docs.ropensci.org/spatsoc/articles/intro-spatsoc.html) vignettes adding entries for edge list generating functions. 
* added edge list generating function `edge_nn` ([PR 11](https://github.com/ropensci/spatsoc/pull/12))
* added edge list generating function `edge_dist` ([PR 11](https://github.com/ropensci/spatsoc/pull/11))


# v 0.1.7 (2019-03-26)
* fix inconsistent blocks across years ([PR 10](https://github.com/ropensci/spatsoc/pull/10))
* update FAQ: remove old randomizations notes, clarify group_times block


# v 0.1.6 (2019-01-10)
* fix bug 'group_times misses nearest hour with mins threshold' ([#5](https://github.com/ropensci/spatsoc/issues/5) and [PR 6](https://github.com/ropensci/spatsoc/pull/6))

# v 0.1.5 (2018-12-04)
* update issue labels and contributing
* change over issue board location from GitLab to rOpenSci repository on GitHub
* added preprint CITATION
* added "https://" to `pkgdown` URL ([PR 1](https://github.com/ropensci/spatsoc/pull/1))

# v 0.1.4 (2018-10-26)
* fin [rOpenSci onboarding process](https://github.com/ropensci/software-review/issues/237)
* fixed bug couldn't provide percent to kernel type `build_polys` or `group_polys`([!3](https://gitlab.com/robit.a/spatsoc/-/merge_requests/3))


# v 0.1.3 
* added `get_gbi` to generate group by individual matrices for better integrating `spatsoc` in social network analysis workflows ([!2](https://gitlab.com/robit.a/spatsoc/-/merge_requests/2))


# v 0.1.2

* **major change to randomizations**: when `iterations = 1`, `randomizations` no longer returns the DT with appended columns. Regardless of the value of iterations, `randomizations` always returns observed rows followed by randomized rows in a long `data.table` ([!1](https://gitlab.com/robit.a/spatsoc/-/merge_requests/1)). 

# v 0.1.1 (2018-09-17)

* improvements to package, function documentation
* [FAQ](https://docs.ropensci.org/spatsoc/articles/faq.html) vignette added
* fixed `build_lines` ordering bug to ensure rows are ordered by date time when building lines
* added CODE_OF_CONDUCT.md and CONTRIBUTING.md
* [Using spatsoc in social network analysis](https://docs.ropensci.org/spatsoc/articles/using-in-sna.html) vignette added

# v 0.1.0 (2018-07-20)

## Initial release

* temporal grouping function: `group_times`
* spatial grouping functions: `group_pts`, `group_lines`, `group_polys`
* data-stream randomization function: `randomizations`
* spatial build functions: `build_lines`, `build_polys`
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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
# Contributing

Development of `spatsoc` welcomes contribution of feature requests, bug reports and suggested improvements through the [issue board](https://github.com/ropensci/spatsoc/issues). 

### Issues

Use the labels available on the Issue Board:

#### priority
 
| priority  | description | colour                                          	| 
|-----------|-------------|---------------------------------------------------|
| critical 	|             | <span style="background:#e50000">#e50000</span> 	|
| high 	    |             | <span style="background:#fc7436">#fc7436</span> 	|
| medium 	  |             | <span style="background:#ffbb5c">#ffbb5c</span> 	|
| low 	    |             | <span style="background:#dfdfdf">#dfdfdf</span> 	|


#### status

| status      	| description                                	| colour                                            |
|-------------	|--------------------------------------------	|--------------------------------------------------	|
| completed   	|                                            	| <span style="background:#34495E">#34495E</span> 	|
| confirmed   	| bug is reproducible                        	| <span style="background:#e50000">#e50000</span> 	|
| in progress 	| [workinonit](https://youtu.be/5nO7IA1DeeI) 	| <span style="background:#aefd3d">#aefd3d</span> 	|
| duplicate   	|                                            	| <span style="background:#dfdfdf">#dfdfdf</span> 	|

#### type

|  type         	| description                        | colour                                           |
|---------------	|----------------------------------- |-------------------------------------------------	|
| documentation 	|             	                     | <span style="background:#639f3b">#639f3b</span> 	|
| support       	| usage or implementation troubles   | <span style="background:#e48ec0">#e48ec0</span>  |
| bug           	| not working as described, intended | <span style="background:#e50000">#e50000</span> 	|
| discussion    	|             	                     | <span style="background:#47b6a3">#47b6a3</span> 	|
| enhancement   	|             	                     | <span style="background:#8935b8">#8935b8</span> 	|
| typo          	|             	                     | <span style="background:#dfdfdf">#dfdfdf</span> 	|


### Prefer to Email? 

Email Alec Robitaille (see DESCRIPTION). 

Though note that private discussions over email don't help others - of course email is totally warranted if it's a sensitive problem of any kind.

### Code of Conduct

Please note that the `spatsoc` project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### Thanks for contributing!

This contributing guide is adapted from the tidyverse contributing guide available at https://raw.githubusercontent.com/r-lib/usethis/master/inst/templates/tidy-contributing.md
## v 0.1.16
See NEWS.md for all updates. 

* added an option for `edge_dist` to handle threshold = NULL. If NULL, `edge_dist` will return all neighbours observed (eg. useful if one wanted to calculated mean nearest neighbour distance at each timegroup). 
* updated EPSG argument according to newest recommendations in tests, man and vignettes ([PR 38](https://github.com/ropensci/spatsoc/pull/38)
* removed expect_silent tests ([PR 37](https://github.com/ropensci/spatsoc/pull/37))
* switched CI for tests and code coverage to GitHub Actions ([PR 36](https://github.com/ropensci/spatsoc/pull/36))

## Test environments
* windows-latest (release) 
* macOS-latest (release)
* ubuntu-20.04 (release)
* ubuntu-20.04 (devel)
* ubuntu 16.04 (3.5)

## R CMD check results

There were no ERRORs, WARNING or NOTES
---
output: github_document
---

```{r opts, include = FALSE}
knitr::opts_chunk$set(message = FALSE, 
                      warning = FALSE,
                      eval = FALSE, 
                      echo = TRUE)
```

<!-- badges: start -->
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
[![peer-review](https://badges.ropensci.org/237_status.svg)](https://github.com/ropensci/software-review/issues/237)
[![CRAN](https://www.r-pkg.org/badges/version/spatsoc)](https://cran.r-project.org/package=spatsoc)
`r badger::badge_devel("robitalec/spatsoc", "blue")`
[![cran checks](https://cranchecks.info/badges/summary/spatsoc)](https://cran.r-project.org/web/checks/check_results_spatsoc.html)
[![codecov](https://codecov.io/gh/ropensci/spatsoc/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/spatsoc)
[![R-CMD-check](https://github.com/ropensci/spatsoc/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/spatsoc/actions)
<!-- badges: end -->

# spatsoc 

### [News](#news) - [Installation](#installation) - [Usage](#usage) - [Contributing](#contributing)

`spatsoc` is an R package for detecting spatial and temporal groups in GPS relocations. It can be used to convert GPS relocations to gambit-of-the-group format to build proximity-based social networks with grouping and edge-list generating functions. In addition, the `randomizations` function provides data-stream randomization methods suitable for GPS data and the `get_gbi` function generates group by individual matrices useful for building networks with `asnipe::get_network`. 

See below for [installation](#installation) and basic [usage](#usage). 

For more details, see the [blog post](https://ropensci.org/blog/2018/12/04/spatsoc/) and vignettes:

* [Introduction to spatsoc](https://docs.ropensci.org/spatsoc/articles/intro-spatsoc.html)
* [Frequently asked questions](https://docs.ropensci.org/spatsoc/articles/faq.html)
* [Using spatsoc in social network analysis](https://docs.ropensci.org/spatsoc/articles/using-in-sna.html)
* [Using edge list and dyad id functions](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html)


## News
We wrote a [`targets`](https://github.com/ropensci/targets) workflow, available at [github.com/robitalec/targets-spatsoc-networks](https://github.com/ropensci/targets). `targets` is an incredible package for designing workflows in R and, with it, we can reproducibly run all steps from raw telemetry data to output networks and metrics. Check it out and let us know how it works for you!


Edge-list generating functions added:

  * `edge_nn`
  * `edge_dist`
  
and dyad id function:

  * `dyad_id`
  
(feedback welcome as always!)

Both documented further in a vignette: [Using edge list and dyad id functions](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html).


Also, our article describing `spatsoc` is published at Methods in Ecology and Evolution. [Link here](https://doi.org/10.1111/2041-210X.13215). Thanks to reviewers and editors at [rOpenSci](https://github.com/ropensci/software-review/issues/237) and at [MEE](https://besjournals.onlinelibrary.wiley.com/journal/2041210x). 


More detailed news [here](https://docs.ropensci.org/spatsoc/news/index.html).


## Installation

```{r install}
# Stable release
install.packages('spatsoc')

# Development version
remotes::install_github('ropensci/spatsoc')
```

`spatsoc` depends on `rgeos` and requires [GEOS](https://trac.osgeo.org/geos/) installed on the system. 

* Debian/Ubuntu: `apt-get install libgeos-dev`
* Arch: `pacman -S geos`
* Fedora: `dnf install geos geos-devel`
* Mac: `brew install geos`
* Windows: see [here](https://trac.osgeo.org/osgeo4w/)


## Usage
### Load package, import data
`spatsoc` expects a `data.table` for all of its functions. If you have a `data.frame`, you can use `data.table::setDT()` to convert it by reference. If your data is a text file (e.g.: CSV), you can use `data.table::fread()` to import it as a `data.table`. 

```{r library}
library(spatsoc)
library(data.table)
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]
```


### Temporal grouping
`group_times` groups rows temporally using a threshold defined in units of minutes (B), hours (C) or days (D). 

<img src="man/figures/fig1.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

### Spatial grouping 
`group_pts` groups points spatially using a distance matrix (B) and a spatial threshold defined by the user (50m in this case). Combined with `group_times`, the returned 'group' column represents spatiotemporal, point based groups (D). 


<img src="man/figures/fig2.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

`group_lines` groups sequences of points (forming a line) spatially by buffering each line (A) by the user defined spatial threshold. Combined with `group_times`, the returned 'group' column represents spatiotemporal, line overlap based groups (B). 

<img src="man/figures/fig3.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>



`group_polys` groups home ranges by spatial and proportional overlap. Combined with `group_times`, the returned 'group' column represents spatiotemporal, polygon overlap based groups. 


<img src="man/figures/fig4.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>


### Edge-list generating functions
`edge_dist` and `edge_nn` generate edge-lists. `edge_dist` measures the spatial distance between individuals (A) and returns all pairs within the user specified distance threshold (B). `edge_nn` measures the distance between individuals (C) and returns the nearest neighbour to each individual (D). 

<img src="man/figures/fig5.png" style="max-height:400px; display:block; margin-left: auto; margin-right: auto;"/>

### Social network analysis functions
`randomizations` for data-stream randomization and `get_gbi` for generating group by individual matrices. 




# Contributing
Please note that this project is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
  
Development of `spatsoc` welcomes contribution of feature requests, bug reports and suggested improvements through the [issue board](https://github.com/ropensci/spatsoc/issues). 

See details in [CONTRIBUTING.md](CONTRIBUTING.md). 

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Using spatsoc in social network analysis - grouping functions"
author: "Alec Robitaille, Quinn Webber and Eric Vander Wal"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    number_sections: false
    toc: false
vignette: >
  %\VignetteIndexEntry{Using spatsoc for social network analysis}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---

```{r knitropts, include = FALSE}
knitr::opts_chunk$set(message = TRUE, 
                      warning = FALSE,
                      eval = FALSE, 
                      echo = TRUE)
```


`spatsoc` can be used in social network analysis to generate gambit of the group format data from GPS relocation data, perform data stream randomization and generate group by individual matrices. 

Gambit of the group format data is generated using the grouping functions:

* `group_times`
* `group_pts`
* `group_lines`
* `group_polys`

Data stream randomization is performed using the `randomizations` function. 

Group by individual matrices are generated using the `get_gbi` function. 


**Note**: edge list generating functions are also available and are described in the vignette [Using edge list generating functions and dyad_id](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html). 


# Generate gambit of the group data
spatsoc provides users with one temporal (`group_times`) and three spatial (`group_pts`, `group_lines`, `group_polys`) functions to generate gambit of the group data from GPS relocations. Users can consider spatial grouping at three different scales combined with an appropriate temporal grouping threshold. The gambit of the group data is then used to generate a group by individual matrix and build the network. 


## 1. Load packages and prepare data
`spatsoc` expects a `data.table` for all `DT` arguments and date time columns to be formatted `POSIXct`. 

```{r}
## Load packages
library(spatsoc)
library(data.table)
library(asnipe)
library(igraph)

## Read data as a data.table
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

## Cast datetime column to POSIXct
DT[, datetime := as.POSIXct(datetime)]

## Calculate the year of the relocation 
DT[, yr := year(datetime)]
```


Next, we will group relocations temporally with `group_times` and spatially with one of `group_pts`, `group_lines`, `group_polys`. Note: these are mutually exclusive, only select one spatial grouping function at a time. 

## 2. a) `group_pts` 

Point based grouping by calculating distance between relocations in each timegroup. Depending on species and study system, relevant temporal and spatial grouping thresholds are used. In this case, relocations within 5 minutes and 50 meters are grouped together. 

```{r}
## Temporal groups
group_times(DT, datetime = 'datetime', threshold = '5 minutes')

## Spatial groups
group_pts(
  DT,
  threshold = 50,
  id = 'ID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup'
)

```

## 2. b) `group_lines`

Line based grouping by measuring intersection of, optionally buffered, trajectories for each individual in each timegroup. Longer temporal thresholds are used to measure, for example, intersecting daily trajectories. 

```{r, eval = FALSE}
# EPSG code for relocations
utm <- 'EPSG:32736'

## Group relocations by julian day
group_times(DT, datetime = 'datetime', threshold = '1 day')

## Group lines for each individual and julian day
group_lines(
  DT,
  threshold = 50,
  projection = utm,
  id = 'ID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  sortBy = 'datetime'
)
```


## 2. c) `group_polys`

Polygon based grouping by generating home ranges using `adehabitatHR` and measuring intersection or proportional overlap. Longer temporal thresholds are used to create seasonal, monthly, yearly home ranges.

```{r, eval = FALSE}
# EPSG code for relocations
utm <- 'EPSG:32736'

## Option 1: area = FALSE and home range intersection 'group' column added to DT 
group_polys(
  DT,
  area = FALSE,
  hrType = 'mcp',
  hrParams = list(percent = 95),
  projection = utm,
  id = 'ID',
  coords = c('X', 'Y')
)

## Option 2: area = TRUE 
#  results must be assigned to a new variable 
#  data.table returned has ID1, ID2 and proportion and area overlap
areaDT <- group_polys(
  DT,
  area = TRUE,
  hrType = 'mcp',
  hrParams = list(percent = 95),
  projection = utm,
  id = 'ID',
  coords = c('X', 'Y')
)

```

# Build observed network 
Once we've generated groups using `group_times` and one of the spatial grouping functions, we can generate a group by individual matrix. 

The following code chunk showing `get_gbi` can be used for outputs from any of `group_pts`, `group_lines` or `group_polys(area = FALSE)`. For the purpose of this vignette however, we will consider the outputs from `group_pts` ([2. a)](#a-group_pts)) for the following code chunk.

Note: we show this example creating the group by individual matrix and network for only 2016 to illustrate how `spatsoc` can be used for simpler data with no splitting of temporal or spatial subgroups (e.g.: yearly, population). See the random network section for how to use `spatsoc` in social network analysis for multi-year or other complex data. 

## 3. `get_gbi`
```{r}
## Subset DT to only year 2016
subDT <- DT[yr == 2016]

## Generate group by individual matrix
# group column generated by spatsoc::group_pts
gbiMtrx <- get_gbi(DT = subDT, group = 'group', id = 'ID')
```

Note: `spatsoc::get_gbi` is identical in function to `asnipe::get_group_by_individual`, but is more efficient (some benchmarks measuring >10x improvements) thanks to `data.table::dcast`.


## 4. `asnipe::get_network`
Next, we can use `asnipe::get_network` to build the observed social network. Ensure that the argument "data_format" is "GBI". Use other arguments that are relevant to your analysis, here we calculate a Simple ratio index. 

```{r}
## Generate observed network
net <- get_network(gbiMtrx,
                   data_format = "GBI",
                   association_index = "SRI")
```


# Data stream randomization
Three types of data stream randomization are provided by `spatsoc`'s `randomizations` function:

* step: randomizes identities of relocations between individuals within each time step.
* daily: randomizes identities of relocations between individuals within each day.
* trajectory: randomizes daily trajectories within individuals ([Spiegel et al. 2016](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12553)).

The results of `randomizations` must be assigned. The function returns the `id` and `datetime` columns provided (and anything provided to `splitBy`). In addition, columns 'observed' and 'iteration' are returned indicating observed rows and which iteration rows correspond to (where 0 is the observed). 

As with spatial grouping functions, these methods are mutually exclusive. Pick one `type` and rebuild the network after randomization. 

Note: the `coords` argument is only required for trajectory type randomization, since after randomizing with this method, the 'coords' are needed to redo spatial grouping (with `group_pts`, `group_lines` or `group_polys`). 


## 5. a) `type = 'step'`
`'step'` randomizes identities of relocations between individuals within each time step. The `datetime` argument expects an integer group generated by `group_times`. The `group` argument expects the column name of the group generated from the spatial grouping functions. 

Four columns are returned when `type = 'step'` along with `id`, `datetime` and `splitBy` columns:

* 'randomID' - randomly selected ID from IDs within each time step
* 'observed' - observed rows (TRUE/FALSE)
* 'iteration' - which iteration rows correspond to (0 is observed)

```{r}
# Calculate year column to ensure randomization only occurs within years since data spans multiple years
DT[, yr := year(datetime)]

## Step type randomizations
#  providing 'timegroup' (from group_times) as datetime
#  splitBy = 'yr' to force randomization only within year
randStep <- randomizations(
   DT,
   type = 'step',
   id = 'ID',
   group = 'group',
   coords = NULL,
   datetime = 'timegroup',
   iterations = 3,
   splitBy = 'yr'
)
```


## 5. b) `type = 'daily'`
`'daily'` randomizes identities of relocations between individuals within each day. The `datetime` argument expects a datetime `POSIXct` format column. 

Four columns are returned when `type = 'daily'` along with `id`, `datetime` and `splitBy` columns:

* 'randomID' - randomly selected ID for each day
* 'jul' - julian day
* 'observed' - observed rows (TRUE/FALSE)
* 'iteration' - which iteration rows correspond to (0 is observed)

```{r}
# Calculate year column to ensure randomization only occurs within years since data spans multiple years
DT[, yr := year(datetime)]

## Daily type randomizations
# splitBy = 'yr' to force randomization only within year
randDaily <- randomizations(
   DT,
   type = 'daily',
   id = 'ID',
   group = 'group',
   coords = NULL,
   datetime = 'datetime',
   splitBy = 'yr',
   iterations = 20
)
```

## 5. c) `type = 'trajectory'`
`'trajectory'` randomizes daily trajectories within individuals ([Spiegel et al. 2016](https://besjournals.onlinelibrary.wiley.com/doi/abs/10.1111/2041-210X.12553)). The `datetime` argument expects a datetime `POSIXct` format column. 

Five columns are returned when `type = 'trajectory'` along with `id`, `datetime` and `splitBy` columns:

* random date time ("random" prefixed to *datetime* argument)
* 'jul' - observed julian day
* 'observed' - observed rows (TRUE/FALSE)
* 'iteration' - which iteration rows correspond to (0 is observed)
* 'randomJul' - random julian day relocations are swapped to from observed julian day

```{r}
# Calculate year column to ensure randomization only occurs within years since data spans multiple years
DT[, yr := year(datetime)]

## Trajectory type randomization
randTraj <- randomizations(
   DT,
   type = 'trajectory',
   id = 'ID',
   group = NULL,
   coords = c('X', 'Y'),
   datetime = 'datetime',
   splitBy = 'yr',
   iterations = 20
)
```


# Build random network 
Once we've randomized the data stream with `randomizations`, we can build the random network. 

We will use the `get_gbi` function directly when `type` is either 'step' or 'daily'. For `type = 'trajectory'`, we will recalculate spatial groups with one of `group_pts`, `group_lines`, `group_polys` for the randomized data. In this case, the example shows `group_pts`. 

Since we want to create a group by individual matrix for each random iteration (and in this case, each year), we will use `mapply` to work on subsets of the randomized data. 

Note: building the random networks depends on the `type` used and therefore, the following chunks are mutually exclusive. Use the one that corresponds to the randomization type you used above. 

## 6. a) `type = 'step'`
`randomizations` with `type = 'step'` returns a 'randomID' which should be used instead of the observed 'ID' to generate the group by indiviual matrix. 

After `get_gbi`, we use `asnipe::get_network` to build the random network. 

```{r}
## Create a data.table of unique combinations of iteration and year, exluding observed rows
iterYearLs <- unique(randStep[!(observed), .(iteration, yr)])

## Generate group by individual matrix 
# for each combination of iteration number and year
# 'group' generated by spatsoc::group_pts
# 'randomID' used instead of observed ID (type = 'step')
gbiLs <- mapply(FUN = function(i, y) {
  get_gbi(randStep[iteration == i & yr == y],
          'group', 'randomID')
  },
  i = iterYearLs$iter,
  y = iterYearLs$yr,
  SIMPLIFY = FALSE
)

## Generate a list of random networks
netLs <- lapply(gbiLs, FUN = get_network,
                data_format = "GBI", association_index = "SRI")

```


## 6. b) `type = 'daily'`
`randomizations` with `type = 'step'` returns a 'randomID' which should be used instead of the observed 'ID' to generate the group by indiviual matrix. 

After `get_gbi`, we use `asnipe::get_network` to build the random network. 

In this case, we will generate a fake column representing a "population" to show how we can translate the `mapply` chunk above to three (or more variables). 

```{r}
## Generate fake population
randDaily[, population := sample(1:2, .N, replace = TRUE)]

## Create a data.table of unique combinations of iteration, year, and population, exluding observed rows
iterYearLs <- unique(randStep[!(observed), .(iteration, yr, population)])

## Generate group by individual matrix 
# for each combination of iteration number and year
# 'group' generated by spatsoc::group_pts
# 'randomID' used instead of observed ID (type = 'step')
gbiLs <- mapply(FUN = function(i, y, p) {
  get_gbi(randDaily[iteration == i & 
                      yr == y & population == p],
          'group', 'randomID')
  },
  i = iterYearLs$iter,
  y = iterYearLs$yr,
  p = iterYearLs$population,
  SIMPLIFY = FALSE
)

## Generate a list of random networks
netLs <- lapply(gbiLs, FUN = get_network,
                data_format = "GBI", association_index = "SRI")

```



## 6. c) `type = 'trajectory'`
`randomizations` with `type = 'trajectory'` returns a random date time which should be used instead of the observed date time to generate random gambit of the group data. 

First, we pass the randomized data to `group_times` using the random date time for `datetime`. 

After `get_gbi`, we use `asnipe::get_network` to build the random network. 

```{r}
## Randomized temporal groups
# 'datetime' is the randomdatetime produced by randomizations(type = 'trajectory')
group_times(randTraj, datetime = 'randomdatetime', threshold = '5 minutes')

## Randomized spatial groups
# 'iteration' used in splitBy to ensure only points within each iteration are grouped
group_pts(randTraj, threshold = 50, id = 'ID', coords = c('X', 'Y'),
          timegroup = 'timegroup', splitBy = 'iteration')

## Create a data.table of unique combinations of iteration and year, exluding observed rows
iterYearLs <- unique(randStep[!(observed), .(iteration, yr)])

## Generate group by individual matrix 
# for each combination of iteration number and year
# 'group' generated by spatsoc::group_pts
# 'ID' used since datetimes were randomized within individuals
gbiLs <- mapply(FUN = function(i, y) {
  get_gbi(randTraj[iteration == i & yr == y],
          'group', 'ID')
  },
  i = iterYearLs$iter,
  y = iterYearLs$yr,
  SIMPLIFY = FALSE
)

## Generate a list of random networks
netLs <- lapply(gbiLs, FUN = get_network,
                data_format = "GBI", association_index = "SRI")

```


# Network metrics
Finally, we can calculate some network metrics. Please note that there are many ways of interpreting, analyzing and measuring networks, so this will simply show one option. 


## 7. Calculate observed network metrics
To calculate observed network metrics, use the network (`net`) produced in [4.](#asnipeget_network) from 2016 data.

```{r}
## Generate graph
g <- graph.adjacency(net, 'undirected', 
                     diag = FALSE, weighted = TRUE)

## Metrics for all individuals 
observed <- data.table(
  centrality = evcent(g)$vector,
  strength = graph.strength(g),
  degree = degree(g),
  ID = names(degree(g)),
  yr = subDT[, unique(yr)]
)
```


## 8. Calculate random network metrics
With the list of random networks from [6.](#build-random-network), we can generate a list of graphs with `igraph::graph.adjacency` (for example) and calculate random network metrics. 

This example uses the `netLs` generated by [6. a)](#a-type-step-1) which was split by year and iteration. 

```{r}
## Generate graph and calculate network metrics
mets <- lapply(seq_along(netLs), function(n) {
  g <- graph.adjacency(netLs[[n]], 'undirected', 
                       diag = FALSE, weighted = TRUE)
  
  data.table(
    centrality = evcent(g)$vector,
    strength = graph.strength(g),
    degree = degree(g),
    ID = names(degree(g)),
    iteration = iterYearLs$iter[[n]],
    yr = iterYearLs$yr[[n]]
    )
})

## Metrics for all individuals across all iterations and years
random <- rbindlist(mets)

## Mean values for each individual and year
meanMets <- random[, lapply(.SD, mean), by = .(ID, yr),
                .SDcols = c('centrality', 'strength', 'degree')]
```

## 9. Compare observed and random metrics
Instead of calculating observed and random metrics separately (shown in [7.](#calculate-observed-network-metrics) and [8.](#calculate-random-network-metrics)), we can calculate metrics for both at the same time and compare. 

This chunk expects the outputs from [5. a)](#a-type-step), skipping steps 6.-8.

Note: by removing the `!(observed)` subset from `randStep` performed in [6. a)](#a-type-step-1), we will include observed rows where `iteration == 0`. This will return a `gbiLs` where the observed and random rows are included in the same `data.table`. 

```{r}
## Create a data.table of unique combinations of iteration and year, including observed and random rows
iterYearLs <- unique(randStep[, .(iteration, yr)])

## Generate group by individual matrix 
# for each combination of iteration and year
# 'group' generated by spatsoc::group_pts
# 'randomID' used instead of observed ID (type = 'step')
gbiLs <- mapply(FUN = function(i, y) {
  get_gbi(randStep[iteration == i & yr == y],
          'group', 'randomID')
  },
  i = iterYearLs$iter,
  y = iterYearLs$yr,
  SIMPLIFY = FALSE
)

## Generate a list of random networks
netLs <- lapply(gbiLs, FUN = get_network,
                data_format = "GBI", association_index = "SRI")

## Generate graph and calculate network metrics
mets <- lapply(seq_along(netLs), function(n) {
  g <- graph.adjacency(netLs[[n]], 'undirected', 
                       diag = FALSE, weighted = TRUE)
  
  data.table(
    centrality = evcent(g)$vector,
    strength = graph.strength(g),
    ID = names(degree(g)),
    iteration = iterYearLs$iter[[n]],
    yr = iterYearLs$yr[[n]]
    )
})

## Observed and random for all individuals across all iterations and years
out <- rbindlist(mets)

## Split observed and random
out[, observed := ifelse(iteration == 0, TRUE, FALSE)]

## Mean values for each individual and year, by observed/random
meanMets <- out[, lapply(.SD, mean), by = .(ID, yr, observed),
                .SDcols = c('centrality', 'strength')]

```
---
title: "Introduction to spatsoc"
author: "Alec Robitaille, Quinn Webber and Eric Vander Wal"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to spatsoc}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---

```{r knitropts, include = FALSE}
knitr::opts_chunk$set(message = TRUE, 
                      warning = FALSE,
                      eval = FALSE, 
                      echo = FALSE)
```

The `spatsoc` package provides functionality for analyzing animal relocation data in time and space to identify potential interactions among individuals and build gambit-of-the-group data for constructing social networks. 

The package contains grouping and edge list generating functions that are used for identifying spatially and temporally explicit groups from input data. In addition, we provide social network analysis functions for randomizing individual identifiers within groups, designed to test whether social networks generated from animal relocation data were based on non-random social proximity among individuals and for generating group by individual matrices.

The functions were developed for application across animal relocation data, for example, proximity based social network analyses and spatial and temporal clustering of points.

## Data preparation

```{r libs, eval = TRUE, include = FALSE}
library(spatsoc); library(data.table)
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]
```

`spatsoc` expects a `data.table` for all of its functions. If you have a `data.frame`, you can use `data.table::setDT()` to convert it by reference. If your data is a CSV, you can use `data.table::fread()` to import it as a `data.table`. 

The data consist of relocations of `r DT[, uniqueN(ID)]` individuals over `r DT[, max(yday(datetime)) - min(yday(datetime))]` days. Using these data, we can compare the various grouping methods available in `spatsoc`. Note: these examples will use a subset of the data, only individuals H, I and J. 


```{r, echo = TRUE}
# Load packages
library(spatsoc)
library(data.table)

# Read in spatsoc's example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Use subset of individuals
DT <- DT[ID %in% c('H', 'I', 'J')]

# Cast character column 'datetime' as POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]
```

```{r, eval = TRUE}
DT <- DT[ID %chin% c('H', 'I', 'J')]
knitr::kable(DT[, .SD[1:2], ID][order(ID)])
```

## Temporal grouping
The `group_times` function is used to group relocations temporally. It is flexible to a threshold provided in units of minutes, hours or days. Since GPS fixes taken at regular intervals have some level of variability, we will provide a time threshold (`threshold`), to consider all fixes within this threshold taken at the same time. Alternatively, we may want to understand different scales of grouping, perhaps daily movement trajectories or seasonal home range overlap. 

```{r groupmins, echo = TRUE}
group_times(DT, datetime = 'datetime', threshold = '5 minutes')
```

```{r tableSetUp, eval = TRUE}
nRows <- 9
```

```{r tabgroupmins, eval = TRUE}
knitr::kable(
  group_times(DT, threshold = '5 minutes', datetime = 'datetime')[, 
    .(ID, X, Y, datetime, minutes, timegroup)][
      order(datetime)][
        1:nRows])
```


A message is returned when `group_times` is run again on the same `DT`, as the columns already exist in the input `DT` and will be overwritten. 

```{r grouphours, echo = TRUE}
group_times(DT, datetime = 'datetime', threshold = '2 hours')
```

```{r tabgrouphours, eval = TRUE}
knitr::kable(
  group_times(DT, threshold = '2 hours', datetime = 'datetime')[, 
    .(ID, X, Y, datetime, hours, timegroup)][
      order(datetime)][
        1:nRows])
```

```{r groupdays, echo = TRUE}
group_times(DT, datetime = 'datetime', threshold = '5 days')
```

```{r tabgroupdays, eval = TRUE}
knitr::kable(
  group_times(DT, threshold = '5 days', datetime = 'datetime')[, .SD[sample(.N, 3)], by = .(timegroup, block)][order(datetime)][
        1:nRows, .(ID, X, Y, datetime, block, timegroup)])
```

## Spatial grouping
The `group_pts` function compares the relocations of all individuals in each timegroup and groups individuals based on a distance threshold provided by the user. The `group_pts` function uses the "chain rule" where three or more individuals that are all within the defined threshold distance of at least one other individual are considered in the same group. For point based spatial grouping with a distance threshold that does not use the chain rule, see `edge_dist` below.

```{r grouppts, echo = TRUE}
group_times(DT = DT, datetime = 'datetime', threshold = '15 minutes')
group_pts(DT, threshold = 50, id = 'ID', coords = c('X', 'Y'), timegroup = 'timegroup')
```

```{r fakegrouppts, eval = TRUE}
DT <- group_times(DT = DT, datetime = 'datetime', 
                     threshold = '15 minutes')
DT <- group_pts(
    DT = DT,
    threshold = 50, id = 'ID', coords = c('X', 'Y'),
    timegroup = 'timegroup')

knitr::kable(
  DT[
      between(group,  771, 774)][order(timegroup)][
        1:nRows, .(ID, X, Y, timegroup, group)]
)
```


The `group_lines` function groups individuals whose trajectories intersect in a specified time interval. This represents a coarser grouping method than `group_pts` which can help understand shared space at daily, weekly or other temporal resolutions.

```{r fakegrouplines, echo = TRUE}
utm <- 'EPSG:32736'

group_times(DT = DT, datetime = 'datetime', threshold = '1 day')
group_lines(DT, threshold = 50, projection = utm, 
            id = 'ID', coords = c('X', 'Y'),
            timegroup = 'timegroup', sortBy = 'datetime')
```

```{r grouplines, eval = TRUE}
utm <- 'EPSG:32736'

DT <- group_times(DT = DT, datetime = 'datetime', 
                threshold = '1 day')
DT <- group_lines(DT,
                  threshold = 50, projection = utm, 
                  id = 'ID', coords = c('X', 'Y'), 
                  sortBy = 'datetime', timegroup = 'timegroup')
knitr::kable(
  unique(DT[, .(ID, timegroup, group)])[, .SD[1:3], ID][order(timegroup)]
)
```


The `group_polys` function groups individuals whose home ranges intersect. This represents the coarsest grouping method, to provide a measure of overlap across seasons, years or all available relocations. It can either return the proportion of home range area overlapping between individuals or simple groups. Home ranges are calculated using `adehabitatHR::kernelUD` or `adehabitatHR::mcp`. Alternatively, a `SpatialPolygonsDataFrame` can be input to the `spPolys` argument.

```{r fakegrouppolys, echo = TRUE}
utm <- 'EPSG:32736'
group_times(DT = DT, datetime = 'datetime', threshold = '8 days')
group_polys(DT = DT, area = TRUE, hrType = 'mcp',
           hrParams = list('percent' = 95),
           projection = utm,
           coords = c('X', 'Y'), id = 'ID')
```

```{r grouppolys, eval = TRUE}
utm <- 'EPSG:32736'
DT <- group_times(DT = DT, datetime = 'datetime', 
                threshold = '8 days')
knitr::kable(
  data.frame(group_polys(
    DT, 
    area = TRUE, hrType = 'mcp',
           hrParams = list('percent' = 95),
           projection = utm,
           coords = c('X', 'Y'), id = 'ID')[
             , .(ID1, ID2, area, proportion)])
)
```

## Edge list generation
The `edge_dist` function calculates the geographic distance between between individuals within each timegroup and returns all paired relocations within the spatial threshold. `edge_dist` uses a distance matrix like group_pts, but, in contrast, does not use the chain rule to group relocations.

```{r edgedist, echo = TRUE}
group_times(DT = DT, datetime = 'datetime', threshold = '15 minutes')
edge_dist(DT, threshold = 50, id = 'ID', coords = c('X', 'Y'), timegroup = 'timegroup', fillNA = TRUE)
```

```{r fakeedgedist, eval = TRUE}
DT <- group_times(DT = DT, datetime = 'datetime', 
                     threshold = '15 minutes')
edges <- edge_dist(DT, threshold = 50, id = 'ID', coords = c('X', 'Y'), timegroup = 'timegroup', fillNA = TRUE)

knitr::kable(
  edges[between(timegroup, 158, 160)]
)
```

The `edge_nn` function calculates the nearest neighbour to each individual within each time group. If the optional distance threshold is provided, it is used to limit the maximum distance between neighbours. `edge_nn` returns an edge list of each individual and their nearest neighbour.

```{r edgenn, echo = TRUE}
group_times(DT = DT, datetime = 'datetime', threshold = '15 minutes')
edge_nn(DT, id = 'ID', coords = c('X', 'Y'), timegroup = 'timegroup')
```

```{r fakeedgenn, eval = TRUE}
DT <- group_times(DT = DT, datetime = 'datetime', 
                     threshold = '15 minutes')
edges <- edge_nn(DT, id = 'ID', coords = c('X', 'Y'), timegroup = 'timegroup')

knitr::kable(
  edges[1:6]
)
```

## Notes
Package dependencies for `spatsoc` are `sp`, `rgeos`, `igraph`, `adehabitatHR` and `data.table`. `data.table` provides efficient methods for manipulating large (or small) datasets. As a result, input `DT` for all `spatsoc` functions must be a `data.table` and if it isn't, you can simply use `data.table::setDT(df)` to convert it by reference. 

In addition, since the `rgeos` package is used in most functions (`group_lines` and `group_polys`) the input `DT`'s coordinate system is important. `rgeos` expects planar coordinates and this requirement is carried forward for `spatsoc`. Since `rgeos` is used, system dependencies include `GEOS`. 
---
title: "Frequently asked questions about spatsoc"
author: "Alec Robitaille, Quinn Webber and Eric Vander Wal"
date: "`r Sys.Date()`"
output: 
  rmarkdown::html_vignette:
    number_sections: true
    toc: true
vignette: >
  %\VignetteIndexEntry{FAQ}
  %\VignetteEngine{knitr::knitr}
  %\VignetteEncoding{UTF-8}
---

```{r knitropts, include = FALSE}
knitr::opts_chunk$set(message = FALSE, 
                      warning = FALSE,
                      eval = FALSE, 
                      echo = TRUE)
```

spatsoc is an R package for detecting spatial and temporal groups in GPS relocations. It can be used to build proximity-based social networks using gambit-of-the-group format and edge-lists. In addition, the randomization function provides data-stream randomization methods suitable for GPS data.

# Usage
`spatsoc` leverages `data.table` to modify by reference and iteratively work on subsets of the input data. The first input for all functions in `spatsoc` is `DT`, an input `data.table`. If your data is a `data.frame`, you can convert it by reference using `setDT(DF)`. 

## Spatial and temporal grouping

`spatsoc` is designed to work in two steps: temporal followed by either spatial grouping or edge list generating. Considering your specific study species and system, determine a relevant temporal and spatial grouping threshold. This may be 5 minutes and 50 meters or 2 days and 100 meters or any other thresholds - the functions provided by `spatsoc` are flexible to user input. In some cases, the spatial grouping function selected is only relevant with certain temporal grouping thresholds. For example, we wouldn't expect a threshold of 5 minutes with `group_polys`. 

```{r, eval = TRUE, results = 'hide'}
# Load packages
library(spatsoc)
library(data.table)

# Read data as a data.table
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast datetime column to POSIXct
DT[, datetime := as.POSIXct(datetime)]

# Temporal groups
group_times(DT, datetime = 'datetime', threshold = '5 minutes')

# Spatial groups
group_pts(
  DT,
  threshold = 50,
  id = 'ID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup'
)
```

```{r, eval = TRUE, echo = FALSE}
knitr::kable(DT[order(group, timegroup)][1:5, .(ID, X, Y, datetime, timegroup, group)])
```


## Social network analysis
See the vignette about [using spatsoc in social network analysis](https://docs.ropensci.org/spatsoc/articles/using-in-sna.html).


# Installation

## System dependencies

### GEOS
Install `GEOS`:

- Debian/Ubuntu: `apt-get install libgeos-dev`
- Arch: `pacman -S geos`
- Fedora: `dnf install geos geos-devel`
- Mac: `brew install geos`
- Windows: see [here](https://trac.osgeo.org/osgeo4w/)

## Package dependencies

* `data.table`
* `igraph`
* `sp`
* `adehabitatHR`
* `rgeos`



# Functions
## group_times

`group_times(DT, datetime, threshold)`

* `DT`: input `data.table`
* `datetime`: date time column name in input data.table
* `threshold`: threshold for grouping 

### DT
A `data.table` with a date time formatted column. The input `DT` will be returned with columns appended. The `timegroup` column corresponds to the temporal group assigned to each row. Please note that the actual value of the time group is meaningless. Reordered data will return a different time group. What is meaningful, however,  is the contents of each group. Each group will contain all rows nearest to the threshold provided. 


### datetime format
The `group_times` function expects either one column (`POSIXct`) or two columns (`IDate` and `ITime`). 

Given a character column representing the date time, convert it to `POSIXct` or `IDate` and `ITime`:

```{r posixct}
DT[, datetime := as.POSIXct(datetime)]
DT[, c('idate', 'itime') := IDateTime(datetime)]
```


These are then provided to the function using the names of the column in the input data. 

`group_times(DT, datetime = 'datetime', threshold = '5 minutes')`

or

`group_times(DT, datetime = c('idate', 'itime'), threshold = '5 minutes')`


### threshold recommendations
The `threshold` provided to `group_times` should be related to the fix rate of the input dataset and to the specific study system and species. If relocations are recorded every two hours, a `threshold = '2 hours'` will group all rows to the nearest two hour group (10am, 12pm, 2pm, 4pm, ...). This, however, means that the relocations can be up to one hour apart from each other. Picking a smaller threshold, e.g.: `threshold = '15 minutes'` may be more relevant in some cases. The flexibility of `spatsoc`'s threshold argument means the user must carefully consider what threshold is reasonable to their specific system. 

### Limitations of threshold
The `threshold` of `group_times` is considered only within the scope of 24 hours and this poses limitations on it:

* `threshold` must evenly divide into 60 minutes or 24 hours
* multi-day blocks are consistent **across years** and timegroups from these are **by year**. 
* number of minutes cannot exceed 60
* `threshold` cannot be fractional

### Columns returned by group_times
The main column returned by `group_times` is "timegroup". It represents the temporal group of each row, where those nearest (either above or below) within the threshold are grouped. Its actual value does not have any meaning, but the contents of each group do. That means if the data is reordered, a row may have a different time group, but the other rows in that group should not change.  

The extra columns are provided to help the user investigate, troubleshoot and interpret the timegroup. 

| threshold unit | column(s) added |
|----------------|-----------------|
| minute         | "minutes" column added identifying the nearest minute group for each row.      |
| hour           | "hours" column added identifying the nearest hour group for each row.         |
| day            | "block" columns added identifying the multiday block for each row. |


### Warnings and messages

* "columns found in input DT and will be overwritten by this function"

This message is returned to the user when a column matching those returned by `group_times` is found in the input DT. This is commonly the case when `group_times` is run multiple times consecutively. 


* "no threshold provided, using the time field directly to group"

This message is returned to the user when the `threshold` is NULL. This is the default setting of `threshold` and, at times, may be suitable. In this case, the date times in the `datetime` column will be grouped exactly. Usually, a threshold should be provided. 

* "the minimum and maximum days in DT are not evenly divisible by the provided block length"

This warning is returned to the user when the `threshold` with unit days does not divide evenly into the range of days in DT. For example, if DT had data covering 30 days, and a threshold of '7 days' was used, this warning would be returned. Note, this warning is returned for the range of days for the entire data set and not by year. 

## group_pts


`group_pts(DT, threshold, id, coords, timegroup, splitBy)`
  
* `DT`: input `data.table`
* `threshold`: threshold for grouping 
* `id`: column name of IDs in `DT`
* `coords`: column names of x and y coordinates in `DT`
* `timegroup`: (optional) column name of time group
* `splitBy`: (optional) column names of extra variables to group on

### <a name="group_pts DT"></a>DT
The input `data.table`. It will returned with a column named group appended, which represents the spatial (and temporal if `timegroup` is provided) group. 

### threshold
The threshold must be in the units of the coordinates. 

### coords
The coordinates must be planar, such as UTM (of whichever zone your relocations are in). 

## group_lines
`group_lines(DT, threshold, projection, id, coords, timegroup, sortBy, splitBy, spLines)`
  
* `DT`: input `data.table`
* `threshold`: threshold for grouping 
* `projection`: projection of coordinates in `DT`
* `id`: column name of IDs in `DT`
* `coords`: column names of x and y coordinates in `DT`
* `timegroup`: (optional) column name of time group
* `sortBy`: column name of date time to sort rows for building lines
* `splitBy`: (optional) column names of extra variables to group on
* `spLines`: alternatively, provide solely a `SpatialLines` object

### DT
See [3.2.1](#dt-1). 

### threshold
The `threshold` argument represents a buffer area around each line. When `threshold = 0`, the lines are grouped by spatial overlap. If the threshold is greater than 0, the lines buffered, then grouped by spatial overlap. 

### projection
The `projection` argument expects a character string defining the projection to
be passed to `sp::CRS`. For example, for UTM zone 36S (EPSG 32736), the projection
argument is `projection = "EPSG:32736"`. See \url{https://spatialreference.org} 
for a list of EPSG codes. Please note, R spatial has followed updates to GDAL 
and PROJ for handling projections, see more at
\url{https://www.r-spatial.org/r/2020/03/17/wkt.html}. 

Because of changes to projection handling, `group_polys`, `build_polys`, `group_lines` and `build_lines` may return the warning once, or many times: 

```
In proj4string(xy) : CRS object has comment, which is lost in output
Calls: group_polys ... do.call -> <Anonymous> -> proj4string -> proj4string
```

This warning is explicitly verbose, to ensure we are considering the updated use of CRS [in `sp`](https://www.r-spatial.org/r/2020/03/17/wkt.html#crs-objects-in-sp) and other spatial packages.  


### sortBy
The `sortBy` argument expects a date time formatted column name, which is used to order the rows for each individual (and `splitBy`). 

## group_polys
`group_polys(DT, area, hrType, hrParams, projection, id, coords, timegroup, splitBy, spLines)`
  
* `DT`: input `data.table`
* `area`: boolean argument if proportional area should be returned
* `hrType`: type of home range created 
* `hrParams`: parameters relevant to the type of home range created
* `projection`: projection of coordinates in `DT`
* `id`: column name of IDs in `DT`
* `coords`: column names of x and y coordinates in `DT`
* `timegroup`: (optional) column name of time group
* `splitBy`: (optional) column names of extra variables to group on
* `spPolys`: alternatively, provide solely a `SpatialPolygons` object


### DT and area
If `area = FALSE`, see [3.2.1](#dt-1). If `area = TRUE`, the DT will not be appended with a group column instead a `data.table` with IDs and proportional area overlap will be returned. 

The default unit for area overlap is square meters. 

<!-- direction of proportion -->

### projection
The `projection` argument expects a character string defining the projection to be
passed to `sp::CRS`. For example, for UTM zone 36S (EPSG 32736), the projection
argument is `projection = "EPSG:32736"`. See \url{https://spatialreference.org} 
for a list of EPSG codes. Please note, R spatial has followed updates to GDAL 
and PROJ for handling projections, see more at
\url{https://www.r-spatial.org/r/2020/03/17/wkt.html}.

Because of changes to projection handling, `group_polys`, `build_polys`, `group_lines` and `build_lines` may return the warning once, or many times: 

```
In proj4string(xy) : CRS object has comment, which is lost in output
Calls: group_polys ... do.call -> <Anonymous> -> proj4string -> proj4string
```

This warning is explicitly verbose, to ensure we are considering the updated use of CRS [in `sp`](https://www.r-spatial.org/r/2020/03/17/wkt.html#crs-objects-in-sp) and other spatial packages.  


### hrType and hrParams
Currently, `spatsoc` offers two types of home ranges provided by the `adehabitatHR` package: 'mcp' (`mcp`) and 'kernel' (`kernelUD` and `getverticeshr`). The parameters must match the arguments of those functions. 

Internally, we match arguments to the functions allowing the user to provide, for example, both the percent (provided to `getverticeshr`) and grid arguments (provided to `mcp`).

```{r}
group_polys(
  DT,
  area = FALSE,
  projection = utm,
  hrType = 'mcp',
  hrParams = list(grid = 60, percent = 95),
  id = 'ID',
  coords = c('X', 'Y')
)
```

## edge_dist
`edge_dist(DT = NULL, threshold = NULL, id = NULL, coords = NULL, timegroup = NULL, splitBy = NULL, fillNA = TRUE)`

* `DT`: input `data.table`
* `threshold`: threshold for grouping 
* `id`: column name of IDs in `DT`
* `coords`: column names of x and y coordinates in `DT`
* `timegroup`: (optional) column name of time group
* `splitBy`: (optional) column names of extra variables to group on
* `fillNA`: boolean indicating if NAs should be returned for individuals that were not within the threshold distance of any other. If TRUE, NAs are returned. If FALSE, only edges between individuals within the threshold distance are returned.

This is the non-chain rule implementation similar to `group_pts`. Edges are defined by the distance threshold and NAs are returned for individuals within each timegroup if they are not within the threshold distance of any other individual (if `fillNA` is TRUE). 

**See the vignette [Using edge list generating functions and dyad_id](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html) for details about the `edge_dist` function.**

## edge_nn
`edge_nn(DT = NULL, id = NULL, coords = NULL, timegroup = NULL, splitBy = NULL, threshold = NULL)`

* `DT`: input `data.table`
* `id`: column name of IDs in `DT`
* `coords`: column names of x and y coordinates in `DT`
* `timegroup`: (optional) column name of time group
* `splitBy`: (optional) column names of extra variables to group on
* `threshold`: (optional) spatial distance threshold to set maximum distance between an individual and their neighbour.

This function can be used to generate edge lists defined either by nearest neighbour or nearest neighbour with a maximum distance. NAs are returned for nearest neighbour for an individual was alone in a timegroup (and/or splitBy) or if the distance between an individual and it's nearest neighbour is greater than the threshold. 


**See the vignette [Using edge list generating functions and dyad_id](https://docs.ropensci.org/spatsoc/articles/using-edge-and-dyad.html) for details about the `edge_nn` function.**

## randomizations
`randomizations(DT, type, id, datetime, splitBy, iterations)`
  
* `DT`: input `data.table`
* `type`: one of 'daily', 'step' or 'trajectory' 
* `id`: Character string of ID column name
* `datetime`: field used for providing date time or time group - see details
* `splitBy`: List of fields in DT to split the randomization process by
* `iterations`: The number of iterations to randomize


**See the vignette [Using spatsoc in social network analysis](https://docs.ropensci.org/spatsoc/articles/using-in-sna.html) for details about the `randomizations` function (specifically the section 'Data stream randomization')**


# Package design
## Don't I need to reassign to save the output?

(Almost) all functions in `spatsoc` use data.table's modify-by-reference to reduce recopying large datasets and improve performance. The exceptions are `group_polys(area = TRUE)`, `randomizations` and the edge list generating functions `edge_dist` and `edge_nn`.


## Why does a function print the result, but columns aren't added to my DT?

Check that your `data.table` has columns allocated (with `data.table::truelength`) and if not, use `data.table::setDT` or `data.table::alloc.col`. This can happen if you are reading your data from `RDS` or `RData` files.  [See here.](https://cran.r-project.org/package=data.table/vignettes/datatable-faq.html#reading-data.table-from-rds-or-rdata-file)

```{r setdt}
if (truelength(DT) == 0) {
  setDT(DT)
}
# then go to spatsoc
group_times(DT, datetime = 'datetime', threshold = '5 minutes')
```

or simply:

```{r alloc}
DT <- readRDS('path/to/data.Rds')
alloc.col(DT)
```




# Summary information
Here are some useful code chunks for understanding the spatial and temporal extent of your data and the outputs of `spatsoc` functions. 

## Number of individuals
```{r}
# Number of unique individuals
DT[, uniqueN(ID)]

# Number of unique individuals by timegroup
DT[, uniqueN(ID), by = timegroup]
```

## Temporal range

```{r}
# Min, max datetime
DT[, range(datetime)]

# Difference between relocations in hours
DT[order(datetime), 
   .(difHours = as.numeric(difftime(datetime, shift(datetime), units = 'hours'))), 
   by = ID]

# Difference between relocations in hours
DT[order(datetime), 
   .(difMins = as.numeric(difftime(datetime, shift(datetime), units = 'mins'))), 
   by = ID]
```

## Spatial extent
Simple spatial extents can be calculated for all individuals or by individual. 
```{r}
# All individuals
DT[, .(minX = min(X),
       maxX = max(X),
       minY = min(Y),
       maxY = max(Y),)]

# By individual
DT[, .(minX = min(X),
       maxX = max(X),
       minY = min(Y),
       maxY = max(Y),),
   by = ID]
```

## `spatsoc` outputs
After using the grouping functions, we can determine the number of individuals in a temporal or spatial group. 
```{r}
# Number of unique individuals by timegroup
DT[, uniqueN(ID), by = timegroup]

# Number of unique individuals by group
DT[, uniqueN(ID), by = group]
```
---
title: "Using edge list generating functions and dyad_id"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using edge list generating functions and dyad_id}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = FALSE,
  echo = TRUE,
  comment = "#>"
)
```



`spatsoc` can be used in social network analysis to generate edge lists from GPS relocation data. 


Edge lists are generated using either the `edge_dist` or the `edge_nn` function. 


**Note**: The grouping functions and their application in social network analysis are further described in the vignette [Using spatsoc in social network analysis - grouping functions](https://docs.ropensci.org/spatsoc/articles/using-in-sna.html). 


## Generate edge lists
spatsoc provides users with one temporal (`group_times`) and two edge list generating functions (`edge_dist`, `edge_nn`) to generate edge lists from GPS relocations. Users can consider edges defined by either the spatial proximity between individuals (with `edge_dist`), by nearest neighbour (with `edge_nn`) or by nearest neighbour with a maximum distance (with `edge_nn`). The edge lists can be used directly by the animal social network package `asnipe` to generate networks. 

### 1. Load packages and prepare data
`spatsoc` expects a `data.table` for all `DT` arguments and date time columns to be formatted `POSIXct`. 

```{r, message = FALSE, warning = FALSE, eval = TRUE}
## Load packages
library(spatsoc)
library(data.table)

## Read data as a data.table
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

## Cast datetime column to POSIXct
DT[, datetime := as.POSIXct(datetime)]
```


Next, we will group relocations temporally with `group_times` and generate edges lists with one of `edge_dist`, `edge_dist`. Note: these are mutually exclusive, only select one edge list generating function at a time. 

### 2. a) `edge_dist` 

Distance based edge lists where relocations in each timegroup are considered edges if they are within the spatial distance defined by the user with the `threshold` argument. Depending on species and study system, relevant temporal and spatial distance thresholds are used. In this case, relocations within 5 minutes and 50 meters are considered edges. 

This is the non-chain rule implementation similar to `group_pts`. Edges are defined by the distance threshold and NAs are returned for individuals within each timegroup if they are not within the threshold distance of any other individual (if `fillNA` is TRUE). 

Optionally, `edge_dist` can return the distances between individuals (less than the threshold) in a column named 'distance' with argument `returnDist = TRUE`. 

```{r, eval = TRUE}
# Temporal groups
group_times(DT, datetime = 'datetime', threshold = '5 minutes')

# Edge list generation
edges <- edge_dist(
  DT,
  threshold = 100,
  id = 'ID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup',
  returnDist = TRUE,
  fillNA = TRUE
)
```

### 2. b) `edge_nn`

Nearest neighbour based edge lists where each individual is connected to their nearest neighbour. `edge_nn` can be used to generate edge lists defined either by nearest neighbour or nearest neighbour with a maximum distance. As with grouping functions and `edge_dist`, temporal and spatial threshold depend on  species and study system. 

NAs are returned for nearest neighbour for an individual was alone in a timegroup (and/or splitBy) or if the distance between an individual and its nearest neighbour is greater than the threshold. 

Optionally, `edge_nn` can return the distances between individuals (less than the threshold) in a column named 'distance' with argument `returnDist = TRUE`. 

```{r, eval = FALSE}
# Temporal groups
group_times(DT, datetime = 'datetime', threshold = '5 minutes')

# Edge list generation
edges <- edge_nn(
  DT,
  id = 'ID',
  coords = c('X', 'Y'),
  timegroup = 'timegroup'
)

# Edge list generation using maximum distance threshold
edges <- edge_nn(
  DT, 
  id = 'ID', 
  coords = c('X', 'Y'),
  timegroup = 'timegroup', 
  threshold = 100
)

# Edge list generation using maximum distance threshold, returning distances
edges <- edge_nn(
  DT, 
  id = 'ID', 
  coords = c('X', 'Y'),
  timegroup = 'timegroup', 
  threshold = 100,
  returnDist = TRUE
)

```


## Dyads

### 3. `dyad_id`

The function `dyad_id` can be used to generate a unique, undirected dyad identifier for edge lists. 

```{r, eval = TRUE}
# In this case, using the edges generated in 2. a) edge_dist
dyad_id(edges, id1 = 'ID1', id2 = 'ID2')
```


Once we have generated dyad ids, we can measure consecutive relocations, start and end relocation, etc. **Note:** since the edges are duplicated A-B and B-A, you will need to use the unique timegroup*dyadID or divide counts by 2. 


### 4. Dyad stats

```{r, eval = TRUE}
# Get the unique dyads by timegroup
# NOTE: we are explicitly selecting only where dyadID is not NA
dyads <- unique(edges[!is.na(dyadID)], by = c('timegroup', 'dyadID'))

# NOTE: if we wanted to also include where dyadID is NA, we should do it explicitly
# dyadNN <- unique(DT[!is.na(NN)], by = c('timegroup', 'dyadID'))

# Get where NN was NA
# dyadNA <- DT[is.na(NN)]

# Combine where NN is NA
# dyads <- rbindlist(list(dyadNN, dyadNA))


# Set the order of the rows
setorder(dyads, timegroup)

## Count number of timegroups dyads are observed together
dyads[, nObs := .N, by = .(dyadID)]

## Count consecutive relocations together
# Shift the timegroup within dyadIDs
dyads[, shifttimegrp := shift(timegroup, 1), by =  dyadID]

# Difference between consecutive timegroups for each dyadID
# where difftimegrp == 1, the dyads remained together in consecutive timegroups
dyads[, difftimegrp := timegroup - shifttimegrp]


# Run id of diff timegroups
dyads[, runid := rleid(difftimegrp), by = dyadID]

# N consecutive observations of dyadIDs
dyads[, runCount := fifelse(difftimegrp == 1, .N, NA_integer_), by = .(runid, dyadID)]

## Start and end of consecutive relocations for each dyad
# Dont consider where runs aren't more than one relocation
dyads[runCount > 1, start := fifelse(timegroup == min(timegroup), TRUE, FALSE), by = .(runid, dyadID)]

dyads[runCount > 1, end := fifelse(timegroup == max(timegroup), TRUE, FALSE), by = .(runid, dyadID)]

## Example output
dyads[dyadID == 'B-H', 
      .(timegroup, nObs, shifttimegrp, difftimegrp, runid, runCount, start, end)]
```

<!-- mean xy, todo -->

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_times.R
\name{group_times}
\alias{group_times}
\title{Group Times}
\usage{
group_times(DT = NULL, datetime = NULL, threshold = NULL)
}
\arguments{
\item{DT}{input data.table}

\item{datetime}{name of date time column(s). either 1 POSIXct or 2 IDate and
ITime. e.g.: 'datetime' or c('idate', 'itime')}

\item{threshold}{threshold for grouping times. e.g.: '2 hours', '10 minutes',
etc. if not provided, times will be matched exactly. Note that provided
threshold must be in the expected format: '## unit'}
}
\value{
\code{group_times} returns the input \code{DT} appended with a
\code{timegroup} column and additional temporal grouping columns to help
investigate, troubleshoot and interpret the timegroup.

The actual value of \code{timegroup} is arbitrary and represents the
identity of a given \code{timegroup} which 1 or more individuals are
assigned to. If the data was reordered, the group may change, but the
contents of each group would not.

The temporal grouping columns added depend on the \code{threshold}
provided:

\itemize{ \item \code{threshold} with unit minutes: "minutes" column added
identifying the nearest minute group for each row. \item \code{threshold}
with unit hours: "hours" column added identifying the nearest hour group
for each row. \item \code{threshold} with unit days: "block" columns added
identifying the multiday block for each row. }

A message is returned when any of these columns already exist in the input
\code{DT}, because they will be overwritten.
}
\description{
\code{group_times} groups rows into time groups. The function accepts date
time formatted data and a threshold argument. The threshold argument is used
to specify a time window within which rows are grouped.
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

The \code{datetime} argument expects the name of a column in \code{DT} which
is of type \code{POSIXct} or the name of two columns in \code{DT} which are
of type \code{IDate} and \code{ITime}.

\code{threshold} must be provided in units of minutes, hours or days. The
character string should start with an integer followed by a unit, separated
by a space. It is interpreted in terms of 24 hours which poses the following
limitations:

\itemize{ \item minutes, hours and days cannot be fractional \item minutes
must divide evenly into 60 \item minutes must not exceed 60 \item minutes,
hours which are nearer to the next day, are grouped as such \item hours must
divide evenly into 24 \item multi-day blocks should divide into the range of
days, else the blocks may not be the same length }

In addition, the \code{threshold} is considered a fixed window throughout the
time series and the rows are grouped to the nearest interval.

If \code{threshold} is NULL, rows are grouped using the \code{datetime}
column directly.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

group_times(DT, datetime = 'datetime', threshold = '5 minutes')

group_times(DT, datetime = 'datetime', threshold = '2 hours')

group_times(DT, datetime = 'datetime', threshold = '10 days')

}
\seealso{
\code{\link{group_pts}} \code{\link{group_lines}}
\code{\link{group_polys}}
}
\concept{Temporal grouping}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/randomizations.R
\name{randomizations}
\alias{randomizations}
\title{Data-stream randomizations}
\usage{
randomizations(
  DT = NULL,
  type = NULL,
  id = NULL,
  group = NULL,
  coords = NULL,
  datetime = NULL,
  splitBy = NULL,
  iterations = NULL
)
}
\arguments{
\item{DT}{input data.table}

\item{type}{one of 'daily', 'step' or 'trajectory' - see details}

\item{id}{Character string of ID column name}

\item{group}{generated from spatial grouping functions - see details}

\item{coords}{Character vector of X coordinate and Y coordinate column names}

\item{datetime}{field used for providing date time or time group - see
details}

\item{splitBy}{List of fields in DT to split the randomization process by}

\item{iterations}{The number of iterations to randomize}
}
\value{
\code{randomizations} returns the random date time or random id along
with the original \code{DT}, depending on the randomization \code{type}.
The length of the returned \code{data.table} is the original number of rows
multiplied by the number of iterations + 1. For example, 3 iterations will
return 4x - one observed and three randomized.

Two columns are always returned: \itemize{ \item observed - if the rows
represent the observed (TRUE/FALSE) \item iteration - iteration of rows
(where 0 is the observed) }

In addition, depending on the randomization type, random ID or random date
time columns are returned:

\itemize{ \item step - \code{randomID} each time step \item daily -
\code{randomID} for each day and \code{jul} indicating julian day \item
trajectory - a random date time ("random" prefixed to \code{datetime}
argument), observed \code{jul} and \code{randomJul} indicating the random
day relocations are swapped to. }
}
\description{
\code{randomizations} performs data-stream social network randomization. The
function accepts a \code{data.table} with relocation data, individual
identifiers and a randomization \code{type}. The \code{data.table} is
randomized either using \code{step} or \code{daily} between-individual
methods, or within-individual daily \code{trajectory} method described by
Spiegel et al. (2016).
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

Three randomization \code{type}s are provided: \enumerate{ \item step -
randomizes identities of relocations between individuals within each time
step. \item daily - randomizes identities of relocations between individuals
within each day. \item trajectory - randomizes daily trajectories within
individuals (Spiegel et al. 2016). }

Depending on the \code{type}, the \code{datetime} must be a certain format:

\itemize{ \item step - datetime is integer group created by
\code{group_times} \item daily - datetime is \code{POSIXct} format \item
trajectory - datetime is \code{POSIXct} format }

The \code{id}, \code{datetime},  (and optional \code{splitBy}) arguments
expect the names of respective columns in \code{DT} which correspond to the
individual identifier, date time, and additional grouping columns. The
\code{coords} argument is only required when the \code{type} is "trajectory",
since the coordinates are required for recalculating spatial groups with
\code{group_pts}, \code{group_lines} or \code{group_polys}.

Please note that if the data extends over multiple years, a column indicating
the year should be provided to the \code{splitBy} argument. This will ensure
randomizations only occur within each year.

The \code{group} argument is expected only when \code{type} is 'step' or
'daily'.

For example, using \code{\link[data.table:IDateTime]{data.table::year}}:

\preformatted{ DT[, yr := year(datetime)] randomizations(DT, type = 'step',
id = 'ID', datetime = 'timegroup', splitBy = 'yr') }

\code{iterations} is set to 1 if not provided. Take caution with a large
value for \code{iterations} with large input \code{DT}.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Date time columns
DT[, datetime := as.POSIXct(datetime)]
DT[, yr := year(datetime)]

# Temporal grouping
group_times(DT, datetime = 'datetime', threshold = '5 minutes')

# Spatial grouping with timegroup
group_pts(DT, threshold = 5, id = 'ID', coords = c('X', 'Y'), timegroup = 'timegroup')

# Randomization: step
randStep <- randomizations(
    DT,
    type = 'step',
    id = 'ID',
    group = 'group',
    datetime = 'timegroup',
    splitBy = 'yr',
    iterations = 2
)

# Randomization: daily
randDaily <- randomizations(
    DT,
    type = 'daily',
    id = 'ID',
    group = 'group',
    datetime = 'datetime',
    splitBy = 'yr',
    iterations = 2
)

# Randomization: trajectory
randTraj <- randomizations(
    DT,
    type = 'trajectory',
    id = 'ID',
    group = NULL,
    coords = c('X', 'Y'),
    datetime = 'datetime',
    splitBy = 'yr',
    iterations = 2
)

}
\references{
\url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.12553}
}
\seealso{
Other Social network tools: 
\code{\link{get_gbi}()}
}
\concept{Social network tools}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_polys.R
\name{group_polys}
\alias{group_polys}
\title{Group Polygons}
\usage{
group_polys(
  DT = NULL,
  area = NULL,
  hrType = NULL,
  hrParams = NULL,
  projection = NULL,
  id = NULL,
  coords = NULL,
  splitBy = NULL,
  spPolys = NULL
)
}
\arguments{
\item{DT}{input data.table}

\item{area}{boolean indicating either overlap group (when \code{FALSE}) or
area and proportion of overlap (when \code{TRUE})}

\item{hrType}{type of HR estimation, either 'mcp' or 'kernel'}

\item{hrParams}{a named list of parameters for \code{adehabitatHR} functions}

\item{projection}{character string defining the projection to be passed to
\code{sp::CRS}. For example, for UTM zone 36S (EPSG 32736),
the projection argument is 'EPSG:32736'. See details.}

\item{id}{Character string of ID column name}

\item{coords}{Character vector of X coordinate and Y coordinate column names}

\item{splitBy}{(optional) character string or vector of grouping column
name(s) upon which the grouping will be calculated}

\item{spPolys}{Alternatively, provide solely a SpatialPolygons object}
}
\value{
When \code{area} is \code{FALSE}, \code{group_polys} returns the
input \code{DT} appended with a \code{group} column. As with the other
grouping functions,  the actual value of \code{group} is arbitrary and
represents the identity of a given group where 1 or more individuals are
assigned to a group. If the data was reordered, the \code{group} may
change, but the contents of each group would not. When \code{area} is
\code{TRUE}, \code{group_polys} returns a proportional area overlap
\code{data.table}. In this case, ID refers to the focal individual of which
the total area is compared against the overlapping area of ID2.

If \code{area} is \code{FALSE}, a message is returned when a column named
\code{group} already exists in the input \code{DT}, because it will be
overwritten.
}
\description{
\code{group_polys} groups rows into spatial groups by overlapping polygons
(home ranges). The function accepts a \code{data.table} with relocation data,
individual identifiers and an \code{area} argument.  The relocation data is
transformed into home range \code{SpatialPolygons}. If the \code{area}
argument is \code{FALSE}, \code{group_polys} returns grouping calculated by
overlap. If the \code{area} argument is \code{TRUE}, the area and proportion
of overlap is calculated. Relocation data should be in two columns
representing the X and Y coordinates.
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

The \code{id}, \code{coords} (and optional \code{splitBy}) arguments expect
the names of respective columns in \code{DT} which correspond to the
individual identifier, X and Y coordinates, and additional grouping columns.

The \code{projection} argument expects a character string defining the EPSG
code. For example, for UTM zone 36N (EPSG 32736), the projection argument is
'EPSG:32736'. See \url{https://spatialreference.org} for a list of EPSG
codes. Please note, R spatial has followed updates to GDAL and PROJ for
handling projections, see more at
\url{https://www.r-spatial.org/r/2020/03/17/wkt.html}. It is likely
that \code{build_polys} will return "Warning in proj4string(xy) :
CRS object has comment, which is lost in output" due to these changes.

The \code{hrType} must be either one of "kernel" or "mcp". The
\code{hrParams} must be a named list of arguments matching those of
\code{adehabitatHR::kernelUD} or \code{adehabitatHR::mcp}.

The \code{splitBy} argument offers further control over grouping. If within
your \code{DT}, you have multiple populations, subgroups or other distinct
parts, you can provide the name of the column which identifies them to
\code{splitBy}. The grouping performed by \code{group_polys} will only
consider rows within each \code{splitBy} subgroup.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

# EPSG code for example data
utm <- 'EPSG:32736'

group_polys(DT, area = FALSE, hrType = 'mcp',
            hrParams = list(percent = 95), projection = utm,
            id = 'ID', coords = c('X', 'Y'))

areaDT <- group_polys(DT, area = TRUE, hrType = 'mcp',
                      hrParams = list(percent = 95), projection = utm,
                      id = 'ID', coords = c('X', 'Y'))
}
\seealso{
\code{\link{build_polys}} \code{\link{group_times}}

Other Spatial grouping: 
\code{\link{group_lines}()},
\code{\link{group_pts}()}
}
\concept{Spatial grouping}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_polys.R
\name{build_polys}
\alias{build_polys}
\title{Build Polygons}
\usage{
build_polys(
  DT = NULL,
  projection = NULL,
  hrType = NULL,
  hrParams = NULL,
  id = NULL,
  coords = NULL,
  splitBy = NULL,
  spPts = NULL
)
}
\arguments{
\item{DT}{input data.table}

\item{projection}{character string defining the projection to be passed to
\code{sp::CRS}. For example, for UTM zone 36S (EPSG 32736),
the projection argument is 'EPSG:32736'. See details.}

\item{hrType}{type of HR estimation, either 'mcp' or 'kernel'}

\item{hrParams}{a named list of parameters for \code{adehabitatHR} functions}

\item{id}{Character string of ID column name}

\item{coords}{Character vector of X coordinate and Y coordinate column names}

\item{splitBy}{(optional) character string or vector of grouping column
name(s) upon which the grouping will be calculated}

\item{spPts}{alternatively, provide solely a SpatialPointsDataFrame with one
column representing the ID of each point.}
}
\value{
\code{build_polys} returns a \code{SpatialPolygons} object
with a polyon for each individual (and optionally \code{splitBy}
combination).

An error is returned when \code{hrParams} do not match the arguments
of the \code{hrType} \code{adehabitatHR} function.
}
\description{
\code{build_polys} creates a \code{SpatialPolygons} object from a
\code{data.table}. The function accepts a \code{data.table} with
relocation data, individual identifiers, a \code{projection},
\code{hrType} and \code{hrParams}. The relocation data is transformed
into \code{SpatialPolygons} for each individual and optionally, each
\code{splitBy}. Relocation data should be in two columns representing
the X and Y coordinates.
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

The \code{id}, \code{coords} (and optional \code{splitBy}) arguments
expect the names of respective columns in \code{DT} which correspond
to the individual identifier, X and Y coordinates, and additional
grouping columns.

The \code{projection} argument expects a character string defining
the EPSG code. For example, for UTM zone 36N (EPSG 32736), the projection
argument is "EPSG:32736". See \url{https://spatialreference.org}
for a list of EPSG codes. Please note, R spatial has followed updates
to GDAL and PROJ for handling projections, see more at
\url{https://www.r-spatial.org/r/2020/03/17/wkt.html}. It is likely
that \code{build_polys} will return "Warning in proj4string(xy) :
CRS object has comment, which is lost in output" due to these changes.

The \code{hrType} must be either one of "kernel" or "mcp". The
\code{hrParams} must be a named list of arguments matching those
of \code{adehabitatHR::kernelUD} and \code{adehabitatHR::getverticeshr}
or \code{adehabitatHR::mcp}.

The \code{splitBy} argument offers further control building
\code{SpatialPolygons}. If in your \code{DT}, you have multiple
temporal groups (e.g.: years) for example, you can provide the
name of the column which identifies them and build \code{SpatialPolygons}
for each individual in each year.

\code{group_polys} uses \code{build_polys} for grouping overlapping
polygons created from relocations.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

# EPSG code for example data
utm <- 'EPSG:32736'

# Build polygons for each individual using kernelUD and getverticeshr
build_polys(DT, projection = utm, hrType = 'kernel',
            hrParams = list(grid = 60, percent = 95),
            id = 'ID', coords = c('X', 'Y'))

# Build polygons for each individual by year
DT[, yr := year(datetime)]
build_polys(DT, projection = utm, hrType = 'mcp', hrParams = list(percent = 95),
            id = 'ID', coords = c('X', 'Y'), splitBy = 'yr')

# Build polygons from SpatialPointsDataFrame
library(sp)
pts <- SpatialPointsDataFrame(coords = DT[, .(X, Y)],
                              proj4string = CRS(utm),
                              data = DT[, .(ID)]
)

build_polys(spPts = pts, hrType = 'mcp', hrParams = list(percent = 95))

}
\seealso{
\code{\link{group_polys}}

Other Build functions: 
\code{\link{build_lines}()}
}
\concept{Build functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edge_nn.R
\name{edge_nn}
\alias{edge_nn}
\title{Nearest neighbour based edge lists}
\usage{
edge_nn(
  DT = NULL,
  id = NULL,
  coords = NULL,
  timegroup,
  splitBy = NULL,
  threshold = NULL,
  returnDist = FALSE
)
}
\arguments{
\item{DT}{input data.table}

\item{id}{Character string of ID column name}

\item{coords}{Character vector of X coordinate and Y coordinate column names}

\item{timegroup}{timegroup field in the DT upon which the grouping will be
calculated}

\item{splitBy}{(optional) character string or vector of grouping column
name(s) upon which the grouping will be calculated}

\item{threshold}{(optional) spatial distance threshold to set maximum
distance between an individual and their neighbour.}

\item{returnDist}{boolean indicating if the distance between individuals
should be returned. If FALSE (default), only ID, NN columns (and timegroup,
splitBy columns if provided) are returned. If TRUE, another column
"distance" is returned indicating the distance between ID and NN.}
}
\value{
\code{edge_nn} returns a \code{data.table}  with three columns:
timegroup, ID and NN. If 'returnDist' is TRUE, column 'distance' is
returned indicating the distance between ID and NN.

The ID and NN columns represent the edges defined by the nearest neighbours
(and temporal thresholds with \code{group_times}).

If an individual was alone in a timegroup or splitBy, or did not have any
neighbours within the threshold distance, they are assigned NA for nearest
neighbour.
}
\description{
\code{edge_nn} returns edge lists defined by the nearest neighbour. The
function accepts a \code{data.table} with relocation data, individual
identifiers and a threshold argument. The threshold argument is used to
specify the criteria for distance between points which defines a group.
Relocation data should be in two columns representing the X and Y
coordinates.
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

The \code{id}, \code{coords} (and optional \code{timegroup} and
\code{splitBy}) arguments expect the names of a column in \code{DT} which
correspond to the individual identifier, X and Y coordinates, timegroup
(generated by \code{group_times}) and additional grouping columns.

The \code{threshold} must be provided in the units of the coordinates. The
\code{threshold} must be larger than 0. The coordinates must be planar
coordinates (e.g.: UTM). In the case of UTM, a \code{threshold} = 50 would
indicate a 50m distance threshold.

The \code{timegroup} argument is optional, but recommended to pair with
\code{\link{group_times}}. The intended framework is to group rows temporally
with \code{\link{group_times}} then spatially with \code{edge_nn} (or
grouping functions).

The \code{splitBy} argument offers further control over grouping. If within
your \code{DT}, you have multiple populations, subgroups or other distinct
parts, you can provide the name of the column which identifies them to
\code{splitBy}. \code{edge_nn} will only consider rows within each
\code{splitBy} subgroup.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

# Temporal grouping
group_times(DT, datetime = 'datetime', threshold = '20 minutes')

# Edge list generation
edges <- edge_nn(DT, id = 'ID', coords = c('X', 'Y'),
        timegroup = 'timegroup')

# Edge list generation using maximum distance threshold
edges <- edge_nn(DT, id = 'ID', coords = c('X', 'Y'),
        timegroup = 'timegroup', threshold = 100)

# Edge list generation, returning distance between nearest neighbours
edge_nn(DT, id = 'ID', coords = c('X', 'Y'),
        timegroup = 'timegroup', threshold = 100,
        returnDist = TRUE)

}
\seealso{
Other Edge-list generation: 
\code{\link{edge_dist}()}
}
\concept{Edge-list generation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/edge_dist.R
\name{edge_dist}
\alias{edge_dist}
\title{Distance based edge lists}
\usage{
edge_dist(
  DT = NULL,
  threshold,
  id = NULL,
  coords = NULL,
  timegroup,
  splitBy = NULL,
  returnDist = FALSE,
  fillNA = TRUE
)
}
\arguments{
\item{DT}{input data.table}

\item{threshold}{distance for grouping points, in the units of the
coordinates}

\item{id}{Character string of ID column name}

\item{coords}{Character vector of X coordinate and Y coordinate column names}

\item{timegroup}{timegroup field in the DT upon which the grouping will be
calculated}

\item{splitBy}{(optional) character string or vector of grouping column
name(s) upon which the grouping will be calculated}

\item{returnDist}{boolean indicating if the distance between individuals
should be returned. If FALSE (default), only ID1, ID2 columns (and
timegroup, splitBy columns if provided) are returned. If TRUE, another
column "distance" is returned indicating the distance between ID1 and ID2.}

\item{fillNA}{boolean indicating if NAs should be returned for individuals
that were not within the threshold distance of any other. If TRUE, NAs are
returned. If FALSE, only edges between individuals within the threshold
distance are returned.}
}
\value{
\code{edge_dist} returns a \code{data.table} with columns ID1, ID2,
timegroup (if supplied) and any columns provided in splitBy. If
'returnDist' is TRUE, column 'distance' is returned indicating the distance
between ID1 and ID2.

The ID1 and ID2 columns represent the edges defined by the spatial (and
temporal with \code{group_times}) thresholds.
}
\description{
\code{edge_dist} returns edge lists defined by a spatial distance within the
user defined threshold. The function accepts a \code{data.table} with
relocation data, individual identifiers and a threshold argument. The
threshold argument is used to specify the criteria for distance between
points which defines a group. Relocation data should be in two columns
representing the X and Y coordinates.
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

The \code{id}, \code{coords} (and optional \code{timegroup} and
\code{splitBy}) arguments expect the names of a column in \code{DT} which
correspond to the individual identifier, X and Y coordinates, timegroup
(generated by \code{group_times}) and additional grouping columns.

If provided, the \code{threshold} must be provided in the units of the coordinates and must be larger than 0.
If the \code{threshold} is NULL, the distance to all other individuals will be returned. The coordinates must be planar
coordinates (e.g.: UTM). In the case of UTM, a \code{threshold} = 50 would
indicate a 50m distance threshold.

The \code{timegroup} argument is optional, but recommended to pair with
\code{\link{group_times}}. The intended framework is to group rows temporally
with \code{\link{group_times}} then spatially with \code{edge_dist} (or
grouping functions).

The \code{splitBy} argument offers further control over grouping. If within
your \code{DT}, you have multiple populations, subgroups or other distinct
parts, you can provide the name of the column which identifies them to
\code{splitBy}. \code{edge_dist} will only consider rows within each
\code{splitBy} subgroup.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

# Temporal grouping
group_times(DT, datetime = 'datetime', threshold = '20 minutes')

# Edge list generation
edges <- edge_dist(
    DT,
    threshold = 100,
    id = 'ID',
    coords = c('X', 'Y'),
    timegroup = 'timegroup',
    returnDist = TRUE,
    fillNA = TRUE
  )
}
\seealso{
Other Edge-list generation: 
\code{\link{edge_nn}()}
}
\concept{Edge-list generation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dyads.R
\name{dyad_id}
\alias{dyad_id}
\title{Dyad ID}
\usage{
dyad_id(DT = NULL, id1 = NULL, id2 = NULL)
}
\arguments{
\item{DT}{input data.table with columns id1 and id2, as generated by
\code{edge_dist} or \code{edge_nn}}

\item{id1}{ID1 column name generated by \code{edge_dist} or \code{edge_nn}}

\item{id2}{ID2 column name generated by \code{edge_dist} or \code{edge_nn}}
}
\value{
\code{dyad_id} returns the input \code{data.table} with appended "dyadID"
column
}
\description{
Generate a dyad ID for edge list generated by \code{\link{edge_nn}} or
\code{\link{edge_dist}}.
}
\details{
An undirected edge identifier between, for example individuals A and B will
be A-B (and reverse B and A will be A-B). Internally sorts and pastes id
columns.

More details in the edge and dyad vignette (in progress).
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

# Temporal grouping
group_times(DT, datetime = 'datetime', threshold = '20 minutes')

# Edge list generation
edges <- edge_dist(
    DT,
    threshold = 100,
    id = 'ID',
    coords = c('X', 'Y'),
    timegroup = 'timegroup',
    returnDist = TRUE,
    fillNA = TRUE
  )

# Generate dyad IDs
dyad_id(edges, 'ID1', 'ID2')
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_lines.R
\name{group_lines}
\alias{group_lines}
\title{Group Lines}
\usage{
group_lines(
  DT = NULL,
  threshold = NULL,
  projection = NULL,
  id = NULL,
  coords = NULL,
  timegroup = NULL,
  sortBy = NULL,
  splitBy = NULL,
  spLines = NULL
)
}
\arguments{
\item{DT}{input data.table}

\item{threshold}{The width of the buffer around the lines in the units of the
projection. Supply 0 to compare intersection without buffering.}

\item{projection}{character string defining the projection to be passed to
\code{sp::CRS}. For example, for UTM zone 36S (EPSG 32736),
the projection argument is 'EPSG:32736'. See details.}

\item{id}{Character string of ID column name}

\item{coords}{Character vector of X coordinate and Y coordinate column names}

\item{timegroup}{timegroup field in the DT upon which the grouping will be
calculated}

\item{sortBy}{Character string of date time column(s) to sort rows by. Must
be a POSIXct.}

\item{splitBy}{(optional) character string or vector of grouping column
name(s) upon which the grouping will be calculated}

\item{spLines}{Alternatively to providing a DT, provide a SpatialLines object
created with the sp package. If a spLines object is provided, groups cannot
be calculated by a    timegroup or splitBy.}
}
\value{
\code{group_lines} returns the input \code{DT} appended with a
\code{group} column.

This column represents the spatial (and if \code{timegroup} was provided -
spatiotemporal) group calculated by overlapping lines. As with the other
grouping functions,  the actual value of group is arbitrary and represents
the identity of a given group where 1 or more individuals are assigned to a
group. If the data was reordered, the group may change, but the contents of
each group would not.

A message is returned when a column named \code{group} already exists in
the input \code{DT}, because it will be overwritten.
}
\description{
\code{group_lines} groups rows into spatial groups by creating trajectories
and grouping based on spatial overlap. The function accepts a
\code{data.table} with relocation data, individual identifiers and a
\code{threshold}. The relocation data is transformed into \code{SpatialLines}
and overlapping \code{SpatialLines} are grouped. The \code{threshold}
argument is used to specify the criteria for distance between lines.
Relocation data should be in two columns representing the X and Y
coordinates.
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

The \code{id}, \code{coords}, \code{sortBy} (and optional \code{timegroup}
and \code{splitBy}) arguments expect the names of respective columns in
\code{DT} which correspond to the individual identifier, X and Y coordinates,
sorting, timegroup (generated by \code{group_times}) and additional grouping
columns.

The \code{projection} argument expects a character string defining
the EPSG code. For example, for UTM zone 36N (EPSG 32736), the projection
argument is "EPSG:32736". See \url{https://spatialreference.org}
for a list of EPSG codes. Please note, R spatial has followed updates
to GDAL and PROJ for handling projections, see more at
\url{https://www.r-spatial.org/r/2020/03/17/wkt.html}. It is likely
that \code{build_polys} will return "Warning in proj4string(xy) :
CRS object has comment, which is lost in output" due to these changes.

The \code{sortBy} is used to order the input \code{data.table} when creating
\code{SpatialLines}. It must a \code{POSIXct} to ensure the rows are sorted
by date time.

The \code{threshold} must be provided in the units of the coordinates. The
\code{threshold} can be equal to 0 if strict overlap is required, else it
needs to be greater than 0. The coordinates must be planar coordinates (e.g.:
UTM). In the case of UTM, a \code{threshold} = 50 would indicate a 50m
distance threshold.

The \code{timegroup} argument is optional, but recommended to pair with
\code{\link{group_times}}. The intended framework is to group rows temporally
with \code{\link{group_times}} then spatially with \code{group_lines} (or
\code{\link{group_pts}}, \code{\link{group_polys}}). With \code{group_lines},
pick a relevant \code{group_times} \code{threshold} such as \code{'1 day'} or
\code{'7 days'} which is informed by your study species and system.

The \code{splitBy} argument offers further control over grouping. If within
your \code{DT}, you have multiple populations, subgroups or other distinct
parts, you can provide the name of the column which identifies them to
\code{splitBy}. The grouping performed by \code{group_lines} will only
consider rows within each \code{splitBy} subgroup.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Subset only individuals A, B, and C
DT <- DT[ID \%in\% c('A', 'B', 'C')]

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

# EPSG code for example data
utm <- 'EPSG:32736'

\donttest{group_lines(DT, threshold = 50, projection = utm, sortBy = 'datetime',
            id = 'ID', coords = c('X', 'Y'))}

## Daily movement tracks
# Temporal grouping
group_times(DT, datetime = 'datetime', threshold = '1 day')

# Subset only first 50 days
DT <- DT[timegroup < 25]

# Spatial grouping
group_lines(DT, threshold = 50, projection = utm,
            id = 'ID', coords = c('X', 'Y'),
            timegroup = 'timegroup', sortBy = 'datetime')

## Daily movement tracks by population
group_lines(DT, threshold = 50, projection = utm,
            id = 'ID', coords = c('X', 'Y'),
            timegroup = 'timegroup', sortBy = 'datetime',
            splitBy = 'population')
}
\seealso{
\code{\link{build_lines}} \code{\link{group_times}}

Other Spatial grouping: 
\code{\link{group_polys}()},
\code{\link{group_pts}()}
}
\concept{Spatial grouping}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/group_pts.R
\name{group_pts}
\alias{group_pts}
\title{Group Points}
\usage{
group_pts(
  DT = NULL,
  threshold = NULL,
  id = NULL,
  coords = NULL,
  timegroup,
  splitBy = NULL
)
}
\arguments{
\item{DT}{input data.table}

\item{threshold}{distance for grouping points, in the units of the
coordinates}

\item{id}{Character string of ID column name}

\item{coords}{Character vector of X coordinate and Y coordinate column names}

\item{timegroup}{timegroup field in the DT upon which the grouping will be
calculated}

\item{splitBy}{(optional) character string or vector of grouping column
name(s) upon which the grouping will be calculated}
}
\value{
\code{group_pts} returns the input \code{DT} appended with a
\code{group} column.

This column represents the spatial (and if \code{timegroup} was provided -
spatiotemporal) group. As with the other grouping functions,  the actual
value of \code{group} is arbitrary and represents the identity of a given
group where 1 or more individuals are assigned to a group. If the data was
reordered, the \code{group} may change, but the contents of each group
would not.

A message is returned when a column named \code{group} already exists in
the input \code{DT}, because it will be overwritten.
}
\description{
\code{group_pts} groups rows into spatial groups. The function accepts a
\code{data.table} with relocation data, individual identifiers and a
threshold argument. The threshold argument is used to specify the criteria
for distance between points which defines a group. Relocation data should be
in two columns representing the X and Y coordinates.
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

The \code{id}, \code{coords} (and optional \code{timegroup} and
\code{splitBy}) arguments expect the names of a column in \code{DT} which
correspond to the individual identifier, X and Y coordinates, timegroup
(generated by \code{group_times}) and additional grouping columns.

The \code{threshold} must be provided in the units of the coordinates. The
\code{threshold} must be larger than 0. The coordinates must be planar
coordinates (e.g.: UTM). In the case of UTM, a \code{threshold} = 50 would
indicate a 50m distance threshold.

The \code{timegroup} argument is optional, but recommended to pair with
\code{\link{group_times}}. The intended framework is to group rows temporally
with \code{\link{group_times}} then spatially with \code{group_pts} (or
\code{\link{group_lines}}, \code{\link{group_polys}}).

The \code{splitBy} argument offers further control over grouping. If within
your \code{DT}, you have multiple populations, subgroups or other distinct
parts, you can provide the name of the column which identifies them to
\code{splitBy}. The grouping performed by \code{group_pts} will only consider
rows within each \code{splitBy} subgroup.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

# Temporal grouping
group_times(DT, datetime = 'datetime', threshold = '20 minutes')

# Spatial grouping with timegroup
group_pts(DT, threshold = 5, id = 'ID',
          coords = c('X', 'Y'), timegroup = 'timegroup')

# Spatial grouping with timegroup and splitBy on population
group_pts(DT, threshold = 5, id = 'ID', coords = c('X', 'Y'),
         timegroup = 'timegroup', splitBy = 'population')
}
\seealso{
\code{\link{group_times}}

Other Spatial grouping: 
\code{\link{group_lines}()},
\code{\link{group_polys}()}
}
\concept{Spatial grouping}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_gbi.R
\name{get_gbi}
\alias{get_gbi}
\title{Generate group by individual matrix}
\usage{
get_gbi(DT = NULL, group = "group", id = NULL)
}
\arguments{
\item{DT}{input data.table}

\item{group}{Character string of group column (generated from one of
spatsoc's spatial grouping functions)}

\item{id}{Character string of ID column name}
}
\value{
\code{get_gbi} returns a group by individual matrix (columns
represent individuals and rows represent groups).

Note that \code{get_gbi} is identical in function for turning the outputs
of \code{spatsoc} into social networks as
\code{\link[asnipe:get_group_by_individual]{asnipe::get_group_by_individual}}
but is more efficient thanks to
\code{\link[data.table:dcast.data.table]{data.table::dcast}}.
}
\description{
\code{get_gbi} generates a group by individual matrix. The function accepts a
\code{data.table} with individual identifiers and a group column. The group
by individual matrix can then be used to build a network using
\code{\link[asnipe:get_network]{asnipe::get_network}}.
}
\details{
The \code{DT} must be a \code{data.table}. If your data is a
\code{data.frame}, you can convert it by reference using
\code{\link[data.table:setDT]{data.table::setDT}}.

The \code{group} argument expects the name of a column which corresponds to
an integer group identifier (generated by \code{\link{spatsoc}}'s grouping
functions).

The \code{id} argument expects the name of a column which corresponds to the
individual identifier.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]
DT[, yr := year(datetime)]

# EPSG code for example data
utm <- 'EPSG:32736'

group_polys(DT, area = FALSE, hrType = 'mcp',
            hrParams = list(percent = 95),
            projection = utm, id = 'ID', coords = c('X', 'Y'),
            splitBy = 'yr')

gbiMtrx <- get_gbi(DT = DT, group = 'group', id = 'ID')

}
\seealso{
\code{\link{group_pts}} \code{\link{group_lines}}
\code{\link{group_polys}}

Other Social network tools: 
\code{\link{randomizations}()}
}
\concept{Social network tools}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extdata.R
\name{DT}
\alias{DT}
\title{Movement of 10 "Newfoundland Bog Cows"}
\format{
A data.table with 14297 rows and 5 variables: \describe{
\item{ID}{individual identifier} \item{X}{X coordinate of the relocation
(UTM 36N)} \item{Y}{Y coordinate of the relocation (UTM 36N)}
\item{datetime}{character string representing the date time}
\item{population}{sub population within the individuals} }
}
\description{
A dataset containing the GPS relocations of 10 individuals in winter
2016-2017.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/spatsoc.R
\docType{package}
\name{spatsoc}
\alias{spatsoc}
\alias{_PACKAGE}
\alias{spatsoc-package}
\title{spatsoc}
\description{
spatsoc is an R package for detecting spatial and temporal groups in GPS
relocations. It can be used to convert GPS relocations to gambit-of-the-group
format to build proximity-based social networks. In addition, the
randomization function provides data-stream randomization methods suitable
for GPS data.
}
\details{
The spatsoc package provides one temporal grouping function:

\itemize{ \item \code{\link{group_times}} } three spatial grouping functions:
\itemize{ \item \code{\link{group_pts}} \item \code{\link{group_lines}} \item
\code{\link{group_polys}} }

two edge list generating functions:

\itemize{ \item \code{\link{edge_dist}} \item \code{\link{edge_nn}} }

and two social network functions: \itemize{ \item
\code{\link{randomizations}} \item \code{\link{get_gbi}} }
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/spatsoc/}
  \item \url{https://github.com/ropensci/spatsoc}
  \item \url{http://spatsoc.robitalec.ca}
  \item Report bugs at \url{https://github.com/ropensci/spatsoc/issues}
}

}
\author{
\strong{Maintainer}: Alec L. Robitaille \email{robit.alec@gmail.com} (\href{https://orcid.org/0000-0002-4706-1762}{ORCID})

Authors:
\itemize{
  \item Quinn Webber (\href{https://orcid.org/0000-0002-0434-9360}{ORCID})
  \item Eric Vander Wal (\href{https://orcid.org/0000-0002-8534-4317}{ORCID})
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/build_lines.R
\name{build_lines}
\alias{build_lines}
\title{Build Lines}
\usage{
build_lines(
  DT = NULL,
  projection = NULL,
  id = NULL,
  coords = NULL,
  sortBy = NULL,
  splitBy = NULL
)
}
\arguments{
\item{DT}{input data.table}

\item{projection}{character string defining the projection to be passed to
\code{sp::CRS}. For example, for UTM zone 36S (EPSG 32736),
the projection argument is 'EPSG:32736'. See details.}

\item{id}{Character string of ID column name}

\item{coords}{Character vector of X coordinate and Y coordinate column names}

\item{sortBy}{Character string of date time column(s) to sort rows by. Must
be a POSIXct.}

\item{splitBy}{(optional) character string or vector of grouping column
name(s) upon which the grouping will be calculated}
}
\value{
\code{build_lines} returns a \code{SpatialLines} object with a line
for each individual (and optionally \code{splitBy} combination).

An error is returned when an individual has less than 2 relocations, making
it impossible to build a line.
}
\description{
\code{build_lines} creates a \code{SpatialLines} object from a \code{data.table}.
The function accepts a \code{data.table} with relocation data, individual
identifiers a sorting column and a \code{projection}. The relocation data
is transformed into \code{SpatialLines} for each individual and optionally,
each \code{splitBy}. Relocation data should be in two columns representing
the X and Y coordinates.
}
\details{
The \code{projection} argument expects a character string defining the EPSG
code. For example, for UTM zone 36N (EPSG 32736), the projection argument is
'EPSG:32736'. See \url{https://spatialreference.org} for a list of
EPSG codes. Please note, R spatial has followed updates to GDAL and PROJ
for handling projections, see more at
\url{https://www.r-spatial.org/r/2020/03/17/wkt.html}.

The \code{sortBy} is used to order the input \code{data.table} when creating
\code{SpatialLines}. It must a \code{POSIXct} to ensure the rows are sorted
by date time.

The \code{splitBy} argument offers further control building \code{SpatialLines}.
If in your \code{DT}, you have multiple temporal groups (e.g.: years) for
example, you can provide the name of the column which identifies them and
build \code{SpatialLines} for each individual in each year.

\code{build_lines} is used by \code{group_lines} for grouping overlapping
lines created from relocations.
}
\examples{
# Load data.table
library(data.table)

# Read example data
DT <- fread(system.file("extdata", "DT.csv", package = "spatsoc"))

# Cast the character column to POSIXct
DT[, datetime := as.POSIXct(datetime, tz = 'UTC')]

# EPSG code for example data
utm <- 'EPSG:32736'

# Build lines for each individual
lines <- build_lines(DT, projection = utm, id = 'ID', coords = c('X', 'Y'),
            sortBy = 'datetime')

# Build lines for each individual by year
DT[, yr := year(datetime)]
lines <- build_lines(DT, projection = utm, id = 'ID', coords = c('X', 'Y'),
            sortBy = 'datetime', splitBy = 'yr')

}
\seealso{
\code{\link{group_lines}}

Other Build functions: 
\code{\link{build_polys}()}
}
\concept{Build functions}
