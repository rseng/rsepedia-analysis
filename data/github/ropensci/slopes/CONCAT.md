slopes: An R Package to Calculate Slopes of Roads, Rivers and
Trajectories
================
27 September 2021

# Summary

This package provides functions and example data to support research
into the slope (also known as longitudinal gradient or steepness) of
linear geographic entities such as roads (Ariza-López et al. 2019) and
rivers (Cohen et al. 2018). The package was initially developed to
calculate the steepness of street segments but can be used to calculate
steepness of any linear feature that can be represented as LINESTRING
geometries in the ‘sf’ class system (Pebesma 2018). The package takes
two main types of input data for slope calculation: vector geographic
objects representing linear features, and raster geographic objects with
elevation values (which can be downloaded using functionality in the
package) representing a continuous terrain surface. Where no raster
object is provided the package attempts to download elevation data using
the ‘ceramic’ package.

# Statement of need

Although there are several ways to name “slope,” such as “steepness,”
“hilliness,” “inclination,” “aspect,” “gradient,” “declivity,” the
referred `slopes` in this package can be defined as the “longitudinal
gradient” of linear geographic entities, as defined in the context of
rivers by(Cohen et al. 2018).

The package was initially developed to research road slopes to support
evidence-based sustainable transport policies. Accounting for gradient
when planning for new cycling infrastructure and road space reallocation
for walking and cycling can improve outcomes, for example by helping to
identify routes that avoid steep hills. The package can be used to
calculate and visualise slopes of rivers and trajectories representing
movement on roads of the type published as open data by Ariza-López et
al. (2019).

Data on slopes are useful in many fields of research, including
[hydrology](https://en.wikipedia.org/wiki/Stream_gradient), natural
hazards (including
[flooding](https://www.humanitarianresponse.info/fr/operations/afghanistan/infographic/afg-river-gradient-and-flood-hazard)
and [landslide risk
management](https://assets.publishing.service.gov.uk/media/57a08d0740f0b652dd0016f4/R7815-ADD017_col.pdf)),
recreational and competitive sports such as
[cycling](http://theclimbingcyclist.com/gradients-and-cycling-an-introduction/),
[hiking](https://trailism.com/trail-grades/), and
[skiing](https://www.snowplaza.co.uk/blog/16682-skiing-steeps-what-does-gradient-mean-ski-piste/).
Slopes are also also important in some branches of [transport and
emissions
modelling](https://www.sciencedirect.com/science/article/pii/S2352146516302642)
and [ecology](https://doi.org/10.1016/j.ncon.2016.10.001). A growing
number of people working with geospatial data require accurate estimates
of gradient, including:

-   Transport planning practitioners who require accurate estimates of
    roadway gradient for estimating energy consumption, safety and mode
    shift potential in hilly cities (such as Lisbon, the case study city
    used in the examples in the documentation).
-   Vehicle routing software developers, who need to build systems are
    sensitive to going up or down steep hills (e.g. bicycles, trains,
    and large trucks), such as active travel planning, logistics, and
    emergency services.
-   Natural hazard researchers and risk assessors require estimates of
    linear gradient to inform safety and mitigation plans associated
    with project on hilly terrain.
-   Aquatic ecologists, flooding researchers and others, who could
    benefit from estimates of river gradient to support modelling of
    storm hydrographs

There likely other domains where slopes could be useful, such as
agriculture, geology, and civil engineering.

An example of the demand for data provided by the package is a map
showing gradients across Sao Paulo (Brazil, see image below) that has
received more than 300 ‘likes’ on Twitter and generated conversations:
<https://twitter.com/DanielGuth/status/1347270685161304069>

<img src="https://camo.githubusercontent.com/30a3b814dd72aef5b51db635f2ab6e1b6b6c57b856d239822788967a4932d655/68747470733a2f2f7062732e7477696d672e636f6d2f6d656469612f45724a32647238574d414948774d6e3f666f726d61743d6a7067266e616d653d6c61726765" style="width:50.0%" />

# Usage and Key functions

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/slopes")
```

### Installation for DEM downloads

If you do not already have DEM data and want to make use of the
package’s ability to download them using the `ceramic` package, install
the package with suggested dependencies, as follows:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/slopes", dependencies = "Suggests")
```

Furthermore, you will need to add a MapBox API key to be able to get DEM
datasets, by signing up and registering for a key at
<https://account.mapbox.com/access-tokens/> and then following these
steps:

``` r
usethis::edit_r_environ()
MAPBOX_API_KEY=xxxxx # replace XXX with your api key
```

The key functions in the package are `elevation_add()`, which adds a
third ‘Z’ coordinate value for each vertex defining LINESTRING objects,
and `slope_xyz()` which calculates slopes for each linear feature in a
simple features object.

By default, the elevation of each vertex is estimated using [bilinear
interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation)
(`method = "bilinear"`) which calculates point height based on proximity
to the centroids of surrounding cells. The value of the `method`
argument is passed to the `method` argument in
[`raster::extract()`](https://rspatial.github.io/raster/reference/extract.html)
or
[`terra::extract()`](https://rspatial.github.io/terra/reference/extract.html)
depending on the class of the input raster dataset. See Kidner, Dorey,
and Smith (1999) for descriptions of alternative elevation interpolation
and extrapolation algorithms.

<!-- # Calculating slopes on regional transport networks -->
<!-- # Acknowledgements -->

# References

<div id="refs" class="references csl-bib-body hanging-indent">

<div id="ref-ariza-lopez_dataset_2019" class="csl-entry">

Ariza-López, Francisco Javier, Antonio Tomás Mozas-Calvache, Manuel
Antonio Ureña-Cámara, and Paula Gil de la Vega. 2019. “Dataset of
Three-Dimensional Traces of Roads.” *Scientific Data* 6 (1): 1–10.
<https://doi.org/10.1038/s41597-019-0147-x>.

</div>

<div id="ref-cohen_global_2018" class="csl-entry">

Cohen, Sagy, Tong Wan, Md Tazmul Islam, and J. P. M. Syvitski. 2018.
“Global River Slope: A New Geospatial Dataset and Global-Scale
Analysis.” *Journal of Hydrology* 563 (August): 1057–67.
<https://doi.org/10.1016/j.jhydrol.2018.06.066>.

</div>

<div id="ref-kidner_what_1999" class="csl-entry">

Kidner, David, Mark Dorey, and Derek Smith. 1999. “GeoComputation.” In.
Vol. 99. Mary Washington College Fredericksburg, Virginia, USA:
www.geocomputation.org.
<http://www.geocomputation.org/1999/082/gc_082.htm>.

</div>

<div id="ref-pebesma_simple_2018" class="csl-entry">

Pebesma, Edzer. 2018. “Simple Features for R: Standardized Support for
Spatial Vector Data.” *The R Journal*.
<https://journal.r-project.org/archive/2018/RJ-2018-009/index.html>.

</div>

</div>

<!-- README.md is generated from README.Rmd. Please edit that file -->

# slopes package

<!-- badges: start -->

[![R-CMD-check](https://github.com/ropensci/slopes/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/slopes/actions)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/slopes/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/slopes?branch=master)
[![Status at rOpenSci Software Peer
Review](https://badges.ropensci.org/420_status.svg)](https://github.com/ropensci/software-review/issues/420)
<!-- badges: end -->

The **slopes** R package calculates the slope (longitudinal steepness,
also known as gradient) of roads, rivers and other linear (simple)
features, based on two main inputs:

-   [vector](https://geocompr.robinlovelace.net/spatial-class.html#vector-data)
    linestring geometries defined by classes in the
    [`sf`](https://r-spatial.github.io/sf/) package
-   [raster](https://geocompr.robinlovelace.net/spatial-class.html#raster-data)
    objects with pixel values reporting average height, commonly known
    as digital elevation model (**DEM**) datasets, defined by classes in
    the [`raster`](https://cran.r-project.org/package=raster) or more
    recent [`terra`](https://rspatial.org/terra) packages

Data on slopes are useful in many fields of research, including
[hydrology](https://en.wikipedia.org/wiki/Stream_gradient), natural
hazards (including
[flooding](https://www.humanitarianresponse.info/fr/operations/afghanistan/infographic/afg-river-gradient-and-flood-hazard)
and [landslide risk
management](https://assets.publishing.service.gov.uk/media/57a08d0740f0b652dd0016f4/R7815-ADD017_col.pdf)),
recreational and competitive sports such as
[cycling](http://theclimbingcyclist.com/gradients-and-cycling-an-introduction/),
[hiking](https://trailism.com/trail-grades/), and
[skiing](https://www.snowplaza.co.uk/blog/16682-skiing-steeps-what-does-gradient-mean-ski-piste/).
Slopes are also also important in some branches of [transport and
emissions
modelling](https://www.sciencedirect.com/science/article/pii/S2352146516302642)
and [ecology](https://doi.org/10.1016/j.ncon.2016.10.001). See the
[`intro-to-slopes`
vignette](https://ropensci.github.io/slopes/articles/intro-to-slopes.html)
for details on fields using slope data and the need for this package.

This README covers installation and basic usage. For more information
about slopes and how to use the package to calculate them, see the [get
started](https://ropensci.github.io/slopes/) and the [introducion to
slopes](https://ropensci.github.io/intro-to-slopes/) vignette.

## How it works

The package takes two main types of input data for slope calculation: -
vector geographic objects representing **linear features**, and -
**elevation values** from a digital elevation model representing a
continuous terrain surface or which can be downloaded using
functionality in the package

The package can be used with two sources of elevation data: - openly
available elevation data via an interface to the [ceramic
package](https://github.com/hypertidy/ceramic), enabling estimation of
hilliness for routes anywhere worldwide even when local elevation data
is lacking. The package takes geographic lines objects and returns
elevation data per vertex (providing the output as a 3D point geometry
in the sf package by default) and per line feature (providing average
gradient by default). - an elevation model, available on your machine.

## Getting started

### Installation

<!-- You can install the released version of slopes from [CRAN](https://CRAN.R-project.org) with: -->
<!-- ``` r -->
<!-- install.packages("slopes") -->
<!-- ``` -->

Install the development version from [GitHub](https://github.com/) with:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/slopes")
```

#### Installation for DEM downloads

If you do not already have DEM data and want to make use of the
package’s ability to download them using the `ceramic` package, install
the package with suggested dependencies, as follows:

``` r
# install.packages("remotes")
remotes::install_github("ropensci/slopes", dependencies = "Suggests")
```

Furthermore, you will need to add a MapBox API key to be able to get DEM
datasets, by signing up and registering for a key at
<https://account.mapbox.com/access-tokens/> and then following these
steps:

``` r
usethis::edit_r_environ()
MAPBOX_API_KEY=xxxxx # replace XXX with your api key
```

## Basic examples

Load the package in the usual way. We will also load the `sf` library:

``` r
library(slopes)
library(sf)
```

The minimum input data requirement for using the package is an `sf`
object containing LINESTRING geometries, as illustrated below (requires
a MapBox API key):

``` r
sf_linestring = lisbon_route # import or load a linestring object
```

``` r
sf_linestring_xyz = elevation_add(sf_linestring)  # dem = NULL
#> Loading required namespace: ceramic
#> Preparing to download: 12 tiles at zoom = 12 from 
#> https://api.mapbox.com/v4/mapbox.terrain-rgb/
```

With the default argument `dem = NULL`, the function downloads the
necessary elevation information from Mapbox. You can also this use a
local digital elevation model (`dem = ...`), as shown in the example
below:

``` r
sf_linestring_xyz_local = elevation_add(sf_linestring, dem = dem_lisbon_raster)
```

In both cases you can obtain the average gradient of the linestring with
`slope_xyz()` and plot the elevation profile with `plot_slope()` as
follows:

``` r
slope_xyz(sf_linestring_xyz_local)
#>          1 
#> 0.07817098
plot_slope(sf_linestring_xyz_local)
```

<img src="man/figures/README-elevationprofile-1.png" width="100%" />

*See more functions in [Get
started](https://ropensci.github.io/slopes/articles/slopes.html)
vignette.*

## See more in vignettes

-   [Get
    started](https://ropensci.github.io/slopes/articles/slopes.html)
-   [An introduction to
    slopes](https://ropensci.github.io/slopes/articles/intro-to-slopes.html)
-   [Reproducible example: gradients of a road network for a given
    city](https://ropensci.github.io/slopes/articles/roadnetworkcycling.html)
-   [Verification of
    slopes](https://ropensci.github.io/slopes/articles/verification.html)
-   [Benchmarking slopes
    calculation](https://ropensci.github.io/slopes/articles/benchmark.html)

## Code of Conduct

Please note that this package is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project, you agree to abide by its terms.
# slopes 1.0.1

* Package source code now hosted at https://github.com/ropensci/slopes
* New documentation section showing how directed routes work: https://docs.ropensci.org/slopes/articles/slopes.html#splitting-the-network

# slopes 1.0.0

* We have submitted and responded to all comments to [rOpenSci review](https://github.com/ropensci/software-review/issues/420).  
* Many changes, including breaking changes to function names.  
* Added a `NEWS.md` file to track changes to the package.  


# slopes 0.0.1

* Initial version of the package on GitHub.  
# Contributor Covenant Code of Conduct

## Our Pledge

We as members, contributors, and leaders pledge to make participation in our
community a harassment-free experience for everyone, regardless of age, body
size, visible or invisible disability, ethnicity, sex characteristics, gender
identity and expression, level of experience, education, socio-economic status,
nationality, personal appearance, race, religion, or sexual identity and
orientation.

We pledge to act and interact in ways that contribute to an open, welcoming,
diverse, inclusive, and healthy community.

## Our Standards

Examples of behavior that contributes to a positive environment for our
community include:

* Demonstrating empathy and kindness toward other people
* Being respectful of differing opinions, viewpoints, and experiences
* Giving and gracefully accepting constructive feedback
* Accepting responsibility and apologizing to those affected by our mistakes,
and learning from the experience
* Focusing on what is best not just for us as individuals, but for the overall
community

Examples of unacceptable behavior include:

* The use of sexualized language or imagery, and sexual attention or
advances of any kind
* Trolling, insulting or derogatory comments, and personal or political attacks
* Public or private harassment
* Publishing others' private information, such as a physical or email
address, without their explicit permission
* Other conduct which could reasonably be considered inappropriate in a
professional setting

## Enforcement Responsibilities

Community leaders are responsible for clarifying and enforcing our standards
of acceptable behavior and will take appropriate and fair corrective action in
response to any behavior that they deem inappropriate, threatening, offensive,
or harmful.

Community leaders have the right and responsibility to remove, edit, or reject
comments, commits, code, wiki edits, issues, and other contributions that are
not aligned to this Code of Conduct, and will communicate reasons for moderation
decisions when appropriate.

## Scope

This Code of Conduct applies within all community spaces, and also applies
when an individual is officially representing the community in public spaces.
Examples of representing our community include using an official e-mail
address, posting via an official social media account, or acting as an appointed
representative at an online or offline event.

## Enforcement

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported to the community leaders responsible for enforcement at [INSERT CONTACT
METHOD]. All complaints will be reviewed and investigated promptly and fairly.

All community leaders are obligated to respect the privacy and security of the
reporter of any incident.

## Enforcement Guidelines

Community leaders will follow these Community Impact Guidelines in determining
the consequences for any action they deem in violation of this Code of Conduct:

### 1. Correction

**Community Impact**: Use of inappropriate language or other behavior deemed
unprofessional or unwelcome in the community.

**Consequence**: A private, written warning from community leaders, providing
clarity around the nature of the violation and an explanation of why the
behavior was inappropriate. A public apology may be requested.

### 2. Warning

**Community Impact**: A violation through a single incident or series of
actions.

**Consequence**: A warning with consequences for continued behavior. No
interaction with the people involved, including unsolicited interaction with
those enforcing the Code of Conduct, for a specified period of time. This
includes avoiding interactions in community spaces as well as external channels
like social media. Violating these terms may lead to a temporary or permanent
ban.

### 3. Temporary Ban

**Community Impact**: A serious violation of community standards, including
sustained inappropriate behavior.

**Consequence**: A temporary ban from any sort of interaction or public
communication with the community for a specified period of time. No public or
private interaction with the people involved, including unsolicited interaction
with those enforcing the Code of Conduct, is allowed during this period.
Violating these terms may lead to a permanent ban.

### 4. Permanent Ban

**Community Impact**: Demonstrating a pattern of violation of community
standards, including sustained inappropriate behavior, harassment of an
individual, or aggression toward or disparagement of classes of individuals.

**Consequence**: A permanent ban from any sort of public interaction within the
community.

## Attribution

This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 2.0,
available at https://www.contributor-covenant.org/version/2/0/
code_of_conduct.html.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
https://www.contributor-covenant.org/faq. Translations are available at https://
www.contributor-covenant.org/translations.
# Contributing to slopes

This outlines how to propose a change to slopes. 
For more detailed info about contributing to this, and other tidyverse packages, please see the
[**development contributing guide**](https://rstd.io/tidy-contrib). 

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("ITSLeeds/slopes", fork = TRUE)`.

*   Install all development dependences with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). Follow the style described in <https://style.tidyverse.org/news.html>.

### Code style

*   New code should follow the tidyverse [style guide](https://style.tidyverse.org). 
    You can use the [styler](https://CRAN.R-project.org/package=styler) package to apply these styles, but please don't restyle code that has nothing to do with your PR.  

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  

## Code of Conduct

Please note that the slopes project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# slopes package

<!-- badges: start -->
[![R-CMD-check](https://github.com/ropensci/slopes/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/slopes/actions)
[![Codecov test coverage](https://codecov.io/gh/ropensci/slopes/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/slopes?branch=master)
[![Status at rOpenSci Software Peer Review](https://badges.ropensci.org/420_status.svg)](https://github.com/ropensci/software-review/issues/420)
<!-- badges: end -->

The **slopes** R package calculates the slope (longitudinal steepness, also known as gradient) of roads, rivers and other linear (simple) features, based on two main inputs:

- [vector](https://geocompr.robinlovelace.net/spatial-class.html#vector-data) linestring geometries defined by classes in the [`sf`](https://r-spatial.github.io/sf/) package
- [raster](https://geocompr.robinlovelace.net/spatial-class.html#raster-data) objects with pixel values reporting average height, commonly known as digital elevation model (**DEM**) datasets, defined by classes in the [`raster`](https://cran.r-project.org/package=raster) or more recent [`terra`](https://rspatial.org/terra) packages

Data on slopes are useful in many fields of research, including [hydrology](https://en.wikipedia.org/wiki/Stream_gradient), natural hazards (including [flooding](https://www.humanitarianresponse.info/fr/operations/afghanistan/infographic/afg-river-gradient-and-flood-hazard) and [landslide risk management](https://assets.publishing.service.gov.uk/media/57a08d0740f0b652dd0016f4/R7815-ADD017_col.pdf)), recreational and competitive sports such as [cycling](http://theclimbingcyclist.com/gradients-and-cycling-an-introduction/), [hiking](https://trailism.com/trail-grades/), and [skiing](https://www.snowplaza.co.uk/blog/16682-skiing-steeps-what-does-gradient-mean-ski-piste/).
Slopes are also also important in some branches of [transport and emissions modelling](https://www.sciencedirect.com/science/article/pii/S2352146516302642) and [ecology](https://doi.org/10.1016/j.ncon.2016.10.001).
See the [`intro-to-slopes` vignette](https://ropensci.github.io/slopes/articles/intro-to-slopes.html) for details on fields using slope data and the need for this package.

This README covers installation and basic usage. For more information about slopes and how to use the package to calculate them, see the [get started](https://ropensci.github.io/slopes/) and the [introducion to slopes](https://ropensci.github.io/intro-to-slopes/) vignette.


## How it works

The package takes two main types of input data for slope calculation: 
- vector geographic objects representing **linear features**, and 
- **elevation values** from a digital elevation model representing a continuous terrain surface or which can be downloaded using functionality in the package

The package can be used with two sources of elevation data:
- openly available elevation data via an interface to the [ceramic package](https://github.com/hypertidy/ceramic), enabling estimation of hilliness for routes anywhere worldwide even when local elevation data is lacking. The package takes geographic lines objects and returns elevation data per vertex (providing the output as a 3D point geometry in the sf package by default) and per line feature (providing average gradient by default).
- an elevation model, available on your machine.

## Getting started

### Installation

<!-- You can install the released version of slopes from [CRAN](https://CRAN.R-project.org) with: -->

<!-- ``` r -->
<!-- install.packages("slopes") -->
<!-- ``` -->

Install the development version from [GitHub](https://github.com/) with:

```{r, eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/slopes")
```

#### Installation for DEM downloads

If you do not already have DEM data and want to make use of the package's ability to download them using the `ceramic` package, install the package with suggested dependencies, as follows:

```{r, eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/slopes", dependencies = "Suggests")
```

Furthermore, you will need to add a MapBox API key to be able to get DEM datasets, by signing up and registering for a key at https://account.mapbox.com/access-tokens/ and then following these steps:

```{r, eval=FALSE}
usethis::edit_r_environ()
MAPBOX_API_KEY=xxxxx # replace XXX with your api key
```

## Basic examples

Load the package in the usual way. We will also load the `sf` library:

```{r message=FALSE, warning=FALSE}
library(slopes)
library(sf)
```

The minimum input data requirement for using the package is an `sf` object containing LINESTRING geometries, as illustrated below (requires a MapBox API key):

```{r}
sf_linestring = lisbon_route # import or load a linestring object
```

```{r, eval=FALSE}
sf_linestring_xyz = elevation_add(sf_linestring)  # dem = NULL
#> Loading required namespace: ceramic
#> Preparing to download: 12 tiles at zoom = 12 from 
#> https://api.mapbox.com/v4/mapbox.terrain-rgb/
```

```{r, echo=FALSE}
# note: the following should be TRUE
# identical(sf_linestring_xyz, lisbon_route_xyz_mapbox)
sf_linestring_xyz = lisbon_route_xyz_mapbox
```

With the default argument `dem = NULL`, the function downloads the necessary elevation information from Mapbox.
You can also this use a local digital elevation model (`dem = ...`), as shown in the example below:

```{r}
sf_linestring_xyz_local = elevation_add(sf_linestring, dem = dem_lisbon_raster)
```

In both cases you can obtain the average gradient of the linestring with `slope_xyz()` and plot the elevation profile with `plot_slope()` as follows:

```{r elevationprofile}
slope_xyz(sf_linestring_xyz_local)
plot_slope(sf_linestring_xyz_local)
```

_See more functions in [Get started](https://ropensci.github.io/slopes/articles/slopes.html) vignette._

## See more in vignettes

-   [Get started](https://ropensci.github.io/slopes/articles/slopes.html)
-   [An introduction to slopes](https://ropensci.github.io/slopes/articles/intro-to-slopes.html)
-   [Reproducible example: gradients of a road network for a given city](https://ropensci.github.io/slopes/articles/roadnetworkcycling.html)
-   [Verification of slopes](https://ropensci.github.io/slopes/articles/verification.html)
-   [Benchmarking slopes calculation](https://ropensci.github.io/slopes/articles/benchmark.html)

## Code of Conduct

Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). 
By contributing to this project, you agree to abide by its terms.
---
title: 'slopes: An R Package to Calculate Slopes of Roads, Rivers and Trajectories'
tags:
  - R
  - topography
  - slopes
  - gradient
  - steepness
  - roads
  - elevation
authors:
  - name: Robin Lovelace
    orcid: 0000-0001-5679-6536
    affiliation: 1
  - name: Rosa Félix
    orcid: 0000-0002-5642-6006
    affiliation: 2
date: 27 September 2021
bibliography: paper.bib
# We'll later convert to a Markdown Document and remove paper.Rmd:
output: github_document
editor_options:
  markdown:
    wrap: sentence
---

# Summary

This package provides functions and example data to support research into the slope (also known as longitudinal gradient or steepness) of linear geographic entities such as roads [@ariza-lopez_dataset_2019] and rivers [@cohen_global_2018].
The package was initially developed to calculate the steepness of street segments but can be used to calculate steepness of any linear feature that can be represented as LINESTRING geometries in the 'sf' class system [@pebesma_simple_2018].
The package takes two main types of input data for slope calculation: vector geographic objects representing linear features, and raster geographic objects with elevation values (which can be downloaded using functionality in the package) representing a continuous terrain surface.
Where no raster object is provided the package attempts to download elevation data using the 'ceramic' package.

# Statement of need

Although there are several ways to name "slope", such as "steepness", "hilliness", "inclination", "aspect", "gradient", "declivity", the referred `slopes` in this package can be defined as the "longitudinal gradient" of linear geographic entities, as defined in the context of rivers by[@cohen_global_2018].

The package was initially developed to research road slopes to support evidence-based sustainable transport policies.
Accounting for gradient when planning for new cycling infrastructure and road space reallocation for walking and cycling can improve outcomes, for example by helping to identify routes that avoid steep hills.
The package can be used to calculate and visualise slopes of rivers and trajectories representing movement on roads of the type published as open data by @ariza-lopez_dataset_2019.

Data on slopes are useful in many fields of research, including [hydrology](https://en.wikipedia.org/wiki/Stream_gradient), natural hazards (including [flooding](https://www.humanitarianresponse.info/fr/operations/afghanistan/infographic/afg-river-gradient-and-flood-hazard) and [landslide risk management](https://assets.publishing.service.gov.uk/media/57a08d0740f0b652dd0016f4/R7815-ADD017_col.pdf)), recreational and competitive sports such as [cycling](http://theclimbingcyclist.com/gradients-and-cycling-an-introduction/), [hiking](https://trailism.com/trail-grades/), and [skiing](https://www.snowplaza.co.uk/blog/16682-skiing-steeps-what-does-gradient-mean-ski-piste/).
Slopes are also also important in some branches of [transport and emissions modelling](https://www.sciencedirect.com/science/article/pii/S2352146516302642) and [ecology](https://doi.org/10.1016/j.ncon.2016.10.001).
A growing number of people working with geospatial data require accurate estimates of gradient, including:

-   Transport planning practitioners who require accurate estimates of roadway gradient for estimating energy consumption, safety and mode shift potential in hilly cities (such as Lisbon, the case study city used in the examples in the documentation).
-   Vehicle routing software developers, who need to build systems are sensitive to going up or down steep hills (e.g. bicycles, trains, and large trucks), such as active travel planning, logistics, and emergency services.
-   Natural hazard researchers and risk assessors require estimates of linear gradient to inform safety and mitigation plans associated with project on hilly terrain.
-   Aquatic ecologists, flooding researchers and others, who could benefit from estimates of river gradient to support modelling of storm hydrographs

There likely other domains where slopes could be useful, such as agriculture, geology, and civil engineering.

An example of the demand for data provided by the package is a map showing gradients across Sao Paulo (Brazil, see image below) that has received more than 300 'likes' on Twitter and generated conversations: <https://twitter.com/DanielGuth/status/1347270685161304069>

![](https://camo.githubusercontent.com/30a3b814dd72aef5b51db635f2ab6e1b6b6c57b856d239822788967a4932d655/68747470733a2f2f7062732e7477696d672e636f6d2f6d656469612f45724a32647238574d414948774d6e3f666f726d61743d6a7067266e616d653d6c61726765){width="50%"}

# Usage and Key functions

Install the development version from [GitHub](https://github.com/) with:

```{r, eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/slopes")
```

### Installation for DEM downloads

If you do not already have DEM data and want to make use of the package's ability to download them using the `ceramic` package, install the package with suggested dependencies, as follows:

```{r, eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/slopes", dependencies = "Suggests")
```

Furthermore, you will need to add a MapBox API key to be able to get DEM datasets, by signing up and registering for a key at https://account.mapbox.com/access-tokens/ and then following these steps:

```{r, eval=FALSE}
usethis::edit_r_environ()
MAPBOX_API_KEY=xxxxx # replace XXX with your api key
```

The key functions in the package are `elevation_add()`, which adds a third 'Z' coordinate value for each vertex defining LINESTRING objects, and `slope_xyz()` which calculates slopes for each linear feature in a simple features object.

By default, the elevation of each vertex is estimated using [bilinear interpolation](https://en.wikipedia.org/wiki/Bilinear_interpolation) (`method = "bilinear"`) which calculates point height based on proximity to the centroids of surrounding cells.
The value of the `method` argument is passed to the `method` argument in [`raster::extract()`](https://rspatial.github.io/raster/reference/extract.html) or [`terra::extract()`](https://rspatial.github.io/terra/reference/extract.html) depending on the class of the input raster dataset.
See @kidner_what_1999 for descriptions of alternative elevation interpolation and extrapolation algorithms.

<!-- # Calculating slopes on regional transport networks -->

<!-- # Acknowledgements -->

# References
<!-- Brief: https://callforpapers.2021.foss4g.org/foss4g-2021-academic/ -->

## Slopes: a package for reproducible slope calculation, analysis and visualisation

Slopes are important for many purposes, including flood risk, agriculture, geology, and infrastructure constructions.
In transport planning, consideration of gradient is especially important for walking and cycling routes, which provided our initial motivation for developing this package.
Slopes can be calculated using proprietary products such as ArcMap but we wanted to calculate slopes using free and open source software to enable transparency, reproducibility and accessibility of our methods.
We developed the software in R because of prior experience with the language and the mature 'R-spatial' community which has developed mature codebases for working with geographic data in a reproducible command line environment, including `sf` (for working with vector datasets representing roads and other linear features) and `raster` (for representing digital elevation models, DEMs).

Building on these foundations the package is now working and has been used to calculate slopes on hundreds of roads in several cities.
Comparison with ArcMap's 3D analyst show that the approach is competitive with the go-to proprietary produce in terms of computational speed and that we can reproduce ArcMap's results: tests show an R-squared value of 0.99.
We hope the package will be of use and interest to the FOSS4G and in the talk will discuss ideas for taking the work forward, e.g. by implementing the logic into other languages/environments such as Rust, Python and even as a QGIS plugin.
Could there be scope for an inter-disciplinary and language-agnostic community interested in slope analysis?
We would like support efforts to strengthen links between geospatial developers who use R and the wider FOSS4G community, for example by comparing the slopes package with other open source approaches for slope calculation and analysis for mutually beneficial learning.
We will conclude with discussion of possible future directions of travel for the project, including possibilities for 3D visualisation, auto-download of elevation point data sampled across linear features (currently the package automates the download of DEM data but requires a MapBox API key) and using the slope values to generate evidence in support of sustainable transport policies. 

# Description

The package calculates longitudinal steepness of linear features such as roads and rivers, based on two main inputs: vector linestring geometries and raster digital elevation model (DEM) datasets.

After installing R, it can be installed as follows:

```{r, eval=FALSE}
remotes::install_github("ropensci/slopes")
```

The minimum data requirements for using the package are elevation points, either as a vector, a matrix or as a digital elevation model (DEM) encoded as a raster dataset. Typically you will also have a geographic object representing the roads or similar features. These two types of input data are represented in the code output and plot below.

# Notes


# Authors


Rosa Félix (1)

Robin Lovelace (2)

(1) University of Lisbon

(1) University of Leeds



---
title: "Debugging issues in the slopes package"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Debugging issues in the slopes package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(slopes)
```

The aim of this vignette is to document issues that have been identified and (hopefully, eventually) fixed.
It does not run by default.
To run it, change `FALSE` to `TRUE` in the code chunk below.

```{r}
knitr::opts_chunk$set(eval = FALSE)
```


## Error in weighted mean (#12)

This issue showed that the function `slope_raster()` errors in the final line:

```{r}
library(sf)
library(raster)
library(slopes)

#import datasets
network_original = st_read("https://github.com/U-Shift/Declives-RedeViaria/blob/main/shapefiles/RedeViariaLisboa_dadosabertos.gpkg?raw=true")
network = st_zm(network_original, drop = T) # make sure it has no Z values stored
network
u = "https://github.com/U-Shift/Declives-RedeViaria/blob/main/raster/LisboaNASA_clip.tif?raw=true"
download.file(u, "LisboaNASA_clip.tif")
dem = raster("LisboaNASA_clip.tif")
dem

# do they overlap?
raster::plot(dem)
plot(sf::st_geometry(network), add = TRUE)
# why this error?
network$slopeNASA = slope_raster(network, e = dem)
# In addition: Warning messages:
# 1: In max(d, na.rm = TRUE) :
#   no non-missing arguments to max; returning -Inf
# 2: In max(d, na.rm = TRUE) :
#   no non-missing arguments to max; returning -Inf
```

Does it error on a subset of the network?

```{r}
network_subset = network[1:9, ]
network$slopesNASA = slope_raster(network_subset, e = dem) # also fails with the same error
dem_subset = crop(dem, network)
network_subset$slopesNASA = slope_raster(network_subset, e = dem_subset) # also fails with the same error
unique(st_geometry_type(network_subset))
unique(st_geometry_type(lisbon_road_segment))
mapview::mapview(network_subset)
network_subset_ls = sf::st_cast(network_subset, to = "LINESTRING")
mapview::mapview(network_subset_ls)
network_subset_ls$slopesNASA = slope_raster(network_subset_ls, e = dem_subset) 
```

Testing on the full network:


```{r}
# install latest version
remotes::install_github("ropensci/slopes")
library(sf)
library(raster)
library(slopes)

#import datasets
network_original = st_read("https://github.com/U-Shift/Declives-RedeViaria/blob/main/shapefiles/RedeViariaLisboa_dadosabertos.gpkg?raw=true")
network = st_zm(network_original, drop = T) # make sure it has no Z values stored
network
u = "https://github.com/U-Shift/Declives-RedeViaria/blob/main/raster/LisboaNASA_clip.tif?raw=true"
download.file(u, "LisboaNASA_clip.tif")
dem = raster("LisboaNASA_clip.tif")
dem

# do they overlap?
raster::plot(dem)
plot(sf::st_geometry(network), add = TRUE)
# why this error?
network$slopeNASA = slope_raster(network, e = dem)
network_ls = sf::st_cast(network, "LINESTRING")
network_ls$slopeNASA = slope_raster(network_ls, e = dem)
plot(network_ls["slopeNASA"], lwd = 5)
```

# Issues with terra

```{r}
remotes::install_github("ropensci/slopes")
library(slopes)
class(lisbon_road_network)
raster::plot(dem_lisbon_raster)
plot(sf::st_geometry(lisbon_road_network), add = TRUE)
demterra = terra::rast(dem_lisbon_raster)
plot(sf::st_geometry(lisbon_road_network), add = TRUE)
lisbon_road_network$sloperaster = slope_raster(lisbon_road_network, e = dem_lisbon_raster)
library(terra)
plot(demterra)
lisbon_road_network$slopeterra = slope_raster(lisbon_road_network, e = demterra, terra = T)
summary(lisbon_road_network$sloperaster)
summary(lisbon_road_network$slopeterra)

# Looking at the code, this seems to be the culprit:

# res = as.numeric(terra::extract(e, m[, 1:2], method = method))
# vs
# res = as.numeric(raster::extract(e, m[, 1:2], method = method))
# lets try it...
method = "bilinear"
m = sf::st_coordinates(lisbon_road_segment)[1:2, ]
e = demterra
as.numeric(terra::extract(e, m[, 1:2], method = method))
as.numeric(terra::extract(dem_lisbon_raster, m[, 1:2], method = method))
res_terra = terra::extract(e, m[, 1:2], method = method)
class(res_terra)
res_terra
as.numeric(res_terra[, "r1"])
```


---
# title: "slopes: a package for calculating slopes<br>📈🎢🗻🛣️"
title: "slopes: a package for calculating slopes "
subtitle: "📈 of roads, rivers and other linear (simple) features 📉"  
author: "Robin Lovelace & Rosa Félix, UseR 2021<br><br><br><br><br><br><br><br><br><br>"
output:
  xaringan::moon_reader:
    css: xaringan-themer.css
    nature:
      slideNumberFormat: "%current%"
      highlightStyle: github
      highlightLines: true
      ratio: 16:9
      countIncrementalSlides: true
---

```{r, eval=FALSE, echo=FALSE}
# see slides manually uploaded online: https://slopes-slides.netlify.app/slides.html#1
# to run these slides locally:
xaringan::inf_mr("data-raw/slides.Rmd")
```

```{r xaringanExtra, echo=FALSE}
# From https://github.com/gadenbuie/xaringanExtra
xaringanExtra::use_xaringan_extra(c("tile_view", "animate_css", "tachyons"))
```


```{r setup, include=FALSE}
options(htmltools.dir.version = FALSE)
knitr::opts_chunk$set(
  fig.width=9, fig.height=3.5, fig.retina=3,
  out.width = "100%",
  # cache = TRUE,
  echo = TRUE,
  message = FALSE, 
  warning = FALSE,
  fig.show = TRUE,
  hiline = TRUE
)
```

```{r xaringan-themer, include=FALSE, warning=FALSE}
library(xaringanthemer)
style_duo_accent(
  title_slide_background_color = "#FFFFFF",
  title_slide_background_size = "100%",
  title_slide_background_image = "https://user-images.githubusercontent.com/1825120/121391204-04c75c80-c946-11eb-8d46-ab5d8ada55c2.png",
  title_slide_background_position = "bottom",
  title_slide_text_color = "#080808",
  primary_color = "#080808",
  secondary_color = "#FF961C",
  inverse_header_color = "#FFFFFF"
)
```

background-image: url(https://camo.githubusercontent.com/30a3b814dd72aef5b51db635f2ab6e1b6b6c57b856d239822788967a4932d655/68747470733a2f2f7062732e7477696d672e636f6d2f6d656469612f45724a32647238574d414948774d6e3f666f726d61743d6a7067266e616d653d6c61726765)
background-position: center
background-size: 100%

--
<br>

# Contents:

## Why slopes?

--

## Key functions

--

## Future plans

--

---

# Why we developed the slopes package

.left-column[

- Real world problems to solve involving slopes
- Existing tools were not up to the job
  - Expensive and hard to reproduce findings (ESRI's 3D analyst)
  - Hard to scale-up (online services)

- R programming challenge

- Support for route planning in active transportation

]

--

.right-column[

Real world problem: infrastructure prioritisation. Source: [paper](https://www.jtlu.org/index.php/jtlu/article/view/862) and www.pct.bike

![](https://user-images.githubusercontent.com/1825120/123165723-805bfa00-d46c-11eb-8446-969ea69e0287.png)

]

???

- Rosa motivation

---

# Applications

.left-column[

-   Transport planning
  - Active travel planning
  - Logistics/route planning
  - Emergency services
-   River/flooding research
-   Civil engineering

]

.right-column[

![](https://onlinelibrary.wiley.com/cms/asset/d3b27e5b-fb41-4169-94d3-633aa23a5ccc/gean12253-fig-0003-m.png)

Image source: Goodchild ([2020](https://doi.org/10.1111/gean.12253)): Beyond Tobler’s Hiking Function 

]

---

# Installation and set-up

```{r, eval=FALSE}
remotes::install_github("ropensci/slopes")
```


```{r}
library(slopes)
library(tmap)
tmap_mode("view")
```

---

#### How the package works

Key functions:

1. `slope_xyz()`: calculates the slope associated with linestrings that have xyz coordinates
1. `slope_raster()`: Calculate slopes of linestrings based on local raster map
1. `elevation_add()`: Adds a third dimension to linestring coordinates
1. `plot_slope()`: Plots the slope profile associated with a linestring
1. See https://ropensci.github.io/slopes/reference/index.html for more


```{r}
lisbon_route_3d_segments = stplanr::rnet_breakup_vertices(lisbon_route_3d)
lisbon_route_3d_segments$slope = slope_xyz(lisbon_route_3d_segments)
tm_shape(lisbon_route_3d_segments) +
  tm_lines(col = "slope", lwd = 3, palette = "viridis")
```


---

# Package data

The key input datasets are:

--

.pull-left[

linestrings representing roads/rivers/other, and...

```{r maps}
tm_shape(lisbon_road_network) + tm_lines()
```

]

--

.pull-right[

and digital elevations:

```{r}
tm_shape(dem_lisbon_raster) + tm_raster(palette = "BrBG", alpha = 0.3)
```

]


---

### Adding the Z dimension

```{r}
lisbon_route
```

#### `Dimension:     XY`

```{r}
lisbon_route_slopes = elevation_add(routes = lisbon_route, dem = slopes::dem_lisbon_raster)
```

```{r, eval=FALSE}
lisbon_route_slopes
## Simple feature collection with 1 feature and 3 fields
## Geometry type: LINESTRING
## Dimension:     XYZ
```


#### `Dimension:     XYZ`

---

## Plotting the Z dimension

.pull-left[

```{r, error=TRUE, fig.height=6}
plot_slope(lisbon_route)
```

]

.pull-right[


```{r, fig.height=6}
plot_slope(lisbon_route_slopes)
```

]

---

### Find slopes when you don't have a DEM

```{r, eval=FALSE}
usethis::edit_r_environ()
# Type in (register on the mapbox website):
MAPBOX_API_KEY=xxxxx
```

```{r}
library(stplanr)
origin = tmaptools::geocode_OSM("rail station zurich", as.sf = TRUE)
destination = tmaptools::geocode_OSM("eth zurich", as.sf = TRUE)
route = osrm::osrmRoute(src = origin, dst = destination, returnclass = "sf")
library(stplanr)
route = route(origin, destination, route_fun = cyclestreets::journey)
route_3d = elevation_add(route, dem = NULL)
```

```{r}
route_3d$gradient_slopes = slope_xyz(route_3d) # todo: calculate slopes in elevation_add by default?
```

---

.pull-left[

```{r}
library(tmap)
m = tm_shape(route_3d) + 
  tm_lines("gradient_slopes", lwd = 3,
           palette = "viridis")
tmap_mode("plot")
m
# Todo: add slope_map function with default palette?
```

]

.pull-right[

```{r}
tmap_mode("view")
m
```

]

---

.pull-left[

###  Worked example

See [vignette](https://ropensci.github.io/slopes/articles/roadnetworkcycling.html) 

Load packages

```{r}
# Get linear features you want the gradients of
library(slopes)
library(dplyr)
library(sf)
# remotes::install_github("ITSLeeds/osmextract")
library(osmextract) # see UseR talk on osmextract package
library(tmap)
network = oe_get("Isle of Wight", vectortranslate_options = c("-where", "highway IS NOT NULL")) 
```

]

.pull-right[

```{r, fig.height=9}
u = "https://github.com/U-Shift/Declives-RedeViaria/releases/download/0.2/IsleOfWightNASA_clip.tif"
f = basename(u) # Get digital elevation data
download.file(url = u, destfile = f, mode = "wb")
dem = raster::raster(f)
library(raster)
plot(dem)
plot(sf::st_geometry(network), add = TRUE) #check if they overlay
```

]

---

## Calculating slopes

```{r}
sys_time = system.time({
  network$slope = slope_raster(network, dem)
})
sys_time
nrow(network)
nrow(network) / sys_time[3]

network$slope = network$slope * 100 # percentage
summary(network$slope) # check the values
```

---

## Results

```{r, eval=FALSE}
qtm(network, "slope") # with a few extra arguments...
```

See [`roadnetworkcycling`](https://ropensci.github.io/slopes/articles/roadnetworkcycling.html) vignette for details and here for interactive map: http://web.tecnico.ulisboa.pt/~rosamfelix/gis/declives/SlopesIoW.html

![](https://user-images.githubusercontent.com/1825120/121820334-2435f080-cc8a-11eb-962c-79dcba97e459.png)

---

## A smaller example using packaged data

```{r}
routes = lisbon_road_network
dem = dem_lisbon_raster
routes$slope = slope_raster(routes, dem)
plot(dem)
plot(routes["slope"], add = TRUE)
```

---

# Future plans

.pull-left[

-   `slope_get_dem()` to get digital elevation model data
-   Add other elevation sources using api keys (e.g Google, others)?
-   Improve `plot_slope()` visualization
-   Explore accuracy of data vs 'ground truth'
-   Finish review, publish in JOSS and on CRAN
-   See [ropensci/slopes](https://github.com/ITSLeeds/slopes) on github to get involved!

]

--

.pull-right[

![](https://www.thestar.co.uk/jp-ct.co.uk/image/onecms:9b91244b-bef9-44d5-af34-56943ea0a809:1bd497fc-203b-4814-8905-1268a2c6c2a5/NSST-16-01-20-blake-NMSYupload.jpeg?&width=1200)

]


???

invite people to help in github?

--

Blake Street, Sheffield. Source: [thestar.co.uk](https://www.thestar.co.uk/news/people/which-sheffields-steepest-street-these-are-ones-you-think-are-worst-1368220)

---


# Where is the ground truth? 

.pull-left[

- Comparison of gradients from the slopes packages and an online routing service

```{r, echo=FALSE}
# route_3d$elevation_change
# route_3d$distances
plot(route_3d$gradient_smooth, route_3d$gradient_slopes)
```

]

.pull-right[

- How to get ground truth data? Quite hard!

![](https://media.springernature.com/full/springer-static/image/art%3A10.1038%2Fs41597-019-0147-x/MediaObjects/41597_2019_147_Fig5_HTML.png)

]

---

## Estimating slopes of bridges

```{r, out.width="80%", fig.show='hold', echo=FALSE}
knitr::include_graphics(c(
  "slope-edinburgh-bridge.png"
  # "slope-edinburgh-bridge2.png"
))
```


???

- RL

- My research into tools for prioritising cycling investment

- UK not very hilly but models fail slightly in hilly areas

---

```{r, out.width="80%", fig.show='hold', echo=FALSE}
knitr::include_graphics(c(
  # "slope-edinburgh-bridge.png",
  "slope-edinburgh-bridge2.png"
))
```



# Thanks!

Slides created via the R packages:

[**xaringan**](https://github.com/yihui/xaringan)<br>
[gadenbuie/xaringanthemer](https://github.com/gadenbuie/xaringanthemer)
(And the [Sharing Xaringan Slides](https://www.garrickadenbuie.com/blog/sharing-xaringan-slides/) blog post!)

The chakra comes from [remark.js](https://remarkjs.com), [**knitr**](http://yihui.name/knitr), and [R Markdown](https://rmarkdown.rstudio.com).

---
title: "Verification of slopes"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Verification of slopes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: slope-references.bib
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# Introduction

This article aims to provide external verification of the results provided by the package.
So far only one verification dataset has been used, but we hope to find others.
If you know of verification datasets, please let us know --- initially we planned to use a dataset from a paper on river slopes [@cohen_global_2018], but we could find no way of extracting the underlying data to do the calculation.

For this article we primarily used the following packages, although others are loaded in subsequent code chunks.

```{r setup}
library(slopes)
library(sf)
```

The results are reproducible (requires downloading input data manually and installing additional packages).
To keep package build times low, only the results are presented below.

# Comparison with results from ArcMap 3D Analyst




# Three-dimensional traces of roads dataset

<!-- todo: make the segments -->

An input dataset, comprising a 3D linestring recorded using a dual frequency GNSS receiver (a [Leica 1200](https://gef.nerc.ac.uk/equipment/gnss.php)) with a vertical accuracy of 20 mm
<!-- 138 GPS 3D traces of a hilly road from a peer reviewed journal article -->
[@ariza-lopez_dataset_2019] was downloaded from the 
<!-- [figshare website as a .zip file](https://ndownloader.figshare.com/files/14331197) - raw data -->
[figshare website as a .zip file](https://ndownloader.figshare.com/files/14331185)
and unzipped and inflated in the working directory as follows (not evaluated to reduce package build times): 

```{r, eval=FALSE}
download.file("https://ndownloader.figshare.com/files/14331185", "3DGRT_AXIS_EPSG25830_v2.zip")
unzip("3DGRT_AXIS_EPSG25830_v2.zip")
trace = sf::read_sf("3DGRT_AXIS_EPSG25830_v2.shp")
plot(trace)
nrow(trace)
#> 11304
summary(trace$X3DGRT_h)
#>  Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
#>   642.9   690.3   751.4   759.9   834.3   884.9 
```


```{r, eval=FALSE, echo=FALSE}
# original trace dataset
traces = sf::read_sf("vignettes/3DGRT_TRACES_EPSG25830_v2.shp")
traces = sf::read_sf("3DGRT_TRACES_EPSG25830_v2.shp")
nrow(traces)
#> [1] 111113
```

To verify our estimates of hilliness, we generated slope estimates for each segment and compared them with [Table 7](https://www.nature.com/articles/s41597-019-0147-x/tables/7) in @ariza-lopez_dataset_2019.
The absolute gradient measure published in that paper were:

```{r}
res_gps = c(0.00, 4.58, 1136.36, 6.97)
res_final = c(0.00, 4.96, 40.70, 3.41)
res = data.frame(cbind(
  c("GPS", "Dual frequency GNSS receiver"),
  rbind(res_gps, res_final)
))
names(res) = c("Source", "min", " mean", " max", " stdev")
knitr::kable(res, row.names = FALSE)
```

```{r, eval=FALSE}
# mapview::mapview(trace) # check extent: it's above 6km in height
# remotes::install_github("hypertidy/ceramic")
loc = colMeans(sf::st_coordinates(sf::st_transform(trace, 4326)))
e = ceramic::cc_elevation(loc = loc[1:2], buffer = 3000)
trace_projected = sf::st_transform(trace, 3857)
plot(e)
plot(trace_projected$geometry, add = TRUE)
```

```{r, echo=FALSE, eval=FALSE}
# aim: get max distance from centrepoint
bb = sf::st_bbox(sf::st_transform(trace, 4326))
geosphere::distHaversine(c(bb[1], bb[2]), c(bb[3], bb[2]))
geosphere::distHaversine(c(bb[1], bb[2]), c(bb[1], bb[4]))
# max of those 2 and divide by 2
```


```{r, echo=FALSE}
knitr::include_graphics("https://user-images.githubusercontent.com/1825120/81125221-75c06780-8f2f-11ea-8cea-ad6322ef99e7.png")
```

The slopes were estimated as follows:

```{r, eval=FALSE}
# source: https://www.robinlovelace.net/presentations/munster.html#31
points2line_trajectory = function(p) {
  c = st_coordinates(p)
  i = seq(nrow(p) - 2)
  l = purrr::map(i, ~ sf::st_linestring(c[.x:(.x + 1), ]))
  lfc = sf::st_sfc(l)
  a = seq(length(lfc)) + 1 # sequence to subset
  p_data = cbind(sf::st_set_geometry(p[a, ], NULL))
  sf::st_sf(p_data, geometry = lfc)
}
r = points2line_trajectory(trace_projected)
# summary(st_length(r)) # mean distance is 1m! Doesn't make sense, need to create segments
s = slope_raster(r, e = e)
slope_summary = data.frame(min = min(s), mean = mean(s), max = max(s), stdev = sd(s))
slope_summary = slope_summary * 100
knitr::kable(slope_summary, digits = 1)
```

| min| mean|  max| stdev|
|---:|----:|----:|-----:|
|   0|  6.2| 48.2|   5.6|

Combined with the previous table from @ariza-lopez_dataset_2019, these results can be compared with those obtained from mainstream GPS, and an accurate GNSS receiver:

|Source                       |min | mean | max    | stdev |
|:----------------------------|:---|:-----|:-------|:------|
|GPS                          |0   |4.58  |1136.36 |6.97   |
|Dual frequency GNSS receiver |0   |4.96  |40.7    |3.41   |
|Slopes R package             |0   |6.2   |48.2    |5.6    |

It is notable that the package substantially overestimates the gradient, perhaps due to the low resolution of the underlying elevation raster.
However, the slopes package seems to provide less noisy slope estimates than the GPS approach, with lower maximum values and low standard deviation.

# References




```{r, eval=FALSE, echo=FALSE}
# failed tests
raster::extract(e, trace_projected)
raster::writeRaster(e, "e.tif")
e_terra = terra::rast("e.tif")
terra::crs(e_terra)
v = terra::vect("vignettes/3DGRT_TRACES_EPSG25830_v2.shp")
e_wgs = terra::project(e_terra, v)
e_stars = stars::st_as_stars(e)
e_wgs = sf::st_transform(e_stars, 4326)
stars::write_stars(e_wgs, "e_wgs.tif")
e2 = raster::raster("e_wgs.tif")
raster::plot(e)
plot(trace$geometry, add = TRUE)
```


```{r, echo=FALSE, eval=FALSE}
# discarded way:
remotes::install_github("jhollist/elevatr")
sp_bbox = sp::bbox(sf::as_Spatial(sf::st_transform(trace, 4326)))
e = elevatr::get_aws_terrain(locations = sp_bbox, prj = "+init:4326")
```

---
title: "Get started"
output: bookdown::html_vignette2
vignette: >
  %\VignetteIndexEntry{Get started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

Welcome to the slopes vignette, a type of long-form documentation/article that introduces the core functions and functionality of the `slopes` package.

# Installation

You can install the released version of slopes from [CRAN](https://CRAN.R-project.org) with:

```r
install.packages("slopes")
```

Install the development version from [GitHub](https://github.com/) with:

```{r, eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/slopes")
```

### Installation for DEM downloads

If you do not already have DEM data and want to make use of the package's ability to download them using the `ceramic` package, install the package with suggested dependencies, as follows:

```{r, eval=FALSE}
# install.packages("remotes")
remotes::install_github("ropensci/slopes", dependencies = "Suggests")
```

Furthermore, you will need to add a MapBox API key to be able to get DEM datasets, by signing up and registering for a key at https://account.mapbox.com/access-tokens/ and then following these steps:

```{r, eval=FALSE}
usethis::edit_r_environ()
MAPBOX_API_KEY=xxxxx # replace XXX with your api key
```


# Functions 

## Elevation

- `elevation_add()` Take a linestring and add a third dimension (z) to its coordinates
- `elevation_get()` Get elevation data from hosted maptile services (returns a raster)
- `elevation_extract()` Extract elevations from coordinates

- `z_value()` retrieves elevation values for each z (as vector of sequential vertices)
- `z_start()` retrieves the elevation value of the first linestring vertice
- `z_end()` retrieves the elevation value of the last linestring vertice
- `z_mean()` retrieves the elevation mean value 
- `z_max()`  retrieves the elevation max value 
- `z_min()`  retrieves the elevation min value 

## Distance

- `sequential_dist()` Calculate the sequential distances between sequential coordinate pairs

## Slope

- `slope_vector()` calculates the slopes associated with consecutive elements in one dimensional distance and associated elevations.
- `slope_distance()` calculates the slopes associated with consecutive distances and elevations.
- `slope_distance_mean()` calculates the mean average slopes associated with consecutive distances and elevations.
- `slope_distance_weighted()` calculates the slopes associated with consecutive distances and elevations, with the mean value associated with each set of distance/elevation vectors weighted in proportion to the distance between each elevation measurement, so longer sections have proportionally more influence on the resulting gradient estimate.

- `slope_raster()` Calculate slopes of linestrings based on local raster map
- `slope_matrix()` Calculate the gradient of line segments from a 3D matrix of coordinates
- `slope_matrix_weighted()` Calculate the weighted gradient of line segments from a 3D matrix of coordinates

- `slope_xyz()` Calculates the slope associated with linestrings that have xyz coordinates

## Plot

- `plot_dz()` Plot a digital elevation profile based on xyz data
- `plot_slope()` Plots the slope profile associated with a linestring with base R graphics

# Package datasets

The `slopes` package comes with some datasets to play with:

**Linestrings:**

- `lisbon_road_segment`: a single road segment in Lisbon (XY)
- `lisbon_route`: a route with some variation in elevation in Lisbon (XY)
- `cyclestreets_route`: a bike route in Leeds (XY)

**Road network:**

- `lisbon_road_network`: a sample of road segments in downtown Lisbon
- `magnolia_xy`: a sample of road segments in center Seattle, in the Magnolia neighborhood

**Digital elevation model (DEM):**

- `dem_lisbon_raster` a DEM of downtown Lisbon (EPSG:3763)


# Examples

Load the package in the usual way. We will also load the `sf` library:

```{r message=FALSE, warning=FALSE}
library(slopes)
library(sf)
```

The minimum input data requirement for using the package is an `sf` object containing LINESTRING geometries.  

You can also create `sf` objects from a matrix of coordinates, as illustrated below (don't worry about the details for now, you can read up on how all this works in the `sf` package [documentation](https://r-spatial.github.io/sf/articles/sf1.html)):

```{r, eval=FALSE, echo=FALSE}
m = st_coordinates(sf::st_transform(lisbon_road_segment, 4326))
s = seq(from = 1, to = nrow(m), length.out = 4)
round(m[s, 1:2], 5)
dput(round(m[s, 1], 4))
dput(round(m[s, 2], 4))
```

```{r}
m = cbind(
  c(-9.1333, -9.134, -9.13),
  c(38.714, 38.712, 38.710)
)
sf_linestring = sf::st_sf(
  data.frame(id = 1),
  geometry = st_sfc(st_linestring(m)),
  crs = 4326
)
class(sf_linestring)
st_geometry_type(sf_linestring)
```

> maybe remove this? or add step 1 and step 2 again.

## Single road segment + no DEM

You can check your input dataset is suitable with the functions `class()` from base R and `st_geometry_type()` from the `sf` package, as demonstrated below on the example object `lisbon_road_segment` that is contained within the package:

```{r}
sf_linestring = lisbon_road_segment
class(sf_linestring)
st_geometry_type(sf_linestring)
```


A quick way of testing if your object can have slopes calculated for it is to plot it in an interactive map and to check that underneath the object there is indeed terrain that will give the linestrings gradient:

```{r linestringmap, message=FALSE, warning=FALSE}
library(tmap)
tmap_mode("view")
tm_shape(sf_linestring) +
  tm_lines(lwd = 5) +
  tm_basemap(leaflet::providers$Esri.WorldTopoMap)
```

Imagine you want to calculate the gradient of the route shown above. 
You can do this as a two step process as follows.

**Step 1**: add elevations to each coordinate in the linestring (requires a [MapBox API](https://account.mapbox.com/access-tokens/) key):

```{r, eval=FALSE}
sf_linestring_xyz = elevation_add(sf_linestring) # dem = NULL
#> Loading required namespace: ceramic
#> Preparing to download: 9 tiles at zoom = 18 from 
#> https://api.mapbox.com/v4/mapbox.terrain-rgb/
```

```{r, echo=FALSE}
# note: the following should be TRUE
# identical(sf_linestring_xyz, lisbon_road_segment_xyz_mapbox)
sf_linestring_xyz = lisbon_road_segment_xyz_mapbox
```

With the argument `dem = NULL`, the function downloads the necessary elevation information from Mapbox. You can use this argument with a local digital elevation model (`dem = ...`).

You can check the elevations added to the new `sf_linestring_xyz` object by printing its coordinates, as follows (note the new Z column that goes from above 87 m above sea level to only 79 m in a short distance).

```{r}
st_coordinates(sf_linestring_xyz)
```

You can use the `z_` functions to extract such values:

```{r}
z_value(sf_linestring_xyz) # returns all the elevation values between xy coordinates

z_mean(sf_linestring_xyz) # elevation mean value
z_min(sf_linestring_xyz) # elevation min value 
z_max(sf_linestring_xyz) # elevation max value 
z_start(sf_linestring_xyz) # first z
z_end(sf_linestring_xyz) # last z
```

**Step 2**: calculate the average slope of the linestring

```{r}
slope_xyz(sf_linestring_xyz)
```

The result, just over 0.2, tells us that it's quite a steep slope: a 21% gradient *on average*.

## Route + available DEM

Using the slopes package we can estimate the gradient of individual road segments. When these segments are combined into routes, we then need a means of assessing the hilliness of the entire route. A range of indices can be used to represent route hilliness. The choice of which index is most appropriate may be context dependent (see the [introducion to slopes](https://ropensci.github.io/intro-to-slopes/) vignette).
  
Again, let us use the same function with a entire route, `lisbon_route`, also available in the package:

```{r message=FALSE, warning=FALSE}
sf_route = lisbon_route
class(sf_route)
st_geometry_type(sf_route)

tm_shape(sf_route) +
  tm_lines(lwd = 3) +
  tm_basemap(leaflet::providers$Esri.WorldTopoMap)
```

**Step 1**: add elevations to each coordinate in the route:

```{r, eval=FALSE}
sf_route_xyz = elevation_add(sf_route)
#> Loading required namespace: ceramic
#> Preparing to download: 12 tiles at zoom = 15 from 
#> https://api.mapbox.com/v4/mapbox.terrain-rgb/
```

```{r, echo=FALSE}
# note: the following should be TRUE
# identical(sf_route_xyz, lisbon_road_segment_xyz_mapbox)
sf_route_xyz = lisbon_route_xyz_mapbox
```

**Step 2**: calculate the average slope of the route

```{r}
slope_xyz(sf_route_xyz)
```

The result shows a 7.7% gradient *on average*.

Now, if you already have a DEM, you can calculate the slopes directly as follows, with `slope_raster()`:

```{r}
class(dem_lisbon_raster)
slope_raster(routes = sf_route,
             dem = dem_lisbon_raster)
```

The result shows a 7.8% gradient *on average*.
As you can see, the retrieved result from elevation information available in Mapbox and in this Digital Elevation Model, is quite similar. (See more about these differences in [Verification of slopes](https://ropensci.github.io/slopes/articles/verification.html).) 

## Route with xyz coordinates

If your linestring object already has X, Y and Z coordinates (e.g. from a GPS device), you can use the `slope_` functions directly. 

```{r eval=FALSE, include=FALSE}
#not to use like this... it would ge good to have a gps example to demonstrate

slope_vector(sf_route_xyz)
slope_distance(sf_route_xyz)
slope_distance_mean(sf_route_xyz)
slope_distance_weighted(sf_route_xyz)

slope_vector(sf_linestring_xyz)
slope_distance(sf_linestring_xyz)
slope_distance_mean(sf_linestring_xyz)
slope_distance_weighted(sf_linestring_xyz)
```

```{r}
# for a line xz
x = c(0, 2, 3, 4, 5, 9)
elevations = c(1, 2, 2, 4, 3, 1) / 10
slope_vector(x, elevations)

# for a path xyz
xy = st_coordinates(sf_linestring)
dist = sequential_dist(xy, lonlat = FALSE)
elevations = elevation_extract(xy, dem_lisbon_raster)

slope_distance(dist, elevations)
slope_distance_mean(dist, elevations)
slope_distance_weighted(dist, elevations)
```

In any case, to use the `slopes` package you need **elevation points**, either as a vector, a matrix or as a digital elevation model (DEM) encoded as a raster dataset.

# Calculating and plotting gradients

## Road network

Typical use cases for the package are calculating the slopes of geographic objects representing roads or other linear features.
These two types of input data are represented in the code output and plot below.

```{r dem-lisbon}
# A raster dataset included in the package:
class(dem_lisbon_raster) # digital elevation model
summary(raster::values(dem_lisbon_raster)) # heights range from 0 to ~100m
raster::plot(dem_lisbon_raster)

# A vector dataset included in the package:
class(lisbon_road_network)
plot(sf::st_geometry(lisbon_road_network), add = TRUE)
```

Calculate the average gradient of **each road segment** as follows:

```{r}
lisbon_road_network$slope = slope_raster(lisbon_road_network, dem = dem_lisbon_raster)
summary(lisbon_road_network$slope)
```

This created a new column, `slope` that represents the average, distance weighted slope associated with each road segment.
The units represent the percentage incline, that is the change in elevation divided by distance.
The summary of the result tells us that the average gradient of slopes in the example data is just over 5%.
  
This result is equivalent to that returned by ESRI's `Slope_3d()` in the [3D Analyst extension](https://desktop.arcgis.com/en/arcmap/10.3/tools/3d-analyst-toolbox/slope.htm), with a correlation between the ArcMap implementation and our implementation of more than 0.95 on our test dataset (we find higher correlations on larger datasets - see the [verification of slopes](https://ropensci.github.io/slopes/articles/verification.html article):

```{r}
cor(
  lisbon_road_network$slope,    # slopes calculates by the slopes package
  lisbon_road_network$Avg_Slope # slopes calculated by ArcMap's 3D Analyst extension
)
```

We can now visualise the average slopes of each route calculated by the `slopes` package as follows:

```{r slope-vis}
raster::plot(dem_lisbon_raster)
plot(lisbon_road_network["slope"], add = TRUE, lwd = 5)
```

## Elevation profile

Taking the [first route example](#route--available-dem), imagine that we want to go from from the Santa Catarina area in the East of the map to the Castelo de São Jorge in the West.
This route goes down a valley and up the other side:

```{r route}
# library(tmap)
# tmap_mode("view")
qtm(lisbon_route)
```

```{r, echo=FALSE, eval=FALSE}
# Removed because it's not rendering in RMarkdown
mapview::mapview(lisbon_road_network["slope"], map.types = "Esri.WorldStreetMap")
mapview::mapview(lisbon_route)
```

We can convert the `lisbon_route` object into a 3d linestring object with X, Y and Z coordinates, using the elevation values stored in the DEM, as follows:

```{r, eval=TRUE}
lisbon_route_xyz = elevation_add(lisbon_route, dem_lisbon_raster) 
```

We can now visualise the elevation profile of the route as follows:

```{r plot_slope}
plot_slope(lisbon_route_xyz)
```

## Splitting the network

The `lisbon_route_xyz` example is useful but often you will want to calculate the slopes not of an entire route (in this case one that is 2.5 km long) but of segments.
There are various ways to split segements, including using algorithms from other packages or [GIS programs](https://github.com/paleolimbot/qgisprocess/issues/26), but here we'll use the `stplanr` function `rnet_breakup_vertices()` (see  [`vignette("roadnetworkcycling")`](https://ropensci.github.io/slopes/articles/roadnetworkcycling.html) for an example of this function working on a large road network):

```{r}
sf::st_length(lisbon_route_xyz) # check route length: 2.5 km
lisbon_route_segments = stplanr::rnet_breakup_vertices(lisbon_route_xyz)
summary(sf::st_length(lisbon_route_segments)) # mean of 50 m
```

We can now calculate the slope for each of these segments.

```{r}
lisbon_route_segments$slope = slope_xyz(lisbon_route_segments)
summary(lisbon_route_segments$slope)
```

## Directed slopes

The route has a direction that is implicit in the order of the vertices and segments.
From the perspective of someone travelling along the route, the slopes have a direction which is important: it's easier to go uphill than downhill.
To calculate the slopes with direction, add the `directed` argument as follows.

```{r}
lisbon_route_segments$slope_directed = slope_xyz(lisbon_route_segments, directed = TRUE)
summary(lisbon_route_segments$slope_directed)
```

Plotting the directed and undirected slopes side-by-side shows the importance of considering slope direction for route planning, which may want to avoid steep hills going uphill but not downhill for certain types of travel, for example.

```{r, fig.show='hold', out.width="50%"}
breaks = c(0, 3, 5, 8, 10, 20, 50)
breaks_proportion = breaks / 100
breaks_directed = c(-rev(breaks_proportion), (breaks_proportion[-1]))
plot(lisbon_route_segments["slope"], breaks = breaks_proportion)
plot(lisbon_route_segments["slope_directed"], breaks = breaks_directed)
```

```{r, eval=FALSE, echo=FALSE}
# test code
z = sf::st_make_grid(lisbon_route_xyz, cellsize = 100)
sampled_points = sf::st_line_sample(lisbon_route_xyz, n = 30)
points_sf = sf::st_sf(geometry = sf::st_cast(sampled_points, "POINT"))
plot(points_sf)
lisbon_route_segments = stplanr::route_split(lisbon_route_xyz, p = points_sf[2:3, ])
lisbon_route_segments = stplanr::rnet_breakup_vertices(lisbon_route_xyz)
library(tmap)
tmap_mode("view")
qtm(lisbon_route_segments, lines.lwd = 9, lines.col = 1:nrow(lisbon_route_segments))
plot(lisbon_route_segments, col = 1:nrow(lisbon_route_segments))
```


```{r, eval=FALSE, echo=FALSE}
# Test: try using QGIS
remotes::install_github("paleolimbot/qgisprocess")
library(qgisprocess)
qgis_configure()
algorithms = qgis_algorithms()
View(algorithms)
result = qgis_run_algorithm(
  algorithm = "grass7:v.split",
  INPUT = lisbon_route_xyz,
  LENGTH = 500
  )
route_segments = sf::st_read(result$OUTPUT)
route_segments
plot(lisbon_route_xyz$geometry)
plot(route_segments$geom, add = T, lwd = 3)
mapview::mapview(route_segments)
```


## Using `elevation_add()` with and without a `dem =` argument

If you do not have a raster dataset representing elevations, you can automatically download them by omitting the argument `dem = NULL` (a step that is automatically done in the function `elevation_add()` shown in the basic example above, results of the subsequent code chunk not shown):

```{r, message=FALSE, warning=FALSE, eval=FALSE}
dem_mapbox = elevation_get(lisbon_route)
lisbon_road_proj = st_transform(lisbon_route, raster::crs(dem_mapbox))
lisbon_route_xyz_mapbox = elevation_add(lisbon_road_proj, dem = dem_mapbox)
plot_slope(lisbon_route_xyz_mapbox)
```

As outlined in the basic example above this can be done more concisely, as:

```{r, eval=FALSE}
lisbon_route_xyz_auto = elevation_add(lisbon_route) #dem = NULL
```
```{r, echo=FALSE}
lisbon_route_xyz_auto = lisbon_route_xyz_mapbox
```
```{r}
plot_slope(lisbon_route_xyz_auto)
```

Note that the elevations shown in both plots differ, since the first is based on DEM elevation available, and the second is based in _Mapbox_ elevation.

# Commulative elevation change

The following example calculate the elevations of a route in Leeds, and plots its commutative sum along the route (not evaluated).

```{r, eval=FALSE}
cyclestreets_xyz = elevation_add(cyclestreets_route) 
plot_slope(cyclestreets_xyz)
plot(cumsum(cyclestreets_xyz$distances), cumsum(cyclestreets_xyz$elevation_change))
```


# See more in vignettes

-   [slopes package](https://ropensci.github.io/slopes/index.html)
-   [An introduction to slopes](https://ropensci.github.io/slopes/articles/intro-to-slopes.html)
-   [Reproducible example: gradients of a road network for a given city](https://ropensci.github.io/slopes/articles/roadnetworkcycling.html)
-   [Verification of slopes](https://ropensci.github.io/slopes/articles/verification.html)
-   [Benchmarking slopes calculation](https://ropensci.github.io/slopes/articles/benchmark.html)


---
title: "Benchmarking slopes calculation"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Benchmarking slopes calculation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(slopes)
library(bench)
library(raster)
```

# Performance

A benchmark can reveal how many route gradients can be calculated per second:

```{r, results='hide'}
e = dem_lisbon_raster
r = lisbon_road_network
et = terra::rast(e)
res = bench::mark(check = FALSE,
  slope_raster = slope_raster(r, e),
  slope_terra = slope_raster(r, et)
)
```

```{r}
res
```


That is approximately

```{r}
round(res$`itr/sec` * nrow(r))
```

routes per second using the `raster` and `terra` (the default if installed, using `RasterLayer` and native `SpatRaster` objects) packages to extract elevation estimates from the raster datasets, respectively.

The message: use the `terra` package to read-in DEM data for slope extraction if speed is important.

To go faster, you can chose the `simple` method to gain some speed at the expense of accuracy:

```{r, results='hide'}
e = dem_lisbon_raster
r = lisbon_road_network
res = bench::mark(check = FALSE,
  bilinear1 = slope_raster(r, e),
  bilinear2 = slope_raster(r, et),
  simple1 = slope_raster(r, e, method = "simple"),
  simple2 = slope_raster(r, et, method = "simple")
)
```

```{r}
res
```

```{r}
round(res$`itr/sec` * nrow(r))
```
---
title: "Example: gradients of a road network for a given city"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Example: gradients of a road network for a given city}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width=7, fig.height=5, fig.align = "center",
  eval = FALSE
)
```

An example of the demand for data provided by the `slopes` package is a map showing gradients across São Paulo (Brazil, see image below), with a simplistic classification for cycling difficulty.

![using slopes() to create a road network gradient for cycling for São Paulo (Brazil)](https://pbs.twimg.com/media/ErJ2dr8WMAIHwMn?format=jpg&name=small)

This vignette will guide through the production of an interactive slope map for a road network, using `slopes`, `osmextract`, `sf`, `stplanr` and `tmap`.  

For the convenience of sample, we will use [Isle of Wight](https://en.wikipedia.org/wiki/Isle_of_Wight) case, with 384 km^2^. See [Other examples] below.

This will follow three steps:

1.  Download of road network from [OpenStreetMap](https://www.openstreetmap.org/)
2.  Prepare the network
3.  Compute slopes and export the map in html


## Extract the OSM network from geofabrik

For this step you may use `osmextract` [package](https://ropensci.github.io/osmextract/articles/osmextract.html) which downloads the most recent information available at OSM (https://download.geofabrik.de/index.html) and converts to _GeoPackage_ (.gpkg), the equivalent to _shapefile_.

```{r setup1, warning=FALSE, message = FALSE}
library(dplyr)
library(sf)
# remotes::install_github("ITSLeeds/osmextract")
library(osmextract)
```

```{r get_iow, warning=FALSE, message = FALSE}
# get the network
iow_osm = oe_get("Isle of Wight", provider = "geofabrik", stringsAsFactors = FALSE, 
                 quiet = FALSE, force_download = TRUE, force_vectortranslate = TRUE) # 7 MB

# filter the major roads
iow_network = iow_osm %>% 
  dplyr::filter(highway %in% c('primary', "primary_link", 'secondary',"secondary_link", 
                               'tertiary', "tertiary_link", "trunk", "trunk_link", 
                               "residential", "cycleway", "living_street", "unclassified", 
                               "motorway", "motorway_link", "pedestrian", "steps", "track")) #remove: "service"
```

## Clean the road network

These are optional steps that give better results, although they may slow down the process since they increase the number of segments present in the network.

### Filter the unconnected segments

The [`rnet_group()`](https://docs.ropensci.org/stplanr/reference/rnet_group.html) function from `stplanar` package assesses the connectivity of each segment assigns a group number (similar to a clustering process). Then we may filter the main group, the one with more connected segments.

```{r setup2, warning=FALSE, message = FALSE}
# remotes::install_github("ropensci/stplanr")
library(stplanr)
```

```{r filter}
# filter unconnected roads
iow_network$group = rnet_group(iow_network)
iow_network_clean = iow_network %>% filter(group == 1) # the network with more connected segments
```

### Break the segments on vertices

A very long segment will have an assigned average slope, but a very long segment can be broken into its nodes and have its own slope in each part of the segment. On one hand, we want the segments to break at their nodes. On the other hand, we don't want artificial *nodes* to be created where two lines cross, in particular where they have different **z** levels (_eg._ *brunels*: bridges and tunnels).  

The [`rnet_breakup_vertices`](https://docs.ropensci.org/stplanr/reference/rnet_breakup_vertices.html) from `stplanr` breaks the segments at their inner vertices, preserving the **brunels**.

```{r breaking, warning=FALSE, message = FALSE}
iow_network_segments = rnet_breakup_vertices(iow_network_clean)
```

In this case, there are around 1.6 x the number of segments than in the original network. 

<!-- `r # round(nrow(iow_network_segments)/nrow(iow_network_clean),2)`x segments than the original network. -->

## Get slope values for each segment

For this case we will use `slope_raster()` [function](https://ropensci.github.io/slopes/reference/slope_raster.html), to retrieve the z values from a digital elevation model. This raster was obtained from STRM NASA mission.

The **SRTM** (*Shuttle Radar Topography Mission*) NASA’s mission provides [freely available](https://gisgeography.com/srtm-shuttle-radar-topography-mission/) worldwide DEM, with a resolution of 25 to 30m and with a vertical accuracy of 16m - [more](https://www2.jpl.nasa.gov/srtm/). The resolution for USA might be better.

Alternatively, **COPERNICUS** ESA's mission also provides [freely available](https://land.copernicus.eu/imagery-in-situ/eu-dem) DEM for all Europe, with a 25m resolution and a vertical accuracy of 7m - [more](https://land.copernicus.eu/user-corner/publications/eu-dem-flyer/view).

Depending of how large is your road network, you can use `elevation_add()` [function](https://ropensci.github.io/slopes/reference/elevation_add.html) - this will require a valid [Mapbox api key](https://docs.mapbox.com/api/overview/).

```{r import_dem, message=FALSE}
# Import and plot DEM
u = "https://github.com/U-Shift/Declives-RedeViaria/releases/download/0.2/IsleOfWightNASA_clip.tif"
f = basename(u)
download.file(url = u, destfile = f, mode = "wb")
dem = raster::raster(f)
# res(dem) #27m of resolution
network = iow_network_segments

library(raster)
plot(dem)
plot(sf::st_geometry(network), add = TRUE) #check if they overlay
```

All the required data is prepared to estimate the road segments' gradient.

```{r slopes_values}
# Get the slope value for each segment (abs), using slopes package
library(slopes)
library(geodist)
network$slope = slope_raster(network, dem)
network$slope = network$slope*100 #percentage
summary(network$slope) #check the values
```

Half of the road segments in Isle of Wight have a gradient below 3.1%.
<!-- `r # round(median(network$slope),1)`%. -->

We will adopt a simplistic qualitative classification for **cycling effort uphill**, and compare the number of segments in each class.

```{r classify}
# Classify slopes
network$slope_class = network$slope %>%
  cut(
    breaks = c(0, 3, 5, 8, 10, 20, Inf),
    labels = c("0-3: flat", "3-5: mild", "5-8: medium", "8-10: hard", 
               "10-20: extreme", ">20: impossible"),
    right = F
  )
round(prop.table(table(network$slope_class))*100,1)
```
<!-- It means that **`r # round(prop.table(table(network$slope_class))[[1]]*100)`%** of the roads are flat or almost flat (0-3%) and about **`r # round(prop.table(table(network$slope_class))[[1]]*100)+round(prop.table(table(network$slope_class))[[2]]*100)`%** of the roads are easily cyclable (0-5%). -->
It means that **49%** of the roads are flat or almost flat (0-3%) and about **75%** of the roads are easily cyclable (0-5%).

Now let us put this information on a map (see [here](https://rpubs.com/RobinLovelace/781081) for interactive version).

```{r map, message = FALSE, eval=FALSE}
# more useful information
network$length = st_length(network)

# make an interactive map
library(tmap)
palredgreen = c("#267300", "#70A800", "#FFAA00", "#E60000", "#A80000", "#730000") #color palette
# tmap_mode("view")
tmap_options(basemaps = leaflet::providers$CartoDB.Positron) #basemap

slopemap =
  tm_shape(network) +
  tm_lines(
    col = "slope_class",
    palette = palredgreen,
    lwd = 2, #line width
    title.col = "Slope [%]",
    popup.vars = c("Highway" = "highway",
                   "Length" = "length",
                  "Slope: " = "slope",
                  "Class: " = "slope_class"),
    popup.format = list(digits = 1),
    # id = "slope"
    id = "name" #if it gets too memory consuming, delete this line
  )

slopemap
```


```{r export, echo=TRUE, eval=FALSE}
#export to html
tmap_save(slopemap, "html/SlopesIoW.html") 

# export information as geopackage
st_write(network, "shapefiles/SlopesIoW.gpkg", append=F)
```

#### Result:

![](https://user-images.githubusercontent.com/1825120/121820334-2435f080-cc8a-11eb-962c-79dcba97e459.png)

-   [Isle of Wight (UK)](http://web.tecnico.ulisboa.pt/~rosamfelix/gis/declives/SlopesIoW.html)

## Other examples

-   [São Paulo (Brazil)](https://web.tecnico.ulisboa.pt/~rosamfelix/gis/declives/DeclivesSaoPaulo.html)
-   [Lisbon (Portugal)](https://web.tecnico.ulisboa.pt/~rosamfelix/gis/declives/DeclivesLisboa.html)
-   [Oporto (Portugal)](https://web.tecnico.ulisboa.pt/~rosamfelix/gis/declives/DeclivesPorto_EU.html)
-   [Leeds (UK)](https://web.tecnico.ulisboa.pt/~rosamfelix/gis/declives/SlopesLeeds.html)
-   [Zurich (CH)](https://web.tecnico.ulisboa.pt/~rosamfelix/gis/declives/SlopesZurich.html)


```{r tidyup, include=FALSE}
rm(iow_osm,iow_network_clean,iow_network_segments, iow_network, slopemap)
file.remove(f) # remove the file, tidy up
```

---
title: "An introduction to slopes"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{An introduction to slopes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
bibliography: slope-references.bib
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

```{r setup}
library(slopes)
```

## Introduction

Although there are several ways to name "slope", such as "steepness", "hilliness", "inclination", "aspect", "gradient", "declivity", the referred `slopes` in this package can be defined as the "longitudinal gradient" of linear geographic entities, as defined in the context of rivers by[@cohen_global_2018].  

The package was initially developed to research road slopes to support evidence-based sustainable transport policies.
Accounting for gradient when planning for new cycling infrastructure and road space reallocation for walking and cycling can improve outcomes, for example by helping to identify routes that avoid steep hills.
The package can be used to calculate and visualise slopes of rivers and trajectories representing movement on roads of the type published as open data by @ariza-lopez_dataset_2019.

Data on slopes are useful in many fields of research, including [hydrology](https://en.wikipedia.org/wiki/Stream_gradient), natural hazards (including [flooding](https://www.humanitarianresponse.info/fr/operations/afghanistan/infographic/afg-river-gradient-and-flood-hazard) and [landslide risk management](https://assets.publishing.service.gov.uk/media/57a08d0740f0b652dd0016f4/R7815-ADD017_col.pdf)), recreational and competitive sports such as [cycling](http://theclimbingcyclist.com/gradients-and-cycling-an-introduction/), [hiking](https://trailism.com/trail-grades/), and [skiing](https://www.snowplaza.co.uk/blog/16682-skiing-steeps-what-does-gradient-mean-ski-piste/).
Slopes are also also important in some branches of [transport and emissions modelling](https://www.sciencedirect.com/science/article/pii/S2352146516302642) and [ecology](https://doi.org/10.1016/j.ncon.2016.10.001).
A growing number of people working with geospatial data require accurate estimates of gradient, including:

- Transport planning practitioners who require accurate estimates of roadway gradient for estimating energy consumption, safety and mode shift potential in hilly cities (such as Lisbon, the case study city used in the examples in the documentation). 
- Vehicle routing software developers, who need to build systems are sensitive to going up or down steep hills (e.g. bicycles, trains, and large trucks), such as active travel planning, logistics, and emergency services.
- Natural hazard researchers and risk assessors require estimates of linear gradient to inform safety and mitigation plans associated with project on hilly terrain.
- Aquatic ecologists, flooding researchers and others, who could benefit from estimates of river gradient to support modelling of storm hydrographs

There likely other domains where slopes could be useful, such as agriculture, geology, and civil engineering.

An example of the demand for data provided by the package is a map showing gradients across Sao Paulo (Brazil, see image below) that has received more than 300 'likes' on Twitter and generated conversations:  https://twitter.com/DanielGuth/status/1347270685161304069

![](https://camo.githubusercontent.com/30a3b814dd72aef5b51db635f2ab6e1b6b6c57b856d239822788967a4932d655/68747470733a2f2f7062732e7477696d672e636f6d2f6d656469612f45724a32647238574d414948774d6e3f666f726d61743d6a7067266e616d653d6c61726765){ width=50% } 


## Calculating slopes

The most common slope calculation method is defined by the vertical difference of the final and start point or line height (z1 and z0) divided by the horizontal length that separates them.

$$
s = \Delta z/l
$$

Depending on the purpose of application, it might me relevant to understand how hilliness is estimated.  

![Traffic sign](https://sinalnorte.com/wp-content/uploads/2018/02/A3b.jpg){ width=10% }  

### Measures of route hilliness

There are many ways to measure hilliness, mean distance weighted hilliness being perhaps the most common.
These measures, and their implementation (or current lack thereof) in the package is summarised below.

- **Mean distance weighted gradient**. Perhaps the simplest and most widely applicable measure is the mean gradient of the route. This should be weighted by the distance of each segment. Implemented by default in the `slope_raster()` function.

- **Max gradient**. For activities like cycling, where steep hills have a disproportionate impact, it may be useful to consider the maximum gradient. Not yet implemented.

<!-- Todo: update when we have `method = max` -->

- **Xth percentile gradient**. Since the maximum gradient gives no information about the rest of the route segments, other measures such as the 75th percentile gradient could be more informative. Not yet implemented.

- **Inverted harmonic mean**. If we use the following formula we will get an index that (like the arithmetic mean) makes use of the full dataset, but that is weighted towards the higher gradient segments. Whether this index, the formula of which is shown below, is helpful, remains to be tested. Not yet implemented.

$$
H(x) = 1 - distance.weighted.harmonic.mean(1-x)
$$


### Segments in a route: Cumulative slope 

The length of a segment in a route is also a relevant factor to have in consideration. If it is ok to bike through a segment of 8% with xx length, it is not so ok to bike in four segments in a row like that one (8% for 4xx length), as illustrated bellow.

```{r cumulative-slopes, fig.cap="Illustration of the importance of slope length. 4 segments with an 8% gradient is not the same as a single segment with a gradient of 8%.", out.width="40%", echo=FALSE}
knitr::include_graphics("SLOPES-commulative-slope-1.png")
```


This is accounted for in slope calculation methods that take the distance-weighted mean of slopes.

```{r}
x = c(0, 2, 3, 4, 5, 9)
y = c(0, 0, 0, 0, 0, 9)
z = c(1, 2, 2, 4, 3, 1) / 10
m = cbind(x, y, z)
d = sequential_dist(m = m, lonlat = FALSE)

slopes::slope_distance_weighted(d = d, elevations = z)
slopes::slope_distance_mean(d = d, elevations = z)
```

The slope estimate that results from the distance-weighted mean is lower than the simple mean.
This is common: steep slopes tend to be short!

A graphical representation of the scenario demonstrated above is shown in Figure \@ref(fig:weighted), that shows the relatively long and flat final segment reduces the slope by half.

```{r weighted, fig.cap="Illustration of example data that demonstrates distance-weighted mean gradient, used by default in the slopes package."}
plot(x, z, ylim = c(-0.5, 0.5), type = "l")
(gxy = slope_matrix(m, lonlat = FALSE))
abline(h = 0, lty = 2)
points(x[-length(x)], gxy, col = "blue")
title("Distance elevation profile",
  sub = "Points show calculated gradients of subsequent lines")
```




# References

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slopes.R
\name{slope_vector}
\alias{slope_vector}
\alias{slope_distance}
\alias{slope_distance_mean}
\alias{slope_distance_weighted}
\title{Calculate the gradient of line segments from distance and elevation vectors}
\usage{
slope_vector(x, elevations)

slope_distance(d, elevations)

slope_distance_mean(d, elevations, directed = FALSE)

slope_distance_weighted(d, elevations, directed = FALSE)
}
\arguments{
\item{x}{Vector of locations}

\item{elevations}{Elevations in same units as x (assumed to be metres)}

\item{d}{Vector of distances between points}

\item{directed}{Should the value be directed? \code{FALSE} by default.
If \code{TRUE} the result will be negative when it represents a downslope
(when the end point is lower than the start point).}
}
\value{
A vector of slope gradients associated with each linear element
(each line between consecutive vertices) associated with linear features.
Returned values for \code{slope_distance_mean()} and
\code{slope_distance_mean_weighted()} are summary statistics for all
linear elements in the linestring.
The output value is a proportion representing the change in elevation
for a given change in horizontal movement along the linestring.
0.02, for example, represents a low gradient of 2\% while 0.08 represents
a steep gradient of 8\%.
}
\description{
\code{slope_vector()} calculates the slopes associated with consecutive elements
in one dimensional distance and associated elevations (see examples).

\code{slope_distance()} calculates the slopes associated with consecutive
distances and elevations.

\code{slope_distance_mean()} calculates the mean average slopes associated with
consecutive distances and elevations.

\code{slope_distance_weighted()} calculates the slopes associated with
consecutive distances and elevations,
with the mean value associated with each set of distance/elevation
vectors weighted in proportion to the distance between each elevation
measurement, so longer sections have proportionally more influence
on the resulting gradient estimate (see examples).
}
\examples{
x = c(0, 2, 3, 4, 5, 9)
elevations = c(1, 2, 2, 4, 3, 0) / 10 # downward slope overall
slope_vector(x, elevations)
library(sf)
m = st_coordinates(lisbon_road_segment)
d = sequential_dist(m, lonlat = FALSE)
elevations = elevation_extract(m, dem_lisbon_raster)
slope_distance(d, elevations)
slope_distance_mean(d, elevations)
slope_distance_mean(d, elevations, directed = TRUE)
slope_distance_mean(rev(d), rev(elevations), directed = TRUE)
slope_distance_weighted(d, elevations)
slope_distance_weighted(d, elevations, directed = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slopes.R
\name{slope_matrix}
\alias{slope_matrix}
\alias{slope_matrix_mean}
\alias{slope_matrix_weighted}
\title{Calculate the gradient of line segments from a 3D matrix of coordinates}
\usage{
slope_matrix(m, elevations = m[, 3], lonlat = TRUE)

slope_matrix_mean(m, elevations = m[, 3], lonlat = TRUE, directed = FALSE)

slope_matrix_weighted(m, elevations = m[, 3], lonlat = TRUE, directed = FALSE)
}
\arguments{
\item{m}{Matrix containing coordinates and elevations.
The matrix should have three columns: x, y, and z, in that order. Typically
these correspond to location in the West-East, South-North, and vertical
elevation axes respectively.
In data with geographic coordinates, Z values are assumed to be in
metres. In data with projected coordinates, Z values are assumed to have
the same units as the X and Y coordinates.}

\item{elevations}{Elevations in same units as x (assumed to be metres).
Default value: \code{m[, 3]}, meaning the 'z' coordinate in a matrix of
coordinates.}

\item{lonlat}{Are the coordinates in lon/lat (geographic) coordinates? TRUE by default.}

\item{directed}{Should the value be directed? \code{FALSE} by default.
If \code{TRUE} the result will be negative when it represents a downslope
(when the end point is lower than the start point).}
}
\value{
A vector of slope gradients associated with each linear element
(each line between consecutive vertices) associated with linear features.
Returned values for \code{slope_matrix_mean()} and
\code{slope_matrix_weighted()} are summary statistics for all
linear elements in the linestring.
The output value is a proportion representing the change in elevation
for a given change in horizontal movement along the linestring.
0.02, for example, represents a low gradient of 2\% while 0.08 represents
a steep gradient of 8\%.
}
\description{
Calculate the gradient of line segments from a 3D matrix of coordinates
}
\examples{
x = c(0, 2, 3, 4, 5, 9)
y = c(0, 0, 0, 0, 0, 9)
z = c(1, 2, 2, 4, 3, 0) / 10
m = cbind(x, y, z)
slope_matrix_weighted(m, lonlat = FALSE)
slope_matrix_weighted(m, lonlat = FALSE, directed = TRUE)
# 0 value returned if no change in elevation:
slope_matrix_weighted(m,lonlat = FALSE, directed = TRUE,
  elevations = c(1, 2, 2, 4, 3, 1))
slope_matrix_mean(m, lonlat = FALSE)
slope_matrix_mean(m, lonlat = FALSE, directed = TRUE)
plot(x, z, ylim = c(-0.5, 0.5), type = "l")
(gx = slope_vector(x, z))
(gxy = slope_matrix(m, lonlat = FALSE))
abline(h = 0, lty = 2)
points(x[-length(x)], gx, col = "red")
points(x[-length(x)], gxy, col = "blue")
title("Distance (in x coordinates) elevation profile",
  sub = "Points show calculated gradients of subsequent lines")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slope_get.R
\name{elevation_get}
\alias{elevation_get}
\title{Get elevation data from hosted maptile services}
\usage{
elevation_get(routes, ..., output_format = "raster")
}
\arguments{
\item{routes}{Routes, the gradients of which are to be calculated.
The object must be of class \code{sf} or \code{sfc} with \code{LINESTRING} geometries.}

\item{...}{Options passed to \code{cc_elevation()}}

\item{output_format}{What format to return the data in?
Accepts \code{"raster"} (the default) and \code{"terra"}.}
}
\value{
A raster object with cell values representing elevations in the
bounding box of the input \code{routes} object.
}
\description{
\code{elevation_get()} uses the
\href{https://hypertidy.github.io/ceramic/reference/cc_location.html}{\code{cc_elevation()}}
function from the \code{ceramic} package to get
DEM data in raster format anywhere worldwide.
It requires an API that can be added by following guidance in the package's
\href{https://github.com/ropensci/slopes}{README}
and in the
\href{https://ropensci.github.io/slopes/articles/slopes.html}{\code{slopes} vignette}.
}
\details{
Note: if you use the \code{cc_elevation()} function directly to get DEM data,
you can cache the data, as described in the package's
\href{https://github.com/hypertidy/ceramic#local-caching-of-tiles}{README}.
}
\examples{
# Time-consuming examples that require an internet connection and API key:
\donttest{
library(sf)
library(raster)
routes = cyclestreets_route
e = elevation_get(routes)
class(e)
crs(e)
e
plot(e)
plot(st_geometry(routes), add = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slopes.R
\name{sequential_dist}
\alias{sequential_dist}
\title{Calculate the sequential distances between sequential coordinate pairs}
\usage{
sequential_dist(m, lonlat = TRUE)
}
\arguments{
\item{m}{Matrix containing coordinates and elevations.
The matrix should have three columns: x, y, and z, in that order. Typically
these correspond to location in the West-East, South-North, and vertical
elevation axes respectively.
In data with geographic coordinates, Z values are assumed to be in
metres. In data with projected coordinates, Z values are assumed to have
the same units as the X and Y coordinates.}

\item{lonlat}{Are the coordinates in lon/lat (geographic) coordinates? TRUE by default.}
}
\value{
A vector of distance values in meters if \code{lonlat = TRUE}
or the map units of the input data if \code{lonlat = FALSE} between
consecutive vertices.
}
\description{
Set \code{lonlat} to \code{FALSE} if you have projected data, e.g. with coordinates
representing distance in meters, not degrees. Lonlat coodinates are assumed
(\code{lonlat = TRUE} is the default).
}
\examples{
x = c(0, 2, 3, 4, 5, 9)
y = c(0, 0, 0, 0, 0, 1)
m = cbind(x, y)
d = sequential_dist(m, lonlat = FALSE)
d
nrow(m)
length(d)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slopes-package.R
\docType{package}
\name{slopes-package}
\alias{slopes}
\alias{slopes-package}
\title{slopes: Calculate Slopes of Roads, Rivers and Trajectories}
\description{
Calculate the slope (also known as steepness and gradient)
of linear features such as roads, rivers and ski routes.

Created by Robin Lovelace and Rosa Félix, who developed the package
to support their research in active travel, and the prioritisation of
investment in cycling, an activity which is sensitive to slopes,
in particular.

See \href{https://github.com/ropensci/slopes}{github.com/ITSLeeds/slopes}
for the source code and \href{https://ropensci.github.io/slopes/}{ropensci.github.io/slopes}
for the website.

All package functions can be found in the
\href{https://ropensci.github.io/slopes/reference/index.html}{Reference page}
on the package website.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/slopes/}
  \item \url{https://docs.ropensci.org/slopes/}
  \item Report bugs at \url{https://github.com/ropensci/slopes/issues}
}

}
\author{
\strong{Maintainer}: Robin Lovelace \email{rob00x@gmail.com} (\href{https://orcid.org/0000-0001-5679-6536}{ORCID})

Authors:
\itemize{
  \item Rosa Félix \email{rosamfelix@tecnico.ulisboa.pt} (\href{https://orcid.org/0000-0002-5642-6006}{ORCID})
  \item Joey Talbot \email{j.d.talbot@leeds.ac.uk} (\href{https://orcid.org/0000-0002-6520-4560}{ORCID})
}

Other contributors:
\itemize{
  \item Dan Olner (Dan reviewed the package for rOpenSci, see https://github.com/ropensci/software-review/issues/420#issuecomment-857662657 ) [reviewer]
  \item Andy Teucher (Andy reviewed the package for rOpenSci, see https://github.com/ropensci/software-review/issues/420#issuecomment-858231647 ) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slopes.R
\name{slope_raster}
\alias{slope_raster}
\title{Calculate the gradient of line segments from a raster dataset}
\usage{
slope_raster(
  routes,
  dem,
  lonlat = sf::st_is_longlat(routes),
  method = "bilinear",
  fun = slope_matrix_weighted,
  terra = has_terra() && methods::is(dem, "SpatRaster"),
  directed = FALSE
)
}
\arguments{
\item{routes}{Routes, the gradients of which are to be calculated.
The object must be of class \code{sf} or \code{sfc} with \code{LINESTRING} geometries.}

\item{dem}{Raster overlapping with \code{routes} and values representing elevations}

\item{lonlat}{Are the routes provided in longitude/latitude coordinates?
By default, value is from the CRS of the routes (\code{sf::st_is_longlat(routes)}).}

\item{method}{The method of estimating elevation at points,
passed to the \code{extract} function for extracting values from raster
datasets. Default: \code{"bilinear"}.}

\item{fun}{The slope function to calculate per route,
\code{slope_matrix_weighted} by default.}

\item{terra}{Should the \code{terra} package be used?
\code{TRUE} by default if the package is installed \emph{and}
if \code{dem} is of class \code{SpatRast}}

\item{directed}{Should the value be directed? \code{FALSE} by default.
If \code{TRUE} the result will be negative when it represents a downslope
(when the end point is lower than the start point).}
}
\value{
A vector of slopes equal in length to the number simple features
(rows representing linestrings) in the input object.
}
\description{
This function takes an \code{sf} representing routes over geographical space
and a raster dataset representing the terrain as inputs.
It returns the average gradient of each route feature.
}
\details{
If calculating slopes associated with OSM data, the results may be better
if the network is first split-up, e.g. using the function
\code{stplanr::rnet_breakup_vertices()} from the
\href{https://docs.ropensci.org/stplanr/reference/}{\code{stplanr}} package.
\strong{Note:} The \code{routes} object must have a geometry type of \code{LINESTRING}.
The \code{sf::st_cast()} function can convert from \code{MULTILINESTRING} (and other)
geometries to \code{LINESTRING}s as follows:
\code{r_linestring = sf::st_cast(routes, "LINESTRING")}.
}
\examples{
library(sf)
routes = lisbon_road_network[1:3, ]
dem = dem_lisbon_raster
(s = slope_raster(routes, dem))
cor(routes$Avg_Slope, s)
slope_raster(routes, dem, directed = TRUE)
# Demonstrate that reverse routes have the opposite directed slope
slope_raster(st_reverse(routes), dem, directed = TRUE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{lisbon_road_segment}
\alias{lisbon_road_segment}
\alias{lisbon_road_segment_3d}
\alias{lisbon_road_segment_xyz_mapbox}
\title{A road segment in Lisbon, Portugal}
\format{
An object of class \code{sf}
}
\source{
Produced by ESRI's
\href{http://pro.arcgis.com/en/pro-app/latest/help/analysis/}{3D Analyst extension}
}
\usage{
lisbon_road_segment
}
\description{
A single road segment and a 3d version.
Different versions of this dataset are provided.
}
\details{
The \code{lisbon_road_segment} has 23 columns and 1 row.

The \code{lisbon_road_segment_xyz_mapbox} was created with:
\code{lisbon_road_segment_xyz_mapbox = elevation_add(lisbon_road_segment)}.
}
\examples{
lisbon_road_segment
lisbon_road_segment_3d
lisbon_road_segment_xyz_mapbox
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slopes.R
\name{elevation_add}
\alias{elevation_add}
\title{Take a linestring and add a third (z) dimension to its coordinates}
\usage{
elevation_add(
  routes,
  dem = NULL,
  method = "bilinear",
  terra = has_terra() && methods::is(dem, "SpatRaster")
)
}
\arguments{
\item{routes}{Routes, the gradients of which are to be calculated.
The object must be of class \code{sf} or \code{sfc} with \code{LINESTRING} geometries.}

\item{dem}{Raster overlapping with \code{routes} and values representing elevations}

\item{method}{The method of estimating elevation at points,
passed to the \code{extract} function for extracting values from raster
datasets. Default: \code{"bilinear"}.}

\item{terra}{Should the \code{terra} package be used?
\code{TRUE} by default if the package is installed \emph{and}
if \code{dem} is of class \code{SpatRast}}
}
\value{
An sf object that is identical to the input \code{routes}, except that
the coordinate values in the ouput has a third \code{z} dimension representing
the elevation of each vertex that defines a linear feature such as a road.
}
\description{
Take a linestring and add a third (z) dimension to its coordinates
}
\examples{
library(sf)
routes = lisbon_road_network[204, ]
dem = dem_lisbon_raster
(r3d = elevation_add(routes, dem))
library(sf)
st_z_range(routes)
st_z_range(r3d)
plot(st_coordinates(r3d)[, 3])
plot_slope(r3d)
\donttest{
# Get elevation data (requires internet connection and API key):
r3d_get = elevation_add(cyclestreets_route)
plot_slope(r3d_get)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{magnolia_xy}
\alias{magnolia_xy}
\title{Road segments in Magnolia, Seattle}
\format{
An object of class \code{sf}
}
\source{
Accessed in early 2021 from the \code{seattle-streets} layer from the
\href{https://data-seattlecitygis.opendata.arcgis.com/}{data-seattlecitygis}
website.
}
\usage{
magnolia_xy
}
\description{
A dataset representing road segments in the Magnolia area of Seattle
with X, Y and Z (elevation) dimensions for each coordinate.
}
\examples{
names(magnolia_xy)
plot(magnolia_xy["SLOPE_PCT"])
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slopes.R
\name{slope_xyz}
\alias{slope_xyz}
\title{Extract slopes from xyz data frame or sf objects}
\usage{
slope_xyz(
  route_xyz,
  fun = slope_matrix_weighted,
  lonlat = TRUE,
  directed = FALSE
)
}
\arguments{
\item{route_xyz}{An \code{sf} or \code{sfc} object with \code{XYZ} coordinate dimensions}

\item{fun}{The slope function to calculate per route,
\code{slope_matrix_weighted} by default.}

\item{lonlat}{Are the coordinates in lon/lat order? TRUE by default}

\item{directed}{Should the value be directed? \code{FALSE} by default.
If \code{TRUE} the result will be negative when it represents a downslope
(when the end point is lower than the start point).}
}
\value{
A vector of slopes equal in length to the number simple features
(rows representing linestrings) in the input object.
}
\description{
The function takes a sf object with 'XYZ' coordinates and returns a vector
of numeric values representing the average slope of each linestring in the
sf data frame input.
}
\details{
The default function to calculate the mean slope is \code{slope_matrix_weighted()}.
You can also use \code{slope_matrix_mean()} from the package or any other
function that takes the same inputs as these functions not in the package.
}
\examples{
route_xyz = lisbon_road_segment_3d
slope_xyz(route_xyz, lonlat = FALSE)
slope_xyz(route_xyz$geom, lonlat = FALSE)
slope_xyz(route_xyz, lonlat = FALSE, directed = TRUE)
slope_xyz(route_xyz, lonlat = FALSE, fun = slope_matrix_mean)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/z.R
\name{z_value}
\alias{z_value}
\alias{z_start}
\alias{z_end}
\alias{z_mean}
\alias{z_max}
\alias{z_min}
\alias{z_elevation_change_start_end}
\alias{z_direction}
\alias{z_cumulative_difference}
\title{Calculate summary values for 'Z' elevation attributes}
\usage{
z_value(x)

z_start(x)

z_end(x)

z_mean(x)

z_max(x)

z_min(x)

z_elevation_change_start_end(x)

z_direction(x)

z_cumulative_difference(x)
}
\arguments{
\item{x}{An \code{sfc} object with 'XYZ' coordinates}
}
\value{
A vector of values representing elevations associated with
simple feature geometries that have elevations (XYZ coordinates).
}
\description{
The \verb{slope_z*()} functions calculate summary values for the Z axis
in \code{sfc} objects with \code{XYZ} geometries.
}
\examples{
x = slopes::lisbon_route_3d
x
z_value(x)[1:5]
xy = slopes::lisbon_route
try(z_value(xy)) # error message
z_start(x)
z_end(x)
z_direction(x)
z_elevation_change_start_end(x)
z_direction(x)
z_cumulative_difference(x)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_slope.R
\name{plot_dz}
\alias{plot_dz}
\title{Plot a digital elevation profile based on xyz data}
\usage{
plot_dz(
  d,
  z,
  fill = TRUE,
  horiz = FALSE,
  pal = colorspace::diverging_hcl,
  ...,
  legend_position = "top",
  col = "black",
  cex = 0.9,
  bg = grDevices::rgb(1, 1, 1, 0.8),
  title = "Slope colors (percentage gradient)",
  brks = NULL,
  seq_brks = NULL,
  ncol = 4
)
}
\arguments{
\item{d}{Cumulative distance}

\item{z}{Elevations at points across a linestring}

\item{fill}{Should the profile be filled? \code{TRUE} by default}

\item{horiz}{Should the legend be horizontal (\code{FALSE} by default)}

\item{pal}{Color palette to use, \code{colorspace::diverging_hcl} by default.}

\item{...}{Additional parameters to pass to legend}

\item{legend_position}{The legend position. One of "bottomright", "bottom",
"bottomleft", "left", "topleft", "top" (the default), "topright", "right"
and "center".}

\item{col}{Line colour, black by default}

\item{cex}{Legend size, 0.9 by default}

\item{bg}{Legend background colour, \code{grDevices::rgb(1, 1, 1, 0.8)} by default.}

\item{title}{Title of the legend, \code{NULL} by default.}

\item{brks}{Breaks in colour palette to show.
\code{c(1, 3, 6, 10, 20, 40, 100)} by default.}

\item{seq_brks}{Sequence of breaks to show in legend.
Includes negative numbers and omits zero by default}

\item{ncol}{Number of columns in legend, 4 by default.}
}
\value{
A plot showing the elevation profile associated with a linestring.
}
\description{
Plot a digital elevation profile based on xyz data
}
\examples{
library(sf)
route_xyz = lisbon_road_segment_3d
m = st_coordinates(route_xyz)
d = cumsum(sequential_dist(m, lonlat = FALSE))
d = c(0, d)
z = m[, 3]
slopes:::plot_dz(d, z, brks = c(3, 6, 10, 20, 40, 100))
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{lisbon_route}
\alias{lisbon_route}
\alias{lisbon_route_3d}
\alias{lisbon_route_xyz_mapbox}
\title{A route composed of a single linestring in Lisbon, Portugal}
\format{
An object of class \code{sf}
}
\source{
See the \code{lisbon_route.R} script in \code{data-raw}
}
\usage{
lisbon_route
}
\description{
A route representing a trip from the Santa Catarina area
in the East of central Lisbon the map to the Castelo de São Jorge
in the West of central Lisbon.
}
\details{
Different versions of this dataset are provided.

The \code{lisbon_route} object has 1 row and 4 columns: geometry, ID,
length and whether or not a path was found.

The \code{lisbon_route_xyz_mapbox} was created with:
\code{lisbon_route_xyz_mapbox = elevation_add(lisbon_route)}.
}
\examples{
lisbon_route
lisbon_route_3d
lisbon_route_xyz_mapbox
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{dem_lisbon_raster}
\alias{dem_lisbon_raster}
\title{Elevation in central Lisbon, Portugal}
\format{
A raster dataset containing elevation above sea level
in a 1km bounding box in Lisbon, Portugal.
}
\source{
\url{https://github.com/rspatial/terra/issues/29}
}
\usage{
dem_lisbon_raster
}
\description{
A dataset containing elevation in and around Lisbon
with a geographic resolution of 10m.
The dataset is 200 pixels wide by 133 pixels high, covering
2.7 square kilometres of central Lisbon.
}
\details{
The dataset was acquired by Instituto Superior
Técnico (University of Lisbon) in 2012, covers all the Northern
Metropolitan Area of Lisbon, and has a 10m cell resolution,
when projected at the official Portuguese EPSG: 3763 - TM06/ETRS89.
The dataset was released as an open access dataset with permission from the
University of Lisbon to support this project.
}
\examples{
library(sf)
library(raster)
dim(dem_lisbon_raster)
res(dem_lisbon_raster)
names(dem_lisbon_raster)
plot(dem_lisbon_raster)
plot(lisbon_road_network["Avg_Slope"], add = TRUE)
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/slopes.R
\name{elevation_extract}
\alias{elevation_extract}
\title{Extract elevations from coordinates}
\usage{
elevation_extract(
  m,
  dem,
  method = "bilinear",
  terra = has_terra() && methods::is(dem, "SpatRaster")
)
}
\arguments{
\item{m}{Matrix containing coordinates and elevations or an sf
object representing a linear feature.}

\item{dem}{Raster overlapping with \code{routes} and values representing elevations}

\item{method}{The method of estimating elevation at points,
passed to the \code{extract} function for extracting values from raster
datasets. Default: \code{"bilinear"}.}

\item{terra}{Should the \code{terra} package be used?
\code{TRUE} by default if the package is installed \emph{and}
if \code{dem} is of class \code{SpatRast}}
}
\value{
A vector of elevation values.
}
\description{
This function takes a series of points located in geographical space
and a digital elevation model as inputs and returns a vector of
elevation estimates associated with each point.
The function takes locations
represented as a matrix of XY (or longitude latitude) coordinates
and a digital elevation model (DEM) with class \code{raster} or \code{terra}.
It returns a vector of values representing estimates of elevation
associated with each of the points.
}
\details{
By default, the elevations are estimated using
\href{https://en.wikipedia.org/wiki/Bilinear_interpolation}{bilinear interpolation}
(\code{method = "bilinear"})
which calculates point height based on proximity to the centroids of
surrounding cells.
The value of the \code{method} argument is passed to the \code{method} argument in
\href{https://rspatial.github.io/raster/reference/extract.html}{\code{raster::extract()}}
or
\href{https://rspatial.github.io/terra/reference/extract.html}{\code{terra::extract()}}
depending on the class of the input raster dataset.

See Kidner et al. (1999)
for descriptions of alternative elevation interpolation and extrapolation
algorithms.
}
\examples{
dem = dem_lisbon_raster
elevation_extract(lisbon_road_network[1, ], dem)
m = sf::st_coordinates(lisbon_road_network[1, ])
elevation_extract(m, dem)
elevation_extract(m, dem, method = "simple")
# Test with terra (requires internet connection):
\donttest{
if(slopes:::has_terra()) {
et = terra::rast(dem_lisbon_raster)
elevation_extract(m, et)
}
}
}
\references{
Kidner, David, Mark Dorey, and Derek Smith.
"What’s the point? Interpolation and extrapolation with a regular grid DEM."
Fourth International Conference on GeoComputation, Fredericksburg,
VA, USA. 1999.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_slope.R
\name{plot_slope}
\alias{plot_slope}
\title{Plot slope data for a 3d linestring with base R graphics}
\usage{
plot_slope(
  route_xyz,
  lonlat = sf::st_is_longlat(route_xyz),
  fill = TRUE,
  horiz = FALSE,
  pal = colorspace::diverging_hcl,
  legend_position = "top",
  col = "black",
  cex = 0.9,
  bg = grDevices::rgb(1, 1, 1, 0.8),
  title = "Slope colors (percentage gradient)",
  brks = c(3, 6, 10, 20, 40, 100),
  seq_brks = seq(from = 3, to = length(brks) * 2 - 2),
  ncol = 4,
  ...
)
}
\arguments{
\item{route_xyz}{An sf linestring with x, y and z coordinates,
representing a route or other linear object.}

\item{lonlat}{Are the routes provided in longitude/latitude coordinates?
By default, value is from the CRS of the routes (\code{sf::st_is_longlat(routes)}).}

\item{fill}{Should the profile be filled? \code{TRUE} by default}

\item{horiz}{Should the legend be horizontal (\code{FALSE} by default)}

\item{pal}{Color palette to use, \code{colorspace::diverging_hcl} by default.}

\item{legend_position}{The legend position. One of "bottomright", "bottom",
"bottomleft", "left", "topleft", "top" (the default), "topright", "right"
and "center".}

\item{col}{Line colour, black by default}

\item{cex}{Legend size, 0.9 by default}

\item{bg}{Legend background colour, \code{grDevices::rgb(1, 1, 1, 0.8)} by default.}

\item{title}{Title of the legend, \code{NULL} by default.}

\item{brks}{Breaks in colour palette to show.
\code{c(1, 3, 6, 10, 20, 40, 100)} by default.}

\item{seq_brks}{Sequence of breaks to show in legend.
Includes negative numbers and omits zero by default}

\item{ncol}{Number of columns in legend, 4 by default.}

\item{...}{Additional parameters to pass to legend}
}
\value{
A plot showing the elevation profile associated with a linestring.
}
\description{
Plot slope data for a 3d linestring with base R graphics
}
\examples{
plot_slope(lisbon_route_3d)
route_xyz = lisbon_road_segment_3d
plot_slope(route_xyz)
plot_slope(route_xyz, brks = c(1, 2, 4, 8, 16, 30))
plot_slope(route_xyz, s = 5:8)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{cyclestreets_route}
\alias{cyclestreets_route}
\title{A journey from CycleStreets.net}
\format{
An object of class \code{sf} with 18 rows and 14 columns on route
characteristics. See https://rpackage.cyclestreets.net/reference/journey.html
for details.
}
\source{
CycleStreets.net
}
\usage{
cyclestreets_route
}
\description{
Road segments representing suggested route to cycle
in Leeds, UK.
}
\details{
Simple feature collection with 30 features and 32 fields

See \code{data-raw/cyclestreets_route.R} in the package's github repo for details.
}
\examples{
library(sf)
class(cyclestreets_route)
plot(cyclestreets_route$geometry)
cyclestreets_route
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{lisbon_road_network}
\alias{lisbon_road_network}
\title{Road segments in Lisbon}
\format{
An object of class \code{sf}, key variables of which include
\describe{
\item{OBJECTID}{ID of the object}
\item{Z_Min}{The minimum elevation on the linear feature from ArcMAP}
\item{Z_Max}{The max elevation on the linear feature from ArcMAP}
\item{Z_Mean}{The mean elevation on the linear feature from ArcMAP}
\item{Slope_Min}{The minimum slope on the linear feature from ArcMAP}
\item{Slope_Max}{The max slope on the linear feature from ArcMAP}
\item{Slope_Mean}{The mean slope on the linear feature from ArcMAP}
\item{geom}{The geometry defining the LINESTRING component of the segment}
}
}
\source{
Produced by ESRI's
\href{http://pro.arcgis.com/en/pro-app/latest/help/analysis/}{3D Analyst extension}
}
\usage{
lisbon_road_network
}
\description{
A dataset representing road segments in Lisbon,
with X, Y and Z (elevation) dimensions for each coordinate.
}
\details{
The dataset covers 32 km of roads in central Lisbon, overlapping with the
area covered by the \code{dem_lisbon_raster} dataset.
}
\examples{
library(sf)
names(lisbon_road_network)
sum(st_length(lisbon_road_network))
plot(lisbon_road_network["Avg_Slope"])
}
\keyword{datasets}
