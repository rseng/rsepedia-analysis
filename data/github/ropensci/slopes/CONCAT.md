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
