nsink
================

<!-- README.md is generated from README.Rmd. Please edit that file -->
<!-- badges: start -->

[![R build
status](https://github.com/jhollist/nsink/workflows/R-CMD-check/badge.svg)](https://github.com/jhollist/nsink/actions)
[![Codecov test
coverage](https://codecov.io/gh/jhollist/nsink/branch/main/graph/badge.svg)](https://codecov.io/gh/jhollist/nsink?branch=main)
[![Lifecycle:
stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6341565.svg)](https://doi.org/10.5281/zenodo.6341565)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04039/status.svg)](https://doi.org/10.21105/joss.04039)
<!-- badges: end -->

# Statement of need

The `nsink` package is an R implementation of the methods described in
[Kellogg et. al (2010)](https://doi.org/10.1016/j.ecoleng.2010.02.006).
Previous implementation of this approach relied on a manual, vector
based approach that was time consuming to prepare. This approach uses a
hybrid raster-vector approach that takes relatively little time to set
up for each new watershed and relies on readily available data. Total
run times vary, but range from minutes up to 5 hours depending on
options selected. Previous versions took weeks of manual data
manipulation. Thus, `nsink` was developed to satisfy the need for
quicker implementation of the NSink method as described in [Kellogg et.
al (2010)](https://doi.org/10.1016/j.ecoleng.2010.02.006).

# `nsink` functionality

As of 2022-03-14 user functions for the `nsink` package are:

-   `nsink_get_huc_id()`: A function for searching the name of a USGS
    Watershed Boundary Dataset Hydrologic Unit
    (<https://www.usgs.gov/core-science-systems/ngp/national-hydrography/watershed-boundary-dataset>)
    and retrieving its 12-digit Hydrologic Unit Code (HUC).  
-   `nsink_get_data()`: Using any acceptable HUC ID (e.g. 2-digit to
    12-digit), this function downloads the NHDPlus, SSURGO, NLCD Land
    Cover, and the NLCD Impervious for that HUC.  
-   `nsink_prep_data()`: `nsink` needs data in a common coordinate
    reference system, from mutliple NHDPlus tables, and from different
    portions of SSURGO. This function completes these data preparation
    steps and outputs all data, clipped to the HUC boundary.
-   `nsink_calc_removal()`: Quantifying relative N removal across a
    landscape is a key aspects of an `nsink` analysis. The
    `nsink_calc_removal()` function takes the object returned from
    `nsink_prep_data()` and calculates relative N removal for each
    landscape sink. See Kellogg et al \[-@kellogg2010geospatial\] for
    details on relative N removal estimation for each sink.
-   `nsink_generate_flowpath()`: This function uses a combination of
    flow determined by topography, via a flow-direction raster, for the
    land-based portions of a flow path and of downstream flow along the
    NHDPlus stream network.  
-   `nsink_summarize_flowpath()`: Summarizing removal along a specified
    flow path requires relative N removal and a generated flow path.
    This function uses these and returns a summary of relative N removal
    along a flow path for each sink.
-   `nsink_generate_static_maps()`: This function analyzes N removal at
    the watershed scale by summarizing the results of multiple flow
    paths. Four static maps are returned: 1)removal efficiency;
    2)loading index; 3)transport index; 4)delivery index. Removal
    efficiency is a rasterized version of the `nsink_calc_removal()`
    output. Loading index is N sources based on NLCD categories.
    Transport index is a heat map with the cumulative relative N removal
    along flow paths originating from a grid of points, density set by
    the user, across a watershed, highlighting the gradient of
    downstream N retention. Delivery index is the result of multiplying
    the loading index and the transport index, and shows potential N
    delivery from different sources, taking into account the relative N
    removal as water moves downstream.
-   `nsink_plot()`: A function that plots each raster in the list
    returned from `nsink_generate_static_maps()`.  
-   `nsink_build()`: One of the drivers behind the development of the
    `nsink` package was to provide `n-sink` analysis output that could
    be used more broadly (e.g. within a GIS). The `nsink_build()` runs a
    complete `nsink` analysis and outputs R objects, shapefiles and/or
    TIFFs.
-   `nsink_load()`: Essentially the inverse of the `nsink_build()`
    function, this function takes a folder of files, likely created by
    `nsink_build()`, and reads them into R.

# Installation instructions

At this time we plan on maintaining the `nsink` package as a GitHub only
package and thus it won’t be available directly from CRAN. You may use
the `install_github()` function from the `remotes` package to install
it. The code below will take care of installing `remotes` and installing
`nsink` from the GitHub repository.

``` r
install.packages("remotes")
remotes::install_github("usepa/nsink", dependencies = TRUE, build_vignettes = TRUE)
```

And then to load up the package:

``` r
library(nsink)
```

# Documentation and examples

All functions are documented, with examples, and that documentation may
be accessed, in R, via the usual help functions. Additionally, an
introduction to the `nsink` package with a more detailed workflow is
documented in a vignette.

``` r
# Load up package
library(nsink)

# Access package level help
help(package = "nsink")

# Access the Introduction to nsink vignette
vignette("intro", package = "nsink")
```

# Contributing

If you would like to contribute to the `nsink` package, please first
read the [CONTRIBUTING](.github/CONTRIBUTING.md). In short,
contributions are happily accepted either via suggestions in the
[Issues](https://github.com/USEPA/nsink/issues) or via pull request.
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
reported to the community leaders responsible for enforcement at hollister.jeff@epa.gov. 
All complaints will be reviewed and investigated promptly and fairly.

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
available at <https://www.contributor-covenant.org/version/2/0/code_of_conduct.html>.

Community Impact Guidelines were inspired by [Mozilla's code of conduct
enforcement ladder](https://github.com/mozilla/diversity).

[homepage]: https://www.contributor-covenant.org

For answers to common questions about this code of conduct, see the FAQ at
<https://www.contributor-covenant.org/faq>. Translations are available at <https://www.contributor-covenant.org/translations>.
---
title: 'nsink: An R package for flow path nitrogen removal estimation'
tags:
- R
- nitrogen
- nitrogen sinks
- landscape
- gis
date: "2021-12-20"
output: github_document
authors:
- name: Jeffrey W. Hollister
  orcid: 0000-0002-9254-9740
  affiliation: 1
- name: Dorothy Q. Kellogg
  orcid: 0000-0002-9509-4606
  affiliation: 2
- name: Qian Lei-Parent
  orcid: 0000-0002-1904-2513
  affiliation: 3
- name: Emily Wilson
  orcid: 0000-0003-0035-5752
  affiliation: 3
- name: Cary Chadwick
  orcid: 0000-0002-6952-7535
  affiliation: 3
- name: David Dickson
  orcid: 0000-0002-2660-6460
  affiliation: 3
- name: Arthur Gold
  orcid: 0000-0002-0290-1377
  affiliation: 2
- name: Chester Arnold
  affiliation: 3

bibliography: paper.bib
affiliations:
- name: U. S. Environmental Protection Agency, Atlantic Coastal Environmental Sciences
    Division, Narragansett, RI 02882
  index: 1
- name: University of Rhode Island, Department of Natural Resources Science, Kingston, RI 02881
  index: 2
- name: University of Connecticut, Center for Land Use Education and Research, Storrs, CT 06268
  index: 3
---
  
  
# Summary

The `nsink` package estimates cumulative nitrogen (N) removal along a specified flow path and is based on methodologies outlined in Kellogg et al. [ -@kellogg2010geospatial]. For a user-specified watershed (i.e., hydrologic unit code (HUC)), `nsink` downloads all required datasets from public datasets in the United States, prepares data for use, summarizes N removal along a flow path and creates several static maps.  The results of an `nsink` analysis may be exported to standard geospatial files for use in other applications.  

# Statement of need

Excess N delivery via surface water to downstream aquatic resources contributes to impaired water quality and impacts ecosystem services including harmful algal blooms (HABs) and hypoxia [@rabalais2002beyond]. Identifying landscape N sinks (i.e., areas where N is effectively removed from the aquatic system) and analyzing N delivery at the watershed scale is helpful to watershed managers, land use planners and conservation organizations.  The theoretical underpinnings for identifying N sinks rely on decades of research and are explained in Kellogg et al. [-@kellogg2010geospatial]. 

Prior N-sink implementations were done case-by-case.  Data acquisition and manipulation were mostly manual and took weeks to months to complete for a single 12-digit HUC.  The effort required for the analysis limited it's application as scaling beyond a few pilot studies was not feasible.  The goal of `nsink` was to address this limitation and provide an open source solution that could be run on a single small watershed (e.g., 12-digit HUC) in minutes to hours with minimal manual input.

# The `nsink` package

## Package Installation
The `nsink` package is available from <https://github.com/usepa/nsink> and may be installed in R with the following:

```r
# If not installed, install remotes
install.packages("remotes")

# Install nsink from GitHub
remotes::install_github("USEPA/nsink", dependencies = TRUE, build_vignettes = TRUE)
```

## Package Details

The `nsink` package is designed around the major steps in running an N-Sink analysis and includes functions for the following tasks:

1. Prepare for analysis
    - Get data
    - Prepare data for analysis
    - Calculate relative N removal layer for hydric soils, lakes and streams.
2. Run a point-based analysis 
    - Calculate a flow path 
    - Summarize relative N removal along a flow path
3. Run a HUC-based analysis
    - Develop static maps
    - Generate output datasets

### Required Data

The ability to run an `nsink` analysis relies on several datasets for the conterminous United States.  By limiting our approach to these national datasets we are ensuring scalability of `nsink` because the datasets will be available for most locations in the United States.  The datasets that `nsink` uses are the National Hydrography Dataset Plus version 2 (NHDPlus), Soil Survey Geographic Database (SSURGO), the National Land Cover Dataset (NLCD) land cover, and the National Land Cover Dataset (NLCD) impervious surface [@moore2019user; @soil2017web; @jin2019overall]. These datasets are all available via an Application Programming Interface (API) or via direct download.   

### Dependencies

The `nsink` package depends on several existing R packages to facilitate spatial data handling, data acquisition, data management, data analysis and data processing.  These are detailed in Table 1.  

Table 1. R package dependencies for the `nsink` package

|Package|Task|Citation|
|-------|----|--------|
|`sf`|Spatial Data Handling and Analysis|@sfpaper; @sf|
|`raster`|Spatial Data Handling and Analysis|@raster|
|`stars`|Spatial Data Handling and Analysis|@stars|
|`fasterize`|Spatial Data Handling and Analysis|@fasterize|
|`lwgeom`|Spatial Data Handling and Analysis|@lwgeom|
|`gstat`|Spatial Data Handling and Analysis|@gstatpaper2004; @gstatpaper2016; @gstat|
|`sp`|Spatial Data Handling and Analysis|@sppaper; @spbook; @sp|
|`units`|Unit Transformations|@unitspaper; @units|
|`FedData`|Data Acquisition|@feddata|
|`httr`|Data Acquisition|@httr|
|`dplyr`|Data Management and Analysis|@dplyr|
|`zoo`|Data Management and Analysis|@zoopaper; @zoo|
|`igraph`|Data Management and Analysis|@igraphpaper; @igraph|
|`readr`|Data Management and Analysis|@readr|
|`foreign`|Data Management and Analysis|@foreign|
|`rlang`|Data Management and Analysis|@rlang|
|`furrr`|Parallel Processing|@furrr|
|`future`|Parallel Processing|@futurepaper; @future|


### Functionality

Currently, `nsink` provides 10 exported functions to facilitate a flow path analysis of relative N removal. The `nsink` repository (<https://github.com/usepa/nsink>) and R package documentation contain detailed documentation of each function.  The pacakge also has a vignette that outlines a typical workflow for running an N-Sink analysis with the `nsink` package.  Upon install, the vignette is accessed in R with `vignette("intro", package = "nsink")`. 

# Acknowledgements
    
Many people have contributed in various ways to the development of the N-Sink concept.  In particular, we would like to thank, Chet Arnold, Cary Chadwick, David Dickson, and Emily Wilson of the University of Connecticut's Center for Land Use Education and Research as well as Peter August, Chris Damon, and Art Gold of the University of Rhode Island's Department of Natural Resources Science.  Both the UCONN and URI crews have contributed tremendously to the development of the N-Sink concept.  Additionally, we are grateful to Stephen Shivers, Michael Dumelle, Justin Bousquin, Joe LiVolsi, Tim Gleason, and Wayne Munns for constructive early reviews of this paper. Lastly, Ken Forshay from the US EPA's Center for Environmental Solutions and Emergency Response deserves our thanks for shepherding the development of N-Sink for many years. The views expressed in this article are those of the authors and do not necessarily represent the views or policies of the U.S. Environmental Protection Agency. Any mention of trade names, products, or services does not imply an endorsement by the U.S. Government or the U.S. Environmental Protection Agency. The EPA does not endorse any commercial products, services, or enterprises. This contribution is identified by the tracking number ORD-044618 of the Atlantic Coastal Environmental Sciences Division, Office of Research and Development, Center for Environmental Measurement and Modeling, US Environmental Protection Agency.
    
# References
v1.2.0 (2022-03-09)
===================

# Changes
- Finalized edits for JOSS publication
- 7 zip is no longer required for use.  The archive pacakge is used to extract
the NHDPlus 7z files.
- Added contributing info to README

# Bug Fixes
- Impervious used to be all NA and impervious percentages.  No zero.  Now it is no NA and ranges from zero to 100.  That broke the removal estimates.  Now fixed


v1.1.0 (2021-12-20)
===================

- First real release
- JOSS Submission
# Contributing to nsink

This outlines how to propose a change to nsink. 

## Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation directly using the GitHub web interface, as long as the changes are made in the _source_ file. 
This generally means you'll need to edit [roxygen2 comments](https://roxygen2.r-lib.org/articles/roxygen2.html) in an `.R`, not a `.Rd` file. 
You can find the `.R` file that generates the `.Rd` by reading the comment in the first line.

## Bigger changes

If you want to make a bigger change, it's a good idea to first file an issue and make sure someone from the team agrees that it’s needed. 
If you’ve found a bug, please file an issue that illustrates the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex) (this will also help you write a unit test, if needed).

### Pull request process

*   Fork the package and clone onto your computer. If you haven't done this before, we recommend using `usethis::create_from_github("jhollist/nsink", fork = TRUE)`.

*   Install all development dependencies with `devtools::install_dev_deps()`, and then make sure the package passes R CMD check by running `devtools::check()`. 
    If R CMD check doesn't pass cleanly, it's a good idea to ask for help before continuing. 
*   Create a Git branch for your pull request (PR). We recommend using `usethis::pr_init("brief-description-of-change")`.

*   Make your changes, commit to git, and then create a PR by running `usethis::pr_push()`, and following the prompts in your browser.
    The title of your PR should briefly describe the change.
    The body of your PR should contain `Fixes #issue-number`.

*  For user-facing changes, add a bullet to the top of `NEWS.md` (i.e. just below the first header). 

### Code style

*  We use [roxygen2](https://cran.r-project.org/package=roxygen2), with [Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/rd-formatting.html), for documentation.  

*  We use [testthat](https://cran.r-project.org/package=testthat) for unit tests. 
   Contributions with test cases included are easier to accept.  

## Code of Conduct

Please note that the nsink project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.
---
title: "nsink"
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

<!-- badges: start -->
[![R build status](https://github.com/jhollist/nsink/workflows/R-CMD-check/badge.svg)](https://github.com/jhollist/nsink/actions)
[![Codecov test coverage](https://codecov.io/gh/jhollist/nsink/branch/main/graph/badge.svg)](https://codecov.io/gh/jhollist/nsink?branch=main)
[![Lifecycle: stable](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html#stable)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.6341565.svg)](https://doi.org/10.5281/zenodo.6341565)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.04039/status.svg)](https://doi.org/10.21105/joss.04039)
<!-- badges: end -->

# Statement of need

The `nsink` package is an R implementation of the methods described in [Kellogg et. al (2010)](https://doi.org/10.1016/j.ecoleng.2010.02.006).  Previous implementation of this approach relied on a manual, vector based approach that was time consuming to prepare.  This approach uses a hybrid raster-vector approach that takes relatively little time to set up for each new watershed and relies on readily available data. Total run times vary, but range from minutes up to 5 hours depending on options selected.  Previous versions took weeks of manual data manipulation.  Thus, `nsink` was developed to satisfy the need for quicker implementation of the NSink method as described in [Kellogg et. al (2010)](https://doi.org/10.1016/j.ecoleng.2010.02.006).  

# `nsink` functionality
As of `r lubridate::today()` user functions for the `nsink` package are:

- `nsink_get_huc_id()`: A function for searching the name of a USGS Watershed Boundary Dataset Hydrologic Unit (<https://www.usgs.gov/core-science-systems/ngp/national-hydrography/watershed-boundary-dataset>) and retrieving its 12-digit Hydrologic Unit Code (HUC).  
- `nsink_get_data()`: Using any acceptable HUC ID (e.g. 2-digit to 12-digit), this function downloads the NHDPlus, SSURGO, NLCD Land Cover, and the NLCD Impervious for that HUC.   
- `nsink_prep_data()`: `nsink` needs data in a common coordinate reference system, from mutliple NHDPlus tables, and from different portions of SSURGO.  This function completes these data preparation steps and outputs all data, clipped to the HUC boundary.
- `nsink_calc_removal()`: Quantifying relative N removal across a landscape is a key aspects of an `nsink` analysis.  The `nsink_calc_removal()` function takes the object returned from `nsink_prep_data()` and calculates relative N removal for each landscape sink.  See Kellogg et al [-@kellogg2010geospatial] for details on relative N removal estimation for each sink.
- `nsink_generate_flowpath()`: This function uses a combination of flow determined by topography, via a flow-direction raster, for the land-based portions of a flow path and of downstream flow along the NHDPlus stream network.   
- `nsink_summarize_flowpath()`: Summarizing removal along a specified flow path requires relative N removal and a generated flow path.  This function uses these and returns a  summary of relative N removal along a flow path for each sink. 
- `nsink_generate_static_maps()`: This function analyzes N removal at the watershed scale by summarizing the results of multiple flow paths. Four static maps are returned: 1)removal efficiency; 2)loading index; 3)transport index; 4)delivery index.  Removal efficiency is a rasterized version of the `nsink_calc_removal()` output.  Loading index is N sources based on NLCD categories.   Transport index is a heat map with the cumulative relative N removal along flow paths originating from a grid of points, density set by the user, across a watershed, highlighting the gradient of downstream N retention. Delivery index is the result of multiplying the loading index and the transport index, and shows potential N delivery from different sources, taking into account the relative N removal as water moves downstream. 
- `nsink_plot()`: A function that plots each raster in the list returned from `nsink_generate_static_maps()`.   
- `nsink_build()`: One of the drivers behind the development of the `nsink` package was to provide `n-sink` analysis output that could be used more broadly (e.g. within a GIS).  The `nsink_build()` runs a complete `nsink` analysis and outputs R objects, shapefiles and/or TIFFs.
- `nsink_load()`: Essentially the inverse of the `nsink_build()` function, this function takes a folder of files, likely created by `nsink_build()`, and reads them into R.

# Installation instructions

At this time we plan on maintaining the `nsink` package as a GitHub only package and thus it won't be available directly from CRAN.  You may use the `install_github()` function from the `remotes` package to install it.  The code below will take care of installing `remotes` and installing `nsink` from the GitHub repository.

```{r, eval=FALSE}
install.packages("remotes")
remotes::install_github("usepa/nsink", dependencies = TRUE, build_vignettes = TRUE)
```

And then to load up the package:

```{r eval=FALSE}
library(nsink)
```

# Documentation and examples

All functions are documented, with examples, and that documentation may be accessed, in R, via the usual help functions.  Additionally, an introduction to the `nsink` package with a more detailed workflow is documented in a vignette.  

```{r eval=FALSE}
# Load up package
library(nsink)

# Access package level help
help(package = "nsink")

# Access the Introduction to nsink vignette
vignette("intro", package = "nsink")
```

# Contributing

If you would like to contribute to the `nsink` package, please first read the [CONTRIBUTING](.github/CONTRIBUTING.md).  In short, contributions are happily 
accepted either via suggestions in the 
[Issues](https://github.com/USEPA/nsink/issues) or via pull request.
---
title: "Introduction to `nsink`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to `nsink`}
  %\VignetteEncoding{UTF-8}
  %\VignetteEngine{knitr::rmarkdown}
bibliography: nsink.bib
csl: ecology.csl
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  eval = TRUE,
  cache = FALSE,
  tidy = FALSE,
  comment = "#>",
  warning = FALSE, 
  message = FALSE
)
load(system.file("testdata.rda", package="nsink"))
```

```{r setup, echo=FALSE}
library(nsink)
```

The package `nsink` implements an approach to estimate relative nitrogen removal
along a flow path.  This approach is detailed in Kellogg et al. 
[-@kellogg2010geospatial] and builds on peer-reviewed literature in the form of 
reviews and meta-analyses [i.e., @mayer2007meta; @alexander2007role; 
@seitzinger2006denitrification] to estimate nitrogen (N) removal within three 
types of landscape sinks -- wetlands, streams and lakes -- along any given flow 
path within a HUC12 basin. The `nsink` package implements this approach, using 
publicly available spatial data to identify flow paths and estimate N removal in
landscape sinks. Removal rates depend on retention time, which is influenced by 
physical characteristics identified using publicly available spatial data -- 
National Hydrography Dataset Plus (NHDPlus), Soil Survey Geographic Database 
(SSURGO), the National Land Cover Dataset (NLCD) land cover, and the National 
Land Cover Dataset (NLCD) impervious surface. Static maps of a specified HUC-12 
basin are generated -- N Removal Efficiency, N Loading Index, N Transport Index, 
and N Delivery Index. These maps may be used to inform local decision-making by 
highlighting areas that are more prone to N "leakiness" and areas that 
contribute to N removal.

The `nsink` package provides several functions to set up and run an N-Sink 
analysis for a specified 12-digit HUC code.  All required data are downloaded, 
prepared for the analysis, HUC-wide nitrogen removal calculated, and flow paths 
summarized.  Additionally, a convenience function that will run all of the 
required functions for a specified HUC is included.  Details on each of the steps are 
outlined in this vignette.

# Get data
The first step in the N-Sink process is to acquire data needed for the analysis.
The `nsink` package utilizes openly available data from several U.S. Federal 
Government sources.  Each dataset uses a 12-digit HUC ID number to select the 
data for download.  The first step is to identify the HUC ID and then download 
the data.

To identify the HUC ID you can use the `nsink_get_huc_id()` function which will 
use a 12-digit HUC name to search all HUCs.  Matches are returned as a data 
frame with an option to return partial or exact matches.

```{r get_huc}
# Get HUC ID - Palmer showing multiple matches
nsink_get_huc_id("Palmer")

# The Niantic
niantic_huc_id <- nsink_get_huc_id("Niantic River")$huc_12
niantic_huc_id
```

With the HUC ID in hand we can now use the `nsink_get_data()` function to 
download the required data.  All data are from publicly available sources and
as of `r lubridate::today()` no authentication is required to access these 
sources.  The HUC ID is required and users may specify a path for storing the 
data as well as indicate whether or not to download the data again if they 
already exist in the data directory.  Also note, the file archiver [7-zip](https://www.7-zip.org/) is required by `nsink_get_data()` to extract the 
NHD Plus files.


```{r get_data_norun, eval=FALSE}
# Get data for selected HUC
niantic_download <- nsink_get_data(niantic_huc_id, 
                                   data_dir = "nsink_niantic_data")
```

In addition to download the data, the function returns the basic information about your download: HUC ID and download location.


```{r get_data_return}
niantic_download
```

# Prepare the data
Once the data is downloaded there are several additional data processing steps 
that are required to subset the data just to the HUC and set all data to a 
common coordinate reference system (CRS).  

These include:

- filter out the HUC boundary
- mask all other data sets to the HUC boundary
- convert all columns names to lower case
- create new columns
- harmonize raster extents
- set all data to common CRS

The `nsink_prep_data()` function will complete all of these processing steps.  
It requires a HUC ID, a specified CRS, and a path to a data directory.  It 
returns a list with all required data for subsequent N-Sink analyses.

A quick note on the CRS.  In the near future, the preferred way to specify the CRS values will either be with Well-Known Text (WKT) or [EPSG Codes](http://spatialreference.org/ref/epsg/). Proj.4 strings will eventually be deprecated.  Currently the packages that `nsink` relies on are at different stages in implementing the changes to PROJ.  `nsink` currently works with all options, but Proj.4 strings are not recomended.  This vignette shows examples with EPSG codes.

```{r prep_data, eval=FALSE}
# EPSG for CONUS Albers Equal Area
aea <- 5072

# Prep data for selected HUC
niantic_data <- nsink_prep_data(niantic_huc_id, projection = aea, 
                data_dir = "nsink_niantic_data")
```

# Calculate removal

The next step in the N-Sink process is to calculate relative nitrogen removal.  
Details on how the nitrogen removal estimates are calculated are available in 
Kellogg et al. [-@kellogg2010geospatial]. The `nsink_calc_removal()` function 
takes the prepared data as an input and returns a list with three items:

- `raster_method`: This item contains a raster based approach to calculating 
removal.  Used for the static maps of removal.
- `land_removal`: This represents land based nitrogen removal which is hydric 
soils with areas of impervious surface removed.
- `network_removal`: This contains removal along the NHD Plus flow network.  
Removal is calculated separately for streams and waterbodies (e.g. lakes and 
reservoirs).

```{r calc_removal, eval=FALSE}
# Calculate removal from prepped data
niantic_removal <- nsink_calc_removal(niantic_data)
```

# Generate and summarize flowpaths

A useful part of the N-Sink approach is the ability to summarize that removal 
along the length of a specified flowpath.  The `nsink` package provides two 
functions that facilitate this process.  The `nsink_generate_flowpath()` 
function takes a point location as an `sf` object and the prepped data 
(generated by `nsink_prep_data()`) as input and returns an `sf` LINESTRING of 
the flowpath starting from the input point location and terminating at the 
furthest downstream location in the input NHD Plus.  The flowpath on land is 
generated from a flow direction grid.  Once that flowpath intersects the stream 
network, flow is determined by flow along the NHD Plus stream network.  First, 
create the `sf` POINT object.         

```{r flowpath_start}
# Load up the sf package
library(sf)
# Starting point
pt <- c(1948121, 2295822)
start_loc <- st_sf(st_sfc(st_point(c(pt)), crs = aea))
```

You may also determine your point location interactively by plotting your data 
and using the `locator()` function .  First create a simple plot.

```{r flowpath_plot, eval=FALSE}
# Create a simple plot
plot(st_geometry(niantic_data$huc))
plot(st_geometry(niantic_data$lakes), add = T, col = "darkblue")
plot(st_geometry(niantic_data$streams), add = T, col = "blue")
```

With the map made, you can use that to interactively select a location and use 
the x and y to create the `sf` POINT object.

```{r interactive_loc, eval=FALSE}
# Select location on map for starting point
pt <- unlist(locator(n = 1))
# convert to sf POINT
start_loc_inter <- st_sf(st_sfc(st_point(pt), crs = aea))
```

With a point identified, we can use that as the starting location for our 
flowpath.

```{r generate_flowpath, eval=FALSE}
niantic_fp <- nsink_generate_flowpath(start_loc, niantic_data)
```

The returned value has both the `flowpath_ends`, the portion of the flowpath on
the land which is created using the flow direction grid, and the 
`flowpath_network` which is the portion of the flowpath from the NHD Plus 
network that occur after the upstream `flowpath_ends` intersect the network.

```{r thefp}
niantic_fp
```

With a flowpath generated we can summarize the relative nitrogen removal along 
that flowpath with the `nsink_summarize_flowpath()` function. It takes the 
flowpath and removal as input.  A data frame is returned with each segment 
identified by type, the percent removal associated with that segment, and 
relative removal.  Total relative removal is 100 - minimum of the `n_out` 
column.

```{r summarize_it, eval=FALSE}
niantic_fp_removal <- nsink_summarize_flowpath(niantic_fp, niantic_removal)
niantic_fp_removal
100-min(niantic_fp_removal$n_out)
```
```{r summarize_show, echo=FALSE}
niantic_fp_removal
100-min(niantic_fp_removal$n_out)
```

# Static maps

Individual flow paths are useful for specific applications, but often it is more
useful to look at removal patterns across the landscape.  The 
`nsink_generate_static_maps()` function provides these HUC wide rasters.  
Required inputs are the prepped data, removal raster, and sampling density.   
The function returns four separate rasters.

- `removal_effic`: Landscape wide estimate of relative nitrogen removal 
percentage.
- `loading_idx`: An index of relative nitrogen loads by land cover class derived 
from published sources
- `transport_idx`: Relative nitrogen transport for a sample of all 
possible flowpaths in a given HUC.  This is an expensive computational task, so 
`nsink` generates a removal hotspot map based on a sample of flowpaths and the 
final hotspot map is interpolated from these samples and referred to as the 
nitrogen transport index.  The `samp_density` argument controls the number of 
sample flowpaths generated.  
- `delivery_idx`: The delivery index is the combination of the loading index and 
the transport index  It represents which areas of the landscape are 
delivering the most nitrogen to the outflow of the watershed.

```{r static_maps, eval=FALSE}
niantic_static_maps <- nsink_generate_static_maps(niantic_data, niantic_removal, 
                                                  900)
```

And with these static maps made, you can plot them quickly with 
`nsink_plot()` and sepcifying which plot you would like to see with the `map` 
argument which can be "removal", "transport", or "delivery".  
An example of `nsink_plot()` is below.

```{r plots_fake, eval=FALSE}
nsink_plot(niantic_static_maps, "transport")
```

```{r plots, echo=FALSE, warning=FALSE}
nsink_plot(niantic_static_maps, "transport")
```

# Convenience function: Build it all!

The workflow described above includes all the basic functionality.  Some users 
may wish to use `nsink` to calculate the base layers for an N-Sink analysis and 
then build an application outside of R.  A convenience function that downloads 
all data, prepares, calculates removal, and generates static maps has been 
included to facilitate this type of analysis.  The `nsink_build()` function 
requires a HUC ID, coordinate reference system, and sampling density.  An output
folder is also needed but has a default location.  Optional arguments for 
forcing a new download and playing a sound to signal when the build has 
finished. Nothing returns to R, but all prepped data files and `.tif` files are 
written into the output folder for use in other applications.


```{r build, eval=FALSE}
niantic_huc_id <- nsink_get_huc_id("Niantic River")$huc_12
aea <- 5072
nsink_build(niantic_huc_id, aea, samp_dens = 900)
```

# References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_calc_removal.R
\name{nsink_calc_removal}
\alias{nsink_calc_removal}
\title{Calculates N-Sink nitrogen removal percent}
\usage{
nsink_calc_removal(
  input_data,
  off_network_lakes = NULL,
  off_network_streams = NULL,
  off_network_canalsditches = NULL
)
}
\arguments{
\item{input_data}{A list of input datasets created with
\code{\link{nsink_prep_data}}.}

\item{off_network_lakes}{Optional argument to set removal for waterbodies
that are not part of the NHDPlus hydrologic network. Default value is to
use the 75th percentile of removal from other lakes in the HUC.  If another
value is desired provide a single numeric ranging from 0 to 1.}

\item{off_network_streams}{Optional argument to set removal for streams that
are not part of the NHDPlus hydrologic network. Default value is to use
the median removal from first order streams in the HUC.  If another value
is desired provide a single numeric ranging from 0 to 1.}

\item{off_network_canalsditches}{Optional argument to set removal for canals
and ditches that are not part of the NHDPlus hydrologic network. Default
value is to use the 25th percentile of removal from third order streams in
the HUC. If another value is desired provide a single numeric ranging from
0 to 1.}
}
\value{
A list with three items, 1) a raster stack with one layer with
        the nitrogen removal, a second layer with the type of removal (e.g.
        hydric soils, lakes, streams), 2) a polygon representing removal from
        land, and 3) a polygon representing removal from the stream network,
        including stream removal, and lake removal.
}
\description{
Starting with base data layers of NHDPlus, SSURGO, impervious surface, flow
velocity, and time of travel, this function calculates percentage of Nitrogen
removal.  Nitrogen removal methods and calculations are from
\href{https://doi.org/10.1016/j.ecoleng.2010.02.006}{Kellogg et al. (2010)}.
This function assumes data has been downloaded with
\code{\link{nsink_get_data}} and has been prepared with
\code{\link{nsink_prep_data}}.
}
\examples{
\dontrun{
library(nsink)
niantic_huc <- nsink_get_huc_id("Niantic River")$huc_12
niantic_data <- nsink_get_data(niantic_huc, data_dir = "nsink_data")
aea <- 5072
niantic_nsink_data <- nsink_prep_data(niantic_huc, projection = aea ,
                                      data_dir = "nsink_data")
removal <- nsink_calc_removal(niantic_nsink_data)
}
}
\references{
Kellogg, D. Q., Gold, A. J., Cox, S., Addy, K., & August, P. V.
            (2010). A geospatial approach for assessing denitrification sinks
            within lower-order catchments. Ecological Engineering, 36(11),
            1596-1606.
            \href{https://doi.org/10.1016/j.ecoleng.2010.02.006}{Link}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_generate_flowpath.R
\name{nsink_generate_flowpath}
\alias{nsink_generate_flowpath}
\title{Generate and clean a flowpath for N-Sink}
\usage{
nsink_generate_flowpath(starting_location, input_data)
}
\arguments{
\item{starting_location}{An \code{\link{sf}} point location as a starting
point for the flowpath.  Projection must match
projection in input_data.}

\item{input_data}{A list of input data with (at least) "fdr", "streams",
"tot", and "raster_template". These may be generated with
\code{\link{nsink_prep_data}}.}
}
\value{
An \code{\link{sf}} LINESTRING object of the flowpath that starts at
        the \code{starting_location} and ends at the ouflow of the HUC.
}
\description{
This function takes an XY location as a starting point and generates a
flowpath for use in the N-Sink nitrogen removal analysis. The flowpath is a
combination of a flow direction derived flowpath on land plus NHDPlus derived
stream-reach flowpath.
}
\examples{
\dontrun{
library(nsink)
niantic_huc <- nsink_get_huc_id("Niantic River")$huc_12
niantic_data <- nsink_get_data(niantic_huc, data_dir = "nsink_data")
aea <- 5072
niantic_nsink_data <- nsink_prep_data(niantic_huc, projection = aea,
                                      data_dir = "nsink_niantic_data")
pt <- c(1948121, 2295822)
start_loc <- st_sf(st_sfc(st_point(c(pt)), crs = aea))
fp <- nsink_generate_flowpath(start_loc, niantic_nsink_data)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_fdr}
\alias{nsink_prep_fdr}
\title{Prepare flow direction data for N-Sink}
\usage{
nsink_prep_fdr(huc_sf, huc_raster, data_dir)
}
\arguments{
\item{huc_sf}{An sf object of the Watershed Boundaries Dataset HUC12}

\item{huc_raster}{A raster object of the Watershed Boundaries Dataset HUC12}

\item{data_dir}{Base directory that contains N-Sink data folders.  Data may be
downloaded with the \code{\link{nsink_get_data}} function.}
}
\value{
A raster object of the flow direction for the huc_sf but in
        the original fdr projection
}
\description{
Standardizes flow direction data by transforming data, and clipping to HUC.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_build.R
\name{nsink_write_static_maps}
\alias{nsink_write_static_maps}
\title{Write static maps to files}
\usage{
nsink_write_static_maps(static_maps, output_dir)
}
\arguments{
\item{static_maps}{A list of static maps, as output by
\code{\link{nsink_generate_static_maps}}}

\item{output_dir}{Output folder to save .tif static maps to}
}
\description{
Writes out static maps as tiffs
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_ssurgo}
\alias{nsink_prep_ssurgo}
\title{Prepare SSURGO data for N-Sink}
\usage{
nsink_prep_ssurgo(huc_sf, data_dir)
}
\arguments{
\item{huc_sf}{An sf object of the Watershed Boundaries Dataset HUC12}

\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}
}
\value{
An sf object of the SSURGO data for the huc_sf with hydric data added.
}
\description{
Standardizes impervious data by transforming data, and reducing columns to
what is needed.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{nsink_fix_data_directory}
\alias{nsink_fix_data_directory}
\title{Fix the data directory}
\usage{
nsink_fix_data_directory(data_dir)
}
\arguments{
\item{data_dir}{The data directory}
}
\value{
A string with the normalized path
}
\description{
This function takes the data directory and checks for existence, creates it
if it doesn't exist, then adds a trailing slash and normalizes the path for
the operating system
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_plot.R
\name{nsink_plot_transport}
\alias{nsink_plot_transport}
\title{N-Sink plot function for the nitrogen transport index}
\usage{
nsink_plot_transport(
  transport_idx,
  breaks = c(0, 10, 20, 30, 40, 50, 60, 70, 80, 90, 100),
  colors = c("#38A1D0", "#7AB4C0", "#A2C8B0", "#C3DC9D", "#E2F088", "#F5E871",
    "#F9C159", "#F99844", "#F56A2E", "#EF2820")
)
}
\arguments{
\item{transport_idx}{A transport index \code{raster} that can be created
via \code{\link{nsink_generate_static_maps}}.}

\item{breaks}{A vector of values specifying breakpoints to break the
transport \code{raster}.  Must be one more than number of
colors.}

\item{colors}{A vector of hexcodes for the colors}
}
\description{
This function creates a simple plot with a pre-selected color palette,
although different breaks and colors can be provided by the user.  This is
meant as a quick means to visualize removal efficiency.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{nsink_get_plus_remotepath}
\alias{nsink_get_plus_remotepath}
\title{Get remote path for NHDPlus components}
\usage{
nsink_get_plus_remotepath(
  rpu,
  component = c("NHDSnapshot", "FdrFac", "EROMExtension", "WBDSnapshot",
    "NHDPlusAttributes")
)
}
\arguments{
\item{rpu}{raster processing unit for NHDPlus. available in nsink:::wbd_lookup}

\item{component}{which component to download}
}
\description{
Code modified from https://github.com/jsta/nhdR/blob/master/R/utils.R.  Still
need to figure out best way to acknowledge jsta as author and include GPL
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_streams}
\alias{nsink_prep_streams}
\title{Prepare streams data for N-Sink}
\usage{
nsink_prep_streams(huc_sf, data_dir)
}
\arguments{
\item{huc_sf}{An sf object of the Watershed Boundaries Dataset HUC12}

\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}
}
\value{
returns an sf object of the NHDPlus streams for the huc_sf
}
\description{
Standardizes streams data by transforming data, clipping to HUC, ...
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_plot.R
\name{nsink_plot_removal}
\alias{nsink_plot_removal}
\title{N-Sink plot function for removal efficiency}
\usage{
nsink_plot_removal(
  removal_effic,
  breaks = c(0.2, 0.4, 0.6, 0.8),
  colors = c("#D3FFBE", "#70A800", "#267300")
)
}
\arguments{
\item{removal_effic}{A removal efficiency \code{raster} that can be created
via \code{\link{nsink_generate_static_maps}}.}

\item{breaks}{A vector of values specifying breakpoints to break the
removal \code{raster}.  Must be one more than number of
colors.}

\item{colors}{A vector of hexcodes for the colors}
}
\description{
This function creates a simple plot with a pre-selected color palette,
although different breaks and colors can be provided by the user.  This is
meant as a quick means to visualize removal efficiency.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_q}
\alias{nsink_prep_q}
\title{Prepare flow data for N-Sink}
\usage{
nsink_prep_q(data_dir)
}
\arguments{
\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}
}
\value{
A tibble of the flow data
}
\description{
Standardizes flow data from the EROM tables.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_plot.R
\name{nsink_plot}
\alias{nsink_plot}
\title{N-Sink plot function for the a list of nsink static maps}
\usage{
nsink_plot(static_maps, map = c("removal", "transport", "delivery"))
}
\arguments{
\item{static_maps}{A list of \code{raster} objects that can be created
via \code{\link{nsink_generate_static_maps}}.  The list
should have removal_effic, delivery_idx, and
transport_idx.}

\item{map}{A character of either, "removal", "transport", or "delivery."}
}
\description{
This function creates a simple plot with pre-selected color palettes from a
list of static maps created by \code{\link{nsink_generate_static_maps}}.
This is meant as a quick means to visualize the various static maps.
}
\examples{
\dontrun{
library(nsink)
niantic_huc <- nsink_get_huc_id("Niantic River")$huc_12
niantic_data <- nsink_get_data(niantic_huc, data_dir = "nsink_data")
aea <- 5072
niantic_nsink_data <- nsink_prep_data(niantic_huc, projection = aea,
                                      data_dir = "nsink_data")
removal <- nsink_calc_removal(niantic_nsink_data)
static_maps <- nsink_generate_static_maps(niantic_nsink_data, removal,
samp_dens = 900)
nsink_plot(static_maps, "transport")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_plot.R
\name{nsink_plot_delivery}
\alias{nsink_plot_delivery}
\title{N-Sink plot function for the nitrogen delivery index}
\usage{
nsink_plot_delivery(
  delivery_idx,
  breaks = c(20, 40, 60, 80, 100),
  colors = c("#FFBEBE", "#F57A7A", "#A80405", "#652600")
)
}
\arguments{
\item{delivery_idx}{A delivery index \code{raster} that can be created
via \code{\link{nsink_generate_static_maps}}.}

\item{breaks}{A vector of values specifying breakpoints to break the
delivery \code{raster}.  Must be one more than number of
colors.}

\item{colors}{A vector of hexcodes for the colors}
}
\description{
This function creates a simple plot with a pre-selected color palette,
although different breaks and colors can be provided by the user.  This is
meant as a quick means to visualize delivery index.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_summarize_flowpath.R
\name{nsink_generate_from_to_nodes}
\alias{nsink_generate_from_to_nodes}
\title{nsink generate from land to nodes}
\usage{
nsink_generate_from_to_nodes(land_off_network)
}
\arguments{
\item{land_off_network}{The land off the network path}
}
\description{
Code borrowed from https://www.r-spatial.org/r/2019/09/26/spatial-networks.html
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_lakes}
\alias{nsink_prep_lakes}
\title{Prepare lakes data for N-Sink}
\usage{
nsink_prep_lakes(huc_sf, data_dir)
}
\arguments{
\item{huc_sf}{An sf object of the Watershed Boundaries Dataset HUC12}

\item{data_dir}{Base directory that contains N-Sink data folders.  Data may be
downloaded with the \code{\link{nsink_get_data}} function.}
}
\value{
An sf object of the NHDPlus lakes for the huc_sf
}
\description{
Standardizes lakes data by transforming data, clipping to HUC, and
renaming columns.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_generate_static_maps.R
\name{nsink_generate_n_removal_heatmap}
\alias{nsink_generate_n_removal_heatmap}
\title{Generate Nitrogen Removal Heatmap}
\usage{
nsink_generate_n_removal_heatmap(
  input_data,
  removal,
  samp_dens,
  ncpu = future::availableCores() - 1,
  seed = 23
)
}
\arguments{
\item{input_data}{A list of input datasets created with
\code{\link{nsink_prep_data}}.}

\item{removal}{The removal raster stack or removal list, generated by
\code{\link{nsink_calc_removal}}}

\item{samp_dens}{A value, in the units of the input data, divided by total area of
the input HUC.  It is used to determine the number of points,
determined through a regular sample, to calculate removal.  For
instance, a value of 90 would roughly equate to a point per every
90 meters.}

\item{ncpu}{Number of CPUs to use for calculating flowpath removal for larger
(i.e. greater than 50) number of flowpaths.  Default is the number
of cores available minus one.}

\item{seed}{Random seed to ensure reproducibility of sample point creation
across runs.  Default set to 23.}
}
\description{
Generates the heatmap.  Only flowpaths that have one segment or more
completely contained within the HUC are included. Practically, this will only
remove flowpaths that do not intersect the stream network and extend beyond
the HUC.  This is rare.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_remove_openwater}
\alias{nsink_remove_openwater}
\title{Remove open water portions of the HUC}
\usage{
nsink_remove_openwater(huc_sf, data_dir)
}
\arguments{
\item{huc_sf}{An sf object of the Watershed Boundaries Dataset HUC12}

\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}
}
\value{
An sf object of the HUC without salt water/open water area.
}
\description{
Uses SSURGO polys and the SSURGO musym = Ws to check for and remove any
portion of a HUC that is actually salt water.  This should account for the
coastal HUCs with large areas of open water.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_build.R
\name{nsink_build}
\alias{nsink_build}
\title{Build out required datasets for N-Sink}
\usage{
nsink_build(
  huc,
  projection,
  output_dir = normalizePath("nsink_output", winslash = "/"),
  data_dir = normalizePath("nsink_data", winslash = "/"),
  force = FALSE,
  samp_dens = 300,
  year = "2016",
  ...
)
}
\arguments{
\item{huc}{A character with the 12 digit HUC ID.  May be searched with
\code{\link{nsink_get_huc_id}}}

\item{projection}{Projection to use for all spatial data, specified as either an
EPSG code (as numeric) or WKT (as string).}

\item{output_dir}{Folder to write processed nsink files to.
Currently, the processed files will be overwritten if
the same output folder is used.  To run different
HUC12's specify separate output folders.}

\item{data_dir}{Folder to hold downloaded data.  The same data
directory can be used to hold data for multiple HUCs.  Data
will not be downloaded again if it already exists in this
folder.}

\item{force}{Logical value used to force a new download if data already
exists on file system.}

\item{samp_dens}{The \code{samp_dens} controls the density of points to use when
creating the nitrogen removal heat map.  The area of the
watershed is sampled with points that are separated by the
\code{samp_dens} value, in the units of the input data.
The larger the value, the fewer the points.}

\item{year}{Year argument to be passed to FedData's \code{\link{get_nlcd}}
function. Defaults to 2016.}

\item{...}{Passes to \code{\link{nsink_calc_removal}} for the off network
arguments: \code{off_network_lakes}, \code{off_network_streams},
and \code{off_network_canalsditches}.}
}
\value{
A list providing details on the huc used and the output location of
        the dataset.
}
\description{
This function is a wrapper around the other functions and runs all of those
required to build out the full dataset needed for a huc and develops the four
static N-Sink maps: the nitrogen loading index, nitrogen removal effeciency,
nitrogen transport index, and the nitrogen delivery index.  The primary
purpose of this is to use the nsink package to develop the required datasets
for an nsink application to be built outside of R (e.g. ArcGIS).  This will
take some time to complete as it is downloading 500-600 Mb of data,
processing that data and then creating output files.
}
\examples{
\dontrun{
library(nsink)
aea <- 5072
nsink_build(nsink_get_huc_id("Niantic River")$huc_12, aea,
            output_dir = "nsink_output", data_dir = "nsink_data",
             samp_dens = 600)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{get_nhd_plus}
\alias{get_nhd_plus}
\title{Function to download nhd files}
\usage{
get_nhd_plus(
  download_url,
  data_dir = normalizePath("nsink_data/", winslash = "/"),
  download_again = FALSE
)
}
\arguments{
\item{download_url}{url to download}

\item{data_dir}{The data dir}

\item{force}{Force new download}
}
\description{
Function to download nhd files
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_calc_removal.R
\name{nsink_calc_off_network_removal}
\alias{nsink_calc_off_network_removal}
\title{Calculates off network waterbody and stream nitrogen removal}
\usage{
nsink_calc_off_network_removal(
  input_data,
  off_network_lakes,
  off_network_streams,
  off_network_canalsditches
)
}
\arguments{
\item{input_data}{A named list with "streams", "lakes",
"network_removal", "tot", and "raster_template".}

\item{off_network_lakes}{Optional argument to set removal for waterbodies
that are not part of the NHDPlus hydrologic network. Default value is to
use the 75th percentile of removal from other lakes in the HUC.  If another
value is desired provide a single numeric ranging from 0 to 1.}

\item{off_network_streams}{Optional argument to set removal for streams that
are not part of the NHDPlus hydrologic network. Default value is to use
the median removal from first order streams in the HUC.  If another value
is desired provide a single numeric ranging from 0 to 1.}

\item{off_network_canalsditches}{Optional argument to set removal for canals
and ditches that are not part of the NHDPlus hydrologic network. Default
value is to use the 25th percentile of removal from third order streams in
the HUC. If another value is desired provide a single numeric ranging from
0 to 1.}
}
\value{
Raster and vectors of off network nitrogen removal
}
\description{
Calculates off network waterbody and stream nitrogen removal
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_calc_removal.R
\name{nsink_calc_stream_removal}
\alias{nsink_calc_stream_removal}
\title{Calculates stream-based nitrogen removal}
\usage{
nsink_calc_stream_removal(input_data)
}
\arguments{
\item{input_data}{A named list with "streams", "q", "tot", and
"raster_template".}
}
\value{
Raster and vector versions of stream based nitrogen removal
}
\description{
Calculates stream-based nitrogen removal
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_load.R
\name{nsink_load}
\alias{nsink_load}
\title{Load an existing N-Sink analysis folder}
\usage{
nsink_load(input_folder, base_name = "nsink_", projection = NULL, ...)
}
\arguments{
\item{input_folder}{Folder that contains nsink files produced by
\code{\link{nsink_build}}}

\item{base_name}{A base name used to assign objects to the global environment.}

\item{projection}{An optional CRS specified as a either an
EPSG code (as numeric) or WKT (as string).
Useful if projection is returned as unknown.}

\item{...}{Passes to \code{\link{nsink_calc_removal}} for the off network
arguments: \code{off_network_lakes}, \code{off_network_streams},
and \code{off_network_canalsditches}.}
}
\value{
Creates several lists in the global environment that would normally
        be created when running an N-Sink analysis.  These include:
        a \code{\link{nsink_prep_data}} object,
        a \code{\link{nsink_calc_removal}} object, and a
        \code{\link{nsink_generate_static_maps}} object
}
\description{
Load an existing N-Sink analysis folder
}
\examples{
\dontrun{
library(nsink)

aea <- 5072
nsink_build(nsink_get_huc_id("Niantic River")$huc_12, aea,
            output_folder = "nsink_output", samp_dens = 300)
nsink_load(input_folder = "nsink_output",
           base_name = "nsink_")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_data}
\alias{nsink_prep_data}
\title{Prepares N-Sink data for a given HUC}
\usage{
nsink_prep_data(
  huc,
  projection,
  data_dir = normalizePath("nsink_data/", winslash = "/"),
  year = "2016"
)
}
\arguments{
\item{huc}{A character string of the HUC12 ID.  Use
\code{\link{nsink_get_huc_id}} to look up ID by name.}

\item{projection}{CRS to use, passed as ethier EPSG code (as numeric)
or WKT (as character).
This must be a projected CRS and not geographic as many of
the measurements required for the nsink analysis require
reliable length and area measurments.}

\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}

\item{year}{The year of the nlcd and impervious data that was retrieved with
FedData's, \code{\link{get_nlcd}} function.}
}
\value{
returns a list of sf, raster, or tabular objects for each of the
        required datasets plus the huc.
}
\description{
In addition to having local access to the required dataset, those datasets
need to have some preparation.  This function standardizes projections and
extents and clips all datasets to the boundary of the specified HUC.
Additionally, any tabular datasets (e.g. flow, time of travel, etc.) are
included in the output as well.
}
\examples{
\dontrun{
library(nsink)
aea <- 5072
niantic_huc <- nsink_get_huc_id("Niantic River")$huc_12
niantic_nsink_data <- nsink_prep_data(huc = niantic_huc, projection = aea,
data_dir = "nsink_data")
# Example using EPSG code for projection
epsg <- 3748L
niantic_nsink_data <- nsink_prep_data(huc = niantic_huc, projection = epsg,
                data_dir = "nsink_data")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_lakemorpho}
\alias{nsink_prep_lakemorpho}
\title{Prepare lake morphology data for N-Sink}
\usage{
nsink_prep_lakemorpho(data_dir)
}
\arguments{
\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}
}
\value{
A tibble of the lake morphology data
}
\description{
Standardizes lake morphology from the lake morphology tables.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_get_data.R
\name{nsink_get_huc_id}
\alias{nsink_get_huc_id}
\title{Look up HUC 12 ID from a HUC name}
\usage{
nsink_get_huc_id(huc_name, exact = FALSE)
}
\arguments{
\item{huc_name}{Character string of a HUC Name or partial HUC name}

\item{exact}{Logical indicating whether or not to do an exact match}
}
\value{
A data frame with HUC_12 and HU_12_NAME that match the huc_name
}
\description{
This function takes a HUC Name and returns matching HUC 12 IDs.  The default
behavior is to select all possible matching IDs without matching the case of
the string.  If an exact match is required, use the  \code{exact} argument.
}
\examples{
nsink_get_huc_id(huc_name = "Niantic River")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_generate_flowpath.R
\name{nsink_get_flowline}
\alias{nsink_get_flowline}
\title{Get flowlines that intersect with a flowpath}
\usage{
nsink_get_flowline(flowpath_ends, streams, tot)
}
\arguments{
\item{flowpath_ends}{An \code{sf} LINESTRING of the flowpath ends, generated
with \code{\link{nsink_get_flowpath_ends}}}

\item{streams}{NHDPlus streams from \code{\link{nsink_prep_data}}}

\item{tot}{NHDPlus time of travel from \code{\link{nsink_prep_data}} which
provides the from and to nodes.}
}
\value{
An \code{sf} object of the NHDPlus flowlines that occur after a
        raster flowpath intersects the stream network.
}
\description{
Extract flowlines that intersect with flowpath ends.  This uses the actual
flowlines as a part of the flowpath instead of simply using the raster
derived flowpaths which do not follow the flowlines exactly.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_summarize_flowpath.R
\name{nsink_create_summary}
\alias{nsink_create_summary}
\title{Create a nitrogen removal summary}
\usage{
nsink_create_summary(land_removal, network_removal)
}
\arguments{
\item{land_removal}{A data frame of land based nitrogen removal via the
generated flowpath and hydric removal raster.}

\item{network_removal}{A data frame of stream network nitrogen removal via
calculated stream and lake removal.}
}
\value{
A data frame summarizing nitrogen removal along a flowpath
}
\description{
This functions takes a nitrogen removal and removal type data frame for land
and the flowpath network and creates a data frame that summarizes the
removal along that flowpath.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{nsink_run_7z}
\alias{nsink_run_7z}
\title{Finds 7-zip}
\usage{
nsink_run_7z(zipfile, destdir, extract_again = FALSE)
}
\arguments{
\item{zipfile}{The zipfile to be extracted}

\item{destdir}{Where to put the extracted files}

\item{force}{Whether or not to extract again if the destination files
already exist}
}
\description{
This code is modified from https://github.com/jsta/nhdR/blob/master/R/utils.R
and https://github.com/jsta/nhdR/blob/master/R/get.R to determine if 7 zip is
available.  If available it unzips a 7z zipfile to a destination directory.
This avoids needing to use archive package which is only available via
GitHub.
}
\author{
Joseph Stachelek, \email{stachel2@msu.edu}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_nlcd}
\alias{nsink_prep_nlcd}
\title{Prepare NLCD data for N-Sink}
\usage{
nsink_prep_nlcd(huc_sf, huc_raster, data_dir, year)
}
\arguments{
\item{huc_sf}{An sf object of the Watershed Boundaries Dataset HUC12}

\item{huc_raster}{A raster object of the Watershed Boundaries Dataset HUC12}

\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}

\item{year}{The year of the nlcd and impervious data that was retrieved with
FedData's, \code{\link{get_nlcd}} function.}
}
\value{
A raster object of the NLCD for the huc_sf
}
\description{
Standardizes NLCD data by by projecting to the coorect CRS,
and standardizing file name.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_calc_removal.R
\name{nsink_calc_land_removal}
\alias{nsink_calc_land_removal}
\title{Calculates land-based nitrogen removal}
\usage{
nsink_calc_land_removal(input_data)
}
\arguments{
\item{input_data}{A named list with "ssurgo", "impervious", "lakes", and
"raster_template".}
}
\value{
List with raster and vector versions of land based nitrogen removal
}
\description{
Calculates land-based nitrogen removal
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_impervious}
\alias{nsink_prep_impervious}
\title{Prepare impervious cover data for N-Sink}
\usage{
nsink_prep_impervious(huc_sf, huc_raster, data_dir, year)
}
\arguments{
\item{huc_sf}{An sf object of the Watershed Boundaries Dataset HUC12}

\item{huc_raster}{A raster object of the Watershed Boundaries Dataset HUC12}

\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}

\item{year}{The year of the nlcd and impervious data that was retrieved with
FedData's, \code{\link{get_nlcd}} function.}
}
\value{
A raster object of the impervious cover for the huc_sf
}
\description{
Standardizes impervious data by projecting to the coorect CRS,
and standardizing file name.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink.R
\name{nsink}
\alias{nsink}
\title{nsink: A package that implements flow path analysis of nitrogen removal}
\description{
The N-Sink approach is based off of research outlined in
\href{https://doi.org/10.1016/j.ecoleng.2010.02.006}{Kellogg et al (2010)}.
This approach builds on peer-reviewed literature in the form of reviews and
meta-analyses (i.e., \href{https://doi.org/doi:10.2134/jeq2006.0462}{Mayer et
al (2007)}, \href{https://doi.org/10.1111/j.1752-1688.2007.00005.x}{Alexander
et al (2007)}, and
\href{https://doi.org/10.1890/1051-0761(2006)016[2064:DALAWA]2.0.CO;2}{Seitzinger
et al (2006)}) to estimate nitrogen (N) removal within three types of
landscape sinks -- wetlands, streams and lakes -- along any given flow path
within a HUC12 basin. The \code{nsink} package implements this approach,
using publicly available spatial data to identify flow paths and estimate N
removal in landscape sinks. Removal rates depend on retention time, which is
influenced by physical characteristics identified using publicly available
spatial data -- National Hydrography Dataset (NHD), Watershed Boundary
Dataset (WBD), National Land Cover Dataset (NLCD), and Soil Survey Geographic
Dataset (SSURGO). Static maps of a specified HUC-12 basin are generated -- N
Removal Efficiency, N Loading Index, N Transport Index, and N Delivery Index.
These maps may be used to inform local decision-making by highlighting areas
that are more prone to N "leakiness" and areas that contribute to N removal.
}
\references{
Kellogg, D. Q., Gold, A. J., Cox, S., Addy, K., & August, P. V.
            (2010). A geospatial approach for assessing denitrification sinks
            within lower-order catchments. Ecological Engineering, 36(11),
            1596-1606.
            \href{https://doi.org/10.1016/j.ecoleng.2010.02.006}{Link}

            Mayer, P. M., Reynolds Jr., S. K., McCutchen, M. D., & Canfield, T. J.
            (2007). Meta-analysis of nitrogen removal in riparian buffers.
            J. Environ. Qual. 36, 1172-1180.
            \href{https://doi.org/doi:10.2134/jeq2006.0462}{Link}

            Alexander, R. B., Boyer, E. W., Smith, R.A., Schwarz, G.E. & Moore, R. B.
            (2007). The role of headwater streams in downstream water quality.
            J. Am. Water Resou.Assoc. 43, 41-59.
            \href{https://doi.org/10.1111/j.1752-1688.2007.00005.x}{Link}

            Seitzinger, S. P., Harrison, J.A., Hohlke, J.K., Bouwman, A.F., Lowrance, R.,
            Peterson, B., Tobias, C., & Van Drecht, G. (2006). Denitrification across
            landscapes and waterscapes: a synthesis. Ecol. Appl. 16, 1064-2090.
            \href{https://doi.org/10.1890/1051-0761(2006)016[2064:DALAWA]2.0.CO;2}{Link}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_calc_removal.R
\name{nsink_merge_removal}
\alias{nsink_merge_removal}
\title{Merges removal rasters into single raster}
\usage{
nsink_merge_removal(removal_rasters)
}
\arguments{
\item{removal_rasters}{A named list of "land_removal", "stream_removal",
"lake_removal", and "raster_template" rasters plus a
sf object "huc".}
}
\value{
Raster of landscape nitrogen removal
}
\description{
Merges removal rasters into single raster
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_generate_static_maps.R
\name{nsink_generate_static_maps}
\alias{nsink_generate_static_maps}
\title{Generate N-Sink Static Maps for a given HUC}
\usage{
nsink_generate_static_maps(
  input_data,
  removal,
  samp_dens,
  ncpu = future::availableCores() - 1,
  seed = 23
)
}
\arguments{
\item{input_data}{A list of input datasets created with
\code{\link{nsink_prep_data}}.}

\item{removal}{The removal raster stack or removal list, generated by
\code{\link{nsink_calc_removal}}}

\item{samp_dens}{A value, in the units of the input data, divided by total area of
the input HUC.  It is used to determine the number of points,
determined through a regular sample, to calculate removal.  For
instance, a value of 90 would roughly equate to a point per every
90 meters.}

\item{ncpu}{Number of CPUs to use for calculating flowpath removal for larger
(i.e. greater than 50) number of flowpaths.  Default is the number
of cores available minus one.}

\item{seed}{Random seed to ensure reproducibility of sample point creation
for transport maps. Default set to 23.}
}
\value{
This function returns a list of rasters: nitrogen removal efficiency,
        nitrogen loading index, nitrogen transport index, and the
        nitrogen delivery index.
}
\description{
Generate N-Sink Static Maps for a given HUC
}
\examples{
\dontrun{
library(nsink)
niantic_huc <- nsink_get_huc_id("Niantic River")$huc_12
niantic_data <- nsink_get_data(niantic_huc, data_dir = "nsink_data")
aea <- 5072
niantic_nsink_data <- nsink_prep_data(niantic_huc, projection = aea,
                                      data_dir = "nsink_data")
removal <- nsink_calc_removal(niantic_nsink_data)
static_maps <- nsink_generate_static_maps(niantic_nsink_data, removal,samp_dens = 900)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_calc_removal.R
\name{nsink_calc_lake_removal}
\alias{nsink_calc_lake_removal}
\title{Calculates lake-based nitrogen removal}
\usage{
nsink_calc_lake_removal(input_data)
}
\arguments{
\item{input_data}{A named list with "streams", "lakes", "tot", "lakemorpho",
and "raster_template".}
}
\value{
Raster and vector versions of lake based nitrogen removal
}
\description{
Calculates lake-based nitrogen removal
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{nsink_get_closest}
\alias{nsink_get_closest}
\title{Get closest}
\usage{
nsink_get_closest(v1, v2)
}
\arguments{
\item{v1}{The first vector of values, likely length or area of a feature}

\item{v2}{The second vector of values, also likely length or area of a
feature.}
}
\value{
Returns a vector, of length v1, of the index from v2 that is closest
        in absolute value to each value in v1.
}
\description{
This function will return the index of one vector that is closest, by
absolute values, to the values in another vector.  Usually used to
identify the closest lake or stream, in area or length respectively, to
another lake or stream.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{nsink_get_closest_lt}
\alias{nsink_get_closest_lt}
\title{Get closest, but less than}
\usage{
nsink_get_closest_lt(v1, v2)
}
\arguments{
\item{v1}{The first vector of values, likely length or area of a feature}

\item{v2}{The second vector of values, also likely length or area of a
feature.}
}
\value{
Returns a vector, of length v1, of the index from v2 that is less
        than and  closest in absolute value to each value in v1.
}
\description{
This function will return the index of one vector that is closest and less
than, by absolute values, to the values in another vector.  If the value in
v1 is less than all values in v2, then just the next closest is returned.
Usually used to identify the closest lake or stream, in area or length
respectively, to another lake or stream.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_prep_data.R
\name{nsink_prep_tot}
\alias{nsink_prep_tot}
\title{Prepare time of travel data for N-Sink}
\usage{
nsink_prep_tot(data_dir)
}
\arguments{
\item{data_dir}{Base directory that contains N-Sink data folders.  Data may
be downloaded with the \code{\link{nsink_get_data}} function.}
}
\value{
A tibble of the time of travel data
}
\description{
Standardizes time of travel from the NHDPlus VAA tables.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_summarize_flowpath.R
\name{nsink_group_land_off_network}
\alias{nsink_group_land_off_network}
\title{nsink group land off network}
\usage{
nsink_group_land_off_network(land_removal_df)
}
\arguments{
\item{land_removal_df}{Land and off network summary df}
}
\description{
nsink group land off network
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_build.R
\name{nsink_write_prepped_data}
\alias{nsink_write_prepped_data}
\title{Write prepped data to files}
\usage{
nsink_write_prepped_data(prepped_data, output_dir)
}
\arguments{
\item{prepped_data}{A list of prepped data, as output by
\code{\link{nsink_prep_data}}}

\item{output_dir}{Output folder to save processed nsink files to}
}
\description{
Writes out data either as shapefiles, for vector data, or tiffs for raster
data.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_summarize_flowpath.R
\name{nsink_create_segment_ids}
\alias{nsink_create_segment_ids}
\title{Create ID's for unique values in a vector}
\usage{
nsink_create_segment_ids(x)
}
\arguments{
\item{x}{A vector of values to create unique ID's}
}
\value{
A vector of unique ID's
}
\description{
This functions takes a vector of values and creates a unique ID for adjacent
and identical values.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_generate_flowpath.R
\name{nsink_get_flowpath_ends}
\alias{nsink_get_flowpath_ends}
\title{Get flowpath beginning and ends}
\usage{
nsink_get_flowpath_ends(flowpath, streams, tot)
}
\arguments{
\item{flowpath}{An \code{sf} LINESTRING of the flowpath, generated with \code{\link{nsink_generate_flowpath}}}

\item{streams}{NHDPlus streams from \code{\link{nsink_prep_data}}}

\item{tot}{NHDPlus time of travel from \code{\link{nsink_prep_data}} which
provides the from and to nodes.}
}
\value{
An \code{sf} object of the portions of the flowpath that are not
        represented by the NHDPlus flowlines
}
\description{
Flowpath from land is only portion that needs to be generated from the flow
direction grid.  This function extracts those portions of the generated
flowpath.  There may be issues with off network waterbodies...  May need to
find all sections without connected flowpaths...
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_generate_static_maps.R
\name{nsink_generate_n_loading_index}
\alias{nsink_generate_n_loading_index}
\title{Generates the Nitrogen Loading Index}
\usage{
nsink_generate_n_loading_index(input_data)
}
\arguments{
\item{input_data}{List of input data}
}
\description{
This function reclassifies the NLCD data to an index of Nitrogen loading.
The index ranges from 0 to 1.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_generate_flowpath.R
\name{nsink_split_flowline}
\alias{nsink_split_flowline}
\title{Split flowlines where they intersect with a flowpath}
\usage{
nsink_split_flowline(flowpath_ends, flowpath_network)
}
\arguments{
\item{flowpath_ends}{The ends of the flowpath that are not a part of the
network}

\item{flowpath_network}{The flowpath network}
}
\value{
An \code{sf} object of the NHDPlus flowlines split where the fp ends
        intersect with the flowlines.
}
\description{
Takes a flowpath input
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_summarize_flowpath.R
\name{nsink_summarize_flowpath}
\alias{nsink_summarize_flowpath}
\title{Summarize nitrogen removal along a flowpath}
\usage{
nsink_summarize_flowpath(flowpath, removal)
}
\arguments{
\item{flowpath}{A flowpath to summarize nitrogen removal}

\item{removal}{The removal raster stack or removal list, generated by
\code{\link{nsink_calc_removal}}}
}
\value{
A data frame is returned with a summary of nitrogen removal.  The
}
\description{
Nitrogen removal varies along a flowpath as it may include different land
cover and waterbody types that have different nitrogen reduction capabilities.
This function requires a flowpath generated by
\code{\link{nsink_generate_flowpath}} as input and returns an
estimate of total flow path removal as well as removal by type.
}
\examples{
\dontrun{
library(nsink)
niantic_huc <- nsink_get_huc_id("Niantic River")$huc_12
niantic_data <- nsink_get_data(niantic_huc, data_dir = "nsink_data")
aea <- 5072
niantic_nsink_data <- nsink_prep_data(niantic_huc, projection = aea,
                                      data_dir = "nsink_data")
removal <- nsink_calc_removal(niantic_nsink_data)
pt <- c(1948121, 2295822)
start_loc <- st_sf(st_sfc(st_point(c(pt)), crs = aea))
fp <- nsink_generate_flowpath(start_loc, niantic_nsink_data)
flow_summary <- nsink_summarize_flowpath(fp, removal)
flow_summary
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_calc_removal.R
\name{nsink_calc_removal_type}
\alias{nsink_calc_removal_type}
\title{Create removal type raster}
\usage{
nsink_calc_removal_type(removal_rasters)
}
\arguments{
\item{removal_rasters}{A named list of "land_removal", "stream_removal",
"lake_removal", and "raster_template" rasters plus a
sf object "huc".}
}
\value{
Raster of landscape nitrogen removal
}
\description{
Create removal type raster
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nsink_get_data.R
\name{nsink_get_data}
\alias{nsink_get_data}
\title{Gets N-Sink data for a given HUC}
\usage{
nsink_get_data(
  huc,
  data_dir = normalizePath("nsink_data", winslash = "/"),
  force = FALSE,
  year = "2016"
)
}
\arguments{
\item{huc}{A character string of a HUC identifier.  Currently can run only
a single HUC at a time.  This package was developed using 12-digit
HUCS, but has been (lightly) tested with larger HUCs and appears
to work, but is not certain for all cases.}

\item{data_dir}{A directory to store N-Sink data downloads.  Defaults to
"nsink_data" inside of the current working directory.
Created if it doesn't exist.  May be used for multiple HUCs
and only data that doesn't currently exist will be
downloaded.}

\item{force}{Logical to determine if files should be downloaded
again if they already exist locally.}

\item{year}{An argument to be passed to FedData's \code{\link{get_nlcd}}
function. Default is 2016.}
}
\value{
Returns a list with the huc used and the directory where the data is
        stored.
}
\description{
The required datasets for the N-sink analysis are available from multiple,
online resources.  This function takes a HUC as input and downloads local
copies of those datasets.  It will place these in a local directory and
inside that directory a new folder for each NHD Plus raster processing unit
will be created.  If you use the same folder for the data directory, it will
act as a cache (imperfect, though) and only download any new datasets that
have not already been downloaded.
}
\examples{
\dontrun{
library(nsink)
niantic_huc <- nsink_get_huc_id("Niantic River")$huc_12
nsink_get_data(huc = niantic_huc, data_dir = "nsink_data", force = TRUE)
}
}
