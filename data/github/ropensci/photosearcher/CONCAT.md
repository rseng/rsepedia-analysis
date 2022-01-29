
[![Travis build
status](https://travis-ci.org/ropensci/photosearcher.svg?branch=master)](https://travis-ci.org/ropensci/photosearcher)
[![Codecov test
coverage](https://codecov.io/gh/nfox29/photosearcher/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/photosearcher?branch=master)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- README.md is generated from README.Rmd. Please edit that file -->

# photosearcher <a href="https://docs.ropensci.org/photosearcher"><img src="man/figures/logo.png" align="right" height="132" /></a>

<p align="left">

<img src=".//man/figures/photosearcher_hex.png" height="200">

</p>

The goal of photosearcher is to provide a repeatable methodology for
accessing the Flickr API. More information can be found on the [package
website](https://docs.ropensci.org/photosearcher/). For more information
and examples of the functions check out the [package
vignette](https://docs.ropensci.org/photosearcher/articles/photosearcher.html).

We have produced this package to help facilitate reproducible code for
answering research questions. In the PLoS journals alone there are over
180 articles with Flickr in their keywords. [Click
here](https://docs.ropensci.org/photosearcher/articles/flickr_in_research.html)
for an overview of the use of Flickr in PLoS journals. Articles using
Flickr are published in a wide range of journals for a wide range of
research fields including [biology and life
sciences](https://www.nature.com/articles/s41598-017-18007-4), [computer
and information
sciences](https://www.inderscience.com/info/inarticle.php?artid=99808),
[medicine and health
sciences](https://www.sciencedirect.com/science/article/abs/pii/S0272494418303086)
and
[politics](https://journals.sagepub.com/doi/full/10.1177/1470357218780530?casa_token=UubfU8-MbuAAAAAA%3AAQBSE3ipGOdMi33J6ISalSySECPxvmxvmgDys3-ni7Z5EuQHNGlPMhOxjq6hyfPLmo1tFEIJYCiR).

## Terms of use

This package should be used within accordance of the [Flickr API terms
of use](https://www.flickr.com/help/terms/api).

## Installation

You can install the released version of photosearcher from github with:

``` r
devtools::install_github("nfox29/photosearcher")
```

## Getting an API key

The package requires a valid API key from the [Flickr development
page](https://www.flickr.com/services/apps/create/). The first time you
call a function from the package, you will be prompted to create and
enter your API key. The API key will then be saved as
photosearcher\_key.sysdata in your working directory and is used for all
function.

## Package functions

The package currently focuses on the ability to use the Flickr API to
search for images and their metadata through the `photo_search` function
([see the flickr.photos.search
method](https://www.flickr.com/services/api/flickr.photos.search.html))
for more information. These photographs can be downloaded using
`download_images` function which saves the images as a .jpeg file.

The package also allows users to find the top tags for a given location
(`location_tags`) and the tags most commonly associated with a given tag
(`related_tags`). The Flickr website offers full [API
Documentation](https://www.flickr.com/services/api/) for all of its call
methods.

## Searching for photo metadata

In this example, we demonstrate how to search for metadata on all images
labelled with the text or keywords *rock climbing* between 2010 and 2019
in the UK and Ireland.

``` r
library(photosearcher)

rock_climbing <- photo_search(
  mindate_taken = "2010-01-01",
  maxdate_taken = "2018-01-01",
  text = "rock climbing",
  bbox = "-12.875977,49.210420,2.636719,59.977005",
  has_geo = TRUE
)  
```

When `has_geo == TRUE` only metadata about images with latitude and
longitude information will be retrieved.

These can be plotted using other packages at user preference. In the
below example, we convert these to an `sf` object and plot using
`ggplot2`.

``` r
library(sf)
rock_climbing <- st_as_sf(rock_climbing, coords = c("longitude", "latitude"))

library(ggplot2)
ggplot() +
  geom_polygon(data = map_data("world", region = c("Ireland", "UK")), 
                               aes(x=long, y = lat, group = group),
                               fill = "lightgrey") + 
  geom_sf(data = rock_climbing) + 
  theme_bw()
```

## Issues and bugs

This package requires an internet connection as well as a connection to
the Flickr API, which may not be constantly available.

If you discover a bug not associated with connection to the API that is
not already a [reported
issue](https://github.com/ropensci/photosearcher/issues), please [open a
new issue](https://github.com/ropensci/photosearcher/issues/new)
providing a reproducible example.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
# How to contribute

Thank you for reading the contribution guidelines. The Flickr API has a whole range of methods currently not included here and we welcome any volunteer developers who wish to help increase the number of these included as functions in this package.
For a full list of the methods available through the Flickr API check out the [Flickr API website]( https://www.flickr.com/services/api/).

## Submitting changes

When sending a [GitHub Pull Request]( https://github.com/nfox29/photosearcher/pull/new/master
) stating clearly what you have changes you have made. For small changes single lines  messages are suffice, but for larger changes please provide a longer detailed message. 

## Coding conventions

We have a few simple coding conventions to follow:
* The package is designed following the [rOpenSci Guidelines]( http://devguide.ropensci.org/building.html). 
* Code style should be automated using the [styler package](https://github.com/r-lib/styler).
* As the package is open source, please ensure your code is clear to read for other. Where possible comment your code throughout.

Thank you,
Nathan
<!-- IF THIS INVOLVES AUTHENTICATION: DO NOT SHARE YOUR USERNAME/PASSWORD, OR API KEYS/TOKENS IN THIS ISSUE - MOST LIKELY THE MAINTAINER WILL HAVE THEIR OWN EQUIVALENT KEY -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
