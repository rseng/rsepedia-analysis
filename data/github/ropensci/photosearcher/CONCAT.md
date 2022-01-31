
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
---
output: github_document
---
[![Travis build status](https://travis-ci.org/ropensci/photosearcher.svg?branch=master)](https://travis-ci.org/ropensci/photosearcher) [![Codecov test coverage](https://codecov.io/gh/nfox29/photosearcher/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/photosearcher?branch=master) [![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
<!-- README.md is generated from README.Rmd. Please edit that file -->
```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# photosearcher <a href="https://docs.ropensci.org/photosearcher"><img src="man/figures/logo.png" align="right" height="132" /></a>

The goal of photosearcher is to provide a repeatable methodology for accessing the Flickr API. More information can be found on the [package website](https://docs.ropensci.org/photosearcher/). For more information and examples of the functions check out the [package vignette](https://docs.ropensci.org/photosearcher/articles/photosearcher.html).

We have produced this package to help facilitate reproducible code for answering research questions. In the PLoS journals alone there are over 180 articles with Flickr in their keywords. [Click here]( https://docs.ropensci.org/photosearcher/articles/flickr_in_research.html) for an overview of the use of Flickr in PLoS journals. Articles using Flickr are published in a wide range of journals for a wide range of research fields including [biology and life sciences](https://www.nature.com/articles/s41598-017-18007-4), [computer and information sciences](https://www.inderscience.com/info/inarticle.php?artid=99808), [medicine and health sciences](https://www.sciencedirect.com/science/article/abs/pii/S0272494418303086) and [politics](https://journals.sagepub.com/doi/full/10.1177/1470357218780530?casa_token=UubfU8-MbuAAAAAA%3AAQBSE3ipGOdMi33J6ISalSySECPxvmxvmgDys3-ni7Z5EuQHNGlPMhOxjq6hyfPLmo1tFEIJYCiR).

## Terms of use

This package should be used within accordance of the [Flickr API terms of use](https://www.flickr.com/help/terms/api).

## Installation

You can install the released version of photosearcher from github with:

``` r
devtools::install_github("nfox29/photosearcher")
```

## Getting an API key

The package requires a valid API key from the [Flickr development page](https://www.flickr.com/services/apps/create/). The first time you call a function from the package, you will be prompted to create and enter your API key. The API key will then be saved as photosearcher_key.sysdata in your working directory and is used for all function.

## Package functions

The package currently focuses on the ability to use the Flickr API to search for images and their metadata through the `photo_search` function ([ see the flickr.photos.search method](https://www.flickr.com/services/api/flickr.photos.search.html)) for more information. These photographs can be downloaded using `download_images` function which saves the images as a .jpeg file. 

The package also allows users to find the top tags for a given location (`location_tags`) and the tags most commonly associated with a given tag (`related_tags`). The Flickr website offers full [API Documentation](https://www.flickr.com/services/api/) for all of its call methods.

## Searching for photo metadata

In this example, we demonstrate how to search for metadata on all images labelled with the text or keywords *rock climbing* between 2010 and 2019 in the UK and Ireland. 

```{r eval=FALSE}
library(photosearcher)

rock_climbing <- photo_search(
  mindate_taken = "2010-01-01",
  maxdate_taken = "2018-01-01",
  text = "rock climbing",
  bbox = "-12.875977,49.210420,2.636719,59.977005",
  has_geo = TRUE
)  
```

When `has_geo == TRUE` only metadata about images with latitude and longitude information will be retrieved. 

These can be plotted using other packages at user preference. In the below example, we convert these to an `sf` object and plot using `ggplot2`. 

```{r eval=FALSE}
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

This package requires an internet connection as well as a connection to the Flickr API, which may not be constantly available. 

If you discover a bug not associated with connection to the API that is not already a [reported issue](https://github.com/ropensci/photosearcher/issues), please [open a new issue](https://github.com/ropensci/photosearcher/issues/new) providing a reproducible example.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "photosearcher"
author: "Nathan Fox"
date: "2019-29-04"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{photosearcher}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

This vignette provides an overview of the functions in R package `photosearcher`.

## Searching for photographs
The `photo_search()` function can be used to search Flickr for photographs that meet specific search criteria. For example, you could search for sightings of animals to better understand the species distribution:

```{r eval=FALSE}
#Search for photos of foxes in the UK for the year 2017
foxes <- photo_search(mindate_taken = "2017-01-01",
                      maxdate_taken = "2018-01-01",
                      text = "foxes",
                      bbox = "-12.875977,49.210420,2.636719,59.977005",
                      has_geo = TRUE)  

#Search for images of trees globally for the 1st of January 2019
trees <- photo_search(mindate_taken = "2019-01-01",
                      maxdate_taken = "2019-01-01",
                      text = "tree",
                      has_geo = TRUE)

#Search for images with text mountain and tagged lake
mountain_lake <- photo_search(mindate_taken = "2019-01-01",
                             maxdate_taken = "2019-02-01",
                             text = "mountain",
                             tags = "lake",
                             has_geo = FALSE)
```

Using the `text` parameter will search only for the phrase given. The `tag` parameter can be search for photographs that contain any of the given tags using `tags_any = TRUE`, or for only photographs that contain all of the given tags using `tags_any = FALSE`.

```{r eval=FALSE}
#Search for photographs with any of the tags
tags_any <- photo_search(mindate_taken = "2017-01-01",
                      maxdate_taken = "2018-01-01",
                      tags = c("lake", "mountain"),
                      tags_any = TRUE,
                      bbox = "-12.875977,49.210420,2.636719,59.977005",
                      has_geo = TRUE)

#Search for photographs with all of the tags
tags_all <- photo_search(mindate_taken = "2017-01-01",
                      maxdate_taken = "2018-01-01",
                      tags = c("lake", "mountain"),
                      tags_any = FALSE,
                      bbox = "-12.875977,49.210420,2.636719,59.977005",
                      has_geo = TRUE)

```


Users can also use the `sf` package to search via a simple features layer, this can be used to find images of in a specific location for example to aid in resource management for recreational activities in English national parks. When using an parameter that searches for a location or when `has_geo == TRUE` only metadata about images with latitude and longitude information will be retrieved. 

```{r eval=FALSE}
#load shape file
national_parks <- sf::st_read(system.file("shape/National_Parks_England.shp", package="photosearcher"))

#Search for photos of foxes in the UK for the year 2017
parks_hiking <- photo_search(mindate_taken = "2012-01-01",
                             maxdate_taken = "2012-04-01",
                             text = "hiking",
                             sf_layer = national_parks,
                             has_geo = TRUE) 
```

## Downloading images
The package can also download the images based of of their unique id number. You can direct the download to any save directory, if it does not exists it will be automatically created. Downloading images allows for validation of their contents. Here, we demonstrate how to download the first five images that were found using the foxes `photo_search` example, as well as for a single image. When downloading image, if users need certain pixel width and height dimensions (i.e. for Google Vision limits), these values can be specified.

```{r eval=FALSE}
#Downloading fox images from photosearch
downloaded_foxes <- download_images(photo_id = foxes$id[1:5], 
                                    save_dir = ".")

#Download a specific images based off its id
single_download <- download_images(photo_id = 48704764812, 
                                   save_dir = ".",
                                   max_image_height = 1200,
                                   max_image_width = 1200)
```

The package will only download images that the owner on Flickr has granted download permission. Those without permission will not be downloaded. The function also provides dataframe outlining which photographs were able to be downloaded. Permission to download does not automatically provide permission to use and distribute, check the photographs licence before use. The `photo_search` function provides the licence information for each image.

## Finding information on a user
The Flickr API can also return publicly available data on its users. The `user_info` function returns non-identifying information including hometown and occupation. This information could be used in studies aiming to calculate distance people travel to take photographs of certain places or things.  

```{r eval=FALSE}
#Find a users hometown
user <- user_info(user_id = "155421853@N05")

user$hometown

user$occupation
```

## Get the top tags for a location
Users interested in finding popular attractions or activities in \ given location can request the top tags for a given location. 

```{r eval=FALSE}
#Find the tags for Southampton, UK
southampton <- location_tags(woe_id = 35356)

#Find the tags for New York state, US
new_york <- location_tags(woe_id = 	2347591)
```

The woe_id argument is a Flickr specific place identifier. A places woeID can be found using the `find_place` function or using [online tools](http://woeid.rosselliot.co.nz/). 

## Find a locations woeID
Users can also call on the Flickr API to find information about a given location.
```{r eval=FALSE}
  #Find woeID for Kuala Lumpur
  kuala_lumpur <- find_place(place = "Kuala Lumpur")

  kuala_lumpaur$woeid
```

## Finding tags most associated with another tag
The Flickr API can perform clustered usage analysis on a given tag, returning the tags that are most often tagged alongside the given tag. For more information on how this is done check the [Flickrs getRelated documentation](https://www.flickr.com/services/api/flickr.tags.getRelated.html)

```{r eval=FALSE}
  #Find tags associated with zebra
  zebra <- related_tags(tag = "zebra")

  #Find tags most assoicated with hiking
  hiking <- related_tags(tag = "hiking")
  
  #Find tags most associated with the word church
  church <- related_tags(tag = "church")

```

## Finding Exif information about a photo
The Flickr API can also return Exchangeable image file format (Exif) data. For more information on Flickrs Exif data see the [Flickrs getExif documentation](https://www.flickr.com/services/api/flickr.photos.getExif.html) 
```{r eval=FALSE}
  #Find Exif information about a single photo
  exif_photo <- get_exif(photo_id = 48704764812)
```

## Finding information on a single photo
As well as searching for metadata on photographs based up on search terms, users can also find information on a single photograph using its unique photograph number.
```{r eval=FALSE}
  #Find information about a single photo
  single_photo <- get_photoinfo(photo_id = 48704764812)
```

## Find interesting photographs 

Flickr has an [algorithm](http://www.steves-digicams.com/knowledge-center/how-tos/online-sharing-social-networking/what-is-flickr-interestingness.html#b) for selecting a list of 500 interesting images for each day. Users can get this list using the `interesting_list` function.

```{r eval=FALSE}
  #Find information about a single photo
  interesting <- interesting_list(date = "2019-01-01")
```
---
title: "Flickr in research"
author: "Nathan Fox"
date: "2019-29-04"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Flickr in research}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE)

library(rplos)
library(tidyverse)
library(wordcloud)

get_kw <- function(x) {
  x %>% str_split(",") %>% 
    sapply(function(y) sub(".*/", "", y)) %>% 
    unique()  %>% 
    as.vector() %>% 
    paste(collapse = ",")
}
```

# Flickr papers

To understand how Flickr is currently being use in scientific research we investigated the discipline and scope of articles with the keyword `flickr` in the [PLoS journals](https://www.plos.org/).

```{r}
plos <- rplos::searchplos(q='flickr', 
                          fl=c('id','publication_date','subject', 'subject_level_1'),
                          limit = 200)$data %>% 
  
  mutate(keywords = map_chr(subject, get_kw),
         subjects = subject_level_1)
```

There are a total of `r nrow(plos)` papers with the keyword `flickr` in the PLoS journals. 

# Subject frequency

To understand which disciplines of science are using Flickr we calculated the frequency of papers in each research field:

```{r}
plos %>% 
  pull(subjects) %>% 
  paste(collapse = ",") %>% 
  str_split(",") %>% 
  table %>% 
  as_tibble()
```


# Keyword Wordcloud

To understand what topics are covered by Flickr researches we found the frequency of other keywords attributed to Flickr papers:

```{r}
kws <- plos %>% pull(keywords) %>% paste(collapse = ",") %>% str_split(",") %>% table %>% as.data.frame



wordcloud(words = kws$., freq = kws$Freq, min.freq = 1,
          max.words=200, random.order=FALSE, rot.per=0.35, 
          colors=brewer.pal(8, "Dark2"))

```





% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_photoinfo.R
\name{get_photoinfo}
\alias{get_photoinfo}
\title{Get metadata for a single photo}
\usage{
get_photoinfo(photo_id = NULL)
}
\arguments{
\item{photo_id}{Character, required. The id of the photo to get information
for.}
}
\value{
Dataframe of information on given image
}
\description{
Returns image metadata for a single photograph.
}
\details{
When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{
get_photoinfo(photo_id = 48704764812)
}
}
\seealso{
Other Get data for known photograph: \code{\link{download_images}},
  \code{\link{get_exif}}
}
\concept{Get data for known photograph}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/location_tags.R
\name{location_tags}
\alias{location_tags}
\title{Get top tags for a location}
\usage{
location_tags(woe_id)
}
\arguments{
\item{woe_id}{numeric. a "Where on Earth" location tag.}
}
\value{
character. List of the top 100 tags associated with the woe_id.
}
\description{
Takes user defined location and returns the top tags related to the location.
Uses the flickr.places.tagsForPlace API method from the Flickr API. See
\url{https://www.flickr.com/services/api/flickr.places.tagsForPlace.html} for
more information on the API method.
}
\details{
When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{
location_tags(woe_id = 35356)
}
}
\seealso{
Other Information about places: \code{\link{findPlaces}}
}
\concept{Information about places}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interesting_list.R
\name{interesting_list}
\alias{interesting_list}
\title{Interesting list}
\usage{
interesting_list(date = "2019-01-01")
}
\arguments{
\item{date}{Character, required. Interestingness list for the date provided.
The date should be in the form of "YYYY-MM-DD".}
}
\value{
data.frame. Output consists of 57 variables including; latitude and
longitude of photograph, date and time it was taken, associated tags and
image urls.

Full list of variables returned:

\itemize{
\item id: photograph's unique id number
\item owner: the unique id of the Flickr user
\item secret: photograph unique secret number
\item server: Flickr server data
\item farm: Flickr server data
\item title: photograph title
\item ispublic: whether photograph is public; 1 = yes, 0 = no
\item isfriend whether user is friend; 1 = yes, 0 = no
\item isfamily whether user is family; 1 = yes, 0 = no
\item license: use licence of the image see \url{https://www.flickr.com/services/api/flickr.photos.licenses.getInfo.html} for details
\item datetaken: date and time of image capture
\item datetakengranularity: accuracy of image date see \url{https://www.flickr.com/services/api/misc.dates.html} for more information on dates
\item datetakenunknown: whether date is unknown see \url{https://www.flickr.com/services/api/misc.dates.html} for more information on dates
\item count_views: number of view the photograph has had,
\item count_comments: number of comments on the photograph
\item count_faves: number of times the photograph has been favourited
\item tags: user defined tags on the photograph
\item latitude: latitude of where the image was taken
\item longitude: longitude of where the image was taken
\item accuracy: accuracy of spatial reference see \url{https://www.flickr.com/services/api/flickr.photos.search.html } for more information
\item context: a numeric value representing the photo's geotagginess beyond latitude and longitude \url{https://www.flickr.com/services/api/flickr.photos.search.html } for more information
\item place_id: unique numeric number representing the location of the photograph
\item woeid: unique numeric number representing the location of the photograph
\item geo_is_family: whether only friends can see geo; 1 = yes, 0 = no
\item geo_is_friend: whether only family can see geo; 1 = yes, 0 = no
\item geo_is_contact: whether only contact can see geo; 1 = yes, 0 = no
\item geo_is_public whether geo is public; 1 = yes, 0 = no
\item url_sq: URL for square image
\item height_sq: height for square image
\item width_sq: width for square image
\item url_t: URL for square image thumbnail image 100 on longest side
\item height_t: height for thumbnail image 100 on longest side,
\item width_t: width for thumbnail image 100 on longest side
\item url_s: URL for small square image 75x75
\item height_s: height for small square image 75x75
\item width_s: width for small square image 75x75
\item url_q: URL for large square image 150x150
\item height_q: height for large square image 150x150
\item width_q: width for large square image 150x150
\item url_m: URL for small image 240 on longest side
\item height_m: height for small image 240 on longest side
\item width_m: width for small image 240 on longest side
\item url_n: URL for small image 320 on longest side
\item height_n: height for small image 320 on longest side
\item width_n: width for small image 320 on longest side
\item url_z: URL for medium image 640 on longest side
\item height_z: height for medium image 640 on longest side
\item width_z: width for medium image 640 on longest side
\item url_c: URL for medium image 800 on longest side
\item height_c: height for medium image 800 on longest side
\item width_c: width for medium image 800 on longest side
\item url_l: URL for large image 1024 on longest side
\item height_l: height for large image 1024 on longest side
\item width_l: width for large image 1024 on longest side
\item url_o: URL for original image, either a jpg, gif or png, depending on source format
\item height_o: height for original image, either a jpg, gif or png, depending on source format
\item width_o: width for original image, either a jpg, gif or png, depending on source format
}
}
\description{
Returns a Flickr generated list of photographs deemed interesting. For
information on how this list is calculated see:
\url{http://www.steves-digicams.com/knowledge-center/how-tos/online-sharing-social-networking/what-is-flickr-interestingness.html#b}
}
\details{
When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{

interesting_list(date = "2019-01-01")

interesting_list(date = "2017-05-05")

interesting_list(date = "2011-11-25")

}
}
\seealso{
Other Search for images: \code{\link{photo_search}}
}
\concept{Search for images}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/photosearcher-package.R
\docType{package}
\name{photosearcher-package}
\alias{photosearcher}
\alias{photosearcher-package}
\title{photosearcher: Photo Searcher}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Queries the Flick API (https://www.flickr.com/services/api/) to return photograph metadata as well as the ability to download the images as jpegs.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/photosearcher}
  \item Report bugs at \url{https://github.com/ropensci/photosearcher/issues}
}

}
\author{
\strong{Maintainer}: Nathan Fox \email{nf2g13@soton.ac.uk}

Authors:
\itemize{
  \item Francesca Mancini
  \item Laura Graham
  \item Louis Sutter
  \item Tom August
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_exif.R
\name{get_exif}
\alias{get_exif}
\title{Get Exif data}
\usage{
get_exif(photo_id = NULL)
}
\arguments{
\item{photo_id}{Id of photograph}
}
\value{
A dataframe of "exchangeable image file format" information for the
given photograph
}
\description{
Returns Exchangeable image file format data for a single photograph. For more
information on how Exif data differs from metadata see:
\url{https://www.flickr.com/services/api/flickr.photos.getExif.html}
}
\details{
When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{
get_exif(photo_id = 48704764812)
}
}
\seealso{
Other Get data for known photograph: \code{\link{download_images}},
  \code{\link{get_photoinfo}}
}
\concept{Get data for known photograph}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_info.R
\name{user_info}
\alias{user_info}
\title{Get user information}
\usage{
user_info(user_id)
}
\arguments{
\item{user_id}{character. The id of the user you wish to obtain information
for.}
}
\value{
data.frame. Dataframe of 5 variables from the searched users publicly
available information: id, occupation, hometown, city, country.
}
\description{
Takes Flickr user ID and returns the profile data.
}
\details{
Uses the flickr.profile.getProfile API method from the Flickr API. See
\url{https://www.flickr.com/services/api/flickr.profile.getProfile.html} for
more information on the API method.

See \url{https://www.pixsy.com/academy/flickr-id/} for a guide on finding
your Flickr ID.

When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{
user_info(user_id = "155421853@N05")
}
}
\concept{User information}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_place.R
\name{findPlaces}
\alias{findPlaces}
\alias{find_place}
\title{Find place location data}
\usage{
find_place(place)
}
\arguments{
\item{place}{character. The place for the query}
}
\value{
data.frame. Information on locations that share the name with the
search location. Nine variables are returned: place_id, woe_id, latitude,
longitude, place_url, place_type, place_type_id, timezone, woe_name.
}
\description{
Find the WoeID and other location data for a given place
}
\details{
Takes user defined location and returns location data for the search. Uses
the flickr.places.find API method from the Flickr API. See
\url{https://www.flickr.com/services/api/flickr.places.find.html} for more
information on the API method.

When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{
find_place(place = "New York")

find_place(place = "England")
}
}
\seealso{
Other Information about places: \code{\link{location_tags}}
}
\concept{Information about places}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_images.R
\name{download_images}
\alias{download_images}
\title{Download images}
\usage{
download_images(photo_id, save_dir = NULL, max_image_height = NULL,
  max_image_width = NULL, overwrite_file = FALSE, quiet = FALSE)
}
\arguments{
\item{photo_id}{numeric or character vector. id of photo to download, can be
single id, list or column for photo_search outputs}

\item{save_dir}{character. name of directory for photos to be saved in.}

\item{max_image_height}{numeric. maximum number of pixels for images height}

\item{max_image_width}{numeric. maximum number of pixels for images width}

\item{overwrite_file}{logical. Whether to overwritten existing files. if
TRUE, files will be overwritten and you will be warned in the output.
Default is FALSE.}

\item{quiet}{logical. If TRUE, suppress status messages (if any), and the
progress bar.}
}
\value{
character. A vector of the images attempted to be downloaded and
whether they were. If an image was not downloaded, information on why is
provided. Images will be saved to \code{save_dir}.
}
\description{
Downloads images based on their Flickr ID. Uses the flickr.photos.getSizes
API method from the Flickr API to test whether you have permission to
download an image. See
\url{https://www.flickr.com/services/api/flickr.photos.getSizes.html} for
more information on the API method. If permission is available the image is
downloaded and saved as a .jpeg in a given save directory.
}
\details{
Please be aware that download times will vary depending on number of
photographs, size of photographs, internet speed and other factors.
Downloading a large amount of photographs may take some time.

When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{
download_images(photo_id = 48704764812, save_dir = ".")

download_images(photo_id = 48704764812, max_image_height = 1200,
max_image_width = 1200, save_dir = ".") }
}
\seealso{
Other Get data for known photograph: \code{\link{get_exif}},
  \code{\link{get_photoinfo}}
}
\concept{Get data for known photograph}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/related_tags.R
\name{related_tags}
\alias{related_tags}
\title{Get related tags}
\usage{
related_tags(tag)
}
\arguments{
\item{tag}{character. tag to search.}
}
\value{
character. Tags most associated with input tag.
}
\description{
Takes a tag and returns the top tags related to that tag.
}
\details{
Uses the flickr.tags.getRelated API method from the Flickr API. See
\url{https://www.flickr.com/services/api/flickr.tags.getRelated.html} for
more information on the API method.

When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{
related_tags(tag = "car")

related_tags(tag = "monkey")

related_tags(tag = "river")
}
}
\concept{Tag information}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/photo_search.R
\name{photo_search}
\alias{photo_search}
\title{Search for photo metadata}
\usage{
photo_search(mindate_taken = NULL, maxdate_taken = NULL,
  mindate_uploaded = NULL, maxdate_uploaded = NULL, user_id = NULL,
  text = NULL, tags = NULL, tags_any = TRUE, bbox = NULL,
  woe_id = NULL, sf_layer = NULL, has_geo = TRUE)
}
\arguments{
\item{mindate_taken}{Character, or date required. Minimum taken date. Photos
with an taken date greater than or equal to this value will be returned. The
date should be in the form of "YYYY-MM-DD".}

\item{maxdate_taken}{Character, or date required. Maximum taken date. Photos
with an taken date less than or equal to this value will be returned. The
date should be in the form of "YYYY-MM-DD".}

\item{mindate_uploaded}{Character or date, optional. Minimum upload date.
Photos with an upload date greater than or equal to this value will be
returned. The date can be in the form of a unix timestamp or mysql datetime.}

\item{maxdate_uploaded}{Character or date, optional. Maximum upload date.
Photos with an upload date less than or equal to this value will be
returned. The date can be in the form of a unix timestamp or mysql datetime.}

\item{user_id}{Character, optional. The Flickr ID of the user who's photo to
search. If this parameter isn't passed then everybody's public photos will
be searched.}

\item{text}{Character, optional. A free text search. Photos who's title,
description or tags contain the text will be returned. You can exclude
results that match a term by prepending it with a - character. Free text
searches for words in order provided, for example a search for "climbing
rock" will be different to "rock climbing".}

\item{tags}{Character vector, optional. A comma-delimited list of tags. Photos
with one or more of the tags listed will be returned. You can exclude
results that match a term by prepending it with a - character.}

\item{tags_any}{Logical, optional. If TRUE, photos containing any of the tags
will be returned. If FALSE, only photos containing all tags will be
returned. Defaulted to return any tags.}

\item{bbox}{String, optional bounding box of search area provide as:
"minimum_longitude,minimum_latitude,maximum_longitude,maximum_latitude".}

\item{woe_id}{String, optional "where on earth identifier" can be supplied
instead of bbox. Use \code{\link{find_place}} to obtain woe_id for a place.}

\item{sf_layer}{A simple features layer, optional, area to search for photos,
can be supplied instead of a bbox or woeID.}

\item{has_geo}{Logical, optional parameter for whether returned photos need
associated spatial data.}
}
\value{
data.frame. Output consists of 57 variables including; latitude and
longitude of photograph, date and time it was taken, associated tags and
image urls.

Full list of variables returned:

\itemize{ \item id: photograph's unique id number \item owner: the unique id
of the Flickr user \item secret: photograph unique secret number \item
server: Flickr server data \item farm: Flickr server data \item title:
photograph title \item ispublic: whether photograph is public; 1 = yes, 0 =
no \item isfriend whether user is friend; 1 = yes, 0 = no \item isfamily
whether user is family; 1 = yes, 0 = no \item license: use licence of the
image see
\url{https://www.flickr.com/services/api/flickr.photos.licenses.getInfo.html}
for details \item datetaken: date and time of image capture \item
datetakengranularity: accuracy of image date see
\url{https://www.flickr.com/services/api/misc.dates.html} for more
information on dates \item datetakenunknown: whether date is unknown see
\url{https://www.flickr.com/services/api/misc.dates.html} for more
information on dates \item count_views: number of view the photograph has
had, \item count_comments: number of comments on the photograph \item
count_faves: number of times the photograph has been favourited \item tags:
user defined tags on the photograph \item latitude: latitude of where the
image was taken \item longitude: longitude of where the image was taken
\item accuracy: accuracy of spatial reference see
\url{https://www.flickr.com/services/api/flickr.photos.search.html } for
more information \item context: a numeric value representing the photo's
geotagginess beyond latitude and longitude
\url{https://www.flickr.com/services/api/flickr.photos.search.html } for
more information \item place_id: unique numeric number representing the
location of the photograph \item woeid: unique numeric number representing
the location of the photograph \item geo_is_family: whether only friends can
see geo; 1 = yes, 0 = no \item geo_is_friend: whether only family can see
geo; 1 = yes, 0 = no \item geo_is_contact: whether only contact can see geo;
1 = yes, 0 = no \item geo_is_public whether geo is public; 1 = yes, 0 = no
\item url_sq: URL for square image \item height_sq: height for square image
\item width_sq: width for square image \item url_t: URL for square image
thumbnail image 100 on longest side \item height_t: height for thumbnail
image 100 on longest side, \item width_t: width for thumbnail image 100 on
longest side \item url_s: URL for small square image 75x75 \item height_s:
height for small square image 75x75 \item width_s: width for small square
image 75x75 \item url_q: URL for large square image 150x150 \item height_q:
height for large square image 150x150 \item width_q: width for large square
image 150x150 \item url_m: URL for small image 240 on longest side \item
height_m: height for small image 240 on longest side \item width_m: width
for small image 240 on longest side \item url_n: URL for small image 320 on
longest side \item height_n: height for small image 320 on longest side
\item width_n: width for small image 320 on longest side \item url_z: URL
for medium image 640 on longest side \item height_z: height for medium image
640 on longest side \item width_z: width for medium image 640 on longest
side \item url_c: URL for medium image 800 on longest side \item height_c:
height for medium image 800 on longest side \item width_c: width for medium
image 800 on longest side \item url_l: URL for large image 1024 on longest
side \item height_l: height for large image 1024 on longest side \item
width_l: width for large image 1024 on longest side \item url_o: URL for
original image, either a jpg, gif or png, depending on source format \item
height_o: height for original image, either a jpg, gif or png, depending on
source format \item width_o: width for original image, either a jpg, gif or
png, depending on source format }
}
\description{
Returns image metadata for photos matching the search terms.
}
\details{
Uses the flickr.photos.search API method from the Flickr API. This search
method requires a limiting factor to prevent parameterless searches - to
ensure this is met the function requires both a minimum and a maximum date
that searched photographs were taken on. See
\url{https://www.flickr.com/services/api/flickr.photos.search.html} for more
information on the API method.

When running the function you need an API key saved as
photosearcher_key.sysdata in your working directory. If this is the first
function you run you will be prompted to create and enter your API key. The
API key will then be saved as photosearcher_key.sysdata in your working
directory and is used for all function.
}
\examples{
\dontrun{
photo_search(
  mindate_taken = "2019-01-01",
  maxdate_taken = "2019-01-02",
  text = "tree",
  bbox = "-13.623047,47.279229,3.251953,60.630102",
  has_geo = TRUE
  )

photo_search(
  mindate_taken = "2017-01-01",
  maxdate_taken = "2017-01-02",
  mindate_uploaded = "2017-03-04",
  maxdate_uploaded = "2017-05-05",
  tags = "lake"
  )

photo_search(
  mindate_taken = "2018-01-01",
  maxdate_taken = "2018-01-03",
  tags = c("mountain", "lake"),
  tags_any = TRUE
)

photo_search(
  mindate_taken = "2018-01-01",
  maxdate_taken = "2018-01-03",
  tags = c("mountain", "lake"),
  tags_any = FALSE
)

}
}
\seealso{
Other Search for images: \code{\link{interesting_list}}
}
\concept{Search for images}
