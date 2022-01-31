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
# AntWeb
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)
![CRAN/GitHub 0.7_/_0.7.4.99](https://img.shields.io/badge/CRAN/GitHub-0.7_/_0.7.4.99-blue.svg) 

[__AntWeb__](http://www.antweb.org/) is a repository of ant specimen records maintained by the [California Academy of Sciences](http://www.calacademy.org/). From the website's description:  

> AntWeb is the world's largest online database of images, specimen records, and natural history information on ants. It is community driven and open to contribution from anyone with specimen records, natural history comments, or images.

__Resources__  
* [AntWeb](http://www.antweb.org/)   
* [AntWeb API](http://www.antweb.org/api/)
* [API version 2](http://www.antweb.org/api/v2/)

## Package Status and Installation

[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/antweb?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/antweb)
[![Travis-CI Build Status](https://travis-ci.org/ropensci/antweb.svg?branch=master)](https://travis-ci.org/)
 [![codecov](https://codecov.io/gh/RMHogervorst/antweb/branch/master/graph/badge.svg)](https://codecov.io/gh/RMHogervorst/antweb)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/AntWeb?color=blue)](https://github.com/metacran/cranlogs.app)

__Installation Instructions__

__Stable version__  

```coffee
install.packages("AntWeb", dependencies = TRUE)
# version 0.6
```

__Development version__  

```coffee
# If you don't already have the devtools package installed, run
# install.packages("devtools")
# unlike most packages, devtools requires additional non-R dependencies depending on your OS. 
# See → https://github.com/ropensci/rOpenSci/wiki/Installing-devtools
library(devtools)
install_github("ropensci/AntWeb", ref = "dev")
```

## Usage

| Function name | Description | Example | 
| ------------- | ----------- | ------- |
| `aw_data`  | Search for data by taxonomic level, full species name, a bounding box, habitat, elevation or type.     |    __Search by a species name__ <br> `aw_data(scientific_name = "acanthognathus brevicornis")` <br> __or by a genus__ <br> `crem <- aw_data(genus = "crematogaster")`  <br> __Search by a bounding box__ <br> `aw_data(bbox = '37.77,-122.46,37.76,-122.47')` <br> __Search by an elevation band__ <br> `aw_data(min_elevation = 1500, max_elevation = 2000)` |
| `aw_unique` | Obtain a list of unique levels by various taxonomic ranks    | `aw_unique(rank = "subfamily")` <br>`genus_list <- aw_unique(rank = "genus")`<br>`aw_unique(rank = "species")` |
| `aw_images` | Search photos by type or time since added.     |    ` aw_images(since = 5)`<br> `aw_images(since = 5, type = "h")` |
| `aw_coords` | Search for specimens by location and radius     |    `aw_coords(coord = "37.76,-122.45", r = 5)` |
| `aw_code` | Search for a specimen by record number   |  `aw_code(occurrenceid = "CAS:ANTWEB:alas188691")` |
| `aw_map` | Map georeferenced data | `adf <- aw_data(genus = "acanthognathus", georeferenced = TRUE)`<br>`aw_map(adf)` |

## Citation

```r
To cite package ‘AntWeb’ in publications use:

  'Karthik Ram' (2014). AntWeb: programmatic interface
  to the AntWeb. R package version 0.7.2.99.
  https://github.com/ropensci/AntWeb

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {AntWeb: programmatic interface to the AntWeb},
    author = {'Karthik Ram'},
    year = {2014},
    note = {R package version 0.7.2.99},
    url = {https://github.com/ropensci/AntWeb},
  }
```

---
  
This package is part of a richer suite called [SPOCC Species Occurrence Data](https://github.com/ropensci/spocc), along with several other packages, that provide access to occurrence records from multiple databases. We recommend using SPOCC as the primary R interface to AntWeb unless your needs are limited to this single source.    

---

### Questions, bugs, and suggestions

Please file any bugs or questions as [issues](https://github.com/ropensci/AntWeb/issues/new) or send in a pull request.

---

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md).
By participating in this project you agree to abide by its terms.


[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)

 
<!--
%\VignetteEngine{knitr::knitr}
%\VignetteIndexEntry{An Introduction to the EML package}
-->

# Guide to using the AntWeb R interface

![](http://www.antweb.org/images/casent0104669/casent0104669_d_1_high.jpg)

AntWeb is the world's largest online database of images, specimen records, and natural history information on ants. The data repository is maintained by the [California Academy of Sciences](http://www.calacademy.org/). It is community driven and open to contribution from anyone with specimen records, natural history comments, or images. This package provides a programmatic interface to the data. This package is part of the rOpenSci suite of tools.

```{r, echo = FALSE, message = FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  error = FALSE,
  tidy = TRUE,
  fig.width=12, 
  fig.height=8, 
  fig.path='Figs/'
  )
```

## Obtaining specimen data  

### Searching the database  

1. Searching by Genus

```{r, genus}
leaf_cutter_ants  <- aw_data(genus = "acromyrmex")
head(leaf_cutter_ants)
``

### Retrieving records by specimen id  


## Searching records by location  

## Searching images  






% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_data.R
\name{aw_data_all}
\alias{aw_data_all}
\title{Download all aw_data available for any request}
\usage{
aw_data_all(..., progress = "text")
}
\arguments{
\item{...}{All the same arguments that get passed to \code{aw_data}}

\item{progress}{Default is on and set to \code{text}. Set to \code{none} to suppress}
}
\description{
This is a thin wrapper around aw_data
}
\examples{
\dontrun{
# crem <- aw_data_all(genus = "crematogaster", georeferenced = TRUE)
}
}
\seealso{
aw_data
}
\keyword{data}
\keyword{download}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_data.R
\name{aw_data}
\alias{aw_data}
\title{Retrieve data from the AntWeb}
\usage{
aw_data(genus = NULL, species = NULL, scientific_name = NULL,
  georeferenced = NULL, min_elevation = NULL, max_elevation = NULL,
  type = NULL, habitat = NULL, country = NULL, min_date = NULL,
  max_date = NULL, bbox = NULL, limit = NULL, offset = NULL,
  quiet = FALSE)
}
\arguments{
\item{genus}{An ant genus name}

\item{species}{a species name}

\item{scientific_name}{An easier way to pass the Genus and species name together, especially when the data are derived from other packages.}

\item{georeferenced}{Default is \code{FALSE}. Set to \code{TRUE} to return only data with lat/long information. Note that this filtering takes place on the client-side, not server side.}

\item{min_elevation}{A lower elevation bound}

\item{max_elevation}{An upper elevation bound}

\item{type}{A holotype}

\item{habitat}{A fuzzy search by any habitat}

\item{country}{A country name}

\item{min_date}{A lower date bound in the format \code{yyyy-mm-dd}}

\item{max_date}{An upper date bound in the format \code{yyyy-mm-dd}}

\item{bbox}{A lat long bounding box. Format is \code{lat,long,lat,long}. Use this website: http://boundingbox.klokantech.com/ to quickly grab a bbox (set format on bottom left to csv and be sure to switch the order from long, lat, long, lat to lat, long, lat, long)
Just set the format on the bottom left to CSV.}

\item{limit}{A numeric value to limit number of records}

\item{offset}{An offset best used with limit as a way to paginate records}

\item{quiet}{If true, any informative messages will be suppressed}
}
\value{
data.frame
}
\description{
This function allows a user to query the AntWeb database by any taxonomic rank or full species name.
}
\examples{
  
# data <- aw_data(genus = "acanthognathus", species = "brevicornis")
# data3 <- aw_data(genus = "acanthognathus", species = "brevicornis", georeferenced = TRUE)
# data2 <- aw_data(scientific_name = "acanthognathus brevicornis")
# sandstone <- aw_data(genus = "Aphaenogaster", habitat = "sandstone")
# data_genus_only <- aw_data(genus = "acanthognathus", limit = 25)
# leaf_cutter_ants  <- aw_data(genus = "acromyrmex")
# data  <- aw_data(genus = "Technomyrmex", bbox = '37.77,-122.46,37.76,-122.47')
# Search just using a bounding box
# data  <- aw_data(bbox = '37.77,-122.46,37.76,-122.47')
# Search by a elevation band
# aw_data(min_elevation = 1500, max_elevation = 2000)
# When you throw a really specimen rich band like below, you'll get a huge number of requests. 
# Only the first 1000 records will download first. 
# aw_data(min_elevation = 200, max_elevation = 400)
# aw_data(min_date = '1980-01-01', max_date = '1981-01-01')
# fail <- aw_data(scientific_name = "auberti levithorax") # This should fail gracefully
}
\keyword{data}
\keyword{download}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_map.R
\name{aw_map}
\alias{aw_map}
\title{LeafletJS Map}
\usage{
aw_map(aw_obj, dest = tempdir(), title = "AntWeb species map",
  incl.data = TRUE)
}
\arguments{
\item{aw_obj}{Result from a search on AntWeb}

\item{dest}{Location where the html file and geojson file should be stored. Default is the temp directory}

\item{title}{Title of the map.}

\item{incl.data}{Default is \code{TRUE}. Writes geoJSON data into the html file to get around security restrictions in browsers like Google Chrome. Set to \code{FALSE} to read from a separate local geoJSON file.}
}
\description{
Builds an interactive map of locations for any list of species
}
\examples{
\dontrun{
 acanthognathus_df <- aw_data(genus = "acanthognathus", georeferenced = TRUE)
 aw_map(acanthognathus_df)
# Or just plot data by habitat. So for e.g. using sandstone as a substrate
sandstone <- aw_data(habitat = "sandstone")
aw_map(sandstone)
}
}
\keyword{map}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_package.R
\docType{package}
\name{AntWeb}
\alias{AntWeb}
\alias{AntWeb-package}
\title{AntWeb}
\description{
AntWeb
}
\details{
\href{http://www.antweb.org/}{The AntWeb} world's largest online database of images, specimen records, and natural history information on ants. The database is maintained and hosted by the \href{http://www.calacademy.org/}{California Academy of Sciences}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_images.R
\name{aw_images}
\alias{aw_images}
\title{aw_images}
\usage{
aw_images(since = NULL, img_type = NULL)
}
\arguments{
\item{since}{number of days in the past to query}

\item{img_type}{h for head, d for dorsal, p for profile, and l for label. If a img_type is not specified, all images are retrieved.}
}
\value{
data.frame
}
\description{
Download ant images based on time elapsed and/or type.
}
\examples{
\dontrun{
z <- aw_images(since = 5)
z1 <- aw_images(since = 5, img_type = "d")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_coords.R
\name{aw_coords}
\alias{aw_coords}
\title{aw_coords}
\usage{
aw_coords(coord = NULL, r = NULL)
}
\arguments{
\item{coord}{Latitude and Longitude. Should be supplied as \code{lat,long}. Example: \code{37.76,-122.45}}

\item{r}{A radius in kilometers. For 2 km add \code{r = 2}}
}
\value{
\code{\link{aw_data}}
}
\description{
Retrieve AntWeb data by location. A radius argument can be supplied as a search radius around a point on th emap.
}
\examples{
 
# data_by_loc <- aw_coords(coord = "37.76,-122.45", r = 2)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_code.R
\name{aw_code}
\alias{aw_code}
\title{aw_code}
\usage{
aw_code(occurrenceid = NULL, catalogNumber = NULL)
}
\arguments{
\item{occurrenceid}{A unique id in the AntWeb database identifying a particular specimen}

\item{catalogNumber}{Specimen catalogue number}
}
\value{
list
}
\description{
Retrieve data by specimen id
}
\examples{
# data_by_code <- aw_code(occurrenceid = "CAS:ANTWEB:alas188691") 
# data_by_code <- aw_code(catalognumber="inb0003695883")
}
\seealso{
\code{\link{aw_data}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_utils.R
\name{aw_cbind}
\alias{aw_cbind}
\title{aw_cbind}
\usage{
aw_cbind(results)
}
\arguments{
\item{results}{A list of objects of class \code{antweb}}
}
\description{
Allows for combining split AntWeb calls (e.g. paginated calls) back into one single result object
}
\examples{
\dontrun{
x1 <- aw_data(genus = "crematogaster", georeferenced = TRUE)
x2 <- aw_data(genus = "crematogaster", georeferenced = TRUE, offset = 1000)
x12 <- aw_cbind(list(x1, x2))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_data.R
\name{aw_unique}
\alias{aw_unique}
\title{aw_unique}
\usage{
aw_unique(rank = NULL, name = NULL)
}
\arguments{
\item{rank}{A taxonomic rank. Allowed values are  \code{subfamily}, \code{genus} or \code{species}}

\item{name}{Optional. If left blank, the query will return a list of all unique names inside the supplied rank.}
}
\value{
data.frame
}
\description{
Get a list of unique names within any taxonomic rank
}
\examples{
 \dontrun{
subfamily_list <- aw_unique(rank = "subfamily")
# genus_list <- aw_unique(rank = "genus")
# species_list <- aw_unique(rank = "species")
}
}
\seealso{
\code{\link{aw_data}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_utils.R
\name{print.antweb}
\alias{print.antweb}
\title{Print a summary for an antweb object}
\usage{
\method{print}{antweb}(x, ...)
}
\arguments{
\item{x}{An object of class \code{antweb}}

\item{...}{additional arguments}
}
\description{
Print a summary for an antweb object
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/aw_distinct.R
\name{aw_distinct}
\alias{aw_distinct}
\title{aw_distinct}
\usage{
aw_distinct(rank = "genus", habitat = NULL, country = NULL,
  min_elevation = NULL, max_elevation = NULL, limit = 1000,
  offset = NULL)
}
\arguments{
\item{rank}{= "genus" Default is genus. But you can also use phylum, sub-phylum etc}

\item{habitat}{The habitat type}

\item{country}{Country name}

\item{min_elevation}{Min elevation recorded for specimen}

\item{max_elevation}{Max elevation recorded for specimen}

\item{limit}{= 1000 Default limit. Set higher if necessary}

\item{offset}{To be used in conjunction with limit}
}
\description{
Retrieves a data.frame of distinct ranks based on various restrictions
}
\examples{
\dontrun{
aw_distinct(rank = "genus", country = "Madagascar")
}
}
