<!-- README.md is generated from README.Rmd. Please edit that file -->
[![Travis-CI Build Status](https://travis-ci.org/ropensci/opencontext.png?branch=master)](https://travis-ci.org/ropensci/opencontext)

opencontext: An R API client for the Open Context archaeological data repository
--------------------------------------------------------------------------------

This packages enables browsing and downloading data from [Open Context](http://opencontext.org/) using R. Open Context reviews, edits, and publishes archaeological research data and archives data with university-backed repositories, including the California Digital Library.

Installation
------------

Install `opencontext`

``` r
install.packages("devtools")
devtools::install_github("ropensci/opencontext")
```

``` r
library("opencontext")
```

Browse countries
----------------

To browse the countries that Open Context has data on:

``` r
countries <- oc_browse("countries")
```

The result is a data frame that include the names of the countries in `countries$label`. URLs that we can use to get more information about what projects, etc. are available for each country in `countries$id`

Browse locations
----------------

To browse the locations for one country, for example, Turkey:

``` r
library("dplyr", warn.conflicts = FALSE)
locations <- oc_browse(type = "countries") %>%
   filter(label == "Turkey") %>%
   oc_get_countries(type = "location")
#> Getting data for Turkey
```

To browse the names of locations that have archaeological data in Turkey, run `locations$label`. We can see that the first location in this example is Çatalhöyük.

Browse projects
---------------

To inspect the projects available for a location in a country, for example, for Çatalhöyük in Turkey:

``` r
projects_at_Çatalhöyük_Turkey <- oc_get_locations("Turkey", "Çatalhöyük")
#> Getting data for Turkey
#> Getting data for Çatalhöyük
```

Once again, the `label` column has the names of the projects: `projects_at_Çatalhöyük_Turkey$label`.

With a little further effort we can browse excavation/survey areas within the project, and get datasets of measurements of objects collected from these areas (along with chronological and spatial data for these objects).

Get data from a specific project
--------------------------------

Now that we've identified a specific project, we can ingest data from that project into our R session.

------------------------------------------------------------------------

[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
---
output:
  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

[![Travis-CI Build Status](https://travis-ci.org/ropensci/opencontext.png?branch=master)](https://travis-ci.org/ropensci/opencontext)

## opencontext: An R API client for the Open Context archaeological data repository

This packages enables browsing and downloading data from [Open Context](http://opencontext.org/) using R. Open Context reviews, edits, and publishes archaeological research data and archives data with university-backed repositories, including the California Digital Library.

## Installation

Install `opencontext`

```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/opencontext")
```

```{r}
library("opencontext")
```

## Browse countries 

To browse the countries that Open Context has data on:

```{r}
countries <- oc_browse("countries")
```

The result is a data frame that include the names of the countries in `countries$label`. URLs that we can use to get more information about what projects, etc. are available for each country in `countries$id`

## Browse locations

To browse the locations for one country, for example, Turkey:

```{r}
library("dplyr", warn.conflicts = FALSE)
locations <- oc_browse(type = "countries") %>%
   filter(label == "Turkey") %>%
   oc_get_countries(type = "location")
```

To browse the names of locations that have archaeological data in Turkey, run `locations$label`. We can see that the first location in this example is Çatalhöyük. 

## Browse projects

To inspect the projects available for a location in a country, for example, for Çatalhöyük in Turkey:

```{r}
projects_at_Çatalhöyük_Turkey <- oc_get_locations("Turkey", "Çatalhöyük")
```

Once again, the `label` column has the names of the projects: `projects_at_Çatalhöyük_Turkey$label`. 

With a little further effort we can browse excavation/survey areas within the project, and get datasets of measurements of objects collected from these areas (along with chronological and spatial data for these objects).

## Get data from a specific project

Now that we've identified a specific project, we can ingest data from that project into our R session. 



---
[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/oc_get.R
\name{oc_get_countries}
\alias{oc_get_countries}
\title{Retrieve data in Open Context by country}
\usage{
oc_get_countries(country, type = c("projects", "locations", "descriptions"))
}
\arguments{
\item{country}{A character vectory of country names or a data frame
returned from \code{\link{oc_browse}}.}

\item{type}{The type of data to return: \code{"projects"},
\code{"locations"}, \code{"descriptions"}. The default is
\code{"projects"}.}
}
\value{
A data frame with the additional class \code{oc_dataframe}.
}
\description{
Given a character vector of one or more countries, or a data frame of
countries returned from \code{\link{oc_browse}}, this function retrieves data
related to those countries. The function can return projects, locations, or
descriptions.
}
\examples{
oc_get_countries("germany", type= "projects")

library(dplyr)
oc_browse(type = "countries") \%>\%
  filter(label == "Turkey") \%>\%
  oc_get_countries(type = "locations")
}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/oc_browse.R
\name{oc_browse}
\alias{oc_browse}
\title{Browse the Open Context archeological database}
\usage{
oc_browse(type = c("countries", "projects", "descriptions"),
  print_url = FALSE, ...)
}
\arguments{
\item{type}{The kind of to be returned. You can chose either
\code{'countries'} to get a data frame of names of countries that have Open
Context datasets, or \code{'projects'} to get a data frame project names,
or \code{'descriptions'} to get a data frame of data attributes that are
widely used in Open Context data sets.}

\item{print_url}{Whether or not to display a message with the URL of the
query. You can navigate to this URL to see the web interface's version of
the data returned by the API.}

\item{...}{Additional arguments passed to \code{\link[httr]{GET}}.}
}
\value{
A data frame with additional class \code{oc_dataframe}.
}
\description{
This function returns a data frame of certain types of top level data from
Open Context. You can get either a data frame of countries for which Open
Context has data, project names that have data on Open Context, or a list of
descriptions (Common Standards) of data attributes that are widely used in
Open Context datasets.
}
\examples{
oc_browse("countries")
oc_browse("projects")
}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/opencontext-package.r
\docType{package}
\name{opencontext}
\alias{opencontext}
\alias{opencontext-package}
\title{An API client for the Open Context archaeological database}
\description{
An API client for the Open Context archaeological database
}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/oc_get.R
\name{oc_get_locations}
\alias{oc_get_locations}
\title{Retrieve the names of projects in a given Country and Location}
\usage{
oc_get_locations(country, locations, type = c("projects", "descriptions"))
}
\arguments{
\item{country}{A country name}

\item{locations}{A character vector of locations in that country}

\item{type}{The type of records to return.}
}
\description{
Retrieve the names of projects in a given Country and Location
}
\examples{
oc_get_locations("Turkey", "Ulucak")
}

% Generated by roxygen2 (4.1.0): do not edit by hand
% Please edit documentation in R/oc_get.R
\name{oc_get_projects}
\alias{oc_get_projects}
\title{Retrieve data given an Open Context project name}
\usage{
oc_get_projects(project)
}
\arguments{
\item{project}{A character vector of project names}
}
\value{
A data frame of resources associated with the projects.
}
\description{
Retrieve data given an Open Context project name
}
\examples{
oc_get_projects("Kenan Tepe")
}

