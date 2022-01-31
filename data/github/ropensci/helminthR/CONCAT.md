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
helminthR
=======

[![R build status](https://github.com/ropensci/helminthR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/helminthR/actions)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/rmq9euldm5gy9qup?svg=true)](https://ci.appveyor.com/project/taddallas/helminthr)
[![codecov.io](https://codecov.io/github/ropensci/helminthR/coverage.svg?branch=master)](https://codecov.io/github/ropensci/helminthR?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/helminthR)](https://github.com/r-hub/cranlogs.app)


> Programmatically access the London Natural History Museum's [helminth database](https://www.nhm.ac.uk/research-curation/scientific-resources/taxonomy-systematics/host-parasites/index.html).

See software note in _Ecography_ ([available here](https://onlinelibrary.wiley.com/doi/full/10.1111/ecog.02131))


### Installation

From GitHub


```r
# install.packages("devtools")
devtools::install_github("rOpenSci/helminthR")
library("helminthR")
```

From CRAN


```r
install.packages("helminthR")
```



### Main functions

#### `findHost()`

Given a host genus and (optionally) species and location, this function returns all host-parasite associations of a given host species. The example below determines all parasite records for helminth infections of _Gorilla gorilla_.


```r
gorillaParasites <- findHost('Gorilla', 'gorilla')
head(gorillaParasites)
```

#### `findParasite()`

Given a helminth parasite genus (and optionally species, and location), this function returns a list of host-parasite records for that parasite. In the example below, I query the database for occurrences of the genus _Strongyloides_.


```r
strongHosts <- findParasite(genus='Strongyloides')
str(strongHosts)
```



#### `listLocations()` and `findLocation()`

List all location names (`listLocations()`). These names can be given to the `findLocation()` function, which finds all host-parasite associations that have occurred in the given location. Below, I look at host-parasite associations recorded in France.



```r
FrenchHostPars <- findLocation(location='France')
str(FrenchHostPars)
```




### Contribute!

Feel free to fork it and contribute some functionality.  



## Meta

* Please [report any issues or bugs](https://github.com/ropensci/helminthR/issues).
* License: GPL-3
* Get citation information for `helminthR` in R doing `citation(package = 'helminthR')`
* Please note that this project is released with a [Contributor Code of Conduct](https://www.contributor-covenant.org/).
By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
helminthR 1.0.9
==============

Fixing some NOTES and WARNINGS in the vignette creation process for the package. Switched from Travis CI to GitHub Actions for CI on GitHub. 



helminthR 1.0.8
==============

Small fix to the cleanData function to now use taxize to get host species taxonomic data.


helminthR 1.0.7
==============

Removed geocoding functionality previously present in listLocations(), as this now requires an API key. A cached version of the geographic coordinates of locations is provided as package data (`data(locations)`). 

Added extra catch in `cleanDat.R` to remove species who are identified as "something spp." instead of just removing those identified as "something sp.". 




helminthR 1.0.6
==============
* bug fix that was causing null results for some location specifications. 





helminthR 1.0.5
==============

* Released to CRAN.
## Test environments

* ubuntu 20.04-release
* ubuntu 20.04-devel
* windows-latest
* macOS-latest


## R CMD check results

R CMD check results
0 errors | 0 warnings | 0 notes


## Reverse dependencies

There are no reverse dependencies


---

I have read and agree to the the CRAN
policies at https://cran.r-project.org/web/packages/policies.html

This is an update to handle some NOTES and a few WARNINGS that were being thrown by certain operating systems when compiling the vignette. 

Thanks!
Tad Dallas
helminthR
=======

  [![R build status](https://github.com/ropensci/helminthR/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/helminthR/actions)
[![Windows Build Status](https://ci.appveyor.com/api/projects/status/rmq9euldm5gy9qup?svg=true)](https://ci.appveyor.com/project/taddallas/helminthr)
[![codecov.io](https://codecov.io/github/ropensci/helminthR/coverage.svg?branch=master)](https://codecov.io/github/ropensci/helminthR?branch=master)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/helminthR)](https://github.com/metacran/cranlogs.app)


> Programmatically access the London Natural History Museum's [helminth database](http://www.nhm.ac.uk/research-curation/scientific-resources/taxonomy-systematics/host-parasites/index.html).

See software note in _Ecography_ ([available here](http://onlinelibrary.wiley.com/doi/10.1111/ecog.02131/full))


### Installation

From GitHub

```{r eval=FALSE}

# install.packages("devtools")
devtools::install_github("rOpenSci/helminthR")
library("helminthR")

```

From CRAN

```{r eval=FALSE}

install.packages("helminthR")

```



### Main functions

#### `findHost()`

Given a host genus and (optionally) species and location, this function returns all host-parasite associations of a given host species. The example below determines all parasite records for helminth infections of _Gorilla gorilla_.

```{r eval=FALSE}

gorillaParasites <- findHost('Gorilla', 'gorilla')
head(gorillaParasites)

```

#### `findParasite()`

Given a helminth parasite genus (and optionally species, and location), this function returns a list of host-parasite records for that parasite. In the example below, I query the database for occurrences of the genus _Strongyloides_.

```{r eval=FALSE}

strongHosts <- findParasite(genus='Strongyloides')
str(strongHosts)

```



### `data(locations)` and `findLocation()`

A data file containing all the location names that can be queried, along with putative latitude and longitude coordinates for the centroid of each location can be found in `data(locations)`. Note that this will replace any object in the global environment named `locations`. These names can be given to the `findLocation()` function, which finds all host-parasite associations that have occurred in the given location. Below, I look at host-parasite associations recorded in France.


```{r eval=FALSE}

FrenchHostPars <- findLocation(location='France')
str(FrenchHostPars)

```




### Contribute!

Feel free to fork it and contribute some functionality.  



## Meta

* Please [report any issues or bugs](https://github.com/ropensci/helminthR/issues).
* License: GPL-3
* Get citation information for `helminthR` in R doing `citation(package = 'helminthR')`
* Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md).
By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "Introduction to the helminthR package"
author: "Tad Dallas"
date: ""
output: 
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteIndexEntry{Introduction to the helminthR package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


```{r echo=FALSE}
library(knitr)
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(if (abs(lines[1])>1) more else NULL,
            x[lines],
            if (length(x)>lines[abs(length(lines))]) more else NULL
           )
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
 })

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
library(helminthR)
```



This is an introduction to the `helminthR` package, whicih allows for the programmatic access to the London Natural History Museum's helminth parasite database ([available here](https://www.nhm.ac.uk/research-curation/scientific-resources/taxonomy-systematics/host-parasites/index.html)). This database represents the largest known database of host-helminth interactions, containing host-helminth occurrence records for over 350 locations, including aquatic, marine, and terrestrial locations. 



See software note in _Ecography_ ([available here](https://onlinelibrary.wiley.com/doi/10.1111/ecog.02131))



### Installation

From GitHub

```{r eval=FALSE}

# install.packages("devtools")
devtools::install_github("rOpenSci/helminthR")
library("helminthR")

```

From CRAN

```{r eval=FALSE}

install.packages("helminthR")

```



## Package functionality

The package allows for the acquisition of host-helminth interaction records based on host name (genus or species), parasite name (genus, species, or group), and/or location (accepted region name as provided in `data(locations)`. Parasite groups include "Acanthocephalans", "Cestodes", "Monogeans", "Nematodes", "Trematodes", or "Turbs" (short for Turbellarians). The user can further define host species as occuring 

1. "In the wild"
2. "Zoo captivity" 
3. "Domesticated"
4. "Experimental"
5. "Commercial source"
6. "Accidental infestation"


by inputting the corresponding number above in the `hostState` argument. 

The package itself has three main functions; `findHost`, `findParasite`, and `findLocation`. 



### Find all helminth parasites of a given host species

Given a host genus and (optionally) species and location, this function returns all host-parasite associations of a given host species. The example below determines all parasite records for helminth infections of _Gorilla gorilla_. We also use the `citation` argument here to obtain information on the citations which the host-helminth occurrences are based on. 

```{r eval=FALSE}

gorillaParasites <- findHost(genus='Gorilla', species='gorilla', 
	hostState=1, speciesOnly=TRUE, citation=TRUE)

```


The above function will query the database for helminth parasites of _Gorilla gorilla_ that were captured in the wild, and will remove helminth parasites not identified to species. If the user wishes to query multiple host species at the same time, the user can do the following 

```{r eval=FALSE}

hosts <- c('Gorilla gorilla', 'Peromyscus leucopus')
plyr::ldply(hosts, 
	function(x){
		findHost(unlist(strsplit(x, ' '))[1], 
			unlist(strsplit(x,' '))[2])})
	}

```



### Find all hosts of a given helminth parasite

Given a helminth parasite genus (and optionally species, and location), this function returns a list of host-parasite records for that parasite. In the example below, I query the database for occurrences of the genus _Strongyloides_.

```{r eval=TRUE}

strongHosts <- findParasite(genus='Strongyloides')
dim(strongHosts)

```

```{r}

head(strongHosts)

```



### `data(locations)` and `findLocation()`

A data file containing all the location names that can be queried, along with putative latitude and longitude coordinates for the centroid of each location can be found in `data(locations)`. Note that this will replace any object in the global environment named `locations`. These names can be given to the `findLocation()` function, which finds all host-parasite associations that have occurred in the given location. Below, I look at host-parasite associations recorded in France.


```{r eval=FALSE}

montanaOcc <- findLocation(location='Montana')
dim(montanaOcc)

```

```{r eval=FALSE}

head(montanaOcc)

```


### plotting host-helminth networks

Below, I provide an example of code for plotting the bipartite network of host-helminth interactions found in the state of Montana. 

```{r eval=FALSE}

g <- igraph::graph.incidence(table(montanaOcc[,1:2]))
igraph::V(g)$name <- NA
igraph::E(g)$color <- 'black'

plot(g, 
	vertex.color=c("black","dodgerblue")[igraph::V(g)$type+1],
	vertex.size=5
)
 
```








% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cleanData.R
\name{cleanData}
\alias{cleanData}
\title{Clean helminth parasite occurrence data}
\usage{
cleanData(edge, speciesOnly = FALSE, validateHosts = FALSE)
}
\arguments{
\item{edge}{Host-parasite edgelist obtained from \code{\link{findLocation}},
\code{\link{findHost}}, or \code{\link{findParasite}}}

\item{speciesOnly}{boolean flag to remove host and parasite species
where data are only available at genus level (default = FALSE)}

\item{validateHosts}{boolean flag to check host species names
against Catalogue of Life information and output taxonomic
information (default = FALSE)}
}
\value{
cleanEdge Host-parasite edgelist, but cleaned
}
\description{
Given a host-parasite edgelist, this function can validate species names,
provide further taxonomic information (thanks to \code{taxize}), 
and remove records only to genus level.
}
\details{
Use \code{data(locations)} for a list of possible locations.
}
\author{
Tad Dallas
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/locations.R
\docType{data}
\name{locations}
\alias{locations}
\title{Table of geographic location names, and associated coordinates}
\format{
\describe{
   \item{Location}{Name of geographic location}
   \item{Latitude}{Latitude of location centroid}
   \item{Longitude}{Longitude of location centroid}
 }
}
\usage{
data(locations)
}
\description{
Lists geographic locations that can be input to \code{\link{findHost}} or
\code{\link{findParasite}} and the corresponding latitude and longitude coordinates 
of the country's centroid. The georeferencing was performed dynamically using the 
Google Maps API, but they have since restricted access. The data on locations is now
provided in this data file called \code{locations} -- \code{data(locations)} -- and is based on 
an earlier usage of \code{ggmap}. The geographic coordinates may not be accurate, and users
should check for accuracy (and feel free to file an issue or PR on Github with corrections).
}
\references{
Gibson, D. I., Bray, R. A., & Harris, E. A. (Compilers) (2005).
Host-Parasite Database of the Natural History Museum, London.
}
\keyword{data}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findLocation.R
\name{findLocation}
\alias{findLocation}
\title{Find host-parasite interactions for a given location}
\usage{
findLocation(
  location = NULL,
  group = NULL,
  citation = FALSE,
  hostState = NULL,
  speciesOnly = FALSE,
  validateHosts = FALSE,
  removeDuplicates = FALSE
)
}
\arguments{
\item{location}{Location of host-parasite interaction.}

\item{group}{Parasite group - Cestodes, Acanthocephalans, Monogeneans, 
Nematodes, Trematodes, or Turbellarian etc. (Turb)}

\item{citation}{Boolean. Should the output include the citation link and 
the number of supporting citations? default is FALSE}

\item{hostState}{number corresponding to one of six different host states. 
The default value is NULL and includes all host states.}

\item{speciesOnly}{boolean flag to remove host and parasite species
where data are only available at genus level (default = FALSE)}

\item{validateHosts}{boolean flag to check host species names
against Catalogue of Life information and output taxonomic
information (default = FALSE)}

\item{removeDuplicates}{(boolean) should duplicate host-parasite 
combinations be removed? (default is FALSE)}
}
\value{
Three (or five) column data.frame containing host species, 
	parasite species (shortened name and full name), and citation link and 
	number of citations (if \code{citation = TRUE}), with each row corresponding 
	to an occurrence of a parasite species on a host species.
}
\description{
Given a location (available from \code{data{locations}}) this function 
returns all host-parasite associations in that location.
}
\details{
\code{hostState} can take values 1-6 corresponding to if the recorded 
host was found 
\itemize{ 
	\item (1) "In the wild"
	\item (2) "Zoo captivity" 
	\item (3) "Domesticated"
	\item (4) "Experimental"
	\item (5) "Commercial source"
	\item (6) "Accidental infestation"
 }
}
\examples{
\donttest{ FrenchHostPars <- helminthR::findLocation(location="France")}
}
\references{
Gibson, D. I., Bray, R. A., & Harris, E. A. (Compilers) (2005).
Host-Parasite Database of the Natural History Museum, London.
<http://www.nhm.ac.uk/research-curation/scientific-resources/taxonomy-systematics/host-parasites/>
}
\seealso{
\code{\link{findHost}}
}
\author{
Tad Dallas
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helminthR-package.R
\docType{package}
\name{helminthR-package}
\alias{helminthR-package}
\alias{helminthR}
\title{Access London Natural History Museum host-helminth record database}
\description{
'helminthR': A programmatic interface to the London 
  Natural History Museum's host-parasite database.

The package currently allows you to query by host species, parasite species, 
and geographic location. No information is provided on parasite prevalence or intensity.
}
\references{
Gibson, D. I., Bray, R. A., & Harris, E. A. (Compilers) (2005).
Host-Parasite Database of the Natural History Museum, London. 
<http://www.nhm.ac.uk/research-curation/scientific-resources/taxonomy-systematics/host-parasites/>
}
\author{
Tad Dallas \email{tad.a.dallas@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findHost.R
\name{findHost}
\alias{findHost}
\title{Find parasite occurrence data for given host.}
\usage{
findHost(
  genus = NULL,
  species = NULL,
  location = NULL,
  citation = FALSE,
  hostState = NULL,
  speciesOnly = FALSE,
  validateHosts = FALSE,
  parGroup = NULL,
  removeDuplicates = FALSE
)
}
\arguments{
\item{genus}{Host genus}

\item{species}{Host species}

\item{location}{Geographic location.}

\item{citation}{Boolean. Should the output include the citation link and 
the number of supporting citations? default is FALSE}

\item{hostState}{number corresponding to one of six different host states. 
The default value is NULL and includes all host states}

\item{speciesOnly}{boolean flag to remove host and parasite species
where data are only available at genus level (default = FALSE)}

\item{validateHosts}{boolean flag to check host species names
against Catalogue of Life information and output taxonomic
information (default = FALSE)}

\item{parGroup}{name of parasite group to query (default queries all groups)}

\item{removeDuplicates}{(boolean) should duplicate host-parasite 
combinations be removed? (default is FALSE)}
}
\value{
Three (or five) column data.frame containing host species, 
	parasite species (shortened name and full name), and citation link and 
	number of citations (if `citation`=TRUE), with each row corresponding 
	to an occurrence of a parasite species on a host species.
}
\description{
Given a host genus, species, and/or location, returns a list of parasite
occurrences on that host or for that location. 
Use \code{data(locations)} for a list of possible locations.
}
\details{
\code{hostState} can take values 1-6 corresponding to if the recorded host 
	was found 
\itemize{ 
	\item (1) "In the wild"
	\item (2) "Zoo captivity" 
	\item (3) "Domesticated"
	\item (4) "Experimental"
	\item (5) "Commercial source"
	\item (6) "Accidental infestation"
 }

 A value of NULL should be entered if you would like to include 
	all hostStates.

\code{parGroup} can be specified as "Acanthocephalans", "Cestodes",
	"Monogeans", "Nematodes", "Trematodes", or "Turbs" (Turbellarians etc.). 
	The default is to query all helminth parasite taxa.
}
\examples{

\donttest{gorillaParasites <- helminthR::findHost("Gorilla", "gorilla")}

# An example of how to query multiple hosts when you have a 
# vector of host species names

hosts <- c("Gorilla gorilla", "Peromyscus leucopus")
\donttest{plyr::ldply(hosts, function(x)
    {helminthR::findHost(unlist(strsplit(x, " "))[1], unlist(strsplit(x," "))[2])})}

}
\references{
Gibson, D. I., Bray, R. A., & Harris, E. A. (Compilers) (2005).
Host-Parasite Database of the Natural History Museum, London.
<http://www.nhm.ac.uk/research-curation/scientific-resources/taxonomy-systematics/host-parasites/>
}
\seealso{
\code{\link{findParasite}}
}
\author{
Tad Dallas
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/findParasite.R
\name{findParasite}
\alias{findParasite}
\title{Find host-parasite interactions for a given parasite species.}
\usage{
findParasite(
  genus = NULL,
  species = NULL,
  group = NULL,
  subgroup = NULL,
  location = NULL,
  citation = FALSE,
  hostState = NULL,
  speciesOnly = FALSE,
  validateHosts = FALSE,
  removeDuplicates = FALSE
)
}
\arguments{
\item{genus}{Parasite genus}

\item{species}{Parasite species}

\item{group}{Parasite group - Cestodes, Acanthocephalans, Monogeneans, 
Nematodes, Trematodes, or Turbellarian etc. (Turb)}

\item{subgroup}{Parasite subgroup (family names largely)}

\item{location}{Location of host-parasite interaction.}

\item{citation}{Boolean. Should the output include the citation link and 
the number of supporting citations? default is FALSE}

\item{hostState}{number corresponding to one of six different host states. 
The default value is NULL
     includes all host states}

\item{speciesOnly}{boolean flag to remove host and parasite species
where data are only available at genus level (default = FALSE)}

\item{validateHosts}{boolean flag to check host species names
against Catalogue of Life information and output taxonomic
information (default = FALSE)}

\item{removeDuplicates}{(boolean) should duplicate host-parasite 
combinations be removed? (default is FALSE)}
}
\value{
Three (or five) column data.frame containing host species, 
  parasite species (shortened name and full name), and citation link and 
  number of citations (if \code{citation = TRUE}), with each row corresponding 
  to an occurrence of a parasite species on a host species.
}
\description{
Given a host genus and/or species, this function returns a matrix containing
host-parasite interaction data. Search available locations using 
\code{data(locations)}.
}
\details{
\code{hostState} can take values 1-6 corresponding to if the recorded host 
was found 
\itemize{ 
	\item (1) "In the wild"
	\item (2) "Zoo captivity" 
	\item (3) "Domesticated"
	\item (4) "Experimental"
	\item (5) "Commercial source"
	\item (6) "Accidental infestation"
 }
}
\examples{

\donttest{strongHosts <- helminthR::findParasite(genus = "Strongyloides")}

# An example of how to query multiple parasite species when 
# you have a vector of parasite species names

parasites <- c("Ascaris aculeati", "Oxyuris flagellum")
\donttest{
 plyr::ldply(parasites, 
   function(x){
     helminthR::findParasite(unlist(strsplit(x, " "))[1], 
       unlist(strsplit(x," "))[2])
   }
 )
}
}
\references{
Gibson, D. I., Bray, R. A., & Harris, E. A. (Compilers) (2005).
Host-Parasite Database of the Natural History Museum, London. 
<http://www.nhm.ac.uk/research-curation/scientific-resources/taxonomy-systematics/host-parasites/>
}
\seealso{
\code{\link{findHost}}
}
\author{
Tad Dallas
}
