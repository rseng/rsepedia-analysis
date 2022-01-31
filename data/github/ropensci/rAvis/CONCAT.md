
[![Build Status](https://travis-ci.org/ropensci/rAvis.svg)](https://travis-ci.org/ropensci/rAvis)

[![Coverage Status](https://coveralls.io/repos/ropensci/rAvis/badge.svg)](https://coveralls.io/r/ropensci/rAvis)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/rAvis)](https://github.com/metacran/cranlogs.app)
[![cran version](http://www.r-pkg.org/badges/version/rAvis)](http://cran.rstudio.com/web/packages/rAvis)


rAvis
=====





`rAvis`: an R-package to download the information stored in “proyectoavis”, a citizen science bird project.

[Proyectoavis web site](http://proyectoavis.com/cgi-bin/portada.cgi)

## Installation

### Stable version from CRAN:


```r
install.packages("rAvis")
library("rAvis")
```


### Install with devtools

In the meantime you can install this development version with devtools package


```r
install.packages("devtools")
library("devtools")
install_github("ropensci/rAvis")
library("rAvis")
```


Load the library


```r
library("rAvis")
```



## Try some command

### Get species observation data


```r
Bubo <- avisQuerySpecies("Bubo bubo")
head(Bubo)
```

```
##   Id..Obs. Nombre.comun   Especie     Fecha Numero              Com.Aut
## 1    88274    Búho Real Bubo bubo 20-4-2014      2               Aragón
## 2    87902    Búho Real Bubo bubo 29-3-2014      1               Aragón
## 3    87516    Búho Real Bubo bubo 12-3-2014      1               Aragón
## 4    86733    Búho Real Bubo bubo  3-9-2013      1 Comunidad Valenciana
## 5    86732    Búho Real Bubo bubo  3-9-2013      1 Comunidad Valenciana
## 6    86731    Búho Real Bubo bubo  3-9-2013      1     Región de Murcia
##   Provincia     UTM Observador Periodo     Hora    Edad          Sexo
## 1    Teruel 30TXL41      torri      NA 08:05:00                      
## 2    Teruel 30TXL52      torri      NA 17:30:00                      
## 3    Teruel 30TXL41      torri      NA 17:30:00                      
## 4  Alicante 30SXH80  p.perales      NA 21:00:00  adulto         macho
## 5  Alicante 30SXH80  p.perales      NA 21:03:00 juvenil indeterminado
## 6    Murcia 30SXH80  p.perales      NA 21:00:00  adulto         macho
##                  Interes Grado.Reprod. Categ.Fenol.              Habitat
## 1                                    0            0                     
## 2                                    0            0                     
## 3    reproducción segura            11            0                     
## 4                                   16            3 roquedos de interior
## 5 comportamiento inusual             0            3   terrenos agrícolas
## 6                                   16            3 roquedos de interior
##   Codigo.de.Habitat
## 1                NA
## 2                NA
## 3                NA
## 4                NA
## 5                NA
## 6                NA
##                                                                                    Notas
## 1 barranco la Salobre - RUBIELOS de la CERIDA. DDAK 250m. S.O. 1 adulto y 1 pll al menos
## 2                                  Arroyo de los Calderones - TORRE los NEGROS. Con JMPS
## 3           Rambla del Ramblon  fuente la Salobre - RUBIELOS de la CERIDA. Con SMO y CAL
## 4                                                                                       
## 5                                      Posado sobre poste eléctrico a plena luz del día-
## 6                                                                                       
##         x     y
## 1 -1.2818 40.77
## 2 -1.1609 40.86
## 3 -1.2818 40.77
## 4 -0.8937 37.97
## 5 -0.8937 37.97
## 6 -0.8937 37.97
```


### Render a map


```r
avisMapSpecies("Pica pica")
```

![plot of chunk unnamed-chunk-6](inst/assets/figureunnamed-chunk-6.png) 


or with a physical map behind:


```r
avisMapSpecies("Pica pica", "phys")
```

![plot of chunk unnamed-chunk-7](inst/assets/figureunnamed-chunk-7.png) 


### Disable INFO messages


```r
avisSetup(verbose = FALSE)
```



## See help

Package help for more info

```r
?rAvis
```

## Meta

Please report any issues or bugs](https://github.com/ropensci/rAvis/issues).

License: GPL-2

This package is part of the [rOpenSci](http://ropensci.org/packages) project.

To cite package `rAvis` in publications use:

```coffee
To cite package ‘rAvis’ in publications use:

Javier González Hernández and Sara Varela (2014). rAvis: Interface to the bird-watching datasets at proyectoavis.com. R
package version 0.1.2.

A BibTeX entry for LaTeX users is

@Manual{,
title = {rAvis: Interface to the bird-watching datasets at proyectoavis.com},
author = {Javier González Hernández and Sara Varela},
year = {2014},
note = {R package version 0.1.2},
}

ATTENTION: This citation information has been auto-generated from the package DESCRIPTION file and may need manual editing, see
‘help("citation")’.
```

Get citation information for `rAvis` in R by `citation(package = 'rAvis')`

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
rAvis 0.1.4
===============

### BUG FIXES

* The layout of Proyecto Avis webpage changed, so we changed some functions of the package in order to being able to download the data again.  
rAvis
=====

```{r include=FALSE}
opts_chunk$set(fig.path = "inst/assets/figure")
```


`rAvis`: an R-package to download the information stored in “proyectoavis”, 
a citizen science bird project.

[Proyectoavis web site](http://proyectoavis.com/cgi-bin/portada.cgi)

## Installation

### Stable version from CRAN:

```{r eval=FALSE}
install.packages ("rAvis")
library("rAvis")
```

### Install with devtools

In the meantime you can install this development version with devtools package

```{r eval=FALSE}
install.packages("devtools")
library("devtools")
install_github("ropensci/rAvis")
library("rAvis")
```

Load the library

```{r}
library("rAvis")
```


## Try some command

### Get species observation data

```{r}
Bubo <- avisQuerySpecies("Bubo bubo")
head(Bubo)
```

### Render a map

```{r}
avisMapSpecies("Pica pica")
```

or with a physical map behind:

```{r}
avisMapSpecies("Pica pica", "phys")
```

### Disable INFO messages

```{r}
avisSetup(verbose=FALSE)
```


## See help

Package help for more info

```r
?rAvis
```

## Meta

Please report any issues or bugs](https://github.com/ropensci/rAvis/issues).

License: GPL-2

This package is part of the [rOpenSci](http://ropensci.org/packages) project.

To cite package `rAvis` in publications use:

```coffee
To cite package ‘rAvis’ in publications use:

Javier González Hernández and Sara Varela (2015). rAvis: Interface to the bird-watching datasets at proyectoavis.com. R
package version 0.1.3

A BibTeX entry for LaTeX users is

@Manual{,
title = {rAvis: Interface to the bird-watching datasets at proyectoavis.com},
author = {Javier González Hernández and Sara Varela},
year = {2015},
note = {R package version 0.1.3},
}

```

Get citation information for `rAvis` in R by `citation(package = 'rAvis')`

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/searchInterfaceFunctions.R
\name{avisQueryContributor}
\alias{avisQueryContributor}
\title{avisQueryContributor}
\usage{
avisQueryContributor(contributor_ids, args = list())
}
\arguments{
\item{contributor_ids}{must be either an integer or a list of contributors ids (integers)}

\item{args}{A list of normalized parameters to add filters to the query.
Currently in Spanish, but this might become outdated. See avisQuery.}
}
\value{
a dataframe with the results of your specific query to Proyecto AVIS database
}
\description{
Is a wrapper for avisQuery that allows you to perform
a query for more than one contributor at once.
}
\examples{
\dontrun{
avisQueryContributor(370)
avisQueryContributor(list(370, 399), args = list(year = 2002))
}
}
\seealso{
avisContributorsSummary
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/searchInterfaceFunctions.R
\name{avisQuery}
\alias{avisQuery}
\title{avisQuery}
\usage{
avisQuery(id_species = "", species = "", family = "", order = "",
age = "", sex = "", breeding = "", habitat = "", month = "", year = "",
args = list())
}
\arguments{
\item{id_species}{a number setting the id of the species according to proyectoavis.com database.
You may get the id of a species with \code{\link{avisSpeciesId}}}

\item{species}{scientific name of the species (one single species): e.g. "Passer domesticus"}

\item{family}{To filter the data by family: e.g. "Passeridae", "Falconidae", etc.}

\item{order}{To filter the data by Order: e.g. "Passeriformes", "Falconiformes", etc.}

\item{age}{To filter the data by age: "pollo", "juvenil", "adulto", "indeterminado".}

\item{sex}{To filter the data by sex: "macho", "hembra", "indeterminado", "pareja", "machos y hembras"}

\item{breeding}{To filter the data by breeding-migratory status: "reproducción posible", "reproducción probable",
"reproducción segura", "migración", "invernada"}

\item{habitat}{Filter by habitat: "bosque", "matorral", "pastizales", "terrenos agrícolas", "zonas humanizadas",
"zonas húmedas interiores", "roquedos de interior", "costas", "otros"}

\item{month}{Filter by month: 1 to 12}

\item{year}{Filter by year: e.g. 2001}

\item{args}{List of arguments accepted by www.proyectoavis.com endpoint. You may use
this list to set the arguments of the function (species, sex, breeding...), or
you may also set all the parameters supported by the endpoint, but not normalized for its use in this package.
These arguments are: id_ca, id_provincia, dia_ini, mes_ini, ano_ini, dia_fin, mes_fin, ano_fin, usu,
plazo, hora_ini, minuto_ini, hora_fin, minuto_fin, codigo_habitat, gr, cf, utm_10, utm_1 (see www.proyectoavis.com)}
}
\value{
a dataframe with the results of your specific query to Proyecto AVIS database.
}
\description{
General function for querying the database using several filters, like order,
family, species, age, sex, habitat, etc.
}
\details{
In case you set a query parameter by its name (eg: avisQuery (species="Bubo bubo"))
and also you set it inside the 'args' parameter (eg: avisQuery (species="Bubo bubo", args=list(species="Tyto alba")),
the value setted by its name will prevail (in the example, "Bubo bubo" will apply).
}
\examples{
\dontrun{
# get all the observations of the species of the Order Falconiformes
avisQuery (order = "Falconiformes")
# get all the observations of the species of the Family Falconidae
avisQuery(family = "Falconidae")
# get the observations of immatures of Iberian Imperial Eagle
avisQuery (species= "Aquila adalberti", age = "juvenil")
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/searchInterfaceFunctions.R
\name{avisQuerySpecies}
\alias{avisQuerySpecies}
\title{avisQuerySpecies}
\usage{
avisQuerySpecies(names, args = list())
}
\arguments{
\item{names}{Must be either a string or a list of scientific names}

\item{args}{A list of normalized parameters to add filters to the query.
Currently in Spanish, but this might become outdated. See avisQuery.}
}
\value{
a dataframe with the results of your specific query to Proyecto AVIS database
}
\description{
Is a wrapper for avisQuery that allows you to perform
a query for more than one species at once.'names' must
be either a string or a list of species names, 'args'
is a list of query parameters
(see avisQuery) that adds further filters to the query.
}
\examples{
\dontrun{
avisQuerySpecies("Bubo bubo")
avisQuerySpecies(list("Bubo bubo", "Tyto alba"), args = list(year = 2012))
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{peninsula}
\alias{peninsula}
\title{A physical map of the Iberian Peninsula}
\format{tif image}
\source{
http://www.openstreetmap.org/
}
\description{
A tif image downloaded from http://www.openstreetmap.org/ using the R- library OpenStreetMap
}
\keyword{datasets}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/rAvis_package.R
\docType{package}
\name{rAvis}
\alias{rAvis}
\alias{rAvis-package}
\title{rAvis: An R-package to download the information stored in Proyecto AVIS,
a citizen science bird project.}
\description{
We developed several functions to explore and donwload
the information stored in ProyectoAVIS database (www.proyectoavis.com),
in an easy and visual way.
}
\details{
We programmed two main functions to set flexible queries
about the species occurrences and the birdwatcher
observations: avisQuerySpecies and avisQueryContributor.
Besides, there are also general functions
to explore the database, like avisMapSpecies.

\tabular{ll}{
Package: \tab rAvis \cr
Type: \tab Package\cr
Version: \tab 0.1.3\cr
Date: \tab 2013-11-24\cr
License: \tab GPL-2 \cr
}
}
\examples{
\dontrun{
avisSpeciesSummary()

avisMapSpecies ("Pica pica", maptype="phys")

avisQuerySpecies(list("Bubo bubo", "Tyto alba"), args = list(year = 2012))
}
}
\author{
Javier Gonzalez Hernandez \email{javigzz@yahoo.es}

Sara Varela \email{svarela@paleobiogeography.org}
}
\references{
Varela S, Gonzalez-Hernandez J, Casabella E, Barrientos R (2014)
rAvis: An R-Package for Downloading Information Stored in Proyecto AVIS,
a Citizen Science Bird Project. PLoS ONE 9(3): e91650. doi: 10.1371/journal.pone.0091650
}
\seealso{
{
http://proyectoavis.com
}
}
\keyword{package}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/remoteUsersDataFunctions.R
\name{avisContributorsSummary}
\alias{avisContributorsSummary}
\title{avisContributorsSummary}
\usage{
avisContributorsSummary()
}
\value{
This function returns a matrix
}
\description{
Returns a table with the
observations aggregated by birdwatcher.
}
\note{
This function does not allow arguments
}
\examples{
\dontrun{
birdwatchers<- avisContributorsSummary()
par (mfrow =c(2,2))
plot (birdwatchers[,2],birdwatchers[,3], xlab=colnames (birdwatchers)[2],
ylab=colnames (birdwatchers)[3], pch=19)
plot (birdwatchers[,2],birdwatchers[,4], xlab=colnames (birdwatchers)[2],
ylab=colnames (birdwatchers)[4], pch=19)
plot (birdwatchers[,2],birdwatchers[,5], xlab=colnames (birdwatchers)[2],
ylab=colnames (birdwatchers)[5], pch=19)
plot (birdwatchers[,2],birdwatchers[,6], xlab=colnames (birdwatchers)[2],
ylab=colnames (birdwatchers)[6], pch=19)
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/remoteSpeciesDataFunctions.R
\name{avisSpeciesSummary}
\alias{avisSpeciesSummary}
\title{avisSpeciesSummary}
\usage{
avisSpeciesSummary()
}
\value{
returns a dataframe
}
\description{
Download a table with a summary of the records
stored in Proyecto AVIS (http://proyectoavis.com)
aggregated by species; number of observations of
each species, number of individuals recorded,
number of different UTMs (10x10km) with observations,
number of birdwatchers that recorded the species
}
\note{
This functions does not allow arguments
}
\examples{
\dontrun{
avis_summary<- avisSpeciesSummary()
# general overview of the data aggregated by species
par (mfrow =c(2,2))
hist (avis_summary$Observations, col="red", border=FALSE, main=NULL)
hist (avis_summary$Individuals, col="red", border=FALSE, main=NULL)
hist (avis_summary$UTM.10x10, col="red", border=FALSE, main=NULL)
hist (avis_summary$Birdwatchers, col="red", border=FALSE, main=NULL)
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/remoteSpeciesDataFunctions.R
\name{avisHasSpecies}
\alias{avisHasSpecies}
\title{avisHasSpecies}
\usage{
avisHasSpecies(nameraw)
}
\arguments{
\item{nameraw}{scientific name of the species (e.g. "Pica pica")}
}
\value{
Logical: returns TRUE for species with observations in the database and
FALSE otherwise
}
\description{
check if a species name exists in Proyecto AVIS.
}
\examples{
\dontrun{
avisHasSpecies("Pica pica")
avisHasSpecies("Pica pic")
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/remoteSpeciesDataFunctions.R
\name{avisSpeciesId}
\alias{avisSpeciesId}
\title{avisSpeciesId}
\usage{
avisSpeciesId(nameraw)
}
\arguments{
\item{nameraw}{scientific name of the species (e.g. "Pica pica")}
}
\value{
an integer
}
\description{
Returns the id of the selected species
}
\examples{
\dontrun{
avisSpeciesId("Pica pica")
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{canarias}
\alias{canarias}
\title{A physical map of the Canary Islands}
\format{tif image}
\source{
http://www.openstreetmap.org/
}
\description{
A tif image downloaded from http://www.openstreetmap.org/ using the R- library OpenStreetMap
}
\keyword{datasets}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{avisMap}
\alias{avisMap}
\title{Renders a map for the observations provided in 'obs'}
\usage{
avisMap(obs, label = "", maptype = "admin")
}
\arguments{
\item{obs}{set of observations returned by any of the avisQueryXXX functions}

\item{label}{label for the map. E.g. "Occurrences of Pica pica in Proyecto AVIS"}

\item{maptype}{Available types of map are 'admin',
administrative provinces of Spain (by default)
or 'phys, physical map of Spain.}
}
\value{
a plot with the occurrences of the species in the Iberian Peninsula. Maps have high resolution, so they could be printed.
}
\description{
This function should be used with avisQuerySpecies, to set a particular
query (with or without filters) and get the observations that we want to map.
It just allow to map one species. See avisMapSpecies for multiple maps.
}
\examples{
\dontrun{
obs<- avisQuerySpecies ("Pica pica", args = list(habitat = "bosque"))
avisMap(obs, label = "Pica pica")
avisMap(obs, label = "Pica pica", maptype = "phys")
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ravis_credentials}
\alias{ravis_credentials}
\title{ravis_credentials}
\format{string}
\description{
Credentials for API
}
\keyword{datasets}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/setupFunctions.R
\name{avisSetup}
\alias{avisSetup}
\title{avisSetup}
\usage{
avisSetup (...)
}
\arguments{
\item{...}{Package settings parameters. Available params: verbose = TRUE/FALSE}
}
\description{
Sets up settings that apply to the behabiour of the package
Allow users to turn off the information messages of the functions.
}
\examples{
\dontrun{
avisSetup(verbose=FALSE)
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ravis_shape_spain}
\alias{ravis_shape_spain}
\title{A Spanish administrative map}
\format{shapefile}
\source{
http://www.diva-gis.org/
}
\description{
A shapefile downloaded from  http://www.diva-gis.org/
}
\keyword{datasets}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/remoteUsersDataFunctions.R
\name{avisContributorAggregatedObservations}
\alias{avisContributorAggregatedObservations}
\title{avisContributorAggregatedObservations}
\usage{
avisContributorAggregatedObservations(contributor_id)
}
\arguments{
\item{contributor_id}{a number setting the id of the birdwatcher (see avisContributorSummary)}
}
\value{
This function returns a dataframe
}
\description{
A function to download the information about the
 observations of a birdwatcher.
}
\examples{
# Explore the contributions of Colectivo Ornitologico Ciguena Negra
\dontrun{
avisContributorAggregatedObservations (370)
}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{ravisUTMLatLong}
\alias{ravisUTMLatLong}
\title{UTM-Latlong}
\format{matrix}
\description{
Geographic coordinates (lat-long) of the centroids of the Spanish UTM squares
}
\keyword{datasets}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/plotFunctions.R
\name{avisMapSpecies}
\alias{avisMapSpecies}
\title{Renders a map for each of the species provided in names}
\usage{
avisMapSpecies(names, maptype = "admin", ...)
}
\arguments{
\item{names}{scientific name of the species
(it could be a list of scientific names). E.g. "Pica pica"}

\item{maptype}{Available types of map are 'admin',
administrative provinces of Spain (by default)
or 'phys', physical map of Spain.}

\item{...}{other filters passed to the observations query with avisQuerySpecies}
}
\value{
a plot with the occurrences of the species in the Iberian Peninsula. Maps have high resolution, so they could be printed.
}
\description{
This function map the species occurrences in the Iberian Peninsula.
}
\details{
For constructing these maps we used free online map repositories.
We downloaded the Spanish administrative map from  http://www.diva-gis.org/
and the Spanish physical map of http://www.openstreetmap.org/
using the R- library OpenStreetMap.
}
\examples{
\dontrun{

avisMapSpecies("Bubo bubo", "phys")

# if interested in several species, you can explore the database using avisMapSpecies
avisMapSpecies (list("Tyto alba", "Athene noctua", "Bubo bubo", "Strix aluco"),
               maptype="phys")

# and you can save those maps individually using the tiff function

directory<- "C:/your_directory"
species<- list("Tyto alba", "Athene noctua", "Bubo bubo", "Strix aluco")
for (x in species){
 tiff (file.path (directory, paste ("/", x, ".tiff", sep="")))
 avisMapSpecies (x)
 dev.off()
}

}
}

% Generated by roxygen2 (4.1.1): do not edit by hand
% Please edit documentation in R/remoteSpeciesDataFunctions.R
\name{avisAllSpecies}
\alias{avisAllSpecies}
\title{avisAllSpecies}
\usage{
avisAllSpecies()
}
\value{
returns a vector
}
\description{
Returns a vector with the ids of the species in Proyecto AVIS
}
\note{
This functions does not allow arguments
}
\examples{
\dontrun{
avisAllSpecies()
}
}

