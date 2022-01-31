[![Build Status](https://travis-ci.org/ropensci/paleobioDB.svg?branch=master)](https://travis-ci.org/ropensci/paleobioDB)
[![codecov.io](https://codecov.io/github/ropensci/paleobioDB/coverage.svg?branch=master)](https://codecov.io/github/ropensci/paleobioDB?branch=master)
[![rstudio mirror downloads](http://cranlogs.r-pkg.org/badges/paleobioDB)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/paleobioDB)](https://CRAN.R-project.org/package=paleobioDB)


paleobioDB
=======

### About

`paleobioDB` is a package for downloading, visualizing and processing data from [Paleobiology Database](http://paleobiodb.org/).


### Quick start

**Install**

Install paleobioDB from CRAN

```coffee
install.packages("paleobioDB")
library(paleobioDB)
```

Install paleobioDB developing version from github

```coffee
install.packages("devtools")
library(devtools)
install_github("ropensci/paleobioDB")
library(paleobioDB)
```

**General overview**

`paleobioDB` has 19 functions to wrap each endpoint of the PaleobioDB API, plus 8 functions to visualize and process the fossil data. The API documentation for the Paleobiology Database can be found [here](http://paleobiodb.org/data1.1/).

## Download fossil occurrences from the PaleobioDB

**pbdb_occurrences** 

e.g., to download all the fossil data that belongs to the family Canidae, set base_name = "Canidae".  

```coffee
> canidae<-  pbdb_occurrences (limit="all",
                             base_name="canidae", vocab="pbdb",
                             interval="Quaternary",             
                             show=c("coords", "phylo", "ident"))
head(canidae)
```

```coffee
## occurrence_no record_type collection_no            taxon_name taxon_rank taxon_no
## 150070  occurrence         13293              Cuon sp.      genus    41204
## 186572  occurrence         18320         Canis cf. sp.      genus    41198
## 186573  occurrence         18320            Vulpes sp.      genus    41248
## 186574  occurrence         18320        Borophagus sp.      genus    41196
## 192926  occurrence         19617        Canis edwardii    species    44838
## 192927  occurrence         19617 Canis armbrusteri cf.    species    44827
## matched_rank     early_interval    late_interval early_age late_age reference_no
## 5 Middle Pleistocene Late Pleistocene     0.81   0.0117         4412
## 5   Late Hemphillian          Blancan    10.300   1.8000         6086
## 5   Late Hemphillian          Blancan    10.300   1.8000         6086
## 5   Late Hemphillian          Blancan    10.300   1.8000         6086
## 3            Blancan     Irvingtonian     4.900   0.3000         2673
## 3            Blancan     Irvingtonian     4.900   0.3000         2673
## lng      lat  family family_no     order order_no    class class_no   phylum
## 111.56667 22.76667 Canidae     41189 Carnivora    36905 Mammalia    36651 Chordata
## -85.79195 40.45444 Canidae     41189 Carnivora    36905 Mammalia    36651 Chordata
## -85.79195 40.45444 Canidae     41189 Carnivora    36905 Mammalia    36651 Chordata
## -85.79195 40.45444 Canidae     41189 Carnivora    36905 Mammalia    36651 Chordata
## -112.40000 35.70000 Canidae     41189 Carnivora    36905 Mammalia    36651 Chordata
## -112.40000 35.70000 Canidae     41189 Carnivora    36905 Mammalia    36651 Chordata
## phylum_no genus_name species_name genus_reso reid_no species_reso matched_name
## 33815       Cuon          sp.       <NA>      NA         <NA>         <NA>
## 33815      Canis          sp.        cf.      NA         <NA>         <NA>
## 33815     Vulpes          sp.       <NA>      NA         <NA>         <NA>
## 33815 Borophagus          sp.       <NA>      NA         <NA>         <NA>
## 33815      Canis     edwardii       <NA>    8376         <NA>         <NA>
## 33815      Canis  armbrusteri       <NA>    8377          cf.         <NA>
## subgenus_name subgenus_reso
## <NA>          <NA>
## <NA>          <NA>
## <NA>          <NA>
## <NA>          <NA>
## <NA>          <NA>
## <NA>          <NA>
```

**CAUTION WITH THE RAW DATA**

Beware of synonyms and errors, they could twist your estimations about species richness, evolutionary and extinction rates, etc. paleobioDB users should be critical about the raw data downloaded from the database and filter the data before analyzing it.

For instance, when using "base_name" for downloading the information with the function pbdb_occurrences, check out the synonyms and errors that could appear in "taxon_name", "genus_name", etc. In our example, in canidae$genus_name there are errors: "Canidae" and "Caninae" appeared as genus names. If not eliminated, they will increase the richness of genera. 


## Map the fossil records

**pbdb_map**

Returns a map with the species occurrences.

```coffee
> pbdb_map(canidae)
``` 
![plot of chunk map](man/figure/pbdb_map.png)


**pbdb_map_occur**
Returns a map and a raster object with the sampling effort (number of fossil records per cell).

```coffee
> pbdb_map_occur (canidae, res= 5)
``` 
```coffee
## class       : RasterLayer 
## dimensions  : 34, 74, 2516  (nrow, ncol, ncell)
## resolution  : 5, 5  (x, y)
## extent      : -179.9572, 190.0428, -86.42609, 83.57391  (xmin, ## xmax, ymin, ymax)
## coord. ref. : NA 
## data source : in memory
## names       : layer 
## values      : 1, 40  (min, max)
``` 

![plot of chunk map](man/figure/pbdb_map_occur.png)


**pbdb_map_richness**
Returns a map and a raster object with the number of different species, genera, family, etc. per cell. The user can change the resolution of the cells. 

```coffee
> pbdb_map_richness (canidae, res= 5, rank="species")
```
```coffee
## class       : RasterLayer 
## dimensions  : 34, 74, 2516  (nrow, ncol, ncell)
## resolution  : 5, 5  (x, y)
## extent      : -179.9572, 190.0428, -86.42609, 83.57391  (xmin, xmax, ymin, ymax)
## coord. ref. : NA 
## data source : in memory
## names       : layer 
## values      : 1, 12  (min, max)
```

![plot of chunk map](man/figure/pbdb_map_occur.png)


## Explore your fossil data 


**pbdb_temporal_range**

Returns a dataframe and a plot with the time span of the species, genera, families, etc. in your query.

```coffee
> pbdb_temp_range (canidae, rank="species")
``` 
```coffee
                           max    min
## Canis brevirostris        5.3330 0.0000
## Canis mesomelas           5.3330 0.
## Alopex praeglacialis      5.3330 0.0117
## Nyctereutes megamastoides 5.3330 0.0117
## Vulpes atlantica          5.3330 0.0117
## Canis latrans             4.9000 0.0000
...

``` 
![plot temprange](man/figure/pbdb_temporal_range.png)



**pbdb_richness**

Returns a dataframe and a plot with the number of species (or genera, families, etc.) across time. You should set the temporal extent and the temporal resolution for the steps.

```coffee
> pbdb_richness (canidae, rank="species", temporal_extent=c(0,10), res=1)
```

```coffee
## labels2 richness
## <=1       23
## 1-2       56
## 2-3       53
## 3-4       19
## 4-5       18
## 5-6        5
## 6-7        0
## 7-8        0
## 8-9        0
## 9-10       0
## >10        0
``` 
![plot richness](man/figure/pbdb_richness.png)


**pbdb_orig_ext**

Returns a dataframe and a plot with the number of new appearances and last appearances of species, genera, families, etc. in your query across the time. You should set the temporal extent and the resolution of the steps. 

```coffee
# evolutionary rates= orig_ext=1
> pbdb_orig_ext (canidae, rank="species", orig_ext=1, temporal_extent=c(0,10), res=1)
```
```coffee
##              new ext
## 1-2 to 0-1    0  28
## 2-3 to 1-2   34   6
## 3-4 to 2-3    1   0
## 4-5 to 3-4   13   0
## 5-6 to 4-5    5   0
## 6-7 to 5-6    0   0
## 7-8 to 6-7    0   0
## 8-9 to 7-8    0   0
## 9-10 to 8-9   0   0
```

![plot of chunk map](man/figure/pbdb_orig_ext_1.png)


```coffee
# extinction rates= orig_ext=2
pbdb_orig_ext(canidae, rank="species", orig_ext=2, temporal_extent=c(0,10), res=1)
``` 
```coffee
##             new ext
## 1-2 to 0-1    0  28
## 2-3 to 1-2   34   6
## 3-4 to 2-3    1   0
## 4-5 to 3-4   13   0
## 5-6 to 4-5    5   0
## 6-7 to 5-6    0   0
## 7-8 to 6-7    0   0
## 8-9 to 7-8    0   0
## 9-10 to 8-9   0   0
``` 

![plot of chunk map](man/figure/pbdb_orig_ext_2.png)

**pbdb_subtaxa**

Returns a plot and a dataframe with the number of species, genera, families, etc. in your dataset.
  
```coffee
> pbdb_subtaxa (canidae, do.plot=TRUE)         
```
```coffee
## species genera families orders classes phyla
## 75     24        1      1       1     1
```
![plot subtaxa](man/figure/pbdb_subtaxa.png)


**pbdb_temporal_resolution**

Returns a plot and a dataframe with a main summary of the temporal resolution of the fossil records

```coffee
> pbdb_temporal_resolution (canidae)
```   

```coffee
## $summary
## Min. 1st Qu.  Median    Mean 3rd Qu.    Max. 
## 0.0117  0.1143  1.5000  1.5360  2.5760 23.0200 
## 
## $temporal_resolution
## [1]  0.7693  8.5000  8.5000  8.5000  4.6000
## [6]  4.6000  4.6000  3.1000  3.1000  3.1000
## [11]  3.1000  4.6000  3.1000  3.1000  3.1000
## [16]  3.1000  3.1000  3.1000  3.1000  3.1000
## [21]  3.1000  3.1000  3.1000  3.1000  3.1000
## [26]  3.1000  3.1000  3.1000  3.1000  3.1000
...
  
```
![plot tempres](man/figure/pbdb_temporal_resolution.png)

## Docker

We are including a Dockerfile to ease working on the package as it fulfills all the system dependencies of the
package.

How to load the package with Docker:

1. Install Docker. Reference here: https://docs.docker.com/get-started
2. Build the *docker image*: from the root folder of this repository. Type:
```coffee
docker build -t rpbdb Docker
```
This command will create a *docker image* in your system based on some of the [rocker/tidyverse](https://hub.docker.com/r/rocker/tidyverse/) images.
You can see the new image with ```docker image ls```.
3. Start a container for this image. Type the following command picking some *<password>* of your choice.
```coffee
docker run -d --rm -p 8787:8787 -e PASSWORD=<password> -v $PWD:/home/rstudio rpbdb
```
This will start a container with access to your current folder where all the code of the package is.
Inside the container, the code will be located in */home/rstudio*. It also exposes the port 8787 of the container so you may access
to the RStudio web application which is bundled in the *rocker* base image.
4. Navigate to http://localhost:8787. Enter with user=*rstudio* and the password you used in the command above.
5. You may enter to the container via console with:
```coffeef
docker exec -ti ravis bash
```
Either from RStudio or from within the container you can install the package out of the code with:
```coffee
cd /home/rstudio
R
library(devtools)
install.packages(".", repos=NULL, type = "source")
```

## Meta

Please report any [issues or bugs](https://github.com/ropensci/pbdb/issues).

License: GPL-2

To cite package `paleobioDB` in publications use:

```coffee
To cite package `paleobioDB` in publications use:

Sara Varela, Javier Gonzalez-Hernandez and Luciano Fabris Sgarbi (2016). paleobioDB: an R-package for downloading, visualizing and processing data from the Paleobiology Database. R package version 0.5. https://github.com/ropensci/paleobioDB

A BibTeX entry for LaTeX users is

  @Manual{,
    title = {paleobioDB: an R-package for downloading, visualizing and processing data from the Paleobiology Database},
    author = {{Sara Varela} and {Javier Gonzalez-Hernandez} and {Luciano Fabris Sgarbi}},
    year = {2014},
    note = {R package version 0.7},
    base = {https://github.com/ropensci/paleobioDB},
  }
```

---

This package is part of the [rOpenSci](http://ropensci.org/packages) project.

[![](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
paleobioDB 0.7.0
================

MINOR IMPROVEMENTS

* Ensures travis fails on either API or internal tests failure

BUG FIXES

* Fix [issue 45](https://github.com/ropensci/paleobioDB/issues/45) about wrong API endpoint.
* Fix errors in several functions in `rpbdb_temporal_functions.R`


paleobioDB 0.6.0
================

MINOR IMPROVEMENTS

* Integrates codecov
* Adds tests execution to travis. Tests are not part of the CRAN check for this
package since they rely on calling to the PaleobioDB API.
* Improves error reporting in pbdb_temp_range. [issue 28](https://github.com/ropensci/paleobioDB/issues/28)

BUG FIXES

* Fix error in pbdb_occurrences under R 3.5.


paleobioDB 0.5.0
===============

BUG FIXES

* Fix bug after JSON responses from the API started including the "elapsed_time" before the "records" element.


paleobioDB 0.4.0
===============

MINOR IMPROVEMENTS

* From now our functions for exploring the data of the paleobioDB use "matched_name" and "matched_number" instead of raw names and raw taxon numbers of the former versions of the package. This means that the default output of richness and first and last occurrences of species/genus/etc. changed. Before we used the original taxonomic identification of the fossil records, now we use the revised taxonomic identification of the fossil records.

BUG FIXES

* Fix duplicated records bug in pbdb_occurrences caused by API returning and array of two identical values in field "reference_no". From now on, if a filed is returned by the API as an array it will be mapped to a semicolon separated string.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_strata_auto}
\alias{pbdb_strata_auto}
\title{pbdb_strata_auto}
\usage{
pbdb_strata_auto (...)
}
\arguments{
\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/strata/auto}. Eg:
\itemize{
  \item \emph{name}: a full or partial name. You can
    use \% and _ as wildcards, but the
    query will be very slow if you put a wildcard at the beginning
  \item \emph{rank}: returns only strata of the specified 
    rank: formation, group or member.
  \item \emph{lngmin}: numeric. The longitude boundaries will be normalized to fall
    between -180 and 180. Note that if you specify lngmin then you must also
    specify lngmax. Returns only records whose geographic location falls within
    the given bounding box (defined by lngmin, lngmax, latmin, latmax). It
    generate two adjacent bounding boxes if the range crosses the antimeridian. 
  \item \emph{lngmax}: numeric. The longitude boundaries will be normalized to fall
    between -180 and 180.
  \item \emph{latmin}: numeric. between -90 and 90. Note that if you specify latmin
    then you must also specify latmax.
  \item \emph{latmax}: numeric. between -90 and 90.
  \item \emph{loc}: Return only strata associated with some occurrence whose
    geographic location falls within the specified geometry, specified in WKT
    format.
  \item \emph{vocab}: set vocab="pbdb" to show the complete name of the variables (by
    default variables have short 3-letter names)
  \item ...
}}
}
\value{
a dataframe with information from the strata that matches our letters.
}
\description{
Returns a list of strata matching the given prefix or partial name. 
This can be used to implement auto-completion for strata names, 
and can be limited by geographic location if desired.
}
\examples{
\dontrun{
pbdb_strata_auto (name= "Pin", vocab="pbdb") 
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_strata}
\alias{pbdb_strata}
\title{pbdb_strata}
\usage{
pbdb_strata (...)
}
\arguments{
\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/strata/list}. Eg:
\itemize{
  \item \emph{name}: a full or partial name. You can
    use \% and _ as wildcards, but the
    query will be very slow if you put a wildcard at the beginning
  \item \emph{rank}: returns only strata of the specified 
    rank: formation, group or member.
  \item \emph{lngmin}: numeric. The longitude boundaries will be normalized to fall
    between -180 and 180. Note that if you specify lngmin then you must also
    specify lngmax. Returns only records whose geographic location falls within
    the given bounding box (defined by lngmin, lngmax, latmin, latmax). It
    generate two adjacent bounding boxes if the range crosses the antimeridian. 
  \item \emph{lngmax}: numeric. The longitude boundaries will be normalized to fall
    between -180 and 180.
  \item \emph{latmin}: numeric. between -90 and 90. Note that if you specify latmin
    then you must also specify latmax.
  \item \emph{latmax}: numeric. between -90 and 90.
  \item \emph{loc}: Return only strata associated with some occurrence whose
    geographic location falls within the specified geometry, specified in WKT
    format.
  \item \emph{vocab}: set vocab="pbdb" to show the complete name of the variables (by
    default variables have short 3-letter names)
  \item ...
}}
}
\value{
a dataframe with information from the selected strata
}
\description{
Returns information about geological strata, 
selected by name, rank, and/or geographic location.
}
\examples{
\dontrun{
pbdb_strata (lngmin=0, lngmax=15, latmin=0, latmax=15, rank="formation", vocab="pbdb") 
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_temporal_functions.R
\name{pbdb_orig_ext}
\alias{pbdb_orig_ext}
\title{pbdb_orig_ext}
\usage{
pbdb_orig_ext (data, rank, 
temporal_extent, res, orig_ext,  
colour="#0000FF30", bord="#0000FF", do.plot=TRUE)
}
\arguments{
\item{data}{dataframe with our query to the paleoBD \code{\link{pbdb_occurrences}}. 
Important, it is required to show the name of the families, orders, etc. in the dataframe, 
to do that set: show=c("phylo", "ident") (see example).}

\item{rank}{to set which taxon rank you are interested. By default rank= "species"}

\item{temporal_extent}{vector to set the temporal extent (min, max)}

\item{res}{numeric. to set the intervals of the temporal extent}

\item{orig_ext}{1= origination, 2=extinction.}

\item{colour}{to change the colour of the bars in the plot, skyblue2 by default.}

\item{bord}{to set the colour of the border of the polygon}

\item{do.plot}{TRUE/FALSE (TRUE by default).}
}
\value{
a  dataframe with the 
number of first appearances and extinctions of the selected taxon rank across time, 
and a plot with the first appearances or extinctions of the selected taxon rank across time.
}
\description{
Plots the appearance of new taxa across time.
}
\examples{
\dontrun{
canidae<-  pbdb_occurrences (limit="all", vocab="pbdb",
base_name="Canidae", show=c("phylo", "ident"))

# plot of the evolutive rates.
pbdb_orig_ext (canidae, rank="genus", temporal_extent=c(0, 10), 
res=1, orig_ext=1) 

# plot of the extinction rates.
pbdb_orig_ext (canidae, rank="species", temporal_extent=c(0, 10), 
res=1, orig_ext=2) 
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_temporal_functions.R
\name{pbdb_richness}
\alias{pbdb_richness}
\title{pbdb_richness}
\usage{
pbdb_richness (data, rank, res, temporal_extent, colour, bord, do.plot)
}
\arguments{
\item{data}{dataframe with our query to the paleoBD \code{\link{pbdb_occurrences}}. 
Important, it is required to show the name of the families, orders, etc. in the dataframe, 
to do that
set: show=c("phylo", "ident") (see example).}

\item{rank}{to set which taxon rank you are interested. By default rank= "species"}

\item{res}{numeric. to set the intervals of the temporal extent}

\item{temporal_extent}{vector to set the temporal extent (min, max)}

\item{colour}{to change the colour of the bars in the plot, skyblue2 by default.}

\item{bord}{to set the colour of the border of the polygon}

\item{do.plot}{TRUE/FALSE (TRUE by default).}
}
\value{
a plot and a dataframe with the richness aggregated by the taxon rank in the specified temporal extent and resolution.
}
\description{
Plots the number of the interested.
}
\examples{
\dontrun{
data<-  pbdb_occurrences (limit="all", vocab="pbdb",
base_name="Canidae", show=c("phylo", "ident"))
pbdb_richness (data, rank="species", res=1, temporal_extent=c(0,3))
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_geographic_functions.R
\name{pbdb_map}
\alias{pbdb_map}
\title{pbdb_map}
\usage{
pbdb_map (data, col.int="white" ,pch=19, col.ocean="black", main=NULL,
col.point=c("light blue","blue"), ...)
}
\arguments{
\item{data}{Input dataframe. This dataframe is the output of \code{\link{pbdb_occurrences}} function using the
argument: \code{show = "coords"}. See too: \strong{Details} and \strong{Examples}}

\item{col.int}{The colour of the mainland.}

\item{pch}{See: \code{\link{par}}}

\item{col.ocean}{The colour of the ocean.}

\item{main}{To set the title of the map. See: \code{\link{par}}}

\item{col.point}{Two or more colours. To generate the colour gradient used to show the number of occurrences per cell in map}

\item{...}{Others parameters. See \code{\link{par}} and \code{\link{map}}}
}
\value{
A map showing the distribution of the fossil records, with the points with a color gradient, according to the number of occurrences per cell.
}
\description{
Maps the fossil records
}
\details{
The function opens a new window for the map

\strong{CAUTION!} The argument \code{show = "coords"} in \code{\link{pbdb_occurrences}} function is required. 
We recommend the use of a cairo device (\code{\link{X11}}) for better visualization of the graphs. See \strong{Examples}
}
\examples{
\dontrun{
data<- pbdb_occurrences (limit="all", vocab= "pbdb",
base_name="Canis", show="coords")
X11(width=12, height=8)
pbdb_map(data)
pbdb_map(data,pch=1)
pbdb_map(data,pch=19,col.point=c("pink","red"), col.ocean="light blue",
main="canis")
}

}
\seealso{
See \code{\link{pbdb_occurrences}}, \code{\link{map}}, \code{\link{par}} and \code{\link{colors}} help pages
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_ref_occurrences}
\alias{pbdb_ref_occurrences}
\title{pbdb_ref_occurrences}
\usage{
pbdb_ref_occurrences (...)
}
\arguments{
\item{...}{arguments passed to the API. See all available arguments in
\url{http://paleobiodb.org/data1.1/occs/refs}
\itemize{
  \item \emph{author} select only references for which any of the authors 
    matches the specified name
  \item \emph{year} select only references published in the specified year
  \item \emph{pubtitle} select only references that involve the specified 
    publication
  \item \emph{order} specifies the order in which the results are returned. You can
    specify multiple values separated by commas, and each value may be appended
    with .asc or .desc. Accepted values are: author, year, pubtitle, created,
    modified, rank.
}}
}
\value{
a dataframe with the information about the references 
that match the query
}
\description{
Returns information about the bibliographic references 
associated with fossil occurrences from the database.
}
\details{
Go to \code{\link{pbdb_occurrences}} to see an explanation about the main 
filtering parameters.
}
\examples{
\dontrun{
pbdb_ref_occurrences (vocab="pbdb", 
base_name="Canis", year=2000)
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_ref_collections}
\alias{pbdb_ref_collections}
\title{pbdb_ref_collections}
\usage{
pbdb_ref_collections (...)
}
\arguments{
\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/colls/refs}. Eg:
\itemize{
  \item \emph{id} comma-separated list of collection identifiers
  \item \emph{author} select only references for which any of the authors 
    matches the specified name
  \item \emph{year} select only references published in the specified year
  \item \emph{pubtitle} select only references that involve the specified 
    publication
  \item \emph{order} specifies the order in which the results are returned. You can
    specify multiple values separated by commas, and each value may be appended
    with .asc or .desc. Accepted values are: author, year, pubtitle, created,
    modified, rank.
  \item ...
}}
}
\value{
a dataframe with the information about the references that match the query
}
\description{
Returns information about the references from which the selected collection data were entered.
}
\examples{
\dontrun{
  pbdb_ref_collections (id=1)
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_temporal_functions.R
\name{pbdb_temporal_resolution}
\alias{pbdb_temporal_resolution}
\title{pbdb_temporal_resolution}
\usage{
pbdb_temporal_resolution (data, do.plot=TRUE)
}
\arguments{
\item{data}{dataframe with our query to the paleoBD \code{\link{pbdb_occurrences}}}

\item{do.plot}{TRUE/FALSE. To show a frequency plot of the data (TRUE by default).}
}
\value{
a plot and a list with a summary of the temporal resolution of the data 
(min, max, 1st and 3rd quartils, median and mean), and the temporal resolution of each fossil record (Ma).
}
\description{
to show the temporal resolution of the fossil data
}
\examples{
\dontrun{
data<- pbdb_occurrences (taxon_name= "Canidae", interval= "Quaternary")
pbdb_temporal_resolution (data)
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_taxa_auto}
\alias{pbdb_taxa_auto}
\title{pbdb_taxa_auto}
\usage{
pbdb_taxa_auto (...)
}
\arguments{
\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/taxa/auto_doc.html}. Eg:
\itemize{
  \item \emph{name}: a partial name or prefix. 
    It must have at least 3 significant characters,
    and may include both a genus
    (possibly abbreviated) and a species.
  \item \emph{limit}: set the limit to the number of matches
  \item ...
}}
}
\value{
a dataframe with information about the matches 
(taxon rank and number of occurrences in the database)
}
\description{
Returns a list of names matching the given prefix or partial name.
}
\examples{
\dontrun{
pbdb_taxa_auto (name="Cani", limit=10) 
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_references}
\alias{pbdb_references}
\title{pbdb_references}
\usage{
pbdb_references (...)
}
\arguments{
\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/refs/list}. Eg:
\itemize{
  \item \emph{author} select only references for which any of the authors 
    matches the specified name
  \item \emph{year} select only references published in the specified year
  \item \emph{pubtitle} select only references that involve the specified 
    publication
  \item \emph{order} specifies the order in which the results are returned. You can
    specify multiple values separated by commas, and each value may be appended
    with .asc or .desc. Accepted values are: author, year, pubtitle, created,
    modified, rank.
  \item ...
}}
}
\value{
a dataframe with the information about the references that match the query
}
\description{
Returns information about multiple references, selected according to the parameters you provide.
}
\examples{
\dontrun{
  pbdb_references (author= "Polly")
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_collection}
\alias{pbdb_collection}
\title{pbdb_collection}
\usage{
pbdb_collection (id, ...)
}
\arguments{
\item{id}{identifier of the collection. This parameter is required.}

\item{...}{additional arguments passed to the API. See all available arguments in
\url{http://paleobiodb.org/data1.1/colls/single}. Eg:
\itemize{
  \item \emph{vocab}: set vocab="pbdb" to show the complete name of the variables (by
    default variables have short 3-letter names)
  \item \emph{show}: show extra variables
  \item ...
}}
}
\value{
a dataframe with a single occurrence
}
\description{
Returns information about a single collection record from 
the Paleobiology Database.
}
\details{
Go to \code{\link{pbdb_occurrences}} to see an explanation about 
the main filtering parameters.
}
\examples{
\dontrun{
pbdb_collection (id=1003, vocab="pbdb", show="loc")

}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_reference}
\alias{pbdb_reference}
\title{pbdb_reference}
\usage{
pbdb_reference (id, ...)
}
\arguments{
\item{id}{identifier of the reference. This parameter is required.}

\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/refs/single}. Eg:
\itemize{
  \item \emph{vocab}: set vocab="pbdb" to show the complete name of the variables (by
    default variables have short 3-letter names)
  \item ...
}}
}
\value{
a dataframe with a single reference
}
\description{
Returns information about a single reference, selected by identifier.
Go to \code{\link{pbdb_occurrences}} to see an explanation about the main filtering parameters
}
\examples{
\dontrun{
pbdb_collection (id=1003, vocab="pbdb", show="loc")
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_taxon}
\alias{pbdb_taxon}
\title{pbdb_taxon}
\usage{
pbdb_taxon (...)
}
\arguments{
\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/taxa/single}. Eg:
\itemize{
  \item \emph{name}: returns information about the most fundamental 
    taxonomic name matching this string. 
    The \% and _ characters may be used as wildcards.
  \item ...
}}
}
\value{
a dataframe with information from a single taxon
}
\description{
Returns information about a single taxonomic name, 
identified either by name or by identifier.
}
\examples{
\dontrun{
pbdb_taxon (name="Canis", vocab="pbdb", 
show=c("attr", "app", "size"))

}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_ref_taxa}
\alias{pbdb_ref_taxa}
\title{pbdb_ref_taxa}
\usage{
pbdb_ref_taxa (...)
}
\arguments{
\item{...}{arguments passed to the API. See all available arguments in
\url{http://paleobiodb.org/data1.1/taxa/refs}
\itemize{
  \item \emph{name}: returns information about the most fundamental 
    taxonomic name matching this string. 
    The \% and _ characters may be used as wildcards.
  \item \emph{id}: returns information about the taxonomic name 
    corresponding to this identifier. You may not specify both 
    name and id in the same query.
  \item \emph{exact}: if this parameter is specified, then the taxon exactly 
    matching the specified name or identifier is selected, 
    rather than the senior synonym which is the default.
  \item \emph{show}: show extra variables
  \item \emph{rel}: set rel="synonyms" to select all synonyms of 
    the base taxon or taxa; rel="children" to select the 
    taxa immediately contained within the base taxon or taxa; 
    rel="common_ancestor" to select the most specific taxon 
    that contains all of the base taxa.
  \item \emph{extant}: TRUE/FALSE to select extant/extinct taxa.
}}
}
\value{
a dataframe with references from a list of taxa
}
\description{
This URL path returns information about the source references associated
with taxa in the Paleobiology Database. You can use the same parameters 
that are available with pbdb_taxa, but Reference records are returned 
instead of Taxon records. One record is returned per reference, 
even if it is associated with multiple taxa.
}
\examples{
\dontrun{
  pbdb_ref_taxa (name="Canidae", vocab="pbdb", show=c("attr", "app", "size", "nav")) 
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_collections}
\alias{pbdb_collections}
\title{pbdb_collections}
\usage{
pbdb_collections (...)
}
\arguments{
\item{...}{documentation for all the parameters is available 
in http://paleobiodb.org/data1.1/colls/list
go to \code{\link{pbdb_occurrences}} to see an explanation about 
the main filtering parameters}
}
\value{
a dataframe with the collections that match the query
}
\description{
Returns information about multiple collections, selected 
according to the parameters you provide.
}
\examples{
\dontrun{
pbdb_collections (base_name="Cetacea", interval="Miocene")
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/paleobioDB-package.R
\docType{package}
\name{paleobioDB}
\alias{paleobioDB}
\alias{paleobioDB-package}
\title{paleobioDB: An R-package for downloading, visualizing and processing data from the Paleobiology Database}
\description{
We have developed paleobioDB, an R-package designed to make easy and flexible queries 
of the Paleobiology Database, as well as to visualize and download selected data. 
This package will make it easy to access paleontological data in a way that should 
allow those data to be further analyzed, 
including via packages and libraries available in R.
}
\details{
We programmed two different groups of functions. 
First, we developed a set of general and flexible functions to wrap the 
\href{http://paleobiodb.org/data1.1/}{PaleobioDB API}.
These functions connect R with each of the endpoints of the PaleobioDB API. 
Second, based on these base functions, we programmed a second set of functions 
intended to explore and visualize the fossil 
occurrences in their geographic, temporal and taxonomic dimensions.

\tabular{ll}{
Package: \tab paleobioDB\cr
Type: \tab Package\cr
Version: \tab 0.4\cr
Date: \tab 2015-07-16\cr
License: \tab GPL-2\cr
}
}
\examples{
\dontrun{

canidae<-  pbdb_occurrences (limit="all", base_name="canidae", 
interval="Quaternary", show=c("coords", "phylo", "ident"))

## to explore the number of subtaxa
pbdb_subtaxa (canidae)

## to explore the temporal resolution of the fossil records
pbdb_temporal_resolution (canidae)

## returns a dataframe and a plot with the temporal span 
##  of the species, genera, etc.
pbdb_temp_range (canidae, rank= "genus", names=FALSE)

## returns a dataframe and a plot showing the species, genera, etc. 
richness across time
pbdb_richness (canidae, rank= "species", 
temporal_extent= c (0,10), res= 1)

## returns a dataframe and a plot showing the evolutionary 
and extinction rates across time

## evolutionary rates= evo_ext=1
pbdb_orig_ext (canidae, rank="species", temporal_extent=c(0, 10),
res=1, orig_ext=1)

## extinction rates= evo_ext=2
pbdb_orig_ext (canidae, rank="species", temporal_extent=c(0, 10),
              res=1, orig_ext=2)

## maps the fossil occurrences
pbdb_map (canidae, main = "Canidae", pch= 19, cex=0.7)

## maps the sampling effort
pbdb_map_occur (canidae, res= 5)

## maps the species, genera, etc. richness
pbdb_map_richness (canidae, rank="species", res= 5)

}

}
\author{
Sara Varela \email{svarela@paleobiogeography.org}

Javier Gonzalez \email{javigzz@yahoo.es}

Luciano Fabris Sgarbi \email{luciano.f.sgarbi@gmail.com}
}
\references{
Sara Varela, Javier Gonzalez-Hernandez,
Luciano Fabris Sgarbi, Charles Marshall, Mark D. Uhen, 
Shanan Peters, Michael McClennen, 2014. paleobioDB: 
an R-package for downloading, visualizing and processing 
data from the Paleobiology Database (under review)
}
\seealso{
{
\url{http://paleobiodb.org}
}
}
\keyword{package}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_scales}
\alias{pbdb_scales}
\title{pbdb_scales}
\usage{
pbdb_scales(...)
}
\arguments{
\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/scales/list}. Eg:
\itemize{
  \item \emph{vocab}: set vocab="pbdb" to show the complete name of the variables (by
    default variables have short 3-letter names)
  \item ...
}}
}
\value{
a dataframe with information from the selected scales
}
\description{
Returns information about multiple time scales.
}
\examples{
\dontrun{
## Get a dataframe with all the scales available in PBDB 
## by setting no ids
pbdb_scales ()
}

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_occurrence}
\alias{pbdb_occurrence}
\title{pbdb_occurrence}
\usage{
pbdb_occurrence (id, ...)
}
\arguments{
\item{id}{identifier of the occurrence. This parameter is required}

\item{...}{arguments passed to the API. See all available arguments in
\url{http://paleobiodb.org/data1.1/occs/single}. Eg:
\itemize{
 \item \emph{vocab}: set vocab="pbdb" to show the complete name of the variables (by
   default variables have short 3-letter names)
}}
}
\value{
a dataframe with a single occurrence
}
\description{
Returns information about a single occurrence record from the Paleobiology 
Database.
}
\details{
Documentation for all the parameters is available at 
http://paleobiodb.org/data1.1/occs/single. In the parameter list above, we
describe the most common filters that paleontologists and ecologists might
use.
}
\examples{
\dontrun{
pbdb_occurrence (id=1001)
pbdb_occurrence (id=1001, vocab="pbdb", show="coords")
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_interval}
\alias{pbdb_interval}
\title{pbdb_interval}
\usage{
pbdb_interval (id, ...)
}
\arguments{
\item{id}{identifier of the temporal interval. This parameter is required.}

\item{...}{additional arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/intervals/single}. Eg:
\itemize{
  \item \emph{vocab}: set vocab="pbdb" to show the complete name of the variables (by
    default variables have short 3-letter names)
}}
}
\value{
a dataframe with information from a single 
temporal interval
}
\description{
Returns information about a single interval, selected by identifier.
}
\examples{
\dontrun{
pbdb_interval (id=1, vocab="pbdb")
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_occurrences}
\alias{pbdb_occurrences}
\title{pbdb_occurrences}
\usage{
pbdb_occurrences(...)
}
\arguments{
\item{...}{arguments passed to the API. See all available arguments in
\url{http://paleobiodb.org/data1.1/occs/list}
\itemize{
  \item \emph{limit}: sets limit to "all" to download all the occurrences. 
    By default the limit is 500. 
  \item \emph{taxon_name}: Return only records associated with the 
    specified taxonomic name(s). 
    You may specify multiple names, separated by commas.
  \item \emph{base_name}:  Return records associated with the specified 
    taxonomic name(s) and any of their children (e.g. base_name="Canis" will 
  \item \emph{lngmin}: numeric. The longitude boundaries will be normalized 
    to fall between -180 and 180. Note that if you specify 
    lngmin then you must also specify lngmax. 
    Returns only records whose geographic location falls 
    within the given bounding box (defined by lngmin, lngmax, 
    latmin, latmax).
    It generates two adjacent bounding boxes if the range crosses
    the antimeridian. 
  \item \emph{lngmax}: numeric. The longitude boundaries will be normalized 
    to fall between -180 and 180.
  \item \emph{latmin}: numeric. between -90 and 90. 
    Note that if you specify latmin then you must also specify latmax.
  \item \emph{latmax}: numeric. between -90 and 90.
  \item \emph{min_ma}: return only records whose temporal 
    locality is at least this old, specified in Ma.
  \item \emph{max_ma}: return only records whose temporal 
    locality is at most this old, specified in Ma.
  \item \emph{interval}: return only records whose temporal 
    locality falls within the named geologic time interval 
    (e.g. "Miocene").
  \item \emph{continent}: return only records whose geographic 
    location falls within the specified continent(s). 
  \item \emph{show}: to show extra variables (e.g. coords, phylo, ident)
}}
}
\value{
a dataframe with the species occurrences
}
\description{
Returns information about species occurrence records stored in the
Paleobiology Database.
}
\details{
Documentation for all the parameters is available 
at \url{http://paleobiodb.org/data1.1/occs/list}. We describe the most common
filters that paleontologists and ecologists might use in the parameter list above.
}
\examples{
\dontrun{
pbdb_occurrences (id=c(10, 11), show=c("coords", "phylo", "ident")) 
pbdb_occurrences (limit="all", vocab= "pbdb", 
taxon_name="Canis", show=c("coords", "phylo", "ident"))
pbdb_occurrences (limit="all", vocab= "pbdb", 
base_name="Canidae", show=c("coords", "phylo", "ident"))
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_intervals}
\alias{pbdb_intervals}
\title{pbdb_intervals}
\usage{
pbdb_intervals (...)
}
\arguments{
\item{...}{arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/intervals/list}. Eg:
\itemize{
  \item \emph{min_ma}: return only intervals that are at least this old
  \item \emph{max_ma}: return only intervals that are at most this old
  \item \emph{order}: return the intervals in order starting as specified. 
    Possible values include older, younger. Defaults to younger
  \item \emph{vocab}: set vocab="pbdb" to show 
    the complete name of the variables (by
    default variables have short 3-letter names)
  \item ...
}}
}
\value{
a dataframe with information from several temporal intervals
}
\description{
Returns information about multiple intervals, 
selected according to the parameters you provide.
}
\examples{
\dontrun{
pbdb_intervals (min_ma= 0, max_ma=2, vocab="pbdb") 
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_geographic_functions.R
\name{pbdb_map_richness}
\alias{pbdb_map_richness}
\title{pbdb_map_richness}
\usage{
pbdb_map_richness (data, rank="species", do.plot=TRUE, res=5,
col.int="white", col.ocean="black",
col.rich=c("light blue","blue"),...)
}
\arguments{
\item{data}{Input dataframe. This dataframe is the output of \code{\link{pbdb_occurrences}} function using the argument: \code{show = c("phylo", "coords", "ident")}. See too: \strong{Details} and \strong{Examples}}

\item{rank}{To set which taxon rank you are interested for calculate richness. The options are: "species", "genus", "family", "order", "class" or "phylum")}

\item{do.plot}{Logical; \code{TRUE} the function returns a RasterLayer and a plot.}

\item{res}{The resolution of the RasterLayer object (in decimal degrees). See: \code{\link{raster}}}

\item{col.int}{The colour of the mainland}

\item{col.ocean}{The colour of the ocean}

\item{col.rich}{Two or more colours. To generate the colour gradient used to show the richness per cell in map}

\item{...}{Others parameters. See \code{\link{par}} and \code{\link{map}}}
}
\value{
A RasterLayer object and a plot with richness of species, genera, families, etc. per cell. This RasterLayer object have the resolution controlled by
the argument \code{res}. The default is \code{res=1}.
}
\description{
Creates a RasterLayer object and a plot with richness of species, genera, families, etc. per cell.
}
\details{
\strong{CAUTION!} The argument \code{show = "coords"} in \code{\link{pbdb_occurrences}} function is required. 
We recommend the use of a cairo device (\code{\link{X11}}) for better visualization of the graphs. See \strong{Examples}
}
\examples{
\dontrun{
data<- pbdb_occurrences (limit=1000, vocab= "pbdb", base_name="mammalia",
show=c("phylo","coords","ident"))
X11(width=13, height=7.8)
pbdb_map_richness (data,res=8,rank="genus")
pbdb_map_richness (data,res=8,rank="family")
## to obtain the raster file and not plot the map
pbdb_map_richness (data,res=8,rank="family",do.plot=F)
}

}
\seealso{
See \code{\link{pbdb_occurrences}}, \code{\link{map}}, \code{\link{par}} and \code{\link{colors}} help pages
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_scale}
\alias{pbdb_scale}
\title{pbdb_scale}
\usage{
pbdb_scale (id, ...)
}
\arguments{
\item{id}{identifier of the temporal interval. This parameter is required.}

\item{...}{additional arguments passed to the API. See
documentation for accepted parameters in
\url{http://paleobiodb.org/data1.1/scales/single}. Eg:
\itemize{
  \item \emph{vocab}: set vocab="pbdb" to show the complete name of the variables (by
    default variables have short 3-letter names)
  \item ...
}}
}
\value{
a dataframe with information from a single scale
}
\description{
Returns information about a single time scale, selected by 
identifier.
}
\examples{
\dontrun{
  pbdb_scale (id=1, vocab="pbdb")
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_geographic_functions.R
\name{pbdb_map_occur}
\alias{pbdb_map_occur}
\title{pbdb_map_occur}
\usage{
pbdb_map_occur (data, res=5, col.int="white", col.ocean="black",
col.eff=c("light blue","blue"), do.plot=TRUE, ...)
}
\arguments{
\item{data}{Input dataframe. This dataframe is the output of \code{\link{pbdb_occurrences}} function using the argument: \code{show="coords"}. See too: \strong{Details} and \strong{Examples}}

\item{res}{the resolution of the RasterLayer object (in decimal degrees). See: \code{\link{raster}}}

\item{col.int}{The colour of the mainland}

\item{col.ocean}{The colour of the ocean}

\item{col.eff}{Two or more colours. To generate the colour gradient used to show the number of occurrences per cell in map}

\item{do.plot}{Logical; \code{TRUE} the function returns a RasterLayer and a plot.}

\item{...}{Others parameters. See \code{\link{par}} and \code{\link{map}}}
}
\value{
A RasterLayer object and a plot with the sampling effort (number of fossil records per cell). This RasterLayer object have the resolution controlled by the argument \code{res}. The deflaut is \code{res=1}.
}
\description{
Creates a RasterLayer object and a plot of the sampling effort (number of fossil records per cell).
}
\details{
\strong{CAUTION!} The argument \code{show = "coords"} in \code{\link{pbdb_occurrences}} function is required. 
We recommend the use of a cairo device (\code{\link{X11}}) for better visualization of the graphs. See \strong{Examples}
}
\examples{
\dontrun{
data<- pbdb_occurrences (limit="all", vocab= "pbdb", base_name="Canis",
show="coords")
X11(width=13, height=7.8)
pbdb_map_occur (data,res=2)
## to obtain the raster file without plotting it
pbdb_map_occur (data,res=3,do.plot=F)
}

}
\seealso{
See \code{\link{pbdb_occurrences}}, \code{\link{map}}, \code{\link{par}} and \code{\link{colors}} help pages
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_taxa}
\alias{pbdb_taxa}
\title{pbdb_taxa}
\usage{
pbdb_taxa (...)
}
\arguments{
\item{...}{arguments passed to the API. See all available arguments in
\url{http://paleobiodb.org/data1.1/taxa/list}
\itemize{
  \item \emph{name}: returns information about the most fundamental 
    taxonomic name matching this string. 
    The \% and _ characters may be used as wildcards.
  \item \emph{id}: return information about the taxonomic name 
    corresponding to this identifier. You may not specify both 
    name and id in the same query.
  \item \emph{exact}: if this parameter is specified, then the taxon exactly 
    matching the specified name or identifier is selected, 
    rather than the senior synonym which is the default.
  \item \emph{show}: to show extra variables: \emph{attr} 
    the attribution of this taxon (author and year); 
    \emph{app} the age of first and last appearance of this taxon 
    from the occurrences recorded in this database; 
    \emph{size} the number of subtaxa appearing in this database; 
    \emph{nav} additional information for the PBDB Navigator taxon browser
  \item \emph{rel}: set rel="synonyms" to select all synonyms of 
    the base taxon or taxa; rel="children" to select the 
    taxa immediately contained within the base taxon or taxa; 
    rel="common_ancestor" to select the most specific taxon 
    that contains all of the base taxa.
  \item \emph{extant}: TRUE/FALSE to select extant/extinct taxa.
}}
}
\value{
a dataframe with information from a list of taxa
}
\description{
Returns information about multiple taxonomic names.  This function can be
used to query for all of the children or parents of a given taxon, among
other operations.
}
\examples{
\dontrun{
pbdb_taxa (name="Canidae", vocab="pbdb", 
show=c("attr", "app", "size", "nav"))
pbdb_taxa (id =c(10, 11), vocab="pbdb", 
show=c("attr", "app", "size", "nav"))
pbdb_taxa (id =c(10, 11), vocab="pbdb", 
show=c("attr", "app", "size", "nav"), rel="common_ancestor")
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_temporal_functions.R
\name{pbdb_temp_range}
\alias{pbdb_temp_range}
\title{pbdb_temp_range}
\usage{
pbdb_temp_range (data, rank, col = "#0000FF", 
names = TRUE, do.plot =TRUE)
}
\arguments{
\item{data}{dataframe with our query to the paleoBD \code{\link{pbdb_occurrences}}. 
Important, it is required to show the name of the families, orders, etc. in the dataframe, 
to do that
set: show=c("phylo", "ident") (see example).}

\item{rank}{to set which taxon rank you are interested.}

\item{col}{to change the colour of the bars in the plot, skyblue2 by default.}

\item{names}{TRUE/FALSE (TRUE by default). To include or not the name of the taxa in the plot}

\item{do.plot}{TRUE/FALSE (TRUE by default).}
}
\value{
a plot and a dataframe with the time span of the taxa selected (species, genus, etc.)
}
\description{
constructs a plot and a dataframe with the temporal range of the taxa (species, genera, families, etc.) within in a selected higher taxon.
}
\examples{
\dontrun{
canis_quaternary<- pbdb_occurrences (limit="all", base_name="Canis", 
                 interval="Quaternary", show=c("coords", "phylo", "ident"))
pbdb_temp_range (canis_quaternary, rank="species", names=FALSE)
}
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_taxonomic_functions.R
\name{pbdb_subtaxa}
\alias{pbdb_subtaxa}
\title{pbdb_subtaxa}
\usage{
pbdb_subtaxa (data, do.plot, col)
}
\arguments{
\item{data}{dataframe with our query to the 
paleoBD \code{\link{pbdb_occurrences}}}

\item{do.plot}{by default this function make a plot to 
visualize the distribution of taxa. Set to FALSE to skip the plot.}

\item{col}{set the colour of the histogram. skyblue2 by default.}
}
\value{
a plot and a dataframe with the number of subtaxa in the data.
}
\description{
count the number of subtaxa within a given taxa. 
e.g. number of species within a genus.
}
\examples{
\dontrun{
canidae_quat<-  pbdb_occurrences (limit="all", 
base_name="Canidae",  interval="Quaternary", 
show=c("coords", "phylo", "ident"))
pbdb_subtaxa (canidae_quat)
}

}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pbdb_querys.R
\name{pbdb_collections_geo}
\alias{pbdb_collections_geo}
\title{pbdb_collections_geo}
\usage{
pbdb_collections_geo (...)
}
\arguments{
\item{...}{documentation for all the parameters is 
available in http://paleobiodb.org/data1.1/colls/summary
go to \code{\link{pbdb_occurrences}} to see an explanation about 
the main filtering parameters}
}
\value{
a dataframe with the collections that match the query
}
\description{
This path returns information about geographic clusters 
of collections from the Paleobiology Database. 
These clusters are defined in order to facilitate the 
generation of maps at low resolutions. 
You can make a config request via 
http://paleobiodb.org/data1.1/config
in order to get a list of the available summary levels.
}
\examples{
\dontrun{
pbdb_collections_geo (vocab="pbdb", lngmin=0.0, 
lngmax=15.0, latmin=0.0, latmax=15.0, level=2)
}
}

