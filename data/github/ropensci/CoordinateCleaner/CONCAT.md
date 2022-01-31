# CoordinateCleaner v2.0-20
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/CoordinateCleaner)](https://cranlogs.r-pkg.org:443/badges/CoordinateCleaner)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/CoordinateCleaner)](https://cranlogs.r-pkg.org:443/badges/grand-total/CoordinateCleaner)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/CoordinateCleaner)](https://cranlogs.r-pkg.org:443/badges/CoordinateCleaner)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2539408.svg)](https://doi.org/10.5281/zenodo.2539408)
[![rOpenSci peer-review](https://badges.ropensci.org/210_status.svg)](https://github.com/ropensci/software-review/issues/210)

Automated flagging of common spatial and temporal errors in biological and palaeontological collection data, for the use in conservation, ecology and palaeontology. Specifically includes tests for

* General coordinate validity
* Country and province centroids
* Capital coordinates
* Coordinates of biodiversity institutions
* Spatial outliers
* Temporal outliers
* Coordinate-country discordance
* Duplicated coordinates per species
* Assignment to the location of the GBIF headquarters
* Urban areas
* Seas
* Plain zeros
* Equal longitude and latitude
* Rounded coordinates
* DDMM to DD.DD coordinate conversion errors
* Large temporal uncertainty (fossils)
* Equal minimum and maximum ages (fossils)
* Spatio-temporal outliers (fossils)

CoordinateCleaner can be particularly useful to improve data quality when using data from GBIF (e.g. obtained with [rgbif]( https://github.com/ropensci/rgbif)) or the Paleobiology database (e.g. obtained with [paleobioDB](https://github.com/ropensci/paleobioDB)) for historical biogeography (e.g. with [BioGeoBEARS](https://CRAN.R-project.org/package=BioGeoBEARS) or [phytools](https://CRAN.R-project.org/package=phytools)), automated conservation assessment (e.g. with [speciesgeocodeR](https://github.com/azizka/speciesgeocodeR/wiki) or [conR](https://CRAN.R-project.org/package=ConR)) or species distribution modelling (e.g. with [dismo](https://CRAN.R-project.org/package=dismo) or [sdm](https://CRAN.R-project.org/package=sdm)). See [scrubr](https://github.com/ropensci/scrubr) and [taxize](https://github.com/ropensci/taxize) for complementary taxonomic cleaning or [biogeo](https://github.com/cran/biogeo) for correcting spatial coordinate errors.

See [News](https://github.com/ropensci/CoordinateCleaner/blob/master/NEWS.md) for update information.

# Installation
## Stable from CRAN

```r
install.packages("CoordinateCleaner")
library(CoordinateCleaner)
```

## Developmental from GitHub
```r
devtools::install_github("ropensci/CoordinateCleaner")
library(CoordinateCleaner)
```

# Usage
A simple example:

```r
# Simulate example data
minages <- runif(250, 0, 65)
exmpl <- data.frame(species = sample(letters, size = 250, replace = TRUE),
                    decimallongitude = runif(250, min = 42, max = 51),
                    decimallatitude = runif(250, min = -26, max = -11),
                    min_ma = minages,
                    max_ma = minages + runif(250, 0.1, 65),
                    dataset = "clean")

# Run record-level tests
rl <- clean_coordinates(x = exmpl)
summary(rl)
plot(rl)

# Dataset level 
dsl <- clean_dataset(exmpl)

# For fossils
fl <- clean_fossils(x = exmpl,
                          taxon = "species",
                          lon = "decimallongitude", 
                          lat = "decimallatitude")
summary(fl)

# Alternative example using the pipe
library(tidyverse)

cl <- exmpl %>%
  cc_val()%>%
  cc_cap()%>%
  cd_ddmm()%>%
  cf_range(lon = "decimallongitude", 
           lat = "decimallatitude", 
           taxon  ="species")
```

# Documentation
Pipelines for cleaning data from the Global Biodiversity Information Facility (GBIF) and the Paleobiology Database (PaleobioDB) are available in [here](https://ropensci.github.io/CoordinateCleaner/articles/).


# Contributing
See the [CONTRIBUTING](https://github.com/ropensci/CoordinateCleaner/blob/master/CONTRIBUTING.md) document.

# Citation
Zizka A, Silvestro D, Andermann T, Azevedo J, Duarte Ritter C, Edler D, Farooq H, Herdean A, Ariza M, Scharn R, Svanteson S, Wengtrom N, Zizka V & Antonelli A (2019) CoordinateCleaner: standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution, 10(5):744-751, doi:10.1111/2041-210X.13152, https://github.com/ropensci/CoordinateCleaner

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

CoordinateCleaner 2.0-20 (2021-10-08)
=========================

### MINOR IMPROVEMENTS
  * Fixed typos in vignette
  * Fixed CRAN warning with rebuilding vignettes
  * Removed defunct functions
  * Added wordlist

CoordinateCleaner 2.0-19 (2020-10-13)
=========================

### MINOR IMPROVEMENTS
  * Adaption of the Description file style.
  
CoordinateCleaner 2.0-18 (2020-10-09)
=========================

### MINOR IMPROVEMENTS
  * Precomputed vignettes that need internet access
  
  
CoordinateCleaner 2.0-17 (2020-08-18)
=========================

### MINOR IMPROVEMENTS
  * Fixed a bug that occasionally caused invalid polygons for cc_coun and cc_sea
  * Added references for country centroids from other sources

CoordinateCleaner 2.0-16 (2020-06-19)
=========================

### MINOR IMPROVEMENTS
  * added an option to adapt the geographic extent to the plotting method for objects of class "spatialvalid"
  * fixed a bug with the handling of the outlier test in clean_coordinates
  * improved the handling of large data.frames in the plotting method for objects of class "spatialvalid"
  * readapted the default country columns name to the rnaturalearth column names


CoordinateCleaner 2.0-15 (2020-05-04)
=========================

### MINOR IMPROVEMENTS
  * adapted the url format in the description files


CoordinateCleaner 2.0-14 (2020-05-04)
=========================

### MINOR IMPROVEMENTS
  * fixed a bug in cc_cen when setting an alternative reference
  * added the ref_col argument to cc_coun to customize the column with ISO codes in the reference data
  * adapted code to changes in sp and rgdal
  * defunct CleanCoordinates, CleanCoordinatesDS, and CleanCoordinatesFOS
  * fixed issue with input data.frame with unordered rownames in cc_outl
  * fixed the 'ras not found' bug in cc_outl
  
  

CoordinateCleaner 2.0-13 (2019-06-18)
=========================

### MINOR IMPROVEMENTS
  * addressed the "ras not found" bug in cc_outl
  
CoordinateCleaner 2.0-12 (2019-05-2)
=========================

### MINOR IMPROVEMENTS
  * improved documentation of cc_outl
  * improved handling of rownames in cc_outl
  
  
CoordinateCleaner 2.0-11 (2019-04-24)
=========================

### MINOR IMPROVEMENTS
  * changes to the description file

CoordinateCleaner 2.0-10 (2019-04-23)
=========================

### MINOR IMPROVEMENTS
  * improved error handling by cc_sea and cc_urb, in case the default reference cannot be obtained from the web
  * added a reference for the methodology to the description file


CoordinateCleaner 2.0-9 (2019-04-02)
=========================

### MINOR IMPROVEMENTS
  * recoded cc_outl, and added a thinning argument to account for sampling bias
  * fixed a bug with the cc_outl test, that produced erroneous flags under some settings of mltpl
  * extended the example dataset for the coordinate level-test suite to be more realistic


CoordinateCleaner 2.0-8 (2019-03-21)
=========================

### MINOR IMPROVEMENTS
  * moved vignettes to online documentation
  * added an area column to the countryref dataset
  * fixed some minor spelling issues in the documentation


CoordinateCleaner 2.0-7 (2019-01-22)
=========================

### MINOR IMPROVEMENTS
  * added citation
  * reduced testing time on CRAN
  * improved documentation of the cc_outl function

CoordinateCleaner 2.0-6 (2019-01-16)
=========================

### MINOR IMPROVEMENTS
  * further url fixes

CoordinateCleaner 2.0-5 (2019-01-15)
=========================

### MINOR IMPROVEMENTS
  * fixed broken url to the CIA factbook
  

CoordinateCleaner 2.0-4 (2019-01-14)
=========================

### MINOR IMPROVEMENTS
  * minor bugfix with cc_cap
  * corrected duplicated vignette index entries 
  * updated maintainer email
  
  
CoordinateCleaner 2.0-3 (2018-10-22)
=========================

### MINOR IMPROVEMENTS
  * removed convenience functionality to only download data from rnaturalearth at first use, to comply with CRAN guidelines


CoordinateCleaner 2.0-2 (2018-10-12)
=========================

### MAJOR IMPROVEMENTS

  * tutorial on outlier detection on the bookdown documentation
  * tutorial on using custom gazetteers
  * rasterisation heuristic in cc_outl
  * added sampling correction to cc_outl  
  * added verify option to cc_inst
  * transfer to rOpenSci
  
  
### MINOR IMPROVEMENTS

  * reduced packages size, by switching to data download from rnaturalearth for urbanareas and landmass
  * fixed issue with names of plot.spatialvalid
  * grouped functions on documentation webpage
  * fixed broken links in the help pages
  ' improved documentation structure
  

CoordinateCleaner 2.0-1 (2018-06-08)
=========================

### MAJOR IMPROVEMENTS

  * changed and more consistent naming scheme for the functions
  
### MINOR IMPROVEMENTS
  * fixed typos in Readme
  * set a download from naturalearth as default for cc_urb
  * reduced vignette memory use and size
  * enables sf format for custom references
  * added speedup option for cc_sea
  * added webpage (https://azizka.github.io/CoordinateCleaner/)


CoordinateCleaner 1.2-1 (2018-06-08)
=========================

### MAJOR IMPROVEMENTS

  * Adapted function and argument names consistently to underscore_case
  * Simplified internal code structure of wrapper functions
  
### MINOR IMPROVEMENTS
  * adapted package to rOpenSci reviews
  
### DEPRECATED AND DEFUNCT
  * CleanCoordinates deprecated, replaced by clean_coordinates
  * CleanCoordinatesDS deprecated, replaced by clean_dataset
  * CleanCoordinatesFOS deprecated, replaced by clean_fossils
  * WritePyrate deprecated, replaced by write_pyrate

CoordinateCleaner 1.1-1 (2018-05-15)
=========================

### MINOR IMPROVEMENTS
  * Switched documentation and NAMESPACE generation to roxygen2
  * Switched from sapply to vapply
  * Improved code readability


CoordinateCleaner 1.1-0 (2018-04-08)
=========================

### NEW FEATURES

### MINOR IMPROVEMENTS

  * Adaption of code to rOpenSci guidelines

### BUG FIXES

### DEPRECATED AND DEFUNCT
# CONTRIBUTING #

## Bugs, suggestions or feature requests?

* Submit an issue on the [Issues page](https://github.com/azizka/CoordinateCleaner/issues) - be sure to include R session information and a reproducible example.

## Code contribution

* If you want to contribute to the package - awesome. Please get in touch with [zizka.alexander@gmail.com](mailto:zizka.alexander@gmail.com).# Version 2.0-20

Fixing a warning resulting from unbalanced code chunk delimiters in one of the vignettes in the previous version (2.0-18) on CRAN. 

Additionally, change of maintainer email due to a change in institution.

The NOTE on spelling errors in the DESCRIPTION is spurious in my opinion, since it flags my last name and Latin abbreviation "et al"" from the literature reference.

The package has been removed from CRAN due to slow response time. Response time was slow because I am on parental leave.# CoordinateCleaner v2.0-19
[![Build Status](https://travis-ci.org/ropensci/CoordinateCleaner.svg?branch=master)](https://travis-ci.org/ropensci/CoordinateCleaner)
[![codecov.io](https://codecov.io/github/ropensci/coordinatecleaner/graphs/badge.svg?branch=master)](https://codecov.io/github/ropensci/CoordinateCleaner)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/CoordinateCleaner)](https://cran.r-project.org/package=CoordinateCleaner)
[![downloads](https://cranlogs.r-pkg.org/badges/grand-total/CoordinateCleaner)](http://cranlogs.r-pkg.org/badges/grand-total/CoordinateCleaner)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/CoordinateCleaner)](http://cranlogs.r-pkg.org/badges/CoordinateCleaner)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.2539408.svg)](https://doi.org/10.5281/zenodo.2539408)
[![rOpenSci peer-review](https://badges.ropensci.org/210_status.svg)](https://github.com/ropensci/software-review/issues/210)

**There were some changes to the documentation of CoordinateCleaner with version 2.0-17. You can find three vignettes on cleaning point occurrence records (for instance from www.gbif.org), fossil data and using custom gazetteers with the package now. Additional information on dataset-level cleaning and identifying geographic outliers can be found [here](https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13152).**

Automated flagging of common spatial and temporal errors in biological and palaeontological collection data, for the use in conservation, ecology and palaeontology. Specifically includes tests for

* General coordinate validity
* Country and province centroids
* Capital coordinates
* Coordinates of biodiversity institutions
* Spatial outliers
* Temporal outliers
* Coordinate-country discordance
* Duplicated coordinates per species
* Assignment to the location of the GBIF headquarters
* Urban areas
* Seas
* Plain zeros
* Equal longitude and latitude
* Rounded coordinates
* DDMM to DD.DD coordinate conversion errors
* Large temporal uncertainty (fossils)
* Equal minimum and maximum ages (fossils)
* Spatio-temporal outliers (fossils)

CoordinateCleaner can be particularly useful to improve data quality when using data from GBIF (e.g. obtained with [rgbif]( https://github.com/ropensci/rgbif)) or the Paleobiology database (e.g. obtained with [paleobioDB](https://github.com/ropensci/paleobioDB)) for historical biogeography (e.g. with [BioGeoBEARS](https://CRAN.R-project.org/package=BioGeoBEARS) or [phytools](https://CRAN.R-project.org/package=phytools)), automated conservation assessment (e.g. with [speciesgeocodeR](https://github.com/azizka/speciesgeocodeR/wiki) or [conR](https://CRAN.R-project.org/package=ConR)) or species distribution modelling (e.g. with [dismo](https://CRAN.R-project.org/package=dismo) or [sdm](https://CRAN.R-project.org/package=sdm)). See [scrubr](https://github.com/ropensci/scrubr) and [taxize](https://github.com/ropensci/taxize) for complementary taxonomic cleaning or [biogeo](https://github.com/cran/biogeo) for correcting spatial coordinate errors. You can find a detailed comparison of the functionality of `CoordinateCleaner`, `scrubr`, and `biogeo` [here](https://ropensci.github.io/CoordinateCleaner/articles/comparison_other_software.html).

See [News](https://github.com/ropensci/CoordinateCleaner/blob/master/NEWS.md) for update information.

# Installation
## Stable from CRAN

```r
install.packages("CoordinateCleaner")
library(CoordinateCleaner)
```

## Developmental using devtools
```r
devtools::install_github("ropensci/CoordinateCleaner")
library(CoordinateCleaner)
```

# Usage
A simple example:

```r
# Simulate example data
minages <- runif(250, 0, 65)
exmpl <- data.frame(species = sample(letters, size = 250, replace = TRUE),
                    decimallongitude = runif(250, min = 42, max = 51),
                    decimallatitude = runif(250, min = -26, max = -11),
                    min_ma = minages,
                    max_ma = minages + runif(250, 0.1, 65),
                    dataset = "clean")

# Run record-level tests
rl <- clean_coordinates(x = exmpl)
summary(rl)
plot(rl)

# Dataset level 
dsl <- clean_dataset(exmpl)

# For fossils
fl <- clean_fossils(x = exmpl,
                          taxon = "species",
                          lon = "decimallongitude", 
                          lat = "decimallatitude")
summary(fl)

# Alternative example using the pipe
library(tidyverse)

cl <- exmpl %>%
  cc_val()%>%
  cc_cap()%>%
  cd_ddmm()%>%
  cf_range(lon = "decimallongitude", 
           lat = "decimallatitude", 
           taxon  ="species")
```

# Documentation
Pipelines for cleaning data from the Global Biodiversity Information Facility (GBIF) and the Paleobiology Database (PaleobioDB) are available in [here](https://ropensci.github.io/CoordinateCleaner/articles/).


# Contributing
See the [CONTRIBUTING](https://github.com/ropensci/CoordinateCleaner/blob/master/CONTRIBUTING.md) document.

# Citation
Zizka A, Silvestro D, Andermann T, Azevedo J, Duarte Ritter C, Edler D, Farooq H, Herdean A, Ariza M, Scharn R, Svanteson S, Wengtrom N, Zizka V & Antonelli A (2019) CoordinateCleaner: standardized cleaning of occurrence records from biological collection databases. Methods in Ecology and Evolution, 10(5):744-751, doi:10.1111/2041-210X.13152, https://github.com/ropensci/CoordinateCleaner

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "Cleaning GBIF data for the use in biogeography"
output: 
  rmarkdown::html_vignette:
    fig_caption: true
    number_sections: true
    self_contained: no
bibliography: CoordinateCleaner.bib
csl: apa.csl
vignette: >
  %\VignetteIndexEntry{Cleaning GBIF data for the use in biogeography}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteSuggests{rgbif}
  \usepackage[utf8]{inputenc}
---



# Background
Big data aggregators such as the Global Biodiversity Information Facility (GBIF, www.gbif.org) have vastly increased the public availability of species occurrence records, with GBIF alone comprising more than 800 million records across all taxonomic groups. The data provided via these sources have revolutionized scientific biogeography and are highly valuable for research. However, some issues exist concerning data quality, mostly because these data are comprised from a variety of different collection methods (museum specimens, scientific surveys, citizen science, population counts for conservation purposes and genetic barcoding among others) and different sources (museums, herbaria, collections of individual researchers, citizen science, photo apps) and digitized and edited by various people and algorithms at different points in time and space.

In this tutorial we provide a pipeline on how to clean occurrence records retrieved from GBIF (or any other database) using *CoordinateCleaner* and meta data. The tutorial includes major steps we consider necessary, but by no means is complete and we explicitly encourage you to explore your data further before use. For the tutorial we will use a data set of occurrence records of a single species (lion, *Panthera leo*) downloaded from GBIF. On this example we can gauge the quality of cleaning steps, because we already have a good idea where we expect lions to occur. Of course, usually for multi-species data sets we do not have this kind of information, and that is the whole point of the automated cleaning. You can easily follow the tutorial using your own data instead. For the tutorial we will assume a global macroecological analysis with a resolution of about 100km as downstream analyses. Remember to adjust test sensitivity, if your analyses have a coarser or finer resolution.

With this tutorial you will be able to:

1. Visualize the data and identify potential problems 
2. Use  *CoordinateCleaner* to automatically flag problematic records
3. Use GBIF provided meta-data to improve coordinate quality, tailored to your downstream analyses
4. Use automated cleaning algorithms of *CoordinateCleaner* to identify problematic contributing datasets

# Identifying erroneous coordinates with *CoordinateCleaner*

The `clean_coordinates` function is a wrapper function around all record-level tests of *CoordinateCleaner*. The idea behind these tests is to use geographic gazetteers to identify records that are most likely erroneous (or very imprecise). We based the choice of tests on common problems observed in biological collection databases [@Maldonado2015], including assignment to country centroids, sea coordinate and outliers among others. You can get an overview over the individual tests using `?clean_coordinates` or via the [package vignettes](https://ropensci.github.io/CoordinateCleaner/). This tutorial assumes occurrence data in the format as downloaded from GBIF, for other formats you might need to adapt the column names. You might need to install some of the required packages for the tutorial using `install.packages`.

## Install `CoordinateCleaner`
You can install the latest stable version of CoordinateCleaner from CRAN using `install.packages("CoordinateCleaner")`. Alternatively you can install the latest development version from GitHub using the devtools package. We recommend the latter, to stay up-to-date. Also, make sure to have the latest R version installed.


```r
install.packages("devtools")
library(devtools)

install_github("ropensci/CoordinateCleaner")
```

## Set up libraries and data
You might need to confirm to install the rnaturalearth package when loading `CoordinateCleaner`






```r
library(countrycode)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(rgbif)
library(sp)

#obtain data from GBIF via rgbif
dat <- occ_search(scientificName = "Panthera leo", limit = 5000, hasCoordinate = T)

dat <- dat$data

# names(dat) #a lot of columns

#select columns of interest
dat <- dat %>%
  dplyr::select(species, decimalLongitude, decimalLatitude, countryCode, individualCount,
         gbifID, family, taxonRank, coordinateUncertaintyInMeters, year,
         basisOfRecord, institutionCode, datasetName)

# remove records without coordinates
dat <- dat%>%
  filter(!is.na(decimalLongitude))%>%
  filter(!is.na(decimalLatitude))
```

## Visualize the data on a map

```r
#plot data to get an overview
wm <- borders("world", colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = dat, aes(x = decimalLongitude, y = decimalLatitude),
             colour = "darkred", size = 0.5)+
  theme_bw()
```

![\label{fig:al}Occurrence records for Panthera leo obtained from GBIF.](gbif-clgbif5-1.png)

This map clearly indicates, that we need to prepare the data further, if we want them to represent the current day (or historic) distribution of lions.

## Use *CoordinateCleaner* to automatically flag problematic records

### Option A) Using the `clean_coordinates` wrapper function
As a first step we will run the automatic cleaning algorithm of CoordinateCleaner. The `clean_coordinates` function is a wrapper around a large set of automated cleaning steps to flag errors that are common to biological collections, including: sea coordinates, zero coordinates, coordinate - country mismatches, coordinates assigned to country and province centroids, coordinates within city areas, outlier coordinates and coordinates assigned to biodiversity institutions. You can switch on each test individually using logical flags, modify the sensitivity of most individual tests using the ".rad" arguments, and provide custom gazetteers using the ".ref" arguments. See `?clean_coordinates` for help. To use the country - coordinate mismatch test we need to convert the country from ISO2 to ISO3 format.


```r
#convert country code from ISO2c to ISO3c
dat$countryCode <-  countrycode(dat$countryCode, origin =  'iso2c', destination = 'iso3c')


#flag problems
dat <- data.frame(dat)
flags <- clean_coordinates(x = dat, 
                           lon = "decimalLongitude", 
                           lat = "decimalLatitude",
                           countries = "countryCode",
                           species = "species",
                          tests = c("capitals", "centroids", "equal","gbif", "institutions",
                                    "zeros", "countries")) # most test are on by default
## Testing coordinate validity
## Flagged 0 records.
## Testing equal lat/lon
## Flagged 0 records.
## Testing zero coordinates
## Flagged 1 records.
## Testing country capitals
## Flagged 22 records.
## Testing country centroids
## Flagged 0 records.
## Testing country identity
## Flagged 284 records.
## Testing GBIF headquarters, flagging records around Copenhagen
## Flagged 0 records.
## Testing biodiversity institutions
## Flagged 0 records.
## Flagged 305 of 5000 records, EQ = 0.06.
summary(flags)
plot(flags, lon = "decimalLongitude", lat = "decimalLatitude")
```

![\label{fig:automated}Records flagged by the automated cleaning.](gbif-clgbif6-1.png)

```
##     .val     .equ     .zer     .cap     .cen     .con     .gbf    .inst .summary 
##        0        0        1       22        0      284        0        0      305
```

The automatic test flagged 6.1% of the records. For the purpose of this tutorial we will exclude the flagged records, but in general it is recommendable to explore them further.


```r
#Exclude problematic records
dat_cl <- dat[flags$.summary,]

#The flagged records
dat_fl <- dat[!flags$.summary,]
```

### Option B) Using the magrittr pipe (%>%)
Alternatively, you can run all tests implemented in *CoordinateCleaner* with a individual function and connect them using the magrittr pipe operator, which will directly result in a `data.frame` comprising only cleaned records.


```r
 #to avoid specifying it in each function
names(dat)[2:3] <- c("decimallongitude", "decimallatitude")

clean <- dat%>%
  cc_val()%>%
  cc_equ()%>%
  cc_cap()%>%
  cc_cen()%>%
  cc_coun(iso3 = "countryCode")%>%
  cc_gbif()%>%
  cc_inst()%>%
  cc_sea()%>%
  cc_zero()%>%
  cc_outl()%>%
  cc_dupl()
```

In this way, you can also add the individual test results as columns to your initial data.frame:


```r
dat %>%
    as_tibble() %>% 
    mutate(val = cc_val(., value = "flagged"),
           sea = cc_sea(., value = "flagged"))
```

### Temporal outliers
While the `cc_outl` function identifies geographic outliers, record in GBIF migh also have doubtful temporal information, i.e. for the time of collection, which can be problematic for example for analyses of range dynamics. The `cf_age` function used for fossil cleaning can also be used to check GBIF records for temporal outliers.


```r
flags <- cf_age(x = dat_cl,
                lon = "decimalLongitude",
                lat = "decimalLatitude",
                taxon = "species", 
                min_age = "year", 
                max_age = "year", 
                value = "flagged")
## Testing temporal outliers on taxon level
## Flagged 340 records.

dat_cl[!flags, "year"]
##   [1] 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1996 1995 1995 1995
##  [26] 1995 1995 1995 1994 1994 1994 1994 1994 1994 1994 1994 1994 1994 1994 1994 1994 1993 1992 1992 1992 1992 1992 1992 1992 1992
##  [51] 1992 1992 1992 1992 1992 1991 1991 1991 1991 1991 1990 1990 1990 1990 1990 1989 1989 1989 1989 1989 1989 1989 1989 1989 1989
##  [76] 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988 1988
## [101] 1988 1988 1988 1988 1988 1988 1988 1988 1987 1987 1987 1987 1987 1986 1986 1986 1986 1986 1985 1985 1985 1985 1985 1985 1985
## [126] 1985 1985 1985 1984 1984 1984 1984 1984 1984 1984 1983 1983 1983 1983 1983 1983 1983 1983 1983 1983 1983 1983 1982 1982 1982
## [151] 1982 1981 1981 1981 1980 1980 1980 1980 1980 1980 1980 1980 1980 1980 1979 1978 1978 1974 1974 1972 1972 1972 1970 1970 1970
## [176] 1970 1970 1970 1970 1970 1969 1969 1969 1969 1969 1969 1969 1969 1969 1969 1969 1969 1968 1968 1968 1968 1967 1967 1967 1967
## [201] 1967 1967 1966 1966 1966 1966 1966 1966 1966 1966 1965 1964 1964 1963 1960 1959 1959 1959 1958 1958 1949 1949 1948 1948 1948
## [226] 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948 1948
## [251] 1948 1941 1941 1941 1941 1941 1941 1941 1941 1939 1938 1938 1938 1938 1937 1937 1937 1937 1937 1937 1937 1936 1936 1931 1931
## [276] 1931 1931 1931 1931 1930 1930 1930 1930 1930 1930 1930 1930 1930 1929 1929 1929 1929 1929 1929 1929 1929 1928 1928 1928 1928
## [301] 1928 1928 1928 1928 1928 1927 1927 1927 1927 1927 1927 1927 1927 1927 1927 1927 1927 1927 1927 1927 1927 1927 1923 1923 1922
## [326] 1920 1920 1913 1912 1912 1912 1912 1912 1912 1912 1911 1911 1911 1911 1911

dat_cl <- dat_cl[flags, ]
```

# Improving data quality using GBIF meta-data
That helped a lot, but unfortunately some unwanted records remain, especially within Europe (Fig. \ref{fig:automated}). This is mostly because we have used the occurrence records uncritically and ignored the meta-data. GBIF offers a whole lot of useful meta-data which we will use now to further refine quality of our dataset. First we'll remove coordinates with very low precision and from unsuitable data sources. We will remove all records with a precision below 100 km as this represent the grain size of our downstream analysis, but we recommend you to chose it based on your downstream analyses. We also exclude fossils as we are interested in recent distributions; and records from unknown sources, as we deem them not reliable enough.


```r
#Remove records with low coordinate precision
hist(dat_cl$coordinateUncertaintyInMeters / 1000, breaks = 20)
```

![\label{fig:automated}A histogram of the coordinate precision in the dataset..](gbif-clgbif11-1.png)

```r

dat_cl <- dat_cl %>%
  filter(coordinateUncertaintyInMeters / 1000 <= 100 | is.na(coordinateUncertaintyInMeters))

#Remove unsuitable data sources, especially fossils 
#which are responsible for the majority of problems in this case
table(dat$basisOfRecord)
## 
##     FOSSIL_SPECIMEN   HUMAN_OBSERVATION MACHINE_OBSERVATION     MATERIAL_SAMPLE         OBSERVATION  PRESERVED_SPECIMEN 
##                  10                4706                   3                  37                   1                 241 
##             UNKNOWN 
##                   2

dat_cl <- filter(dat_cl, basisOfRecord == "HUMAN_OBSERVATION" | 
                         basisOfRecord == "OBSERVATION" |
                         basisOfRecord == "PRESERVED_SPECIMEN")
```

In the next step we will remove records with suspicious individual counts. GBIF includes few records of absence (individual count = 0) and suspiciously high occurrence counts, which might indicate inappropriate data or data entry problems. 

```r
#Individual count
table(dat_cl$individualCount)
```

```
## 
##   1   2   3   4   5   6   9  11  14  15  20 
## 108  35   4   5   5   5   1   3   1   1   1
```

```r
dat_cl <- dat_cl%>%
  filter(individualCount > 0 | is.na(individualCount))%>%
  filter(individualCount < 99 | is.na(individualCount)) # high counts are not a problem
```

We might also want to exclude very old records, as they are more likely to be unreliable. For instance, records from before the second world war are often very imprecise, especially if they were geo-referenced based on political entities. Additionally old records might be likely from areas where species went extinct (for example due to land-use change).


```r
#Age of records
table(dat_cl$year)
```

```
## 
## 1997 1998 1999 2000 2001 2002 2003 2004 2005 2006 2007 2008 2009 2010 2011 2012 2013 2014 2015 2016 2017 2018 2019 2020 2021 
##    6   13   82   34   14   17   30   28   39   48   69   78   76  103  149  156  147  283  419  284  501  492  610  290  202
```

```r
dat_cl <- dat_cl%>%
  filter(year > 1945) # remove records from before second world war
```

On top of the geographic cleaning, we also want to make sure to only include species level records and records from the right taxon. The latter is not a problem in this case, as we only have one species, but it can be helpful for large datasets. Taxonomic problems such as spelling mistakes in the names or synonyms can be a severe problem. We'll not treat taxonomic cleaning here, but if you need to, check out the [taxize R package](https://docs.ropensci.org/taxize/) or the [taxonomic name resolution service](https://tnrs.biendata.org) (plants only).


```r
table(dat_cl$family) #that looks good
## 
## Felidae 
##    4170
dat_cl <- dat_cl%>%
  filter(family == 'Felidae')

table(dat_cl$taxonRank) # this is also good
## 
##    SPECIES SUBSPECIES 
##       1065       3105
```

We excluded almost 50% of the initial data points with the data cleaning, and the general picture has improved considerably. We confined the records mostly to what can be considered current day distribution of the species of interest (Fig. \ref{fig:final}).

We have, however, also lost quite a number of records. In general, there is no "one-size-fits-it-all" for data quality of geographic species occurrence records. Of course highest coordinate precision is desirable, but what is acceptable will strongly depend on the downstream analyses. For species distribution modelling, usually high precision is necessary e.g. 1-10 km, but for other analyses such as biogeographic reconstructions using tectonic plates, a record might be considered good enough quality, as long as it is on the right continent. As another example for conservation purposes it might be sufficient to know that a species is present within a certain country.

# Improving data quality using external information
Figure \ref{fig:final} shows the success of automated cleaning. However, three records within Europe remain. A short inspection of the data suggests that these are a dubious human observation and five specimens, potentially assigned to their specimen location, or fossils with misclassified meta-data. One option to automatically flag these records is to rerun the outlier test on the cleaned data. However, this would most likely also flag the isolated Indian population (which is a true presence) as problematic. 

## Flag records based on fixed longitude and latitude
The first option alternative is to exclude records outside a certain study extent. In our example this is the easiest solution because we know that lions do not occur in high latitudes any more.


```r
#exclude based on study area
dat_fin <- filter(dat_cl, decimalLatitude < 40)
```

## Flag records based on species natural ranges
In cases where simple latitudinal or longitudinal borders are not useful, an alternative is to use species ranges from external source as reference and flag all records falling outside these ranges. For amphibians, birds, mammals and reptiles the International Union for the conservation of nature (IUCN) provides detailed shape files of species' natural distribution ranges. These can be downloaded for free at https://www.iucnredlist.org/resources/spatial-data-download. *CoordinateCleaner* implements a straight forward way to use these, or any other, ranges to flag records in the `cc_iucn` function. Since downloading the IUCN shapes requires log-in we will approximate lion's natural range from scratch for our example. For plants check out the botanical countries of the [World Checklist of selected plant families](http://wcsp.science.kew.org/home.do).



```r

#create simple natural range for lions
range <- Polygon(cbind(c(-23, -7, 31, 71, 83, 42, 41, 24, -23), c(14, 37, 32, 27, 18, 0, -16, -38, 14)))
range <- Polygons(list(range), ID = c("A"))
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"
range <- SpatialPolygons(list(range), proj4string = CRS(wgs84))
df <- data.frame(species = c("Panthera leo"), row.names = c("A"))
nat_range <- SpatialPolygonsDataFrame(range, data = as.data.frame(df))

# Visualize range
plo <- fortify(nat_range)
## Regions defined for each Polygons

ggplot() +
  borders("world", colour="gray50", fill="gray50")+
  geom_polygon(data = plo, aes(x = long, y = lat, group = group))+
  theme_bw()+
  theme(legend.position = "none",
        axis.title = element_blank())
```

![plot of chunk clgbif16](gbif-clgbif16-1.png)

```r

# run cc_iucn()
range_flags <- cc_iucn(x = dat_cl,
                       range = nat_range,
                       lon = "decimalLongitude",
                       lat = "decimalLatitude",
                       value = "flagged")
## Testing natural ranges
## Warning in sp::proj4string(range): CRS object has comment, which is lost in output

## Warning in sp::proj4string(range): CRS object has comment, which is lost in output
## Flagged 139 records.

dat_fin <- dat_cl[range_flags, ]
```

![\label{fig:final}The dataset of occurrence of lions after different cleaning phases.](gbif-clgbif17-1.png)

# Identifying problematic data sets
Some types of potentially problematic coordinates can cause bias, but are not identifiable on record-level if the relevant meta-data are missing. This is especially the case if the erroneous records have been combined with precise GPS-based point occurrences into datasets of mixed precision. Two important cases are: (A) coordinate conversion errors based on the misinterpretation of the degree sign as decimal delimiter and (B) data derived from rasterized data collection designs (e.g. presence in a 50x50 km grid cell). *CoordinateCleaner* implements two algorithms to identify these problems on a dataset level.

## Identify dataset with ddmm to dd.dd conversion error
We will first run the test for erroneous data conversion due to the misinterpretation of the degree sign as decimal delimiter. We will use the `cd_ddmm` function, alternatively, you can use the `clean_dataset` wrapper. See supplementary material S1 for a detailed description of the algorithm and implementation of the test. You can control the output of the function via the `value` argument.


```r
out.ddmm <- cd_ddmm(dat_cl, lon = "decimalLongitude", lat = "decimalLatitude", 
                    ds = "species", diagnostic = T, diff = 1,
                    value = "dataset")
```

![plot of chunk clgbif18](gbif-clgbif18-1.png)

This looks good. The test indicates a slightly higher fraction of records with decimals below .60 than expected at random, but this is within the expected range and thus the test indicates no bias, which is confirmed by the diagnostic plot. In the case of a strong bias, the green points would be clustered in the bottom left quarter of the plot.

## Test for rasterized sampling
As a second step we will use the `cd_round` function to identify datasets with a significant proportion of coordinates that have been collected in large scale lattice designs. These records might have a low precision and might therefore be problematic for some analyses. For instance presence derived from a 1 degree grid of a national atlas might be to coarse for small scale species distribution models. 


```r
par(mfrow = c(2,2), mar = rep(2, 4))
out.round <- cd_round(dat_fin, lon = "decimalLongitude", 
                      lat = "decimalLatitude", 
                      ds = "species",
                      value = "dataset",
                      T1 = 7,
                      graphs = T)
## Testing for rasterized collection
```

![\label{fig:final}Diagnostic plots testing for rasterized sampling or excessive rounding. The left panel shows histograms of the record distribution, the right panel shows the autocorrelation plots. The upper panel shows longitude, the lower panel shows latitude. The logical flag in the heading of the right panel indicates the binary flag.](gbif-clgbif19-1.png)

These results look good. The dataset does not show rasterized collection schemes (see Supplementary material S1 for examples of biased datasets). The test has detected and flagged some small scale and low intensity periodicity in the longitude coordinates, however, the entire dataset is only flagged if both longitude and latitude show a pattern (as expected from rasterized sampling). You can modify the test sensitivity using various arguments. See `?cd_round` for more information.

The lion dataset is relatively small and consistent, at least in the way that it only comprises on species. For larger scale analyses you might need to deal with larger datasets, composed from a larger variety of sources.

## References
---
title: "Comparison of CoordinateCleaner to other tools"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Comparison of CoordinateCleaner to other tools}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

##Background
Erroneous database entries and problematic geographic coordinates are a central issue in biogeography and there is a set of tools available to address different dimensions of the problem. CoordinateCleaner focuses on the fast and reproducible flagging of large amounts of records, and additional functions to detect dataset-level and fossil-specific biases. In the R-environment the [scrubr](https://github.com/ropensci/scrubr) and [biogeo](https://CRAN.R-project.org/package=biogeo) offer cleaning approaches complementary to CoordinateCleaner. The scrubr package combines basic geographic cleaning (comparable to cc_dupl, cc_zero and cc_count in CoordinateCleaner) but adds options to clean taxonomic names (See also [taxize](https://github.com/ropensci/taxize)) and date information. biogeo includes some basic automated geographic cleaning (similar to `cc_val`, `cc_count` and `cc_outl`) but rather focusses on correcting suspicious coordinates on a manual basis using environmental information.


Table 1. Function by function comparison of CoordinateCleaner, scrubr and biogeo.

| Functionality                                                 | CoordinateCleaner 2.0-2                                   | scrubr 0.1.1     |biogeo 1.0| Percent overlap                      |
|--------------------------------|--------------|------------------|------------------|--------------------|
| Missing coordinates                                           | cc_val                                                    | coord_incomplete |missingvalsexclude| 100%                                 |
| Coordinates outside CRS                                       | cc_val                                                    | coord_impossible |-| 100%                                 |
| Duplicated records                                            | cc_dupl                                                   | dedup            |duplicatesexclude| The aim is identical, methods differ |
| 0/0 coordinates                                               | cc_zero                                                   | coord_unlikely   |-| 100%                                 |
| Identical lon/lat                                             | cc_equ                                                    | -                |-| 0%                                   |
| Country capitals                                              | cc_cap                                                    | -                |-| 0%                                   |
| Political unit centroids                                      | cc_cen                                                    | -  |-| 0%                                   |
| Coordinates  in-congruent with additional location information | cc_count                                                  | coord_within    |errorcheck, quickclean| 100%                                 |
| Coordinates assigned to GBIF headquarters                      | cc_gbif                                                   | -                |-| 0%                                   |
| Coordinates assigned to the location of biodiversity institutions  |    cc_inst                                     | -                |-| 0%                                   |
| Coordinates outside natural range                             | cc_iucn                                                  |-                 |-| 0%                                   |
| Spatial outliers                                              | cc_outl                                                   | -                |outliers| 50%, biogeo uses environmental distance                                   |
| Coordinates within the ocean                                  | cc_sea                                                    | -                |-| 0%                                   |
| Coordinates in urban area                                     | cc_urb                                                    | -                |-| 0%                                   |
| Coordinate conversion error                                   | dc_ddmm                                                   | -                |-| 0%                                   |
| Rounded coordinates/rasterized collection                    | dc_round                                                  | -                |precisioncheck| 20%, biogeo test for predefined rasters                                   |
| Fossils: invalid age range                                    | tc_equal                                                  | -                |-| 0%                                   |
| Fossils: excessive age range                                  | tc_range                                                  | -                |-| 0%                                   |
| Fossils: temporal outlier                                     | tc_outl                                                   | -                |-| 0%                                   |
| Fossils: PyRate interface                                     | WritePyrate                                               | -                |-| 0%                                   |
| Wrapper functions to run all test                             | CleanCoordinates, CleanCoordinatesDS, CleanCoordinatesFOS | -                |-| 0%                                   |
| Database of biodiversity institutions                          | institutions                                              | -                |-| 0%                                   |
| Taxonomic cleaning                                            | -                                                         | tax_no_epithet   |-| 0%                                   |
| Missing date                                                  | -                                                         | date_missing     |-| 0%                                   |
| Add date                                                      | -                                                         | date_create      |-| 0%                                   |
| Date format                                                   | -                                                         | date_standardize |-| 0%
| Reformatting coordinate annotation |-|-|a large set of functions|0 %|
| Correcting coordinates using guessing and environmental distance |-|-|a large set of functions|0 %|---
title: "Cleaning fossil data for the use in biogeography and palaeontology"
output: 
  rmarkdown::html_vignette:
    fig_caption: true
    number_sections: true
    citation_package: natbib
vignette: >
  %\VignetteIndexEntry{Cleaning fossil data for the use in biogeography and palaeontology}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
  %\VignetteSuggests{paleobioDB}
---



# Background
The public availability of fossils for large-scale analyses is rapidly increasing, mainly due to increased databasing efforts and data aggregators such as the paleobiology database (www.paleobiodb.org) or Neotoma (www.neotomadb.org), among others. However, data quality is an issue, in particular, for old collections or collections with uncertain taxonomy and/or bad preservation. Similar problems as known from biological collection databases (See supplementary material S2) are relevant for fossils, but in addition fossils might be dated wrongly or with very low precision. This tutorial presents a pipeline to clean fossil data from the paleobiology database (or any other) before using it in biogeographic or evolutionary analyses. We focus on identifying overly imprecisely geo-referenced and/or dated records by combining automated cleaning using *CoordinateCleaner* with cleaning based on meta-data. The proposed steps are by no means exhaustive, and keep in mind that what is "good data" depends entirely on your downstream analyses! We wrote this tutorial (and the whole *CoordinateCleaner* package) to help you to identify potential problems quicker and more reproducibly to improve data quality in large datasets. For the sake of this tutorial we will mostly remove flagged records from the dataset, however, we recommend to double check them individually.

# Install `CoordinateCleaner`
You can install the latest stable version of CoordinateCleaner (Currently 2.0-1) from CRAN using `install.packages("CoordinateCleaner")`. Alternatively you can install the latest development version from GitHub using the devtools package. We recommend the latter, to stay up-to-date. Also, make sure to have the latest R version installed. In this tutorial, relevant R-code is shown in grey boxes, the resulting output lines are marked by ##. 


```r
install.packages("devtools")
library(devtools)

install_github("ropensci/CoordinateCleaner")
```

# Load required libraries
As a first step we will load the R libraries required for the tutorial. You might need to install some of them using `install.packages`. 


```r
library(dplyr)
library(ggplot2)
library(CoordinateCleaner)
library(countrycode)
library(paleobioDB)
```

# Load test dataset
For this tutorial we will use a dataset of vascular plant fossils from the last 65 million years, downloaded from the paleobiology database using the plaeobioDB package. For the tutorial we'll limit the data to maximum 5,000 records, to keep the downloading time reasonable. If you obtained your data from the web mask of the paleobiology database, or use an entirely different database, you will have to adapt the column names in the script.


```r
#load data
dat <- paleobioDB::pbdb_occurrences(base_name = "Magnoliopsida", vocab = "pbdb", limit = 5000,
                        show = c("coords", "phylo", "attr", "loc", "time", "rem"))
rownames(dat) <- NULL
```

Alternatively, CoordinateCleaner includes an example dataset, downloaded from paleobioDB as specified above. We will use this one for the rest of this tutorial.






# Visualize the records on a map
As a first step we will visualize the records on a map, to get a general overview.


```r
#plot data to get an overview
wm <- borders("world", colour="gray50", fill="gray50")
ggplot()+ coord_fixed()+ wm +
  geom_point(data = dat, aes(x = lng, y = lat),
             colour = "darkred", size = 0.5)+
  theme_bw()
```

![plot of chunk unnamed-chunk-15](pbdb-unnamed-chunk-15-1.png)

# CoordinateCleaner
 CoordinateCleaner includes a suite of automated tests to identify problems common to biological and palaebiological databases.

## Spatial issues
We'll first check coordinate validity to check if all coordinates are numeric and part of a lat/lon coordinate reference system using `cc_val`.


```r
cl <- cc_val(dat, lat = "lat", lon = "lng")
```

```
## Testing coordinate validity
```

```
## Removed 0 records.
```

Looks good, then we will test for coordinates with equal longitude and latitude. You can use the `test` argument to specify if coordinates should be flagged if their absolute values are identical (e.g. 56,-56).


```r
cl <- cc_equ(cl, lat = "lat", lon = "lng")
```

```
## Testing equal lat/lon
```

```
## Removed 0 records.
```

For the purpose of the tutorial, we will always exclude flagged records. If you want to further explore them, and you should if by any means possible, use the `value = "flagged"` argument, valid for all functions. In that case the output value will be a vector of logical values with the same length as `dat`, where TRUE = valid record, FALSE = flagged record. It is generally advisable to check flagged records whenever possible, to avoid data-loss and false flags.


```r
fl <- cc_equ(dat, value = "flagged", lat = "lat", lon = "lng")
## Testing equal lat/lon
## Flagged 0 records.

# extract and check the flagged records
fl_rec <- dat[!fl,] 
head(fl_rec)
##  [1] occurrence_no  record_type    collection_no  taxon_name     taxon_rank     taxon_no       matched_name   matched_rank  
##  [9] matched_no     early_interval late_interval  early_age      late_age       reference_no   lng            lat           
## [17] class          class_no       phylum         phylum_no      cc             state          geogscale      early_age.1   
## [25] late_age.1     cx_int_no      early_int_no   late_int_no    genus          genus_no       family         family_no     
## [33] order          order_no       county         reid_no       
## <0 rows> (or 0-length row.names)
```

We'll also test if the records are identical, or in close vicinity to the centroids of political units. You can modify the buffer around each centroid using the `buffer` argument and the level of testing (country centroids, province centroids, or both) using the `test argument`. In case you have a list of geographic coordinates you consider problematic, for instance a list of cities you can provide them as custom gazetteer using the `ref` argument.


```r
fl <- cc_cen(cl, lat = "lat", lon = "lng", value = "flagged")
## Testing country centroids
## Flagged 6 records.
fl_rec <- cl[!fl, ]
unique(fl_rec$cc)
## [1] JP
## Levels: NZ US CN IR KP AU UK FR JP DE CA RU KE ZM EG CD ZA TZ UG ET MX IT
cl <- cl[fl, ]
```

Next we will test if the coordinates are within the country they are assigned to. This test is a bit more tricky, as it will also flag records, if the country name in the country column is not following ISO3 or if the records have been assigned during a different political landscape. For instance records from former Western and Eastern Germany. Here we need to convert the country annotation in column cc from ISO2 to ISO3; it is advisable to double check which records have be flagged, to avoid unnecessary data loss (see above).


```r
#adapt country code to ISO3, for country test
cs_ma <- "GBR"
names(cs_ma) <- "UK"
cl$cc_iso3 <- countrycode(cl$cc, origin = "iso2c", destination = "iso3c", custom_match = cs_ma)

cl <- cc_coun(cl, lat = "lat", lon = "lng", iso3 = "cc_iso3")
## Testing country identity
## Removed 234 records.
```

Next we will test if any of the records bear the coordinates of a hosting biodiversity institution or the GBIF headquarters, using the `institutions` database of *CoordinateCleaner*. As for the country centroid test you can change the buffer around the institutions with the `buffer` argument.


```r
cl <- cc_inst(cl, lat = "lat", lon = "lng")
## Testing biodiversity institutions
## Removed 0 records.
cl <- cc_gbif(cl, lat = "lat", lon = "lng")
## Testing GBIF headquarters, flagging records around Copenhagen
## Removed 0 records.
```

Finally, we will test for plain zero coordinates (e.g. 0/0).


```r
cl <- cc_zero(cl, lat = "lat", lon = "lng")
## Testing zero coordinates
## Removed 0 records.
```

## Temporal issues
The spatial cleaning above is mostly identical with steps from recent geographic records. Additionally \emph{CoordinateCleaner} includes three functions to test the temporal dimension of fossils. Fossil ages are usually defined with a maximum and a minimum range, based on geological strata. First we will exclude records without dating information (NA) and then test for records with equal minimum and maximum range. Unless your data includes absolutely dated fossils, this will most likely be an data entry error.



```r
cl <- cl[!is.na(cl$late_age),]
cl <- cl[!is.na(cl$early_age),]
cl <- cf_equal(cl, min_age = "late_age", max_age = "early_age")
## Testing age validity
## Removed 0 records.
```

Next we will look at the age range (= max age - min age) of each record. The age range is the dating precision and can vary considerably, depending on the data available for dating. For many analyses, for instance in PyRate, very imprecisely dated records are not suitable. Lets first have a look at the age ranges in our test dataset.


```r
rang <- cl$early_age - cl$late_age
hist(rang, breaks = 40, xlab = "Date range [max age - min age]", main = "")
```

![plot of chunk unnamed-chunk-24](pbdb-unnamed-chunk-24-1.png)

Some individual records are dated with a precision of more than 60 million years! \emph{CoordinateCleaner} offers two ways to flag records based on their age range (1) based on absolute age, e.g. age range > 35 million years or (2) based on age range outlier detection in the entire dataset (e.g. if few records are much less precisely dated than the rest of all records) and (3) based on age range outlier detection on taxon level (e.g. all \emph{Quercus} records that are much less precisely dated than the other \emph{Quercus} records. The second and third approach can be combined and offer some more flexibility over the absolute age limit, but need some consideration on the desired sensitivity. Here, we will run all three variants for illustration, if you use your own data you should decide which one is more suitable depending on your downstream analyses. In the case of (2) and (3) you can tweak the test sensitivity using the `mltpl` argument.


```r
# Outlier dataset
cl <- cf_range(cl, taxon = "", min_age = "late_age", max_age = "early_age")
## Testing temporal range outliers on dataset level
## Removed 57 records.

# Outlier per taxon
cl <- cf_range(cl, taxon = "taxon_name", min_age = "late_age", max_age = "early_age")
## Testing temporal range outliers on taxon level
## Removed 86 records.

# Absolute age limit
cl <- cf_range(cl, taxon = "taxon_name", min_age = "late_age", 
               max_age = "early_age", method = "time", max_range = 35)
## Testing temporal range outliers on taxon level
## Removed 1 records.

rang <- cl$early_age - cl$late_age
hist(rang, breaks = 40, xlab = "Date range [max age - min age]", main = "")
```

![plot of chunk unnamed-chunk-25](pbdb-unnamed-chunk-25-1.png)

Finally we will test for outliers in space-time, that is records that are either very distant in space or in time from all other records (1) in the dataset (2) per taxon. The test is again based on quantile outlier detection and can be modified using  various arguments. Here it is important to carefully consider the desired test sensitivity. See `?cf_outl` for help.


```r
# Outlier dataset
cl <- cf_outl(cl, taxon = "", lat = "lat", lon = "lng",
              min_age = "late_age", max_age = "early_age")
## Testing spatio-temporal outliers on dataset level
## Removed 256 records.

# Outlier taxon
cl <- cf_outl(cl, taxon = "taxon_name", lat = "lat", lon = "lng",
              min_age = "late_age", max_age = "early_age")
## Testing spatio-temporal outliers on taxon level
## Removed 30 records.
```

Done! To check how many records have been flagged in total, you can compare the two datasets.


```r
nrow(dat) - nrow(cl)
```

```
## [1] 689
```


All test combined have removed about 689 records (13.8%). If you want to identify all flagged records, to double check or correct them, take a look at the `clean_fossils` wrapper function below. 

So far so good, we have significantly refined the data  for our needs. In section 6 we will have a look at the meta-data for further refinement, but before, note that there are two different ways to run *CoordinateCleaner*. You can  connect all functions directly in a row using the magrittr pipe (%>%) operator.


```r
cl <- dat%>%
  cc_val(lat = "lat", lon = "lng")%>%
  cc_equ(lat = "lat", lon = "lng")%>%
  cc_cen(lat = "lat", lon = "lng")%>%
  cc_coun(lat = "lat", lon = "lng", iso3 = "cc")%>%
  cc_gbif(lat = "lat", lon = "lng")%>%
  cc_inst(lat = "lat", lon = "lng")%>%
  cc_zero(lat = "lat", lon = "lng")%>%
  cf_equal(min_age = "late_age", max_age = "early_age")%>%
  cf_range(taxon = "taxon_name",
           min_age = "late_age", max_age = "early_age")%>%
  cf_outl(taxon = "taxon_name", 
          lat = "lat", lon = "lng",
          min_age = "late_age", max_age = "early_age")
```

Alternatively you can use `clean_fossils`, a wrapper around all quality tests provided by *CoordinateCleaner* relevant for fossil data. See `?CleanCoordiantesFOS` for help.


```r
#adapt country code to ISO3, for country test
cs_ma <- "GBR"
names(cs_ma) <- "UK"
dat$cc <- countrycode(dat$cc, origin = "iso2c", destination = "iso3c", custom_match = cs_ma)

#run automated testing
flags <- clean_fossils(x = dat, 
                             taxon = "taxon_name",
                             min_age = "late_age", max_age = "early_age", 
                             value = "spatialvalid")

head(flags)
cl <- dat[flags$.summary,] #the cleaned records
fl_rec <- dat[!flags$.summary,] # the flagged records for verification
```

# Improving data quality using meta-data
Usually, at least some type of meta-data are provided with fossil occurrences, as is the case in the paleobiology database. We'll now explore these and see if we can identify further problems.  

## Basic taxonomy
First we'll take a short look at taxonomy. Fossil taxonomy is very complex and composite databases often have taxonomic issues that are extremely difficult to resolve. Here we will only do some very basic checks to test if: 
1. all taxa in our dataset are plants, 2. they are at least identified to genus level.


```r
#1. This looks OK
table(cl$phylum)
## 
## Spermatophyta 
##          4311

#2. Taxonomic level of identification
table(cl$taxon_rank)
## 
##  subclass   species     genus     class    family     order subfamily 
##         7      3340       536       271       124        20         0
```

The required taxonomic level of course depends on the downstream analyses, but here we will exclude everything other than genus or species, which is a reasonable approach for most PyRate analyses.


```r
cl <- cl %>%
  filter(taxon_rank %in% c("species", "genus"))
```

##  Spatial coordinates
The Paleobiology database includes some information on the basis of the geographic data for many records.


```r
table(cl$geogscale)
```

```
## 
## small collection          outcrop       local area            basin 
##             1953             1750               79                0
```

As expected most records are only roughly geo-referenced, but the precision is still relatively high for many records. 

## Time 
We have checked for potentially problematic records in time and space above, but it is definitively advisable to check again.


```r
#minimum ages
tail(table(cl$late_age))
## 
##  63.3    66  70.6  93.5  93.9 100.5 
##    53   138     8     3     3    21

ggplot(cl)+
  geom_histogram(aes(x = late_age))
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-33](pbdb-unnamed-chunk-33-1.png)

```r

#maximum ages
tail(table(cl$early_age))
## 
##  70.6  83.5  99.6 100.5 105.3   113 
##   126     9     3    11     3    21

ggplot(cl)+
  geom_histogram(aes(x = early_age))
## `stat_bin()` using `bins = 30`. Pick better value with `binwidth`.
```

![plot of chunk unnamed-chunk-33](pbdb-unnamed-chunk-33-2.png)

The minimum and maximum ages look unproblematic, but there are still some records with very large temporal uncertainties, and at least one case where the minimum and maximum age seem reversed. This might be informative in some cases, but for most analysis this might be problematic, so here we exclude all records with temporal uncertainty above 20.442 million years, which will retain 95% of the data. This is an arbitrary choice, and you'll have to choose a more suitable value based on your planned analyses.

## Additional meta-data
If you download data using the paleobiology database web portal, there are some additional fields you can check to generally gauge data precision, namely basis and precision of the geo-reference, the year of publications (old publications are more likely to be imprecise), preservation quality (bad preservation increases taxonomic uncertainty) and plant organ (also taxonomic uncertainty) and collection type. You might need to adapt the column names, and not all columns will be relevant for all taxa (e.g. "plant organ").


```r

table(cl$latlng_basis)
table(cl$latlng_precision)


table(cl$ref_pubyr)
table(cl$collection_type)
table(cl$preservation_quality)
table(cl$plant_organ)

cl <- filter(cl, preservation_quality != "very poor")
```

# Conclusions

Through the various cleaning steps outline above, we have identified some potential major caveats and hopefully increased the quality of the dataset. We have excluded a significant fraction of all records (-76.52 %). Data quality is a delicate issue, especially for fossils from compound data bases and the usefulness of individual records will depend on your downstream analyses. We hope that you find this tutorial useful in exploring data downloaded from the Paleobiology database and to explore the quality of any fossil dataset. 


![plot of chunk unnamed-chunk-35](pbdb-unnamed-chunk-35-1.png)![plot of chunk unnamed-chunk-35](pbdb-unnamed-chunk-35-2.png)

#Writing the result to disk in PyRate format

If you want to use fossil data to estimate speciation and extinction rates, their correlation with environmental factors or to do biogeography, [PyRate](https://github.com/dsilvestro/PyRate) might be useful for you. You can use the `write_pyrate` function of *CoordinateCleaner* to write ready-to-use PyRate input files from your now cleaned dataset to disk. To do so, you additionally need an assessment of which species are extinct and which are extant today, which you can provide via the `status` argument of `write_pyrate`. If you want to specify a path were to save the file, use the `path` argument instead of `fname`.


```r
# replace  blanks in taxon names
cl$taxon_name <- gsub("[[:blank:]]{1,}","_", cl$taxon_name)

#simulated current status, soley for demonstration purposes, replace with your own data
mock_status <- data.frame(taxon_name = unique(cl$taxon_name),
                          status = sample(c("extinct", "extant"), 
                                          size = length(unique(cl$taxon_name)), 
                                          replace = TRUE))

#add current status to fossils
cl2 <- inner_join(cl, mock_status, by = "taxon_name")

#Write PyRate input to disk
write_pyrate(cl, fname = "paleobioDB_angiosperms", status = cl2$status,
            taxon = "taxon_name", min_age = "late_age", max_age = "early_age")
```
---
title: "Using customized gazetteers"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using customized gazetteers}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteSuggests{rgbif}
  \usepackage[utf8]{inputenc}
---




CoordinateCleaner identifies potentially erroneous geographic records with coordinates assigned to the sea, countr coordinate, country capitals, urban areas, institutions, the GBIF headquarters and countries based on the comparison with geographic gazetteers (i.e. reference databases). All of these functions include default reference databases compiled from various sources. These default references have been selected suitable for regional to global analyses. They will also work for smaller scale analyses, but in some case different references might be desirable and available. this could be for instance centroids of small scale political units, a different set of urban areas, or a different coastline when working with coastal species. To account for this, each *CoordinateCleaner* function using a gazetteer has a `ref` argument to specify custom gazetteers.

We will use the case of coastlines and a coastal species to demonstrate the application of custom gazetteers. The purpose of `cc_sea` is to flag records in the sea, since these often represent erroneous and undesired records for terrestrial organisms. The standard gazetteer for this function is fetched from naturalearthdata.com at a 1:50m scale. However, often coordinates available from public databases are only precise at the scale of kilometres, which might lead to an overly critical flagging of coordinates close to the coastline, which is a problem especially for coastal or intertidal species. WE illustrate the issue on for the mangrove tree genus *Avicennia*.


```r
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(rgbif)
library(sp)
library(viridis)


#download data from GBIF
dat <- rgbif::occ_search(scientificName = "Avicennia", limit = 1000,
         hasCoordinate = T)

dat <- dat$data

dat <-  dat %>% 
  dplyr::select(species = name, decimallongitude = decimalLongitude, 
         decimallatitude = decimalLatitude, countryCode)

# run with default gazetteer
outl <- cc_sea(dat, value = "flagged")
## OGR data source with driver: ESRI Shapefile 
## Source: "C:\Users\az64mycy\AppData\Local\Temp\Rtmp4SRhHV", layer: "ne_110m_land"
## with 127 features
## It has 3 fields

plo <- data.frame(dat, outlier =  as.factor(!outl))

#plot results
ggplot()+
  borders(fill = "grey60")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")
```

![plot of chunk cusgaz1](cusgaz-cusgaz1-1.png)

A large number of the coastal records gets flagged, which in this case is undesirable, because it is not a function of the records being wrong, but rather of the precision of the coordinates and the resolution of the reference. To avoid this problem you can use a buffered reference, which avoids flagging records close to the coast line and only flags records from the open ocean. *CoordinateCleaner* comes with a one degree buffered reference (`buffland`). In case a narrower or distance true buffer is necessary, you can provide any SpatialPolygonsDataFrame similar in structure to `buffland` via the `ref` argument.


```r
# The buffered custom gazetteer
data("buffland")
plot(buffland)
```

![plot of chunk cusgaz2](cusgaz-cusgaz2-1.png)

```r

# run with custom gazetteer
outl <- cc_sea(dat, value = "flagged", ref = buffland)

plo <- data.frame(dat, outlier =  as.factor(!outl))

#plot results
ggplot()+
  borders(fill = "grey60")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")
```

![plot of chunk cusgaz2](cusgaz-cusgaz2-2.png)

---
title: "A global gazetteer of biodiversity institutions"
output: rmarkdown::html_vignette
bibliography: CoordinateCleaner.bib
vignette: >
  %\VignetteIndexEntry{A global gazetteer of biodiversity institutions}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteSuggests{caret}
  %\VignetteSuggests{viridis}
  %\VignetteSuggests{countrycode}
  %\VignetteSuggests{tidyverse}
---

```{r setup, include=FALSE, message = F, warning = FALSE}
library(caret)
library(CoordinateCleaner)
library(countrycode)
library(dplyr)
library(ggplot2)
library(magrittr)
library(raster)
library(viridis)
```

```{r, echo = FALSE}
knitr::opts_chunk$set(
   fig.height=8, 
   fig.width=8
)
```


```{r, echo = F, message = FALSE, warning = FALSE}
data(institutions)
institutions <- filter(institutions, !is.na(decimallongitude))
institutions$source <- trimws(institutions$source)
```

## Background
Most of the geographic species occurrence records publicly available from aggregated databases such as the Global Biodiversity Information Facility (GBIF), are either based on collected specimens stored in a museum, university, botanical garden, herbarium or zoo, or on human observations, e.g. vegetation surveys or citizen science projects. A relatively common error in the geographic information of these records are coordinates assigned to the physical location of the institution hosting the specimen. The reasons for these errors may include among others individuals escaped from horticulture, specimens erroneously geo-reference to their physical location as well as records based on pictures taken by laymen in zoos or botanical gardens. These records are problematic as the conditions at these locations do not represent the species' natural habitat and might in fact differ considerably from them.

To identify these records, CoordinateCleaner includes a novel geo-referenced global database of biodiversity institutions - defined here as institutions that generally are concerned with biodiversity research and/or hosting collections of living or mounted biological specimens. We implement a cleaning check using this database as gazetteer in the `cc_inst` function and the `institutions` argument of the `clean_coordinates` function of the *CoordinateCleaner* R-package. Furthermore, we hope that this database can prove useful beyond cleaning geographic records, for instance to assess sampling biases in biological collections.


## Data compilation
We compiled names of biodiversity institutions from six different sources [@BGCI-BotanicGardensConservationInternational2017; @IndexHerbariorum2017;@TheGlobalRegistryofBiodiversityRepositories2017; @Wikipedia2017; @GlobalBiodiveristyInformationFacility2017; @GeoNames2017] and geo-referenced them using the Google maps API via the ggmap package in R [@Kahle2013] using institution names and, if this yielded no results the institutions address. For those records that did not yield any results we used opencage via the opencage R-package [@Salmon2017] for geo-referencing. We manually geo-referenced those institutions that could not be geo-referenced automatically (c. 50%) using the WWW and Google Earth [@GoogleInc2017]. In total the database comprises almost 9700 geo-referenced institutions (and another 2500 entries for which geo-referencing was not possible, either to problems with non-English names or geographic ambiguities). The spatial extent of the database is global, but we acknowledge that there is a focus on English-speaking countries and countries using the Roman alphabet. This is partly a bias due to the data compilation process. We hope that this bias can be overcome by future contributions to the database from researchers in non-English speaking and non-Roman alphabet countries. In general, we acknowledge that the database may not be complete and created a webmask at (http://biodiversity-institutions.surge.sh/) were researchers can easily submit their institution or a comment on an existing institution. The webpage also includes an overview on the institutions included in the dataset.


```{r fig1, echo = F, evaluate = T, warning = F, fig.show = T}
plo <- institutions

plo$source <- factor(plo$source, levels = names(sort(table(plo$source))))

ggplot(data = plo)+
  geom_bar(aes(x = source))+
  xlab("Source")+
  ylab("Count")+
  theme_bw()
```



## Data structure
In addition to  the name and geographic coordinates for each institution, the database includes information on the type of the institutions ("type", e.g. "herbarium" or "university"), the source from where we obtained the name of the institution ("source"), the precision of the coordinates ("geocoding.precision.m" and "geocoding.issue") as well as the city and address (when available, "city" and "address"). The quality of the meta-data might vary among different sources). Furthermore, the database includes a column identifying if the respective institution is located within a protected area [@UNEP-WCMCandIUCN2017], and if this is the case, the World Database of Protected Areas ID of the respective protected area (WDPA, shape file available at: https://www.protectedplanet.net/). We include this flag, as biodiversity institutions within protected areas might or might not be relevant for coordinate cleaning, depending on downstream analyses.




```{r fig2, echo = F, evaluate = T, warning = F, fig.show = T}
#number per 100km x 100 km grid cell
#reference raster
ras <- raster::raster("C:/Users/az64mycy/Dropbox (iDiv)/research_projects/00_CoordinateCleaner/CoordinateCleaner/articles/inst/ANNUAL_NDVI.tif")

#projections
wgs1984 <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs")
behr <- CRS('+proj=cea +lon_0=0 +lat_ts=30 +x_0=0 +y_0=0 +datum=WGS84 +ellps=WGS84 +units=m +no_defs')

#select and reproject institutions
abu <- institutions%>%
  dplyr::select(species = type, decimallongitude, decimallatitude) %>% 
  filter(!is.na(decimallongitude))

abu.b <- abu%>%
  dplyr::select(decimallongitude, decimallatitude)%>%
  sp::SpatialPoints(proj4string = wgs1984)%>%
  spTransform(behr)

abu <- abu.b %>% 
  rasterize(ras, fun ="count") %>% 
  raster::rasterToPoints()%>%
  data.frame()

ggplot()+
  geom_raster(data = abu, aes(x = x, y = y, fill = log(layer)))+
  scale_fill_gradient2(low = "blue", mid = "yellow", high = "red", midpoint = 2, name = "Number of\ninstitutions\n[log]")+
  xlab("Longitude")+
  ylab("Latitude")+
  theme_bw()+
  coord_fixed()+
  theme(axis.text = element_blank(),
        axis.ticks = element_blank(),
        strip.background=element_blank(),
        legend.position = c(0.09, 0.35))

```



```{r fig3, echo = F, evaluate = T, warning = F, fig.show = T}
#Institutions per continent
cont <- institutions%>%
  mutate(continent = countrycode(country, origin = "iso3c", destination = "continent"))


sor <- cont%>%
  group_by(continent)%>%
  summarize(numb = n())%>%
  arrange(numb)


plo <- cont%>%
  mutate(continent = factor(cont$continent, levels = sor$continent))%>%
  filter(!is.na(continent))%>%
  filter(!is.na(type))
  

ggplot(data = plo)+
  geom_bar(aes(x = continent))+
  theme(axis.text= element_text(angle = 90, size = 5))+
  facet_wrap(~type, ncol = 2)+
  ylab("Count")+
  xlab("Institution type")+
  theme_bw()

```



```{r fig4, echo = F, evaluate = T, warning = F, fig.show = T}
#institutions per country
sor <- institutions%>%
  group_by(country)%>%
  summarize(numb = n())%>%
  arrange(numb)

sor <- sor[(nrow(sor) - 10):nrow(sor),]

sor2 <- sor %>%
  mutate(country = countrycode::countrycode(country, origin = "iso3c", destination = "country.name.en"))

plo <- institutions%>%
  filter(!is.na(country))%>%
  filter(country %in% sor$country)%>%
  mutate(country = countrycode::countrycode(country, origin = "iso3c", destination = "country.name.en"))


plo <- plo%>%
  mutate(country = factor(plo$country, levels = sor2$country))


plo$country <- gsub("United Kingdom of Great Britain and Northern Ireland", "UK", plo$country)
plo$country <- gsub("United States of America", "USA", plo$country)

plo <- plo%>%
  mutate(country = factor(plo$country, levels = c(sor2$country[1:8], "UK", "USA")))

ggplot(data = plo)+
  geom_bar(aes(x = country))+
  theme(axis.text= element_text(angle = 90, size = 5))+
  ylab("Count")+
  theme_bw()+
  theme(axis.title.x = element_blank())
```

## Data accessability
The database is open-source and available as R data file (.rda) as part of the *CoordinateCleaner* package either from [CRAN](https://cran.r-project.org/web/packages/CoordinateCleaner/index.html) or [GitHub](https://github.com/ropensci/CoordinateCleaner) under a CC-BY license. We acknowledge, that this database is not complete and can constantly be improved, any feedback can be provided via the GitHub page of \emph{CoordinateCleaner} (https://github.com/ropensci/CoordinateCleaner/). 

# References
---
title: "Dataset-level cleaning"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Dataset-level cleaning}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteSuggests{caret}
  %\VignetteSuggests{readr}
---

```{r setup, include=FALSE, message = F}
knitr::opts_chunk$set(echo = TRUE)
library(caret)
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(knitr)
library(magrittr)
library(raster)
library(readr)
library(tidyr)
library(viridis)

```


```{r, echo = FALSE}
# A function for generating captions and cross-references

fig <- local({
    i <- 0
    list(
        cap=function(refName, text, center=FALSE, col="black", inline=FALSE) {
            i <<- i + 1
            ref[[refName]] <<- i
            css_ctr <- ""
            if (center) css_ctr <- "text-align:center; display:inline-block; width:100%;"
            cap_txt <- paste0("<span style=\"color:", col, "; ", css_ctr, "\">Figure ", i, ": ", text , "</span>")
            anchor <- paste0("<a name=\"", refName, "\"></a>")
            if (inline) {
                paste0(anchor, cap_txt)    
            } else {
                list(anchor=anchor, cap_txt=cap_txt)
            }
        },
        
        ref=function(refName, link=FALSE, checkRef=TRUE) {
            
            ## This function puts in a cross reference to a caption. You refer to the
            ## caption with the refName that was passed to fig$cap() (not the code chunk name).
            ## The cross reference can be hyperlinked.
            
            if (checkRef && !refName %in% names(ref)) stop(paste0("fig$ref() error: ", refName, " not found"))
            if (link) {
                paste0("<A HREF=\"#", refName, "\">Figure ", ref[[refName]], "</A>")
            } else {
                paste0("Figure ", ref[[refName]])
            }
        },
        
        ref_all=function(){
            ## For debugging
            ref
        })
})

## This chunk replaces the default hook for processing plots. It achieves the purposes,
## of laying out auto-numbered captions, but other functionality may be gone.


knit_hooks$set(plot = function(x, options) {
    sty <- ""
    if (options$fig.align == 'default') {
        sty <- ""
    } else {
        sty <- paste0(" style=\"text-align:", options$fig.align, ";\"")
    }
    
    if (is.list(options$fig.cap)) {
        ## options$fig.cap is a list returned by the function fig$cap()
        str_caption <- options$fig.cap$cap_txt
        str_anchr <- options$fig.cap$anchor
    } else {
        ## options$fig.cap is a character object (hard coded, no anchor)
        str_caption <- options$fig.cap
        str_anchr <- ""
    }
    
    paste('<figure', sty, '>', str_anchr, '<img src="',
        opts_knit$get('base.url'), paste(x, collapse = '.'),
        '"><figcaption>', str_caption, '</figcaption></figure>',
        sep = '')
    
})

## This chucnk will read through *this* Rmd file, and attempt to extract all of the 
## labels (not caption text) used for Figure captions. These labels are used
## as anchors, so scanning through the document now will allow us to create cross references
## before the caption actually appears. 

## Get the name of this Rmd file
rmdFn <- knitr::current_input()  # filename of input document

## Read lines and close connection
rmdCon <- file(rmdFn, open = "r")
rmdLines <- readLines(rmdCon)
close(rmdCon)

## Pull out all occurences of at least one back tick, followed 
## by any number of characters, followed by fig$cap (all on one line)
figscap_idx <- grep("`+(.*)fig\\$cap", rmdLines)
rmdLines <- rmdLines[figscap_idx]

## Get rid of everything up until the start of the caption label
## This presumes the caption label is the first argument of fig$cap()
## E.g., fig.cap = fig$cap("my_label", ...)
rmdLinesSansPre <- sub("(.*)fig\\$cap(.*?)[\"']", "", rmdLines)

## Identify everything up until the first quote
match_data <- regexpr("(.*?)[\"']", rmdLinesSansPre)

## Reduce the length by one, because we're not interested in the final quote
attr(match_data, "match.length") <- attr(match_data, "match.length") - 1

## Extract
fig_labels <- regmatches(rmdLinesSansPre, match_data, invert=FALSE)

if (length(fig_labels) > 0) {

    ## Test for duplicates
    if (anyDuplicated(fig_labels) > 0) stop("Duplicate caption labels detected")
    
    ## Create a named list of Figure numbers
    ref <- as.list(1:length(fig_labels))
    names(ref) <- fig_labels
}    


```


## Background
Some problems with biological collection data are not apparent from individual records, but rather linked to properties of an entire data set. The CleanCoordinatesDS function can use dataset properties to flag three types of potential problems, under the assumption that these problems will effect many records in a dataset from the same source (but not necessarily all): 

1. An erroneous conversion of coordinates in degree minute annotation into decimal degrees, where the decimal sign is erroneously translated into the decimal delimiter, e.g. 10°30' to 10.3 °. This is a problem that has been observed in particular for older data sets.

2. A periodicity in the decimals of a data set, as will arise when coordinates are either rounded or recorded in a raster format and coordinates represent raster cell centres. This problem represents low precision rather than errors, but can also be fatal, if undetected and taken as actual localities, for example in distribution modelling.

# 1. Conversion errors between degree-minute and decimal-degree annotation (`cd_ddmm`)
## Background
Geographic coordinates in a longitude/latitude coordinate reference system can be noted in different ways. The most common notations are degree minutes seconds (ddmm, e.g. 38°54'22'') and decimal degrees (dd.dd, e.g. 38.90611°). A hybrid annotation of degrees with decimal minutes is also sometimes used (e.g. 38°54.367). However, analyses using distribution data almost exclusively rely on the machine readable decimal-degree format. Therefore the diversity of formats is challenging for databases comprised of coordinate records from different sources which potentially use different annotation formats. Systematic errors can arise, especially if old data are combined and digitized automatically, without appropriate conversion. A particular problem reported repeatedly is the misinterpretation of the degree sign (°) as the decimal delimiter (e.g. 38°54' converted to 38.54°), which leads to biased geographic occurrence information. A particular caveat in identifying these problems post-hoc in a database compiled from many sources is that biased records might be mixed with unproblematic records. For instance, specimens from a certain herbarium might have been digitized in different instances in a way, that part of the records are biased whereas others are not.

## Algorithm and implementation
As part of this study we present a novel algorithm to identify data sets potentially biased by erroneous conversion from ddmm to dd.dd due to the misinterpretation of the degree sign as decimal delimiter. 

```{r fig5, echo = F, evaluate = T, warning = F, fig.show = T, fig.height=8, fig.width=8, fig.cap = fig$cap("fig_analmat", "Examples of the analysis matrices for data sets with varying fractions of ddmm to dd.dd conversion errors. Error fraction is the percentage of records with conversion error in the data set.")}
load(file = "inst/analyses_matrix_Clean")
dat0 <- raster(dat.t1)%>%
  rasterToPoints()%>%
  as.data.frame()%>%
  mutate(layer = ifelse(layer == 0, NA, layer))%>%
  mutate(layer = as.character(layer)) %>% 
  mutate(layer = parse_factor(layer, levels = c(NA,"1")))%>%
  mutate(ds = "Error fraction 0")

load(file = "inst/analyses_matrix_Error_fraction_0.2")
dat02 <- raster(dat.t1)%>%
  rasterToPoints()%>%
  as.data.frame()%>%
  mutate(layer = ifelse(layer == 0, NA, layer))%>%
  mutate(layer = as.character(layer)) %>% 
  mutate(layer = parse_factor(layer, levels = c(NA,"1")))%>%
  mutate(ds = "Error fraction 0.2")

load("inst/analyses_matrix_Error_fraction_0.5")
dat08 <- raster(dat.t1)%>%
  rasterToPoints()%>%
  as.data.frame()%>%
  mutate(layer = ifelse(layer == 0, NA, layer))%>%
  mutate(layer = as.character(layer)) %>% 
  mutate(layer = parse_factor(layer, levels = c(NA,"1")))%>%
  mutate(ds = "Error fraction 0.5")

load("inst/analyses_matrix_Error_fraction_1")
dat1 <- raster(dat.t1)%>%
  rasterToPoints()%>%
  as.data.frame()%>%
  mutate(layer = ifelse(layer == 0, NA, layer))%>%
  mutate(layer = as.character(layer)) %>% 
  mutate(layer = parse_factor(layer, levels = c(NA,"1")))%>%
  mutate(ds = "Error fraction 1")

dat <- bind_rows(dat0, dat02, dat08, dat1)
dat$ds <- factor(dat$ds, levels = c("Error fraction 0", "Error fraction 0.2", "Error fraction 0.5", "Error fraction 1"))

ggplot(data = dat)+
  geom_raster(aes(x = x, y = y, fill = layer))+
  scale_fill_manual(values = c("darkgreen", "white"))+
  scale_y_continuous(limits=c(0, 1), expand = c(0, 0))+
  scale_x_continuous(limits=c(0, 1), expand = c(0, 0))+
  xlab ("Longitude decimals")+
  ylab("Latitude decimals")+
  theme_bw()+
  theme(legend.position = "none")+
  facet_wrap(~ds)
```

The binomial test and frequency comparison described in the main text are based on an analysis matrix which is recording the distribution of coordinates decimals in the longitude/latitude space (`r fig$ref("fig_analmat", link = TRUE)`). The analysis matrix does not record the number of records in a cell, but only presence/absence, to account for clustered sampling (i.e. a large number of records with similar decimals are most likely not related to conversion error, but rather to multiple samples from the same location, or coordinate rounding, see section 3). 

The test is implemented in the `cd_ddmm` function and the `ddmm` argument of the `clean_dataset` wrapper function. The input follows the general package design, and is a `data.frame` with at least three columns including the decimal longitude and latitude and a data set identifier for each record. The names of the columns can be specified via the `lon`, `lat` and `ds` arguments. There are three additional arguments to customize the test, in particular to modify its sensitivity to the fraction of a data set that is biased: The `pvalue` argument controls the cut-off p-value for significance of the one-sided t-test. The `diff` argument controls the threshold difference for the `ddmm test`, and indicates by which fraction the records with decimals below 0.6 must outnumber the records with decimals above 0.6. The size of the analysis matrix can be adapted using the `mat.size` argument.

### Simulations
We used simulations to asses (A) the effect of varying the diff parameter and (B) the effect of data set size on the performance of the `ddmm` test.

We simulated 100,000 datasets of species occurrences with varying number of records and degree of sample clustering. For each iteration we first draw a random number, N as $\Gamma(\alpha = 2, \beta = 1) *500$ for the number of records. We the simulated N latitude and longitude coordinates between 0° an 90° using $K \in [1,5)$ truncated normal distributions with $\mu_i \sim \mathcal{U}(0,90)$ and $\sigma \sim \mathcal{U}(0.1,5)$. We then added a bias by replacing between 0 - 80% of the samples by records with decimals sampled from $\mathcal{U}(0, 0.599)$. We then analysed the simulated data using `cd_ddmm` using diff parameters between 0.1 and 1.

`r fig$ref("fig_ddmm_diff", link = TRUE)` shows the effect of the `diff` parameter on the sensitivity of the `cd_ddmm` test and `r fig$ref("fig_dataset_size", link = TRUE)` the effect of dataset size. In general, the test is identifying the presence of a bias with rate > 0.1 (10\% bias) well, for datasets with more than 100 individual occurrence records. However, for empirical data we recommend a more conservative `diff` threshold of 0.5-1 to identify datasets with more than 30\% bias. This is because smaller bias rates might be caused by irregular sampling rather than conversion errors. However, the advisable diff threshold depends on downstream analyses and higher `diff` values might be necessary (`r fig$ref("fig_ddmm_diff", link = TRUE)`). We suggest to manually check the decimal distribution in flagged datasets to avoid data loss. This can be done easily by visually inspecting the analyses matrix, by setting the `diagnostic` argument of `cd_ddmm` to TRUE. 

```{r fig7, echo = F, evaluate = T, warning = F, fig.show = T, message = F, fig.height=8, fig.width=8, fig.cap = fig$cap("fig_ddmm_diff", "The effect of the `diff` parameter on the `cd_ddmm` test. The sensitivity of the test decreases with higher `diff` values, meaning that a higher  fraction of records with erroneous coordinate is tolerated. The error fraction in data set is the fraction of records subject to conversion error.")}
# False positive/negative rates - diff parameter
load("inst/simulations_cd_ddmm.Rds")
dat <- plo <- simulations_cd_ddmm

#plot summary stats with non-bias = 0
plo <- plo[!is.na(plo$pass),]
plo <- plo%>%
  mutate(biased = ifelse(bias.rate == 0, F, T))%>%
  mutate(correct = pass != biased)
#remember: pass = T means not biased/test passed

summ <- plo%>%
  group_by(diff, bias.rate)%>%
  summarize(true.positive = sum(biased & !pass),
            true.negative = sum(!biased & pass),
            false.positive = sum(!biased & !pass),
            false.negative = sum(biased & pass),
            sensitivity = true.positive / (true.positive + false.negative),
            specificity = true.negative / (true.negative + false.negative),
            precision = true.positive /(true.positive + false.positive),
            detection.rate = true.positive / (true.positive + false.positive + true.negative + false.negative),
            n = n(),
            false.negative.rate = false.negative / n, 
            true.positive.rate = true.positive / n,
            true.negative.rate = true.negative / n,
            false.positive.rate = false.positive / n)

suma <- summ%>%
  filter(bias.rate == 0)%>%
  dplyr::select(diff, false.positives = false.positive, n, false.positive.rate)%>%
  mutate(false.positive.rate = round(false.positive.rate, 2))

summ <- summ%>%
  dplyr::select(-true.positive, -true.negative, -false.positive, -false.negative, -n, -precision, 
                -specificity, -detection.rate, -sensitivity)%>%
  mutate(false.positive.rate = ifelse(bias.rate > 0, NA, false.positive.rate))%>%
  mutate(true.negative.rate = ifelse(bias.rate > 0, NA, true.negative.rate))

## plot summary stats
plo.sum <- gather(summ, index, value,- diff, -bias.rate)

## diff threshold
ggplot(plo.sum)+
  geom_point(aes(x = bias.rate, y = value, group = as.factor(diff), col = as.factor(diff)))+
  geom_line(aes(x = bias.rate, y = value, group = as.factor(diff), col = as.factor(diff)))+
  scale_colour_viridis(discrete = T, name = "Diff threshold")+
  facet_wrap(~index, scale = "fixed")+
  ylim(0,1)+
  xlab("Bias rate")+
  ylab("Rate")+
  theme_bw()+
  theme(legend.position = "bottom")

```

```{r fig8, echo = F, evaluate = T, warning = F, fig.show = T, message = F, fig.height=8, fig.width=8, fig.cap = fig$cap("fig_dataset_size", "The effect of dataset size on the `cd_ddmm` test. Datasets with less than 100 records are flagged poorly.")}
# False positive/negative rates - dataset size
dat$ds <- cut(dat$dataset.size, breaks = c(0,100, 500, 1000,  10000, 20000))
plo <- dat

#plot summary stats with non-bias = 0
plo <- plo[!is.na(plo$pass),]
plo <- plo%>%
  mutate(biased = ifelse(bias.rate == 0, F, T))%>%
  mutate(correct = pass != biased)
#remember: pass = T means not biased/test passed

summ <- plo%>%
  group_by(ds, bias.rate)%>%
  summarize(true.positive = sum(biased & !pass),
            true.negative = sum(!biased & pass),
            false.positive = sum(!biased & !pass),
            false.negative = sum(biased & pass),
            sensitivity = true.positive / (true.positive + false.negative),
            specificity = true.negative / (true.negative + false.negative),
            precision = true.positive /(true.positive + false.positive),
            detection.rate = true.positive / (true.positive + false.positive + true.negative + false.negative),
            n = n(),
            false.negative.rate = false.negative / n, 
            true.positive.rate = true.positive / n,
            true.negative.rate = true.negative / n,
            false.positive.rate = false.positive / n)

summ <- summ%>%
  dplyr::select(-true.positive, -true.negative, -false.positive, -false.negative, -n, -precision, 
                -specificity, -detection.rate, -sensitivity)%>%
  mutate(false.positive.rate = ifelse(bias.rate > 0, NA, false.positive.rate))%>%
  mutate(true.negative.rate = ifelse(bias.rate > 0, NA, true.negative.rate))


## plot summary stats
plo.sum <- gather(summ, index, value,-ds, -bias.rate)

## diff threshold
ggplot(plo.sum)+
  geom_point(aes(x = bias.rate, y = value, group = as.factor(ds), col = as.factor(ds)))+
  geom_line(aes(x = bias.rate, y = value, group = as.factor(ds), col = as.factor(ds)))+
  scale_colour_viridis(discrete = T, name = "Number of records in dataset")+
  facet_wrap(~index, scale = "fixed")+
  ylim(0,1)+
  xlab("Bias rate")+
  ylab("Rate")+
  theme_bw()+
  theme(legend.position = "bottom")
```

##Expected failure rate
The simulations suggest the expectable false positives for a given diff (p-value constant at 0.025) are low. 

```{r, results='asis', evaluate = T, echo = F}
knitr::kable(suma, caption = "The rate of false positives at various rates of diff.")
```

## A practical example
The conversion error test is implemented in the cd_ddmm function. You can easily run the test with few lines of code. By default the test requires a two degree span for each tested dataset to avoid false flags due to habitat restrictions, as for instance on islands or patchy species distributions, for instance islands of forests in grassland. Nevertheless, we consider the `cd_ddmm` test a tool to identify potentially problematic datasets rather than for automatic filtering and recommend to double check flagged datasets using summary statistics and the diagnostic plots.

```{r echo = T, evaluate = T, warning = F, message = F, fig.height=8, fig.width=8, fig.show = T}
clean <- data.frame(species = letters[1:10],
                    decimallongitude = runif(100, -180, 180),
                    decimallatitude = runif(100, -90,90),
                    dataset = "clean")

#problematic dataset
lon <- sample(0:180, size = 100, replace = TRUE) + runif(100, 0,0.59)
lat <- sample(0:90, size = 100, replace = TRUE) + runif(100, 0,0.59)

biased <-  data.frame(species = letters[1:10],
                      decimallongitude = lon,
                      decimallatitude = lat,
                      dataset = "biased")

dat <- rbind(clean, biased)

# with diagnostic plots and small matrix due to small dataset size and for visualization
par(mfrow = c(1,2))
cd_ddmm(x = dat, diagnostic = TRUE, value = "dataset", mat_size = 100)

#inspect geographic extent of the flagged dataset to exclude island or patchy habitat
min(biased$decimallongitude)
max(biased$decimallongitude)

```

The diagnostic plots clearly show the biased distribution of decimals in the biased dataset and the large geographic extent excludes islands or patchy habitats as cause of this distribution.


# Low coordinate precision due to rasterized sampling or decimal rounding (`cd_round`)
## Background
Species occurrence records with coordinates often do not represent point occurrences but are either derived from rasterized collection designs (e.g. presence/absence in a 100x100 km grid cell) or have been subject to strong decimal rounding. If the relevant meta-data are missing, these issues cannot be identified on the record level, especially if the records have been combined with precise GPS-based point occurrences into data sets of mixed precision. However, knowledge of coordinate precision and awareness of a large number of imprecise records in a data set can be crucial for downstream analyses. For instance, coordinates collected in a 100km x 100km raster might be unsuitable for species distribution modelling on the local and regional scale.

## Algorithm and implementation
The sensitivity of the algorithms presented in the main text test can be customized based on the raster regularity and the fraction of biased records in a dataset (`T1`), the geographic extent of the raster(`reg.dist.min` and `reg.dist.max`) and the number of raster nodes (`reg.out.thresh`). 

```{r, echo = F, evaluate = T, warning = F, fig.show = T, fig.height=8, fig.width=8, fig.cap = fig$cap("fig_patternCoordinates", "The diagnostic plots retrieved from `cd_round` for two datasets with 1000 records each. The upper panel shows an unbias dataset, the lower panel shows a dataset with 10% records from a rasterized sampling (a two degree raster in this case). The left panel shows the distribution of geographic coordinates, the right panel shows the autocorrelation plot. The green vertical lines mark flagged outliers, the red horizontal lines show the flagged most common distance. The logical heading indicates the flag: TRUE = not flagged and FALSE = flagged. The test identified the biased dataset.")}
simDS <- function(n.rec,
                  K = 3,
                  data.range.min = -40, data.range.max = 40,
                  sig.min = 0.1, sig.max = 5, E.res = 0.1,
                  E = 0, nrep){
  
  
  #simulate bias grid
  bi.rang <- E.res * 5
  bi.stpt <- sample(data.range.min:(data.range.max - bi.rang), size = 1)
  data.rang <- (data.range.max + 180) - (data.range.min + 180)
  if(data.rang < bi.rang){
    data.range.max = data.range.min + bi.rang
    bi.stpt <- data.range.min
  }
  
  bi.lon <- bi.stpt + rep(seq(0, E.res * 4,by = E.res), times = round(n.rec * E / 5, 0))
  
  bi.rang <- E.res * 5
  bi.stpt <- sample(data.range.min:(data.range.max - bi.rang), size = 1)
  data.rang <- (data.range.max + 90) - (data.range.min + 90)
  if(data.rang < bi.rang){
    data.range.max = data.range.min + bi.rang
    bi.stpt <- data.range.min
  }
  
  bi.lat <- bi.stpt + rep(seq(0, E.res * 4,by = E.res), times = round(n.rec * E / 5, 0))

  #longitude
  #simulate the data according to specifications
  n.dis <- K
  mu <- runif(n.dis, data.range.min, data.range.max)
  sig <- runif(n.dis, sig.min, sig.max)
  
  #clean data
  cl.lon <- msm::rtnorm(n = round(n.rec * (1 - E), 0), mean = mu, sd = sig, lower = data.range.min, upper = data.range.max)
  
  #latitude
  n.dis <- K
  mu <- runif(n.dis, data.range.min, data.range.max)
  sig <- runif(n.dis, sig.min, sig.max)
  
  #clean data
  cl.lat<- msm::rtnorm(n = round(n.rec * (1 - E), 0), mean = mu, sd = sig, lower = data.range.min, upper = data.range.max)
  
  
  #add biased data
  inp.lon <- c(cl.lon, bi.lon)
  inp.lat <- c(cl.lat, bi.lat)
  inp <- data.frame(decimallongitude = inp.lon,
                    decimallatitude = inp.lat,
                    dataset = nrep)
  
  results <- data.frame(dataset.size = n.rec,
                        data.range = paste(data.range.min, data.range.max, sep = "_"),
                        number.of.seeds = n.dis,
                        seed.means = paste(mu, collapse = "_"),
                        seed.sigma = paste(sig, collapse = "_"),
                        bias.rate = E, 
                        bias.res = E.res,
                        bias.range = E.res * 5,
                        bias.start = bi.stpt,
                        nrep = nrep,
                        stringsAsFactors = F)
  
  
  out <- list(inp, results)
  return(out)

}

# Run simulation plus analyses
#simulate data
  sims.cl <- simDS(n.rec = 1000, K = 5, data.range.min = 0,
              data.range.max = 30, sig.min = 0.1, sig.max = 2,
              E = 0, E.res = 2,
              nrep = "No bias")

  sims.bia <- simDS(n.rec = 1000, K = 5, data.range.min = 0,
              data.range.max = 30, sig.min = 0.1, sig.max = 2,
              E = 0.05, E.res = 2,
              nrep = "10% bias")
  
    #run test
  par(mfrow = c(2,2))
  
  res <- cd_round(x = sims.cl[[1]], lon = "decimallongitude", ds = "dataset", test = "lon", value = "dataset")
  res <- cd_round(x = sims.bia[[1]], lon = "decimallongitude", ds = "dataset", test = "lon", value = "dataset")
  

```


## Simulations
We simulated 100,000 datasets of species occurrences with varying number of records and degree of sample clustering. For each iteration we first draw a random number, $N \sim \Gamma(\alpha = 2, \beta = 1) *500$ for the number of records. We the simulated N latitude and longitude coordinates between 0° an 90° using $K \in [1,5)$ truncated normal distributions with $\mu_i \sim \mathcal{U}(0,90)$ and $\sigma \sim \mathcal{U}(0.1,5)$. We then added a fraction $\rho$ of biased records, where $\rho \in[0,0.6]$ (=0-60%). We assumed a rasterized bias with five nodes (i.e. a raster with five rows and five columns). We first sampled a random coordinate as origin for the raster and then sampled the four remaining nodes at a resolution $\tau \in  [0.1, 2]$. We then analysed the simulated data using the `cd_round` function with `T1` parameters between 3 and 13. `r fig$ref("fig_clustering_t7_results", link = TRUE)` shows the effect of the `T1` parameter on the sensitivity of the `ds_ddmm` test.



```{r fig10, echo = F, evaluate = T, warning = F, message = F, fig.show = T, fig.height=8, fig.width=8, fig.cap = fig$cap("fig_clustering_t7_results", "The effect of the `T1` parameter on the sensitivity of the `ds_round` test.")}

load("inst/autocorrelation_simulations.Rds")

#plot summary stats with non-bias = 0
plo <- plo[!is.na(plo$flag),]
plo <- plo%>%
  # mutate(biased = ifelse(bias.rate < 0.1, F, T))%>%
  mutate(biased = ifelse(bias.rate == 0, F, T))%>%
  mutate(correct = flag != biased)
#remember: flag = T means not biased/test passed
summ <- plo%>%
  group_by(outlier.threshold, bias.rate)%>%
  summarize(true.positive = sum(biased & !flag),
            true.negative = sum(!biased & flag),
            false.positive = sum(!biased & !flag),
            false.negative = sum(biased & flag),
            sensitivity = true.positive / (true.positive + false.negative),
            specificity = true.negative / (true.negative + false.negative),
            precision = true.positive /(true.positive + false.positive),
            detection.rate = true.positive / (true.positive + false.positive + true.negative + false.negative),
            n = n(),
            false.negative.rate = false.negative / n, 
            true.positive.rate = true.positive / n,
            true.negative.rate = true.negative / n,
            false.positive.rate = false.positive / n)

summ <- summ%>%
  dplyr::select(-true.positive, -true.negative, -false.positive, -false.negative, -n, -precision, 
         -specificity, -detection.rate, -sensitivity)%>%
  mutate(false.positive.rate = ifelse(bias.rate > 0, NA, false.positive.rate))%>%
  mutate(true.negative.rate = ifelse(bias.rate > 0, NA, true.negative.rate))


## plot summary stats
plo.sum <- gather(summ, index, value,- outlier.threshold, -bias.rate)

## outlier threshold
ggplot(plo.sum)+
  geom_point(aes(x = bias.rate, y = value, group = as.factor(outlier.threshold), col = as.factor(outlier.threshold)))+
  scale_colour_viridis(discrete = T, name = "Outlier threshold (T1)")+
  facet_wrap(~index, scale = "fixed")+
  ylim(0,1)+
  xlab("Bias rate")+
  ylab("Rate")+
  theme_bw()+
  theme(legend.position = "bottom")
```


## A practical example
The following example illustrates the use of `cd_round`.

```{r, echo = F, evaluate = T, warning = F, fig.show = T, fig.height=8, fig.width=8,  message = F}

#simulate bias grid, one degree resolution, 10% error on a 1000 records dataset
#simulate biased fraction of the data, grid resolution = 1 degree
#simulate non-biased fraction of the data
bi <- sample(3 + 0:5, size = 100, replace = TRUE)
mu <- runif(3, 0, 15)
sig <- runif(3, 0.1, 5)
cl <- rnorm(n = 900, mean = mu, sd = sig)
lon <- c(cl, bi)

bi <- sample(9:13, size = 100, replace = TRUE)
mu <- runif(3, 0, 15)
sig <- runif(3, 0.1, 5)
cl <- rnorm(n = 900, mean = mu, sd = sig)
lat <- c(cl, bi)

biased <- data.frame(decimallongitude = lon,
                  decimallatitude = lat,
                  dataset = "biased")

# simulate unbiased data
lon <- runif(n = 1000, min = -30, max = 30)
lat <- runif(n = 1000, min = -30, max = 30)

clean <- data.frame(decimallongitude = lon,
                    decimallatitude = lat, dataset = "clean")

dat <- rbind(biased, clean)

#run test, only longitude for better visualization
par(mfrow = c(2,2))
cd_round(dat, value = "dataset", test = "lon")
```

If `value = "dataset"` the output of cd_round is a table indicating the number of regular outliers (i.e. the raster nodes) and indicating which datasets have been flagged. Additionally, the diagnostic plots clearly show the rasterized pattern in the biased dataset and confirm the automated flag.
---
title: "Identifying geographic outliers"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Identifying geographic outliers}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteSuggests{rgbif}
  \usepackage[utf8]{inputenc}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE)
```

# The problem of geographic outliers
Isolated occurrence records, distant to all other recordings of a taxon--i.e. geographic outliers--often are erroneous or problematic, because might present 1) data entry errors such as switched coordinates or switched decimal signs, 2) individuals in horticulture or captivity far from natural conditions, 3) alien records outside the natural range of the species, or 4) "freak" individuals that are remnants of historical contingencies. For many analyses such records are undesirable and need to be removed. However, the automatic treatment of such records is difficult, because records from cases 3) and 4) might be relevant for some analyses, distance-based outlier detection is a "majority vote" and hence susceptible to uneven sampling, and because some species (albeit few in our experience) might have geographically disjunct distributions (for instance on oceanic islands) prone to create false positives. Furthermore, it is generally difficult to define where a species range ends and when a record should be considered and outlier. Environmental factors can be used to inform a decision on geographic outliers, however, if records are used for environmental based distribution modelling after cleaning the is a risk of circularity. 

*CoordinateCleaner* implements the `cc_outl` test to automatically flag outlier records based on geographic position by using interquartile range outlier detection ($x > IQR(x) + Q_{75} * mltpl$) on a distance matrix of all records of a species. In our experience this test can largely improve datasets and is unproblematic for most species, but needs more careful thought than the other tests implemented in *CoordinateCleaner*. The default settings of the test are conservative in flagging records, but if peripheral or satellite distributions are common in a dataset and relevant to the analyses outcome, we strongly recommend to inspect flagged records carefully.

In this tutorial we illustrate the potential and caveats of  `cc_outl` using two boundary case examples of species with strong sampling bias-- the Red squirrel and the Eurasian lynx. In the end, we present `cc_iucn` as an alternative for outlier detection for datasets with large ranges and strong sampling bias.

```{r, echo = F, warning = F, message = F}
library(CoordinateCleaner)
library(dplyr)
library(ggplot2)
library(magrittr)
library(readr)
library(rgdal)
library(rgbif)
library(sp)
library(viridis)

# load example data from local to speed up vignette building
squ <- read_delim(file = "inst/0005313-180730143533302.csv", delim = "\t", guess_max = 100000)%>%
  dplyr::select(species = specieskey, decimallongitude, decimallatitude,countrycode)%>%
  filter(!is.na(decimallongitude))%>%
  filter(!is.na(decimallatitude))

lnx <- read_delim(file = "inst/0006556-180730143533302.csv", delim = "\t", guess_max = 100000)%>%
  dplyr::select(species = specieskey, decimallongitude, decimallatitude, countrycode)%>%
  filter(!is.na(decimallongitude))%>%
  filter(!is.na(decimallatitude))
```


# 1. The case of the European lynx
The Eurasian lynx (*Lynx lynx*) is a wide-spread Eurasian carnivore. The Global Biodiversity Information Facility (www.gbif.org) comprises more than `r nrow(lnx)` records for this species, the vast majority in Europe, few records from Asia, and some records from North America and South America. While the European and Asian records correspond to the recent natural range of the species (compare for instance [here](http://www.iucnredlist.org/details/12519/0)), the American records are  doubtful and might correspond to swapped coordinates, individuals in captivity or misidentification of the Canada lynx (*Lynx canadensis*). The Eurasian lynx is an extreme case for outlier detection, because it includes true erroneous outliers (in the Americas), but also very isolated but potentially valid records (in northern China). `cc_outl` can still help to improve the data set. Since this is a large dataset, the test will automatically use a raster-based distance approximation for the calculation of the distance matrix, as for all species with more than 10,000 records.

## Testing for geographic outliers, defaults
We will first download occurrence records for *Lynx lynx* from www.gbif.org using the rgbif package.

```{r, warning = F, message = F, collapse = T, eval = F}
# Load data from GBIF
lnx <- rgbif::occ_search(scientificName = "Lynx lynx", limit = 200000,
         hasCoordinate = T, return = "data")%>%
  dplyr::select(species = key, decimallongitude = decimalLongitude, 
         decimallatitude = decimalLatitude, countrycode)

```

You can then run the default `cc_outl` test.

```{r, warning = F, message = F, collapse = T, eval = T}
# Run basic coordinate cleaning
flags <- clean_coordinates(x = lnx, tests = c("capitals",
  "centroids", "equal", "gbif", "institutions","seas","zeros"))

lnx <- lnx[flags$.summary,]

# Default outlier detection
outl <- cc_outl(lnx, value = "flagged")

plo <- data.frame(lnx, outlier =  as.factor(!outl))

# visualize occurrence records
ggplot()+
  borders(fill = "grey60")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")
```

The default settings of the outlier test flagged `r sum(!outl)` outliers out of `r nrow(lnx)` records. A comparison with the natural range of the species shows that the test performed quite well with the default settings: all American records were flagged correctly as outliers and most Eurasian records were not flagged. However there are some potential issues, since three isolated records in Asia where also flagged as potential outliers, although they likely are valid records. In this case the isolated position of those records is not an error, but probably rather due the extremely low sampling and data availability across large parts of Asia. `cc_outl` offers two options to improve the results 1) tweaking test sensitivity and 2) accounting for sampling intensity.

## Tweaking test sensitivity
We can make the outlier detection more conservative using the `mltpl` argument of `cc_outl`. We'll increase mltpl to 6, meaning that the mean distance of a record to all other records must be more than 6 times the interquantile range of the mean distance of all points.

```{r, warning = F, message = F, collapse = T}
outl_6 <- cc_outl(lnx, value = "flagged", mltpl = 6)

plo <- data.frame(lnx, outlier =  as.factor(!outl_6))

# visualize occurrence records
ggplot()+
  borders(fill = "grey60")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")
```

This improved the results and the dataset might be acceptable for most analyses. Tweaking the test sensitivity is an easy measure and only requires a general idea of how disjunct sampling in a dataset is and how important falsely flagged records (false positives for the test) compared to falsely retained records (false negatives for the test) are for downstream analyses.  Only the two very isolated records in eastern Asia remain problematic. These records represent an extreme case, due to the rare combination of the extremely large distribution range of the species and the extremely sparse sampling over huge areas. 

## Correcting for sampling intensity
`cc_outl` enables to account for uneven sampling, on a country level using the `sampling_thresh` argument. Based on the assumption that in areas (i.e. countries) of high sampling species are more likely to be recorded and hence the distance of any correct occurrence record to the next record of the same species will be small, `sampling_thresh` can avoid the flagging of records from countries with very low general availability of occurrence records. If `sampling_thresh` is larger than zero, `cc_outl` uses the density of biological occurrence records available for a country from the Global Biodiversity Information Facility as proxy of sampling in this country. Based on the statistical distribution of the record density across all countries, flagged records from countries in the lower `sampling_thres` percentile, are never considered outliers. For instance, when `sampling_thresh` equals 0.25, records from the 25% least densely sampled countries cannot be flagged as outliers.

```{r, warning = F, message = F, collapse = T}
# Run outlier test
outl_samp <- cc_outl(lnx, value = "flagged", mltpl = 6,  sampling_thresh = 0.25)

# visualize occurrence records
plo <- data.frame(lnx, outlier =  as.factor(!outl_samp))

ggplot()+
  borders(fill = "grey60")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")
```

Correcting for sampling the records in Eastern China are now correctly assigned and the entire natural range of *Lynx lynx* is included in the dataset. However a record in Bolivia is now falsely included in the dataset. The extreme case of the lynx illustrates inherent limits of automated geographic outlier testing, since a true error (the record in Bolivia) is less distant than some of the correct records (in northern China) under similar sampling conditions. 

# 2. The case of the red squirrel
The Red squirrel (*Sciurus vulgaris*) is widespread throughout Eurasia and has `r nrow(squ)` records in the Global Biodiversity Information Facility (www.gbif.org). The sampling in the natural range is highly clustered in Europe with few records across Asia. Additionally, the dataset includes multiple records that are likely erroneous in the Americas and South Africa. The biased sampling effort towards Europe, together with the large distribution range of the species make the use of geographic outlier detection for data cleaning challenging. Still, `cc_outl` can help to improve the data set. Since this is a large dataset, the test will automatically use a raster-based distance approximation for the calculation of the distance matrix, as for all species with more than 10,000 records.

Again, we will download occurrence records from www.gbif.org and then run `cc_outl` with default settings.

```{r, warning = F, message = F, collapse = T, eval = F}
# Load data from GBIF
squ <- rgbif::occ_search(scientificName = "Lynx lynx", limit = 200000,
         hasCoordinate = T)

squ <- squ$data%>%
  dplyr::select(species = key, decimallongitude, decimallatitude, countrycode,
         hasCoordinate = T)

```

```{r, warning = F, message = F, collapse = T, eval = T}
# Run basic coordinate cleaning
flags <- clean_coordinates(x = squ, tests = c("capitals",
  "centroids", "equal", "gbif", "institutions","seas","zeros"))

squ <- squ[flags$.summary,]

# Default outlier detection
outl <- cc_outl(squ, value = "flagged")

plo <- data.frame(squ, outlier =  as.factor(!outl))

# visualize occurrence records
ggplot()+
  borders(fill = "grey60")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")
```

The default test flagged `r sum(!flags$.summary)` of `r nrow(squ)` records. A comparison with the natural range shows that all erroneous records in the Americas and South Africa where identified correctly. Unfortunately the strong sampling bias towards Europe has also let to a flagging of potentially correct East Asian records. We will again use the sampling thresh argument to address this and avoid the false flagging of records in Eastern Asia.

```{r, warning = F, message = F, collapse = T}
# Run outlier test
outl_samp<- cc_outl(squ, value = "flagged", sampling_thresh = 0.25)

# visualize occurrence records
plo <- data.frame(squ, outlier =  as.factor(!outl_samp))

ggplot()+
  borders(fill = "grey60")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")
```

This improved the results significantly and the cleaned dataset seems a good representation of the natural distribution of *Sciurus vulgaris*. However, similar to the example of the Eurasian lynx, the result is not perfect, since the records in Korea and Japan are flagged as outliers. Hence additional or alternative data cleaning methods might be necessary.

# 3. Natural ranges instead of outlier detection
The examples above showed the potential and caveats of `cc_outl`. If the inclusion of all records from the natural range, also in the periphery (as for the most eastern Asian records for the squirrel) and the exclusion of all records outside there of is crucial for downstream analyses and no reasonable balance between Type I and Type II error can be found, other approaches to identify "outlier" records might be more appropriate. If additional information on species natural ranges is available, for instance from the [International Union for the Conservation of Nature] (http://www.iucnredlist.org/technical-documents/spatial-data) (as is the case for all amphibians, birds, mammals, and reptiles), you can use the `cc_iucn` function of *CoordinateCleaner* to exclude records outside this range. `cc_iucn` directly accepts the IUCN format:

```{r, warning = F, message = F, collapse = T}
#load the IUCN range for the Red squirrel. 
#These are not provided with CoordinateCleaner and 
#need to be downloaded seperately from www.iucn.org
sq_range <- readOGR(dsn = "inst", layer = "species_20025")
sq_range@data$species <- unique(squ$species) #replace species name by GBIF ID to synchronize betwee records and range polygon

# run natural range test
rang <- cc_iucn(x = squ, range = sq_range, value = "flagged")

# plot results
plo <- data.frame(squ, outlier =  as.factor(!rang))

nat_range <- fortify(sq_range)

ggplot()+
  borders(fill = "grey60")+
  geom_polygon(data = nat_range, aes(x = long, y = lat, group = group), fill = "green", alpha = 0.5, col = "grey50")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")
```

If no such information is available, ranges can be generated ad-hoc. `cc_iucn` can be applied for multiple species, assuming a `SpatialPolygonDataFrame` similar in structure to the data provided by the IUCN.

```{r, warning = F, message = F, collapse = T}
# Create custom range polygon for the lynx
## define the polygon shape
lx_range <- Polygon(cbind(c(-10, -10, 50, 170, 138, 83, 36, 14, -10), 
                       c(35, 67, 80, 69, 32, 21, 25, 35, 35)))
lx_range <- Polygons(list(lx_range), ID = c("A"))

## define projection
wgs84 <- "+proj=longlat +datum=WGS84 +no_defs"
lx_range <- SpatialPolygons(list(lx_range), proj4string = CRS(wgs84))

## creat dataset in same format as for the quirrel; species = species name, in this case the GBIF ID

df <- data.frame(species = 2435240, row.names = "A")
lx_range <- SpatialPolygonsDataFrame(lx_range, data = as.data.frame(df))

# run natural range test
rang <- cc_iucn(x = lnx, range = lx_range, value = "flagged")

# plot results
plo <- data.frame(lnx, outlier =  as.factor(!rang))

nat_range <- fortify(lx_range)

ggplot()+
  borders(fill = "grey60")+
  geom_polygon(data = nat_range, aes(x = long, y = lat), fill = "green", alpha = 0.5, col = "grey50")+
  geom_point(data = plo, 
             aes(x = decimallongitude, y = decimallatitude, col = outlier))+
  scale_color_viridis(discrete = T, name = "Flagged outlier")+
  coord_fixed()+
  theme_bw()+
  theme(legend.position = "bottom")

```

`cc_iucn` supports occurrence record datasets and range polygons containing multiple species, so that for instance a set of occurrence records from all mammals can clean using the IUCN shape file containing all mammal natural ranges. IUCN and custom polygons can also be combined, if their structure is adapted.

```{r, warning = F, message = F, collapse = T}
## adapt the structure of the lynx polygon
lx_range@data <- data.frame(t(rep(NA, ncol(sq_range@data))), row.names = "A")
names(lx_range@data) <- names(sq_range@data)

# Combine ranges
nat_range <- rbind(sq_range, lx_range) 

# Combine records
dat <- rbind(squ, lnx)

# run natural range test
rang <- cc_iucn(x = dat, range = nat_range, value = "flagged")
```% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write_pyrate.R
\name{write_pyrate}
\alias{write_pyrate}
\title{Create Input Files for PyRate}
\usage{
write_pyrate(
  x,
  status,
  fname,
  taxon = "accepted_name",
  min_age = "min_ma",
  max_age = "max_ma",
  trait = NULL,
  path = getwd(),
  replicates = 1,
  cutoff = NULL,
  random = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing fossil records with taxon names, ages, 
and geographic coordinates.}

\item{status}{a vector of character strings of length \code{nrow(x)}.
Indicating for each record \dQuote{extinct} or \dQuote{extant}.}

\item{fname}{a character string. The prefix to use for the output files.}

\item{taxon}{character string. The column with the taxon name. 
Default = \dQuote{accepted_name}.}

\item{min_age}{character string. The column with the minimum age. Default
= \dQuote{min_ma}.}

\item{max_age}{character string. The column with the maximum age. Default
= \dQuote{max_ma}.}

\item{trait}{a numeric vector of length \code{nrow(x)}. Indicating trait
values for each record. Optional.  Default = NULL.}

\item{path}{a character string. giving the absolute path to write the output
files. Default is the working directory.}

\item{replicates}{a numerical. The number of replicates for the randomized
age generation. See details. Default = 1.}

\item{cutoff}{a numerical. Specify a threshold to exclude fossil occurrences
with a high temporal uncertainty, i.e. with a wide temporal range between
min_age and max_age. Examples: cutoff=NULL (default; all occurrences are
kept in the data set) cutoff=5 (all occurrences with a temporal range of 5
Myr or higher are excluded from the data set)}

\item{random}{logical. Specify whether to take a random age (between MinT
and MaxT) for each occurrence or the midpoint age. Note that this option
defaults to TRUE if several replicates are generated (i.e. replicates > 1).
Examples: random = TRUE (default) random = FALSE (use midpoint ages)}
}
\value{
PyRate input files in the working directory.
}
\description{
Creates the input necessary to run Pyrate, based on a data.frame with fossil
ages (as derived e.g. from clean_fossils) and a vector of the
extinction status for each sample. Creates files in the working directory!
}
\details{
The replicate option allows the user to generate several replicates of the
data set in a single input file, each time re-drawing the ages of the
occurrences at random from uniform distributions with boundaries MinT and
MaxT. The replicates can be analysed in different runs (see PyRate command
-j) and combining the results of these replicates is a way to account for
the uncertainty of the true ages of the fossil occurrences. Examples:
replicates=1 (default, generates 1 data set), replicates=10 (generates 10
random replicates of the data set).
}
\note{
See \url{https://github.com/dsilvestro/PyRate/wiki} for more details
and tutorials on PyRate and PyRate input.
}
\examples{

minages <- runif(250, 0, 65)
exmpl <- data.frame(accepted_name = sample(letters, size = 250, replace = TRUE),
                    lng = runif(250, min = 42, max = 51),
                    lat = runif(250, min = -26, max = -11),
                    min_ma = minages,
                    max_ma = minages + runif(250, 0.1, 65))

#a vector with the status for each record, 
#make sure species are only classified as either extinct or extant, 
#otherwise the function will drop an error

status <- sample(c("extinct", "extant"), size = nrow(exmpl), replace = TRUE)

#or from a list of species
status <- sample(c("extinct", "extant"), size = length(letters), replace = TRUE)
names(status) <- letters
status <- status[exmpl$accepted_name]

\dontrun{
write_pyrate(x = exmpl,fname = "test", status = status)
}

}
\keyword{Fossil}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoordinateCleaner-package.R
\docType{data}
\name{buffsea}
\alias{buffsea}
\title{Global Coastlines buffered by -1 degree}
\source{
\url{http://www.naturalearthdata.com/downloads/10m-physical-vectors/}
}
\description{
A \code{SpatialPolygonsDataFrame} with global coastlines, with a -1 degree buffer to extent coastlines as alternative reference for \code{\link{cc_sea}}. Can be useful to identify marine species on land without flagging records in estuaries, etc.
}
\examples{

data("buffsea")
}
\keyword{gazetteers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_gbif.R
\name{cc_gbif}
\alias{cc_gbif}
\title{Identify Records Assigned to GBIF Headquarters}
\usage{
cc_gbif(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species",
  buffer = 1000,
  geod = TRUE,
  verify = FALSE,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{species}{character string. The column with the species identity. Only
required if verify = TRUE.}

\item{buffer}{numerical. The buffer around the GBIF headquarters,
where records should be flagged as problematic. Units depend on geod. Default = 100 m.}

\item{geod}{logical. If TRUE the radius is calculated
based on a sphere, buffer is in meters. If FALSE
the radius is calculated in degrees. Default = T.}

\item{verify}{logical. If TRUE records are only flagged if they are the
only record in a given species flagged close to a given reference.
If FALSE, the distance is the only criterion}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records within 0.5 degree radius around the GBIF headquarters in
Copenhagen, DK.
}
\details{
Not recommended if working with records from Denmark or the Copenhagen area.
}
\examples{

x <- data.frame(species = "A", 
                decimallongitude = c(12.58, 12.58), 
                decimallatitude = c(55.67, 30.00))
                
cc_gbif(x)
cc_gbif(x, value = "flagged")

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_inst.R
\name{cc_inst}
\alias{cc_inst}
\title{Identify Records in the Vicinity of Biodiversity Institutions}
\usage{
cc_inst(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species",
  buffer = 100,
  geod = TRUE,
  ref = NULL,
  verify = FALSE,
  verify_mltpl = 10,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{species}{character string. The column with the species identity. Only
required if verify = TRUE.}

\item{buffer}{numerical. The buffer around each institution,
where records should be flagged as problematic, in decimal
degrees.  Default = 100m.}

\item{geod}{logical. If TRUE the radius around each capital is calculated
based on a sphere, buffer is in meters and independent of latitude. If FALSE
the radius is calculated assuming planar coordinates and varies slightly with latitude,
in this case buffer is in degrees. Default = TRUE. See https://seethedatablog.wordpress.com/ 
for detail and credits.}

\item{ref}{SpatialPointsDataFrame. Providing the geographic gazetteer. Can
be any SpatialPointsDataFrame, but the structure must be identical to
\code{\link{institutions}}.  Default = \code{\link{institutions}}}

\item{verify}{logical. If TRUE, records close to institutions are only flagged,
if there are no other records of the same species in the greater vicinity 
(a radius of buffer * verify_mltpl).}

\item{verify_mltpl}{numerical. indicates the factor by which the radius for verify
exceeds the radius of the initial test. Default = 10, which might be suitable if 
geod is TRUE, but might be too large otherwise.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records assigned to the location of zoos, botanical gardens, herbaria,
universities and museums, based on a global database of ~10,000 such
biodiversity institutions. Coordinates from these locations can be related
to data-entry errors, false automated geo-reference or individuals in
captivity/horticulture.
}
\details{
Note: the buffer radius is in degrees, thus will differ slightly between
different latitudes.
}
\examples{

x <- data.frame(species = letters[1:10], 
                decimallongitude = runif(100, -180, 180), 
                decimallatitude = runif(100, -90,90))

#large buffer for demonstration, using geod = FALSE for shorter runtime              
cc_inst(x, value = "flagged", buffer = 10, geod = FALSE) 

\dontrun{
#' cc_inst(x, value = "flagged", buffer = 50000) #geod = T
}

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoordinateCleaner-package.R
\docType{data}
\name{countryref}
\alias{countryref}
\title{Country Centroids and Country Capitals}
\format{
A data frame with 5,305 observations on 13 variables.
#' \describe{ 
\item{iso3}{ISO-3 code for each country, in case of provinces also referring to the country.}
\item{iso2}{ISO-2 code for each country, in case of provinces also referring to the country.} 
\item{adm1_code}{adm code for countries and provinces.} 
\item{name}{a factor; name of the country or province.} 
\item{type}{identifying if the entry refers to a country or province level.} 
\item{centroid.lon}{Longitude of the country centroid.}
\item{centroid.lat}{Latitude of the country centroid.}
\item{capital}{Name of the country capital, empty for provinces.}
\item{capital.lon}{Longitude of the country capital.}
\item{capital.lat}{Latitude of the country capital.}
\item{area_sqkm}{The area of the country or province.}
\item{uncertaintyRadiusMeters}{The uncertainty of the country centroid.}
\item{source}{The data source. Currently only available for \url{http://geo-locate.org}}}
}
\source{
CENTRAL INTELLIGENCE AGENCY (2014) \emph{The World Factbook},
Washington, DC.

\url{https://www.cia.gov/the-world-factbook/}
\url{http://thematicmapping.org/downloads/world_borders.php}
\url{http://geo-locate.org}
}
\description{
A \code{data.frame} with coordinates of country and province centroids and country
capitals as reference for the \code{\link{clean_coordinates}}, \code{\link{cc_cen}} and \code{\link{cc_cap}} functions.
Coordinates are based on the Central Intelligence Agency World Factbook  \url{https://www.cia.gov/the-world-factbook/},
\url{http://thematicmapping.org/downloads/world_borders.php} and geolocate \url{http://geo-locate.org}.
}
\examples{

data(countryref)
head(countryref)
}
\keyword{gazetteers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_outl.R
\name{cc_outl}
\alias{cc_outl}
\title{Identify Geographic Outliers in Species Distributions}
\usage{
cc_outl(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species",
  method = "quantile",
  mltpl = 5,
  tdi = 1000,
  value = "clean",
  sampling_thresh = 0,
  verbose = TRUE,
  min_occs = 7,
  thinning = FALSE,
  thinning_res = 0.5
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{species}{character string. The column with the species name. Default
= \dQuote{species}.}

\item{method}{character string.  Defining the method for outlier
selection.  See details. One of \dQuote{distance}, \dQuote{quantile},
\dQuote{mad}.  Default = \dQuote{quantile}.}

\item{mltpl}{numeric. The multiplier of the interquartile range
(\code{method == 'quantile'}) or median absolute deviation (\code{method ==
'mad'})to identify outliers. See details.  Default = 5.}

\item{tdi}{numeric.  The minimum absolute distance (\code{method ==
'distance'}) of a record to all other records of a species to be identified
as outlier, in km. See details. Default = 1000.}

\item{value}{character string.  Defining the output value. See value.}

\item{sampling_thresh}{numeric. Cut off threshold for the sampling correction.
Indicates the quantile of sampling in which outliers should be ignored. For instance, 
if \code{sampling_thresh} == 0.25, records in the 25% worst sampled countries will 
not be flagged as outliers. Default = 0 (no sampling correction).}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}

\item{min_occs}{Minimum number of geographically unique datapoints needed for a species to be tested. 
This is necessary for reliable outlier estimation.
Species with fewer than min_occs records will not be tested and the output value will be 'TRUE'.
Default is to 7. If \code{method == 'distance'}, consider a lower threshold.}

\item{thinning}{forces a raster approximation for the distance calculation. 
This is routinely used for species with more than 10,000 records for computational reasons, 
but can be enforced for smaller datasets, which is recommended when sampling is very uneven.}

\item{thinning_res}{The resolution for the spatial thinning in decimal degrees. Default = 0.5.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes out or flags records that are outliers in geographic space according to the method
defined via the \code{method} argument. Geographic outliers often represent
erroneous coordinates, for example due to data entry errors, imprecise
geo-references, individuals in horticulture/captivity.
}
\details{
The method for outlier identification depends on the \code{method} argument.
If \dQuote{outlier}: a boxplot method is used and records are flagged as
outliers if their \emph{mean} distance to all other records of the same
species is larger than mltpl * the interquartile range of the mean distance
of all records of this species. If \dQuote{mad}: the median absolute
deviation is used. In this case a record is flagged as outlier, if the
\emph{mean} distance to all other records of the same species is larger than
the median of the mean distance of all points plus/minus the mad of the mean
distances of all records of the species * mltpl. If \dQuote{distance}:
records are flagged as outliers, if the \emph{minimum} distance to the next
record of the species is > \code{tdi}. For species with records from > 10000
unique locations a random sample of 1000 records is used for 
the distance matrix calculation. The test skips species with fewer than \code{min_occs},
 geographically unique records.

The likelihood of occurrence records being erroneous outliers is linked to the sampling effort
in any given location. To account for this, the sampling_cor option fetches 
the number of occurrence records available 
from www.gbif.org, per country as a proxy of sampling effort. The outlier test 
(the mean distance) for each records is than weighted by the log transformed 
number of records per square kilometre in this country. 
See for \url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13152}
an example and further explanation of the outlier test.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

x <- data.frame(species = letters[1:10], 
                decimallongitude = runif(100, -180, 180), 
                decimallatitude = runif(100, -90,90))
                
cc_outl(x)
cc_outl(x, method = "quantile", value = "flagged")
cc_outl(x, method = "distance", value = "flagged", tdi = 10000)
cc_outl(x, method = "distance", value = "flagged", tdi = 1000)

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoordinateCleaner-package.R
\docType{data}
\name{buffland}
\alias{buffland}
\title{Global Coastlines buffered by 1 degree}
\source{
\url{http://www.naturalearthdata.com/downloads/10m-physical-vectors/}
}
\description{
A \code{SpatialPolygonsDataFrame} with global coastlines, with a 1 degree buffer to extent coastlines as alternative reference for \code{\link{cc_sea}}. Can be useful to identify species in the sea, without flagging records in mangroves, marshes, etc.
}
\examples{

data("buffland")
}
\keyword{gazetteers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cf_range.R
\name{cf_range}
\alias{cf_range}
\title{Identify Fossils with Extreme Age Ranges}
\usage{
cf_range(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  min_age = "min_ma",
  max_age = "max_ma",
  taxon = "accepted_name",
  method = "quantile",
  mltpl = 5,
  size_thresh = 7,
  max_range = 500,
  uniq_loc = FALSE,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing fossil records with taxon names, ages, 
and geographic coordinates.}

\item{lon}{character string. The column with the longitude coordinates.
To identify unique records if \code{uniq_loc  = TRUE}.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallatitude}. To identify unique records if \code{uniq_loc  = T}.}

\item{min_age}{character string. The column with the minimum age. Default
= \dQuote{min_ma}.}

\item{max_age}{character string. The column with the maximum age. Default
= \dQuote{max_ma}.}

\item{taxon}{character string. The column with the taxon name. If
\dQuote{}, searches for outliers over the entire dataset, otherwise per
specified taxon. Default = \dQuote{accepted_name}.}

\item{method}{character string.  Defining the method for outlier
selection.  See details. Either \dQuote{quantile} or \dQuote{mad}.  Default
= \dQuote{quantile}.}

\item{mltpl}{numeric. The multiplier of the interquartile range
(\code{method == 'quantile'}) or median absolute deviation (\code{method ==
'mad'}) to identify outliers. See details.  Default = 5.}

\item{size_thresh}{numeric.  The minimum number of records needed for a
dataset to be tested. Default = 10.}

\item{max_range}{numeric. A absolute maximum time interval between min age
and max age. Only relevant for \code{method} = \dQuote{time}.}

\item{uniq_loc}{logical.  If TRUE only single records per location and time
point (and taxon if \code{taxon} != "") are used for the outlier testing.
Default = T.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records with an unexpectedly large temporal range, based on a quantile
outlier test.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

minages <- runif(n = 11, min = 0.1, max = 25)
x <- data.frame(species = c(letters[1:10], "z"),
                lng = c(runif(n = 9, min = 4, max = 16), 75, 7),
                lat = c(runif(n = 11, min = -5, max = 5)),
                min_ma = minages, 
                max_ma = minages + c(runif(n = 10, min = 0, max = 5), 25))

cf_range(x, value = "flagged", taxon = "")

}
\seealso{
Other fossils: 
\code{\link{cf_age}()},
\code{\link{cf_equal}()},
\code{\link{cf_outl}()}
}
\concept{fossils}
\keyword{Fossil}
\keyword{Temporal}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cf_age.R
\name{cf_age}
\alias{cf_age}
\title{Identify Fossils with Outlier Age}
\usage{
cf_age(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  min_age = "min_ma",
  max_age = "max_ma",
  taxon = "accepted_name",
  method = "quantile",
  size_thresh = 7,
  mltpl = 5,
  replicates = 5,
  flag_thresh = 0.5,
  uniq_loc = FALSE,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing fossil records with taxon names, ages, 
and geographic coordinates.}

\item{lon}{character string. The column with the longitude coordinates.
To identify unique records if \code{uniq_loc  = TRUE}.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallatitude}. To identify unique records if \code{uniq_loc  = T}.}

\item{min_age}{character string. The column with the minimum age. Default
= \dQuote{min_ma}.}

\item{max_age}{character string. The column with the maximum age. Default
= \dQuote{max_ma}.}

\item{taxon}{character string. The column with the taxon name. If
\dQuote{}, searches for outliers over the entire dataset, otherwise per
specified taxon. Default = \dQuote{accepted_name}.}

\item{method}{character string.  Defining the method for outlier
selection.  See details. Either \dQuote{quantile} or \dQuote{mad}.  Default
= \dQuote{quantile}.}

\item{size_thresh}{numeric.  The minimum number of records needed for a
dataset to be tested. Default = 10.}

\item{mltpl}{numeric. The multiplier of the interquartile range
(\code{method == 'quantile'}) or median absolute deviation (\code{method ==
'mad'}) to identify outliers. See details.  Default = 5.}

\item{replicates}{numeric. The number of replications for the distance
matrix calculation. See details.  Default = 5.}

\item{flag_thresh}{numeric.  The fraction of passed replicates necessary to pass the test. 
See details. Default = 0.5.}

\item{uniq_loc}{logical.  If TRUE only single records per location and time
point (and taxon if \code{taxon} != "") are used for the outlier testing.
Default = T.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records that are temporal outliers based on
interquantile ranges.
}
\details{
The outlier detection is based on an interquantile range test. A temporal
distance matrix among all records is calculated based on a single point selected by random
between the minimum and maximum age for each record. The mean distance for
each point to all neighbours is calculated and the sum of these distances
is then tested against the interquantile range and flagged as an outlier if
\eqn{x > IQR(x) + q_75 * mltpl}. The test is replicated \sQuote{replicates}
times, to account for dating uncertainty. Records are flagged as outliers
if they are flagged by a fraction of more than \sQuote{flag.thresh}
replicates. Only datasets/taxa comprising more than \sQuote{size_thresh}
records are tested. Distance are calculated as Euclidean distance.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

minages <- c(runif(n = 11, min = 10, max = 25), 62.5)
x <- data.frame(species = c(letters[1:10], rep("z", 2)),
                min_ma = minages,
                max_ma = c(minages[1:11] + runif(n = 11, min = 0, max = 5), 65))

cf_age(x, value = "flagged", taxon = "")

# unique locations only
x <- data.frame(species = c(letters[1:10], rep("z", 2)),
                decimallongitude = c(runif(n = 10, min = 4, max = 16), 75, 7),
                decimallatitude = c(runif(n = 12, min = -5, max = 5)),
                min_ma = minages, 
                max_ma = c(minages[1:11] + runif(n = 11, min = 0, max = 5), 65))

cf_age(x, value = "flagged", taxon = "", uniq_loc = TRUE)

}
\seealso{
Other fossils: 
\code{\link{cf_equal}()},
\code{\link{cf_outl}()},
\code{\link{cf_range}()}
}
\concept{fossils}
\keyword{Coordinate}
\keyword{Fossil}
\keyword{Temporal}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_cap.R
\name{cc_cap}
\alias{cc_cap}
\title{Identify Coordinates in Vicinity of Country Capitals.}
\usage{
cc_cap(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species",
  buffer = 10000,
  geod = TRUE,
  ref = NULL,
  verify = FALSE,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{species}{character string. The column with the species identity. Only
required if verify = TRUE.}

\item{buffer}{The buffer around each capital coordinate (the centre of the
city), where records should be flagged as problematic. Units depend on geod.
Default = 10 kilometres.}

\item{geod}{logical. If TRUE the radius around each capital is calculated
based on a sphere, buffer is in meters and independent of latitude. If FALSE
the radius is calculated assuming planar coordinates and varies slightly with latitude,
in this case buffer is in degrees. Default = TRUE. See https://seethedatablog.wordpress.com/ 
for detail and credits.}

\item{ref}{SpatialPointsDataFrame. Providing the geographic gazetteer. Can
be any SpatialPointsDataFrame, but the structure must be identical to
\code{\link{countryref}}.  Default = \code{\link{countryref}}.}

\item{verify}{logical. If TRUE records are only flagged if they are the
only record in a given species flagged close to a given reference.
If FALSE, the distance is the only criterion}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records within a certain radius around country capitals. Poorly
geo-referenced occurrence records in biological databases are often
erroneously geo-referenced to capitals.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

x <- data.frame(species = letters[1:10], 
                decimallongitude = runif(100, -180, 180), 
                decimallatitude = runif(100, -90,90))

cc_cap(x)
cc_cap(x, value = "flagged")

}
\seealso{
Other Coordinates: 
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_dupl.R
\name{cc_dupl}
\alias{cc_dupl}
\title{Identify Duplicated Records}
\usage{
cc_dupl(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species",
  additions = NULL,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{species}{a character string. The column with the species name. Default
= \dQuote{species}.}

\item{additions}{a vector of character strings. Additional columns to be
included in the test for duplication. For example as below, collector name
and collector number.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags duplicated records based on species name and coordinates, as well as
user-defined additional columns. True (specimen) duplicates or duplicates
from the same species can make up the bulk of records in a biological
collection database, but are undesirable for many analyses. Both can be
flagged with this function, the former given enough additional information.
}
\examples{

x <- data.frame(species = letters[1:10], 
                decimallongitude = sample(x = 0:10, size = 100, replace = TRUE), 
                decimallatitude = sample(x = 0:10, size = 100, replace = TRUE),
                collector = "Bonpl",
                collector.number = c(1001, 354),
                collection = rep(c("K", "WAG","FR", "P", "S"), 20))

cc_dupl(x, value = "flagged")
cc_dupl(x, additions = c("collector", "collector.number"))

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cd_round.R
\name{cd_round}
\alias{cd_round}
\title{Identify Datasets with Rasterized Coordinates}
\usage{
cd_round(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  ds = "dataset",
  T1 = 7,
  reg_out_thresh = 2,
  reg_dist_min = 0.1,
  reg_dist_max = 2,
  min_unique_ds_size = 4,
  graphs = TRUE,
  test = "both",
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{ds}{a character string. The column with the dataset of each record. In
case \code{x} should be treated as a single dataset, identical for all
records.  Default = \dQuote{dataset}.}

\item{T1}{numeric.  The threshold for outlier detection in a in an
interquantile range based test. This is the major parameter to specify the
sensitivity of the test: lower values, equal higher detection rate. Values
between 7-11 are recommended. Default = 7.}

\item{reg_out_thresh}{numeric. Threshold on the number of equal distances
between outlier points.  See details.  Default = 2.}

\item{reg_dist_min}{numeric.  The minimum detection distance between
outliers in degrees (the minimum resolution of grids that will be flagged).
Default = 0.1.}

\item{reg_dist_max}{numeric.  The maximum detection distance between
outliers in degrees (the maximum resolution of grids that will be flagged).
Default = 2.}

\item{min_unique_ds_size}{numeric.  The minimum number of unique locations
(values in the tested column) for datasets to be included in the test.
Default = 4.}

\item{graphs}{logical. If TRUE, diagnostic plots are produced.  Default =
TRUE.}

\item{test}{character string.  Indicates which column to test. Either
\dQuote{lat} for latitude, \dQuote{lon} for longitude, or \dQuote{both} for
both.  In the latter case datasets are only flagged if both test are failed.
Default = \dQuote{both}}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
with summary statistics and flags for each dataset (\dQuote{dataset}) or a
\code{data.frame} containing the records considered correct by the test
(\dQuote{clean}) or a logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE =
test failed/potentially problematic. Default =
\dQuote{clean}.
}
\description{
Flags datasets with periodicity patterns indicative of a rasterized
(lattice) collection scheme, as often obtain from e.g. atlas data. Using a
combination of autocorrelation and sliding-window outlier detection to
identify periodicity patterns in the data. See 
\url{https://besjournals.onlinelibrary.wiley.com/doi/full/10.1111/2041-210X.13152}
for further details and 
a description of the algorithm
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

#simulate bias grid, one degree resolution, 10\% error on a 1000 records dataset
#simulate biased fraction of the data, grid resolution = 1 degree
#simulate non-biased fraction of the data
  bi <- sample(3 + 0:5, size = 100, replace = TRUE)
  mu <- runif(3, 0, 15)
  sig <- runif(3, 0.1, 5)
  cl <- rnorm(n = 900, mean = mu, sd = sig)
  lon <- c(cl, bi)
  
  bi <- sample(9:13, size = 100, replace = TRUE)
  mu <- runif(3, 0, 15)
  sig <- runif(3, 0.1, 5)
  cl <- rnorm(n = 900, mean = mu, sd = sig)
  lat <- c(cl, bi)
  
  #add biased data
  
  inp <- data.frame(decimallongitude = lon,
                    decimallatitude = lat,
                    dataset = "test")
            
          
  #run test
  \dontrun{
  cd_round(inp, value = "dataset")
  }
  

}
\seealso{
Other Datasets: 
\code{\link{cd_ddmm}()}
}
\concept{Datasets}
\keyword{"Coordinate}
\keyword{"Dataset}
\keyword{cleaning"}
\keyword{level}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_sea.R
\name{cc_sea}
\alias{cc_sea}
\title{Identify Non-terrestrial Coordinates}
\usage{
cc_sea(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  ref = NULL,
  scale = 110,
  value = "clean",
  speedup = TRUE,
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{ref}{a SpatialPolygonsDataFrame. Providing the geographic gazetteer.
Can be any SpatialPolygonsDataFrame, but the structure must be identical to
rnaturalearth::ne_download(scale = 110, type = 'land', category = 'physical').  
Default = rnaturalearth::ne_download(scale = 110, type = 'land', category = 'physical')}

\item{scale}{the scale of the default reference, as downloaded from natural earth. 
Must be one of 10, 50, 110. Higher numbers equal higher detail. Default = 110.}

\item{value}{character string.  Defining the output value. See value.}

\item{speedup}{logical. Using heuristic to speed up the analysis for large data sets
with many records per location.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags coordinates outside the reference landmass. Can be used to restrict
datasets to terrestrial taxa, or exclude records from the open ocean, when
depending on the reference (see details). Often records of terrestrial taxa
can be found in the open ocean, mostly due to switched latitude and
longitude.
}
\details{
In some cases flagging records close of the coastline is not recommendable,
because of the low precision of the reference dataset, minor GPS imprecision
or because a dataset might include coast or marshland species. If you only
want to flag records in the open ocean, consider using a buffered landmass
reference, e.g.: \code{\link{buffland}}.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{
x <- data.frame(species = letters[1:10], 
                decimallongitude = runif(10, -30, 30), 
                decimallatitude = runif(10, -30, 30))
                
cc_sea(x, value = "flagged")

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_zero.R
\name{cc_zero}
\alias{cc_zero}
\title{Identify Zero Coordinates}
\usage{
cc_zero(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  buffer = 0.5,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{buffer}{numerical. The buffer around the 0/0 point,
where records should be flagged as problematic, in decimal
degrees.  Default = 0.1.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records with either zero longitude or latitude and a radius around the
point at zero longitude and zero latitude. These problems are often due to
erroneous data-entry or geo-referencing and can lead to typical patterns of
high diversity around the equator.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

x <- data.frame(species = "A", 
                decimallongitude = c(0,34.84, 0, 33.98), 
                decimallatitude = c(23.08, 0, 0, 15.98))
                
cc_zero(x)
cc_zero(x, value = "flagged")

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_dataset.R
\name{clean_dataset}
\alias{clean_dataset}
\title{Coordinate Cleaning using Dataset Properties}
\usage{
clean_dataset(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  ds = "dataset",
  tests = c("ddmm", "periodicity"),
  value = "dataset",
  verbose = TRUE,
  ...
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{ds}{a character string. The column with the dataset of each record. In
case \code{x} should be treated as a single dataset, identical for all
records. Default = \dQuote{dataset}.}

\item{tests}{a vector of character strings, indicating which tests to run.
See details for all tests available. Default = c("ddmm", "periodicity")}

\item{value}{a character string.  Defining the output value. See value.
Default = \dQuote{dataset}.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}

\item{...}{additional arguments to be passed to \code{\link{cd_ddmm}} and
\code{\link{cd_round}} to customize test sensitivity.}
}
\value{
Depending on the \sQuote{value} argument:
\describe{
\item{\dQuote{dataset}}{a \code{data.frame} with the
the test summary statistics for each dataset in \code{x}}
\item{\dQuote{clean}}{a \code{data.frame} containing only
records from datasets in \code{x} that passed the tests}
\item{\dQuote{flagged}}{a logical vector of the same length as
rows in \code{x}, with TRUE = test passed and
FALSE = test failed/potentially problematic.}
}
}
\description{
Tests for problems associated with coordinate conversions and rounding,
based on dataset properties. Includes test to identify contributing datasets with
potential errors with converting ddmm to dd.dd, and
periodicity in the data decimals indicating rounding or a raster basis
linked to low coordinate precision. Specifically:
\itemize{
\item ddmm  tests for erroneous conversion from a degree
minute format (ddmm) to a decimal degree (dd.dd) format
\item periodicity test for periodicity in the data,
which can indicate imprecise coordinates, due to rounding or rasterization.
}
}
\details{
These tests are based on the statistical distribution of coordinates and
their decimals within
datasets of geographic distribution records to identify datasets with
potential errors/biases. Three potential error sources can be identified.
The ddmm flag tests for the particular pattern that emerges if geographical
coordinates in a degree minute annotation are transferred into decimal
degrees, simply replacing the degree symbol with the decimal point. This
kind of problem has been observed by in older datasets first recorded on
paper using typewriters, where e.g. a floating point was used as symbol for
degrees. The function uses a binomial test to check if more records than
expected have decimals below 0.6 (which is the maximum that can be obtained
in minutes, as one degree has 60 minutes) and if the number of these records
is higher than those above 0.59 by a certain proportion. The periodicity
test uses rate estimation in a Poisson process to estimate if there is
periodicity in the decimals of a dataset (as would be expected by for
example rounding or data that was collected in a raster format) and if there
is an over proportional number of records with the decimal 0 (full degrees)
which indicates rounding and thus low precision. The default values are
empirically optimized by with GBIF data, but should probably be adapted.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more details
and tutorials.
}
\examples{
#Create test dataset
clean <- data.frame(dataset = rep("clean", 1000),
                    decimallongitude = runif(min = -43, max = -40, n = 1000),
                    decimallatitude = runif(min = -13, max = -10, n = 1000))
                    
bias.long <- c(round(runif(min = -42, max = -40, n = 500), 1),
               round(runif(min = -42, max = -40, n = 300), 0),
               runif(min = -42, max = -40, n = 200))
bias.lat <- c(round(runif(min = -12, max = -10, n = 500), 1),
              round(runif(min = -12, max = -10, n = 300), 0),
              runif(min = -12, max = -10, n = 200))
bias <- data.frame(dataset = rep("biased", 1000),
                   decimallongitude = bias.long,
                   decimallatitude = bias.lat)
test <- rbind(clean, bias)

\dontrun{                  
#run clean_dataset
flags <- clean_dataset(test)

#check problems
#clean
hist(test[test$dataset == rownames(flags[flags$summary,]), "decimallongitude"])
#biased
hist(test[test$dataset == rownames(flags[!flags$summary,]), "decimallongitude"])

}
}
\seealso{
\code{\link{cd_ddmm}} \code{\link{cd_round}}

Other Wrapper functions: 
\code{\link{clean_coordinates}()},
\code{\link{clean_fossils}()}
}
\concept{Wrapper functions}
\keyword{Coordinate}
\keyword{cleaning}
\keyword{wrapper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_coordinates.R
\name{clean_coordinates}
\alias{clean_coordinates}
\alias{summary.spatialvalid}
\alias{is.spatialvalid}
\title{Geographic Cleaning of Coordinates from Biologic Collections}
\usage{
clean_coordinates(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species",
  countries = NULL,
  tests = c("capitals", "centroids", "equal", "gbif", "institutions", "outliers",
    "seas", "zeros"),
  capitals_rad = 10000,
  centroids_rad = 1000,
  centroids_detail = "both",
  inst_rad = 100,
  outliers_method = "quantile",
  outliers_mtp = 5,
  outliers_td = 1000,
  outliers_size = 7,
  range_rad = 0,
  zeros_rad = 0.5,
  capitals_ref = NULL,
  centroids_ref = NULL,
  country_ref = NULL,
  country_refcol = "iso_a3",
  inst_ref = NULL,
  range_ref = NULL,
  seas_ref = NULL,
  seas_scale = 50,
  urban_ref = NULL,
  value = "spatialvalid",
  verbose = TRUE,
  report = FALSE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{species}{a character string. A vector of the same length as rows in x,
with the species identity for each record.  If NULL, \code{tests} must not
include the "outliers" or "duplicates" tests.}

\item{countries}{a character string. The column with the country assignment of
each record in three letter ISO code. Default = \dQuote{countrycode}. If missing, the
countries test is skipped.}

\item{tests}{a vector of character strings, indicating which tests to run.
See details for all tests available. Default = c("capitals", "centroids",
"equal", "gbif", "institutions", "outliers",
"seas", "zeros")}

\item{capitals_rad}{numeric. The radius around capital coordinates in
meters. Default = 10000.}

\item{centroids_rad}{numeric. The radius around centroid coordinates in
meters. Default = 1000.}

\item{centroids_detail}{a \code{character string}. If set to
\sQuote{country} only country (adm-0) centroids are tested, if set to
\sQuote{provinces} only province (adm-1) centroids are tested.  Default =
\sQuote{both}.}

\item{inst_rad}{numeric. The radius around biodiversity institutions
coordinates in metres. Default = 100.}

\item{outliers_method}{The method used for outlier testing. See details.}

\item{outliers_mtp}{numeric. The multiplier for the interquartile range of
the outlier test.  If NULL \code{outliers.td} is used.  Default = 5.}

\item{outliers_td}{numeric.  The minimum distance of a record to all other
records of a species to be identified as outlier, in km. Default = 1000.}

\item{outliers_size}{numerical.  The minimum number of records in a dataset
to run the taxon-specific outlier test.  Default = 7.}

\item{range_rad}{buffer around natural ranges. Default = 0.}

\item{zeros_rad}{numeric. The radius around 0/0 in degrees. Default = 0.5.}

\item{capitals_ref}{a \code{data.frame} with alternative reference data for
the country capitals test. If missing, the \code{countryref} dataset is used.
Alternatives must be identical in structure.}

\item{centroids_ref}{a \code{data.frame} with alternative reference data for
the centroid test. If NULL, the \code{countryref} dataset is used.
Alternatives must be identical in structure.}

\item{country_ref}{a \code{SpatialPolygonsDataFrame} as alternative
reference for the countries test. If NULL, the
\code{rnaturalearth:ne_countries('medium')} dataset is used.}

\item{country_refcol}{the column name in the reference dataset, containing the relevant
ISO codes for matching. Default is to "iso_a3_eh" which referes to the ISO-3
codes in the reference dataset. See notes.}

\item{inst_ref}{a \code{data.frame} with alternative reference data for the
biodiversity institution test. If NULL, the \code{institutions} dataset
is used.  Alternatives must be identical in structure.}

\item{range_ref}{a \code{SpatialPolygonsDataFrame} of species natural ranges.
Required to include the 'ranges' test. See \code{\link{cc_iucn}} for details.}

\item{seas_ref}{a \code{SpatialPolygonsDataFrame} as alternative reference
for the seas test. If NULL, the
rnaturalearth::ne_download(=scale = 110, type = 'land', category = 'physical')
dataset is used.}

\item{seas_scale}{The scale of the default landmass reference. Must be one of 10, 50, 110.
Higher numbers equal higher detail. Default = 50.}

\item{urban_ref}{a \code{SpatialPolygonsDataFrame} as alternative reference
for the urban test. If NULL, the test is skipped. See details for a
reference gazetteers.}

\item{value}{a character string defining the output value. See the value
section for details. one of \sQuote{spatialvalid}, \sQuote{summary},
\sQuote{clean}. Default = \sQuote{\code{spatialvalid}}.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}

\item{report}{logical or character.  If TRUE a report file is written to the
working directory, summarizing the cleaning results. If a character, the
path to which the file should be written.  Default = FALSE.}
}
\value{
Depending on the output argument:
\describe{
\item{\dQuote{spatialvalid}}{an object of class \code{spatialvalid} similar to x
with one column added for each test. TRUE = clean coordinate entry, FALSE = potentially
problematic coordinate entries.  The .summary column is FALSE if any test flagged
the respective coordinate.}
\item{\dQuote{flagged}}{a logical vector with the
same order as the input data summarizing the results of all test. TRUE =
clean coordinate, FALSE = potentially problematic (= at least one test
failed).}
\item{\dQuote{clean}}{a \code{data.frame} similar to x
with potentially problematic records removed}
}
}
\description{
Cleaning geographic coordinates by multiple empirical tests to flag
potentially erroneous coordinates, addressing issues common in biological
collection databases.
}
\details{
The function needs all coordinates to be formally valid according to WGS84.
If the data contains invalid coordinates, the function will stop and return
a vector flagging the invalid records. TRUE = non-problematic coordinate,
FALSE = potentially problematic coordinates.
\itemize{
\item capitals tests a radius around adm-0 capitals. The
radius is \code{capitals_rad}.
\item centroids tests a radius around country centroids.
The radius is \code{centroids_rad}.
\item countries tests if coordinates are from the
country indicated in the country column.  \emph{Switched off by default.}
\item duplicates tests for duplicate records. This
checks for identical coordinates or if a species vector is provided for
identical coordinates within a species. All but the first records are
flagged as duplicates. \emph{Switched off by default.}
\item equal tests for equal absolute longitude and latitude.
\item gbif tests a one-degree radius around the GBIF
headquarters in Copenhagen, Denmark.
\item institutions tests a radius around known
biodiversity institutions from \code{instiutions}. The radius is
\code{inst_rad}.
\item outliers tests each species for outlier records.
Depending on the \code{outliers_mtp} and \code{outliers.td} arguments either
flags records that are a minimum distance away from all other records of
this species (\code{outliers_td}) or records that are outside a multiple of
the interquartile range of minimum distances to the next neighbour of this
species (\code{outliers_mtp}). Three different methods are available
for the outlier test: "If
\dQuote{outlier} a boxplot method is used and records are flagged as
outliers if their \emph{mean} distance to all other records of the same
species is larger than mltpl * the interquartile range of the mean distance
of all records of this species. If \dQuote{mad} the median absolute
deviation is used. In this case a record is flagged as outlier, if the
\emph{mean} distance to all other records of the same species is larger than
the median of the mean distance of all points plus/minus the mad of the mean
distances of all records of the species * mltpl. If \dQuote{distance}
records are flagged as outliers, if the \emph{minimum} distance to the next
record of the species is > \code{tdi}.
\item ranges tests if records fall within provided natural range polygons on
a per species basis. See \code{\link{cc_iucn}} for details.
\item seas tests if coordinates fall into the ocean.
\item urban tests if coordinates are from urban areas.
\emph{Switched off by default}
\item validity checks if coordinates correspond to a lat/lon coordinate reference system.
This test is always on, since all records need to pass for any other test to run.
\item zeros tests for plain zeros, equal latitude and
longitude and a radius around the point 0/0. The radius is \code{zeros.rad}.
}
}
\note{
Always tests for coordinate validity: non-numeric or missing
coordinates and coordinates exceeding the global extent (lon/lat, WGS84).
See \url{https://ropensci.github.io/CoordinateCleaner/} for more details
and tutorials.

The country_refcol argument allows to adapt the function to the structure of
alternative reference datasets. For instance, for
\code{rnaturalearth::ne_countries(scale = "small")}, the default will fail,
but country_refcol = "iso_a3" will work.
}
\examples{


exmpl <- data.frame(species = sample(letters, size = 250, replace = TRUE),
                    decimallongitude = runif(250, min = 42, max = 51),
                    decimallatitude = runif(250, min = -26, max = -11))

test <- clean_coordinates(x = exmpl, 
                          tests = c("equal"))
                                    
\dontrun{
#run more tests
test <- clean_coordinates(x = exmpl, 
                          tests = c("capitals", 
                          "centroids","equal", 
                          "gbif", "institutions", 
                          "outliers", "seas", 
                          "zeros"))
}
                                 
                                    
summary(test)

}
\seealso{
Other Wrapper functions: 
\code{\link{clean_dataset}()},
\code{\link{clean_fossils}()}
}
\concept{Wrapper functions}
\keyword{Coordinate}
\keyword{cleaning}
\keyword{wrapper}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_urb.R
\name{cc_urb}
\alias{cc_urb}
\title{Identify Records Inside Urban Areas}
\usage{
cc_urb(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  ref = NULL,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{ref}{a SpatialPolygonsDataFrame. Providing the geographic gazetteer
with the urban areas. See details. By default 
rnaturalearth::ne_download(scale = 'medium', type = 'urban_areas').
Can be any \code{SpatialPolygonsDataframe}, but the structure must be
identical to \code{rnaturalearth::ne_download()}.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records from inside urban areas, based on a geographic gazetteer.
Often records from large databases span substantial time periods (centuries)
and old records might represent habitats which today are replaced by city
area.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

\dontrun{
x <- data.frame(species = letters[1:10],
                decimallongitude = runif(100, -180, 180),
                decimallatitude = runif(100, -90,90))

cc_urb(x)
cc_urb(x, value = "flagged")
}

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_equ.R
\name{cc_equ}
\alias{cc_equ}
\title{Identify Records with Identical lat/lon}
\usage{
cc_equ(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  test = "absolute",
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{test}{character string. Defines if coordinates are compared exactly
(\dQuote{identical}) or on the absolute scale (i.e. -1 = 1,
\dQuote{absolute}). Default is to \dQuote{absolute}.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records with equal latitude and longitude coordinates, either exact or
absolute. Equal coordinates can often indicate data entry errors.
}
\examples{

x <- data.frame(species = letters[1:10], 
                decimallongitude = runif(100, -180, 180), 
                decimallatitude = runif(100, -90,90))

cc_equ(x)
cc_equ(x, value = "flagged")

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_coun.R
\name{cc_coun}
\alias{cc_coun}
\title{Identify Coordinates Outside their Reported Country}
\usage{
cc_coun(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  iso3 = "countrycode",
  value = "clean",
  ref = NULL,
  ref_col = "iso_a3",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{iso3}{a character string. The column with the country assignment of
each record in three letter ISO code. Default = \dQuote{countrycode}.}

\item{value}{character string.  Defining the output value. See value.}

\item{ref}{a SpatialPolygonsDataFrame. Providing the geographic gazetteer.
Can be any SpatialPolygonsDataFrame, but the structure must be identical to
\code{rnaturalearth::ne_countries(scale = "medium")}.  
Default = \code{rnaturalearth::ne_countries(scale = "medium")}}

\item{ref_col}{the column name in the reference dataset, containing the relevant
ISO codes for matching. Default is to "iso_a3_eh" which refers to the ISO-3
codes in the reference dataset. See notes.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags mismatches between geographic coordinates and additional country
information (usually this information is reliably reported with specimens).
Such a mismatch can occur for example, if latitude and longitude are
switched.
}
\note{
The ref_col argument allows to adapt the function to the structure of
alternative reference datasets. For instance, for 
\code{rnaturalearth::ne_countries(scale = "small")}, the default will fail, 
but ref_col = "iso_a3" will work.

With the default reference, records are flagged if they fall 
outside the terrestrial territory of countries, hence records in territorial waters might be flagged. 
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

\dontrun{
x <- data.frame(species = letters[1:10], 
                decimallongitude = runif(100, -20, 30), 
                decimallatitude = runif(100, 35,60),
                countrycode = "RUS")

cc_coun(x, value = "flagged")#non-terrestrial records are flagged as wrong. 
}

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_val.R
\name{cc_val}
\alias{cc_val}
\title{Identify Invalid lat/lon Coordinates}
\usage{
cc_val(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags non-numeric and not available coordinates
as well as lat >90, la <-90, lon > 180 and lon < -180 are flagged.
}
\details{
This test is obligatory before running any further tests of
CoordinateCleaner, as additional tests only run with valid coordinates.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

x <- data.frame(species = letters[1:10], 
                decimallongitude = c(runif(106, -180, 180), NA, "13W33'", "67,09", 305), 
                decimallatitude = runif(110, -90,90))
                
cc_val(x)
cc_val(x, value = "flagged")

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cd_ddmm.R
\name{cd_ddmm}
\alias{cd_ddmm}
\title{Identify Datasets with a Degree Conversion Error}
\usage{
cd_ddmm(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  ds = "dataset",
  pvalue = 0.025,
  diff = 1,
  mat_size = 1000,
  min_span = 2,
  value = "clean",
  verbose = TRUE,
  diagnostic = FALSE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{ds}{a character string. The column with the dataset of each record. In
case \code{x} should be treated as a single dataset, identical for all
records.  Default = \dQuote{dataset}.}

\item{pvalue}{numeric. The p-value for the one-sided t-test to flag the test
as passed or not. Both ddmm.pvalue and diff must be met. Default = 0.025.}

\item{diff}{numeric. The threshold difference for the ddmm test. Indicates
by which fraction the records with decimals below 0.6 must outnumber the
records with decimals above 0.6. Default = 1}

\item{mat_size}{numeric. The size of the matrix for the binomial test. Must
be changed in decimals (e.g. 100, 1000, 10000). Adapt to dataset size,
generally 100 is better for datasets < 10000 records, 1000 is better for
datasets with 10000 - 1M records. Higher values also work reasonably well
for smaller datasets, therefore, default = 1000. For large datasets try
10000.}

\item{min_span}{numeric. The minimum geographic extent of datasets to be
tested. Default = 2.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}

\item{diagnostic}{logical. If TRUE plots the analyses matrix for each
dataset.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
with summary statistics and flags for each dataset (\dQuote{dataset}) or a
\code{data.frame} containing the records considered correct by the test
(\dQuote{clean}) or a logical vector (\dQuote{flags}), with TRUE = test passed and FALSE =
test failed/potentially problematic. Default =
\dQuote{clean}.
}
\description{
This test flags datasets where a significant fraction of records has
been subject to a common degree minute to decimal degree conversion error,
where the degree sign is recognized as decimal delimiter.
}
\details{
If the degree sign is recognized as decimal delimiter during coordinate
conversion, no coordinate decimals above 0.59 (59') are possible. The test
here uses a binomial test to test if a significant proportion of records in
a dataset have been subject to this problem. The test is best adjusted via
the diff argument. The lower \code{diff}, the stricter the test. Also scales
with dataset size. Empirically, for datasets with < 5,000 unique coordinate
records \code{diff = 0.1} has proven reasonable flagging most datasets with
>25\% problematic records and all dataset with >50\% problematic records.
For datasets between 5,000 and 100,000 geographic unique records \code{diff
= 0.01} is recommended, for datasets between 100,000 and 1 M records diff =
0.001, and so on.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

clean <- data.frame(species = letters[1:10], 
                decimallongitude = runif(100, -180, 180), 
                decimallatitude = runif(100, -90,90),
                dataset = "FR")
                
cd_ddmm(x = clean, value = "flagged")

#problematic dataset
lon <- sample(0:180, size = 100, replace = TRUE) + runif(100, 0,0.59)
lat <- sample(0:90, size = 100, replace = TRUE) + runif(100, 0,0.59)

prob <-  data.frame(species = letters[1:10], 
                decimallongitude = lon, 
                decimallatitude = lat,
                dataset = "FR")
                
cd_ddmm(x = prob, value = "flagged")

}
\seealso{
Other Datasets: 
\code{\link{cd_round}()}
}
\concept{Datasets}
\keyword{"Coordinate}
\keyword{"Dataset}
\keyword{cleaning"}
\keyword{level}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_iucn.R
\name{cc_iucn}
\alias{cc_iucn}
\title{Identify Records Outside Natural Ranges}
\usage{
cc_iucn(
  x,
  range,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species",
  buffer = 0,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{range}{a SpatialPolygonsDataFrame of natural ranges for species in x. 
Must contain a column named as indicated by \code{species}. See details.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{species}{a character string. The column with the species name. 
Default = \dQuote{species}.}

\item{buffer}{numerical. The buffer around each species' range,
from where records should be flagged as problematic, in decimal
degrees. Default = 0.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records outside of the provided natural range polygon, on a per species basis. 
Expects one entry per species. See the example or 
\url{https://www.iucnredlist.org/resources/spatial-data-download} for 
the required polygon structure.
}
\details{
Download natural range maps in suitable format for amphibians, birds,
mammals and reptiles
from \url{https://www.iucnredlist.org/resources/spatial-data-download}.
Note: the buffer radius is in degrees, thus will differ slightly between
different latitudes.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{
require(sp)

x <- data.frame(species = c("A", "B"),
decimallongitude = runif(100, -170, 170),
decimallatitude = runif(100, -80,80))

range_species_A <- Polygon(cbind(c(-45,-45,-60,-60,-45),c(-10,-25,-25,-10,-10)))
range_species_B <- Polygon(cbind(c(15,15,32,32,15),c(10,-10,-10,10,10)))
range_A <- Polygons(list(range_species_A), ID = c("A"))
range_B <- Polygons(list(range_species_B), ID = c("B"))
range <- SpatialPolygons(list(range_A, range_B))
df <- data.frame(species = c("A", "B"), row.names = c("A", "B"))
range <- SpatialPolygonsDataFrame(range, data = as.data.frame(df))

cc_iucn(x = x, range = range, buffer = 10)

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_cen}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/clean_fossils.R
\name{clean_fossils}
\alias{clean_fossils}
\title{Geographic and Temporal Cleaning of Records from Fossil Collections}
\usage{
clean_fossils(
  x,
  lon = "lng",
  lat = "lat",
  min_age = "min_ma",
  max_age = "max_ma",
  taxon = "accepted_name",
  tests = c("agesequal", "centroids", "equal", "gbif", "institutions", "spatiotemp",
    "temprange", "validity", "zeros"),
  countries = NULL,
  centroids_rad = 0.05,
  centroids_detail = "both",
  inst_rad = 0.001,
  outliers_method = "quantile",
  outliers_threshold = 5,
  outliers_size = 7,
  outliers_replicates = 5,
  zeros_rad = 0.5,
  centroids_ref = NULL,
  country_ref = NULL,
  inst_ref = NULL,
  value = "spatialvalid",
  verbose = TRUE,
  report = FALSE
)
}
\arguments{
\item{x}{data.frame. Containing fossil records, containing taxon names, ages,
and geographic coordinates..}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{min_age}{character string. The column with the minimum age. Default
= \dQuote{min_ma}.}

\item{max_age}{character string. The column with the maximum age. Default
= \dQuote{max_ma}.}

\item{taxon}{character string. The column with the taxon name. If
\dQuote{}, searches for outliers over the entire dataset, otherwise per
specified taxon. Default = \dQuote{accepted_name}.}

\item{tests}{vector of character strings, indicating which tests to run.
See details for all tests available. Default = c("centroids",
"equal", "gbif", "institutions", "temprange", "spatiotemp", "agesequal", "zeros")}

\item{countries}{a character string. The column with the country assignment of
each record in three letter ISO code. Default = \dQuote{countrycode}. If missing, the
countries test is skipped.}

\item{centroids_rad}{numeric. The radius around centroid coordinates in
meters. Default = 1000.}

\item{centroids_detail}{a \code{character string}. If set to
\sQuote{country} only country (adm-0) centroids are tested, if set to
\sQuote{provinces} only province (adm-1) centroids are tested.  Default =
\sQuote{both}.}

\item{inst_rad}{numeric. The radius around biodiversity institutions
coordinates in metres. Default = 100.}

\item{outliers_method}{The method used for outlier testing. See details.}

\item{outliers_threshold}{numerical.  The multiplier for the interquantile
range for outlier detection. The higher the number, the more conservative
the outlier tests.  See \code{\link{cf_outl}} for details. Default = 3.}

\item{outliers_size}{numerical.  The minimum number of records in a dataset
to run the taxon-specific outlier test.  Default = 7.}

\item{outliers_replicates}{numeric. The number of replications for the
distance matrix calculation. See details.  Default = 5.}

\item{zeros_rad}{numeric. The radius around 0/0 in degrees. Default = 0.5.}

\item{centroids_ref}{a \code{data.frame} with alternative reference data for
the centroid test. If NULL, the \code{countryref} dataset is used.
Alternatives must be identical in structure.}

\item{country_ref}{a \code{SpatialPolygonsDataFrame} as alternative
reference for the countries test. If NULL, the
\code{rnaturalearth:ne_countries('medium')} dataset is used.}

\item{inst_ref}{a \code{data.frame} with alternative reference data for the
biodiversity institution test. If NULL, the \code{institutions} dataset
is used.  Alternatives must be identical in structure.}

\item{value}{a character string defining the output value. See the value
section for details. one of \sQuote{spatialvalid}, \sQuote{summary},
\sQuote{clean}. Default = \sQuote{\code{spatialvalid}}.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}

\item{report}{logical or character.  If TRUE a report file is written to the
working directory, summarizing the cleaning results. If a character, the
path to which the file should be written.  Default = FALSE.}
}
\value{
Depending on the output argument:
\describe{
\item{\dQuote{spatialvalid}}{an object of class \code{spatialvalid} similar to x
with one column added for each test. TRUE = clean coordinate entry, FALSE = potentially
problematic coordinate entries.  The .summary column is FALSE if any test flagged
the respective coordinate.}
\item{\dQuote{flagged}}{a logical vector with the
same order as the input data summarizing the results of all test. TRUE =
clean coordinate, FALSE = potentially problematic (= at least one test
failed).}
\item{\dQuote{clean}}{a \code{data.frame} similar to x
with potentially problematic records removed}
}
}
\description{
Cleaning records by multiple empirical tests to flag potentially erroneous
coordinates and time-spans, addressing issues common in fossil collection
databases. Individual tests can be activated via the tests argument:
}
\details{
\itemize{
\item agesequal tests for equal minimum and maximum age.
\item centroids tests a radius around country centroids.
The radius is \code{centroids_rad}.
\item countries tests if coordinates are from the
country indicated in the country column.  \emph{Switched off by default.}
\item equal tests for equal absolute longitude and latitude.
\item gbif tests a one-degree radius around the GBIF
headquarters in Copenhagen, Denmark.
\item institutions tests a radius around known
biodiversity institutions from \code{instiutions}. The radius is
\code{inst_rad}.
\item spatiotemp test for records which are outlier in time and space. See below for details.
\item temprange tests for records with unexpectedly large temporal ranges,
using a quantile-based outlier test.
\item validity checks if coordinates correspond to a lat/lon coordinate reference system.
This test is always on, since all records need to pass for any other test to run.
\item zeros tests for plain zeros, equal latitude and
longitude and a radius around the point 0/0. The radius is \code{zeros_rad}.
The outlier detection in \sQuote{spatiotemp} is based on an interquantile range test. In a first
step a distance matrix of geographic distances among all records is
calculate. Subsequently a similar distance matrix of temporal distances
among all records is calculated based on a single point selected by random
between the minimum and maximum age for each record. The mean distance for
each point to all neighbours is calculated for both matrices and spatial and
temporal distances are scaled to the same range. The sum of these distanced
is then tested against the interquantile range and flagged as an outlier if
\eqn{x > IQR(x) + q_75 * mltpl}. The test is replicated \sQuote{replicates}
times, to account for temporal uncertainty. Records are flagged as outliers
if they are flagged by a fraction of more than \sQuote{flag_thresh}
replicates. Only datasets/taxa comprising more than \sQuote{size.thresh}
records are tested. Note that geographic distances are calculated as
geospheric distances for datasets (or taxa) with fewer than 10,000 records
and approximated as Euclidean distances for datasets/taxa with 10,000 to
25,000 records. Datasets/taxa comprising more than 25,000 records are
skipped.
}
}
\note{
Always tests for coordinate validity: non-numeric or missing
coordinates and coordinates exceeding the global extent (lon/lat, WGS84).

See \url{https://ropensci.github.io/CoordinateCleaner/} for more details
and tutorials.
}
\examples{

minages <- runif(250, 0, 65)
exmpl <- data.frame(accepted_name = sample(letters, size = 250, replace = TRUE),
                    lng = runif(250, min = 42, max = 51),
                    lat = runif(250, min = -26, max = -11),
                    min_ma = minages,
                    max_ma = minages + runif(250, 0.1, 65))

test <- clean_fossils(x = exmpl)

summary(test)

}
\seealso{
Other Wrapper functions: 
\code{\link{clean_coordinates}()},
\code{\link{clean_dataset}()}
}
\concept{Wrapper functions}
\keyword{Coordinate}
\keyword{Fossil}
\keyword{Temporal}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoordinateCleaner-package.R
\docType{data}
\name{pbdb_example}
\alias{pbdb_example}
\title{Example data from the Paleobiologydatabase}
\format{
A data frame with 5000 observations on 36 variables.
}
\source{
\itemize{ 
\item The Paleobiology database \url{https://paleobiodb.org/} 
\item Sara Varela, Javier Gonzalez Hernandez and Luciano Fabris Sgarbi (2016). 
paleobioDB: Download and Process Data from the Paleobiology Database. 
R package version 0.5.0. \url{https://CRAN.R-project.org/package=paleobioDB}.
}
}
\description{
A dataset of 5000 flowering plant fossil occurrences as example for data of the paleobiology Database, downloaded using the paleobioDB packages as specified in the vignette \dQuote{Cleaning_PBDB_fossils_with_CoordinateCleaner}.
}
\examples{

data(institutions)
str(institutions)

}
\keyword{gazetteers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cf_equal.R
\name{cf_equal}
\alias{cf_equal}
\title{Identify Fossils with equal min and max age}
\usage{
cf_equal(
  x,
  min_age = "min_ma",
  max_age = "max_ma",
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing fossil records with taxon names, ages, 
and geographic coordinates.}

\item{min_age}{character string. The column with the minimum age. Default
= \dQuote{min_ma}.}

\item{max_age}{character string. The column with the maximum age. Default
= \dQuote{max_ma}.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records with equal minimum and maximum age.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

minages <- runif(n = 10, min = 0.1, max = 25)
x <- data.frame(species = letters[1:10], 
                min_ma = minages, 
                max_ma = minages + runif(n = 10, min = 0, max = 10))
x <- rbind(x, data.frame(species = "z", 
                min_ma = 5, 
                max_ma = 5))
                
cf_equal(x, value = "flagged")

}
\seealso{
Other fossils: 
\code{\link{cf_age}()},
\code{\link{cf_outl}()},
\code{\link{cf_range}()}
}
\concept{fossils}
\keyword{Fossils}
\keyword{Temporal}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cc_cen.R
\name{cc_cen}
\alias{cc_cen}
\title{Identify Coordinates in Vicinity of Country and Province Centroids}
\usage{
cc_cen(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  species = "species",
  buffer = 1000,
  geod = TRUE,
  test = "both",
  ref = NULL,
  verify = FALSE,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing geographical coordinates and species
names.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{species}{character string. The column with the species identity. Only
required if verify = TRUE.}

\item{buffer}{numerical. The buffer around each province or country
centroid, where records should be flagged as problematic. Units depend on geod.  
Default = 1 kilometre.}

\item{geod}{logical. If TRUE the radius around each capital is calculated
based on a sphere, buffer is in meters and independent of latitude. If FALSE
the radius is calculated assuming planar coordinates and varies slightly with latitude,
in this case buffer is in degrees. Default = TRUE. See https://seethedatablog.wordpress.com/ 
for detail and credits.}

\item{test}{a character string. Specifying the details of the test. One of
c(\dQuote{both}, \dQuote{country}, \dQuote{provinces}).  If both tests for
country and province centroids.}

\item{ref}{SpatialPointsDataFrame. Providing the geographic gazetteer. Can
be any SpatialPointsDataFrame, but the structure must be identical to
\code{\link{countryref}}.  Default = \code{\link{countryref}}.}

\item{verify}{logical. If TRUE records are only flagged if they are the
only record in a given species flagged close to a given reference.
If FALSE, the distance is the only criterion}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records within a radius around the geographic centroids of political
countries and provinces. Poorly geo-referenced occurrence records in
biological databases are often erroneously geo-referenced to centroids.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

x <- data.frame(species = letters[1:10], 
                decimallongitude = runif(100, -180, 180), 
                decimallatitude = runif(100, -90,90))
                
cc_cen(x, geod = FALSE)

\dontrun{
#' cc_inst(x, value = "flagged", buffer = 50000) #geod = T
}

}
\seealso{
Other Coordinates: 
\code{\link{cc_cap}()},
\code{\link{cc_coun}()},
\code{\link{cc_dupl}()},
\code{\link{cc_equ}()},
\code{\link{cc_gbif}()},
\code{\link{cc_inst}()},
\code{\link{cc_iucn}()},
\code{\link{cc_outl}()},
\code{\link{cc_sea}()},
\code{\link{cc_urb}()},
\code{\link{cc_val}()},
\code{\link{cc_zero}()}
}
\concept{Coordinates}
\keyword{Coordinate}
\keyword{cleaning}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoordinateCleaner-package.R
\docType{package}
\name{CoordinateCleaner-package}
\alias{CoordinateCleaner-package}
\alias{CoordinateCleaner}
\title{CoordinateCleaner}
\description{
Automated Cleaning of Occurrence Records from Biological Collections
}
\details{
Automated flagging of common spatial and temporal errors in biological and paleontological collection data, for the use in conservation, ecology and paleontology. Includes automated tests to easily flag (and exclude) records assigned to country or province centroid, the open ocean, the headquarters of the Global Biodiversity Information Facility, urban areas or the location of biodiversity institutions (museums, zoos, botanical gardens, universities). Furthermore identifies per species outlier coordinates, zero coordinates, identical latitude/longitude and invalid coordinates. Also implements an algorithm to identify data sets with a significant proportion of rounded coordinates. Especially suited for large data sets. See <https://ropensci.github.io/CoordinateCleaner/> for more details and tutorials.
}
\author{
Alexander Zizka, Daniele Silvestro, Tobias Andermann, Josue Azevedo, 
Camila Duarte Ritter, Daniel Edler, Harith Farooq, Andrei Herdean, Maria Ariza, 
Ruud Scharn, Sten Svantesson, Niklas Wengstrom, Vera Zizka
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.spatialvalid.R
\name{plot.spatialvalid}
\alias{plot.spatialvalid}
\title{Plot Method for Class Spatialvalid}
\usage{
\method{plot}{spatialvalid}(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  bgmap = NULL,
  clean = TRUE,
  details = FALSE,
  pts_size = 1,
  font_size = 10,
  zoom_f = 0.1,
  ...
)
}
\arguments{
\item{x}{an object of the class \code{spatialvalid} as from
\code{\link{clean_coordinates}}.}

\item{lon}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the latitude coordinates.
Default = \dQuote{decimallatitude}.}

\item{bgmap}{an object of the class \code{SpatialPolygonsDataFrame} used as
background map. Default = ggplot::borders()}

\item{clean}{logical.  If TRUE, non-flagged coordinates are included in the
map.}

\item{details}{logical. If TRUE, occurrences are color-coded by the type of
flag.}

\item{pts_size}{numeric. The point size for the plot.}

\item{font_size}{numeric. The font size for the legend and axes}

\item{zoom_f}{numeric. the fraction by which to expand the plotting area 
from the occurrence records. Increase, if countries do not show 
up on the background map.}

\item{\dots}{arguments to be passed to methods.}
}
\value{
A plot of the records flagged as potentially erroneous by
\code{\link{clean_coordinates}}.
}
\description{
A set of plots to explore objects of the class \code{spatialvalid}. A plot
to visualize the flags from clean_coordinates
}
\examples{


exmpl <- data.frame(species = sample(letters, size = 250, replace = TRUE),
                    decimallongitude = runif(250, min = 42, max = 51),
                    decimallatitude = runif(250, min = -26, max = -11))

test <- clean_coordinates(exmpl, species = "species", 
                          tests = c("sea", "gbif", "zeros"),
                          verbose = FALSE)

summary(test)
plot(test)
}
\seealso{
\code{\link{clean_coordinates}}
}
\keyword{Visualisation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/CoordinateCleaner-package.R
\docType{data}
\name{institutions}
\alias{institutions}
\title{Global Locations of Biodiversity Institutions}
\format{
A data frame with 12170 observations on 12 variables.
}
\source{
Compiled from various sources: \itemize{ \item Global Biodiversity
Information Facility \url{https://www.gbif.org/} \item Wikipedia
\url{https://www.wikipedia.org/} \item Geonames \url{https://www.geonames.org/} \item The Global
Registry of Biodiversity Repositories \item Index
Herbariorum \url{http://sweetgum.nybg.org/science/ih/}
\item Botanic Gardens Conservation International \url{https://www.bgci.org/}
}
}
\description{
A global gazetteer for biodiversity institutions from various sources,
including zoos, museums, botanical gardens, GBIF contributors, herbaria,
university collections.
}
\examples{

data(institutions)
str(institutions)

}
\keyword{gazetteers}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cf_outl.R
\name{cf_outl}
\alias{cf_outl}
\title{Identify Outlier Records in Space and Time}
\usage{
cf_outl(
  x,
  lon = "decimallongitude",
  lat = "decimallatitude",
  min_age = "min_ma",
  max_age = "max_ma",
  taxon = "accepted_name",
  method = "quantile",
  size_thresh = 7,
  mltpl = 5,
  replicates = 5,
  flag_thresh = 0.5,
  uniq_loc = FALSE,
  value = "clean",
  verbose = TRUE
)
}
\arguments{
\item{x}{data.frame. Containing fossil records with taxon names, ages, 
and geographic coordinates.}

\item{lon}{character string. The column with the longitude coordinates.
To identify unique records if \code{uniq_loc  = TRUE}.
Default = \dQuote{decimallongitude}.}

\item{lat}{character string. The column with the longitude coordinates.
Default = \dQuote{decimallatitude}. To identify unique records if \code{uniq_loc  = T}.}

\item{min_age}{character string. The column with the minimum age. Default
= \dQuote{min_ma}.}

\item{max_age}{character string. The column with the maximum age. Default
= \dQuote{max_ma}.}

\item{taxon}{character string. The column with the taxon name. If
\dQuote{}, searches for outliers over the entire dataset, otherwise per
specified taxon. Default = \dQuote{accepted_name}.}

\item{method}{character string.  Defining the method for outlier
selection.  See details. Either \dQuote{quantile} or \dQuote{mad}.  Default
= \dQuote{quantile}.}

\item{size_thresh}{numeric.  The minimum number of records needed for a
dataset to be tested. Default = 10.}

\item{mltpl}{numeric. The multiplier of the interquartile range
(\code{method == 'quantile'}) or median absolute deviation (\code{method ==
'mad'}) to identify outliers. See details.  Default = 5.}

\item{replicates}{numeric. The number of replications for the distance
matrix calculation. See details.  Default = 5.}

\item{flag_thresh}{numeric.  The fraction of passed replicates necessary to pass the test. 
See details. Default = 0.5.}

\item{uniq_loc}{logical.  If TRUE only single records per location and time
point (and taxon if \code{taxon} != "") are used for the outlier testing.
Default = T.}

\item{value}{character string.  Defining the output value. See value.}

\item{verbose}{logical. If TRUE reports the name of the test and the number
of records flagged.}
}
\value{
Depending on the \sQuote{value} argument, either a \code{data.frame}
containing the records considered correct by the test (\dQuote{clean}) or a
logical vector (\dQuote{flagged}), with TRUE = test passed and FALSE = test failed/potentially
problematic . Default = \dQuote{clean}.
}
\description{
Removes or flags records of fossils that are spatio-temporal outliers based on
interquantile ranges. Records are flagged if they are either extreme in time
or space, or both.
}
\details{
The outlier detection is based on an interquantile range test. In a first
step a distance matrix of geographic distances among all records is
calculate. Subsequently a similar distance matrix of temporal distances
among all records is calculated based on a single point selected by random
between the minimum and maximum age for each record. The mean distance for
each point to all neighbours is calculated for both matrices and spatial and
temporal distances are scaled to the same range. The sum of these distanced
is then tested against the interquantile range and flagged as an outlier if
\eqn{x > IQR(x) + q_75 * mltpl}. The test is replicated \sQuote{replicates}
times, to account for temporal uncertainty. Records are flagged as outliers
if they are flagged by a fraction of more than \sQuote{flag.thres}
replicates. Only datasets/taxa comprising more than \sQuote{size_thresh}
records are tested. Note that geographic distances are calculated as
geospheric distances for datasets (or taxa) with fewer than 10,000 records
and approximated as Euclidean distances for datasets/taxa with 10,000 to
25,000 records. Datasets/taxa comprising more than 25,000 records are
skipped.
}
\note{
See \url{https://ropensci.github.io/CoordinateCleaner/} for more
details and tutorials.
}
\examples{

minages <- c(runif(n = 11, min = 10, max = 25), 62.5)
x <- data.frame(species = c(letters[1:10], rep("z", 2)),
                lng = c(runif(n = 10, min = 4, max = 16), 75, 7),
                lat = c(runif(n = 12, min = -5, max = 5)),
                min_ma = minages, 
                max_ma = c(minages[1:11] + runif(n = 11, min = 0, max = 5), 65))

cf_outl(x, value = "flagged", taxon = "")

}
\seealso{
Other fossils: 
\code{\link{cf_age}()},
\code{\link{cf_equal}()},
\code{\link{cf_range}()}
}
\concept{fossils}
\keyword{Coordinate}
\keyword{Fossil}
\keyword{Temporal}
\keyword{cleaning}
