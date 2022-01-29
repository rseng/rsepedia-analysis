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

The package has been removed from CRAN due to slow response time. Response time was slow because I am on parental leave.