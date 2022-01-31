
<!-- README.md is generated from README.Rmd. Please edit that file -->

# rgugik <img src="man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->

[![CRAN](http://www.r-pkg.org/badges/version/rgugik)](https://cran.r-project.org/package=rgugik)
[![R build
status](https://github.com/kadyb/rgugik/workflows/rcmdcheck/badge.svg)](https://github.com/kadyb/rgugik/actions)
[![codecov](https://codecov.io/gh/kadyb/rgugik/branch/master/graph/badge.svg)](https://app.codecov.io/gh/kadyb/rgugik)
[![License:
MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02948/status.svg)](https://doi.org/10.21105/joss.02948)
<!-- badges: end -->

**rgugik** is an R package for downloading open data from resources of
[Polish Head Office of Geodesy and
Cartography](https://www.gov.pl/web/gugik) including:

-   [Orthophotomaps](http://www.gugik.gov.pl/pzgik/zamow-dane/ortofotomapa)
-   [General Geographic
    Database](http://www.gugik.gov.pl/pzgik/zamow-dane/baza-danych-obiektow-ogolnogeograficznych)
-   [Topographic
    Database](http://www.gugik.gov.pl/pzgik/zamow-dane/baza-danych-obiektow-topograficznych-bdot-10k)
-   [Register of Towns, Streets and
    Addresses](https://emuia.gugik.gov.pl)
-   [State Register of Geographical
    Names](https://www.geoportal.gov.pl/dane/panstwowy-rejestr-nazw-geograficznych)
-   [State Register of
    Borders](http://www.gugik.gov.pl/pzgik/zamow-dane/panstwowy-rejestr-granic-i-powierzchni-jednostek-podzialow-terytorialnych-kraju)
-   Location (geometry) of cadastral parcels using TERYT (parcel ID) or
    coordinates
-   3D models of buildings (LOD1, LOD2)
-   Various digital elevation models as:
    -   [Digital terrain
        model](http://www.gugik.gov.pl/pzgik/zamow-dane/numeryczny-model-terenu)
    -   [Digital surface
        model](http://www.gugik.gov.pl/pzgik/zamow-dane/numeryczny-model-pokrycia-terenu)
    -   [Point
        cloud](http://www.gugik.gov.pl/pzgik/zamow-dane/dane-pomiarowe)

It is also possible to geocode addresses or objects using the
`geocodePL_get()` function.

**Corresponding functions**

| Function                              | Input                  | Dastaset EN                              | Dataset PL                                |
|:--------------------------------------|:-----------------------|:-----------------------------------------|:------------------------------------------|
| `ortho_request()`, `tile_download()`  | geometry               | Orthophotomap                            | Ortofotomapa                              |
| `geodb_download()`                    | voivodeship            | General Geographic Database              | Baza Danych Obiektów Ogólnogeograficznych |
| `topodb_download()`                   | county                 | Topographic Database                     | Baza Danych Obiektów Topograficznych      |
| `emuia_download()`                    | commune                | Register of Towns, Streets and Addresses | Ewidencja Miejscowości, Ulic i Adresów    |
| `geonames_download()`                 | type                   | State Register of Geographical Names     | Państwowy Rejestr Nazw Geograficznych     |
| `borders_get()`, `borders_download()` | type                   | State Register of Borders                | Państwowy Rejestr Granic                  |
| `parcel_get()`                        | parcel ID, coordinates | Location of cadastral parcels            | Lokalizacja działek katastralnych         |
| `models3D_download()`                 | county                 | 3D models of buildings                   | Modele 3D budynków                        |
| `DEM_request()`, `tile_download()`    | geometry               | Digital elevation models                 | Cyfrowe modele wysokościowe               |

There are the additional functions for obtaining digital terrain model:

-   `pointDTM_get()` for small areas (high resolution grid)
-   `pointDTM100_download()` for voivodeships areas (low resolution
    grid)
-   `minmaxDTM_get()` to find the minimum and maximum elevation (small
    areas)

The names of administrative units and their IDs can be obtained using
these functions:

-   `voivodeship_names()` (16)
-   `county_names()` (380)
-   `commune_names()` (2477)

## Installation

You can install the released version from
[CRAN](https://cran.r-project.org/) with:

``` r
install.packages("rgugik")
```

You can install the development version from
[GitHub](https://github.com) with:

``` r
# install.packages("remotes")
remotes::install_github("kadyb/rgugik")
```

## Usage

### Orthophotomap

-   `ortho_request()` - returns a data frame with metadata and links to
    the orthoimages for a given geometry (point, line or polygon)
-   `tile_download()` - downloads orthoimages based on the data frame
    obtained using the `ortho_request()` function

``` r
library(rgugik)
library(sf)
library(raster)

polygon_path = system.file("datasets/search_area.gpkg", package = "rgugik")
polygon = read_sf(polygon_path)

req_df = ortho_request(polygon)

# select the oldest image
req_df = req_df[req_df$year == 2001, ]

# print metadata
t(req_df)
#>             5                                                                               
#> sheetID     "N-33-130-D-b-2-3"                                                              
#> year        "2001"                                                                          
#> resolution  "1"                                                                             
#> composition "RGB"                                                                           
#> sensor      "Satellite"                                                                     
#> CRS         "PL-1992"                                                                       
#> isFilled    "TRUE"                                                                          
#> URL         "https://opendata.geoportal.gov.pl/ortofotomapa/41/41_3756_N-33-130-D-b-2-3.tif"
#> seriesID    "41"                                                                            
#> sha1        "312c81963a31e268fc20c442733c48e1aa33838f"                                      
#> date        "2001-01-01"                                                                    
#> filename    "41_3756_N-33-130-D-b-2-3"

# download image
tile_download(req_df)
#> 1/1

img = brick("41_3756_N-33-130-D-b-2-3.tif")
plotRGB(img)
```

<img src="man/figures/README-f1-1.png" width="100%" />

### Administrative boundaries

``` r
library(rgugik)
library(sf)

# get counties from opolskie voivodeship (TERYT 16)
counties = county_names
counties = counties[substr(counties$TERYT, 1, 2) == "16", "TERYT"]
counties_geom = borders_get(TERYT = counties)
plot(st_geometry(counties_geom), main = "Opolskie")
```

<img src="man/figures/README-f2-1.png" width="100%" />

### Vignettes

More advanced examples of the practical (step by step) use of this
package can be found in the vignettes:

-   [Orthophotomap](https://kadyb.github.io/rgugik/articles/orthophotomap.html)
-   [Digital elevation
    model](https://kadyb.github.io/rgugik/articles/DEM.html)
-   [Topographic
    Database](https://kadyb.github.io/rgugik/articles/topodb.html)

## Acknowledgment

[Head Office of Geodesy and Cartography in
Poland](https://www.gov.pl/web/gugik) is the main source of the provided
data. The data is made available in accordance with the [Act of May 17,
1989 Geodetic and Cartographic
Law](http://isap.sejm.gov.pl/isap.nsf/DocDetails.xsp?id=WDU19890300163)
(amended on 16 April 2020).

All datasets can be explored interactively using the
[Geoportal](https://mapy.geoportal.gov.pl).

## Contribution

Contributions to this package are welcome. The preferred method of
contribution is through a GitHub pull request. Feel also free to contact
us by creating [an issue](https://github.com/kadyb/rgugik/issues). More
detailed information can be found in the
[CONTRIBUTING](https://github.com/kadyb/rgugik/blob/master/CONTRIBUTING.md)
document.

Maintainers and contributors must follow this repository’s [CODE OF
CONDUCT](https://github.com/kadyb/rgugik/blob/master/CODE_OF_CONDUCT.md).

## Citation

To cite **rgugik** in publications, please use the following
[article](https://doi.org/10.21105/joss.02948):

    Dyba, K. and Nowosad, J. (2021). rgugik: Search and Retrieve Spatial Data from the Polish Head Office of Geodesy and Cartography in R. Journal of Open Source Software, 6(59), 2948, https://doi.org/10.21105/joss.02948

BibTeX version can be obtained with `citation("rgugik")`.

## Related projects

If you don’t feel familiar with R, there is a similar
[QGIS](https://www.qgis.org/en/site/) tool in the
[EnviroSolutions](https://github.com/envirosolutionspl) repository.

# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who contribute through reporting issues, posting feature requests, updating documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for everyone, regardless of level of experience, gender, gender identity and expression, sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or imagery, derogatory comments or personal attacks, trolling, public or private harassment, insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments, commits, code, wiki edits, issues, and other contributions that are not aligned to this Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the [Contributor Covenant](http:contributor-covenant.org), version 1.0.0, available at https://www.contributor-covenant.org/version/1/0/0/code-of-conduct.html


# CONTRIBUTING

No matter your current skills, it’s possible to contribute to the `rgugik` package.
We appreciate any contribution no matter the amount.

## Bugs

If you’ve found a bug, please create a minimal reproducible example using the [reprex](https://www.tidyverse.org/help#reprex) package first.
Spend some time trying to make it as minimal as possible, this will facilitate the task and speed up the entire process.
Next, submit an issue on the [Issues page](https://github.com/kadyb/rgugik/issues).

## Contributions

### Fixing typos

You can fix typos, spelling mistakes, or grammatical errors in the documentation.
We use [roxygen2](https://roxygen2.r-lib.org/), so the documentation should be generated using `.R` files, not by editing the `.Rd` files directly.

### Greater changes

If you want to make a greater change, it's a good idea to file an issue first and make sure someone from the team agrees that it’s needed.
We don’t want you to spend a bunch of time on something that we don’t think is a suitable idea.

Once accepted, you can follow the pull request process:

1. Fork this repo to your GitHub account.
2. Clone your version to your machine, e.g., `git clone https://github.com/kadyb/rgugik.git`.
3. Make sure to track progress upstream (i.e., our version of `rgugik` at `kadyb/rgugik`) by doing `git remote add upstream https://github.com/kadyb/rgugik.git`. 
Before making any changes, make sure to pull changes in from upstream by either doing `git fetch upstream` then merge later, or `git pull upstream` to fetch and merge in one step.
4. Make your changes (make changes to a new branch).
5. If you alter package functionality at all (e.g., the code itself, not just documentation) please do write some tests to cover the new functionality.
6. Push changes to your account.
7. Submit a pull request to the master branch at `kadyb/rgugik`.

We use [testthat](https://testthat.r-lib.org/) for unit tests. Contributions with test cases included are prioritized to accept.

Please make sure that your new code and documentation match the existing style.
We use [lintr](https://github.com/jimhester/lintr) for static code analysis (i.e., code style).

## Questions

Questions are welcomed on the [Issues page](https://github.com/kadyb/rgugik/issues).
Adding a reproducible example may make it easier for us to answer.

## Thanks for contributing!
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r include=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```

# rgugik <img src="man/figures/logo.png" align="right" width="150"/>

<!-- badges: start -->
[![CRAN](http://www.r-pkg.org/badges/version/rgugik)](https://cran.r-project.org/package=rgugik)
[![R build status](https://github.com/kadyb/rgugik/workflows/rcmdcheck/badge.svg)](https://github.com/kadyb/rgugik/actions)
[![codecov](https://codecov.io/gh/kadyb/rgugik/branch/master/graph/badge.svg)](https://app.codecov.io/gh/kadyb/rgugik)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02948/status.svg)](https://doi.org/10.21105/joss.02948)
<!-- badges: end -->

**rgugik** is an R package for downloading open data from resources of [Polish Head Office of Geodesy and Cartography](https://www.gov.pl/web/gugik) including:

  - [Orthophotomaps](http://www.gugik.gov.pl/pzgik/zamow-dane/ortofotomapa)
  - [General Geographic Database](http://www.gugik.gov.pl/pzgik/zamow-dane/baza-danych-obiektow-ogolnogeograficznych)
  - [Topographic Database](http://www.gugik.gov.pl/pzgik/zamow-dane/baza-danych-obiektow-topograficznych-bdot-10k)
  - [Register of Towns, Streets and Addresses](https://emuia.gugik.gov.pl)
  - [State Register of Geographical Names](https://www.geoportal.gov.pl/dane/panstwowy-rejestr-nazw-geograficznych)
  - [State Register of Borders](http://www.gugik.gov.pl/pzgik/zamow-dane/panstwowy-rejestr-granic-i-powierzchni-jednostek-podzialow-terytorialnych-kraju)
  - Location (geometry) of cadastral parcels using TERYT (parcel ID) or coordinates
  - 3D models of buildings (LOD1, LOD2)
  - Various digital elevation models as:
    - [Digital terrain model](http://www.gugik.gov.pl/pzgik/zamow-dane/numeryczny-model-terenu)
    - [Digital surface model](http://www.gugik.gov.pl/pzgik/zamow-dane/numeryczny-model-pokrycia-terenu)
    - [Point cloud](http://www.gugik.gov.pl/pzgik/zamow-dane/dane-pomiarowe)

It is also possible to geocode addresses or objects using the `geocodePL_get()` function.

**Corresponding functions**

```{r echo=FALSE}
ds_pl = c("Ortofotomapa",
          "Baza Danych Obiektów Ogólnogeograficznych",
          "Baza Danych Obiektów Topograficznych",
          "Ewidencja Miejscowości, Ulic i Adresów",
          "Państwowy Rejestr Nazw Geograficznych",
          "Państwowy Rejestr Granic",
          "Lokalizacja działek katastralnych",
          "Modele 3D budynków",
          "Cyfrowe modele wysokościowe")

ds_en = c("Orthophotomap",
          "General Geographic Database",
          "Topographic Database",
          "Register of Towns, Streets and Addresses",
          "State Register of Geographical Names",
          "State Register of Borders",
          "Location of cadastral parcels",
          "3D models of buildings",
          "Digital elevation models")

fun = c("`ortho_request()`, `tile_download()`",
        "`geodb_download()`",
        "`topodb_download()`",
        "`emuia_download()`",
        "`geonames_download()`",
        "`borders_get()`, `borders_download()`",
        "`parcel_get()`",
        "`models3D_download()`",
        "`DEM_request()`, `tile_download()`")

input = c("geometry",
          "voivodeship",
          "county",
          "commune",
          "type",
          "type",
          "parcel ID, coordinates",
          "county",
          "geometry")

df = data.frame(fun, input, ds_en, ds_pl)
colnames(df) = c("Function", "Input", "Dastaset EN", "Dataset PL")

knitr::kable(df)
```

There are the additional functions for obtaining digital terrain model:

  - `pointDTM_get()` for small areas (high resolution grid)
  - `pointDTM100_download()` for voivodeships areas (low resolution grid)
  - `minmaxDTM_get()` to find the minimum and maximum elevation (small areas)

The names of administrative units and their IDs can be obtained using these functions:

  - `voivodeship_names()` (16)
  - `county_names()` (380)
  - `commune_names()` (2477)

## Installation

You can install the released version from [CRAN](https://cran.r-project.org/) with:

```{r eval=FALSE}
install.packages("rgugik")
```

You can install the development version from [GitHub](https://github.com) with:

```{r message=FALSE, warning=FALSE, eval=FALSE}
# install.packages("remotes")
remotes::install_github("kadyb/rgugik")
```

## Usage

### Orthophotomap
  
- `ortho_request()` - returns a data frame with metadata and links to the orthoimages for a given geometry (point, line or polygon)
- `tile_download()` - downloads orthoimages based on the data frame obtained using the `ortho_request()` function
  
```{r f1, message=FALSE, warning=FALSE}
library(rgugik)
library(sf)
library(raster)

polygon_path = system.file("datasets/search_area.gpkg", package = "rgugik")
polygon = read_sf(polygon_path)

req_df = ortho_request(polygon)

# select the oldest image
req_df = req_df[req_df$year == 2001, ]

# print metadata
t(req_df)

# download image
tile_download(req_df)

img = brick("41_3756_N-33-130-D-b-2-3.tif")
plotRGB(img)
```

```{r echo=FALSE, message=FALSE}
invisible(file.remove("41_3756_N-33-130-D-b-2-3.tif"))
```

### Administrative boundaries

```{r f2}
library(rgugik)
library(sf)

# get counties from opolskie voivodeship (TERYT 16)
counties = county_names
counties = counties[substr(counties$TERYT, 1, 2) == "16", "TERYT"]
counties_geom = borders_get(TERYT = counties)
plot(st_geometry(counties_geom), main = "Opolskie")
```

### Vignettes

More advanced examples of the practical (step by step) use of this package can be found in the vignettes:

- [Orthophotomap](https://kadyb.github.io/rgugik/articles/orthophotomap.html)
- [Digital elevation model](https://kadyb.github.io/rgugik/articles/DEM.html)
- [Topographic Database](https://kadyb.github.io/rgugik/articles/topodb.html)

## Acknowledgment

[Head Office of Geodesy and Cartography in Poland](https://www.gov.pl/web/gugik) is the main source of the provided data. The data is made available in accordance with the [Act of May 17, 1989 Geodetic and Cartographic Law](http://isap.sejm.gov.pl/isap.nsf/DocDetails.xsp?id=WDU19890300163) (amended on 16 April 2020).

All datasets can be explored interactively using the [Geoportal](https://mapy.geoportal.gov.pl).

## Contribution

Contributions to this package are welcome. 
The preferred method of contribution is through a GitHub pull request. 
Feel also free to contact us by creating [an issue](https://github.com/kadyb/rgugik/issues).
More detailed information can be found in the [CONTRIBUTING](https://github.com/kadyb/rgugik/blob/master/CONTRIBUTING.md) document.

Maintainers and contributors must follow this repository’s [CODE OF CONDUCT](https://github.com/kadyb/rgugik/blob/master/CODE_OF_CONDUCT.md).

## Citation

To cite **rgugik** in publications, please use the following [article](https://doi.org/10.21105/joss.02948):

```
Dyba, K. and Nowosad, J. (2021). rgugik: Search and Retrieve Spatial Data from the Polish Head Office of Geodesy and Cartography in R. Journal of Open Source Software, 6(59), 2948, https://doi.org/10.21105/joss.02948
```

BibTeX version can be obtained with `citation("rgugik")`.

## Related projects

If you don't feel familiar with R, there is a similar [QGIS](https://www.qgis.org/en/site/) tool in the [EnviroSolutions](https://github.com/envirosolutionspl) repository.
---
title: "Digital elevation model"
author: "Krzysztof Dyba, Jakub Nowosad"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Digital elevation model}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



<style>
body {text-align: justify}
</style>

## Definition

**Digital elevation models** are models describing the terrain surface.
They are created as a result of the processing of aerial photos, laser scanning (LiDAR), geodetic surveying, or radar interferometry (InSAR).
DEMs are one of the key datasets in Geographic Information Systems (GIS) and constitute the basis for many environmental spatial analyses.
In addition, they are a source for derived products such as terrain slope and aspect.
DEM is the general name for a group of models with different characteristics, including:

1. **Digital terrain model** (DTM) - numerical representation of the terrain with its morphological forms.
This model is devoid of any objects above the ground, such as buildings or trees.
2. **Digital surface model** (DSM) - numerical representation of the terrain, including objects on its surface, such as buildings and trees.

<center>
![DTM-DSM](https://upload.wikimedia.org/wikipedia/commons/thumb/6/6c/DTM_DSM.svg/1024px-DTM_DSM.svg.png){ width=50% }
<p style="font-size:8px;">
Wikimedia Commons, the free media repository, https://commons.wikimedia.org/w/index.php?title=File:DTM_DSM.svg&oldid=475779479
(accessed October 7, 2020).
</p>
</center>
<br>

The properties of the DEMs:

1. Format - we can distinguish three main structures: **GRID** (point / cell), **TIN** (irregular topological triangle network) and **contour lines** (vector).
Currently, the most used format is GRID.
2. Accuracy - is related to the vertical measurement error.
3. Spatial resolution - is related to the size of the smallest object that can be detected by the sensor and is determined by the size of the image cell (pixel).
The larger the cell, the more generalized form of the terrain it presents.

## Purpose

The purpose of this vignette is to determine the elevation of the ground surface and objects in the selected area.
The source of the data will be Airborne Laser Scanning already processed to the GRID format.

## Analysis


```r
# attach packages
library(sf)
library(stars)
library(rgugik)
```

Our analysis area is the Morasko Meteorite nature reserve located in the Greater Poland voivodeship.
It was established in 1976 in order to protect the area of impact craters, which, according to researchers, were formed in the fall of the Morasko meteorite about 5,000 years ago.
In addition, the oak-hornbeam forest with rare species of plants (*lilium martagon*, *ceratophyllum submersum*) and birds (*european nightjar*, *black woodpecker*) is protected.

### Vector data

The centroid (geometric center) of the Morasko Meteorite nature reserve has X = 16.895 and Y = 52.487 coordinates in World Geodetic System 1984 (EPSG 4326).
Let's start by creating this point with the **sf** package.


```r
morasko = st_point(c(16.895, 52.489))
morasko = st_sfc(morasko, crs = 4326) # set coordinate system
morasko
```

```
## Geometry set for 1 feature 
## Geometry type: POINT
## Dimension:     XY
## Bounding box:  xmin: 16.895 ymin: 52.489 xmax: 16.895 ymax: 52.489
## Geodetic CRS:  WGS 84
```

```
## POINT (16.895 52.489)
```

Now our point is embedded in space (has a spatial reference).
In the next step, let's create an approximate zone that will include the area of the reserve.
The function `st_buffer()` will be used for this.
Before this operation, we need to transform the coordinate system to a system with metric units, e.g. Poland CS92 (EPSG 2180), using `st_transform()` function.


```r
morasko = st_transform(morasko, crs = 2180)
morasko_buffer = st_buffer(morasko, dist = 400)
```

We have created a buffer with a radius of 400 meters.
Let's visualize it.


```r
plot(morasko_buffer, axes = TRUE, main = "Morasko reserve buffer")
plot(morasko, add = TRUE)
```

![](DEM-GEOMETRY-1.png)

Of course, the area shown above is not exactly the reserve area.
The exact area can be determined from the polygon layer as in [orthophotomap](https://kadyb.github.io/rgugik/articles/orthophotomap.html) example using the General Geographic Database.

### Raster data

Now we can search for available elevation data for this area using the `DEM_request()` function (it is analogous to the `ortho_request()` function).
The only argument of the function is our reserve buffer.


```r
req_df = DEM_request(morasko_buffer)
```

Let's check the obtained results.


```r
# display the first 10 rows and the first 5 columns
req_df[1:10, 1:5]
```

```
##             sheetID year              format resolution avgElevErr
## 1    N-33-130-D-b-1 2007      Intergraph TTN       <NA>        1.5
## 2  N-33-130-D-b-1-1 2020 ARC/INFO ASCII GRID      5.0 m        0.5
## 3  N-33-130-D-b-1-1 2017      ASCII XYZ GRID      1.0 m        0.9
## 4  N-33-130-D-b-1-1 2012 ARC/INFO ASCII GRID      0.5 m        0.1
## 5     6.179.11.13.4 2018           ASCII TBD      1.0 m        0.1
## 6     6.179.11.14.1 2018           ASCII TBD      1.0 m        0.1
## 7     6.179.11.08.4 2018           ASCII TBD      1.0 m        0.1
## 8     6.179.11.14.3 2018           ASCII TBD      1.0 m        0.1
## 9     6.179.11.09.3 2018           ASCII TBD      1.0 m        0.1
## 10    6.179.11.13.2 2018           ASCII TBD      1.0 m        0.1
```

We have received metadata with many types of data of different formats, timeliness, resolution, and accuracy.
For our analysis, we need digital terrain model (DTM) and digital surface model (DSM) in the "ARC/INFO ASCII GRID" format.
Let's make data selection by creating two tables and combining them together.


```r
req_df_DTM = req_df[req_df$format == "ARC/INFO ASCII GRID" &
                    req_df$product == "DTM" &
                    req_df$year == 2019, ]
req_df_DSM = req_df[req_df$format == "ARC/INFO ASCII GRID" &
                    req_df$product == "DSM" &
                    req_df$year == 2019, ]

# combine data tables
req_df = rbind(req_df_DTM, req_df_DSM)
req_df[, 1:5]
```

```
##             sheetID year              format resolution avgElevErr
## 39 N-33-130-D-b-1-1 2019 ARC/INFO ASCII GRID      1.0 m        0.1
## 11 N-33-130-D-b-1-1 2019 ARC/INFO ASCII GRID      0.5 m        0.1
```

Now we can download the data using the `tile_download()` function with our filtered data frame as input.


```r
# 168.7 MB
tile_download(req_df, outdir = "./data")
```

```
## 1/2
## 2/2
```

If you run into any problem with the download, remember that you can pass another download method from `download.file()` as a function argument.


```r
tile_download(req_df, outdir = "./data", method = "wget")
```

### Processing

Let’s load the downloaded numerical models using the `read_stars()` function from the **stars** package, which allows working on spatiotemporal arrays.
We have two files, one represents DTM and second represents DSM.


```r
# load data
DTM = read_stars("data/73044_917579_N-33-130-D-b-1-1.asc", proxy = FALSE)
DSM = read_stars("data/73043_917495_N-33-130-D-b-1-1.asc", proxy = FALSE)

# name raster
names(DTM) = "DTM"
names(DSM) = "DSM"

# set coordinate system
st_crs(DTM) = 2180
st_crs(DSM) = 2180
```

You probably noticed the four-fold difference in their sizes.
It is due to the difference between their cells resolutions.
We need to unify them to a common resolution to be able to combine them into one stack.
It is much better to use a lower resolution than to increase it, because we cannot get more information and the processing will be faster.
Let's use the `st_warp()` function to do this.


```r
DSM = st_warp(DSM, dest = DTM, cellsize = 1)
```

Now, both models have the same dimensions (the number of rows and columns) and spatial resolution.
Thus, we can combine them into one object (`DEM`).


```r
DEM = c(DTM, DSM)
length(DEM)
```

```
## [1] 2
```



Now we have a DEM object that consists of two attributes (DTM and DSM).
In fact, both attributes contains same type of data as they are representing elevation.
Therefore, we can collapse the attributes into a new dimension.
Let's do that using `st_redimension()`.


```r
DEM = st_redimension(DEM)
names(st_dimensions(DEM))[3] = "elev" # name new data dim
st_dimensions(DEM)
```

```
##      from   to offset delta               refsys point   values x/y
## x       1 2188 355733     1 ETRS89 / Poland CS92  TRUE     NULL [x]
## y       1 2379 517029    -1 ETRS89 / Poland CS92  TRUE     NULL [y]
## elev    1    2     NA    NA                   NA    NA DTM, DSM
```

After this operation, our elevation attribute consists of the DTM and DSM layers (dimensions).
Then let's crop the rasters to our buffer.


```r
DEM = st_crop(DEM, morasko_buffer)
```

Let's check what the result looks like.


```r
plot(DEM, col = terrain.colors(99, alpha = NULL))
```

![](DEM-PLOT-1.png)

In the first quadrant of the circle, we can see five smaller circles.
These are the craters formed after the impact of the Morasko meteorite.
The largest fragment found weighs 272 kg and it is the largest meteorite found in Poland.
The collection of found meteorites can be seen at the [Earth Museum](https://muzeumziemi.amu.edu.pl/) in Poznań.

Let's calculate the crater width using the terrain transverse profile.
We can use our centroid and add a second example point 30 degrees towards N.
Next, we connect these points into a line (`st_linestring()`) and then sample this line every 1 m (`st_line_sample()`), because our DEM has this resolution.
As a result, we get one complex geometry (*MULTIPOINT*), which we have to convert into a simple geometry (*POINT*) consisting of many points.
The function `st_cast()` is used for this.


```r
pts_matrix = matrix(c(357121.7, 515765.5,
                      357321.2, 516017.9),
                    ncol = 2, byrow = TRUE)
line = st_sfc(st_linestring(pts_matrix), crs = 2180)
line = st_line_sample(line, density = 1)
line = st_cast(line, "POINT")
```


```r
# plot DTM (first layer)
plot(DEM[, , , 1], main = "DTM [m]", col = terrain.colors(99, alpha = NULL),
     reset = FALSE)
plot(line, col = "red", add = TRUE)
```

![](DEM-LINE-1.png)

In the last step, we extract the elevation values for these points using `st_extract()`.


```r
# take elevation from DTM and DSM layers
elev_line = st_extract(DEM, line)[[1]]
colnames(elev_line) = c("DTM", "DSM")
```

Now we can see how our transverse profile looks like.


```r
# use 'dev.off()' to reset previous plot
plot(elev_line[, "DTM"], type = "l", main = "Digital terrain model",
     ylab = "Elevation [m]", xlab = "Distance [m]", col = "red")
abline(v = c(126, 219), col = "blue")
```

![](DEM-PROFILE-1.png)

The largest width of the impact crater is about 90 m.

Okay, we checked the terrain.
In the last step, let's examine the height of the objects on it.
For this purpose, we calculate the height of the trees by subtracting the DTM from the DSM.
The product of this difference is called normalized DSM, because it takes the terrain elevation as a reference.


```r
calc = function(DEM) (DEM[2] - DEM[1])
nDSM = st_apply(DEM, MARGIN = c("x", "y"), FUN = calc)
plot(nDSM, main = "Trees height [m]",
     col = hcl.colors(9, palette = "Greens", rev = TRUE))
```

![](DEM-TREES-1.png)


---
title: "Orthophotomap"
author: "Krzysztof Dyba, Jakub Nowosad"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Orthophotomap}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



<style>
body {text-align: justify}
</style>

## Definition

**Orthophotomap** is a raster, orthogonal and cartometric representation of the terrain surface created by digital processing of aerial or satellite images.
During the orthorectification, geometric distortions resulting from the land relief are removed by using digital elevation models (DEM).
An orthophotomap is georeferenced, and therefore, allows to determine geographic coordinates for each of its cells.

Orthophotomaps' properties:

1. **Spatial resolution** - is related to the size of the smallest object that can be detected by the sensor and is determined by the size of the image cell (pixel).
The smaller the cell, the more detail it represents.
Too large a pixel means that individual objects in the scene are no longer recognizable.
2. **Composition** - analog images are in shades of gray, while digital images can be in natural colors (RGB) or near infrared (NIR).

## Purpose

The purpose of this vignette is to assess the vegetation condition of the selected area.
It can be done based on remote sensing data (multispectral orthophotomap) and a simple vegetation index.

**NDVI** (Normalized Difference Vegetation Index) is a simple indicator of vegetation that uses the red and the near infrared bands.
Its main application is monitoring and forecasting of agricultural production.
It is calculated using the following formula:

$$NDVI = \frac {NIR - RED} {NIR + RED}$$

Its value ranges from -1 to 1.
The higher the value, the higher the biomass level.
Values close to 0 and below are related to water, bare soil surfaces or buildings.

## Analysis


```r
# attach packages
library(sf)
library(stars)
library(rgugik)
```

The analysis area is the Krajkowo nature reserve located in the Greater Poland voivodeship.
It was established in 1958 in order to protect the breeding places of birds, especially the *grey heron* and the *great black cormorant*, and to protect the landscape of the Warta oxbow.

### Vector data

Data on nature reserves can be found in General Geographic Databases.
We can obtain them using the `geodb_download()` function.
Let's do that.


```r
# 17.6 MB
geodb_download("wielkopolskie", outdir = "./data")
```

If you run into any problem with the download, remember that you can pass another download method from `download.file()` as a function argument.


```r
geodb_download(req_df, outdir = "./data", method = "wget")
```

The downloaded database consists of many files in the GML (*Geography Markup Language*) format.
A brief description of the structure of this database can be found [here](https://kadyb.github.io/rgugik/articles/articles/spatialdb_description.html).
The table with the nature reserves is in the "PL.PZGIK.201.30__OT_TCRZ_A.xml" file.
We can use the **sf** package and its `read_sf()` function to load it.


```r
reserves = read_sf("data/PL.PZGiK.201.30/BDOO/PL.PZGIK.201.30__OT_TCRZ_A.xml")
```

Let's check the structure of our data.


```r
ncol(reserves)
## [1] 28
nrow(reserves)
## [1] 110
```

In simple terms, it is a spatial table consisting of 110 observations (rows) and 28 variables (columns).
The names of the objects are located in the **nazwa** column, which allow us to select the Krajkowo reserve only.


```r
# selection by attribute
krajkowo = reserves[reserves$nazwa == "Krajkowo", ]
```



We can display it in two basic ways:

1. Using the `plot()` function and directly specifying the column with the object geometry: `plot(krajkowo$geometry)`
2. Using the `plot()` and `st_geometry()` functions that obtain geometry from the vector layer.
In the first case, we need to know the name of the column with geometries (e.g. `geometry`, `geom`, etc.), while in the second case, the geometry is selected automatically (it is the safer and preferable way).


```r
plot(st_geometry(krajkowo), axes = TRUE, main = "Krajkowo reserve")
```

![](orto-GEOMETRY-1.png)

We can also calculate the area of this polygon.


```r
krajkowo_area = st_area(krajkowo) # [m^2]
units::set_units(krajkowo_area, "ha") # convert to [ha]
```

```
## 165.2744 [ha]
```

The `st_area()` function returned the area in m^2, and after the conversion we got the result of 165 ha.

### Raster data

Now let's move on to the stage of downloading the orthophotomap.
We use the `ortho_request()` function that show us which images are available for the analyzed area.
We need to provide our Krajkowo polygon as the argument of this function.


```r
req_df = ortho_request(krajkowo)
```

We can display the resulting table using the code below.


```r
# display the first 10 rows and the first 6 columns
req_df[1:10, 1:6]
```

```
##             sheetID year resolution composition  sensor     CRS
## 1  N-33-142-B-d-4-2 2004       0.50         B/W  Analog PL-1992
## 2  N-33-142-B-d-4-4 2004       0.50         B/W  Analog PL-1992
## 3  N-33-142-B-d-4-4 2010       0.25         RGB Digital PL-1992
## 4  N-33-142-B-d-4-2 2010       0.25         RGB Digital PL-1992
## 5  N-33-142-B-d-4-4 2010       0.25         CIR Digital PL-1992
## 6  N-33-142-B-d-4-2 2010       0.25         CIR Digital PL-1992
## 7  N-33-142-B-d-4-2 2016       0.25         CIR Digital PL-1992
## 8  N-33-142-B-d-4-2 2016       0.25         RGB Digital PL-1992
## 9  N-33-142-B-d-4-4 2016       0.25         CIR Digital PL-1992
## 10 N-33-142-B-d-4-4 2016       0.25         RGB Digital PL-1992
```

To complete our task, we need to obtain a near infrared data.
So in the next step, we select those rows for which the `composition` column has the value of "CIR".


```r
# select IR images and overwrite the req_df object
req_df = req_df[req_df$composition == "CIR", ]
```

Then let's sort the table according to the year the photo was taken, with the most recent images at the top of the table.


```r
req_df = req_df[order(-req_df$year), ]
```

Let's display the table again and select the newest compositions.


```r
req_df[, c(1:5, 9)]
```

```
##             sheetID year resolution composition  sensor seriesID
## 7  N-33-142-B-d-4-2 2016       0.25         CIR Digital    69837
## 9  N-33-142-B-d-4-4 2016       0.25         CIR Digital    69837
## 22 N-33-142-B-d-4-2 2013       0.25         CIR Digital    69903
## 24 N-33-142-B-d-4-4 2013       0.25         CIR Digital    69903
## 5  N-33-142-B-d-4-4 2010       0.25         CIR Digital    69763
## 6  N-33-142-B-d-4-2 2010       0.25         CIR Digital    69763
```


```r
req_df = req_df[req_df$year == 2016, ]
```

Note that the result has a pair of objects (images).
This means that our Krajkowo reserve is depicted in two photos within one series.
Therefore, the `seriesID` column is used to combine smaller images into a larger mosaic.


```r
req_df[, c(1:5, 9)]
```

```
##            sheetID year resolution composition  sensor seriesID
## 7 N-33-142-B-d-4-2 2016       0.25         CIR Digital    69837
## 9 N-33-142-B-d-4-4 2016       0.25         CIR Digital    69837
```

The `tile_download()` function is used to download orthophotomaps by taking our selected table as the main argument.
We can also specify the output folder with the `outdir` argument.


```r
# 61.9 MB
tile_download(req_df, outdir = "./data")
```

```
## 1/2
## 2/2
```

### Processing

Let's load the downloaded orthophotomaps using the `read_stars()` function from the **stars** package, which allows working with spatiotemporal arrays.
In our case, it is a raster consisting of three bands (NIR, R, G) for only one point in time.
We don't need to load the entire data array into memory - we can read the file's metadata instead by using the `proxy` argument.


```r
img1 = read_stars("data/69837_329609_N-33-142-B-d-4-2.TIF", proxy = TRUE)
img2 = read_stars("data/69837_329613_N-33-142-B-d-4-4.TIF", proxy = TRUE)
```

Now we can perform two operations: rasters merging and cropping to the reserve area.
The use of a `proxy` allows to get the result almost immediately, while processing the entire image (`proxy = FALSE`) would take several minutes.
Images have their own specific Coordinate Reference Systems, so let's make sure it is the correct one after merging.
It should be EPSG 2180 in this case.


```r
img = st_mosaic(img1, img2)
st_crs(img) = 2180 # overwrite CRS to be sure
img = st_crop(img, krajkowo)
```



Let's display the effect using the `plot()` function, and define the input bands with the `rgb` argument.
It creates a composition consisting of three bands: NIR, R and G in our case.
The composition below is shown in infrared, not in natural colors, which may be misinterpreted from the `rgb` argument name.


```r
plot(img, rgb = c(1, 2, 3), main = NULL)
```

![](orto-CIR-1.png)

In the last step, we calculate the NDVI using the near infrared (1) and red (2) bands.


```r
calc_ndvi = function(img) (img[1] - img[2]) / (img[1] + img[2])
ndvi = st_apply(img, MARGIN = c("x", "y"), FUN = calc_ndvi)
plot(ndvi, main = "NDVI", col = hcl.colors(10, palette = "RdYlGn"))
```

![](orto-NDVI-1.png)

A surprising observation is the relatively low NDVI values for the forest area.
There are two reasons for this, i.e. the photos are taken in mid-March (before the start of the growing season) and probably have not been calibrated.
For this reason, a better source of data for analysis may be satellite images, which are calibrated spectrally and obtained continuously (if no cloudiness occurs).


---
title: "Topographic Database"
author: "Krzysztof Dyba, Jakub Nowosad"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Topographic Database}
  %\VignetteEngine{knitr::rmarkdown}
  \usepackage[utf8]{inputenc}
---



<style>
body {text-align: justify}
</style>

## Definition

**Topographic Database** (pl. *Baza Danych Obiektów Topograficznych*) is a vector (object) database containing the spatial location of topographic objects with their characteristics for Poland.
The content and detail of the database correspond to the topographic map in the scale 1:10000.
The thematic scope includes information on water network, communication network, land cover, buildings and technical structures, utility infrastructure, land use, protected areas, territorial division units, and other objects.
The database is available in the Geography Markup Language (GML) format.
The source of its data comes from:

- other spatial databases maintained by Polish Head Office of Geodesy and Cartography (e.g. *Register of Towns, Streets and Addresses*, *Register of Land and Buildings*, *State Register of Borders*),
- databases maintained by other ministries or institutions (e.g., Ministry of Infrastructure, State Water Management, General Directorate for Environmental Protection),
- fieldworks.

## Purpose

The purpose of this vignette is to perform spatial operations on vector data from the Topographic Database.
We focus on four cases, taking into account different types of geometry, i.e. point, line, and polygon, and their attributes.
Also, we show how they can be visualized.

## Analysis


```r
# attach packages
library(sf)
library(rgugik)
```

Our analysis area is the *bieszczadzki* county located in the Subcarpathian (*podkarpackie*) voivodeship.
It is the farthest south area in Poland that also has the lowest population density (19 people on km^2^).

### Database

We start by downloading the topographic database for our county using the `topodb_download()` function.


```r
# 22.4 MB
topodb_download("bieszczadzki", outdir = "./data")
```

If you run into any problem with the download, remember that you can pass another download method from `download.file()` as a function argument.


```r
topodb_download(req_df, outdir = "./data", method = "wget")
```

The downloaded database consists of many files in the *GML* format.
All the data necessary for analyzes can be found in the `data/PL.PZGiK.332.1801/BDOT10k/` location.
A brief description of the structure of this database can be found [here](https://kadyb.github.io/rgugik/articles/articles/spatialdb_description.html).

First, let's load the file with administrative units ("PL.PZGiK.332.1801__OT_ADJA_A.xml") using the **sf** package and its `read_sf()` function.


```r
territory = read_sf("data/PL.PZGiK.332.1801/BDOT10k/PL.PZGiK.332.1801__OT_ADJA_A.xml")
```

The file contains administrative units at various levels, let's choose the lowest level, i.e. communes.
There are three types of communes in this dataset, specified in the `rodzaj` column: urban (*GM*), rural (*GW*), and urban-rural (*Gmw*).
Let's select them.


```r
communes = territory[territory$rodzaj %in% c("GM", "GW", "Gmw"), "rodzaj"]
table(communes$rodzaj)
```

```
## 
## Gmw  GW 
##   1   2
```

We can see that the *bieszczadzki* county consists of two rural communes and one urban-rural commune.
Let's visualize it.


```r
plot(communes, axes = TRUE, main = "Bieszczadzki county")
```

![](topodb-communes-1.png)

### Lengths and categories of roads

In the first task, we calculate the lengths of roads, taking into account their categories.
Road data can be found in the "PL.PZGiK.332.1801__OT_SKDR_L.xml" file.


```r
roads = read_sf("data/PL.PZGiK.332.1801/BDOT10k/PL.PZGiK.332.1801__OT_SKDR_L.xml")
```

Let's plot them.
We use the `plot()` function again, but this time we combine the two layers into one image.
The first layer (in the background) contains roads and to add another layer, we have to set the argument `reset` to `"FALSE"`.
Then we can add a second layer with the territory borders by setting the `add` argument to `"TRUE"`.


```r
plot(roads["katZarzadzania"], main = "Road category", reset = FALSE)
plot(st_geometry(territory), add = TRUE)
```

![](topodb-roads-1.png)

We have six road categories related to the managing entity.
Those are: national (K), voivodeship (W), county (P), communal (G), institutional (Z), and other (I) roads.

We use the `st_length()` function to find the length of each object in the table.
Next, we create a data frame consisting of the road category and its length.
Then we aggregate this data frame and calculate the sum of the lengths for each category.


```r
length_roads = st_length(roads)
length_roads = data.frame(length = length_roads,
                          class = as.factor(roads$katZarzadzania))
length_roads = aggregate(length ~ class, data = length_roads, FUN = sum)
```

The results are given in meters - let's convert them into kilometers.


```r
# convert to [km]
length_roads$length = units::set_units(length_roads$length, "km")
```

Let's also change the names of the categories.


```r
road_class = c("communal", "other", "national", "county", "voivodeship",
               "institutional")
levels(length_roads$class) = road_class
```

Now we can see the results.
The `other` type of roads dominates that consists mainly of non-public roads.


```r
length_roads
```

```
##           class          length
## 1      communal  186.17329 [km]
## 2         other 2819.55806 [km]
## 3      national   19.24996 [km]
## 4        county  187.17870 [km]
## 5   voivodeship  105.67567 [km]
## 6 institutional   10.44893 [km]
```

We can also calculate the total length of the roads in this area.


```r
sum(length_roads$length)
```

```
## 3328.285 [km]
```

The result is about 3328.285 km.
Another aspect of the data that we can investigate is the density of the road network.
We need to calculate the total area first, and then divide the total road length by the total area.


```r
communes_area = sum(st_area(communes))
communes_area = units::set_units(communes_area, "km2")
density = sum(length_roads$length)/communes_area
density = units::set_units(density, "km/km2")
density
```

```
## 2.921031 [km/km2]
```

The road density is about 2.92 km/km^2^.

### Roads through the rivers

Another dataset included in the Topographic Database contains rivers for this area ("PL.PZGiK.332.1801__OT_SWRS_L.xml").


```r
rivers = read_sf("data/PL.PZGiK.332.1801/BDOT10k/PL.PZGiK.332.1801__OT_SWRS_L.xml")
rivers = rivers[rivers$rodzaj == "Rz", ] # select only rivers
```

Rivers are divided into smaller sections with different parameters, such as river width or data source.
Let's merge sections from the same rivers into a single feature (geometry) - we can use the attribute with an ID (`idMPHP`) for that purpose.
We can also give each river a category number by creating a sequence from 1 to n using the `seq_len` function.


```r
rivers = aggregate(rivers[, c("geometry", "idMPHP")],
                   list(rivers$idMPHP),
                   sum)
rivers$idMPHP = seq_len(length(unique(rivers$idMPHP)))
rivers$idMPHP = as.factor(rivers$idMPHP)
```

Let's visualize the rivers' courses.


```r
plot(rivers["idMPHP"], main = "Rivers", reset = FALSE)
plot(st_geometry(territory), add = TRUE)
```

![](topodb-rivers-1.png)

With rivers and roads, we can designate points of intersection that symbolize bridges and crossings.
We can use the `st_intersection()` function for this.


```r
bridges = st_geometry(st_intersection(rivers, roads))
length(bridges)
```

```
## [1] 81
```

We get 81 such points.
Let's plot them.


```r
# use 'dev.off()' to reset previous plot
plot(st_geometry(rivers), main = "Bridges and crossings", col = "blue")
plot(st_geometry(territory), add = TRUE)
plot(bridges, add = TRUE, pch = 20)
```

![](topodb-bridges-1.png)



### Land cover

Land cover is the physical material at the surface of the earth like grass, trees, bare ground, water, etc.
Let's check the land cover data for our county - it is stored in files with the  `PT` prefix.
We use the `list.files()` function to list them.
The `pattern` argument is important here because it determines what files should be selected.
Our `pattern` should look like this: `PT+.+A\\.xml$` - only files containing areal data (A) of land cover (PT) will be listed.


```r
files = list.files("data/PL.PZGiK.332.1801/BDOT10k",
                   pattern = "PT+.+A\\.xml$",
                   full.names = TRUE)

# print filenames
basename(files)
```

```
##  [1] "PL.PZGiK.332.1801__OT_PTGN_A.xml"
##  [2] "PL.PZGiK.332.1801__OT_PTKM_A.xml"
##  [3] "PL.PZGiK.332.1801__OT_PTLZ_A.xml"
##  [4] "PL.PZGiK.332.1801__OT_PTNZ_A.xml"
##  [5] "PL.PZGiK.332.1801__OT_PTPL_A.xml"
##  [6] "PL.PZGiK.332.1801__OT_PTRK_A.xml"
##  [7] "PL.PZGiK.332.1801__OT_PTSO_A.xml"
##  [8] "PL.PZGiK.332.1801__OT_PTTR_A.xml"
##  [9] "PL.PZGiK.332.1801__OT_PTUT_A.xml"
## [10] "PL.PZGiK.332.1801__OT_PTWP_A.xml"
## [11] "PL.PZGiK.332.1801__OT_PTWZ_A.xml"
## [12] "PL.PZGiK.332.1801__OT_PTZB_A.xml"
```

We also can exclude "PL.PZGiK.332.1801__OT_PTSO_A.xml" file from the list, because this file contains zero objects (features).


```r
# drop "OT_PTSO_A.xml"
files = files[-7]
```

Let's prepare the names of the objects to which the data will be loaded.
The following names are the extended names of the abbreviations stored in the filenames.


```r
layer_names = c("fallowlands", "communication", "forest", "undeveloped",
                "squares", "shrublands", "crops", "grassland",
                "water", "heaps", "buildings")
```

Now we load each GML file, naming it from the list above.
Instead of using a loop, we can use the `lapply()` function, which performs a specific action for each element of the vector.
The action in our case is to load the GML files using `read_sf()`.


```r
layers = lapply(files, read_sf)
names(layers) = layer_names
```

Previously, we used `st_length()` to calculate the line length, now we use corresponding `st_area()` function to calculate the area.
Here we use the `lapply()` function analogously, which will work for each item on the list.
A similar function is `sapply()`, which returns a vector instead of a list.


```r
# calculate areas in each layer
area_landcover = lapply(layers, st_area)
# sum areas for each layer
area_landcover = sapply(area_landcover, sum)
# convert units
area_landcover = units::set_units(area_landcover, "m^2")
area_landcover = units::set_units(area_landcover, "km^2")
names(area_landcover) = layer_names
```

Let's see the results (in kilometers).


```r
area_landcover
```

```
## Units: [km^2]
##   fallowlands communication        forest   undeveloped       squares 
##    0.07561383    1.97014481  860.17248134    0.27817715    0.29038385 
##    shrublands         crops     grassland         water         heaps 
##    1.57266555  255.13528884    0.63447491    8.81181550    0.09220115 
##     buildings 
##   10.38781092
```

Let's make sure that the total land cover is equal to the area of our county.
Some tiny precision differences are possible, so we should set the difference tolerance.
This is possible using the `all.equal()` function.


```r
all.equal(sum(area_landcover), communes_area, tolerance = 0.001)
```

```
## [1] TRUE
```

Everything is correct.
Let's present the results as percentages of the area and sort them in descending order.


```r
landcover_percentage = area_landcover / sum(area_landcover) * 100
units(landcover_percentage) = NULL # drop units
landcover_percentage = sort(landcover_percentage, decreasing = TRUE)
landcover_percentage = round(landcover_percentage, 2)
landcover_percentage
```

```
##        forest         crops     buildings         water communication 
##         75.49         22.39          0.91          0.77          0.17 
##    shrublands     grassland       squares   undeveloped         heaps 
##          0.14          0.06          0.03          0.02          0.01 
##   fallowlands 
##          0.01
```

Over 75% of the county's area is covered by forests and only less than 1% by buildings.



### Buffer

In the last analysis in this vignette, we want to check how many buildings have bus stops within a given distance.
We can apply spatial buffers to solve this question.
Information about bus stops is in the "PL.PZGiK.332.1801__OT_OIKM_P.xml" file, where they are represented by the *OIKM04* value of the `x_kod` attribute.


```r
bus_stop = read_sf("data/PL.PZGiK.332.1801/BDOT10k/PL.PZGiK.332.1801__OT_OIKM_P.xml")
bus_stop = bus_stop[bus_stop$x_kod == "OIKM04", ]
```

Let's prepare a visualization in which the bus stops are marked with blue dots and the buildings are presented as red polygons.


```r
buildings = layers$buildings
plot(st_geometry(communes), main = "Bus stops")
plot(st_geometry(layers$buildings), add = TRUE, border = "red")
plot(st_geometry(bus_stop), add = TRUE, pch = 20, cex = 0.7, col = "blue")
```

![](topodb-busstops-1.png)

Let's create a buffer for each bus stop with a range of 1 km using `st_buffer()`.


```r
bus_buffer = st_buffer(bus_stop, 1000)
```

Now, we can plot it all.


```r
plot(st_geometry(communes), main = "Bus stops buffers")
plot(st_geometry(buildings), add = TRUE, border = "red")
plot(st_geometry(bus_buffer), add = TRUE)
```

![](topodb-buffers-1.png)

To return the buildings within the buffer range, we can perform the `st_within()` operation.


```r
buildings_buffer = st_within(buildings, bus_buffer)
```

The result is a nested list that consists of 2828 buildings and their associated buffers.
Let's count how many buildings are not in any buffer by using `sapply()` as in the previous examples.


```r
buildings_ex = sapply(buildings_buffer, length)
buildings_ex = sum(buildings_ex == 0)
buildings_ex = round(buildings_ex / nrow(buildings) * 100)
buildings_ex
```

```
## [1] 14
```

Answer to our last question: 14% of the buildings in this county do not have access to a bus stop within a 1 km radius.


---
title: Spatial Databases
author: "Krzysztof Dyba, Jakub Nowosad"
output: github_document
---

This article refers to two major spatial databases maintained by Polish Head Office of Geodesy and Cartography:

- **Topographic Database** (*Baza Danych Obiektów Topograficznych*),
- **General Geographic Database** (*Baza Danych Obiektów Ogólnogeograficznych*).

### Main categories of data

```{r include=FALSE}
main_abbre = c(
  "OT",
  "AD",
  "BU",
  "PT",
  "TC",
  "SW",
  "SK",
  "SU",
  "KU",
  "OI"
  )

main_pl = c(
  "Obiekt topograficzny",
  "Administracja",
  "Budowle",
  "Pokrycie terenu",
  "Tereny chronione",
  "Sieć wodna",
  "Sieć komunikacyjna",
  "Sieć uzbrojenia terenu",
  "Kompleksy użytkowania terenu",
  "Obiekty inne"
  )

main_en = c(
  "Topographic object",
  "Administration",
  "Buildings",
  "Land cover",
  "Protected areas",
  "Water network",
  "Communication network",
  "Utility infrastructure",
  "Land use",
  "Other objects")

df_main = data.frame(ABBRE = main_abbre, POLISH = main_pl, ENGLISH = main_en)
```

```{r echo=FALSE}
knitr::kable(df_main)
```

### Second level of categories (names of vector layers)

```{r include=FALSE}
abbre = c(
  "ADJA",
  "ADMS",
  "BUHD",
  "BUIN",
  "BUIT",
  "BUUO",
  "BUWT",
  "BUZM",
  "KUIK",
  "KUKO",
  "KUPG",
  "KUSC",
  "KUSK",
  "OIKM",
  "OIMK",
  "PTGN",
  "PTLZ",
  "PTNZ",
  "PTPL",
  "PTRK",
  "PTSO",
  "PTTR",
  "PTUT",
  "PTWP",
  "PTWZ",
  "PTZB",
  "SKDR",
  "SKPP",
  "SKRW",
  "SKTR",
  "SULN",
  "SWKN",
  "SWRM",
  "SWRS",
  "TCON",
  "TCPK",
  "TCPN",
  "TCRZ"
  )

pl = c(
  "Jednostka podziału administracyjnego",
  "Punkt główny miejscowości",
  "Budowla hydrotechniczna",
  "Budowla inżynierska",
  "Inne urządzenia techniczne",
  "Umocnienie drogowe, kolejowe i wodne",
  "Wysoka budowla techniczna",
  "Budowla ziemna",
  "Inny kompleks użytkowania terenu",
  "Kompleks komunikacyjny",
  "Kompleks przemysłowo-gospodarczy",
  "Kompleks sakralny i cmentarz",
  "Kompleks sportowo rekreacyjny",
  "Obiekt związany z komunikacją",
  "Mokradła",
  "Teren gruntów nieużytkowanych",
  "Teren leśny lub zadrzewiony",
  "Inny teren niezabudowany",
  "Teren placów",
  "Teren roślinności krzewiastej",
  "Teren składowania odpadów",
  "Teren roślinności trawiastej lub upraw rolnych",
  "Teren upraw trwałych",
  "Woda powierzchniowa",
  "Tereny zwałowisk i wyrobisk",
  "Zabudowa",
  "Droga",
  "Przeprawa",
  "Rondo lub węzeł drogowy",
  "Tor lub zespół torów",
  "Linia napowietrzna",
  "Kanał",
  "Rów melioracyjny",
  "Rzeka strumień",
  "Obszar Natura 2000",
  "Parki Krajobrazowe",
  "Parki Narodowe",
  "Rezerwaty"
  )

en = c(
  "Administrative subdivision",
  "Focal point of the town",
  "Hydraulic structure",
  "Engineering construction",
  "Other technical devices",
  "Road, rail and waterway reinforcements",
  "High technical construction",
  "Earth structure",
  "Other land use complex",
  "Communication complex",
  "Industrial and economic complex",
  "Religious complex and cemetery",
  "Sports and recreation complex",
  "Object related to communication",
  "Wetlands",
  "Unused lands area",
  "Forest or wooded area",
  "Other undeveloped area",
  "Squares area",
  "Shrubland area",
  "Landfill area",
  "Grassland or agricultural area",
  "Permanent crops area",
  "Surface water",
  "Heap and excavation areas",
  "Buildings",
  "Road",
  "Passage",
  "Roundabout or road junction",
  "Track or track set",
  "Overhead line",
  "Channel",
  "Drainage ditch",
  "River stream",
  "Natura 2000 area",
  "Landscape parks",
  "National park",
  "Nature reserve"
  )

df = data.frame(ABBRE = abbre, POLISH = pl, ENGLISH = en)
```

```{r echo=FALSE}
knitr::kable(df)
```

The last character in a vector layer name means:
  
  + **P** - point
  + **L** - line
  + **A** - area

### Other files

1. "**UzytkownikXX** (UserXX)" - contact details for the Voivodeship Marshal's Office. "XX" contains TERC (voivodeship ID).
2. Files with names: "**Ciek** (Watercourse)", "**LiniaKolejowa** (Railway line)", "**SzlakDrogowy** (Road trail)", "**WezelKolejowy** (Railway junction)", "**ZbiornikWodny** (Reservoir)" contain metadata for these objects.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tile_download.R
\name{tile_download}
\alias{tile_download}
\title{Download requested tiles}
\usage{
tile_download(
  df_req,
  outdir = ".",
  unzip = TRUE,
  check_SHA = FALSE,
  print_iter = TRUE,
  ...
)
}
\arguments{
\item{df_req}{a data frame obtained using the \code{\link[=ortho_request]{ortho_request()}} and
\code{\link[=DEM_request]{DEM_request()}} functions}

\item{outdir}{(optional) name of the output directory;
by default, files are saved in the working directory}

\item{unzip}{TRUE (default) or FALSE, when TRUE the downloaded archive will
be extracted and removed; only suitable for certain elevation data}

\item{check_SHA}{check the integrity of downloaded files
(logical, FALSE default)}

\item{print_iter}{print the current iteration of all
(logical, TRUE default)}

\item{...}{additional argument for \code{\link[utils:download.file]{utils::download.file()}}}
}
\value{
georeferenced tiles with properties (resolution, year, etc.)
as specified in the input data frame
}
\description{
Download requested tiles
}
\examples{
\dontrun{
library(sf)
polygon_path = system.file("datasets/search_area.gpkg", package = "rgugik")
polygon = read_sf(polygon_path)

req_df = ortho_request(polygon)
tile_download(req_df[1, ]) # download the first image only

req_df = DEM_request(polygon)
tile_download(req_df[1, ]) # download the first DEM only
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/emuia_download.R
\name{emuia_download}
\alias{emuia_download}
\title{Download Register of Towns, Streets and Addresses for communes}
\usage{
emuia_download(commune = NULL, TERYT = NULL, outdir = ".", unzip = TRUE, ...)
}
\arguments{
\item{commune}{commune name in Polish. Check \code{\link[=commune_names]{commune_names()}} function.}

\item{TERYT}{county ID (7 characters)}

\item{outdir}{(optional) name of the output directory;
by default, files are saved in the working directory}

\item{unzip}{TRUE (default) or FALSE, when TRUE the downloaded archive will
be extracted and removed}

\item{...}{additional argument for \code{\link[utils:download.file]{utils::download.file()}}}
}
\value{
a register in SHP format
}
\description{
Download Register of Towns, Streets and Addresses for communes
}
\examples{
\dontrun{
emuia_download(commune = "Kotla") # 38 KB
emuia_download(TERYT = c("0203042", "2412032")) # 75 KB
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geonames_download.R
\name{geonames_download}
\alias{geonames_download}
\title{Download State Register of Geographical Names}
\usage{
geonames_download(type, format = "SHP", outdir = ".", unzip = TRUE, ...)
}
\arguments{
\item{type}{names of places ("place") and/or physiographic objects ("object")}

\item{format}{data format ("GML", "SHP" (default) and/or "XLSX")}

\item{outdir}{(optional) name of the output directory;
by default, files are saved in the working directory}

\item{unzip}{TRUE (default) or FALSE, when TRUE the downloaded archive will
be extracted and removed}

\item{...}{additional argument for \code{\link[utils:download.file]{utils::download.file()}}}
}
\value{
a selected data type in the specified format
}
\description{
Download State Register of Geographical Names
}
\examples{
\dontrun{
geonames_download(type = "place", format = "SHP") # 18.2 MB
}
}
\references{
\url{http://isap.sejm.gov.pl/isap.nsf/download.xsp/WDU20150000219/O/D20150219.pdf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/borders_get.R
\name{borders_get}
\alias{borders_get}
\title{Get the boundaries of administrative units}
\usage{
borders_get(voivodeship = NULL, county = NULL, commune = NULL, TERYT = NULL)
}
\arguments{
\item{voivodeship}{selected voivodeships in Polish.
Check \code{\link[=voivodeship_names]{voivodeship_names()}} function}

\item{county}{county names in Polish.
Check \code{\link[=county_names]{county_names()}} function}

\item{commune}{commune names in Polish.
Check \code{\link[=commune_names]{commune_names()}} function}

\item{TERYT}{voivodeships, counties or communes (2, 4 or 7 characters)}
}
\value{
a sf data.frame (EPSG: 2180)
}
\description{
Get the boundaries of administrative units
}
\examples{
\dontrun{
voivodeship_geom = borders_get(voivodeship = "lubuskie") # 494 KB
county_geom = borders_get(county = "Sopot") # 18 KB
commune_geom = borders_get(commune = c("Hel", "Krynica Morska")) # 11 KB
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/borders_download.R
\name{borders_download}
\alias{borders_download}
\title{Download State Register of Borders}
\usage{
borders_download(type, outdir = ".", unzip = TRUE, ...)
}
\arguments{
\item{type}{"administrative units", "special units" or "addresses"}

\item{outdir}{(optional) name of the output directory;
by default, files are saved in the working directory}

\item{unzip}{TRUE (default) or FALSE, when TRUE the downloaded archive will
be extracted and removed}

\item{...}{additional argument for \code{\link[utils:download.file]{utils::download.file()}}}
}
\value{
a selected data type in SHP format
}
\description{
Download State Register of Borders
}
\examples{
\dontrun{
borders_download("administrative units") # 375 MB
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geocodePL_get.R
\name{geocodePL_get}
\alias{geocodePL_get}
\title{Convert addresses and objects to geographic coordinates}
\usage{
geocodePL_get(
  address = NULL,
  road = NULL,
  rail_crossing = NULL,
  geoname = NULL
)
}
\arguments{
\item{address}{place with or without street and house number}

\item{road}{road number with or without mileage}

\item{rail_crossing}{rail crossing identifier
(11 characters including 2 spaces, format: "XXX XXX XXX")}

\item{geoname}{name of the geographical object from State Register
of Geographical Names (function \code{\link[=geonames_download]{geonames_download()}})}
}
\value{
a sf data.frame (EPSG: 2180) with metadata
}
\description{
Convert addresses and objects to geographic coordinates
}
\examples{
\dontrun{
geocodePL_get(address = "Marki") # place
geocodePL_get(address = "Marki, Andersa") # place and street
geocodePL_get(address = "Marki, Andersa 1") # place, street and house number
geocodePL_get(address = "Królewskie Brzeziny 13") # place and house number

geocodePL_get(road = "632") # road number
geocodePL_get(road = "632 55") # road number and mileage

geocodePL_get(rail_crossing = "001 018 478")

geocodePL_get(geoname = "Las Mierzei") # physiographic object
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/models3D_download.R
\name{models3D_download}
\alias{models3D_download}
\title{Download 3D models of buildings for counties}
\usage{
models3D_download(
  county = NULL,
  TERYT = NULL,
  LOD = "LOD1",
  outdir = ".",
  unzip = TRUE,
  ...
)
}
\arguments{
\item{county}{county name in Polish. Check \code{\link[=county_names]{county_names()}} function.}

\item{TERYT}{county ID (4 characters)}

\item{LOD}{level of detail for building models ("LOD1" or "LOD2").
"LOD1" is default. "LOD2" is only available for ten voivodeships
(TERC: "04", "06", "12", "14", "16", "18", "20", "24", "26", "28").
Check \code{\link[=voivodeship_names]{voivodeship_names()}} function.}

\item{outdir}{(optional) name of the output directory;
by default, files are saved in the working directory}

\item{unzip}{TRUE (default) or FALSE, when TRUE the downloaded archive will
be extracted and removed}

\item{...}{additional argument for \code{\link[utils:download.file]{utils::download.file()}}}
}
\value{
models of buildings in Geography Markup Language format (.GML)
}
\description{
Download 3D models of buildings for counties
}
\examples{
\dontrun{
models3D_download(TERYT = c("2476", "2264")) # 3.6 MB
models3D_download(county = "sejneński", LOD = "LOD2") # 7.0 MB
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parcel_get.R
\name{parcel_get}
\alias{parcel_get}
\title{Get the geometry of cadastral parcels}
\usage{
parcel_get(TERYT = NULL, X = NULL, Y = NULL)
}
\arguments{
\item{TERYT}{parcel ID (18 characters, e.g. "141201_1.0001.6509")}

\item{X}{longitude (EPSG: 2180)}

\item{Y}{latitude (EPSG: 2180)}
}
\value{
a simple feature geometry (in case of TERYT) or data frame with simple
feature geometry and TERYT (in case of coordinates)
}
\description{
Get the geometry of cadastral parcels
}
\examples{
\dontrun{
parcel = parcel_get(TERYT = "141201_1.0001.6509")
parcel = parcel_get(X = 313380.5, Y = 460166.4)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rgugik-package.R
\docType{package}
\name{rgugik-package}
\alias{rgugik}
\alias{rgugik-package}
\title{rgugik: Search and Retrieve Spatial Data from 'GUGiK'}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Automatic open data acquisition from resources of Polish Head Office of Geodesy and Cartography ('Główny Urząd Geodezji i Kartografii') (<www.gugik.gov.pl>). Available datasets include various types of numeric, raster and vector data, such as orthophotomaps, digital elevation models (digital terrain models, digital surface model, point clouds), state register of borders, spatial databases, geometries of cadastral parcels, 3D models of buildings, and more. It is also possible to geocode addresses or objects using the geocodePL_get() function.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://kadyb.github.io/rgugik/}
  \item \url{https://github.com/kadyb/rgugik}
  \item Report bugs at \url{https://github.com/kadyb/rgugik/issues}
}

}
\author{
\strong{Maintainer}: Krzysztof Dyba \email{adres7@gmail.com} (\href{https://orcid.org/0000-0002-8614-3816}{ORCID})

Authors:
\itemize{
  \item Jakub Nowosad \email{nowosad.jakub@gmail.com} (\href{https://orcid.org/0000-0002-1057-3721}{ORCID})
}

Other contributors:
\itemize{
  \item Maciej Beręsewicz \email{maciej.beresewicz@ue.poznan.pl} (\href{https://orcid.org/0000-0002-8281-4301}{ORCID}) [contributor]
  \item GUGiK (source of the data) [contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/voivodeship_names.R
\docType{data}
\name{voivodeship_names}
\alias{voivodeship_names}
\title{Voivodeships in Poland}
\format{
An object of class \code{data.frame} with 16 rows and 3 columns.
}
\usage{
voivodeship_names
}
\description{
The data frame contains Polish and English names of
voivodeships, and their identifiers (TERC, 2 characters).
}
\examples{
voivodeship_names
}
\keyword{dataset}
\keyword{voivodeship}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/minmaxDTM_get.R
\name{minmaxDTM_get}
\alias{minmaxDTM_get}
\title{Get minimum and maximum elevation for a given polygon}
\usage{
minmaxDTM_get(polygon)
}
\arguments{
\item{polygon}{the polygon layer with only one object (area less than 10 ha),
the larger the polygon area, the lower DTM resolution,
the input coordinate system must be EPSG:2180}
}
\value{
a data frame with vector points and min/max terrain elevation
(EPSG:2180)
}
\description{
Get minimum and maximum elevation for a given polygon
}
\examples{
\dontrun{
library(sf)
polygon_path = system.file("datasets/search_area.gpkg", package = "rgugik")
polygon = read_sf(polygon_path)
minmax = minmaxDTM_get(polygon)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ortho_request.R
\name{ortho_request}
\alias{ortho_request}
\alias{orto_request}
\title{Get metadata and links to available orthoimages}
\usage{
ortho_request(x)

orto_request(x)
}
\arguments{
\item{x}{an \code{sf}, \code{sfc} or \code{SpatVector} object with one or more features
(requests are based on the bounding boxes of the provided features)}
}
\value{
a data frame with metadata and links to the orthoimages
}
\description{
Get metadata and links to available orthoimages
}
\examples{
\dontrun{
library(sf)
polygon_path = system.file("datasets/search_area.gpkg", package = "rgugik")
polygon = read_sf(polygon_path)
req_df = ortho_request(polygon)

# simple filtering by attributes
req_df = req_df[req_df$composition == "CIR", ]
req_df = req_df[req_df$resolution <= 0.25 & req_df$year >= 2016, ]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pointDTM_get.R
\name{pointDTM_get}
\alias{pointDTM_get}
\title{Get terrain elevation for a given polygon}
\usage{
pointDTM_get(polygon, distance = 1, print_iter = TRUE)
}
\arguments{
\item{polygon}{the polygon layer with only one object
(its area is limited to the 20 ha * distance parameter),
the input coordinate system must be EPSG:2180}

\item{distance}{distance between points in meters
(must be integer and greater than 1)}

\item{print_iter}{print the current iteration of all
(logical, TRUE default)}
}
\value{
a data frame with vector points and terrain elevation
(EPSG:2180, Vertical Reference System:PL-KRON86-NH)
}
\description{
Get terrain elevation for a given polygon
}
\examples{
\dontrun{
library(sf)
polygon_path = system.file("datasets/search_area.gpkg", package = "rgugik")
polygon = read_sf(polygon_path)
DTM = pointDTM_get(polygon, distance = 2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/topodb_download.R
\name{topodb_download}
\alias{topodb_download}
\title{Download Topographic Databases for counties}
\usage{
topodb_download(county = NULL, TERYT = NULL, outdir = ".", unzip = TRUE, ...)
}
\arguments{
\item{county}{county name in Polish. Check \code{\link[=county_names]{county_names()}} function.}

\item{TERYT}{county ID (4 characters)}

\item{outdir}{(optional) name of the output directory;
by default, files are saved in the working directory}

\item{unzip}{TRUE (default) or FALSE, when TRUE the downloaded archive will
be extracted and removed}

\item{...}{additional argument for \code{\link[utils:download.file]{utils::download.file()}}}
}
\value{
a database in Geography Markup Language format (.GML),
the content and detail level corresponds to the topographic map
in the scale of 1:10000
}
\description{
Download Topographic Databases for counties
}
\examples{
\dontrun{
topodb_download(county = "Świętochłowice") # 2.4 MB
topodb_download(TERYT = c("2476", "2264")) # 4.8 MB
}
}
\references{
description of topographical and general geographical databases,
and technical standards for making maps (in Polish):
\url{https://isap.sejm.gov.pl/isap.nsf/download.xsp/WDU20210001412/O/D20211412.pdf}

brief description of categories and layer names (in English and Polish):
\url{https://kadyb.github.io/rgugik/articles/articles/spatialdb_description.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commune_names.R
\docType{data}
\name{commune_names}
\alias{commune_names}
\title{Communes in Poland}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 2477 rows and 2 columns.
}
\usage{
commune_names
}
\description{
The data frame contains names of communes,
and their identifiers (TERC, 7 characters).
}
\examples{
commune_names
}
\keyword{commune}
\keyword{dataset}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pointDTM100_download.R
\name{pointDTM100_download}
\alias{pointDTM100_download}
\title{Download digital terrain models for voivodeships (100 m resolution)}
\usage{
pointDTM100_download(voivodeships, outdir = ".", unzip = TRUE, ...)
}
\arguments{
\item{voivodeships}{selected voivodeships in Polish or English, or TERC
(function \code{\link[=voivodeship_names]{voivodeship_names()}} can by helpful)}

\item{outdir}{(optional) name of the output directory;
by default, files are saved in the working directory}

\item{unzip}{TRUE (default) or FALSE, when TRUE the downloaded archive will
be extracted and removed}

\item{...}{additional argument for \code{\link[utils:download.file]{utils::download.file()}}}
}
\value{
text files with X, Y, Z columns (EPSG:2180)
}
\description{
Download digital terrain models for voivodeships (100 m resolution)
}
\examples{
\dontrun{
pointDTM100_download(c("opolskie", "świętokrzyskie")) # 8.5 MB
pointDTM100_download(c("Opole", "Swietokrzyskie")) # 8.5 MB
pointDTM100_download(c("16", "26")) # 8.5 MB
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/DEM_request.R
\name{DEM_request}
\alias{DEM_request}
\title{Get metadata and links to available digital elevation models}
\usage{
DEM_request(x)
}
\arguments{
\item{x}{an \code{sf}, \code{sfc} or \code{SpatVector} object with one or more features
(requests are based on the bounding boxes of the provided features)}
}
\value{
a data frame with metadata and links to the digital elevation models
(different formats of digital terrain model, digital surface model and
point clouds)
}
\description{
Get metadata and links to available digital elevation models
}
\examples{
\dontrun{
library(sf)
polygon_path = system.file("datasets/search_area.gpkg", package = "rgugik")
polygon = read_sf(polygon_path)
req_df = DEM_request(polygon)

# simple filtering by attributes
req_df = req_df[req_df$year > 2018, ]
req_df = req_df[req_df$product == "PointCloud" & req_df$format == "LAS", ]
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/county_names.R
\docType{data}
\name{county_names}
\alias{county_names}
\title{Counties in Poland}
\format{
An object of class \code{data.frame} with 380 rows and 3 columns.
}
\usage{
county_names
}
\description{
The data frame contains the names of counties,
their identifiers (TERYT, 4 characters) and the availability of
building models in the LOD2 standard (logical value).
}
\examples{
county_names
}
\keyword{county}
\keyword{dataset}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geodb_download.R
\name{geodb_download}
\alias{geodb_download}
\title{Download General Geographic Databases for entire voivodeships}
\usage{
geodb_download(voivodeships, outdir = ".", unzip = TRUE, ...)
}
\arguments{
\item{voivodeships}{selected voivodeships in Polish or English, or TERC
(object \code{\link{voivodeship_names}} can by helpful)}

\item{outdir}{(optional) name of the output directory;
by default, files are saved in the working directory}

\item{unzip}{TRUE (default) or FALSE, when TRUE the downloaded archive will
be extracted and removed}

\item{...}{additional argument for \code{\link[utils:download.file]{utils::download.file()}}}
}
\value{
a database in Geography Markup Language format (.GML),
the content and detail level corresponds to the general
geographic map in the scale of 1:250000
}
\description{
Download General Geographic Databases for entire voivodeships
}
\examples{
\dontrun{
geodb_download(c("opolskie", "lubuskie")) # 12.7 MB
geodb_download(c("Opole", "Lubusz")) # 12.7 MB
geodb_download(c("16", "08")) # 12.7 MB
}
}
\references{
description of topographical and general geographical databases,
and technical standards for making maps (in Polish):
\url{https://isap.sejm.gov.pl/isap.nsf/download.xsp/WDU20210001412/O/D20211412.pdf}

brief description of categories and layer names (in English and Polish):
\url{https://kadyb.github.io/rgugik/articles/articles/spatialdb_description.html}
}
