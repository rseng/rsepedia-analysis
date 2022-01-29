
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
