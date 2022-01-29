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

<!-- README.md is generated from README.Rmd. Please edit that file -->

# MODISTools <a href='https://github.com/ropensci/MODISTools'><img src='https://raw.githubusercontent.com/ropensci/MODISTools/master/MODISTools-logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/MODISTools/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/MODISTools/actions)
[![codecov](https://codecov.io/gh/ropensci/MODISTools/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/MODISTools)
![Status](https://www.r-pkg.org/badges/version/MODISTools) [![rOpenSci
Peer
Review](https://badges.ropensci.org/246_status.svg)](https://github.com/ropensci/software-review/issues/246)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/MODISTools)
<!-- badges: end -->

Programmatic interface to the [‘MODIS Land Products Subsets’ web
services](https://modis.ornl.gov/data/modis_webservice.html). Allows for
easy downloads of [‘MODIS’](http://modis.gsfc.nasa.gov/) time series
directly to your R workspace or your computer. When using the package
please cite the manuscript as referenced below. Keep in mind that the
original manuscript describes versions prior to release 1.0 of the
package. Functions described in this manuscript do not exist in the
current package, please consult [the
documentation](https://docs.ropensci.org/MODISTools/reference/index.html)
to find matching functionality.

## Installation

### stable release

To install the current stable release use a CRAN repository:

``` r
install.packages("MODISTools")
library("MODISTools")
```

### development release

To install the development releases of the package run the following
commands:

``` r
if(!require(devtools)){install.package("devtools")}
devtools::install_github("khufkens/MODISTools")
library("MODISTools")
```

Vignettes are not rendered by default, if you want to include additional
documentation please use:

``` r
if(!require(devtools)){install.package("devtools")}
devtools::install_github("khufkens/MODISTools", build_vignettes = TRUE)
library("MODISTools")
```

## Use

### Downloading MODIS time series

To extract a time series of modis data for a given location and its
direct environment use the mt_subset() function.

<details>
<summary>
detailed parameter description (click to expand)
</summary>
<p>

| Parameter | Description                                                                                                                     |
|-----------|---------------------------------------------------------------------------------------------------------------------------------|
| product   | a MODIS product                                                                                                                 |
| band      | a MODIS product band (if NULL all bands are downloaded)                                                                         |
| lat       | latitude of the site                                                                                                            |
| lon       | longitude of the site                                                                                                           |
| start     | start year of the time series (data start in 1980)                                                                              |
| end       | end year of the time series (current year - 2 years, use force = TRUE to override)                                              |
| internal  | logical, TRUE or FALSE, if true data is imported into R workspace otherwise it is downloaded into the current working directory |
| out_dir   | path where to store the data when not used internally, defaults to tempdir()                                                    |
| km_lr     | force “out of temporal range” downloads (integer)                                                                               |
| km_ab     | suppress the verbose output (integer)                                                                                           |
| site_name | a site identifier                                                                                                               |
| site_id   | a site_id for predefined locations (not required)                                                                               |
| progress  | logical, TRUE or FALSE (show download progress)                                                                                 |

</p>
</details>

``` r
# load the library
library(MODISTools)

# download data
subset <- mt_subset(product = "MOD11A2",
                    lat = 40,
                    lon = -110,
                    band = "LST_Day_1km",
                    start = "2004-01-01",
                    end = "2004-02-01",
                    km_lr = 1,
                    km_ab = 1,
                    site_name = "testsite",
                    internal = TRUE,
                    progress = FALSE)
print(str(subset))
#> 'data.frame':    36 obs. of  21 variables:
#>  $ xllcorner    : chr  "-9370962.97" "-9370962.97" "-9370962.97" "-9370962.97" ...
#>  $ yllcorner    : chr  "4446875.49" "4446875.49" "4446875.49" "4446875.49" ...
#>  $ cellsize     : chr  "926.625433055834" "926.625433055834" "926.625433055834" "926.625433055834" ...
#>  $ nrows        : int  3 3 3 3 3 3 3 3 3 3 ...
#>  $ ncols        : int  3 3 3 3 3 3 3 3 3 3 ...
#>  $ band         : chr  "LST_Day_1km" "LST_Day_1km" "LST_Day_1km" "LST_Day_1km" ...
#>  $ units        : chr  "Kelvin" "Kelvin" "Kelvin" "Kelvin" ...
#>  $ scale        : chr  "0.02" "0.02" "0.02" "0.02" ...
#>  $ latitude     : num  40 40 40 40 40 40 40 40 40 40 ...
#>  $ longitude    : num  -110 -110 -110 -110 -110 -110 -110 -110 -110 -110 ...
#>  $ site         : chr  "testsite" "testsite" "testsite" "testsite" ...
#>  $ product      : chr  "MOD11A2" "MOD11A2" "MOD11A2" "MOD11A2" ...
#>  $ start        : chr  "2004-01-01" "2004-01-01" "2004-01-01" "2004-01-01" ...
#>  $ end          : chr  "2004-02-01" "2004-02-01" "2004-02-01" "2004-02-01" ...
#>  $ complete     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ modis_date   : chr  "A2004001" "A2004009" "A2004017" "A2004025" ...
#>  $ calendar_date: chr  "2004-01-01" "2004-01-09" "2004-01-17" "2004-01-25" ...
#>  $ tile         : chr  "h09v05" "h09v05" "h09v05" "h09v05" ...
#>  $ proc_date    : chr  "2015212185706" "2015212201022" "2015212213103" "2015213005429" ...
#>  $ pixel        : int  1 1 1 1 2 2 2 2 3 3 ...
#>  $ value        : int  13135 13120 13350 13354 13123 13100 13324 13331 13098 13069 ...
#> NULL
```

The output format is a *tidy* data frame, as shown above. When witten to
a csv with the parameter `internal = FALSE` this will result in a flat
file on disk.

Note that when a a region is defined using km_lr and km_ab multiple
pixels might be returned. These are indexed using the `pixel` column in
the data frame containing the time series data. The remote sensing
values are listed in the `value` column. When no band is specified all
bands of a given product are returned, be mindful of the fact that
different bands might require different multipliers to represent their
true values. To list all available products, bands for particular
products and temporal coverage see function descriptions below.

### Batch downloading MODIS time series

When a large selection of locations is needed you might benefit from
using the batch download function `mt_batch_subset()`, which provides a
wrapper around the `mt_subset()` function in order to speed up large
download batches. This function has a similar syntax to `mt_subset()`
but requires a data frame defining site names (site_name) and locations
(lat / lon) (or a comma delimited file with the same structure) to
specify a list of download locations.

Below an example is provided on how to batch download data for a data
frame of given site names and locations (lat / lon).

``` r
# create data frame with a site_name, lat and lon column
# holding the respective names of sites and their location
df <- data.frame("site_name" = paste("test",1:2))
df$lat <- 40
df$lon <- -110
  
# test batch download
subsets <- mt_batch_subset(df = df,
                     product = "MOD11A2",
                     band = "LST_Day_1km",
                     internal = TRUE,
                     start = "2004-01-01",
                     end = "2004-02-01")

print(str(subsets))
#> 'data.frame':    8 obs. of  21 variables:
#>  $ xllcorner    : chr  "-9370036.39" "-9370036.39" "-9370036.39" "-9370036.39" ...
#>  $ yllcorner    : chr  "4447802.08" "4447802.08" "4447802.08" "4447802.08" ...
#>  $ cellsize     : chr  "926.625433055834" "926.625433055834" "926.625433055834" "926.625433055834" ...
#>  $ nrows        : int  1 1 1 1 1 1 1 1
#>  $ ncols        : int  1 1 1 1 1 1 1 1
#>  $ band         : chr  "LST_Day_1km" "LST_Day_1km" "LST_Day_1km" "LST_Day_1km" ...
#>  $ units        : chr  "Kelvin" "Kelvin" "Kelvin" "Kelvin" ...
#>  $ scale        : chr  "0.02" "0.02" "0.02" "0.02" ...
#>  $ latitude     : num  40 40 40 40 40 40 40 40
#>  $ longitude    : num  -110 -110 -110 -110 -110 -110 -110 -110
#>  $ site         : chr  "test 1" "test 1" "test 1" "test 1" ...
#>  $ product      : chr  "MOD11A2" "MOD11A2" "MOD11A2" "MOD11A2" ...
#>  $ start        : chr  "2004-01-01" "2004-01-01" "2004-01-01" "2004-01-01" ...
#>  $ end          : chr  "2004-02-01" "2004-02-01" "2004-02-01" "2004-02-01" ...
#>  $ complete     : logi  TRUE TRUE TRUE TRUE TRUE TRUE ...
#>  $ modis_date   : chr  "A2004001" "A2004009" "A2004017" "A2004025" ...
#>  $ calendar_date: chr  "2004-01-01" "2004-01-09" "2004-01-17" "2004-01-25" ...
#>  $ tile         : chr  "h09v05" "h09v05" "h09v05" "h09v05" ...
#>  $ proc_date    : chr  "2015212185706" "2015212201022" "2015212213103" "2015213005429" ...
#>  $ pixel        : int  1 1 1 1 1 1 1 1
#>  $ value        : int  13098 13062 13297 13323 13098 13062 13297 13323
#> NULL
```

### Listing products

To list all available products use the mt_products() function.

``` r
products <- mt_products()
head(products)
#>        product
#> 1       Daymet
#> 2 ECO4ESIPTJPL
#> 3      ECO4WUE
#> 4       GEDI03
#> 5      MCD12Q1
#> 6      MCD12Q2
#>                                                                       description
#> 1 Daily Surface Weather Data (Daymet) on a 1-km Grid for North America, Version 4
#> 2            ECOSTRESS Evaporative Stress Index PT-JPL (ESI) Daily L4 Global 70 m
#> 3                       ECOSTRESS Water Use Efficiency (WUE) Daily L4 Global 70 m
#> 4             GEDI Gridded Land Surface Metrics (LSM) L3 1km EASE-Grid, Version 2
#> 5           MODIS/Terra+Aqua Land Cover Type (LC) Yearly L3 Global 500 m SIN Grid
#> 6      MODIS/Terra+Aqua Land Cover Dynamics (LCD) Yearly L3 Global 500 m SIN Grid
#>   frequency resolution_meters
#> 1     1 day              1000
#> 2    Varies                70
#> 3    Varies                70
#> 4  One time              1000
#> 5    1 year               500
#> 6    1 year               500
```

### Listing bands

To list all available bands for a given product use the mt_bands()
function.

``` r
bands <- mt_bands(product = "MOD11A2")
head(bands)
#>               band                          description valid_range fill_value
#> 1   Clear_sky_days               Day clear-sky coverage    1 to 255          0
#> 2 Clear_sky_nights             Night clear-sky coverage    1 to 255          0
#> 3    Day_view_angl View zenith angle of day observation    0 to 130        255
#> 4    Day_view_time        Local time of day observation    0 to 240        255
#> 5          Emis_31                   Band 31 emissivity    1 to 255          0
#> 6          Emis_32                   Band 32 emissivity    1 to 255          0
#>    units scale_factor add_offset
#> 1   <NA>         <NA>       <NA>
#> 2   <NA>         <NA>       <NA>
#> 3 degree            1        -65
#> 4    hrs          0.1          0
#> 5   <NA>        0.002       0.49
#> 6   <NA>        0.002       0.49
```

### listing dates

To list all available dates (temporal coverage) for a given product and
location use the mt_dates() function.

``` r
dates <- mt_dates(product = "MOD11A2", lat = 42, lon = -110)
head(dates)
#>   modis_date calendar_date
#> 1   A2000049    2000-02-18
#> 2   A2000057    2000-02-26
#> 3   A2000065    2000-03-05
#> 4   A2000073    2000-03-13
#> 5   A2000081    2000-03-21
#> 6   A2000089    2000-03-29
```

## References

Tuck et al. (2014). [MODISTools - downloading and processing MODIS
remotely sensed data in R Ecology &
Evolution](https://onlinelibrary.wiley.com/doi/full/10.1002/ece3.1273),
4(24), 4658 - 4668.

## Acknowledgements

Original development was supported by the UK Natural Environment
Research Council (NERC; grants NE/K500811/1 and NE/J011193/1), and the
Hans Rausing Scholarship. Refactoring was supported through the Belgian
Science Policy office COBECORE project (BELSPO; grant
BR/175/A3/COBECORE). Logo design elements are taken from the FontAwesome
library according to [these terms](https://fontawesome.com/license),
where the globe element was inverted and intersected.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# MODISTools 1.1.2

* removing parallel processing in batch tool due to timeout issues

# MODISTools 1.1.1

* added raster conversion function `mt_to_raster()`
* force integer values on buffer values
* removed keywords

# MODISTools 1.1.0

* ropensci code review
* replaced NULL parameters with missing()
* output a tidy data frame
* include coordinate transform functions to sf objects
* pkgdown documentation generation

# MODISTools 1.0.0

* ropensci refactoring
* version increment follows original package
* push to CRAN

# MODISTools 0.1.0

* first draft
I have read and agree to the the CRAN policies at
http://cran.r-project.org/web/packages/policies.html

## test environments, local, CI and r-hub

- local OSX / Ubuntu 16.04 install on R 3.6.3
- Ubuntu 16.04 on Travis-CI (devel / release)
- rhub check on windows OK
- codecove.io code coverage at ~90%

## local / Travis CI R CMD check results

0 errors | 0 warnings | 0 notes
