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
smapr
================

[![codecov](https://codecov.io/gh/ropensci/smapr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/smapr)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/smapr)](https://cran.r-project.org/package=smapr)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](http://cranlogs.r-pkg.org/badges/grand-total/smapr)](http://cran.rstudio.com/web/packages/smapr/index.html)
[![](https://badges.ropensci.org/231_status.svg)](https://github.com/ropensci/onboarding/issues/231)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)

An R package for acquisition and processing of [NASA (Soil Moisture
Active-Passive) SMAP data](http://smap.jpl.nasa.gov/)

## Installation

To install smapr from CRAN:

``` r
install.packages("smapr")
```

To install the development version from GitHub:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/smapr")
```

#### Docker instructions (alternative to a local installation)

If a local installation is not possible for some reason, we have made a
Docker image available with smapr and all its dependencies.

    docker run -d -p 8787:8787 earthlab/smapr

In a web browser, navigate to localhost:8787 and log in with username:
rstudio, password: rstudio.

## Authentication

Access to the NASA SMAP data requires authentication through NASA’s
Earthdata portal. If you do not already have a username and password
through Earthdata, you can register for an account here:
<https://urs.earthdata.nasa.gov/> You cannot use this package without an
Earthdata account.

Once you have an account, you need to pass your Earthdata username
(`ed_un`) and password (`ed_pw`) as environmental variables that can be
read from within your R session. There are a couple of ways to do this:

### Recommended approach

Use `set_smap_credentials('yourusername', 'yourpasswd')`. This will save
your credentials by default, overwriting existing credentials if
`overwrite = TRUE`.

#### Alternative approaches

  - Use `Sys.setenv()` interactively in your R session to set your
    username and password (not including the `<` and `>`):

<!-- end list -->

``` r
Sys.setenv(ed_un = "<your username>", ed_pw = "<your password>")
```

  - Create a text file `.Renviron` in your home directory, which
    contains your username and password. If you don’t know what your
    home directory is, execute `normalizePath("~/")` in the R console
    and it will be printed. Be sure to include a new line at the end of
    the file or R will fail silently when loading it.

Example `.Renviron file` (note the new line at the end\!):

    ed_un=slkdjfsldkjfs
    ed_pw=dlfkjDD124^

Once this file is created, restart your R session and you should now be
able to access these environment variables (e.g., via
`Sys.getenv("ed_un")`).

# SMAP data products

Multiple SMAP data products are provided by the NSIDC, and these
products vary in the amount of processing. Currently, smapr primarily
supports level 3 and level 4 data products, which represent global daily
composite and global three hourly modeled data products, respectively.
There are a wide variety of data layers available in SMAP products,
including surface soil moisture, root zone soil moisture, freeze/thaw
status, surface temperature, vegetation water content, vegetation
opacity, net ecosystem carbon exchange, soil temperature, and
evapotranspiration. NSIDC provides documentation for all SMAP data
products on their [website](https://nsidc.org/data/smap/smap-data.html),
and we provide a summary of data products supported by smapr
below.

| Dataset id  | Description                                         | Resolution |
| ----------- | --------------------------------------------------- | ---------- |
| SPL2SMAP\_S | SMAP/Sentinel-1 Radiometer/Radar Soil Moisture      | 3 km       |
| SPL3FTA     | Radar Northern Hemisphere Daily Freeze/Thaw State   | 3 km       |
| SPL3SMA     | Radar Global Daily Soil Moisture                    | 3 km       |
| SPL3SMP     | Radiometer Global Soil Moisture                     | 36 km      |
| SPL3SMAP    | Radar/Radiometer Global Soil Moisture               | 9 km       |
| SPL4SMAU    | Surface/Rootzone Soil Moisture Analysis Update      | 9 km       |
| SPL4SMGP    | Surface/Rootzone Soil Moisture Geophysical Data     | 9 km       |
| SPL4SMLM    | Surface/Rootzone Soil Moisture Land Model Constants | 9 km       |
| SPL4CMDL    | Carbon Net Ecosystem Exchange                       | 9 km       |

## Typical workflow

At a high level, most workflows follow these steps:

1.  Find SMAP data with `find_smap()`
2.  Download data with `download_smap()`
3.  List data contents with `list_smap()`
4.  Extract data with `extract_smap()`

Each of these steps are outlined below:

### Finding SMAP data

Data are hosted on a server by the National Snow and Ice Data Center.
The `find_smap()` function searches for specific data products and
returns a data frame of available data. As data mature and pass checks,
versions advance. At any specific time, not all versions of all datasets
for all dates may exist. For the most up to date overview of dataset
versions, see the NSIDC SMAP data version
[webpage](https://nsidc.org/data/smap/smap-data.html).

``` r
library(smapr)
library(raster)
#> Loading required package: sp
available_data <- find_smap(id = "SPL3SMAP", date = "2015-05-25", version = 3)
str(available_data)
#> 'data.frame':    1 obs. of  3 variables:
#>  $ name: chr "SMAP_L3_SM_AP_20150525_R13080_001"
#>  $ date: Date, format: "2015-05-25"
#>  $ dir : chr "SPL3SMAP.003/2015.05.25/"
```

### Downloading and inspecting SMAP data

Given a data frame produced by `find_smap`, `download_smap` downloads
the data onto the local file system. Unless a directory is specified as
an argument, the data are stored in the user’s cache.

``` r
downloads <- download_smap(available_data)
#> Downloading https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMAP.003/2015.05.25/SMAP_L3_SM_AP_20150525_R13080_001.h5
#> Downloading https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMAP.003/2015.05.25/SMAP_L3_SM_AP_20150525_R13080_001.qa
#> Downloading https://n5eil01u.ecs.nsidc.org/SMAP/SPL3SMAP.003/2015.05.25/SMAP_L3_SM_AP_20150525_R13080_001.h5.iso.xml
str(downloads)
#> 'data.frame':    1 obs. of  4 variables:
#>  $ name     : chr "SMAP_L3_SM_AP_20150525_R13080_001"
#>  $ date     : Date, format: "2015-05-25"
#>  $ dir      : chr "SPL3SMAP.003/2015.05.25/"
#>  $ local_dir: chr "~/.cache/smap"
```

The SMAP data are provided in HDF5 format, and in any one file there are
actually multiple data sets, including metadata. The `list_smap`
function allows users to inspect the contents of downloaded data at a
high level (`all = FALSE`) or in depth (`all = TRUE`).

``` r
list_smap(downloads, all = FALSE)
#> $SMAP_L3_SM_AP_20150525_R13080_001
#>                           name group     otype dclass  dim
#> 1                     Metadata     . H5I_GROUP   <NA> <NA>
#> 2 Soil_Moisture_Retrieval_Data     . H5I_GROUP   <NA> <NA>
```

To see all of the data fields, set `all = TRUE`.

### Extracting gridded data products

The `extract_smap` function extracts gridded data products (e.g., global
soil moisture) and returns Raster\* objects. If more than one file has
been downloaded and passed into the first argument, `extract_smap`
extracts all of the rasters and returns a
RasterStack.

``` r
sm_raster <- extract_smap(downloads, "Soil_Moisture_Retrieval_Data/soil_moisture")
plot(sm_raster, main = "Level 3 soil moisture: May 25, 2015")
```

<img src="man/figures/extract-data-1.png" style="display: block; margin: auto;" />

The path “Soil\_Moisture\_Retrieval\_Data/soil\_moisture” was determined
from the output of `list_smap(downloads, all = TRUE)`, which lists all
of the data contained in SMAP data files.

### Saving GeoTIFF output

The raster stack can be saved as a GeoTIFF using the `writeRaster`
function from the raster pacakge.

``` r
writeRaster(sm_raster, "sm_raster.tif")
```

## Meta

  - Please [report any issues or
    bugs](https://github.com/ropensci/smapr/issues), after reading our
    contribution [guidelines](CONTRIBUTING.md), and the [Contributor
    Code of Conduct](CONDUCT.md).
  - License: GPL-3
  - See `citation("smapr")` in R to cite this package in
publications.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# smapr 0.2.1

* patch to skip test on cran that requires internet
* updated SMAP data versions

# smapr 0.2.0

* added set_smap_credentials() function for NASA EarthData portal
* expanded vignettes to include cropping, masking, etc.
* added verbose argument to download_smap
* performance improvements to extract_smap()

# smapr 0.1.2

* added support for SMAP/sentinel hybrid soil moisture product
* added a code of conduct
* added a vignette to show a complete workflow

# smapr 0.1.1

* added patch for searching date ranges containing missing collections
* added unit tests for user specified download directories
* updates to examples for new data versions
* adding a CONTRIBUTING.md file

# smapr 0.1.0

* updating remote data location (previous ftp server was removed)
* using NASA Earthdata authentication

# smapr 0.0.1

* first submission to CRAN
# CONTRIBUTING #

### Please contribute!
We love collaboration.

### Found a Bug?

* Submit an issue on our Issues page [here](https://github.com/earthlab/smapr/issues).

### Code contributions?

* **Fork** this repo to your Github account.
* **Clone** your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/smapr.git`.
* Make sure to **track upstream** progress (i.e., on our version of `smapr` at `earthlab/smapr`) by doing `git remote add upstream https://github.com/earthlab/smapr.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your **changes** (bonus points for making changes on a new branch).
* **Push** up to your account.
* Submit a **pull request** to home base at `earthlab/smapr`.

Please follow [this](http://adv-r.had.co.nz/Style.html) styleguide for your contributions.

### Questions?

Get in touch: [maxwell.b.joseph@colorado.edu](mailto:maxwell.b.joseph@colorado.edu)

### Thanks for contributing!
## Test environments
* Ubuntu Linux 14.04 (on travis-ci), R-release, R-devel, R-oldrel
* Windows Server Windows Server 2012 R2 x64 (on appveyor)
* Debian Linux, R-release, GCC (on r-hub)
* Ubuntu Linux 16.04 LTS, R-devel, GCC (on r-hub)
* Ubuntu Linux 16.04 LTS, R-release, GCC (on r-hub)

## R CMD check results
There were no ERRORs or WARNINGs. 

Possibly mis-spelled words in DESCRIPTION: 
  SMAP

SMAP is not a misspelling, it's an acronym for Soil Moisture Active Passive. 

## Downstream dependencies

There are currently no downstream dependencies for this package.
# Before posting

The GitHub issue tracker is intended for bug reports and feature requests. 
Do not post your NASA Earthdata username or password with your issue!

When you post, please include a minimal reproducible example of the problem 
and/or desired behavior, if applicable. 
---
title: "smapr"
output: github_document
---

[![codecov](https://codecov.io/gh/ropensci/smapr/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/smapr)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/smapr)](https://cran.r-project.org/package=smapr)
[![lifecycle](https://img.shields.io/badge/lifecycle-maturing-blue.svg)](https://www.tidyverse.org/lifecycle/#maturing)
[![](http://cranlogs.r-pkg.org/badges/grand-total/smapr)](http://cran.rstudio.com/web/packages/smapr/index.html) 
[![](https://badges.ropensci.org/231_status.svg)](https://github.com/ropensci/onboarding/issues/231)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](http://www.repostatus.org/badges/latest/active.svg)](http://www.repostatus.org/#active)


```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/"
)
```


An R package for acquisition and processing of [NASA (Soil Moisture Active-Passive) SMAP data](http://smap.jpl.nasa.gov/)

## Installation

To install smapr from CRAN: 

```{r cran-installation, eval = FALSE}
install.packages("smapr")
```

To install the development version from GitHub:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/smapr")
```

#### Docker instructions (alternative to a local installation)

If a local installation is not possible for some reason, we have made a Docker 
image available with smapr and all its dependencies.

```
docker run -d -p 8787:8787 earthlab/smapr
```

In a web browser, navigate to localhost:8787 and log in with 
username: rstudio, password: rstudio.


## Authentication

Access to the NASA SMAP data requires authentication through NASA's Earthdata 
portal. 
If you do not already have a username and password through Earthdata, you can 
register for an account here: https://urs.earthdata.nasa.gov/
You cannot use this package without an Earthdata account. 

Once you have an account, you need to pass your Earthdata username (`ed_un`) 
and password (`ed_pw`) as environmental variables that can be read from within 
your R session. 
There are a couple of ways to do this: 

### Recommended approach

Use `set_smap_credentials('yourusername', 'yourpasswd')`. 
This will save your credentials by default, overwriting existing credentials if 
`overwrite = TRUE`. 

#### Alternative approaches

- Use `Sys.setenv()` interactively in your R session to set your username and 
password (not including the `<` and `>`):

```{r, eval = FALSE}
Sys.setenv(ed_un = "<your username>", ed_pw = "<your password>")
```

- Create a text file `.Renviron` in your home directory, which contains your 
username and password. 
If you don't know what your home directory is, execute `normalizePath("~/")` in 
the R console and it will be printed.
Be sure to include a new line at the end of the file or R will fail silently 
when loading it.

Example `.Renviron file` (note the new line at the end!):

```
ed_un=slkdjfsldkjfs
ed_pw=dlfkjDD124^

```

Once this file is created, restart your R session and you should now be able to 
access these environment variables (e.g., via `Sys.getenv("ed_un")`).



# SMAP data products

Multiple SMAP data products are provided by the NSIDC, and these products vary 
in the amount of processing. 
Currently, smapr primarily supports level 3 and level 4 data products, 
which represent global daily composite and global three hourly modeled data 
products, respectively. 
There are a wide variety of data layers available in SMAP products, including surface soil moisture, root zone soil moisture, freeze/thaw status, surface temperature, vegetation water content, vegetation opacity, net ecosystem carbon exchange, soil temperature, and evapotranspiration. 
NSIDC provides documentation for all SMAP data products on their 
[website](https://nsidc.org/data/smap/smap-data.html), and we provide a summary 
of data products supported by smapr below. 

| Dataset id | Description                                         | Resolution |
|------------|-----------------------------------------------------|------------|
| SPL2SMAP_S | SMAP/Sentinel-1 Radiometer/Radar Soil Moisture      | 3 km       |
| SPL3FTA    | Radar Northern Hemisphere Daily Freeze/Thaw State   | 3 km       |
| SPL3SMA    | Radar Global Daily Soil Moisture                    | 3 km       |
| SPL3SMP    | Radiometer Global Soil Moisture                     | 36 km      |
| SPL3SMAP   | Radar/Radiometer Global Soil Moisture               | 9 km       |
| SPL4SMAU   | Surface/Rootzone Soil Moisture Analysis Update      | 9 km       | 
| SPL4SMGP   | Surface/Rootzone Soil Moisture Geophysical Data     | 9 km       |
| SPL4SMLM   | Surface/Rootzone Soil Moisture Land Model Constants | 9 km       |
| SPL4CMDL   | Carbon Net Ecosystem Exchange                       | 9 km       |

## Typical workflow

At a high level, most workflows follow these steps:

1. Find SMAP data with `find_smap()`
2. Download data with `download_smap()`
3. List data contents with `list_smap()`
4. Extract data with `extract_smap()`

Each of these steps are outlined below:

### Finding SMAP data

Data are hosted on a server by the National Snow and Ice Data Center. 
The `find_smap()` function searches for specific data products and returns a 
data frame of available data.
As data mature and pass checks, versions advance. 
At any specific time, not all versions of all datasets for all dates may exist. 
For the most up to date overview of dataset versions, see the NSIDC SMAP data 
version [webpage](https://nsidc.org/data/smap/smap-data.html).

```{r find-data}
library(smapr)
library(raster)
available_data <- find_smap(id = "SPL3SMAP", date = "2015-05-25", version = 3)
str(available_data)
```

### Downloading and inspecting SMAP data

Given a data frame produced by `find_smap`, `download_smap` downloads the data 
onto the local file system. 
Unless a directory is specified as an argument, the data are stored in the 
user's cache. 

```{r download-data}
downloads <- download_smap(available_data)
str(downloads)
```

The SMAP data are provided in HDF5 format, and in any one file there are 
actually multiple data sets, including metadata. 
The `list_smap` function allows users to inspect the contents of downloaded 
data at a high level (`all = FALSE`) or in depth (`all = TRUE`). 

```{r list-data}
list_smap(downloads, all = FALSE)
```

To see all of the data fields, set `all = TRUE`. 

### Extracting gridded data products

The `extract_smap` function extracts gridded data products 
(e.g., global soil moisture) and returns Raster* objects. 
If more than one file has been downloaded and passed into the first argument, `extract_smap` extracts all of the rasters and returns a RasterStack.

```{r extract-data, fig.align='center', fig.width=8, fig.height=7}
sm_raster <- extract_smap(downloads, "Soil_Moisture_Retrieval_Data/soil_moisture")
plot(sm_raster, main = "Level 3 soil moisture: May 25, 2015")
```

The path "Soil_Moisture_Retrieval_Data/soil_moisture" was determined from the 
output of `list_smap(downloads, all = TRUE)`, which lists all of the data 
contained in SMAP data files. 

### Saving GeoTIFF output

The raster stack can be saved as a GeoTIFF using the `writeRaster` function 
from the raster pacakge. 

```{r}
writeRaster(sm_raster, "sm_raster.tif")
```

```{r, echo = FALSE, results='hide'}
# cleanup
file.remove("sm_raster.tif")
```


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/smapr/issues),
after reading our contribution [guidelines](CONTRIBUTING.md), and the 
[Contributor Code of Conduct](CONDUCT.md). 
* License: GPL-3
* See `citation("smapr")` in R to cite this package in publications. 

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: "Introduction to the smapr package"
author: "Maxwell B. Joseph"
date: "2022-01-26"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Introduction to the smapr package}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library(smapr)
library(sp)
library(raster)
```


This vignette outlines a basic use scenario for smapr.
We will acquire and process
[NASA (Soil Moisture Active-Passive) SMAP data](http://smap.jpl.nasa.gov/),
and generate some simple visualizations.


## SMAP data products

Multiple SMAP data products are provided by the NSIDC, and these products vary
in the amount of processing.
Currently, smapr primarily supports level 3 and level 4 data products,
which represent global daily composite and global three hourly modeled data
products, respectively.
NSIDC provides documentation for all SMAP data products on their
[website](https://nsidc.org/data/smap/smap-data.html), and we provide a summary
of data products supported by smapr below.

| Dataset id | Description                                         | Resolution |
|------------|-----------------------------------------------------|------------|
| SPL2SMAP_S | SMAP/Sentinel-1 Radiometer/Radar Soil Moisture      | 3 km       |
| SPL3FTA    | Radar Northern Hemisphere Daily Freeze/Thaw State   | 3 km       |
| SPL3SMA    | Radar Global Daily Soil Moisture                    | 3 km       |
| SPL3SMP    | Radiometer Global Soil Moisture                     | 36 km      |
| SPL3SMAP   | Radar/Radiometer Global Soil Moisture               | 9 km       |
| SPL4SMAU   | Surface/Rootzone Soil Moisture Analysis Update      | 9 km       |
| SPL4SMGP   | Surface/Rootzone Soil Moisture Geophysical Data     | 9 km       |
| SPL4SMLM   | Surface/Rootzone Soil Moisture Land Model Constants | 9 km       |
| SPL4CMDL   | Carbon Net Ecosystem Exchange                       | 9 km       |


This vignette uses the level 4 [SPL4SMAU](https://nsidc.org/data/SPL4SMAU)
(Surface/Rootzone Soil Moisture Analysis Update) data product.

## Preparing to access SMAP data

NASA requires a username and password from their Earthdata portal to access
SMAP data.
You can get these credentials here: https://earthdata.nasa.gov/

Once you have your credentials, you can use the `set_smap_credentials`
function to set them for use by the smapr package:


```r
set_smap_credentials("myusername", "mypassword")
```

This function saves your credentials for later use unless you use the argument
`save = FALSE`.

## Finding data

To find out which SMAP data are available, we'll use the `find_smap` function,
which takes a data set ID, date(s) to search, and a dataset version.


```r
available_data <- find_smap(id = 'SPL4SMAU', dates = '2018-06-01', version = 5)
```

This returns a data frame, where every row is one data file that is available
on NASA's servers.


```r
str(available_data)
#> 'data.frame':	8 obs. of  3 variables:
#>  $ name: chr  "SMAP_L4_SM_aup_20180601T030000_Vv5030_001" "SMAP_L4_SM_aup_20180601T060000_Vv5030_001" "SMAP_L4_SM_aup_20180601T090000_Vv5030_001" "SMAP_L4_SM_aup_20180601T120000_Vv5030_001" ...
#>  $ date: Date, format: "2018-06-01" "2018-06-01" "2018-06-01" "2018-06-01" ...
#>  $ dir : chr  "SPL4SMAU.005/2018.06.01/" "SPL4SMAU.005/2018.06.01/" "SPL4SMAU.005/2018.06.01/" "SPL4SMAU.005/2018.06.01/" ...
```

## Downloading data

To download the data, we can use `download_smap`. Note that this may take a
while, depending on the number of files being downloaded, and the speed of your
internet connection.
Because we're downloading multiple files, we will use the
`verbose = FALSE` argument to avoid printing excessive output to the console.


```r
local_files <- download_smap(available_data, overwrite = FALSE, verbose = FALSE)
```

Each file corresponds to different
times as indicated by the file names:


```r
local_files$name[1:2]
#> [1] "SMAP_L4_SM_aup_20180601T030000_Vv5030_001" "SMAP_L4_SM_aup_20180601T060000_Vv5030_001"
```

## Exploring data

Each file that we downloaded is an HDF5 file with multiple datasets bundled
together.
To list all of the data in a file we can use `list_smap`.
By default, if we give `list_smap` a data frame of local files, it will
return a list of data frames.
Because all of these data files are of the same data product, using `list_smap`
on one file (e.g., the first) will tell us what's available in all of the files:


```r
list_smap(local_files[1, ])
#> $SMAP_L4_SM_aup_20180601T030000_Vv5030_001
#>                                name                              group       otype      dclass         dim
#> 1                                 y                                  . H5I_DATASET   H5T_FLOAT        1624
#> 2                     Forecast_Data                                  .   H5I_GROUP        <NA>        <NA>
#> 3               sm_surface_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 4                     tb_v_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 5             surface_temp_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 6              tb_v_forecast_ensstd                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 7         soil_temp_layer1_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 8               sm_profile_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 9                     tb_h_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 10             sm_rootzone_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 11             tb_h_forecast_ensstd                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 12                             time                                  . H5I_DATASET   H5T_FLOAT           1
#> 13          EASE2_global_projection                                  . H5I_DATASET  H5T_STRING           1
#> 14                    Analysis_Data                                  .   H5I_GROUP        <NA>        <NA>
#> 15       sm_surface_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 16     surface_temp_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 17      sm_rootzone_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 18              sm_surface_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 19       sm_profile_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 20             sm_rootzone_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 21              sm_profile_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 22        soil_temp_layer1_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 23 soil_temp_layer1_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 24            surface_temp_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 25                         cell_lat                                  . H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 26                         cell_row                                  . H5I_DATASET H5T_INTEGER 3856 x 1624
#> 27                         cell_lon                                  . H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 28                         Metadata                                  .   H5I_GROUP        <NA>        <NA>
#> 29                           Source                           Metadata   H5I_GROUP        <NA>        <NA>
#> 30                           L1C_TB                    Metadata/Source   H5I_GROUP        <NA>        <NA>
#> 31           AcquisitionInformation                           Metadata   H5I_GROUP        <NA>        <NA>
#> 32                         platform    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 33                 platformDocument    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 34                    radarDocument    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 35                            radar    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 36               radiometerDocument    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 37                       radiometer    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 38                      DataQuality                           Metadata   H5I_GROUP        <NA>        <NA>
#> 39                              TBH               Metadata/DataQuality   H5I_GROUP        <NA>        <NA>
#> 40                DomainConsistency           Metadata/DataQuality/TBH   H5I_GROUP        <NA>        <NA>
#> 41             CompletenessOmission           Metadata/DataQuality/TBH   H5I_GROUP        <NA>        <NA>
#> 42                              TBV               Metadata/DataQuality   H5I_GROUP        <NA>        <NA>
#> 43                DomainConsistency           Metadata/DataQuality/TBV   H5I_GROUP        <NA>        <NA>
#> 44             CompletenessOmission           Metadata/DataQuality/TBV   H5I_GROUP        <NA>        <NA>
#> 45             SeriesIdentification                           Metadata   H5I_GROUP        <NA>        <NA>
#> 46            DatasetIdentification                           Metadata   H5I_GROUP        <NA>        <NA>
#> 47                           Extent                           Metadata   H5I_GROUP        <NA>        <NA>
#> 48                             CRID                           Metadata   H5I_GROUP        <NA>        <NA>
#> 49                              AUP                      Metadata/CRID   H5I_GROUP        <NA>        <NA>
#> 50                             Root                      Metadata/CRID   H5I_GROUP        <NA>        <NA>
#> 51                           Config                           Metadata   H5I_GROUP        <NA>        <NA>
#> 52        GridSpatialRepresentation                           Metadata   H5I_GROUP        <NA>        <NA>
#> 53                         Latitude Metadata/GridSpatialRepresentation   H5I_GROUP        <NA>        <NA>
#> 54                        Longitude Metadata/GridSpatialRepresentation   H5I_GROUP        <NA>        <NA>
#> 55                      ProcessStep                           Metadata   H5I_GROUP        <NA>        <NA>
#> 56                      cell_column                                  . H5I_DATASET H5T_INTEGER 3856 x 1624
#> 57                Observations_Data                                  .   H5I_GROUP        <NA>        <NA>
#> 58                         tb_v_obs                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 59                tb_v_obs_time_sec                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 60                         tb_h_obs                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 61                  tb_h_orbit_flag                  Observations_Data H5I_DATASET H5T_INTEGER 3856 x 1624
#> 62             tb_h_resolution_flag                  Observations_Data H5I_DATASET H5T_INTEGER 3856 x 1624
#> 63                  tb_h_obs_errstd                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 64                  tb_v_orbit_flag                  Observations_Data H5I_DATASET H5T_INTEGER 3856 x 1624
#> 65             tb_v_resolution_flag                  Observations_Data H5I_DATASET H5T_INTEGER 3856 x 1624
#> 66                   tb_v_obs_assim                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 67                   tb_h_obs_assim                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 68                tb_h_obs_time_sec                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 69                  tb_v_obs_errstd                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 70                                x                                  . H5I_DATASET   H5T_FLOAT        3856
```

To dig deeper, we can use the `all` argument to `list_smap`:


```r
list_smap(local_files[1, ], all = TRUE)
#> $SMAP_L4_SM_aup_20180601T030000_Vv5030_001
#>                                name                              group       otype      dclass         dim
#> 1                                 y                                  . H5I_DATASET   H5T_FLOAT        1624
#> 2                     Forecast_Data                                  .   H5I_GROUP        <NA>        <NA>
#> 3               sm_surface_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 4                     tb_v_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 5             surface_temp_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 6              tb_v_forecast_ensstd                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 7         soil_temp_layer1_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 8               sm_profile_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 9                     tb_h_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 10             sm_rootzone_forecast                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 11             tb_h_forecast_ensstd                      Forecast_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 12                             time                                  . H5I_DATASET   H5T_FLOAT           1
#> 13          EASE2_global_projection                                  . H5I_DATASET  H5T_STRING           1
#> 14                    Analysis_Data                                  .   H5I_GROUP        <NA>        <NA>
#> 15       sm_surface_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 16     surface_temp_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 17      sm_rootzone_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 18              sm_surface_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 19       sm_profile_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 20             sm_rootzone_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 21              sm_profile_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 22        soil_temp_layer1_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 23 soil_temp_layer1_analysis_ensstd                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 24            surface_temp_analysis                      Analysis_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 25                         cell_lat                                  . H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 26                         cell_row                                  . H5I_DATASET H5T_INTEGER 3856 x 1624
#> 27                         cell_lon                                  . H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 28                         Metadata                                  .   H5I_GROUP        <NA>        <NA>
#> 29                           Source                           Metadata   H5I_GROUP        <NA>        <NA>
#> 30                           L1C_TB                    Metadata/Source   H5I_GROUP        <NA>        <NA>
#> 31           AcquisitionInformation                           Metadata   H5I_GROUP        <NA>        <NA>
#> 32                         platform    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 33                 platformDocument    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 34                    radarDocument    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 35                            radar    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 36               radiometerDocument    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 37                       radiometer    Metadata/AcquisitionInformation   H5I_GROUP        <NA>        <NA>
#> 38                      DataQuality                           Metadata   H5I_GROUP        <NA>        <NA>
#> 39                              TBH               Metadata/DataQuality   H5I_GROUP        <NA>        <NA>
#> 40                DomainConsistency           Metadata/DataQuality/TBH   H5I_GROUP        <NA>        <NA>
#> 41             CompletenessOmission           Metadata/DataQuality/TBH   H5I_GROUP        <NA>        <NA>
#> 42                              TBV               Metadata/DataQuality   H5I_GROUP        <NA>        <NA>
#> 43                DomainConsistency           Metadata/DataQuality/TBV   H5I_GROUP        <NA>        <NA>
#> 44             CompletenessOmission           Metadata/DataQuality/TBV   H5I_GROUP        <NA>        <NA>
#> 45             SeriesIdentification                           Metadata   H5I_GROUP        <NA>        <NA>
#> 46            DatasetIdentification                           Metadata   H5I_GROUP        <NA>        <NA>
#> 47                           Extent                           Metadata   H5I_GROUP        <NA>        <NA>
#> 48                             CRID                           Metadata   H5I_GROUP        <NA>        <NA>
#> 49                              AUP                      Metadata/CRID   H5I_GROUP        <NA>        <NA>
#> 50                             Root                      Metadata/CRID   H5I_GROUP        <NA>        <NA>
#> 51                           Config                           Metadata   H5I_GROUP        <NA>        <NA>
#> 52        GridSpatialRepresentation                           Metadata   H5I_GROUP        <NA>        <NA>
#> 53                         Latitude Metadata/GridSpatialRepresentation   H5I_GROUP        <NA>        <NA>
#> 54                        Longitude Metadata/GridSpatialRepresentation   H5I_GROUP        <NA>        <NA>
#> 55                      ProcessStep                           Metadata   H5I_GROUP        <NA>        <NA>
#> 56                      cell_column                                  . H5I_DATASET H5T_INTEGER 3856 x 1624
#> 57                Observations_Data                                  .   H5I_GROUP        <NA>        <NA>
#> 58                         tb_v_obs                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 59                tb_v_obs_time_sec                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 60                         tb_h_obs                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 61                  tb_h_orbit_flag                  Observations_Data H5I_DATASET H5T_INTEGER 3856 x 1624
#> 62             tb_h_resolution_flag                  Observations_Data H5I_DATASET H5T_INTEGER 3856 x 1624
#> 63                  tb_h_obs_errstd                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 64                  tb_v_orbit_flag                  Observations_Data H5I_DATASET H5T_INTEGER 3856 x 1624
#> 65             tb_v_resolution_flag                  Observations_Data H5I_DATASET H5T_INTEGER 3856 x 1624
#> 66                   tb_v_obs_assim                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 67                   tb_h_obs_assim                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 68                tb_h_obs_time_sec                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 69                  tb_v_obs_errstd                  Observations_Data H5I_DATASET   H5T_FLOAT 3856 x 1624
#> 70                                x                                  . H5I_DATASET   H5T_FLOAT        3856
```

Looking at this output, we can conclude that the file contains multiple arrays
(notice the `dim` column).
These arrays correspond to things like estimated root zone soil moisture
(`/Analysis_Data/sm_rootzone_analysis`), estimated surface soil moisture
(`/Analysis_Data/sm_surface_analysis`), and estimated surface temperature
(`/Analysis_Data/surface_temp_analysis`).
See https://nsidc.org/data/smap/spl4sm/data-fields#sm_surface_analysis for more
detailed information on what these datasets represent and how they were
generated.

## Extracting data

The datasets that we are interested in are spatial grids.
The `smapr` package can extract these data into `raster` objects with the
`extract_smap` function, which takes a dataset name as an argument.
These names are paths that can be generated from the output of `list_smap`.
For example, if we want to get rootzone soil moisture, we can see a dataset
with name `sm_rootzone_analysis` in group `/Analysis_Data`, so that the path
to the dataset is `/Analysis_Data/sm_rootzone_analysis`:


```r
sm_raster <- extract_smap(local_files, '/Analysis_Data/sm_rootzone_analysis')
```

This will extract all of the data in the data frame `local_files`, generating
a RasterBrick with one layer per file:


```r
sm_raster
#> class      : RasterBrick 
#> dimensions : 1624, 3856, 6262144, 8  (nrow, ncol, ncell, nlayers)
#> resolution : 9008.055, 9008.055  (x, y)
#> extent     : -17367530, 17367530, -7314541, 7314541  (xmin, xmax, ymin, ymax)
#> crs        : +proj=cea +lat_ts=30 +lon_0=0 +x_0=0 +y_0=0 +datum=WGS84 +units=m +no_defs 
#> source     : tmp.tif 
#> names      : SMAP_L4_S//Vv5030_001, SMAP_L4_S//Vv5030_001, SMAP_L4_S//Vv5030_001, SMAP_L4_S//Vv5030_001, SMAP_L4_S//Vv5030_001, SMAP_L4_S//Vv5030_001, SMAP_L4_S//Vv5030_001, SMAP_L4_S//Vv5030_001 
#> min values :           0.006559700,           0.006590730,           0.006612024,           0.006668247,           0.006737789,           0.006590234,           0.006586422,           0.006500987 
#> max values :                   0.8,                   0.8,                   0.8,                   0.8,                   0.8,                   0.8,                   0.8,                   0.8
```

We can visualize each layer in our RasterBrick:


```r
plot(sm_raster)
```

![plot of chunk plot-raster](vignettes/smapr-intro-plot-raster-1.png)

## Common downstream operations

### Study region cropping and masking

If you want to crop the data to a study region, you can use the `raster::crop`
function.
Let's illustrate by focusing on the state of Colorado, which is approximately
rectangular, so that we can define an `extent` object using latitude and
longitude values that roughly correspond to the state boundaries.
The `raster::extent()` function can be used with any `Spatial*` or `Raster*`
object.


```r
co_extent <- extent(c(-109, -102, 37, 41))
co_extent <- as(co_extent, "SpatialPolygons")
sp::proj4string(co_extent) <- "+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"
co_extent
#> class       : SpatialPolygons 
#> features    : 1 
#> extent      : -109, -102, 37, 41  (xmin, xmax, ymin, ymax)
#> crs         : +proj=longlat +datum=WGS84 +no_defs
```

Now that we have a SpatialPolygons object, we can use this to crop the
soil moisture raster to our study region.
First, we need to ensure that the projections are the same:


```r
proj_co_extent <- spTransform(co_extent, crs(sm_raster))
```

Then, we could crop the soil moisture data to our polygon, which will make the
extents match.


```r
co_soil_moisture <- crop(sm_raster, proj_co_extent)
plot(co_soil_moisture)
```

![plot of chunk crop-raster](vignettes/smapr-intro-crop-raster-1.png)

We might also have a polygon to use as a mask.
For example, to mask out our Colorado polygon, we could do the following
(operating on the first layer of the soil moisture raster for simplicity):


```r
plot(mask(sm_raster[[1]], proj_co_extent))
```

![plot of chunk mask-raster](vignettes/smapr-intro-mask-raster-1.png)

Notice that masking does not crop the raster, it simply sets all values outside
of the polygon to NA.
In some cases, an inverse mask is useful, where values inside the polygon are
set to NA:


```r
plot(mask(sm_raster[[1]], proj_co_extent, inverse = TRUE))
```

![plot of chunk inverse-mask](vignettes/smapr-intro-inverse-mask-1.png)

### Computing summary statistics over layers

We may want to average soil moisture values across layers of a raster brick.
This can be done with the `raster::calc()` function:


```r
mean_sm <- calc(sm_raster, fun = mean)
plot(mean_sm, main = 'Mean soil moisture')
```

![plot of chunk get-mean](vignettes/smapr-intro-get-mean-1.png)

### Comparing surface and soil moisture

Our SPL4SMAU data have estimated surface and rootzone soil moisture layers.
If we want to compare these values, we can load the surface soil moisture data,
compute the mean value over layers as we did for the rootzone soil moisture
raster, and generate a scatterplot.


```r
surface_raster <- extract_smap(local_files,
                               name = '/Analysis_Data/sm_surface_analysis')

# compute mean
mean_surface_sm <- calc(surface_raster, fun = mean)

# compare values
plot(values(mean_sm), values(mean_surface_sm), col = 'dodgerblue', cex = .1,
     xlab = 'Rootzone soil moisture', ylab = 'Surface soil moisture', bty = 'n')
abline(0, 1, lty = 2)
```

![plot of chunk surface-vs-rootzone](vignettes/smapr-intro-surface-vs-rootzone-1.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/list_smap.R
\name{list_smap}
\alias{list_smap}
\title{Lists the contents of SMAP data files}
\usage{
list_smap(files, all = TRUE)
}
\arguments{
\item{files}{A \code{data.frame} produced by \code{download_smap()} that
specifies input data files.}

\item{all}{If TRUE a longer, more detailed list of information on each
entry is provided.}
}
\value{
Returns a list of \code{data.frame} objects that list the contents
of each data file in \code{files}.
}
\description{
This function returns a list of the contents of SMAP data files.
}
\examples{
\dontrun{
files <- find_smap(id = "SPL4SMGP", dates = "2015-03-31", version = 4)
files <- download_smap(files[1, ])
list_smap(files)
list_smap(files, all = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/set_smap_credentials.R
\name{set_smap_credentials}
\alias{set_smap_credentials}
\title{Set credentials for NASA's Earthdata portal}
\usage{
set_smap_credentials(username, password, save = TRUE,
  overwrite = FALSE)
}
\arguments{
\item{username}{A character string of your Earthdata portal username}

\item{password}{A character string of your Earthdata portal password}

\item{save}{Logical: whether to save your credentials to your 
.Renviron file (e.g., ~/.Renviron). Previous Earthdata credentials will not 
be overwritten unless \code{overwrite = TRUE}.}

\item{overwrite}{Logical: whether to overwrite previous Earthdata credentials
in your .Renviron file (only applies when \code{save = TRUE})}
}
\value{
A data.frame with the names of the data files, the remote directory, and
  the date.
}
\description{
To use smapr, users need to provide NASA Earthdata portal credentials. 
This function allows users to interactively set these credentials via the 
user's Earthdata username and password.
}
\details{
If you do not yet have a username and password, register for one here:
https://urs.earthdata.nasa.gov/

A warning: do not commit your username and password to a public repository!
This function is meant to be used interactively, and not embedded within a 
script that you would share.
}
\examples{
\dontrun{
set_smap_credentials('myusername', 'mypassword')
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/extract_smap.R
\name{extract_smap}
\alias{extract_smap}
\title{Extracts contents of SMAP data}
\usage{
extract_smap(data, name, in_memory = FALSE)
}
\arguments{
\item{data}{A data frame produced by \code{download_smap()} that specifies
input files from which to extract data.}

\item{name}{The path in the HDF5 file pointing to data to extract.}

\item{in_memory}{Logical. Should the result be stored in memory? If not, then
raster objects are stored on disk in the cache directory. By default
the result is stored on disk.}
}
\value{
Returns a RasterStack object.
}
\description{
Extracts datasets from SMAP data files.
}
\details{
The arguments \code{group} and \code{dataset} must refer specifically  the
group and name within group for the input file, such as can be obtained with
\code{list_smap()}. This function will extract that particular dataset,
returning a Raster object.
}
\examples{
\dontrun{
files <- find_smap(id = "SPL4SMGP", dates = "2015-03-31", version = 4)
downloads <- download_smap(files[1, ])
sm_raster <- extract_smap(downloads, name = '/Geophysical_Data/sm_surface')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/smapr-package.R
\docType{package}
\name{smapr-package}
\alias{smapr-package}
\title{smapr: A package for acquisition and processing of NASA SMAP data.}
\description{
The smapr package provides a means to discover, acquire, and process
NASA Soil Moisture Active Passive (SMAP) data.
}
\author{
Max Joseph \email{maxwell.b.joseph@colorado.edu}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_smap.R
\name{find_smap}
\alias{find_smap}
\title{Find SMAP data}
\usage{
find_smap(id, dates, version)
}
\arguments{
\item{id}{A character string that refers to a specific SMAP dataset, e.g.,
\code{"SPL4SMGP"} for SMAP L4 Global 3-hourly 9 km Surface and Rootzone Soil
Moisture Geophysical Data. See "Details" for a list of supported data types 
and their associated id codes.}

\item{dates}{An object of class Date or a character string formatted as
%Y-%m-%d (e.g., "2016-04-01") which specifies the date(s) to search.
To search for one specific date, this can be a Date object of length one. To
search over a time interval, it can be a multi-element object of class Date
such as produced by \code{seq.Date}.}

\item{version}{Which data version would you like to search for? Version
information for each data product can be found at
\url{https://nsidc.org/data/smap/data_versions}}
}
\value{
A data.frame with the names of the data files, the remote directory, and
  the date.
}
\description{
This function searches for SMAP data on a specific date, returning a
\code{data.frame} describing available data.
}
\details{
There are many SMAP data products that can be accessed with this function.
Currently, smapr supports level 3 and level 4 data products, each of which
has an associated Data Set ID which is specified by the \code{id} argument,
described at \url{https://nsidc.org/data/smap/smap-data.html} and summarized
below:

\describe{
\item{SPL2SMAP_S}{SMAP/Sentinel-1 Radiometer/Radar Soil Moisture}
\item{SPL3FTA}{Radar Northern Hemisphere Daily Freeze/Thaw State}
\item{SPL3SMA}{Radar Global Daily Soil Moisture}
\item{SPL3SMP}{Radiometer Global Soil Moisture}
\item{SPL3SMAP}{Radar/Radiometer Global Soil Moisture}
\item{SPL4SMAU}{Surface/Rootzone Soil Moisture Analysis Update}
\item{SPL4SMGP}{Surface/Rootzone Soil Moisture Geophysical Data}
\item{SPL4SMLM}{Surface/Rootzone Soil Moisture Land Model Constants}
\item{SPL4CMDL}{Carbon Net Ecosystem Exchange}
}

This function requires a username and password from NASA's Earthdata portal.
If you have an Earthdata username and password, pass them in using the
\code{\link[=set_smap_credentials]{set_smap_credentials()}} function.

If you do not yet have a username and password, register for one here:
\url{https://urs.earthdata.nasa.gov/}
}
\examples{
\dontrun{
# looking for data on one day:
find_smap(id = "SPL4SMGP", dates = "2015-03-31", version = 4)

# searching across a date range
start_date <- as.Date("2015-03-31")
end_date <- as.Date("2015-04-02")
date_sequence <- seq(start_date, end_date, by = 1)
find_smap(id = "SPL4SMGP", dates = date_sequence, version = 4)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_smap.R
\name{download_smap}
\alias{download_smap}
\title{Download SMAP data}
\usage{
download_smap(files, directory = NULL, overwrite = TRUE,
  verbose = TRUE)
}
\arguments{
\item{files}{A \code{data.frame} produced by \code{find_smap()}
that specifies data files to download.}

\item{directory}{A local directory path in which to save data, specified as a
character string. If left as \code{NULL}, data are stored in a user's cache
directory.}

\item{overwrite}{TRUE or FALSE: should existing data files be overwritten?}

\item{verbose}{TRUE or FALSE: should messages be printed to indicate that 
files are being downloaded?}
}
\value{
Returns a \code{data.frame} that appends a column called
\code{local_dir} to the input data frame, which consists of a character
vector specifying the local directory containing the downloaded files.
}
\description{
This function downloads SMAP data in HDF5 format.
}
\details{
This function requires a username and password from NASA's Earthdata portal.
If you have an Earthdata username and password, pass them in using the
\code{\link[=set_smap_credentials]{set_smap_credentials()}} function.

If you do not yet have a username and password, register for one here:
\url{https://urs.earthdata.nasa.gov/}
}
\examples{
\dontrun{
files <- find_smap(id = "SPL4SMGP", dates = "2015-03-31", version = 4)
# files[1, ] refers to the first available data file
downloads <- download_smap(files[1, ])
}
}
