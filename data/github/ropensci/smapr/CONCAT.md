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
