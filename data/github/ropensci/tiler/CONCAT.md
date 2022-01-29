
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tiler <img src="man/figures/logo.png" style="margin-left:10px;margin-bottom:5px;" width="120" align="right">

**Author:** [Matthew Leonawicz](https://github.com/leonawicz)
<a href="https://orcid.org/0000-0001-9452-2771" target="orcid.widget">
<img alt="ORCID logo" src="https://info.orcid.org/wp-content/uploads/2019/11/orcid_16x16.png" width="16" height="16" /></a>
<br/> **License:** [MIT](https://opensource.org/licenses/MIT)<br/>

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![Travis build
status](https://travis-ci.org/ropensci/tiler.svg?branch=master)](https://travis-ci.org/ropensci/tiler)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ropensci/tiler?branch=master&svg=true)](https://ci.appveyor.com/project/leonawicz/tiler)
[![Codecov test
coverage](https://codecov.io/gh/ropensci/tiler/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/tiler?branch=master)

[![](https://badges.ropensci.org/226_status.svg)](https://github.com/ropensci/software-review/issues/226)
[![CRAN
status](http://www.r-pkg.org/badges/version/tiler)](https://cran.r-project.org/package=tiler)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/tiler)](https://cran.r-project.org/package=tiler)
[![Github
Stars](https://img.shields.io/github/stars/ropensci/tiler.svg?style=social&label=Github)](https://github.com/ropensci/tiler)

![](https://github.com/ropensci/tiler/blob/master/data-raw/ne.jpg?raw=true)

## Create geographic and non-geographic map tiles

The `tiler` package provides a tile generator function for creating map
tile sets for use with packages such as `leaflet`. In addition to
generating map tiles based on a common raster layer source, it also
handles the non-geographic edge case, producing map tiles from arbitrary
images. These map tiles, which have a a non-geographic simple coordinate
reference system, can also be used with `leaflet` when applying the
simple CRS option.

Map tiles can be created from an input file with any of the following
extensions: `tif`, `grd` and `nc` for spatial maps and `png`, `jpg` and
`bmp` for basic images.

## Motivation

This package helps R users who wish to create geographic and
non-geographic map tiles easily and seamlessly with only a single line
of R code. The intent is to do this with a package that has

  - minimal heavy package dependencies.
  - minimal extraneous general features and functions that do not have
    to do with tile generation.
  - to create tiles without having to code directly in other software,
    interact directly with Python, or make calls at the command line;
    allowing the R user to remain comfortably within the familiar R
    environment.
  - to support the creation on map tiles from raw images for users who
    wish to create non-standard maps, which may also be followed by
    georeferencing locations of interest in the simplified coordinate
    reference system of the map image.

## Installation

Install `tiler` from CRAN with

``` r
install.packages("tiler")
```

Install the development version from GitHub with

``` r
# install.packages("remotes")
remotes::install_github("ropensci/tiler")
```

For non-geographic tiles, using a `png` file is recommended for quality
and file size. `jpg` may yield a lower quality result, while a large,
high resolution `bmp` file may have an enormous file size compared to
`png`.

`jpg` and `bmp` are *optionally* supported by `tiler`. This means they
are not installed and imported with `tiler`. It is assumed the user will
provide `png` images. If using `jpg` or `bmp` and the packages `jpeg` or
`bmp` are not installed, respectively, `tile` will print a message to
the console notifying of the required package installations.

### System requirements

This package requires Python and the `gdal` library for Python. Windows
users are recommended to install
[OSGeo4W](https://trac.osgeo.org/osgeo4w/) as an easy way to obtain the
required `gdal` support for Python in Windows. See `?tiler_options` or
the package vignette for more information.

-----

Please note that the `tiler` project is released with a [Contributor
Code of
Conduct](https://github.com/ropensci/tiler/blob/master/CODE_OF_CONDUCT.md).
By contributing to this project, you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# tiler 0.2.5

* Documentation updates.

# tiler 0.2.4

* Precompiling of vignette depending on external data.

# tiler 0.2.3

* Minor fixes for CRAN release.

# tiler 0.2.2

* Breaking change: no longer using a `format` argument. All tiles are TMS.
* Updated `gdal2tiles` to version 2.4 release.
* Bug fix.
* Updated documentation.

# tiler 0.2.1

* Improved and simplified instructions and expectations for Windows use. Windows users must add `OSGeo4W.bat` path to system path.
* Added `leaflet` examples of remotely hosted tiles generated by `tiler` to vignette.
* Bug fix related to system path to `OSGeo4W.bat` being ignored on Windows.

# tiler 0.2.0

* All three `gdal2tiles*` scripts have been updated to accept a command line argument when called by R that provides a path for any temporary files, i.e., `tmp.*.vrt` files created by the `gdal2tiles*` scripts. These were previously accumulating in the system temp folder. The new temporary directory is a sub-directory inside `tempdir()`. Therefore, it is cleaned up when exiting R. Nevertheless, `tile` also force deletes this subdirectory immediately after its internal system call to one of the `gdal2tiles*` scripts returns, so the temporary sub-directory does not even exist for the full duration of the `tile` call.
* Added functions `tile_viewer` and `view_tiles` and other supporting functions for generating HTML Leaflet tile preview web page.
* Added arguments to `tile`. `tile` now generates previewer by default.
* Added unit tests.
* Updated vignette.

# tiler 0.1.6

* Made minor formatting changes per CRAN request for resubmission.

# tiler 0.1.5

* Refactored `tile`, added arguments including `resume` and `format`, changed some argument names.
* Added default support for XYZ format tiles in addition to TMS. This brings in another version of `gdal2tiles`.
* Updated documentation.
* Added unit tests.

# tiler 0.1.0

* Created `tile` function for generating map tiles from geographic or non-geographic maps/images.
* Created readme with basic example.
* Added initial introduction vignette.
* Added robust unit tests and other external examples and spot testing of edge cases.

# tiler 0.0.0.9000

* Added initial package scaffolding.
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
## Test environments

* local Windows 10 install, R 4.0.3
* Windows 10 (AppVeyor), R 4.0.3
* Ubuntu 16.04 (Travis CI), R-devel, R-release, R-oldrel
* Mac OSX (Travis CI) R-release (4.0.3)
* win-builder

## Update release

* This update includes a maintainer email address update.

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

Checks pass on all packages. (https://github.com/ropensci/tiler/tree/master/revdep)

## System Requirements details

This package requies Python and python-gdal. These are specified in the DESCRIPTION System Requirements field as well as the Description text.
# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 3.6.1 (2019-07-05) |
|os       |Windows 10 x64               |
|system   |x86_64, mingw32              |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |English_United States.1252   |
|ctype    |English_United States.1252   |
|tz       |America/Denver               |
|date     |2019-11-25                   |

# Dependencies

|package |old   |new   |<U+0394>  |
|:-------|:-----|:-----|:--|
|tiler   |0.2.3 |0.2.4 |*  |
|png     |0.1-7 |0.1-7 |   |
|raster  |3.0-7 |3.0-7 |   |
|Rcpp    |1.0.3 |1.0.3 |   |
|rgdal   |1.4-7 |1.4-7 |   |
|sp      |1.3-2 |1.3-2 |   |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*