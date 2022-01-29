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

PostcodesioR
================

# PostcodesioR <img src='man/figures/logo.png' align="right" height="139" />

[![Travis-CI Build
Status](https://travis-ci.org/ropensci/PostcodesioR.svg?branch=master)](https://travis-ci.org/ropensci/PostcodesioR)
[![Package-License](https://img.shields.io/badge/license-GPL--3-brightgreen.svg?style=flat)](https://www.gnu.org/licenses/gpl-3.0.html)
[![](https://badges.ropensci.org/176_status.svg)](https://github.com/ropensci/software-review/issues/176)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/PostcodesioR)](https://cran.r-project.org/package=PostcodesioR)
[![DOI](https://zenodo.org/badge/64221541.svg)](https://zenodo.org/badge/latestdoi/64221541)
![Downloads](https://cranlogs.r-pkg.org/badges/grand-total/PostcodesioR)

An API wrapper around [postcodes.io](https://postcodes.io/) - free UK
postcode lookup and geocoder. This package helps to find and transform
information about UK administrative geography like postcodes, LSOA,
MSOA, constituencies, counties, wards, districts, CCG or NUTS.

The package is based exclusively on open data provided by postcodes.io.
PostcodesioR can be used by data scientists or social scientists working
with geocoded UK data. A common task when working with such data is
aggregating geocoded data on different administrative levels,
e.g. turning postcode-level data into counties or regions. This package
can help in achieving this and in many other cases when changing the
aggregation of geographic data is required.

## Installation

This package can be installed from GitHub (developmental version) or
CRAN (stable).

In order to install PostcodesioR use one of the following commands:

``` r
# stable version
install.packages("PostcodesioR")
```

or

``` r
# developmental version
if(!require("devtools")) {
  install.packages("devtools")
}
devtools::install_github("ropensci/PostcodesioR")
```

## Loading

Load the package by typing

``` r
library(PostcodesioR)
```

## Examples

Where possible, I tried to return a data frame. Unfortunately, a lot of
API calls return more complex data and in those cases it is safer to use
lists. The API limits the number of returned calls. Check functions’
documentation for more details.

For additional information about the returned data and the function
calls see the original [documentation](https://postcodes.io/docs).

The main function of this package provides information related to a
given postcode

``` r
lookup_result <- postcode_lookup("EC1Y8LX")

#overview
str(lookup_result)
```

    ## 'data.frame':    1 obs. of  35 variables:
    ##  $ postcode                       : chr "EC1Y 8LX"
    ##  $ quality                        : int 1
    ##  $ eastings                       : int 532544
    ##  $ northings                      : int 182128
    ##  $ country                        : chr "England"
    ##  $ nhs_ha                         : chr "London"
    ##  $ longitude                      : num -0.0909
    ##  $ latitude                       : num 51.5
    ##  $ european_electoral_region      : chr "London"
    ##  $ primary_care_trust             : chr "Islington"
    ##  $ region                         : chr "London"
    ##  $ lsoa                           : chr "Islington 023D"
    ##  $ msoa                           : chr "Islington 023"
    ##  $ incode                         : chr "8LX"
    ##  $ outcode                        : chr "EC1Y"
    ##  $ parliamentary_constituency     : chr "Islington South and Finsbury"
    ##  $ admin_district                 : chr "Islington"
    ##  $ parish                         : chr "Islington, unparished area"
    ##  $ admin_county                   : logi NA
    ##  $ admin_ward                     : chr "Bunhill"
    ##  $ ced                            : logi NA
    ##  $ ccg                            : chr "NHS North Central London"
    ##  $ nuts                           : chr "Haringey and Islington"
    ##  $ admin_district_code            : chr "E09000019"
    ##  $ admin_county_code              : chr "E99999999"
    ##  $ admin_ward_code                : chr "E05000367"
    ##  $ parish_code                    : chr "E43000209"
    ##  $ parliamentary_constituency_code: chr "E14000764"
    ##  $ ccg_code                       : chr "E38000240"
    ##  $ ccg_id_code                    : chr "93C"
    ##  $ ced_code                       : chr "E99999999"
    ##  $ nuts_code                      : chr "TLI43"
    ##  $ lsoa_code                      : chr "E01002704"
    ##  $ msoa_code                      : chr "E02000576"
    ##  $ lau2_code                      : chr "E09000019"

Check the
[vignette](https://docs.ropensci.org/PostcodesioR/articles/Introduction.html)
to see all functions in action.

## Notes

Currently, there is a limit to the number of API calls that can be made.
However, [postcodes.io](https://postcodes.io/) provides full list of
geolocation data that can be used locally without limitations. The
original data is sourced from [Office for National Statistics Data
Portal](https://geoportal.statistics.gov.uk/). That
[file](https://github.com/ideal-postcodes/postcodes.io/blob/master/latest)
is rather large so I didn’t include it in the package.

Go to the package’s [website](https://docs.ropensci.org/PostcodesioR/)
or to my [blog](https://walczak.org/tag/postcodesior/) for more
examples.

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/PostcodesioR/blob/master/CONDUCT.md).
By participating in this project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# PostcodesioR 0.3.1

# PostcodesioR 0.3.0

* `scottish_postcode_lookup` added.
* New fields (codes) added.
* README updated (hex logo and downloads).
* New tests.

# PostcodesioR 0.2.0

* `bulk_postcode_lookup` bug fixed.

# PostcodesioR 0.1.1

* Added a `NEWS.md` file to track changes to the package.
## Update (v 0.3.1)

* Minor bus fixes.
* Updated tests.
* Updated documentation.

## Test environments
* local ubuntu 20.04, R 4.0.3
* win-builder (devel and release). No errors, warnings or notes
* R-hub using `rhub::check_for_cran()`. One note related to the change of maintainer's email.
* R-hub tested on:
- Windows Server 2008 R2 SP1, R-devel, 32/64 bit,
- Ubuntu Linux 16.04 LTS, R-release, GCC
- Fedora Linux, R-devel, clang, gfortran

## Update (v 0.3)

* Minor bus fixes.
* `scottish_postcode_lookup` added
* New fields (codes) added
* README updated (hex logo and downloads)
* New tests
* URLs in README updated

## Test environments
* local ubuntu 20.04, R 4.0.3
* win-builder (devel and release). No errors, warnings or notes
* R-hub using `rhub::check_for_cran()`. No errors, warnings or notes
* R-hub tested on:
- Windows Server 2008 R2 SP1, R-devel, 32/64 bit,
- Ubuntu Linux 16.04 LTS, R-release, GCC
- Fedora Linux, R-devel, clang, gfortran

## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

This is a resubmission. In this version I have:

* Wrapped Postcodes.io in DESCRIPTION Title in quotation marks

* Shortened the title

* Replaced `dontrun` with `donttest` in examples

## Test environments
* local ubuntu 16.04, R 3.6.1
* ubuntu 14.04 (on travis-ci), R 3.6.0
* win-builder (devel and release)
* R-hub using `rhub::check_for_cran()`

## R CMD check results

0 errors | 0 warnings | 0 notes

## `rhub::check_for_cran()` results

* Error

Installation failed with PREPERROR on R-hub's Fedora Linux, R-devel, clang, gfortran (and other platforms), on account of failures to install required dependencies. This seems to be an issue outside of my control.

* Note

Windows Server 2008 R2 SP1, R-devel, 32/64 bit and Ubuntu Linux 16.04 LTS, R-release, GCC return a note:

The note regarded potential mis-spelled words. These words were spelled correctly.

Possibly mis-spelled words in DESCRIPTION:
  geocoding (16:22, 18:20)
  io (3:37)

* This is a new release.
