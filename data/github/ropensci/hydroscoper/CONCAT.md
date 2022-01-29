---
title: 'hydroscoper: R interface to the Greek National Data Bank for Hydrological and Meteorological Information'
authors:
- affiliation: 1
  name: Konstantinos Vantas
  orcid: 0000-0001-6387-8791
date: "01 Mar 2017"
output:
  html_document: default
  pdf_document: default
bibliography: paper.bib
tags:
- R
- tidy data
- Hydrology
- Meteorology
- Greece
affiliations:
- index: 1
  name: Faculty of Engineering, Aristotle University of Thessaloniki, Greece
---

# Summary

The Greek National Data bank for Hydrological and Meteorological Information, Hydroscope [@hydroscope2018], is the result of long-standing efforts by numerous Greek scientists in collaboration with various companies and associations. Its main purpose is the formation of the basic infrastructure for the implementation of the European Community Directives 2000/60/EC (i.e. establishment of a framework for Community action in the field of water policy) and 2007/60/EC, (i.e. assessment and management of flood risks).

Hydroscope offers several national data sources in HTML and plain text files via a web interface, using the Enhydris database system [@christofides2011enhydris]. These data are well structured but are in Greek, thus limiting their usefulness. Furthermore, fully reproducible research [@peng2011reproducible] can be tedious and error-prone using Hydroscope's web interface. On the contrary, using the Enhydris API for reproducibility requires external programs and scripting to import the data.

`hydroscoper` [@hydroscoper2018] provides functionality for automatic retrieval and translation of Hydroscope's data to English for use in R [@Rbase]. The main functions that can be utilized is the family of functions, `get_stations`, `get_timeseries`, `get_data`, etc., to easily download JSON and TXT files as tidy data frames [@Wickham2014]. The internal databases of the package can be used to run queries on the available stations and time series, reducing the time needed for downloading and data wrangling [@kandel2011research], as these data are rarely modified.

The data have many applications. In general, availability of meteorological and hydrological information is essential for water resources management, water quality assessment and global change studies [@vafiadis1994hydroscope]. Also, these data can be used by researchers to forward their studies when specific requirements of time series are required, such as the estimation of rainfall erosivity [@vantas_sid2017]. Finally, the data are crucial to Greek organizations for the implementation of the Water Framework Directive 2000-2027 [@WFD].


# References

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
hydroscoper
================

<!-- README.md is generated from README.Rmd. Please edit that file -->

[![Travis-CI Build
Status](https://travis-ci.org/ropensci/hydroscoper.svg?branch=master)](https://travis-ci.org/ropensci/hydroscoper)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ropensci/hydroscoper?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/hydroscoper)
[![codecov](https://codecov.io/github/ropensci/hydroscoper/branch/master/graphs/badge.svg)](https://codecov.io/gh/ropensci/hydroscoper)
[![minimal R
version](https://img.shields.io/badge/R%3E%3D-3.4-6666ff.svg)](https://cran.r-project.org/)
[![CRAN\_Status\_Badge](http://www.r-pkg.org/badges/version/hydroscoper)](https://cran.r-project.org/package=hydroscoper)
[![packageversion](https://img.shields.io/badge/Package%20version-1.4.1-orange.svg?style=flat-square)](https://github.com/ropensci/hydroscoper)
[![](https://cranlogs.r-pkg.org/badges/grand-total/hydroscoper)](https://cran.r-project.org/package=hydroscoper)
[![ropensci](https://badges.ropensci.org/185_status.svg)](https://github.com/ropensci/software-review/issues/185)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1196540.svg)](https://doi.org/10.5281/zenodo.1196540)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00625/status.svg)](https://doi.org/10.21105/joss.00625)

<img src="https://github.com/ropensci/hydroscoper/raw/master/man/figures/hydroscoper_hex.png" align = "right" width = 120/>

`hydroscoper` is an R interface to the Greek National Data Bank for
Hydrological and Meteorological Information, *Hydroscope*. For more
details checkout the package’s
[website](https://docs.ropensci.org/hydroscoper/) and the vignettes:

-   [An introduction to
    `hydroscoper`](https://docs.ropensci.org/hydroscoper/articles/intro_hydroscoper.html)
    with details about the Hydroscope project and the package.
-   [Using `hydroscoper`’s data
    sets](https://docs.ropensci.org/hydroscoper/articles/stations_with_data.html)
    with a simple example of how to use the package’s internal data
    sets.

## Installation

Install the stable release from CRAN with:

``` r
install.packages("hydroscoper")
```

You can install the development version from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("ropensci/hydroscoper")
```

## Using hydroscoper

The functions that are provided by `hydroscoper` are:

-   `get_stations, get_timeseries, ..., etc.` family functions, to
    retrieve tibbles with Hydroscope’s data for a given data source.
-   `get_data`, to retrieve a tibble with time series’ values.  
-   `hydro_coords`, to convert Hydroscope’s points’ raw format to a
    tibble.
-   `hydro_translate` to translate various terms and names from Greek to
    English.

The data sets that are provided by `hydroscoper` are:

-   `stations` a tibble with stations’ data from Hydroscope.
-   `timeseries` a tibble with time series’ data from Hydroscope.
-   `greece_borders` a tibble with the borders of Greece.

## Example

This is a minimal example which shows how to get the station’s *200200*
precipitation time series *56* from the *kyy* sub-domain.

Load libraries and get data:

``` r
library(hydroscoper)
library(tibble)
library(ggplot2)

ts_raw <- get_data(subdomain = "kyy", time_id = 56)
ts_raw
#> # A tibble: 147,519 x 3
#>    date                value comment
#>    <dttm>              <dbl> <chr>  
#>  1 1985-05-06 08:00:00     0 1      
#>  2 1985-05-06 08:30:00     0 1      
#>  3 1985-05-06 09:00:00     0 1      
#>  4 1985-05-06 09:30:00     0 1      
#>  5 1985-05-06 10:00:00     0 1      
#>  6 1985-05-06 10:30:00     0 1      
#>  7 1985-05-06 11:00:00     0 1      
#>  8 1985-05-06 11:30:00     0 1      
#>  9 1985-05-06 12:00:00     0 1      
#> 10 1985-05-06 12:30:00     0 1      
#> # … with 147,509 more rows
```

Let’s create a plot:

``` r
ggplot(data = ts_raw, aes(x = date, y = value))+
  geom_line()+
  labs(title= "30 min precipitation for station 200200",
       x="Date", y = "Rain height (mm)")+
  theme_classic()
```

![](man/figures/README-plot_time_series-1.png)<!-- -->

## Meta

-   Bug reports, suggestions, and code are welcome. Please see
    [Contributing](https://github.com/ropensci/hydroscoper/blob/master/CONTRIBUTING.md).
-   License:
    -   All code is licensed MIT.
    -   All data are from the public data sources in
        `http://www.hydroscope.gr/`.
-   To cite `hydroscoper` please use:

<!-- -->

    Vantas Konstantinos, (2018). hydroscoper: R interface to the Greek National Data Bank for
    Hydrological and Meteorological Information. Journal of Open Source Software,
    3(23), 625 DOI:10.21105/joss.00625

or the BibTeX entry:

    @Article{kvantas2018,
    author = {Konstantinos Vantas},
    title = {{hydroscoper}: R interface to the Greek National Data Bank for Hydrological and Meteorological Information},
    doi = {10.21105/joss.00625},
    year = {2018},
    month = {mar},
    publisher = {The Open Journal},
    volume = {2},
    number = {23},
    journal = {The Journal of Open Source Software}
    }

[![ropensci\_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# hydroscoper 1.4.1 (Release date: 2021-05-14)

* remove call to closeAllConnections() command

# hydroscoper 1.4 (Release date: 2021-03-20)

* Fail gracefully with an informative message if the hydroscope site is not 
  available or has changed
* Update dependencies

# hydroscoper 1.3 (Release date: 2020-07-03)

* Update dependencies

# hydroscoper 1.2 (Release date: 2019-06-02)

* This version removes plyr dependency.

# hydroscoper 1.1.1 (Release date: 2018-08-25)

* This version fixes the errors found in CRAN Package Check Results.

# hydroscoper 1.1.0 (Release date: 2018-07-06)

* New functionality:

  - `find_stations()` returns a tibble with the nearest hydroscope's stations' distances using a given point coordinates.


# hydroscoper 1.0.0 (Release date: 2018-03-16)

* `hydroscoper` was transferred to rOpenSci: https://github.com/ropensci/hydroscoper

* General

  - This is a major update. All the functions are rewritten utilizing the Enhydris API.
  - The included data in the package cover all Hydroscope's databases.
  - Add vignettes "An introduction to `hydroscoper`" and  "Using `hydroscoper`'s data".
  - Add package documentation site.
  - Use pingr package to check if a sub-domain is alive.
  - All the functions return tibbles.
  - Add `greece_borders` dataset.

* New functionality:

  - `get_instruments()` returns a tibble with the instruments' data.
  - `get_water_basins()` returns a tibble with the Water Basins' data.
  - `get_water_divisions()` returns a tibble with the Water Divisions' data.
  - `get_political_divisions()` returns a tibble with the Political Divisions' data.
  - `get_variables()` returns a tibble with the Variables' data.
  - `get_units_of_measurement()` returns a tibble with the Units' data.
  - `get_time_steps()` returns a tibble with the Time Steps' data.
  - `get_owners()` returns a tibble with the Owners' data.
  - `get_instruments_type()`returns a tibble with the Instruments' type data.
  - `get_station_type()` returns a tibble with the Water Basins data.
  - `get_database()` returns a tibble with the Water Basins data.
  - `hydro_coords` returns a tibble with the stations' longitudes and latitudes using as input the variable `point` from `get_stations` function.
  - `hydro_translate()` translates various Greek terms to English.

* Changes

  - `get_stations` and `get_timeseries` use the Enhydris API and are considerably faster.
  - `get_data` uses lower case variable naming: `date, value, comment`

* Defuncs

  - `get_coords` has been removed from the package. Please use `hydro_coords` to convert Hydroscope's points' raw format to a tibble.

--------------------------------------------------------------------------------

# hydroscoper 0.1.0 (Release date: 2017-12-22)

* Initial submission to CRAN.



# Suggestions and bug reports

`hydroscoper` thrives on the suggestions and bug reports you submit to the [issue tracker](https://github.com/ropensci/hydroscoper/issues). Before posting, please search both the open and closed issues to help us avoid duplication.

# Patches to code and documentation

Please help us develop `hydroscoper`. It really helps us manage the exciting but accelerating stream of incoming [issues](https://github.com/ropensci/hydroscoper/issues). To submit a patch, please [fork this repository](https://help.github.com/articles/fork-a-repo/) and submit a [pull request](https://github.com/ropensci/hydroscoper/pulls).

# Best ways to help

Multiple [issues](https://github.com/ropensci/hydroscoper) are labeled ["help wanted"](https://github.com/ropensci/hydroscoper/issues?q=is%3Aissue+is%3Aopen+label%3A%22help+wanted%22). We would particularly appreciate your input there.

# Contributing

Maintainers and other contributors must follow this repository's [code of conduct](https://github.com/ropensci/hydroscoper/blob/master/CONDUCT.md).
## Test environments
* local Mac OS X install, R 4.0.4
* travis.ci Linux, x64, R 4.0.2
* ci.appveyor, R 4.0.4
* win-builder (devel and release)
* check-rhub

## R CMD check results
There were no ERRORs or WARNINGs.

## Downstream dependencies
There are currently no downstream dependencies for this package
