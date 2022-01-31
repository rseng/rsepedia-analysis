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
---
title: hydroscoper
output: github_document
editor_options: 
  chunk_output_type: inline
---
<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo=FALSE, message=FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
)

chk_online <- FALSE

library(pingr)

# helper function to check if a sub-domain is online
online <- function(h_url){
  !is.na(pingr::ping(h_url, count = 1))
}

# check if sub-domains are online
chk_online   <- online("kyy.hydroscope.gr")

```
 
[![Travis-CI Build Status](https://travis-ci.org/ropensci/hydroscoper.svg?branch=master)](https://travis-ci.org/ropensci/hydroscoper)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/hydroscoper?branch=master&svg=true)](https://ci.appveyor.com/project/ropensci/hydroscoper)
[![codecov](https://codecov.io/github/ropensci/hydroscoper/branch/master/graphs/badge.svg)](https://codecov.io/gh/ropensci/hydroscoper) 
[![minimal R version](https://img.shields.io/badge/R%3E%3D-3.4-6666ff.svg)](https://cran.r-project.org/)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/hydroscoper)](https://cran.r-project.org/package=hydroscoper)
[![packageversion](https://img.shields.io/badge/Package%20version-1.4.1-orange.svg?style=flat-square)](https://github.com/ropensci/hydroscoper)
[![](https://cranlogs.r-pkg.org/badges/grand-total/hydroscoper)](https://cran.r-project.org/package=hydroscoper)
[![ropensci](https://badges.ropensci.org/185_status.svg)](https://github.com/ropensci/software-review/issues/185)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1196540.svg)](https://doi.org/10.5281/zenodo.1196540)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00625/status.svg)](https://doi.org/10.21105/joss.00625)

<img src="https://github.com/ropensci/hydroscoper/raw/master/man/figures/hydroscoper_hex.png" align = "right" width = 120/>

`hydroscoper` is an R interface to the  Greek National Data Bank for Hydrological and Meteorological Information,
*Hydroscope*. For more details checkout the package's [website](https://docs.ropensci.org/hydroscoper/) and the vignettes:

 * [An introduction to `hydroscoper`](https://docs.ropensci.org/hydroscoper/articles/intro_hydroscoper.html) with details about the Hydroscope project and the package.
 * [Using `hydroscoper`'s data sets](https://docs.ropensci.org/hydroscoper/articles/stations_with_data.html) with a simple example of how to use the package's internal data sets.
 
## Installation

Install the stable release from CRAN with:

```{r cran_installation, eval = FALSE}
install.packages("hydroscoper")
```

You can install the development version from GitHub with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("ropensci/hydroscoper")
```

## Using hydroscoper

The functions that are provided by `hydroscoper` are:

* `get_stations, get_timeseries, ..., etc.` family functions, to retrieve tibbles with Hydroscope's data for a given data source.
* `get_data`, to retrieve a tibble with time series' values.  
* `hydro_coords`, to convert Hydroscope's points' raw format to a tibble.
* `hydro_translate` to translate various terms and names from Greek to English.

The data sets that are provided by `hydroscoper` are:

* `stations` a tibble with stations' data from Hydroscope.
* `timeseries` a tibble with  time series' data from Hydroscope.
* `greece_borders` a tibble with the borders of Greece.

## Example

This is a minimal example which shows how to get the station's *200200* precipitation time series *56*  from the *kyy* sub-domain.

Load libraries and get data:

```{r load_libraries, eval = chk_online}
library(hydroscoper)
library(tibble)
library(ggplot2)

ts_raw <- get_data(subdomain = "kyy", time_id = 56)
ts_raw
```

Let's create a plot:

```{r plot_time_series, eval = chk_online}
ggplot(data = ts_raw, aes(x = date, y = value))+
  geom_line()+
  labs(title= "30 min precipitation for station 200200",
       x="Date", y = "Rain height (mm)")+
  theme_classic()
```

## Meta

* Bug reports, suggestions, and code are welcome. Please see [Contributing](https://github.com/ropensci/hydroscoper/blob/master/CONTRIBUTING.md).
* License:
    + All code is licensed MIT.
    + All data are from the public data sources in `http://www.hydroscope.gr/`.
* To cite `hydroscoper` please use:
```
Vantas Konstantinos, (2018). hydroscoper: R interface to the Greek National Data Bank for
Hydrological and Meteorological Information. Journal of Open Source Software,
3(23), 625 DOI:10.21105/joss.00625
```
or the BibTeX entry:
```
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
```

[![ropensci_footer](http://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
---
title: "Getting Hydroscope's data"
author: "Konstantinos Vantas"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Vignette Title}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  eval = TRUE,
  comment = "#>"
)

chk_online <- FALSE
chk_stations <- FALSE
chk_ts <- FALSE

```

## Introduction

This vignette shows how to create a tidy data frame for the stations' and time series' data from Hydroscope. We will use the Hydroscope's sub-domains `kyy`, `ypaat` and `emy`, because their servers are maintened by the National Technical University Of Athens and work seamlessly.

This vignette requires an internet connection to run and that the Hydroscope's subdomains to be online. We can check if the sub-domains are online using the `pingr` package:

```{r}
library(pingr)

# helper function to check if a sub-domain is online
online <- function(h_url){
  !is.na(pingr::ping(h_url, count = 1))
}

# check if sub-domains are online
kyy   <- online("kyy.hydroscope.gr")
emy   <- online("emy.hydroscope.gr")
ypaat <- online("ypaat.hydroscope.gr")
chk_online <- kyy & emy & ypaat
print(paste("All sub-domains are online:", chk_online))
```

Note that the following chucks will be evaluated only if all sub-domains are online and if the retrieved data have some expected variables. Downloading will take some time, depending on your internet connection and the Hydroscope's servers traffic.

## Download data

With the following code we can download Hydroscope's sub-domains' data to a named list.

```{r download_all, eval = chk_online}
library(hydroscoper)
library(tibble)
library(plyr)

subdomain <- c("kyy", "ypaat", "emy")

# download all databases
hydro_data <- lapply(subdomain, function(x) get_database(x))
names(hydro_data) <- subdomain

```

## Stations

First of all, we must check if there are specific variables at the sub-domain's stations' data, because we will merge these data to a new data frame.

```{r, check_stations, eval = chk_online}

# check if all expected variables exist in all subdomains stations data
vars <- c("id", "name", "owner", "point", "altitude", "water_basin",
          "water_division")
chk1 <- laply(subdomain, function(id) {
    all(vars %in% names(hydro_data[[id]]$stations))
  })

# check water_basin, water_division and owner variables
vars2 <- c("id", "name")
chk2 <- laply(subdomain, function(id) { 
  all((vars2 %in% names(hydro_data[[id]]$water_basin)) &
        (vars2 %in% names(hydro_data[[id]]$water_division)) &
        (vars2 %in% names(hydro_data[[id]]$owner))
  )})

chk_stations <- all(chk1) & all(chk2)
print(paste("All expected variables exist in stations' related data:",
            chk_stations))
```

With the following code we will create a tidy data frame with the stations data. Note that the Greek terms are translated to English.

```{r, merge_stations_data, eval = chk_stations}
stations <- ldply(subdomain, function(id){

  tmp <-  hydro_data[[id]]

  # extract dataframes to join
  wbas <- tmp$water_basins[c("id", "name")]
  wdiv <- tmp$water_divisions[c("id", "name")]
  owners <- tmp$owners[c("id", "name")]

  # remove area from water_basin
  wbas$name <- laply(wbas$name, function(y) gsub("\\([^()]*\\)", "", y))

  # translate names
  wdiv$name <- hydro_translate(wdiv$name, "division")
  owners$name <- hydro_translate(owners$name, "owner")

  # rename dataframes variables
  names(wbas) <- c("water_basin", "water_basin_name")
  names(wdiv) <- c("water_division", "water_division_name")
  names(owners) <- c("owner", "owner_name")

  # merge data
  res <- merge(tmp$stations, wbas, by = "water_basin", all.x = TRUE)
  res <- merge(res, wdiv, by = "water_division", all.x = TRUE)
  res <- merge(res, owners, by = "owner", all.x = TRUE)

  # create coords
  coords <- hydro_coords(res$point)

  # create a data frame with all the data
  data.frame(station_id = as.integer(res$id),
             name = as.character(res$name),
             water_basin = res$water_basin_name,
             water_division = res$water_division_name,
             owner = res$owner_name,
             longitude = as.numeric(coords$long),
             latitude = as.numeric(coords$lat),
             altitude = as.numeric(res$altitude),
             subdomain = rep(id, nrow(res)),
             stringsAsFactors = FALSE)

})
stations <- as_tibble(stations)
stations
```

## Time series

Similarly, we will check if there are specific variables at each sub-domain's time series' data.

```{r check_timeseries_data, eval = chk_online}

# check if variables exists in timeseries data
vart <- c("id", "gentity", "start_date_utc", "end_date_utc", "variable",
          "time_step", "unit_of_measurement")
chk3 <- laply(subdomain, function(id) { 
  all(vart %in% names(hydro_data[[id]]$timeseries))
  })

# check if variables exists in hydrometeorological variables and timesteps
vart2 <- c("id", "descr")
chk4 <- laply(subdomain, function(id) {
  all(vart2 %in% names(hydro_data[[id]]$variables)) &
  all(vart2 %in% names(hydro_data[[id]]$time_steps))
})

# check if variables exists in unit_of_measurement
vart3 <- c("id", "symbol")
chk5 <- laply(subdomain, function(id) {
  all(vart3 %in% names(hydro_data[[id]]$units))
})

chk_ts <- all(chk3) & all(chk4) & all(chk5)
print(paste("All expected variables exist in time series' data:",
            chk_ts))


```

With the following code we will, also, create a tidy data frame with the time series data with translated terms.

```{r merge_timeseries_data, eval = chk_ts}
# create stations' dataframe
timeseries <- ldply(subdomain, function(id){
  
  tmp <- hydro_data[[id]]
  
  # create data frames to join
  var_names <- tmp$variables[c("id", "descr")]
  ts_names <-  tmp$time_steps[c("id", "descr")]
  ts_units <-tmp$units[c("id", "symbol")]
  
  # translate names
  var_names$descr <- hydro_translate(var_names$descr, "variable")
  ts_names$descr <- hydro_translate(ts_names$descr, "timestep")
  
  # rename dataframes variables
  names(var_names) <- c("variable", "variable_name")
  names(ts_names) <- c("time_step", "time_step_name")
  names(ts_units) <- c("unit_of_measurement", "unit_of_measurement_name")
  
  # merge data
  res <- merge(tmp$timeseries, var_names, by = "variable", all.x = TRUE)
  res <- merge(res, ts_names, by = "time_step", all.x = TRUE)
  res <- merge(res, ts_units, by = "unit_of_measurement", all.x = TRUE)
  
  data.frame(
    "time_id" = as.integer(res$id),
    "station_id" = as.integer(res$gentity),
    "variable" = as.character(res$variable_name),
    "timestep" = as.character(res$time_step_name),
    "units" = as.character(res$unit_of_measurement_name),
    "start_date" = as.character(res$start_date_utc),
    "end_date" = as.character(res$end_date_utc),
    "subdomain" = rep(id, nrow(res)), 
    stringsAsFactors = FALSE)
})
timeseries <- as_tibble(timeseries)
timeseries
```

---
title: "Using `hydroscoper`'s data"
author: "Konstantinos Vantas"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Using hydroscoper's data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---
 
```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

This vignette shows how to use the package's internal data sets.

## Load libraries

```{r load_library}
library(hydroscoper)
library(tibble)
library(ggplot2)
```

## Data sets

There are three data sets stored in the package. `stations` is comprised of the stations' id, name, longitude, latitude, etc.

```{r stations_data}
stations
```

`timeseries` of the time series' id, the corresponding station, variable type, time step etc.

```{r timeseries_data}
timeseries
```

`greece_borders` is a data-frame for use with the function `geom_polygon` from the `ggplot2` package.

## Stations location

`stations` and `greece_borders` can be used to create a map with all Hydroscope's stations. Unfortunately, there is a number of them that have erroneous coordinates (over the sea and far from Greece). Also, there are 120 stations with missing coordinates.

```{r all_stations}
ggplot() + 
  geom_polygon(data = greece_borders,
               aes(long, lat, group = group),
               fill = "grey",
               color = NA) +
  geom_point(data = stations,
             aes(x = longitude, y = latitude, color = subdomain)) +
  scale_color_manual(values=c("#E64B35FF", "#4DBBD5FF", "#00A087FF", 
                              "#3C5488FF"))+
  coord_fixed(ratio=1) +
  theme_bw()
```

## Stations with available time series 

The location of the stations with time series available to download are presented at the following map.

```{r stations_with_timeseries}
stations_ts <- subset(stations, station_id %in% timeseries$station_id &
                        subdomain %in% c("kyy", "ypaat"))


ggplot() + 
  geom_polygon(data = greece_borders,
               aes(long, lat, group = group),
               fill = "grey",
               color = NA) +
  geom_point(data = stations_ts,
             aes(x = longitude, y = latitude, color = subdomain)) +
  scale_color_manual(values=c("#00A087FF", "#3C5488FF"))+
  coord_fixed(ratio=1) +
  theme_bw()
```

Although there is a large number of stations with available data, there is heterogeneity in the coverage of the country.
---
title: "An introduction to `hydroscoper`"
author: "Konstantinos Vantas"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{An introduction to hydroscoper}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

chk_online <- FALSE

library(pingr)

# helper function to check if a sub-domain is online
online <- function(h_url){
  !is.na(pingr::ping(h_url, count = 1))
}

# check if sub-domains are online
chk_online   <- online("kyy.hydroscope.gr")

```

## What is Hydroscope?

`Hydroscope` is the  Greek National Data Bank for Hydrological and Meteorological Information, a result of long-standing efforts by numerous Greek scientists in collaboration with various companies and associations. It was implemented in three phases, funded by the Ministry of Development, the Ministry of Environment and Energy and the European Union.

This National Data Bank provides several data sources from various organisations via a web interface. Each participating organisation keeps its data on its own server using a database system for the storage and management of information. These organisations are:

* Ministry of Environment and Energy.
* Ministry of Rural Development and Food.
* National Meteorological Service.
* National Observatory of Athens.
* Greek Prefectures.
* Public Power Corporation.

The above data are structured as tables and space separated text files, but are in Greek, thus limiting their usefulness. Another issue with Hydroscope is the lack of comprehensive look-up tables about the available data, which are spread across many different databases.  

## What does `hydroscoper`?

`hydroscoper` provides functionality for automatic retrieval and translation of Hydroscope's data to English.  The main functions that can be utilized is the family of functions, `get_stations`, `get_timeseries`, `get_data`, etc., to easily download Hydroscope's data as `tibbles`. 

The package covers Hydroscope's data sources using the Enhydris API. The Enhydris database is implemented in PostgreSQL and details about the  about the Web-service API [here](http://enhydris.readthedocs.io).

### Internal datasets

The internal datasets of the package can be used to run queries on the available Hydroscope's stations and time series data, reducing the time needed for downloading and data wrangling, as these data are rarely modified.  These datasets are:

#### `stations`

It is a comprehensive look-up table with geographical and ownership information of the available stations in all Hydroscope's databases. The variables are: 

  1. `station_id` The station's ID.
  2. `name` The station's name.
  3. `water_basin` The station's Water Basin.
  4. `water_division` The station's Water Division.
  5. `owner` The station's owner.
  6. `longitude` The station's longitude in decimal degrees, (ETRS89).
  7. `latitude` The station's latitude in decimal degrees, (ETRS89).
  8. `altitude` The station's altitude, meters above sea level.
  9. `subdomain` The corresponding Hydroscope's database.


#### `timeseries`

It is also a look-up table with all the available measurements for a given station in a  given Hydroscope's  database, with units of measurement and times of those measurements. The variables are:

  1. `time_id` The time series ID.
  2. `station_id` The corresponding station's ID.
  3. `variable` The time series variable type.
  4. `timestep` The timestep of time series.
  5. `units` The units of the time series.
  6. `start_date` The starting date of time series values.
  7. `end_date` The ending date of time series values.
  8. `subdomain` The corresponding Hydroscope's database.


## Data sources

 * Ministry of Environment and Energy, National Observatory of Athens and Greek Prefectures, `http://kyy.hydroscope.gr/`.
 * Ministry of Rural Development and Food, `http://ypaat.hydroscope.gr`.
 * National Meteorological Service, `http://emy.hydroscope.gr`.
 * Greek Public Power Corporation, `http://deh.hydroscope.gr`.

Note that:

1. Only the two Ministries allow to download time series values freely.
2. `ypaat`, `emy` and `kyy` sub-domains are maintained  by the National Technical University Of Athens and these servers work seamlessly.
3. `deh` sub-domain is maintained  by the Greek Public Power Corporation and occasionally the server is down.

## Example

This is a basic example which shows how to get the stations' and time series' data from the Hydroscope's Ministry of Environment and Energy database, `http://kyy.hydroscope.gr/`.

Load libraries:

```{r load_libraries}
library(hydroscoper)
library(ggplot2)
library(tibble)

```

We will use the package's data `stations` and `timeseries`, to reduce the time needed with data munging. We can subset the station's data for the `kyy` sub-domain with:

```{r subset_data}
# load data
data("stations")

# subset stations data
kyy_stations <- subset(stations, subdomain == "kyy")

# view kyy stations
kyy_stations
```

Let's plot these stations using the package's dataset `greece_borders`.

```{r kyy_stations_map}
ggplot() + 
  geom_polygon(data = greece_borders,
               aes(long, lat, group = group),
               fill = "grey",
               color = NA) +
  geom_point(data = kyy_stations,
             aes(x = longitude, y = latitude),
             color = "#E64B35FF") +
  coord_fixed(ratio=1) +
  theme_bw()
```


To get the time series' data for the station `200200`  we can use:

```{r subset_timeseries, eval = chk_online}
station_ts <- subset(timeseries, station_id == 200200)
station_ts
```

We can download the station's  precipitation time series **56**:

``` {r get_timeseries, eval = chk_online}
ts_raw <- get_data(subdomain = "kyy", time_id = 56)
ts_raw
```

Let's create a plot:

``` {r plot_time_series, eval = chk_online}
ggplot(data = ts_raw, aes(x = date, y = value))+
  geom_line()+
  labs(title= "30 min precipitation", 
       subtitle = "station 200200",
       x="Date", y = "Rain height (mm)")+
  theme_classic()
```



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{greece_borders}
\alias{greece_borders}
\title{Greek borders}
\format{
A tibble with 18,474 rows and 8 variables:
\describe{
  \item{long}{Longitude in decimal degrees, ETRS89}
  \item{lat}{Latitude in decimal degrees, ETRS89}
  \item{order}{order, integer}
  \item{hole}{hole, boolean}
  \item{piece}{piece, integer}
  \item{group}{group, numeric}
}
}
\source{
Konstantinos Vantas
}
\usage{
greece_borders
}
\description{
The borders of Greece are taken from Geoadata.gov.gr.  The
variables are created using the function tidy from the broom package. This
data frame was created for use with the geom_polygon from ggplot2 package.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/find_stations.R
\name{find_stations}
\alias{find_stations}
\title{Find nearest stations using a point's coordinates}
\usage{
find_stations(longitude = 24, latitude = 38)
}
\arguments{
\item{longitude}{a numeric value in degrees}

\item{latitude}{a numeric value in degrees}
}
\value{
If the given longitude is in [24, 38] and the latitude is in [34, 42]
(i.e. are valid values for Greece) returns an ordered tibble with the
station_id, name, subdomain and distance values in km. The station's data
that are used  come from the `stations` dataset. Otherwise returns an error
message.
}
\description{
\code{find_stations} returns a tibble with the nearest stations'
distances using a given point's longitude and latitude values. This function
uses the Haversine formula for distance calculation in km.
}
\examples{

# find the five nearest stations to a point near Thessaloniki,
# (lon, lat) = (22.97, 40.60)
head(find_stations(22.97, 40.60), 5)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hydro_translate.R
\name{hydro_translate}
\alias{hydro_translate}
\title{Translate Greek names and terms to English}
\usage{
hydro_translate(x, value = c("owner", "variable", "timestep", "division"))
}
\arguments{
\item{x}{a string vector}

\item{value}{One of the predefined values in
\code{c("owner", "variable", "timestep", "division")}}
}
\value{
If \code{value} is one of:
\itemize{
\item{\code{owner}, organizations' names.}
\item{\code{variable}, hydrometeorological term.}
\item{\code{timestep}, timestep term.}
\item{\code{division}, Water Division.}
}
returns a character vector with translations of various hydrometeorological
terms or organizations' names from Greek (with latin characters) to English.

The organizations' names in \code{owner} are:

\tabular{ll}{
\strong{Code}       \tab \strong{Name} \cr
min_envir_energy    \tab Ministry of Environment and Energy \cr
min_agricult        \tab Ministry of Rural Development and Food\cr
natio_meteo_service \tab National Meteorological Service  \cr
natio_observ_athens \tab National Observatory of Athens \cr
public_power_corp   \tab Public Power Corporation \cr
natio_argic_resear  \tab National Agricultural Research Foundation \cr
greek_perfectures   \tab Greek Prefectures \cr
crete_eng_faculty   \tab Technical University of Crete \cr
crete_natural_museum \tab Natural History Museum of Crete \cr
}


The Greek Water Divisions codes in \code{division} are:

\tabular{ll}{
\strong{Code}  \tab \strong{Name} \cr
GR01 \tab Dytike Peloponnesos \cr
GR02 \tab Boreia Peloponnesos \cr
GR03 \tab Anatolike Peloponnesos \cr
GR04 \tab Dytike Sterea Ellada \cr
GR05 \tab Epeiros  \cr
GR06 \tab Attike  \cr
GR07 \tab Anatolike Sterea Ellada \cr
GR08 \tab Thessalia  \cr
GR09 \tab Dytike Makedonia  \cr
GR10 \tab Kentrike Makedonia  \cr
GR11 \tab Anatolike Makedonia  \cr
GR12 \tab Thrake  \cr
GR13 \tab Krete  \cr
GR14 \tab Nesoi Aigaiou  \cr
}
}
\description{
\code{hydro_translate} translates various Hydroscope's names and
terms to English.
}
\note{
The dictionary used for the Greek to English translation is:

\tabular{ll}{
 \strong{Transliterated term} \tab \strong{English term} \cr
  agnosto               \tab unknown           \cr
  anemos                \tab wind              \cr
  dieuthynse            \tab direction         \cr
  parelthon             \tab past              \cr
  tachyteta             \tab speed             \cr
  mese                  \tab average           \cr
  brochoptose           \tab precipitation     \cr
  diarkeia              \tab duration          \cr
  exatmise              \tab evaporation       \cr
  exatmisodiapnoe       \tab evapotranspiration\cr
  thermokrasia          \tab temperature       \cr
  edaphous              \tab ground            \cr
  bathos                \tab depth             \cr
  elachiste             \tab min               \cr
  megiste               \tab max               \cr
  piese                 \tab pressure          \cr
  semeiake              \tab point             \cr
  chioni                \tab snow              \cr
  ypsometro             \tab elevation         \cr
  stathme               \tab level             \cr
  plemmyra              \tab flood             \cr
  paroche               \tab flow              \cr
  broche                \tab precipitation     \cr
  katastase             \tab condition         \cr
  ektimemene            \tab estimation        \cr
  athroistiko           \tab cumulative        \cr
  stereo                \tab sediment          \cr
  ygrasia               \tab humidity          \cr
  ygro                  \tab wet               \cr
  apolyte               \tab absolute          \cr
  schetike              \tab relative          \cr
  asbestio              \tab calcium           \cr
  wetu                  \tab precipitation     \cr
  chionobrochometro     \tab snow_rain_gauge   \cr
  xero                  \tab dry               \cr
  ydrometrese           \tab flow_gauge        \cr
  thalasses             \tab sea               \cr
  semeio_drosou         \tab dew_point         \cr
  oratoteta             \tab visibility        \cr
  steria                \tab land              \cr
  thalassa              \tab sea               \cr
  barometro             \tab barometer         \cr
  tase_ydratmon         \tab vapour_pressure   \cr
  psychrometro          \tab psychrometer      \cr
  isodynamo_ypsos       \tab water_equivalent  \cr
  agogimoteta           \tab conductance       \cr
  aktinobolia           \tab radiation         \cr
  anthraka              \tab carbon            \cr
  dioxeidio             \tab dioxide           \cr
  ypoloipo              \tab residual          \cr
  argilio               \tab aluminum          \cr
  argilos               \tab clay              \cr
  arseniko              \tab arsenic           \cr
  pyritiou              \tab silicon           \cr
  aera                  \tab air               \cr
  nephokalypse          \tab cloud_cover       \cr
  nephose               \tab clouds            \cr
  axiosemeiota          \tab remarkably        \cr
  nephe                 \tab clouds            \cr
  kairos                \tab weather           \cr
  diafora               \tab difference        \cr
  atmosfairiki          \tab atmospheric       \cr
  stathera              \tab constant          \cr
  parousa               \tab present           \cr
  parelthousa           \tab past              \cr
  kalymeno              \tab cover             \cr
  el.                   \tab min               \cr
  meg.                  \tab max               \cr
  skleroteta            \tab hardness          \cr
  eliophaneia           \tab sunshine          \cr
  eisroe_se_tamieuteres \tab inflow_reservoir   \cr
}
}
\examples{
\dontrun{

# get data from the Ministry of Environment and Energy
kyy_owners <- get_owners("kyy")
kyy_vars <- get_variables("kyy")
owners_names <- hydro_translate(kyy_owners$name, "owner")
vars <- hydro_translate(kyy_vars$descr, "variable")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{timeseries}
\alias{timeseries}
\title{timeseries}
\format{
A tibble with 10,804 rows and 9 variables:
\describe{
    \item{time_id}{The time series ID}
    \item{station_id}{The corresponding station's ID}
    \item{variable}{The time series variable type}
    \item{timestep}{The timestep of time series}
    \item{units}{The units of the time series}
    \item{start_date}{The starting date of time series values}
    \item{end_date}{The ending date of time series values}
    \item{subdomain}{The corresponding Hydroscope's database}
}
}
\usage{
timeseries
}
\description{
Time series' data from the Greek National Data Bank for
Hydrological and Meteorological Information.  This dataset is a comprehensive
look-up table of all of the available measurements for a given station in a
given Hydroscope's  database, with units of measurement and times  of
those measurements.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hydroscoper_package.R
\docType{package}
\name{hydroscoper}
\alias{hydroscoper}
\alias{_PACKAGE}
\alias{hydroscoper-package}
\title{hydroscoper: Interface to Hydroscope}
\description{
\code{hydroscoper} provides an R interface to the  Greek
National Data Bank for Hydrological and Meteorological Information
\code{http://www.hydroscope.gr}.

\code{hydroscoper} covers Hydroscope's data sources using the
\code{Enhydris API} and provides functions to:
\enumerate{
  \item {Transform the available tables and data sets into
       \href{https://tibble.tidyverse.org/}{tibbles}.}
  \item{Transliterate the Greek Unicode names to Latin.}
  \item{Translate various Greek terms to English.}
}
}
\section{Enhydris API}{


The Enhydris database is implemented in PostgreSQL. Details
can be found \href{https://enhydris.readthedocs.io}{here}
}

\section{Data Sources}{


The data are retrieved from the Hydroscope's databases:
\itemize{
  \item Ministry of Environment, Energy and Climate Change.
  \item Ministry of Rural Development and Food.
  \item  National Meteorological Service.
  \item  Greek Public Power Corporation.
}
}

\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/hydroscoper}
  \item \url{https://docs.ropensci.org/hydroscoper/}
  \item Report bugs at \url{https://github.com/ropensci/hydroscoper/issues}
}

}
\author{
\strong{Maintainer}: Konstantinos Vantas \email{kon.vantas@gmail.com} (\href{https://orcid.org/0000-0001-6387-8791}{ORCID})

Other contributors:
\itemize{
  \item Sharla Gelfand (Sharla Gelfand reviewed the package for rOpenSci, see https://github.com/ropensci/onboarding/issues/185) [contributor, reviewer]
  \item Tim Trice (Tim Trice reviewed the package for rOpenSci, see https://github.com/ropensci/onboarding/issues/185) [reviewer]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{stations}
\alias{stations}
\title{stations}
\format{
A tibble with 2,322 rows and 9 variables:
\describe{
    \item{station_id}{The station's ID from the domain's database}
    \item{name}{The station's name}
    \item{water_basin}{The station's Water Basin}
    \item{water_division}{The station's Water Division}
    \item{owner}{The station's owner}
    \item{longitude}{The station's longitude in decimal degrees, ETRS89}
    \item{latitude}{The station's latitude in decimal degrees, ETRS89}
    \item{altitude}{The station's altitude, meters above sea level}
    \item{subdomain}{The corresponding Hydroscope's database}
}
}
\usage{
stations
}
\description{
Stations' data from the Greek National Data Bank for
Hydrological and Meteorological Information. This dataset is a comprehensive
look-up table with geographical and ownership information of the available
stations in all Hydroscope's databases.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hydro_coords.R
\name{hydro_coords}
\alias{hydro_coords}
\title{Convert coordinates from Hydroscope's points to a tibble}
\usage{
hydro_coords(x)
}
\arguments{
\item{x}{a string vector with the points retrieved from Hydroscope}
}
\value{
a tibble with the longitude and latitude values.
}
\description{
\code{hydro_coords} returns a tibble with the stations'
longitude and latitude using as input the variable \code{point}
from \code{get_stations} function.
}
\examples{
\dontrun{
# get stations from the Greek Ministry of Environment and Energy
kyy_stations <- get_stations("kyy")

# create a tibble with stations' coords
coords <- hydro_coords(kyy_stations$point)
}
}
\author{
Konstantinos Vantas, \email{kon.vantas@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{get_data}
\alias{get_data}
\title{Get time series values in a tibble}
\usage{
get_data(subdomain = c("kyy", "ypaat", "emy", "deh"), time_id)
}
\arguments{
\item{subdomain}{One of the subdomains of hydroscope.gr}

\item{time_id}{A time series ID}
}
\value{
If \code{subdomain} is one of:
\itemize{
\item{\code{kyy}, Ministry of Environment and Energy}
\item{\code{ypaat}, Ministry of Rural Development and Food}
\item{\code{deh}, Greek Public Power Corporation}
\item{\code{emy}, National Meteorological Service}
}
and \code{time_id} exists in that \code{subdomain}, returns a tibble with
the time series values. Otherwise returns an error message.

The dataframe columns are:
\describe{
    \item{date}{The time series Dates (POSIXct)}
    \item{value}{The time series values (numeric)}
    \item{comment}{Comments about the values (character)}
}
}
\description{
\code{get_data} returns a tibble from  a Hydroscope's
time-series text file.
}
\note{
Data are not available freely in the sub-domains:  \code{"deh"}
(Greek Public Power Corporation) and
\code{"emy"} (National Meteorological Service).
}
\examples{
\dontrun{
# get time series 912 from the Greek Ministry of Environment and Energy
time_series <- get_data("kyy", 912)
}

}
\references{
Stations' data are retrieved from the Hydroscope's
databases:
\itemize{
\item Ministry of Environment, Energy and Climate Change.
\item Ministry of Rural Development and Food.
}
}
\author{
Konstantinos Vantas, \email{kon.vantas@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_tables.R
\name{get_tables}
\alias{get_tables}
\alias{get_stations}
\alias{get_timeseries}
\alias{get_instruments}
\alias{get_water_basins}
\alias{get_water_divisions}
\alias{get_political_divisions}
\alias{get_variables}
\alias{get_units_of_measurement}
\alias{get_time_steps}
\alias{get_owners}
\alias{get_instruments_type}
\alias{get_station_type}
\alias{get_database}
\title{Get tibbles from Hydroscope}
\usage{
get_stations(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)

get_timeseries(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)

get_instruments(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)

get_water_basins(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)

get_water_divisions(
  subdomain = c("kyy", "ypaat", "emy", "deh"),
  translit = TRUE
)

get_political_divisions(
  subdomain = c("kyy", "ypaat", "emy", "deh"),
  translit = TRUE
)

get_variables(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)

get_units_of_measurement(
  subdomain = c("kyy", "ypaat", "emy", "deh"),
  translit = TRUE
)

get_time_steps(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)

get_owners(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)

get_instruments_type(
  subdomain = c("kyy", "ypaat", "emy", "deh"),
  translit = TRUE
)

get_station_type(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)

get_database(subdomain = c("kyy", "ypaat", "emy", "deh"), translit = TRUE)
}
\arguments{
\item{subdomain}{One of the subdomains of Hydroscope in the vector
\code{c("kyy", "ypaat", "emy", "deh")}.}

\item{translit}{Automatically transliterate Greek to Latin.}
}
\value{
If \code{subdomain} is one of:
\itemize{
\item{\code{kyy}, Ministry of Environment and Energy.}
\item{\code{ypaat}, Ministry of Rural Development and Food.}
\item{\code{deh}, Greek Public Power Corporation.}
\item{\code{emy}, National Meteorological Service.}
}
returns a tibble or a named list with tibbles from the corresponding
database. Otherwise returns an error message.
}
\description{
A family of functions that return a tibble from a specific
 database from Hydroscope using the Enhydris API. \code{get_database} returns
 a named list of tibbles using all the family's functions.
}
\note{
Objects' IDs are not unique among the different Hydroscope databases.
For example, time series' IDs from http://kyy.hydroscope.gr have same values
with time series' from http://ypaat.hydroscope.gr.

The coordinates of the stations are based on the European Terrestrial
Reference System 1989 (ETRS89).
}
\examples{

\dontrun{

# data will be downloaded from Ministry of Environment and Energy (kyy):
subdomain <- "kyy"

# stations
kyy_stations <- get_stations(subdomain)

# time series
kyy_ts <- get_timeseries(subdomain)

# instruments
kyy_inst <- get_instruments(subdomain)

# water basins
kyy_wbas <- get_water_basins(subdomain)

# water divisions
kyy_wdiv <- get_water_divisions(subdomain)

# political divisions
kyy_pol <- get_political_divisions(subdomain)

# variables
kyy_vars <- get_variables(subdomain)

# units of measurement
kyy_units <- get_units_of_measurement(subdomain)

# time steps
kyy_time_steps <- get_time_steps(subdomain)

# owners
kyy_owners <- get_owners(subdomain)

# instruments type
kyy_instr_type <- get_instruments_type(subdomain)

# stations' type
kyy_st_type <- get_station_type(subdomain)

# use all the get_ functions above to create a named list with tibbles
kyy_db <- get_database(subdomain)
}

}
\references{
The data are retrieved from the Hydroscope's site databases:
\itemize{
\item Ministry of Environment, Energy and Climate Change.
\item Ministry of Rural Development and Food.
\item National Meteorological Service.
\item Greek Public Power Corporation.
}
}
\author{
Konstantinos Vantas, \email{kon.vantas@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{hydroscoper_defunct}
\alias{hydroscoper_defunct}
\alias{get_coords}
\title{Defunct functions in hydroscoper}
\usage{
get_coords(...)
}
\arguments{
\item{...}{Defunct function's parameters}
}
\description{
These functions are no longer available in \pkg{hydroscoper}.
}
\details{
Defunct functions:

\itemize{
 \item \code{\link{get_coords}}: This function is defunct. Please use
 \code{\link{hydro_coords}} to convert Hydroscope's points raw format to a
 tidy data frame.
}
}
