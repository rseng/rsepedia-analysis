---
title: bikedata
tags:
    - public hire bicycle
    - open data
    - R
authors:
    - name: Mark Padgham
      affiliation: 1
      orcid: 0000-0003-2172-5265
    - name: Richard Ellison
      affiliation: 2
affiliations:
    - name: Department of Geoinformatics, University of Salzburg, Austria
      index: 1
    - name: Institute of Transport and Logistics Studies, The University of Sydney, Australia
      index: 2
bibliography: paper.bib
date: 28 Nov 2017
---

# Summary

The R package `bikedata` collates and facilitates access to arguably the world's
largest open ongoing dataset on human mobility. All other comparable sources of
data (such public transit data, or mobile phone data) are either not publicly
available, or have been released only at single distinct times for single
distinct purposes. Many public hire bicycle systems in the U.S.A., along with
Santander Cycles in London, U.K., issue ongoing releases of their usage data,
providing a unique source of data for analysing, visualising, and understanding
human movement and urban environments [@Austwick2013; @Borgnat2011;
@Padgham2012].  Such data provide an invaluable resource for urban planners,
geographers, social and health scientists and policy makers, data visualisation
specialists, and data-affine users of the systems themselves.  The `bikedata`
package aims to provide unified access to usage statistics from all public hire
bicycle systems which provide  data. These currently including Santander Cycles
in London, U.K., and from the U.S.A., citibike in New York City NY, Divvy in
Chicago IL, Capital Bikeshare in Washington DC, Hubway in Boston MA, Metro in
Los Angeles LA, and Indego in Philadelphia PA. Additional systems will be added
on an ongoing basis.  The package facilitates the three necessary steps of (1)
downloading data; (2) storing data in a readily accessible form (in this case in
a single SQLite3 database); (3) extracting aggregate statistics.  The two
primary aggregate statistics are matrices of numbers of trips between all pairs
of stations, and daily time series. Both forms of aggregation may be extracted
for specific dates, times, or demographic characteristics of cyclists.

# References
<!-- README.md is generated from README.Rmd. Please edit that file -->

[![R build
status](https://github.com/ropensci/bikedata/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/bikedata/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/bikedata/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/bikedata)
[![Project Status:
Active](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/bikedata)](https://cran.r-project.org/package=bikedata)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/grand-total/bikedata?color=orange)](https://cran.r-project.org/package=bikedata)
[![](http://badges.ropensci.org/116_status.svg)](https://github.com/ropensci/software-review/issues/116)
[![status](https://joss.theoj.org/papers/10.21105/joss.00471/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00471)

The `bikedata` package aims to enable ready importing of historical trip
data from all public bicycle hire systems which provide data, and will
be expanded on an ongoing basis as more systems publish open data.
Cities and names of associated public bicycle systems currently
included, along with numbers of bikes and of docking stations (from
[wikipedia](https://en.wikipedia.org/wiki/List_of_bicycle-sharing_systems#Cities)),
are

| City                           | Hire Bicycle System                                                   | Number of Bicycles | Number of Docking Stations |
|--------------------------------|-----------------------------------------------------------------------|--------------------|----------------------------|
| London, U.K.                   | [Santander Cycles](https://tfl.gov.uk/modes/cycling/santander-cycles) | 13,600             | 839                        |
| San Francisco Bay Area, U.S.A. | [Ford GoBike](https://www.fordgobike.com/)                            | 7,000              | 540                        |
| New York City NY, U.S.A.       | [citibike](https://www.citibikenyc.com/)                              | 7,000              | 458                        |
| Chicago IL, U.S.A.             | [Divvy](https://www.divvybikes.com/)                                  | 5,837              | 576                        |
| Montreal, Canada               | [Bixi](https://bixi.com/)                                             | 5,220              | 452                        |
| Washingon DC, U.S.A.           | [Capital BikeShare](https://www.capitalbikeshare.com/)                | 4,457              | 406                        |
| Guadalajara, Mexico            | [mibici](https://www.mibici.net/)                                     | 2,116              | 242                        |
| Minneapolis/St Paul MN, U.S.A. | [Nice Ride](https://www.niceridemn.com/)                              | 1,833              | 171                        |
| Boston MA, U.S.A.              | [Hubway](https://www.bluebikes.com/)                                  | 1,461              | 158                        |
| Philadelphia PA, U.S.A.        | [Indego](https://www.rideindego.com)                                  | 1,000              | 105                        |
| Los Angeles CA, U.S.A.         | [Metro](https://bikeshare.metro.net/)                                 | 1,000              | 65                         |

These data include the places and times at which all trips start and
end. Some systems provide additional demographic data including years of
birth and genders of cyclists. The list of cities may be obtained with
the `bike_cities()` functions, and details of which include demographic
data with `bike_demographic_data()`.

The following provides a brief overview of package functionality. For
more detail, see the
[vignette](https://docs.ropensci.org/bikedata/articles/bikedata.html).

------------------------------------------------------------------------

## 1 Installation

Currently a development version only which can be installed with the
following command,

``` r
devtools::install_github("ropensci/bikedata")
```

and then loaded the usual way

``` r
library (bikedata)
```

## 2 Usage

Data may downloaded for a particular city and stored in an `SQLite3`
database with the simple command,

``` r
store_bikedata (city = 'nyc', bikedb = 'bikedb', dates = 201601:201603)
# [1] 2019513
```

where the `bikedb` parameter provides the name for the database, and the
optional argument `dates` can be used to specify a particular range of
dates (Jan-March 2016 in this example). The `store_bikedata` function
returns the total number of trips added to the specified database. The
primary objects returned by the `bikedata` packages are ‘trip matrices’
which contain aggregate numbers of trips between each pair of stations.
These are extracted from the database with:

``` r
tm <- bike_tripmat (bikedb = 'bikedb')
dim (tm); format (sum (tm), big.mark = ',')
```

    #> [1] 518 518
    #> [1] "2,019,513"

During the specified time period there were just over 2 million trips
between 518 bicycle docking stations. Note that the associated databases
can be very large, particularly in the absence of `dates` restrictions,
and extracting these data can take quite some time.

Data can also be aggregated as daily time series with

``` r
bike_daily_trips (bikedb = 'bikedb')
```

    #> # A tibble: 87 x 2
    #>    date       numtrips
    #>    <chr>         <dbl>
    #>  1 2016-01-01    11172
    #>  2 2016-01-02    14794
    #>  3 2016-01-03    15775
    #>  4 2016-01-04    19879
    #>  5 2016-01-05    18326
    #>  6 2016-01-06    24922
    #>  7 2016-01-07    28215
    #>  8 2016-01-08    29131
    #>  9 2016-01-08    21140
    #> 10 2016-01-10    14481
    #> # … with 77 more rows

A summary of all data contained in a given database can be produced as

``` r
bike_summary_stats (bikedb = 'bikedb')
#>    num_trips num_stations          first_trip       last_trip latest_files
#> ny  2019513          518 2016-01-01 00:00    2016-03-31 23:59        FALSE
```

The final field, `latest_files`, indicates whether the files in the
database are up to date with the latest published files.

### 2.1 Filtering trips by dates, times, and weekdays

Trip matrices can be constructed for trips filtered by dates, days of
the week, times of day, or any combination of these. The temporal extent
of a `bikedata` database is given in the above `bike_summary_stats()`
function, or can be directly viewed with

``` r
bike_datelimits (bikedb = 'bikedb')
```

    #>              first               last 
    #> "2016-01-01 00:00" "2016-03-31 23:59"

Additional temporal arguments which may be passed to the `bike_tripmat`
function include `start_date`, `end_date`, `start_time`, `end_time`, and
`weekday`. Dates and times may be specified in almost any format, but
larger units must always precede smaller units (so years before months
before days; hours before minutes before seconds). The following
examples illustrate the variety of acceptable formats for these
arguments.

``` r
tm <- bike_tripmat ('bikedb', start_date = "20160102")
tm <- bike_tripmat ('bikedb', start_date = 20160102, end_date = "16/02/28")
tm <- bike_tripmat ('bikedb', start_time = 0, end_time = 1) # 00:00 - 01:00
tm <- bike_tripmat ('bikedb', start_date = 20160101, end_date = "16,02,28",
                 start_time = 6, end_time = 24) # 06:00 - 23:59
tm <- bike_tripmat ('bikedb', weekday = 1) # 1 = Sunday
tm <- bike_tripmat ('bikedb', weekday = c('m', 'Th'))
tm <- bike_tripmat ('bikedb', weekday = 2:6,
                    start_time = "6:30", end_time = "10:15:25")
```

### 2.2 Filtering trips by demographic characteristics

Trip matrices can also be filtered by demographic characteristics
through specifying the three additional arguments of `member`, `gender`,
and `birth_year`. `member = 0` is equivalent to `member = FALSE`, and
`1` equivalent to `TRUE`. `gender` is specified numerically such that
values of `2`, `1`, and `0` respectively translate to female, male, and
unspecified. The following lines demonstrate this functionality

``` r
sum (bike_tripmat ('bikedb', member = 0))
sum (bike_tripmat ('bikedb', gender = 'female'))
sum (bike_tripmat ('bikedb', weekday = 'sat', birth_year = 1980:1990,
                   gender = 'unspecified'))
```

### 3. Citation

``` r
citation ("bikedata")
#> 
#> To cite bikedata in publications use:
#> 
#>   Mark Padgham, Richard Ellison (2017). bikedata Journal of Open Source Software, 2(20). URL
#>   https://doi.org/10.21105/joss.00471
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Article{,
#>     title = {bikedata},
#>     author = {Mark Padgham and Richard Ellison},
#>     journal = {The Journal of Open Source Software},
#>     year = {2017},
#>     volume = {2},
#>     number = {20},
#>     month = {Dec},
#>     publisher = {The Open Journal},
#>     url = {https://doi.org/10.21105/joss.00471},
#>     doi = {10.21105/joss.00471},
#>   }
```

### 4. Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
0.2.5
==================
Minor changes:
- `store_bikedata()` now has additional parameter `latest_lo_stns` that should
  generally just be left at default value, but which can be set to `FALSE` for
  truly offline use.
- Update bundled version of sqlite3 from 3.28 to 3.30
- minor bug fixes

0.2.4
==================
Back on CRAN after being removed due to dependency (dodgr) having been removed

0.2.3
===================
Minor changes:
- Fix dl_bikedata for Philadelphia
- Improve robustness of tests

0.2.2
===================
Minor changes:
- add NEWS & README to CRAN description
- Minor bug fixes

0.2.1
===================
- New helper function `bike_cities` to directly list cities included in current
  package version

Minor changes:
- Bug fix for San Fran thanks to @tbdv (see issue #78)
- Bug fix for LA (see issue #87)

0.2.0
===================
- Major expansion to include new cities of San Francisco, Minneapolis/St Paul,
  Montreal Canada, and Guadalajara Mexico
- Most code restructured to greatly ease the process of adding new cities (see
  github wiki for how to do this).
- New co-author: Tom Buckley

Minor changes:
- Bug fix for LA data which previously caused error due to invisible mac OSX
  system files being bundled with the distributed data
- More accurate date processing for quarterly LA data

0.1.1
===================
- Important bug in `dodgr` package rectified previously bug-ridden
  `bike_distmat()` calculations (thanks Joris Klingen).
- Files for Washington DC Capital Bike Share system changed from quarterly to
  annual dumps
- One rogue .xlsx file from London now processed and read properly (among all
  other well-behaved .csv files).
- Update bundled sqlite3: 3.21 -> 3.22


0.1.0
===================
- New function `bike_distmat()` calculates distance matrices between all pairs
  of stations as routed through street networks for each city.
- Helper function `bike_match_matrices()` matches distance and trip matrices by
  start and end stations, so they can be directly compared in standard
  statistical routines.
- North American Bike Share Association (NABSA) systems (currently LA and
  Philly) now distinguish member versus non-member based on whether usage is
  30-day pass or "Walk-up".

minor changes
- `dl_bikedata()` function also aliased to `download_bikedata()`, so both do the
  same job.
- Repeated runs of `store_bikedata()` on pre-existing databases sometimes
  re-added old data. This has now been fixed so only new data are added with
  each repeated call.
- Dates for NABSA cities (LA and Philadelphia) are given in different formats,
  all of which are now appropriately handled.
- Internally bundled sqlite3 upgraded v3.19 -> v3.21


0.0.4
===================
- Database no longer automatically indexed; rather indexes must be actively
  generated with `index_bikedata_db()`. This makes multiple usages of
  `store_bikedata()` faster and easier.
- `store_bikedata()` fixed so it only unzips files not already in database (it
  used to unzip them all)
- Internal changes to improve consistency (mostly through using the DBI
  package).


0.0.3
===================
- Minor changes only
- More informative messages when data for specified dates not available

0.0.2
===================
- No change to package functionality
- Drop dplyr dependency after dplyr 0.7 upgrade

0.0.1 (31 May 2017)
===================
- Initial CRAN release
# Contributing to bikedata

## Opening issues

The easiest way to note any behavioural curiosities or to request any new
features is by opening a [github issue](https://github.com/ropensci/bikedata/issues).


## Development guidelines

If you'd like to contribute changes to `bikedata`, we use [the GitHub
flow](https://guides.github.com/introduction/flow/index.html) for proposing,
submitting, reviewing, and accepting changes. If you haven't done this before,
there's a nice overview of git [here](http://r-pkgs.had.co.nz/git.html), as well
as best practices for submitting pull requests
[here](http://r-pkgs.had.co.nz/git.html#pr-make).

The `bikedata` coding style diverges somewhat from [this commonly used R style
guide](http://adv-r.had.co.nz/Style.html), primarily in the following two ways,
both of which improve code readability: (1) All curly braces are vertically aligned:
```r
this <- function ()
{
    x <- 1
}
```
and **not**
```r
this <- function(){
    x <- 1
}
```
and (2) Also highlighted in that code is the additional whitespace which
permeates `bikedata` code. Words of text are separated by whitespace, and so
code words should be too:
```r
this <- function1 (function2 (x))
```
and **not**
```r
this <- function1(function2(x))
```
with the natural result that one ends up writing
```r
this <- function ()
```
with a space between `function` and `()`. That's it.


## Code of Conduct

We want to encourage a warm, welcoming, and safe environment for contributing to
this project. See the [code of
conduct](https://github.com/ropensci/bikedata/blob/master/CODE_OF_CONDUCT.md)
for more information.
# CRAN notes for bikedata_0.2.5 submission

This package was previously, "Archived again on 2020-02-12 as check issues were not corrected on re-submission. UBSAN reports integer overflow, valgrind reports use of uninitialized values." My recent resubmission still manifest the integer overflow problem which this submission now rectifies. The problem arose through me overseeing an inline conversion to <int> as part of a <long int> variable. I have confirmed with g++ UBSAN that the present submission fixes the issue. Please accept my apologies for any inconvenience.

This submission has also been tested on:

# Test environments

- CRAN win-builder: R-oldrelease, R-release, R-devel
* Ubuntu 16.04 (on `travis-ci`): R-release, R-devel, R-oldrelease
* OSX: R-release (on `travis-ci`)

---
title: "bikedata"
keywords: "bicycle hire systems, bike hire systems, bike hire, bicycle hire, database, bike data"
output:
  rmarkdown::html_vignette:
    self_contained: no

  md_document:
    variant: markdown_github
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

[![R build status](https://github.com/ropensci/bikedata/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/bikedata/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/bikedata/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/bikedata)
[![Project Status: Active](http://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/bikedata)](https://cran.r-project.org/package=bikedata) 
[![CRAN Downloads](https://cranlogs.r-pkg.org/badges/grand-total/bikedata?color=orange)](https://cran.r-project.org/package=bikedata)
[![](http://badges.ropensci.org/116_status.svg)](https://github.com/ropensci/software-review/issues/116)
[![status](https://joss.theoj.org/papers/10.21105/joss.00471/status.svg)](https://joss.theoj.org/papers/10.21105/joss.00471)


The `bikedata` package aims to enable ready importing of historical trip data
from all public bicycle hire systems which provide data, and will be expanded on
an ongoing basis as more systems publish open data. Cities and names of
associated public bicycle systems currently included, along with numbers of
bikes and of docking stations (from 
[wikipedia](https://en.wikipedia.org/wiki/List_of_bicycle-sharing_systems#Cities)),
are

City | Hire Bicycle System | Number of Bicycles | Number of Docking Stations
--- | --- | --- | ---
London, U.K. | [Santander Cycles](https://tfl.gov.uk/modes/cycling/santander-cycles) | 13,600 | 839
San Francisco Bay Area, U.S.A. | [Ford GoBike](https://www.fordgobike.com/)  | 7,000 | 540 
New York City NY, U.S.A. | [citibike](https://www.citibikenyc.com/) | 7,000 | 458
Chicago IL, U.S.A. | [Divvy](https://www.divvybikes.com/) | 5,837 | 576
Montreal, Canada | [Bixi](https://bixi.com/) | 5,220 | 452
Washingon DC, U.S.A. | [Capital BikeShare](https://www.capitalbikeshare.com/) | 4,457 | 406
Guadalajara, Mexico | [mibici](https://www.mibici.net/) | 2,116 | 242
Minneapolis/St Paul MN, U.S.A. | [Nice Ride](https://www.niceridemn.com/) | 1,833 | 171
Boston MA, U.S.A. | [Hubway](https://www.bluebikes.com/) | 1,461 | 158
Philadelphia PA, U.S.A. | [Indego](https://www.rideindego.com) | 1,000 | 105
Los Angeles CA, U.S.A. | [Metro](https://bikeshare.metro.net/) | 1,000 | 65

These data include the places and times at which all trips start and end. Some
systems provide additional demographic data including years of birth and genders
of cyclists. The list of cities may be obtained with the `bike_cities()`
functions, and details of which include demographic data with
`bike_demographic_data()`.

The following provides a brief overview of package functionality. For more
detail, see the
[vignette](https://docs.ropensci.org/bikedata/articles/bikedata.html).

------


## 1 Installation

Currently a development version only which can be installed with the following
command,
```{r, eval=FALSE}
devtools::install_github("ropensci/bikedata")
```
```{r usage2, echo=FALSE, message=FALSE}
devtools::load_all (".")
#devtools::load_all (".", recompile=TRUE)
#devtools::document (".")
#goodpractice::gp ("bikedata")
#devtools::check (".")
#testthat::test_dir ("./tests/")
#Rcpp::compileAttributes()
```
and then loaded the usual way
```{r, eval = FALSE}
library (bikedata)
```


```{r echo=FALSE, message=FALSE, warning=FALSE, error=FALSE}
options(width = 120)
```


## 2 Usage

Data may downloaded for a particular city and stored in an `SQLite3` database
with the simple command,
```{r, echo = FALSE, eval = FALSE}
dl_bikedata (city = 'ny', data_dir = '/data/data/bikes/nyc-temp/',
             dates = 201601:201603)
store_bikedata (bikedb = 'bikedb', data_dir = '/data/data/bikes/nyc-temp/')
```
```{r, eval = FALSE}
store_bikedata (city = 'nyc', bikedb = 'bikedb', dates = 201601:201603)
# [1] 2019513
```
where the `bikedb` parameter provides the name for the database, and the
optional argument `dates` can be used to specify a particular range of dates
(Jan-March 2016 in this example).  The `store_bikedata` function returns the
total number of trips added to the specified database. The primary objects
returned by the `bikedata` packages are 'trip matrices' which contain aggregate
numbers of trips between each pair of stations. These are extracted from the
database with:
```{r, eval = FALSE}
tm <- bike_tripmat (bikedb = 'bikedb')
dim (tm); format (sum (tm), big.mark = ',')
```
```{r, echo = FALSE}
c (518, 518)
"2,019,513"
```
During the specified time period there were just over 2 million trips between
518 bicycle docking stations.  Note that the associated databases can be very
large, particularly in the absence of `dates` restrictions, and extracting these
data can take quite some time.

Data can also be aggregated as daily time series with
```{r, eval = FALSE}
bike_daily_trips (bikedb = 'bikedb')
```
```{r, echo = FALSE}
n <- 87
dates <- c ('2016-01-01', '2016-01-02', '2016-01-03', '2016-01-04',
            '2016-01-05', '2016-01-06', '2016-01-07', '2016-01-08',
            '2016-01-08', '2016-01-10', rep (NA, n - 10))
nt <- c (11172, 14794, 15775, 19879, 18326, 24922, 28215, 29131, 21140, 14481,
         rep (NA, n - 10))
tibble::tibble (date = dates, numtrips = nt)
```

A summary of all data contained in a given database can be produced as
```{r, eval = FALSE}
bike_summary_stats (bikedb = 'bikedb')
#>    num_trips num_stations          first_trip       last_trip latest_files
#> ny  2019513          518 2016-01-01 00:00    2016-03-31 23:59        FALSE
```
The final field, `latest_files`, indicates whether the files in the database are
up to date with the latest published files.

### 2.1 Filtering trips by dates, times, and weekdays

Trip matrices can be constructed for trips filtered by dates, days of the week,
times of day, or any combination of these.  The temporal extent of a `bikedata`
database is given in the above `bike_summary_stats()` function, or can be
directly viewed with
```{r, eval = FALSE}
bike_datelimits (bikedb = 'bikedb')
```
```{r, echo = FALSE}
res <- c ("2016-01-01 00:00", "2016-03-31 23:59")
names (res) <- c ("first", "last")
res
```
Additional temporal arguments which may be passed to the `bike_tripmat`
function include `start_date`, `end_date`, `start_time`, `end_time`, and
`weekday`. Dates and times may be specified in almost any format, but larger
units must always precede smaller units (so years before months before days;
hours before minutes before seconds). The following examples illustrate the
variety of acceptable formats for these arguments.
```{r, eval = FALSE}
tm <- bike_tripmat ('bikedb', start_date = "20160102")
tm <- bike_tripmat ('bikedb', start_date = 20160102, end_date = "16/02/28")
tm <- bike_tripmat ('bikedb', start_time = 0, end_time = 1) # 00:00 - 01:00
tm <- bike_tripmat ('bikedb', start_date = 20160101, end_date = "16,02,28",
                 start_time = 6, end_time = 24) # 06:00 - 23:59
tm <- bike_tripmat ('bikedb', weekday = 1) # 1 = Sunday
tm <- bike_tripmat ('bikedb', weekday = c('m', 'Th'))
tm <- bike_tripmat ('bikedb', weekday = 2:6,
                    start_time = "6:30", end_time = "10:15:25")
```

### 2.2 Filtering trips by demographic characteristics

Trip matrices can also be filtered by demographic characteristics through
specifying the three additional arguments of `member`, `gender`, and
`birth_year`. `member = 0` is equivalent to `member = FALSE`, and `1` equivalent
to `TRUE`. `gender` is specified numerically such that values of `2`, `1`, and
`0` respectively translate to female, male, and unspecified. The following lines
demonstrate this functionality
```{r, eval = FALSE}
sum (bike_tripmat ('bikedb', member = 0))
sum (bike_tripmat ('bikedb', gender = 'female'))
sum (bike_tripmat ('bikedb', weekday = 'sat', birth_year = 1980:1990,
                   gender = 'unspecified'))
```

### 3. Citation

```{r}
citation ("bikedata")
```

### 4. Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By contributing to this
project you agree to abide by its terms.

[![ropensci\_footer](https://ropensci.org//public_images/github_footer.png)](https://ropensci.org/)
# Script to assemble the `bike_test_data.rda` file

```{r}
devtools::load_all (".", export_all = FALSE)
library (magrittr)
data_dir <- tempdir ()
nrows <- 200 # number of rows to read from each file

names (bike_test_data)
head (bike_test_data$la)

# ----- DC -----
dl_bikedata (city = "dc", data_dir = data_dir, dates = 201701)
f <- list.files (tempdir ())
f <- file.path (data_dir, f [grep ("capitalbikeshare", f)])
fi <- unzip (f, list = TRUE)$Name
unzip (f, files = fi [1], exdir = data_dir, junkpaths = TRUE)
dc <- read.csv (file.path (data_dir, fi [1]), header = TRUE, nrows = nrows)

# ----- LO -----
dl_bikedata (city = "lo", data_dir = data_dir, dates = 201601)
f <- list.files (tempdir ())
f <- file.path (data_dir, f [grep ("JourneyDataExtract", f)])
lo <- read.csv (f [1], header = TRUE, nrows = nrows)

# ----- BO -----
# These 3 different time periods have 3 different formats
dl_bikedata (city = "bo", data_dir = data_dir, dates = 2012)
dl_bikedata (city = "bo", data_dir = data_dir, dates = 201701)
dl_bikedata (city = "bo", data_dir = data_dir, dates = 201801)
f <- list.files (tempdir ())
f <- file.path (data_dir, f [grep ("hubway", f)])
# grepping all at once doesn't put them in this order:
f <- c (f [grep ("2012", f)], f [grep ("2017", f)], f [grep ("2018", f)])
bo12 <- read.csv (f [1], header = TRUE, nrows = nrows)
fi <- unzip (f [2], list = TRUE)$Name
unzip (f [2], files = fi, exdir = data_dir, junkpaths = TRUE)
bo17 <- read.csv (file.path (data_dir, fi), header = TRUE, nrows = nrows)
fi <- unzip (f [3], list = TRUE)$Name
unzip (f [3], files = fi, exdir = data_dir, junkpaths = TRUE)
bo18 <- read.csv (file.path (data_dir, fi), header = TRUE, nrows = nrows)
# stations also need to be downloaded
dl_files <- bikedata:::get_bike_files (city = 'bo')
dl_files <- dl_files [which (grepl ('Stations', dl_files))]
for (f in dl_files)
{
    furl <- gsub (" ", "%20", f)
    f <- gsub (" ", "", f)
    destfile <- file.path (data_dir, basename(f))
    resp <- httr::GET (furl, httr::write_disk (destfile, overwrite = TRUE))
}
f <- list.files (tempdir ())
f <- file.path (data_dir, f [grep ("hubway_stations", f, ignore.case = TRUE)])
bo_st1 <- read.csv (f [1], header = TRUE)
bo_st2 <- read.csv (f [2], header = TRUE)

# ----- NY -----
dl_bikedata (city = "ny", data_dir = data_dir, dates = 201612)
f <- list.files (tempdir ())
f <- file.path (data_dir, f [grep ("^201612-citibike", f)])
fi <- unzip (f, list = TRUE)$Name
unzip (f, files = fi [1], exdir = data_dir, junkpaths = TRUE)
ny <- read.csv (file.path (data_dir, fi [1]), header = TRUE, nrows = nrows)

# ----- CH -----
dl_bikedata (city = "ch", data_dir = data_dir, dates = 201612)
f <- list.files (tempdir ())
f <- file.path (data_dir, f [grep ("Divvy", f)])
fi <- unzip (f, list = TRUE)$Name
fitr <- fi [grep ("Trips_2016_Q4", fi)]
fist <- fi [grep ("Stations_2016_Q4", fi)]
unzip (f, files = c (fitr, fist), exdir = data_dir, junkpaths = TRUE)
ch_tr <- read.csv (file.path (data_dir, fitr), header = TRUE, nrows = nrows)
ch_st <- read.csv (file.path (data_dir, fist), header = TRUE)

# ----- LA -----
dl_bikedata (city = "la", data_dir = data_dir, dates = 201701)
f <- list.files (tempdir ())
f <- file.path (data_dir, f [grep ("la_metro", f)])
fi <- unzip (f, list = TRUE)$Name
unzip (f, files = fi, exdir = data_dir, junkpaths = TRUE)
la <- read.csv (file.path (data_dir, fi), header = TRUE, nrows = nrows)

# ----- MN -----
# data have to be pre-downloaded
mn_dir <- "/data/data/bikes/mn"
f <- list.files (mn_dir, full.names = TRUE) [3] # random file for 2012
fi <- unzip (f, list = TRUE)$Name
fitr <- fi [grepl ("trip", fi)]
fist <- fi [grepl ("station", fi)]
unzip (f, files = c (fitr, fist), exdir = data_dir, junkpaths = TRUE)
mn_tr <- read.csv (file.path (data_dir, basename (fitr)),
                   header = TRUE, nrows = nrows)
mn_st <- read.csv (file.path (data_dir, basename (fist)),
                   header = TRUE)

bike_test_data <- list (dc = dc,
                        lo = lo,
                        bo12 = bo12,
                        bo17 = bo17,
                        bo18 = bo18,
                        bo_st1 = bo_st1,
                        bo_st2 = bo_st2,
                        ny = ny,
                        ch_tr = ch_tr,
                        ch_st = ch_st,
                        la = la,
                        mn_tr = mn_tr,
                        mn_st = mn_st)
save (bike_test_data, file = "./data/bike_test_data.rda", compress = "xz")
```
# The sysdata.rda object

## DC station locations

The function `R/stations.R/bike_get_dc_stations` has code to extract and process
DC stations. The data can be obtained from 
http://opendata.dc.gov/datasets/capital-bike-share-locations/, using
Download->Spreadsheet. The code is reproduced here
```{r}
stations_dc <- read.csv ("Capital_Bike_Share_Locations.csv")
names (stations_dc) <- tolower (names (stations_dc))
name <- noquote (gsub ("'", "", stations_dc$address)) #nolint
name <- trimws (name, which = 'right') # trim terminal white space
stations_dc <- data.frame (id = stations_dc$terminal_number,
                           name = name,
                           lon = stations_dc$longitude,
                           lat = stations_dc$latitude,
                           stringsAsFactors = FALSE)
```

## Bike Header Field Names

The fields stored in the `bikedata` database are:

| number | field                   |
| ----   | ----------------------- |
| 1      | duration                |
| 2      | start_time              |
| 3      | end_time                |
| 4      | start_station_id        |
| 5      | start_station_name      |
| 6      | start_station_latitude  |
| 7      | start_station_longitude |
| 8      | end_station_id          |
| 9      | end_station_name        |
| 10     | end_station_latitude    |
| 11     | end_station_longitude   |
| 12     | bike_id                 |
| 13     | user_type               |
| 14     | birth_year              |
| 15     | gender                  |

Each file has at least some of these fields, but different systems naturally use
different nomenclatures. The `header_names` structure maps different system
names for these fields onto the above names. All names are converted to lower
case and all white space and underscores removed, so entries here should be all
lower case with no white space.

old DC files had "Duration (ms)", but no longer do.
LA has "passholder_type", which can be "Flex Pass" = annual, or "Monthly Pass"
PH has "passholder_type", which can be "IndegoFlex" or "Indego30"

Note that an extra city column is needed because LA has "start_station" and
"end_station" for the ID columns, while MN has these for the station name
columns.

```{r}
fields <- c ("duration", "starttime", "endtime", "startstationid",
            "startstationname", "startstationlatitude",
            "startstationlongitude", "endstationid", "endstationname",
            "endstationlatitude", "endstationlongitude", "bikeid",
            "usertype", "birthyear", "gender")

duration <- c ("duration", "tripduration", "totalduration", "durationsec",
               "durationseconds", "totalduration(ms)")

starttime <- c ("starttime", "startdate", "iniciodelviaje")
endtime <- c ("endtime", "enddate", "stoptime", "findelviaje")

startstationid <- c ("startstationid", "startstationnumber", "fromstationid",
                     "startterminal", "startstation", "startstationcode",
                     "origenid")
startstationname <- c ("startstationname", "fromstationname", "startstation")
startstationlatitude <- c ("startstationlatitude", "startlat")
startstationlongitude <- c ("startstationlongitude", "startlon")

endstationid <- c ("endstationid", "endstationnumber", "tostationid",
                   "endstation", "endterminal", "endstationcode",
                   "destinoid")
endstationname <- c ("endstationname", "tostationname", "endstation")
endstationlatitude <- c ("endstationlatitude", "endlat")
endstationlongitude <- c ("endstationlongitude", "endlon")

bikeid <- c ("bikeid", "bikenumber", "bike#")
usertype <- c ("usertype", "membertype", "type", "subscribertype",
               "subscriptiontype", "accounttype", "passholdertype",
               "ismember", "usuarioid")
birthyear <- c ("birthyear", "birthday", "memberbirthyear",
                "edad","anodenacimento")
gender <- c ("gender", "membergender", "genero")

field_names <- data.frame (matrix (nrow = 0, ncol = 2))
for (f in fields)
{
    field_names <- rbind (field_names,
                          cbind (rep (f, length (get (f))), get (f)))
}
names (field_names) <- c ("field", "variation")
field_names$index <- field_names$field
levels (field_names$index) <- seq (unique (field_names$index))
field_names$index <- as.numeric (field_names$index)

field_names$city <- "all"
field_names$city [field_names$field == "startstationid" &
                  field_names$variation == "startstation"] <- "la"
field_names$city [field_names$field == "endstationname" &
                  field_names$variation == "endstation"] <- "mn"
field_names$city [field_names$field == "startstationname" &
                  field_names$variation == "startstation"] <- "mn"
field_names$city [field_names$field == "endstationid" &
                  field_names$variation == "endstation"] <- "la"
```

And this then saves the correponding `data.frame` to the package data:
```{r}
data_dir <- file.path (here::here (), "R")
f <- file.path (data_dir, "sysdata.rda")
load ("./R/sysdata.rda")
stations_dc <- sysdata$stations_dc # comment out to refresh using above code
sysdata <- list (stations_dc = stations_dc, field_names = field_names)
save (sysdata, file = f, compress = "xz")
```
---
title: "bikedata"
author: 
  - "Mark Padgham"
date: "`r Sys.Date()`"
output: 
    html_document:
        toc: true
        toc_float: true
        theme: flatly
vignette: >
  %\VignetteIndexEntry{bikedata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## 1. Introduction

`bikedata` is an R package for downloading and aggregating data from 
public bicycle hire, or bike share, systems. Although there are very many
public bicycle hire systems in the world ([see this wikipedia
list](https://en.wikipedia.org/wiki/List_of_bicycle-sharing_systems)),
relatively few openly publish data on system usage. The `bikedata` package aims
to enable ready importing of data from all systems which provide it, and will
be expanded on an ongoing basis as more systems publish open data. Cities and
names of associated public bicycle hire systems currently included in the
`bikedata` package, along with numbers of bikes and of docking stations, are:

City | Hire Bicycle System | Number of Bicycles | Number of Docking Stations
--- | --- | --- | ---
London, U.K. | [Santander Cycles](https://tfl.gov.uk/modes/cycling/santander-cycles) | 13,600 | 839
San Francisco Bay Area, U.S.A. | [Ford GoBike](https://www.fordgobike.com/)  | 7,000 | 540 
New York City NY, U.S.A. | [citibike](https://www.citibikenyc.com/) | 7,000 | 458
Chicago IL, U.S.A. | [Divvy](https://www.divvybikes.com/) | 5,837 | 576
Montreal, Canada | [Bixi](https://bixi.com/) | 5,220 | 452
Washingon DC, U.S.A. | [Capital BikeShare](https://www.capitalbikeshare.com/) | 4,457 | 406
Guadalajara, Mexico | [mibici](https://www.mibici.net/) | 2,116 | 242
Minneapolis/St Paul MN, U.S.A. | [Nice Ride](https://www.niceridemn.com/) | 1,833 | 171
Boston MA, U.S.A. | [Hubway](https://www.bluebikes.com/) | 1,461 | 158
Philadelphia PA, U.S.A. | [Indego](https://www.rideindego.com) | 1,000 | 105
Los Angeles CA, U.S.A. | [Metro](https://bikeshare.metro.net/) | 1,000 | 65

All of these systems record and disseminate individual trip data, minimally
including the times and places at which every trip starts and ends. Some provide
additional anonymised individual data, typically including whether or not a user
is registered with the system and if so, additional data including age, gender,
and residential postal code.  The list of cities may be obtained with the
`bike_cities()` functions, and details of which include demographic data with
`bike_demographic_data()`.

Cities with extensively developed systems and cultures of public hire bicycles,
yet which do not provide (publicly available) data include:

City | Number of Bicycles | Number of Docking Stations
--- | --- | ---
Hangzhou, China | 78,000 | 2,965
Paris, France | 14,500 | 1,229
Barcelona, Spain | 6,000 | 424

The current version of the `bikedata` R package can be installed with the
following command:
```{r install1, eval = FALSE}
install.packages ('bikedata')
```
Or the development version with
```{r install2, eval = FALSE}
devtools::install_github ("mpadge/bikedata")
```
Once installed, it can be loaded in the usual way:
```{r, eval = TRUE}
library (bikedata)
```
```{r, eval = FALSE, echo = FALSE, message = FALSE}
devtools::load_all ()
```


## 2. Main Functions

The `bikedata` function `dl_bikedata()` downloads individual trip data from
any or all or the above listed systems, and the function `store_bikedata()`
stores them in an `SQLite3` database.  For example, the following line will
download all data from the Metro system of Los Angeles CA, U.S.A., and store
them in a database named 'bikedb',
```{r store-la-data, eval = FALSE}
bikedb <- file.path (tempdir (), "bikedb.sqlite") # or whatever
dl_bikedata (city = 'la', dates = 2016, quiet = TRUE)
store_bikedata (data_dir = tempdir (), bikedb = bikedb, quiet = TRUE)
```
```{r store-la-data-dummy, echo = FALSE}
98138
```
The `store_bikedata()` function returns the number of trips added to the
database.  Both the downloaded data and the `SQLite3` database are stored by
default in the temporary directory of the current `R` session. Any data
downloaded into `tempdir()` will of course be deleted on termination of the
active **R** session; use of other directories (as described below) will create
enduring data which must be managed by the user.

Successive calls to `store_bikedata()` will append additional data to the same
database. For example, the following line will append all data from Chicago's
Divvy bike system from the year 2017 to the database created with the first
call above. 
```{r, eval = FALSE}
dl_bikedata (city = 'divvy', dates = 2016, quiet = TRUE)
store_bikedata (bikedb = bikedb, data_dir = tempdir (), quiet = TRUE)
```
```{r, echo = FALSE}
3595383
```
The function again returns the number of trips *added* to the database, which
is now less than the total number of trips stored of:
```{r, eval = FALSE}
bike_db_totals (bikedb = bikedb)
```
```{r, echo = FALSE}
3693521
```

Prior to accessing any data from the `SQLite3` database, it is recommended to
create database indexes using the function `index_bikedata_db()`:
```{r, eval = FALSE}
index_bikedata_db (bikedb = bikedb)
```
This will speed up subsequent extraction of aggregated data.

Having stored individual trip data in a database, the primary function of the
`bikedata` package is `bike_tripmat()`, which extracts aggregate numbers of
trips between all pairs of stations. The minimal arguments to this function are
the name of the database, and the name of a city for databases holding data from
multiple cities.
```{r, eval = FALSE}
tm <- bike_tripmat (bikedb = bikedb, city = 'la')
class (tm); dim (tm); sum (tm)
```
```{r}
```{r, echo = FALSE}
x <- "matrix"; y <- c (64, 64); z <- 98138
x; y; z
```
In 2016, the Los Angeles Metro system had 64 docking stations, and there were a
total of 98,138 individual trips during the year. Trip matrices can also be
extracted in long form using
```{r, eval = FALSE}
bike_tripmat (bikedb = bikedb, city = 'la', long = TRUE)
```
```{r, echo = FALSE}
n <- 4096
ss <- c (rep ('la3005', 10), rep (NA, n - 10))
es <- c ('la3005', 'la3006', 'la3007', 'la3008', 'la3010', 'la3011',
         'la3014', 'la3016', 'la3018', 'la3019', rep (NA, n - 10))
nt <- c (252, 93, 23, 153, 5, 63, 40, 10, 31, 36, rep (NA, n - 10))

tm <- tibble::tibble (start_station_id = ss,
                      end_station_id = es, numtrips = nt)
tm
```
Details of the docking stations associated with these trip matrices can be
obtained with
```{r, eval = FALSE}
bike_stations (bikedb = bikedb)
```
```{r, echo = FALSE}
n <- 660
id <- seq (n)
city <- rep ('la', n)
ids <- c ('la3005', 'la3006', 'la3007', 'la3008', 'la3009',
          'la3010', 'la3011', 'la3013', 'la3014', 'la3016', rep (NA, n - 10))
x <- c (34.04855, 34.04554, 34.05048, 34.04661, 33.98738,
        34.03705, 34.04113, 34.05661, 34.05290, 34.04373, rep (NA, n - 10))
y <- c (-118.2590, -118.2567, -118.2546, -118.2627, -118.4728,
        -118.2549, -118.2680, -118.2372, -118.2416, -118.2601, rep (NA, n - 10))
stations1 <- tibble::tibble (id = id, city = city, stn_id = ids,
                             name = rep ('', n), longitude = x, latitude = y)
stations1
```
Stations can also be extracted for particular cities:
```{r, eval = FALSE}
st <- bike_stations (bikedb = bikedb, city = 'ch')
```
For consistency and to avoid potential confusion of function names, most
functions in the `bikedata` package begin with the prefix `bike_` (except for
`store_bikedata()` and `dl_bikedata()`).

Databases generated by the `bikedata` package will generally be very large
(commonly at least several GB), and many functions may take considerable time to
execute.  It is nevertheless possible to explore package functionality quickly
through using the additional helper function, `bike_write_test_data()`. This
function uses the `bike_dat` data set provided with the package, which contains
details of 200 representative trips for each of the cities listed above. The
function writes these data to disk as `.zip` files which can then be read by the
`store_bikedata()` function.
```{r, eval = FALSE}
bike_write_test_data ()
store_bikedata (bikedb = 'testdb')
bike_summary_stats (bikedb = 'testdb')
```
The `.zip` files generated by `bike_write_test_data()` are created by default
in the `tempdir()` of the current `R` session, and so will be deleted on
session termination. Specifying any alternative `bike_dir` will create enduring
copies of those files in that location which ought to be deleted when finished.

The remainder of this vignette provides further detail on these three distinct
functional aspects of downloading, storage, and extraction of data.

## 3. Downloading Data

Data may be downloaded with the `dl_bikedata()` function. In it's simplest form,
this function requires specification only of a city for which data are to be
downloaded, although the directory is usually specified as well:
```{r, eval = FALSE}
dl_bikedata (city = 'chicago', data_dir = '/data/bikedata/')
```

## 3.1 Downloading data for specific date ranges

Both `store_bikedata()` and `dl_bikedata()` accept an additional argument
(`dates`) specifying ranges of dates for which data should be downloaded and
stored.  The format of this argument is quite flexible so that,
```{r, eval = FALSE}
dl_bikedata (city = 'dc', dates = 16)
```
will download data from Washington DC's Capital Bikeshare system for all 12
months of the year 2016, while,
```{r, eval = FALSE}
dl_bikedata (city = 'ny', dates = 201604:201608)
```
will download New York City data from April to August (inclusively) for that
year. (Note that the default `data_dir` is the `tempdir()` of the current `R`
session, with downloaded files being deleted upon session termination.) Dates
can also be entered as character strings, with the following calls producing
results equivalent to the preceding call, 
```{r, eval = FALSE}
dl_bikedata (city = 'ny', dates = '2016/04:2016/08')
dl_bikedata (city = 'new york', dates = '201604:201608')
dl_bikedata (city = 'n.y.c.', dates = '2016-04:2016-08')
dl_bikedata (city = 'new', dates = '2016 Apr-Aug')
```
The only strict requirement for the format of `dates` is that years must be
specified before months, and that some kind of separator must be used between
the two except when formatted as single six-digit numbers or character strings
(YYYYMM).  The arguments `city = 'new'` and `city = 'CI'` in the final call are
sufficient to uniquely identify New York City's citibike system.

If files have been previously downloaded to a nominated directory, then calling
the `dl_bikedata()` function will only download those data files that do not
already exist. This function may thus be used to periodically refresh the
contents of a nominated directory as new data files become available.

Some systems disseminate data on quarterly (Washington DC and Los Angeles) or
bi-annual (Chicago) bases. The `dates` argument in these cases is translated to
the appropriate quarterly or bi-annual files. These are then downloaded as
single files, and thus the following call
```{r, eval = FALSE}
dl_bikedata (city = 'dc', dates = '2016.03-2016.05')
```
will actually download data for the entire first and second quarters of 2016.
Even though the database constructed with `store_bikedata()` will then contain
data beyond the specified date ranges, it is nevertheless possible to obtain a
trip matrix corresponding to specific dates and/or times, as described below.

The `dates` argument can also be passed to `store_bikedata`. This is useful in
cases where data are to be loaded only from a restricted set of files in the
given data directory.

## 3.2 Refreshing data sources

The `dl_bikedata()` function will only download data that do not already exist
in the nominated directory, and so can be use to periodically refresh data. If,
for example, the following function were previously run at the end of 2017:
```{r, eval = FALSE}
dl_bikedata (city = 'sf', data_dir = '/data/stored/here')
```
then running again in, say, April 2018, would download three additional files
corresponding to the first three months of 2018. These data can then be added to
a previously-constructed database with the usual call
```{r, eval = FALSE}
store_bikedata (city = 'sf', data_dir = '/data/stored/here', bikedb = bikedb)
```
If previous data have been stored in a nominated database, yet deleted from
local storage, then any new data can be added by first getting the names of
previously stored files with
```{r, eval = FALSE}
bike_stored_files (bikedb = bikedb = city = 'sf')
```
And then calling `dl_bikedata`, with `dates` specified to add only those files
not previously stored.


## 4. Storing Data

As mentioned above, individual trip data are stored in a single `SQLite3`
database, created by default in the temporary directory of the current `R`
session. Specifying a path for the `bikedb` argument in the `store_bikedata()`
function will create a database that will remain in that location until
explicitly deleted. 

The nominated database is created if it does not already exist, otherwise
additional data are appended to the existing database. As described above, the
same `dates` argument can be passed to both `dl_bikedata()` and
`store_bikedata()` to download data within specified ranges of dates.

Both `dl_bikedata()` and `store_bikedata()` are primarily intended to be used to
download data for specified cities. It is possible to use the latter to store
all data for all cities simply by calling `store_bikedata (bikedb = bikedb)`,
however doing so will request confirmation that data from *all* cities really
ought to be downloaded and/or stored.  Intended general usage of the
`store_bikedata()` function is illustrated in the following line:
```{r, eval = FALSE}
dl_bikedata (bikedb = bikedb, city = 'ny', dates = '2014 aug - 2015 dec')
ntrips <- store_bikedata (bikedb = bikedb, city = 'ny',
                          data_dir = '/data/stored/here')
```
Note that passing with `city` parameter to `store_bikedata()` is not strictly
necessary, but will ensure that only data for the nominated city are loaded from
directories which may contain additional data from other cities.


## 5. Accessing Aggregate Data

### 5.1 Origin-Destination Matrices

As briefly described in the introduction, the primary function for extracting
aggregate data from the `SQLite3` database established with `store_bikedata()`
is `bike_tripmat()`. With the single mandatory argument naming the database,
this function returns a matrix of numbers of trips between all pairs of
stations.  Trip matrices can be returned either in square form (the default),
with both rows and columns named after the bicycle docking stations and matrix
entries tallying numbers of rides between each pair of stations, or in long form
by requesting `bike_tripmat (..., long = TRUE)`. The latter case will return a
[`tibble`](https://cran.r-project.org/package=tibble) with the three columns of
`station_station_id`, `end_station_id`, and `number_trips`, as demonstrated
above.

The data for the individual stations associated with the trip matrix can be
extracted with `bike_stations()`, which returns a `tibble` containing the 6
columns of city, station code, station name, longitude, and latitude. Station
codes are specified by the operators of each system, and pre-pended with a
2-character city identifier (so, for example, the first of the stations shown
above is `la3005`). The `bike_stations()` function will generally return all
operational stations within a given system, which `bike_tripmat()` will return
only those stations in operation during the requested time period. The previous
call stored all data from Chicago's Divvybikes system for the year 2016 only, so
the trip matrix has less entries than the full stations table, which includes
stations added since then.
```{r, eval = FALSE}
dim (bike_tripmat (bikedb = bikedb, city = 'ch'))
```
```{r, echo = FALSE}
c (581, 581)
```
```{r, eval = FALSE}
dim (bike_stations (bikedb = bikedb, city = 'ch'))
```
```{r, echo = FALSE}
c (596, 6)
```


#### 5.1.1. Temporal filtering of trip matrices

Trip matrices can also be extracted for particular dates, times, and days of
the week, through specifying one or more of the optional arguments:

1. `start_date`
2. `end_date`
3. `start_time`
4. `end_time`
5. `weekday`

Arguments may in all cases be specified in a range of possible formats as long
as they are unambiguous, and as long as 'larger' units precede 'smaller' units
(so years before months before days, and hours before minutes before seconds).
Acceptable formats may be illustrated through specifying a list of arguments to
be passed to `bike_tripmat()`. This is done here through passing two lists to
`bike_tripmat()` via `do.call()`, enabling the second list (`args1`) to be
subsequently modified.
```{r, eval = FALSE}
args0 <- list (bikedb = bikedb, city = 'ny', args)
args1 <- list (start_date = 16, end_time = 12, weekday = 1)
tm <- do.call (bike_tripmat, c (args0, args1))
```
In `args1`, a two-digit `start_date` (or `end_date`) is interpreted to represent
a year, while a one- or two-digit `_time` is interpreted to represent an hour.
A value of `end_time = 24` is interpreted as `end_time = '23:59:59'`, while a
value of `_time = 0` is interpreted as `00:00:00`.  The following further
illustrate the variety of acceptable formats,
```{r, eval = FALSE}
args1 <- list (start_date = '2016 May', end_time = '12:39', weekday = 2:6)
args1 <- list (end_date = 20160720, end_time = 123915, weekday = c ('mo', 'we'))
args1 <- list (end_date = '2016-07-20', end_time = '12:39:15', weekday = 2:6)
```
Both `_date` and `_time` arguments may be specified in either `character` or
`numeric` forms; in the former case with arbitrary (or no) separators.
Regardless of format, larger units must precede smaller units as explained
above.

Weekdays may specified as characters, which must simply be unambiguous and (in
admission of currently inadequate internationalisation) correspond to standard
English names. Minimal character specifications are thus `'so', 'm', 'tu', 'w',
'th', 'f', 'sa'`. The value of `weekday = 1` denotes Sunday, so `weekdays =
2:6` denote the traditional working days, Monday to Friday, while weekends may
be denoted with `weekdays = c ('sa', 'so')` or `weekdays = c (1, 7)`.


#### 5.1.2. Demographic filtering of trip matrices

As described at the outset, the bicycle hire systems of several cities provide
additional demographic information including whether or not cyclists are
registered with the system, and if so, additional information including birth
years and genders. Note that the provision of such information is voluntary, and
that no providers can or do guarantee the accuracy of their data.

Those systems which provide demographic information are listed with the
function `bike_demographic_data()`, which also lists the nominal kinds of
demographic data provided by the different systems.
```{r}
bike_demographic_data ()
```
Data can then be filtered by demographic parameters with additional optional
arguments to `bike_tripmat()` of,

1. `registered` (`TRUE/FALSE`, `'yes'/'no'`, 0/1)
2. `birth_year` (as one or more four-digit numbers or character strings)
3. `gender` ('m/f/.', 'male/female/other')

Users are not required to specify genders, and any values of `gender` other than
character strings beginning with either `f` or `m` (case-insensitive) will be
interpreted to request non-specified or alternative values of gender.  Note
further than many systems offer a range of potential birth years starting from a
default value of 1900, and there are consequently a significant number of
cyclists who declare this as their birth year.

It is of course possible to combine all of these optional parameters in a
single query. For example,
```{r, eval = FALSE}
tm <- bike_tripmat (bikedb = bikedb, city = 'ny', start_date = 2016,
        start_time = 9, end_time = 24, weekday = 2:6, gender = 'xx', 
        birth_year = 1900:1950)
```
The value of `gender = 'xx'` will be interpreted to request data from all
members with nominal alternative genders.  As demographic data are only given
for registered users, the `registered` parameter is redundant in this query.

#### 5.1.3. Standardising trip matrices by durations of operation

Most bicycle hire systems have progressively expanded over time through ongoing
addition of new docking stations. Total numbers of counts within a trip matrix
will thus be generally less for more recently installed stations, and more for
older stations. The `bike_tripmat()` function has an option, `standardise =
FALSE`. Setting `standardise = TRUE` allows trip matrices to be standardised for
durations of station operation, so that numbers of trips between any pair of
stations reflect what they would be if all stations had been in operation for
the same duration.

Standardisation implements a linear scaling of total numbers of trips to and
from each station according to total durations of operation, with counts in
the final trip matrix scaled to have the same total number of trips as the
original matrix.  This standardisation has two immediate consequences:

1. Numbers of trips between any pair of stations will not necessarily be integer
   values, but are rounded for the sake of sanity to three digits, corresponding
   to the maximal likely precision attainable for daily differences in operating
   durations;
2. Trip numbers will generally not equal actual observed numbers. Counts for the
   longest operating durations will be lower than actually recorded, while
   counts for more recent stations will be greater than observed values.

The `standardise` option nevertheless enables travel patterns between different
(groups of) stations to be statistically compared in a way that is free of the
potentially confounding influence of differing durations of operation.



### 5.2. Station Data

Data on docking stations may be accessed with the function `bike_stations()`
as demonstrated above:
```{r, eval = FALSE}
bike_stations (bikedb = bikedb)
```
```{r, echo = FALSE}
stations1
```
This function returns a [`tibble`](https://cran.r-project.org/package=tibble)
detailing the names and locations of all bicycle stations present in the
database. Station data for specific cities may be extracted through specifying
an additional `city` argument.
```{r, eval = FALSE}
bike_stations (bikedb = bikedb, city = 'ch')
```
```{r, echo = FALSE}
n <- 596
id <- 66 + seq (n)
city <- rep ('ch', n)
ids <- c ('ch456', 'ch101', 'ch109', 'ch21', 'ch80',
          'ch346', 'ch341', 'ch444', 'ch511', 'ch376', rep (NA, n - 10))
nms <- c ('2112 W Peterson Ave', '63rd St Beach', '900 W Harrison St',
            'Aberdeen St & Jackson Blvd', 'Aberdeen St & Monroe St', 
            'Ada St & Washington Blvd', 'Adler Planetarium', 
            'Albany Ave & 26th St', 'Albany Ave & Bloomingdale Ave', 
            'Artesian Ave & Hubbard St', rep (NA, n - 10))
x <- c (-87.68359, -87.57612, -87.65002, -87.65479, -87.65393,
        -87.66121, -87.60727, -87.70201, -87.70513, -87.68822, rep (NA, n - 10))
y <- c (41.99118, 41.78102, 41.87468, 41.87773, 41.88046,
        41.88283, 41.86610, 41.84448, 41.91403, 41.88949, rep (NA, n - 10))
stations2 <- tibble::tibble (id = id, city = city, stn_id = ids,
                             name = nms, longitude = x, latitude = y)
stations2
```



### 5.3. Summary Statistics

`bikedata` provides a number of helper functions for extracting summary
statistics from the `SQLite3` database. The function `bike_summary_stats
(bikedb)` generates an overview table. (This function may take some time to
execute on large databases.)

```{r, eval = FALSE}
bike_summary_stats (bikedb)
```
```{r, echo = FALSE}
ntr <- c (3693521, 3595383, 98138)
nst <- c (662, 596, 66)
startd <- c ('2016-01-01 00:07:00', '2016-01-01 00:07:00', 
             '2016-07-07 04:17:00')
endd <- rep ('2016-12-31 23:57:52', 3)
fls <- c (NA, TRUE, TRUE)
tbl <- tibble::tibble (city = c ("total", "ch", "la"),
                       num_trips = ntr, num_stations = nst,
                       first_trip = as.factor (startd), 
                       last_trip = as.factor (endd), latest_files = fls)
tbl
```
Additional helper functions provide individual components from this summary
data, and will generally do so notably faster for large databases than the above
function. The primary individual function is `bike_db_totals()`, which can be
used to extract total numbers of either trips (the default) or stations (by
specifying `trips = FALSE`) from the entire database or from specific cities. 
```{r, eval = FALSE}
bike_db_totals (bikedb = bikedb)
```
```{r, echo = FALSE}
3693521
```
```{r, eval = FALSE}
bike_db_totals (bikedb = bikedb, city = "ch")
```
```{r, echo = FALSE}
3595383
```
```{r, eval = FALSE}
bike_db_totals (bikedb = bikedb, city = "la")
```
```{r, echo = FALSE}
93138
```
```{r, eval = FALSE}
bike_db_totals (bikedb = bikedb, trips = FALSE)
```
```{r, echo = FALSE}
660
```
```{r, eval = FALSE}
bike_db_totals (bikedb = bikedb, trips = FALSE, city = "ch")
```
```{r, echo = FALSE}
596
```
```{r, eval = FALSE}
bike_db_totals (bikedb = bikedb, trips = FALSE, city = "la")
```
```{r, echo = FALSE}
64
```
The other primary components of `bike_summary_stats()` are the dates of first
and last trips for the entire database and for individual cities. These dates
can be obtained directly with the function `bike_datelimits()`:
```{r, eval = FALSE}
bike_datelimits (bikedb = bikedb)
```
```{r, echo = FALSE}
c ('first' = "2016-01-01 00:07:00", 'last' = "2016-12-31 23:57:52")
```
```{r, eval = FALSE}
bike_datelimits (bikedb = bikedb, city = 'ch')
```
```{r}
c ('first' = "2016-01-01 00:07:00", 'last' = "2016-12-31 23:57:52")
```
A helper function is also provided to determine whether the files stored in the
database represent the latest available files.
```{r, eval = FALSE}
bike_latest_files (bikedb = bikedb)
```
```{r}
c ('la' = TRUE, 'ch' = FALSE)
```


### 5.4. Time Series of Daily Trips

The `bike_tripmat()` function provides a *spatial* aggregation of data. An
equivalent *temporal* aggregation is provided by the function
`bike_daily_trips()`, which aggregates trips for individual days.

```{r, eval = FALSE}
bike_daily_trips (bikedb = bikedb, city = 'ch')
```
```{r, echo = FALSE}
n <- 366
dates <- c ('2016-01-01', '2016-01-02', '2016-01-03', '2016-01-04',
            '2016-01-05', '2016-01-06', '2016-01-07', '2016-01-08',
            '2016-01-09', '2016-01-10', rep (NA, n - 10))
nt <- c (935, 1421, 1399, 3833, 4189, 4608, 5028, 3425, 1733, 993,
         rep (NA, n - 10))
tibble::tibble (date = dates, numtrips = nt)
```
Daily trip counts can also be standardised to account for differences in numbers
of stations within a system as for trip matrix standardisation described above.
Such standardisation is helpful because daily numbers of trips will generally
increase with increasing numbers of stations. Standardisation returns a time
series of daily trips reflecting what they would be if all system stations had
been in operation throughout the entire time.
```{r, eval = FALSE}
bike_daily_trips (bikedb = bikedb, city = 'ch', standardise = TRUE)
```
```{r, echo = FALSE}
nt <- c (2468.925, 2481.939, 2200.766, 5509.787, 5884.207, 6298.229, 6630.111,
         4476.455, 2265.021, 1297.845, rep (NA, n - 10))
tibble::tibble (date = dates, numtrips = nt)
```
This [`tibble`](https://cran.r-project.org/package=tibble) reveals two points of
immediate note:

1. Trip numbers are no longer integer values, but are rounded to three decimal
   places to reflect the highest plausible numerical accuracy; and
2. Standardised trip numbers are considerably higher for the initial values,
   because of expansion of the Chicago Divvy system throughout the year 2016.


## 6. Direct database access


Although the `bikedata` package aims to circumvent any need to access the
database directly, through providing ready extraction of trip data for most
analytical or visualisation needs, direct access may be achieved either using
the convenient `dplyr` functions, or the more powerful functionality provided
by the `RSQLite` package.

The following code illustrates access using the `dplyr` package:
```{r, eval = FALSE}
db <- dplyr::src_sqlite (bikedb, create=F)
dplyr::src_tbls (db)
```
```{r}
c ("datafiles", "stations", "trips")
```
```{r, eval = FALSE}
dplyr::collect (dplyr::tbl (db, 'datafiles'))
```
```{r, echo = FALSE}
city <- c ('la', 'la', 'la', 'ch', 'ch')
nms <- c ('la_metro_gbfs_trips_Q1_2017.zip', 'MetroBikeShare_2016_Q3_trips.zip',
          'Metro_trips_Q4_2016.zip', 'Divvy_Trips_2016_Q1Q2.zip',
          'Divvy_Trips_2016_Q3Q4.zip')
tibble::tibble (id = 0:4, city = city, name = nms)
```
```{r, eval = FALSE}
dplyr::collect (dplyr::tbl (db, 'stations'))
```
```{r, echo = FALSE}
stations1
```
```{r, eval = FALSE}
dplyr::collect (dplyr::tbl (db, 'trips'))
```
```{r, echo = FALSE}
n <- 3693511
id <- seq (n)
city <- rep ('la', n)
dur <- c (180, 1980, 300, 10860, 420, 780, 600, 600, 2880, 960,
          rep (NA, n - 10))
st <- c ('2016-01-01 00:15:00', '2016-01-01 00:24:00', '2016-01-01 00:28:00',
         '2016-01-01 00:38:00', '2016-01-01 00:38:00', '2016-01-01 00:39:00',
         '2016-01-01 00:43:00', '2016-01-01 00:56:00', '2016-01-01 00:57:00',
         '2016-01-01 01:54:00', rep (NA, n - 10))
et <- c ('2017-01-01 00:23:00', '2017-01-01 00:36:00', '2017-01-01 00:45:00',
         '2017-01-01 00:43:00', '2017-01-01 00:43:00', '2017-01-01 00:59:00',
         '2017-01-01 00:55:00', '2017-01-01 01:44:00', '2017-01-01 01:44:00',
         '2017-01-01 02:19:00', rep (NA, n - 10))
nachr <- rep ('a', n)
naint <- rep (1L, n)
tibble::tibble (id = id, city = city, trip_duration = dur,
                start_time = st, stop_time = et,
                start_station_id = nachr, end_station_id = nachr,
                bike_id = nachr, user_type = nachr, birth_year = naint,
                gender = naint)
```
The [`RSQLite`](https://cran.r-project.org/package=RSQLite) package enables more
complex queries to be constructed. The names of stations, for example, could be
extracted using the following code
```{r, eval = FALSE}
db <- RSQLite::dbConnect(RSQLite::SQLite(), bikedb, create = FALSE)
qry <- "SELECT stn_id, name FROM stations WHERE city = 'ch'"
stns <- RSQLite::dbGetQuery(db, qry)
RSQLite::dbDisconnect(db)
head (stns)
```
```{r, echo = FALSE}
data.frame (stn_id = c ('ch456', 'ch101', 'ch109', 'ch21', 'ch80', 'ch346'),
            name = c ('2112 W Peterson Ave', '63rd St Beach', 
                      '900 W Harrison St', 'Aberdeen St & Jackson Blvd', 
                      'Aberdeen St & Monroe St', 'Ada St & Washington Blvd'))
```
Many of the queries used in the `bikedata` package are constructed in this way
using the `RSQLite` interface.


## 7. Visualisation of bicycle trips

The `bikedata` package does not provide any functions enabling visualisation of
aggregate trip data, both because of the primary focus on enabling access and
aggregation in the simplest practicable way, and because of the myriad
different ways users of the package are likely to want to visualise the data.
This section therefore relies on other packages to illustrate some of the ways
in which trip matrices may be visualised.

### 7.1 Visualisation using R Base functions

The simplest spatial visualisation involves connecting the geographical
coordinates of stations with straight lines, with numbers of trips represented
by some characteristics of the lines connecting pairs of stations, such as
thickness or colours.  This can be achieved with the following code, which also
illustrates that it is generally more useful for visualisation purposes to
extract trip matrices in long rather than square form.
```{r plot-la-out, echo = TRUE, eval = FALSE}
stns <- bike_stations (bikedb = bikedb, city = 'la')
ntrips <- bike_tripmat (bikedb = bikedb, city = 'la', long = TRUE)
x1 <- stns$longitude [match (ntrips$start_station_id, stns$stn_id)]
y1 <- stns$latitude [match (ntrips$start_station_id, stns$stn_id)]
x2 <- stns$longitude [match (ntrips$end_station_id, stns$stn_id)]
y2 <- stns$latitude [match (ntrips$end_station_id, stns$stn_id)]
# Set plot area to central region of bike system
xlims <- c (-118.27, -118.23)
ylims <- c (34.02, 34.07)
plot (stns$longitude, stns$latitude, xlim = xlims, ylim = ylims)
cols <- rainbow (100)
nt <- ceiling (ntrips$numtrips * 100 / max (ntrips$numtrips))
for (i in seq (x1))
    lines (c (x1 [i], x2 [i]), c (y1 [i], y2 [i]), col = cols [nt [i]],
        lwd = ntrips$numtrips [i] * 10 / max (ntrips$numtrips))
```
![](la_map_simple.png)

### 7.2 A More Sophisticated Visualisation

The following code illustrates a more sophisticated approach to plotting such
data, using routines from the packages `osmdata`, `stplanr`, and `tmap`. Begin
by extracting the street network for Los Angeles using the `osmdata` package.
Current `stplanr` routines require spatial objects of class
[`sp`](https://cran.r-project.org/package=sp) rather than
[`sf`](https://cran.r-project.org/package=sf).
```{r get-la-highways, eval = FALSE}
library (magrittr)
xlims_la <- range (stns$longitude, na.rm = TRUE)
ylims_la <- range (stns$latitude, na.rm = TRUE)
# expand those limits slightly
ex <- 0.1
xlims_la <- xlims_la + c (-ex, ex) * diff (xlims_la)
ylims_la <- ylims_la + c (-ex, ex) * diff (ylims_la)
bbox <- c (xlims_la [1], ylims_la [1], xlims_la [2], ylims_la [2])
bbox <- c (xlims [1], xlims [2], ylims [1], ylims [2])
# Then the actual osmdata query to extract all OpenStreetMap highways
highways <- osmdata::opq (bbox = bbox) %>%
    osmdata::add_osm_feature (key = 'highway') %>%
    osmdata::osmdata_sp (quiet = FALSE)
```
For compatibility with current `stplanr` code, the `stns` table also needs to be
converted to a `SpatialPointsDataFrame` and re-projected.
```{r convert-stns-to-spdf, eval = FALSE}
stns_tbl <- bike_stations (bikedb = bikedb)
stns <- sp::SpatialPointsDataFrame (coords = stns_tbl[,c('longitude','latitude')],
                                    proj4string = sp::CRS("+init=epsg:4326"), 
                                    data = stns_tbl)
stns <- sp::spTransform (stns, highways$osm_lines@proj4string)
```
These data can then be used to create an `stplanr::SpatialLinesNetwork` which
can be used to trace the routes between bicycle stations along the street
network.  This first requires mapping the bicycle station locations to the
nearest nodes in the street network, and converting the start and end stations
of the `ntrips` table to corresponding rows in the street network data frame.
```{r map-stns-to-streetnet, eval = FALSE}
la_net <- stplanr::SpatialLinesNetwork (sl = highways$osm_lines)
# Find the closest node to each station
nodeid <- stplanr::find_network_nodes (la_net, stns$longitude, stns$latitude)
# Convert start and end station IDs in trips table to node IDs in `la_net`
startid <- nodeid [match (ntrips$start_station_id, stns$stn_id)]
endid <- nodeid [match (ntrips$end_station_id, stns$stn_id)]
ntrips$start_station_id <- startid
ntrips$end_station_id <- endid
```
The aggregate trips on each part of the network using the `sum_network_lines()`
function which is part of the current development version of `stplanr`.
```{r, eval = FALSE}
bike_usage <- stplanr::sum_network_links (la_net, data.frame (ntrips))
```
Then finally plot it with `tmap`, again trimming the plot using the previous
limits to exclude a very few isolated stations
```{r plot-tmat, eval = FALSE}
tmap::tm_shape (bike_usage, xlim = xlims, ylim = ylims, is.master=TRUE) + 
    tmap::tm_lines (col="numtrips", lwd="numtrips", title.col = "Number of trips",
                    breaks = c(0, 200, 400, 600, 800, 1000, Inf),
                    legend.lwd.show = FALSE, scale = 5) + 
    tmap::tm_layout (bg.color="gray95", legend.position = c ("right", "bottom"),
                     legend.bg.color = "white", legend.bg.alpha = 0.5)
#tmap::save_tmap (filename = "la_map.png")
```
![](la_map.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database-stats.R
\name{bike_db_totals}
\alias{bike_db_totals}
\title{Count number of entries in sqlite3 database tables}
\usage{
bike_db_totals(bikedb, trips = TRUE, city)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.}

\item{trips}{If true, numbers of trips are counted; otherwise numbers of
stations}

\item{city}{Optional city for which numbers of trips are to be counted}
}
\description{
Count number of entries in sqlite3 database tables
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
bikedb <- file.path (data_dir, 'testdb')
# latest_lo_stns is set to FALSE just to avoid download on CRAN; this should
# normally remain at default value of TRUE:
store_bikedata (data_dir = data_dir, bikedb = bikedb, latest_lo_stns = FALSE)
# create database indexes for quicker access:
index_bikedata_db (bikedb = bikedb)

bike_db_totals (bikedb = bikedb, trips = TRUE) # total trips
bike_db_totals (bikedb = bikedb, trips = TRUE, city = 'ch')
bike_db_totals (bikedb = bikedb, trips = TRUE, city = 'ny')
bike_db_totals (bikedb = bikedb, trips = FALSE) # total stations
bike_db_totals (bikedb = bikedb, trips = FALSE, city = 'ch')
bike_db_totals (bikedb = bikedb, trips = FALSE, city = 'ny')
# numbers of stations can also be extracted with
nrow (bike_stations (bikedb = bikedb))
nrow (bike_stations (bikedb = bikedb, city = 'ch'))

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files ('.', pattern = '.zip'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/store-bikedata.R
\name{bike_rm_db}
\alias{bike_rm_db}
\title{Remove SQLite3 database generated with 'store_bikedat()'}
\usage{
bike_rm_db(bikedb)
}
\arguments{
\item{bikedb}{The SQLite3 database containing the bikedata.}
}
\value{
TRUE if \code{bikedb} successfully removed; otherwise FALSE
}
\description{
If no directory is specified the \code{bikedb} argument passed to
\code{store_bikedata}, the database is created in \code{tempdir()}. This
function provides a convenient way to remove the database in such cases by
simply passing the name.
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
# or download some real data!
# dl_bikedata (city = "la", data_dir = data_dir)
bikedb <- file.path (data_dir, "testdb")
store_bikedata (data_dir = data_dir, bikedb = bikedb)

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files (data_dir, pattern = ".zip"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database-stats.R
\name{bike_datelimits}
\alias{bike_datelimits}
\title{Extract date-time limits from trip database}
\usage{
bike_datelimits(bikedb, city)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.
If no directory specified, it is presumed to be in \code{tempdir()}.}

\item{city}{If given, date limits are calculated only for trips in
that city.}
}
\value{
A vector of 2 elements giving the date-time of the first and last
trips
}
\description{
Extract date-time limits from trip database
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
# dl_bikedata (city = 'la', data_dir = data_dir) # or some real data!
# Remove one London file that triggers an API call which may fail tests:
file.remove (file.path (tempdir(),
             "01aJourneyDataExtract10Jan16-23Jan16.csv"))
bikedb <- file.path (data_dir, 'testdb')
store_bikedata (data_dir = data_dir, bikedb = bikedb)
# create database indexes for quicker access:
index_bikedata_db (bikedb = bikedb)

bike_datelimits ('testdb') # overall limits for all cities
bike_datelimits ('testdb', city = 'NYC')
bike_datelimits ('testdb', city = 'los angeles')

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files ('.', pattern = '.zip'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/store-bikedata.R
\name{index_bikedata_db}
\alias{index_bikedata_db}
\title{Add indexes to database created with store_bikedata}
\usage{
index_bikedata_db(bikedb)
}
\arguments{
\item{bikedb}{The SQLite3 database containing the bikedata.}
}
\description{
Add indexes to database created with store_bikedata
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
# or download some real data!
# dl_bikedata (city = "la", data_dir = data_dir)
bikedb <- file.path (data_dir, "testdb")
store_bikedata (data_dir = data_dir, bikedb = bikedb)
# create database indexes for quicker access:
index_bikedata_db (bikedb = bikedb)

trips <- bike_tripmat (bikedb = bikedb, city = "LA") # trip matrix
stations <- bike_stations (bikedb = bikedb) # station data

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files (data_dir, pattern = ".zip"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/store-bikedata.R
\name{store_bikedata}
\alias{store_bikedata}
\title{Store hire bicycle data in SQLite3 database}
\usage{
store_bikedata(
  bikedb,
  city,
  data_dir,
  dates = NULL,
  latest_lo_stns = TRUE,
  quiet = FALSE
)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database to
use. If it doesn't already exist, it will be created, otherwise data
will be appended to existing database.  If no directory specified,
it is presumed to be in \code{tempdir()}.}

\item{city}{One or more cities for which to download and store bike data, or
names of corresponding bike systems (see Details below).}

\item{data_dir}{A character vector giving the directory containing the
data files downloaded with \code{dl_bikedata} for one or more
cities. Only if this parameter is missing will data be downloaded.}

\item{dates}{If specified and no \code{data_dir} is given, data are
downloaded and stored only for these dates specified as vector of YYYYMM
values.}

\item{latest_lo_stns}{If \code{TRUE} (default), download latest version of
London stations; otherwise use potentially obsolete internal version. (This
parameter should not need to be changed, but can be set to \code{FALSE} to
avoid external calls; for example when not online.)}

\item{quiet}{If FALSE, progress is displayed on screen}
}
\value{
Number of trips added to database
}
\description{
Store previously downloaded data (via the \link{dl_bikedata} function) in a
database for subsequent extraction and analysis.
}
\note{
Data for different cities may all be stored in the same database, with
city identifiers automatically established from the names of downloaded data
files. This function can take quite a long time to execute, and may generate
an SQLite3 database file several gigabytes in size.
}
\section{Details}{

City names are not case sensitive, and must only be long enough to
unambiguously designate the desired city. Names of corresponding bike systems
can also be given.  Currently possible cities (with minimal designations in
parentheses) and names of bike hire systems are:
\tabular{lr}{
 Boston (bo)\tab Hubway\cr
 Chicago (ch)\tab Divvy Bikes\cr
 Washington, D.C. (dc)\tab Capital Bike Share\cr
 Los Angeles (la)\tab Metro Bike Share\cr
 London (lo)\tab Santander Cycles\cr
 Minnesota (mn)\tab NiceRide\cr
 New York City (ny)\tab Citibike\cr
 Philadelphia (ph)\tab Indego\cr
 San Francisco Bay Area (sf)\tab Ford GoBike\cr
}
}

\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
# or download some real data!
# dl_bikedata (city = "la", data_dir = data_dir)
bikedb <- file.path (data_dir, "testdb")
store_bikedata (data_dir = data_dir, bikedb = bikedb)
# create database indexes for quicker access:
index_bikedata_db (bikedb = bikedb)

trips <- bike_tripmat (bikedb = bikedb, city = "LA") # trip matrix
stations <- bike_stations (bikedb = bikedb) # station data

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files (data_dir, pattern = ".zip"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bikedata-package.R
\docType{package}
\name{bikedata}
\alias{bikedata}
\title{Download and aggregate data from public bicycle hire systems}
\description{
Download data from all public bicycle hire systems which provide open data,
currently including
\itemize{
\item Santander Cycles London, U.K.
\item citibike New York City NY, U.S.A.
\item Divvy Chicago IL, U.S.A.
\item Capital BikeShare Washingon DC, U.S.A.
\item Hubway Boston MA, U.S.A.
\item Metro Los Angeles CA, U.S.A.
}
}
\section{Download and store data}{

\itemize{
\item \code{dl_bikedata} Download data for particular cities and dates
\item \code{store_bikedata} Store data in \code{SQLite3} database
}
}

\section{Sample data for testing package}{

\itemize{
\item \code{bike_test_data} Description of test data included with package
\item \code{bike_write_test_data} Write test data to disk in form precisely
reflecting data provided by all systems
\item \code{bike_rm_test_data} Remove data written to disk with
\code{bike_write_test_data}
}
}

\section{Functions to aggregate trip data}{

\itemize{
\item \code{bike_daily_trips} Aggregate daily time series of total trips
\item \code{bike_stations} Extract table detailing locations and names of
bicycle docking stations
\item \code{bike_tripmat} Extract aggregate counts of trips between all pairs
of stations within a given city
}
}

\section{Summary Statistics}{

\itemize{
\item \code{bike_summary_stats} Overall quantitative summary of database
contents.  All of the following functions provide individual aspects of this
summary.
\item \code{bike_db_totals} Count total numbers of trips or stations, either
for entire database or a specified city.
\item \code{bike_datelimits} Return dates of first and last trips, either for
entire database or a specified city.
\item \code{bike_demographic_data} Simple table indicating which cities
include demographic parameters with their data
\item \code{bike_latest_files} Check whether files contained in database are
latest published versions
}
}

\author{
Mark Padgham
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stations.R
\name{bike_stations}
\alias{bike_stations}
\title{Extract station matrix from SQLite3 database}
\usage{
bike_stations(bikedb, city)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.
If no directory specified, it is presumed to be in \code{tempdir()}.}

\item{city}{Optional city (or vector of cities) for which stations are to be
extracted}
}
\value{
Matrix containing data for each station
}
\description{
Extract station matrix from SQLite3 database
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
# or download some real data!
# dl_bikedata (city = 'la', data_dir = data_dir)
bikedb <- file.path (data_dir, 'testdb')
store_bikedata (data_dir = data_dir, bikedb = bikedb)
# create database indexes for quicker access:
index_bikedata_db (bikedb = bikedb)

stations <- bike_stations (bikedb)
head (stations)

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files (data_dir, pattern = '.zip'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write-test-data.R
\name{bike_write_test_data}
\alias{bike_write_test_data}
\title{Writes test data bundled with package to zip files}
\usage{
bike_write_test_data(data_dir = tempdir())
}
\arguments{
\item{data_dir}{Directory in which data are to be extracted. Defaults to
\code{tempdir()}. If any other directory is specified, files ought to be
removed with \code{bike_rm_test_data()}.}
}
\description{
Writes very small test files to disk that can be used to test the package.
The entire package works by reading zip-compressed data files provided by the
various hire bicycle systems. This function generates some equivalent data
that can be read into an \code{SQLite} database by the
\code{store_bikedata()} function, so that all other package functionality can
then be tested from the resultant database. This function is also used in the
examples of all other functions.
}
\examples{
\dontrun{
bike_write_test_data ()
list.files (tempdir ())
bike_rm_test_data ()

bike_write_test_data (data_dir = '.')
list.files ()
bike_rm_test_data (data_dir = '.')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tripmat.R
\name{bike_tripmat}
\alias{bike_tripmat}
\title{Extract station-to-station trip matrix or data.frame from SQLite3 database}
\usage{
bike_tripmat(
  bikedb,
  city,
  start_date,
  end_date,
  start_time,
  end_time,
  weekday,
  member,
  birth_year,
  gender,
  standardise = FALSE,
  long = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.
If no directory specified, it is presumed to be in \code{tempdir()}.}

\item{city}{City for which tripmat is to be aggregated}

\item{start_date}{If given (as year, month, day) , extract only those records
from and including this date}

\item{end_date}{If given (as year, month, day), extract only those records to
and including this date}

\item{start_time}{If given, extract only those records starting from and
including this time of each day}

\item{end_time}{If given, extract only those records ending at and including
this time of each day}

\item{weekday}{If given, extract only those records including the nominated
weekdays. This can be a vector of numeric, starting with Sunday=1, or
unambiguous characters, so "sa" and "tu" for Saturday and Tuesday.}

\item{member}{If given, extract only trips by registered members
(\code{member = 1} or \code{TRUE}) or not (\code{member = 0} or
\code{FALSE}).}

\item{birth_year}{If given, extract only trips by registered members whose
declared birth years equal or lie within the specified value or values.}

\item{gender}{If given, extract only records for trips by registered
users declaring the specified genders (\code{f/m/.} or \code{2/1/0}).}

\item{standardise}{If TRUE, numbers of trips are standardised to the
operating durations of each stations, so trip numbers are increased for
stations that have only operated a short time, and vice versa.}

\item{long}{If FALSE, a square tripmat of (num-stations, num_stations) is
returned; if TRUE, a long-format matrix of (stn-from, stn-to, ntrips) is
returned.}

\item{quiet}{If FALSE, progress is displayed on screen}
}
\value{
If \code{long = FALSE}, a square matrix of numbers of trips between
each station, otherwise a long-form \pkg{tibble} with three columns of of
(\code{start_station_id, end_station_id, numtrips}).
}
\description{
Extract station-to-station trip matrix or data.frame from SQLite3 database
}
\note{
The \code{city} parameter should be given for databases containing data
from multiple cities, otherwise most of the resultant trip matrix is likely
to be empty.  Both dates and times may be given either in numeric or
character format, with arbitrary (or no) delimiters between fields. Single
numeric times are interpreted as hours, with 24 interpreted as day's end at
23:59:59.

If \code{standardise = TRUE}, the trip matrix will have the same number
of trips, but they will be re-distributed as described, with more recent
stations having more trips than older stations. Trip number are also
non-integer in this case, whereas they are always integer-valued for
\code{standardise = FALSE}.
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
# or download some real data!
# dl_bikedata (city = "la", data_dir = data_dir)
bikedb <- file.path (data_dir, "testdb")
store_bikedata (data_dir = data_dir, bikedb = bikedb)
# create database indexes for quicker access:
index_bikedata_db (bikedb = bikedb)


tm <- bike_tripmat (bikedb = bikedb, city = "ny") # full trip matrix
tm <- bike_tripmat (bikedb = bikedb, city = "ny",
                    start_date = 20161201, end_date = 20161201)
tm <- bike_tripmat (bikedb = bikedb, city = "ny", start_time = 1)
tm <- bike_tripmat (bikedb = bikedb, city = "ny", start_time = "01:00")
tm <- bike_tripmat (bikedb = bikedb, city = "ny", end_time = "01:00")
tm <- bike_tripmat (bikedb = bikedb, city = "ny",
                    start_date = 20161201, start_time = 1)
tm <- bike_tripmat (bikedb = bikedb, city = "ny", start_date = 20161201,
                    end_date = 20161201, start_time = 1, end_time = 2)
tm <- bike_tripmat (bikedb = bikedb, city = "ny", weekday = 5)
tm <- bike_tripmat (bikedb = bikedb, city = "ny",
                    weekday = c("f", "sa", "th"))
tm <- bike_tripmat (bikedb = bikedb, city = "ny",
                    weekday = c("f", "th", "sa"))
tm <- bike_tripmat (bikedb = bikedb, city = "ny", member = 1)
tm <- bike_tripmat (bikedb = bikedb, city = "ny", birth_year = 1976)
tm <- bike_tripmat (bikedb = bikedb, city = "ny", birth_year = 1976:1990)
tm <- bike_tripmat (bikedb = bikedb, city = "ny", gender = "f")
tm <- bike_tripmat (bikedb = bikedb, city = "ny",
                    gender = "m", birth_year = 1976:1990)

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files (data_dir, pattern = ".zip"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database-stats.R
\name{bike_summary_stats}
\alias{bike_summary_stats}
\title{Extract summary statistics of database}
\usage{
bike_summary_stats(bikedb)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.
If no directory specified, it is presumed to be in \code{tempdir()}.}
}
\value{
A \code{data.frame} containing numbers of trips and stations along
with times and dates of first and last trips for each city in database and a
final column indicating whether the files match the latest published
versions.
}
\description{
Extract summary statistics of database
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
# dl_bikedata (city = "la", data_dir = data_dir) # or some real data!
# Remove one London file that triggers an API call which may fail tests:
file.remove (file.path (tempdir(),
             "01aJourneyDataExtract10Jan16-23Jan16.csv"))
bikedb <- file.path (data_dir, "testdb")
store_bikedata (data_dir = data_dir, bikedb = bikedb)
# create database indexes for quicker access:
index_bikedata_db (bikedb = bikedb)

bike_summary_stats ("testdb")

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files (".", pattern = ".zip"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distmat.R
\name{bike_distmat}
\alias{bike_distmat}
\title{Extract station-to-station distance matrix}
\usage{
bike_distmat(bikedb, city, expand = 0.5, long = FALSE, quiet = TRUE)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.
If no directory specified, it is presumed to be in \code{tempdir()}.}

\item{city}{City for which tripmat is to be aggregated}

\item{expand}{Distances are calculated by routing through the OpenStreetMap
street network surrounding the bike stations, with the street network
expanded by this amount to ensure all stations can be connected.}

\item{long}{If FALSE, a square distance matrix of (num-stations,
num_stations) is returned; if TRUE, a long-format matrix of (stn-from,
stn-to, distance) is returned.}

\item{quiet}{If FALSE, progress is displayed on screen}
}
\value{
If \code{long = FALSE}, a square matrix of numbers of trips between
each station, otherwise a long-form \pkg{tibble} with three columns of of
(start_station_id, end_station_id, distance)
}
\description{
Extract station-to-station distance matrix
}
\note{
Distance matrices returned from \code{bike_distamat} use all stations
listed for a given system, while trip matrices extracted with
\link{bike_tripmat} will often have fewer stations because operational
station numbers commonly vary over time. The two matrices may be reconciled
with the \code{match_trips2dists} function, enabling then to be directly
compared.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/write-test-data.R
\name{bike_rm_test_data}
\alias{bike_rm_test_data}
\title{Removes test data written with 'bike_write_test_data()'}
\usage{
bike_rm_test_data(data_dir = tempdir())
}
\arguments{
\item{data_dir}{Directory in which data were extracted.}
}
\value{
Number of files successfully removed, which should equal six.
}
\description{
The function \code{bike_write_test_data()} writes several small
zip-compressed files to disk. The default location is \code{tempdir()}, in
which case these files will be automatically removed on termination of
current R session. If, however, any other value for \code{data_dir} is passed
to \code{bike_write_test_data()}, then the resultant files ought be deleted
by calling this function.
}
\examples{
\dontrun{
bike_write_test_data ()
list.files (tempdir ())
bike_rm_test_data ()

bike_write_test_data (data_dir = getwd ())
list.files ()
bike_rm_test_data (data_dir = getwd ())
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{lo_stns}
\alias{lo_stns}
\title{Docking stations for London, U.K.}
\format{
A \code{data.frame} of the four columns described above.
}
\usage{
lo_stns
}
\description{
A \code{data.frame} of station id values, names, and geographic coordinates
for 786 stations for London, U.K. These stations are generally (and by
default) downloaded automatically to ensure they are always up to date, but
such downloading can be disabled in the \code{store_bikedata()} function by
setting \code{latest_lo_stns = FALSE}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/dl-bikedata.R
\name{dl_bikedata}
\alias{dl_bikedata}
\alias{download_bikedata}
\title{Download hire bicycle data}
\usage{
dl_bikedata(city, data_dir = tempdir(), dates = NULL, quiet = FALSE)

download_bikedata(city, data_dir = tempdir(), dates = NULL, quiet = FALSE)
}
\arguments{
\item{city}{City for which to download bike data, or name of corresponding
bike system (see Details below).}

\item{data_dir}{Directory to which to download the files}

\item{dates}{Character vector of dates to download data with dates formated
as YYYYMM.}

\item{quiet}{If FALSE, progress is displayed on screen}
}
\description{
Download data for subsequent storage via \link{store_bikedata}.
}
\note{
Only files that don't already exist in \code{data_dir} will be
downloaded, and this function may thus be used to update a directory of files
by downloading more recent files. If a particular file request fails,
downloading will continue regardless. To ensure all files are downloaded,
this function may need to be run several times until a message appears
declaring that 'All data files already exist'
}
\section{Details}{

This function produces (generally) zip-compressed data in R's temporary
directory. City names are not case sensitive, and must only be long enough to
unambiguously designate the desired city. Names of corresponding bike systems
can also be given.  Currently possible cities (with minimal designations in
parentheses) and names of bike hire systems are:
\tabular{lr}{
 Boston (bo)\tab Hubway\cr
 Chicago (ch)\tab Divvy Bikes\cr
 Washington, D.C. (dc)\tab Capital Bike Share\cr
 Los Angeles (la)\tab Metro Bike Share\cr
 London (lo)\tab Santander Cycles\cr
 Minnesota (mn)\tab NiceRide\cr
 New York City (ny)\tab Citibike\cr
 Philadelphia (ph)\tab Indego\cr
 San Francisco Bay Area (sf)\tab Ford GoBike\cr
}

Ensure you have a fast internet connection and at least 100 Mb space
}

\examples{
\dontrun{
dl_bikedata (city = 'New York City USA', dates = 201601:201613)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database-stats.R
\name{bike_demographic_data}
\alias{bike_demographic_data}
\title{Static summary of which systems provide demographic data}
\usage{
bike_demographic_data()
}
\value{
A \code{data.frame} detailing the kinds of demographic data provided
by the different systems
}
\description{
Static summary of which systems provide demographic data
}
\examples{
bike_demographic_data ()
# Examples of filtering data by demographic parameters:
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
bikedb <- file.path (data_dir, "testdb")
store_bikedata (data_dir = data_dir, bikedb = bikedb)
# create database indexes for quicker access:
index_bikedata_db (bikedb = bikedb)

sum (bike_tripmat (bikedb = bikedb, city = "bo")) # 200 trips
sum (bike_tripmat (bikedb = bikedb, city = "bo", birth_year = 1990)) # 9
sum (bike_tripmat (bikedb = bikedb, city = "bo", gender = "f")) # 22
sum (bike_tripmat (bikedb = bikedb, city = "bo", gender = 2)) # 22
sum (bike_tripmat (bikedb = bikedb, city = "bo", gender = 1)) # = m; 68
sum (bike_tripmat (bikedb = bikedb, city = "bo", gender = 0)) # = n; 9
# Sum of gender-filtered trips is less than total because \code{gender = 0}
# extracts all registered users with unspecified genders, while without
# gender filtering extracts all trips for registered and non-registered
# users.

# The following generates an error because Washinton DC's DivvyBike system
# does not provide demographic data
sum (bike_tripmat (bikedb = bikedb, city = "dc", birth_year = 1990))
bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distmat.R
\name{bike_match_matrices}
\alias{bike_match_matrices}
\title{Match rows and columns of distance and trip matrices}
\usage{
bike_match_matrices(mat1, mat2)
}
\arguments{
\item{mat1}{A wide- or long-form trip or distance matrix returned from
\code{\link{bike_tripmat}} or \code{\link{bike_distmat}}.}

\item{mat2}{The corresponding distance or trip matrix.}
}
\value{
A list of the same matrices with matching start and end stations, and
in the same order passed to the routine (that is, \code{mat1} then
\code{mat2}). Each kind of matrix will be identified and named accordingly as
either "trip" or "dist". Matrices are returned in same format (long or wide)
as submitted.
}
\description{
Match rows and columns of distance and trip matrices
}
\note{
Distance matrices returned from \code{bike_distamat} use all stations
listed for a given system, while trip matrices extracted with
\link{bike_tripmat} will often have fewer stations because operational
station numbers commonly vary over time. This function reconciles the two
matrices through matching all row and column names (or just station IDs for
long-form matrices), enabling then to be directly compared.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{bike_test_data}
\alias{bike_test_data}
\title{Test data for all 6 cities}
\format{
A list of one data frame for each of the five cities of (bo, dc, la,
lo, ny), plus two more for chicago stations and trips (ch_st, ch_tr). Each of
these (except "ch_st") contains 200 representative trips.
}
\usage{
bike_test_data
}
\description{
A data set containing for each of the six cities a \code{data.frame} object
of 200 trips.
}
\note{
These data are only used to convert to \code{.zip}-compressed files
using \code{bike_write_test_data()}. These \code{.zip} files can be
subsequently read into an SQLite3 database using \code{store_bikedata}.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{bike_cities}
\alias{bike_cities}
\title{List of cities currently included in bikedata}
\usage{
bike_cities()
}
\value{
A \code{data.frame} of cities, abbreviations, and names of bike
systems currently able to be accessed.
}
\description{
List of cities currently included in bikedata
}
\examples{
bike_cities ()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database-stats.R
\name{bike_daily_trips}
\alias{bike_daily_trips}
\title{Extract daily trip counts for all stations}
\usage{
bike_daily_trips(
  bikedb,
  city,
  station,
  member,
  birth_year,
  gender,
  standardise = FALSE
)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.
If no directory specified, it is presumed to be in \code{tempdir()}.}

\item{city}{City for which trips are to be counted - mandatory if database
contains data for more than one city}

\item{station}{Optional argument specifying bike station for which trips are
to be counted}

\item{member}{If given, extract only trips by registered members
(\code{member = 1} or \code{TRUE}) or not (\code{member = 0} or
\code{FALSE}).}

\item{birth_year}{If given, extract only trips by registered members whose
declared birth years equal or lie within the specified value or values.}

\item{gender}{If given, extract only records for trips by registered
users declaring the specified genders (\code{f/m/.} or \code{2/1/0}).}

\item{standardise}{If TRUE, daily trip counts are standardised to the
relative numbers of bike stations in operation for each day, so daily trip
counts are increased during (generally early) periods with relatively fewer
stations, and decreased during (generally later) periods with more stations.}
}
\value{
A \code{data.frame} containing daily dates and total numbers of trips
}
\description{
Extract daily trip counts for all stations
}
\examples{
\dontrun{
bike_write_test_data () # by default in tempdir ()
# dl_bikedata (city = "la", data_dir = data_dir) # or some real data!
store_bikedata (data_dir = tempdir (), bikedb = "testdb")
# create database indexes for quicker access:
index_bikedata_db (bikedb = "testdb")

bike_daily_trips (bikedb = "testdb", city = "ny")
bike_daily_trips (bikedb = "testdb", city = "ny", member = TRUE)
bike_daily_trips (bikedb = "testdb", city = "ny", gender = "f")
bike_daily_trips (bikedb = "testdb", city = "ny", station = "173",
                  gender = 1)

bike_rm_test_data ()
bike_rm_db ("testdb")
# don't forget to remove real data!
# file.remove (list.files (".", pattern = ".zip"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database-stats.R
\name{bike_latest_files}
\alias{bike_latest_files}
\title{Check whether files in database are the latest published files}
\usage{
bike_latest_files(bikedb)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.
If no directory specified, it is presumed to be in \code{tempdir()}.}
}
\value{
A named vector of binary values: TRUE is files in \code{bikedb} are
the latest versions; otherwise FALSE, in which case \code{store_bikedata}
could be run to update the database.
}
\description{
Check whether files in database are the latest published files
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
# or download some real data!
# dl_bikedata (city = 'la', data_dir = data_dir)
# Remove one London file that triggers an API call which may fail tests:
file.remove (file.path (tempdir(),
             "01aJourneyDataExtract10Jan16-23Jan16.csv"))
bikedb <- file.path (data_dir, 'testdb')
store_bikedata (data_dir = data_dir, bikedb = bikedb)
# bike_latest_files (bikedb)
# All false because test data are not current, but would pass with real data

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files (data_dir, pattern = '.zip'))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/database-stats.R
\name{bike_stored_files}
\alias{bike_stored_files}
\title{Get names of files read into database}
\usage{
bike_stored_files(bikedb, city)
}
\arguments{
\item{bikedb}{A string containing the path to the SQLite3 database.}

\item{city}{Optional city for which filenames are to be obtained}
}
\description{
Get names of files read into database
}
\examples{
\dontrun{
data_dir <- tempdir ()
bike_write_test_data (data_dir = data_dir)
bikedb <- file.path (data_dir, 'testdb')
store_bikedata (data_dir = data_dir, bikedb = bikedb)
files <- bike_stored_files (bikedb = bikedb)
# returns a tibble with names of all stored files

bike_rm_test_data (data_dir = data_dir)
bike_rm_db (bikedb)
# don't forget to remove real data!
# file.remove (list.files ('.', pattern = '.zip'))
}
}
