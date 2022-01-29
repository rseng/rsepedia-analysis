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

