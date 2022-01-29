
# weathercan <img src="https://github.com/ropensci/weathercan/raw/master/inst/assets/weathercan_logo.png" align = "right" width = 110/>

[![:name status
badge](https://ropensci.r-universe.dev/badges/:name)](https://ropensci.r-universe.dev)
[![weathercan status
badge](https://ropensci.r-universe.dev/badges/weathercan)](https://ropensci.r-universe.dev)
[![R-CMD-check](https://github.com/ropensci/weathercan/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/weathercan/actions)
[![codecov](https://codecov.io/gh/ropensci/weathercan/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/weathercan)

[![](https://badges.ropensci.org/160_status.svg)](https://github.com/ropensci/software-review/issues/160)
[![DOI](https://zenodo.org/badge/60650396.svg)](https://zenodo.org/badge/latestdoi/60650396)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00571/status.svg)](https://doi.org/10.21105/joss.00571)

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/weathercan)](https://cran.r-project.org/package=weathercan)
[![CRAN
Downloads](http://cranlogs.r-pkg.org/badges/grand-total/weathercan)](https://CRAN.R-project.org/package=weathercan)

This package makes it easier to search for and download multiple
months/years of historical weather data from [Environment and Climate
Change Canada (ECCC)
website](https://climate.weather.gc.ca/historical_data/search_historic_data_e.html).

Bear in mind that these downloads can be fairly large and performing
multiple downloads may use up ECCC’s bandwidth unnecessarily. Try to
stick to what you need.

For more details and tutorials checkout the [weathercan
website](https://docs.ropensci.org/weathercan/)

> Check out the Demo weathercan shiny dashboard
> ([html](https://steffilazerte.shinyapps.io/weathercan_shiny/);
> [source](https://github.com/steffilazerte/weathercan_shiny))

## Installation

You can install `weathercan` directly from CRAN:

``` r
install.packages("weathercan")
```

Or you can install from the rOpenSci R-Universe:

``` r
install.packages("weathercan", repos = "https://ropensci.r-universe.dev")
```

View the available vignettes with `vignette(package = "weathercan")`

View a particular vignette with, for example,
`vignette("weathercan", package = "weathercan")`

## General usage

To download data, you first need to know the `station_id` associated
with the station you’re interested in.

### Stations

`weathercan` includes the function `stations()` which returns a list of
stations and their details (including `station_id`).

``` r
head(stations())
```

    ## # A tibble: 6 × 16
    ##   prov  station_name        station_id climate_id WMO_id TC_id   lat   lon  elev tz        interval start   end normals normals_1981_2010 normals_1971_2000
    ##   <chr> <chr>                    <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <chr>    <dbl> <dbl> <lgl>   <lgl>             <lgl>            
    ## 1 AB    DAYSLAND                  1795 301AR54        NA <NA>   52.9 -112.  689. Etc/GMT+7 day       1908  1922 FALSE   FALSE             FALSE            
    ## 2 AB    DAYSLAND                  1795 301AR54        NA <NA>   52.9 -112.  689. Etc/GMT+7 hour        NA    NA FALSE   FALSE             FALSE            
    ## 3 AB    DAYSLAND                  1795 301AR54        NA <NA>   52.9 -112.  689. Etc/GMT+7 month     1908  1922 FALSE   FALSE             FALSE            
    ## 4 AB    EDMONTON CORONATION       1796 301BK03        NA <NA>   53.6 -114.  671. Etc/GMT+7 day       1978  1979 FALSE   FALSE             FALSE            
    ## 5 AB    EDMONTON CORONATION       1796 301BK03        NA <NA>   53.6 -114.  671. Etc/GMT+7 hour        NA    NA FALSE   FALSE             FALSE            
    ## 6 AB    EDMONTON CORONATION       1796 301BK03        NA <NA>   53.6 -114.  671. Etc/GMT+7 month     1978  1979 FALSE   FALSE             FALSE

``` r
glimpse(stations())
```

    ## Rows: 26,337
    ## Columns: 16
    ## $ prov              <chr> "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", "AB", …
    ## $ station_name      <chr> "DAYSLAND", "DAYSLAND", "DAYSLAND", "EDMONTON CORONATION", "EDMONTON CORONATION", "EDMONTON CORONATION", "FLEET", "FLEET", "FLEET", …
    ## $ station_id        <dbl> 1795, 1795, 1795, 1796, 1796, 1796, 1797, 1797, 1797, 1798, 1798, 1798, 1799, 1799, 1799, 1800, 1800, 1800, 1801, 1801, 1801, 1802, …
    ## $ climate_id        <chr> "301AR54", "301AR54", "301AR54", "301BK03", "301BK03", "301BK03", "301B6L0", "301B6L0", "301B6L0", "301B8LR", "301B8LR", "301B8LR", …
    ## $ WMO_id            <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ TC_id             <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
    ## $ lat               <dbl> 52.87, 52.87, 52.87, 53.57, 53.57, 53.57, 52.15, 52.15, 52.15, 53.20, 53.20, 53.20, 52.40, 52.40, 52.40, 54.08, 54.08, 54.08, 53.52,…
    ## $ lon               <dbl> -112.28, -112.28, -112.28, -113.57, -113.57, -113.57, -111.73, -111.73, -111.73, -110.15, -110.15, -110.15, -115.20, -115.20, -115.2…
    ## $ elev              <dbl> 688.8, 688.8, 688.8, 670.6, 670.6, 670.6, 838.2, 838.2, 838.2, 640.0, 640.0, 640.0, 1036.0, 1036.0, 1036.0, 585.2, 585.2, 585.2, 668…
    ## $ tz                <chr> "Etc/GMT+7", "Etc/GMT+7", "Etc/GMT+7", "Etc/GMT+7", "Etc/GMT+7", "Etc/GMT+7", "Etc/GMT+7", "Etc/GMT+7", "Etc/GMT+7", "Etc/GMT+7", "E…
    ## $ interval          <chr> "day", "hour", "month", "day", "hour", "month", "day", "hour", "month", "day", "hour", "month", "day", "hour", "month", "day", "hour…
    ## $ start             <dbl> 1908, NA, 1908, 1978, NA, 1978, 1987, NA, 1987, 1987, NA, 1987, 1980, NA, 1980, 1980, NA, 1980, 1986, NA, 1986, 1987, NA, 1987, 1986…
    ## $ end               <dbl> 1922, NA, 1922, 1979, NA, 1979, 1990, NA, 1990, 1998, NA, 1998, 2009, NA, 2007, 1981, NA, 1981, 2019, NA, 2007, 1991, NA, 1991, 1995…
    ## $ normals           <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRU…
    ## $ normals_1981_2010 <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, TRUE, TRUE, TRUE, FALSE, FALSE, FALSE, TRUE, TRU…
    ## $ normals_1971_2000 <lgl> FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE, FALSE,…

You can look through this data frame directly, or you can use the
`stations_search` function:

``` r
stations_search("Kamloops", interval = "hour")
```

    ## # A tibble: 3 × 16
    ##   prov  station_name station_id climate_id WMO_id TC_id   lat   lon  elev tz        interval start   end normals normals_1981_2010 normals_1971_2000
    ##   <chr> <chr>             <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <chr>    <dbl> <dbl> <lgl>   <lgl>             <lgl>            
    ## 1 BC    KAMLOOPS A         1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      1953  2013 TRUE    TRUE              TRUE             
    ## 2 BC    KAMLOOPS A        51423 1163781     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      2013  2021 FALSE   FALSE             FALSE            
    ## 3 BC    KAMLOOPS AUT      42203 1163842     71741 ZKA    50.7 -120.  345  Etc/GMT+8 hour      2006  2021 FALSE   FALSE             FALSE

Time frame must be one of “hour”, “day”, or “month”.

You can also search by proximity:

``` r
stations_search(coords = c(50.667492, -120.329049), dist = 20, interval = "hour")
```

    ## # A tibble: 3 × 17
    ##   prov  station_name station_id climate_id WMO_id TC_id   lat   lon  elev tz        interval start   end normals normals_1981_2010 normals_1971_2000 distance
    ##   <chr> <chr>             <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <chr>    <dbl> <dbl> <lgl>   <lgl>             <lgl>                <dbl>
    ## 1 BC    KAMLOOPS A         1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      1953  2013 TRUE    TRUE              TRUE                  8.64
    ## 2 BC    KAMLOOPS AUT      42203 1163842     71741 ZKA    50.7 -120.  345  Etc/GMT+8 hour      2006  2021 FALSE   FALSE             FALSE                 8.64
    ## 3 BC    KAMLOOPS A        51423 1163781     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      2013  2021 FALSE   FALSE             FALSE                 9.28

You can update this list of stations with

``` r
stations_dl()
```

    ## According to Environment Canada, Modified Date: 2021-10-31 23:34 UTC

    ## Stations data saved...
    ## Use `stations()` to access most recent version and `stations_meta()` to see when this was last updated

And check when it was last updated with

``` r
stations_meta()
```

    ## $ECCC_modified
    ## [1] "2021-10-31 23:34:00 UTC"
    ## 
    ## $weathercan_modified
    ## [1] "2021-11-30"

**Note:** For reproducibility, if you are using the stations list to
gather your data, it can be a good idea to take note of the ECCC date of
modification and include it in your reports/manuscripts.

### Weather

Once you have your `station_id`(s) you can download weather data:

``` r
kam <- weather_dl(station_ids = 51423, start = "2018-02-01", end = "2018-04-15")
```

    ## As of weathercan v0.3.0 time display is either local time or UTC
    ## See Details under ?weather_dl for more information.
    ## This message is shown once per session

``` r
kam
```

    ## # A tibble: 1,776 × 37
    ##    station_name station_id station_operator prov    lat   lon  elev climate_id WMO_id TC_id date       time                year  month day   hour  weather  hmdx
    ##    <chr>             <dbl> <lgl>            <chr> <dbl> <dbl> <dbl> <chr>      <chr>  <chr> <date>     <dttm>              <chr> <chr> <chr> <chr> <chr>   <dbl>
    ##  1 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 00:00:00 2018  02    01    00:00 <NA>       NA
    ##  2 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 01:00:00 2018  02    01    01:00 Snow       NA
    ##  3 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 02:00:00 2018  02    01    02:00 <NA>       NA
    ##  4 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 03:00:00 2018  02    01    03:00 <NA>       NA
    ##  5 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 04:00:00 2018  02    01    04:00 Cloudy     NA
    ##  6 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 05:00:00 2018  02    01    05:00 <NA>       NA
    ##  7 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 06:00:00 2018  02    01    06:00 <NA>       NA
    ##  8 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 07:00:00 2018  02    01    07:00 Cloudy     NA
    ##  9 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 08:00:00 2018  02    01    08:00 <NA>       NA
    ## 10 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2018-02-01 2018-02-01 09:00:00 2018  02    01    09:00 <NA>       NA
    ## # … with 1,766 more rows

You can also download data from multiple stations at once:

``` r
kam_pg <- weather_dl(station_ids = c(48248, 51423), start = "2018-02-01", end = "2018-04-15")
```

## Climate Normals

To access climate normals, you first need to know the `climate_id`
associated with the station you’re interested in.

``` r
stations_search("Winnipeg", normals_years = "current")
```

    ## # A tibble: 1 × 13
    ##   prov  station_name                station_id climate_id WMO_id TC_id   lat   lon  elev tz        normals normals_1981_2010 normals_1971_2000
    ##   <chr> <chr>                            <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <lgl>   <lgl>             <lgl>            
    ## 1 MB    WINNIPEG RICHARDSON INT'L A       3698 5023222     71852 YWG    49.9 -97.2  239. Etc/GMT+6 TRUE    TRUE              TRUE

Then you can download the climate normals with the `normals_dl()`
function.

``` r
n <- normals_dl("5023222")
```

See the [Getting
Started](https://docs.ropensci.org/weathercan/articles/weathercan.html)
vignette for more details.

## Citation

``` r
citation("weathercan")
```

    ## 
    ## To cite 'weathercan' in publications, please use:
    ## 
    ##   LaZerte, Stefanie E and Sam Albers (2018). weathercan: Download and format weather data from Environment and Climate Change Canada. The
    ##   Journal of Open Source Software 3(22):571. doi:10.21105/joss.00571.
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Article{,
    ##     title = {{weathercan}: {D}ownload and format weather data from Environment and Climate Change Canada},
    ##     author = {Stefanie E LaZerte and Sam Albers},
    ##     journal = {The Journal of Open Source Software},
    ##     volume = {3},
    ##     number = {22},
    ##     pages = {571},
    ##     year = {2018},
    ##     url = {https://joss.theoj.org/papers/10.21105/joss.00571},
    ##   }

## License

The data and the code in this repository are licensed under multiple
licences. All code is licensed
[GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html). All weather data
is licensed under the ([Open Government License -
Canada](http://open.canada.ca/en/open-government-licence-canada)).

## `weathercan` in the wild!

-   Browse [`weathercan` use cases](https://ropensci.org/usecases/) on
    rOpenSci.org
-   Checkout the [`weathercan` Shiny
    App](https://nickrongkp.shinyapps.io/WeatherCan/) by Nick Rong
    (@nickyrong) and Nathan Smith (@WraySmith)
-   R package
    [`RavenR`](https://github.com/rchlumsk/RavenR/tree/master/R) has
    functions for converting ECCC data downloaded by `weathercan` to the
    .rvt format for Raven.
-   R package [`meteoland`](https://github.com/emf-creaf/meteoland) has
    functions for converting ECCC data downloaded by `weathercan` to the
    format required for use in `meteoland`.

## Similar packages

**[`rclimateca`](https://github.com/paleolimbot/rclimateca)**

`weathercan` and `rclimateca` were developed at roughly the same time
and as a result, both present up-to-date methods for accessing and
downloading data from ECCC. The largest differences between the two
packages are: a) `weathercan` includes functions for interpolating
weather data and directly integrating it into other data sources. b)
`weathercan` actively seeks to apply tidy data principles in R and
integrates well with the tidyverse including using tibbles and nested
listcols. c) `rclimateca` contains arguments for specifying short
vs. long data formats. d) `rclimateca` has the option of formatting data
in the MUData format using the
[`mudata2`](https://cran.r-project.org/package=mudata2) package by the
same author.

**[`CHCN`](https://cran.r-project.org/package=CHCN)**

`CHCN` is an older package last updated in 2012. Unfortunately, ECCC
updated their services within the last couple of years which caused a
great many of the previous web scrapers to fail. `CHCN` relies on a
decommissioned [older web-scraper](https://quickcode.io/) and so is
currently broken.

## Contributions

We welcome any and all contributions! To make the process as painless as
possible for all involved, please see our [guide to
contributing](CONTRIBUTING.md)

## Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By participating in
this project you agree to abide by its terms.

[![ropensci\_footer](http://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# weathercan 0.6.2

- Create cache dir for stations data recursively
- Fix choice of local vs. package version of stations data frame
- Update to readr v2
- Add flexibility for csv/tsv stations files (fixes #126)
- Update stations url
- Make examples and tests robust to internet issues


# weathercan 0.6.1

## Small changes
- Save `stations()` data to local cache

# weathercan 0.6.0

## Big changes
- Move from data frame `stations` to function `stations()`. Returns same data
  but is updateable with `stations_dl()` and you can check download dates
  version with `stations_meta()` (fixes #10)
- Download climate normals from climate.weather.gc.ca (fixes #88)
  - More stations available (more than 2x as many!)
  - More year ranges available (1981-2010 and 1971-2000; 
    Note that while climate normals from 1961-1990 are available from ECCC, they
    don't have climate ids making it tricky to download reliably)
    
## Small changes
- Remove old deprecated function arguments
- Better test coverage (#81)
- Better handling of http errors (#101, #119; Thanks @KevCaz!)

## Bug fixes
- Download stations data frame from ECCC Google drive rather than ECCC FTP site
- Update dependency versions (#111, #112, #118)


# weathercan 0.5.0 (2020-01-14)

## Small changes
- Internal changes to address change in formatting of historical weather data provided by ECCC 
  (includes new parameters for the amount of precipitation in mm: `precip_amt`, `precip_amt_flag`; fixes #107) 
- Updated stations data frame

## Bug fixes
- Updated normals column values (fixes #106)

# weathercan 0.4.0 (2020-08-26)

## Bug fixes
- Fixed odd bug where some Linux systems failed to download stations data

## Features and potentially breaking changes
- Added caching in memory with memoise (caches for 24hrs, can change this by restarting the R session)
    - Caches individual downloaded files, so you may see a speed up even if you change the parameters of the download.
- Some missing values in meta data previously were "" but are now explicitly NAs

## Internal changes
- Use readr for reading data
- Use vcr for tests

# weathercan 0.3.4 (2020-04-14)

## Small changes
- Internal changes to fix compatibility with tibble v3.0.0
- Internal changes to fix compatibility with dplyr v1.0.0
- Updated internal stations data frame

# weathercan 0.3.3 (2020-01-24)

## Small changes
- Internal changes to address issues with testing
- Remove all reliance on ECCC servers when testing on CRAN
- Update internal datasets

# weathercan 0.3.2 (2020-01-06)

## Small changes
- Internal changes to address expected changes to normals metadata
- Internal changes to address problems with connections on Windows
- Update links to website

# weathercan 0.3.1 (2019-09-27)

## Small changes
- Internal changes to address change in formatting of historical weather data provided by ECCC (fixes #83)


# weathercan 0.3.0 (2019-09-25)

## Big changes
- New function: `normals_dl()` function downloads climate normals. Addresses issue #38.
- New argument: `stations_search()` has `normals_only` to return only stations with climate normals
- Deprecated `url` argument in favour of `weathercan.urls.stations`, `weathercan.urls.weather` and `weathercan.urls.normals` options.
- Deprecated `tz_disp` in favour of `time_disp`. Now all timezones are UTC, but the displayed time is either local time or UTC. When `time_disp` = "none", the time displayed is local time without daylight savings, similar to how ECCC presents the data. This means that data from different time zones will have similar ecological times (i.e. midnights will be comparable), but the actual times are not UTC. When `time_disp` = "UTC', the time displayed is in UTC timezone, meaning that stations from different times zones will have true times (i.e. midnight observation in Toronto will be three hours before midnight observation in Vancouver). Addresses issue #74.

## Small changes
- Add parameter in `station_search()` to restrict by start and end dates. Addresses issue #35.
- Internal change, switching to .data and "" for all non-standard evaluations as opposed to listing global variables
- Tweaks to keep compatibility wit `tidyr`

## Bug fixes
- Fix bug #69 which resulted in daily downloads missing partial years when the date range spanned two calendar years
- Fix bug #70 where internal `stations` data frame references conflicted with local references to `stations`
- Fix bug #72 which was a security vulnerability in an article's json

# weathercan 0.2.8 (2018-10-08)

## Bug fixes
- Add timezones to the `stations` data frame to remove dependency of Google API. Timezones added with the `lutz` package, so updates the the `stations` data frame now require `lutz` and `sf` packages (added to Suggests).

## Changes
- Sort `stations` by `station_id` not by `station_name`

## Other
- Update all internal data frames

# weathercan 0.2.7 (2018-06-27)

## Bug fixes
- Fix bug created when ECCC changed file metadata for dates after April 1st 2018 (only affected downloads which included dates both before AND after April 1st, 2018) - Results in a new column `station_operator` for all data (NA where unavailable for older stations).
- Adjust code flexibility to handle future changes
- Add catch to warn user if end dates earlier than start dates

## Changes
- Update readme/vignettes/internal data sets to include new columns
- Update internal `stations` data frame as well as `flags` and `glossary`
- Remove `tibble` dependency by relying on `dplyr`

# weathercan 0.2.6 (2018-05-25)

## Bug fixes
- Fix bug created when ECCC removed Data Quality Column
- Adjust code flexibility to handle future changes
- Add tests to catch future changes

# weathercan 0.2.5 (2018-03-02)

## Changes
- More sensible messages when missing station data
- Streamline messages from multiple stations
- Accepts older R version
- `stations_dl` fails gracefully on R versions < 3.3.4
- Update `stations` dataframe

## Bug fixes
- Fix error when missing station data from one of several stations

# weathercan 0.2.4 (2018-02-01)

Now part of [ropensci.org](https://ropensci.org)!

## Changes
- `sp` moved to suggests, users are now prompted to install sp if they want to search stations by coordinates
- `weather_dl()` replaces `weather()`
- `weather_interp()` replaces `add_weather()`
- `stations_dl()` replaces `stations_all()`
- `tz_calc()` replaces `get_tz()`
- Internal code modifications to match best practices

# weathercan 0.2.3 (2017-11-22)

## Changes
- Updated `stations` data
- Added `flags` and `glossary` datasets as well as vignettes
- `stations_search()` warns user if name looks like coords
- `stations_search()` with `coord` now returns closest 10 stations
- `add_weather()` warns user if trying to interpolate weather from >1 station
- Updated code to conform with rOpenSci requirements
- Data downloaded from multiple timezones defaults to UTC

## Bug fixes
- `weather(format = FALSE)` properly returns data
- updated `weather()` to work with `lubridate` 1.7.1

# weathercan 0.2.2 (2017-06-16)

## Changes
- Update and expand vignettes (closes #15)
- Data now returned as tibbles
- Added listcol functionality (closes #8)
- Added internal tests for interpolation
- Updated R version
- Standardized reference to stations dataset (`stn`) in all functions

## Major changes
- envirocan renamed to weathercan (closes #17)

## Bug fixes
- Fixed inclusion of New Brunswick stations (closes #9)
- Downloads with no data return empty tibble and an informative message (closes #21)


# envirocan 0.2.1 (2017-03-04)
- Minor bug fixes: correcting encoding information for downloads, updating function calls to dplyr package, updating stations dataset

# envirocan 0.2.0 (2016-07-08)

- Added new function, `add_weather()` which performs a linear interpolation and merges weather data into an existing data frame.
- Added two new hourly datasets with weather data downloaded for Kamloops and Prince George, BC: kamloops, pg
- Added a new daily dataset for Kamloops: kamloops_day
- Fixed a bug when downloading data from multiple stations at the same time
- Changed 'timeframe' arguments to 'interval'
- Minor internal adjustments


# envirocan 0.1.1.1 (2016-06-23)

- quick fix to correct duplicated monthly data downloads


# envirocan 0.1.1 (2016-06-23)

## Functionality
- Allow blank start/end dates to download whole data set
- Add option to trim missing data from start and end of the range

## Bug fixes
- Add messages so functions fail gracefully if timezone doesn't exist
- Correct bugs that prevented downloading of monthly data


# envirocan 0.1.0 (2016-06-21)

This is the initial release for envirocan.

## Include functionality:

Finding stations:

- `stations` data frame to look up station data
- `stations_search()` function to search for a station by name or proximity
- `stations_all()` to download a new stations data set from Environment Canada

Downloading weather:

- `weather()` function to specify station_id(s) start and end dates for downloading data.

# Contributing to `weathercan`

Thank you for any and all contributions! Following these guidelines will help streamline the process of contributing and make sure that we're all on the same page. While we ask that you read this guide and follow it to the best of your abilities, we welcome contributions from all, regardless of your level of experience.

By participating in this project, you agree to abide by the [code of conduct](https://github.com/ropensci/weathercan/blob/master/CONDUCT.md).

# Types of contributions 

Don't feel that you must be a computer whiz to make meaningful contributions. Feel free to:

- Identify areas for future development ([open an Issue](https://github.com/ropensci/weathercan/issues))
- Identify issues/bugs ([open an Issue](https://github.com/ropensci/weathercan/issues))
- Write tutorials/vignettes ([open a Pull Request](https://github.com/ropensci/weathercan/pulls) to contribute to the ones here, or make your own elsewhere and send us a link)
- Add functionality ([open a Pull Request](https://github.com/ropensci/weathercan/pulls))
- Fix bugs ([open a Pull Request](https://github.com/ropensci/weathercan/pulls))

# New to GitHub?

Getting ready to make your first contribution? Here are a couple of tutorials you may wish to check out:

- [Tutorial for first-timers](https://github.com/Roshanjossey/first-contributions)
- [How to contribute (in-depth lessons)](https://egghead.io/series/how-to-contribute-to-an-open-source-project-on-github)
- [GitHub on setup](https://help.github.com/articles/set-up-git)
- [GitHub on pull requests](https://help.github.com/articles/using-pull-requests/).)


# How to contribute code

- Fork the repository
- Clone the repository from GitHub to your computer e.g,. `git clone https://github.com/ropensci/weathercan.git`
- Make sure to track progress upstream (i.e., on our version of `weathercan` at `ropensci/weathercan`)
  - `git remote add upstream https://github.com/ropensci/weathercan.git`
  - Before making changes make sure to pull changes in from upstream with `git pull upstream`
- Make your changes
  - For changes beyond minor typos, add an item to NEWS.md describing the changes and add yourself to the DESCRIPTION file as a contributor
- Push to your GitHub account
- Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/weathercan`

# Code formatting

- In general follow the convention of <http://r-pkgs.had.co.nz/r.html#style> (snake_case functions and argument names, etc.)
- Where there is conflict, default to the style of `weathercan`
- Use explicit package imports (i.e. package_name::package_function) and avoid @import if at all possible
## Release v.0.6.2

* General changes to avoid CRAN policy violations, namely, errors when API is inaccessible:
  * Use @examplesIf to ONLY run examples if resources are available 
     (I believe this was the problem that resulted in archival)
  * Avoid tests where API calls used (even if mocking, or using vcr)
  * Pre-compiled the vignettes which rely on the API
* Fixed minor bugs with new stations data
* Increased robustness to ECCC data changes

## Test environments
As of November 30th, 2021

* ubuntu 20.04 - Local (4.1.1), GitHub Actions (devel, release, old release)
* Windows Server - GitHub Actions (release), winbuilder (devel), rhub (devel, release, old release)
* OSX 11.6.1 - GitHub Actions (release)
* Solaris - rhub (release)
* Fedora GCC - rhub (devel)
* Fedora CLANG - rhub (devel)
* Debian CLANG - rhub (devel)
* Debian GCC - rhub (release, patched, devel)


## R CMD check results

There were no ERRORs and no WARNINGs

In addition to the NOTE that this is a new submission:

1 NOTEs:

Solaris had one NOTEs:
* checking top-level files ... NOTE
  Files ‘README.md’ or ‘NEWS.md’ cannot be checked without ‘pandoc’ being installed.

## Downstream dependencies

* RavenR suggests weathercan, I checked and found no problems
<!--- Provide a general summary of your changes in the Title above -->

## Description
<!--- Describe your changes in detail -->

## Related Issue
<!--- if this closes an issue make sure include e.g., "fix #4"
or similar - or if just relates to an issue make sure to mention
it like "#4" -->

## Example
<!--- if introducing a new feature or changing behavior of existing
methods/functions, include an example if possible to do in brief form -->

## Best Practices
<!--- Did you remember to include documentation, examples andtests? Unless you're just changing
grammar, please include new tests for your change -->
The following have been updated or added as needed:
[ ] Documentation
[ ] Examples in documentation
[ ] Vignettes
[ ] `testthat` Tests
## Expected Behavior
<!--- If you're describing a bug, tell us what should happen -->
<!--- If you're suggesting a change/improvement, tell us how it should work -->

## Current Behavior
<!--- If describing a bug, tell us what happens instead of the expected behavior -->
<!--- If suggesting a change/improvement, explain the difference from current behavior -->

## Steps to Reproduce (for bugs)
<!--- Provide a link to a live example, or an unambiguous set of steps to -->
<!--- reproduce this bug. Include code to reproduce, if relevant -->
1.
2.
3.
4.

## Possible Solution
<!--- Not obligatory, but suggest a fix/reason for the bug, -->
<!--- or ideas how to implement the addition or change -->

## Context
<!--- How has this issue affected you? What are you trying to accomplish? -->
<!--- Providing context helps us come up with a solution that is most useful in the real world -->

## Your Environment
<!--- Include the output of "devtools::session_info()" to help us understand your system environment -->
## Environment and Climate Change Canada Data

The file `stations.csv` contains a dataset containing station information downloaded from Environment and Climate Change Canada, downloaded from [Environment and Climate Change Canada](http://dd.weather.gc.ca/citypage_weather/docs/site_list_en.csv) under the ([Open Government License - Canada](http://open.canada.ca/en/open-government-licence-canada)). 
---
title: 'weathercan: Download and format weather data from Environment and Climate Change Canada'
tags:
  - R
  - data
  - weather
  - Canada
  - meteorology
authors:
  - name: Stefanie E. LaZerte
    orcid: 0000-0002-7690-8360
    affiliation: 1
  - name: Sam Albers
    orcid: 
    affiliation: 2
affiliations:
  - name: steffilazerte.ca
    index: 1
  - name: University of Northern British Columbia
    index: 2
date: 2017-11-24
bibliography: paper.bib
---

# Summary

Environment and Climate Change Canada maintains an online source of historical Canadian weather data in hourly, daily and monthly formats for various stations across Canada [@canada_historical_2011]. This data is freely available and can be accessed directly from their website. However, downloading data from multiple stations and across larger time periods can take significant time and effort. Further, these downloads require processing before they can be used for analysis. `weathercan` [@lazerte_weathercan_2018] is an R [@r_stats] package that automates and simplifies the downloading and formating of this data.

The first step in using `weathercan` is to identify the station ID(s) of the weather station(s) of interest. Stations can be searched for either by name or proximity to a given location. Searches can be conducted on all possible stations, or filtered to include only those recording weather at the desired time interval. Next, weather data can be downloaded for the specified stations, time range and time interval (i.e. hours, days, months). Data downloaded from multiple stations and over several months are automatically combined into one data frame ready for analysis or plotting (Figure 1). Finally, weather data from a single station can be aligned and merged with existing datasets through linear interpolation.

![](paper_files/figure-markdown/unnamed-chunk-2-1.png)
Figure 1. Data downloaded with `weathercan` is formated and ready for ploting.

`weathercan` is available on GitHub at <https://github.com/ropensci/weathercan>

# References
---
title: 'weathercan: Download and format weather data from Environment and Climate Change Canada'
tags:
  - R
  - data
  - weather
  - Canada
  - meteorology
authors:
  - name: Stefanie E. LaZerte
    orcid: 0000-0002-7690-8360
    affiliation: 1
  - name: Sam Albers
    orcid: 
    affiliation: 2
affiliations:
  - name: steffilazerte.ca
    index: 1
  - name: University of Northern British Columbia
    index: 2
date: 2017-11-24
bibliography: paper.bib
---
