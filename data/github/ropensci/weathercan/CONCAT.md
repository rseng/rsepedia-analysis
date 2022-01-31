
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
---
output: github_document
---

```{r, echo=FALSE, message=FALSE, warning=FALSE}
library(weathercan)
library(dplyr)
library(tibble)
knitr::opts_chunk$set(cache = FALSE,
                      fig.path = "tools/readme/")
old <- options(width = 160)
```

# weathercan <img src="https://github.com/ropensci/weathercan/raw/master/inst/assets/weathercan_logo.png" align = "right" width = 110/>

[![:name status badge](https://ropensci.r-universe.dev/badges/:name)](https://ropensci.r-universe.dev)
[![weathercan status badge](https://ropensci.r-universe.dev/badges/weathercan)](https://ropensci.r-universe.dev)
[![R-CMD-check](https://github.com/ropensci/weathercan/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/weathercan/actions)
[![codecov](https://codecov.io/gh/ropensci/weathercan/branch/master/graph/badge.svg)](https://app.codecov.io/gh/ropensci/weathercan)

[![](https://badges.ropensci.org/160_status.svg)](https://github.com/ropensci/software-review/issues/160) [![DOI](https://zenodo.org/badge/60650396.svg)](https://zenodo.org/badge/latestdoi/60650396) [![DOI](http://joss.theoj.org/papers/10.21105/joss.00571/status.svg)](https://doi.org/10.21105/joss.00571)


[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/weathercan)](https://cran.r-project.org/package=weathercan) [![CRAN Downloads](http://cranlogs.r-pkg.org/badges/grand-total/weathercan)](https://CRAN.R-project.org/package=weathercan)


This package makes it easier to search for and download multiple months/years of historical weather data from [Environment and Climate Change Canada (ECCC) website](https://climate.weather.gc.ca/historical_data/search_historic_data_e.html).

Bear in mind that these downloads can be fairly large and performing multiple downloads may use up ECCC's bandwidth unnecessarily. Try to stick to what you need.

For more details and tutorials checkout the [weathercan website](https://docs.ropensci.org/weathercan/)

> Check out the Demo weathercan shiny dashboard ([html](https://steffilazerte.shinyapps.io/weathercan_shiny/); [source](https://github.com/steffilazerte/weathercan_shiny))

## Installation

You can install `weathercan` directly from CRAN:

```{r, eval = FALSE}
install.packages("weathercan")
```


Or you can install from the rOpenSci R-Universe:

```{r, eval = FALSE}
install.packages("weathercan", repos = "https://ropensci.r-universe.dev")
```


View the available vignettes with `vignette(package = "weathercan")`  
 
View a particular vignette with, for example, `vignette("weathercan", package = "weathercan")`

## General usage

To download data, you first need to know the `station_id` associated with the station you're interested in.

### Stations

`weathercan` includes the function `stations()` which returns a list of stations and their details (including `station_id`).

```{r}
head(stations())
glimpse(stations())
```

You can look through this data frame directly, or you can use the `stations_search` function:

```{r}
stations_search("Kamloops", interval = "hour")
```

Time frame must be one of "hour", "day", or "month".

You can also search by proximity:

```{r}
stations_search(coords = c(50.667492, -120.329049), dist = 20, interval = "hour")
```

You can update this list of stations with 

```{r}
stations_dl()
```

And check when it was last updated with
```{r}
stations_meta()
```

**Note:** For reproducibility, if you are using the stations list to gather your
data, it can be a good idea to take note of the ECCC date of modification and 
include it in your reports/manuscripts.

### Weather

Once you have your `station_id`(s) you can download weather data:

```{r, R.options = list(tibble.max_extra_cols = 0)}
kam <- weather_dl(station_ids = 51423, start = "2018-02-01", end = "2018-04-15")
kam
```

You can also download data from multiple stations at once:

```{r, R.options = list(tibble.max_extra_cols = 0)}
kam_pg <- weather_dl(station_ids = c(48248, 51423), start = "2018-02-01", end = "2018-04-15")
```

## Climate Normals

To access climate normals, you first need to know the `climate_id` associated with the station you're interested in.

```{r}
stations_search("Winnipeg", normals_years = "current")
```

Then you can download the climate normals with the `normals_dl()` function.

```{r}
n <- normals_dl("5023222")
```

See the [Getting Started](https://docs.ropensci.org/weathercan/articles/weathercan.html) 
vignette for more details. 


## Citation

```{r, warning = FALSE}
citation("weathercan")
```

## License

The data and the code in this repository are licensed under multiple licences. All code is licensed [GPL-3](https://www.gnu.org/licenses/gpl-3.0.en.html). All weather data is licensed under the ([Open Government License - Canada](http://open.canada.ca/en/open-government-licence-canada)). 

## `weathercan` in the wild!

- Browse [`weathercan` use cases](https://ropensci.org/usecases/) on rOpenSci.org
- Checkout the [`weathercan` Shiny App](https://nickrongkp.shinyapps.io/WeatherCan/) by Nick Rong (@nickyrong) and Nathan Smith (@WraySmith)
- R package [`RavenR`](https://github.com/rchlumsk/RavenR/tree/master/R) has functions for converting 
  ECCC data downloaded by `weathercan` to the .rvt format for Raven.
- R package [`meteoland`](https://github.com/emf-creaf/meteoland) has functions for converting ECCC
  data downloaded by `weathercan` to the format required for use in `meteoland`.

## Similar packages

**[`rclimateca`](https://github.com/paleolimbot/rclimateca)**

`weathercan` and `rclimateca` were developed at roughly the same time and as a result, both present up-to-date methods for accessing and downloading data from ECCC. The largest differences between the two packages are: a) `weathercan` includes functions for interpolating weather data and directly integrating it into other data sources. b) `weathercan` actively seeks to apply tidy data principles in R and integrates well with the tidyverse including using tibbles and nested listcols. c) `rclimateca` contains arguments for specifying short vs. long data formats. d) `rclimateca` has the option of formatting data in the MUData format using the [`mudata2`](https://cran.r-project.org/package=mudata2) package by the same author.

**[`CHCN`](https://cran.r-project.org/package=CHCN)**

`CHCN` is an older package last updated in 2012. Unfortunately, ECCC updated their services within the last couple of years which caused a great many of the previous web scrapers to fail. `CHCN` relies on a decommissioned [older web-scraper](https://quickcode.io/) and so is currently broken. 

## Contributions

We welcome any and all contributions! To make the process as painless as possible for all involved, please see our [guide to contributing](CONTRIBUTING.md)

## Code of Conduct

Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By participating in this project you agree to abide by its terms.

[![ropensci_footer](http://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)


```{r, include = FALSE}
# Reset options
options(old)
```
---
output: 
  md_document:
    pandoc_args: ["--atx-headers","--wrap=preserve"]
    variant: markdown
    includes:
      in_header: yaml.md
---

```{r, include = FALSE}
library(weathercan)
library(tidyverse)
library(viridis)

knitr::opts_chunk$set(cache = FALSE)

options(width = 90, tibble.max_extra_cols = 0)
```

# Summary

Environment and Climate Change Canada maintains an online source of historical Canadian weather data in hourly, daily and monthly formats for various stations across Canada [@canada_historical_2011]. This data is freely available and can be accessed directly from their website. However, downloading data from multiple stations and across larger time periods can take significant time and effort. Further, these downloads require processing before they can be used for analysis. `weathercan` [@lazerte_weathercan_2018] is an R [@r_stats] package that automates and simplifies the downloading and formating of this data. 

The first step in using `weathercan` is to identify the station ID(s) of the weather station(s) of interest. Stations can be searched for either by name or proximity to a given location. Searches can be conducted on all possible stations, or filtered to include only those recording weather at the desired time interval. Next, weather data can be downloaded for the specified stations, time range and time interval (i.e. hours, days, months). Data downloaded from multiple stations and over several months are automatically combined into one data frame ready for analysis or plotting (Figure 1). Finally, weather data from a single station can be aligned and merged with existing datasets through linear interpolation. 

```{r, echo = FALSE, fig.width = 5, fig.asp = 0.8, dpi = 600}
w <- weather_dl(station_ids = c(50821, 51097), 
                start = "2017-01-01", end = "2017-09-01",
                interval = "hour")

ggplot(data = w, aes(x = time, y = temp, colour = station_name)) +
  theme_bw() +
  theme(legend.position = "top") +
  geom_line() +
  labs(x = "Date", y = "Temperature C") +
  scale_colour_viridis(name = "Station", discrete = TRUE, end = 0.7)
```
Figure 1. Data downloaded with `weathercan` is formated and ready for ploting.


`weathercan` is available on GitHub at <https://github.com/ropensci/weathercan>



# References
---
title: "Climate Normals: Terms and Units"
author: "Steffi LaZerte"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Climate Normals: Terms and Units}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(weathercan)
library(dplyr)
library(tidyr)
library(stringr)
```

This table shows details regarding original column (measurement) names and units of all climate normals measurements. It further provides links back to the ECCC glossary for more details.

See the ECCC website on climate normals for more details: <`r paste0("https://www.canada.ca/en/environment-climate-change/services/",
            "climate-change/canadian-centre-climate-services/display-download/",
            "technical-documentation-climate-normals.html")`>

For details on weather measurements, see the `glossary` vignette.

### General descriptions
```{r, asis = TRUE, echo = FALSE}
glossary_normals[1:18,] %>%
  mutate(description = stringr::str_replace_all(description, "\\n", " ")) %>%
  knitr::kable()
```


### Original names and units
These represent the original ECCC measurement names with units and their corresponding measurements in `weathercan`.

```{r, echo = FALSE}
g <- glossary_normals[19:nrow(glossary_normals),] %>%
  select(-description) %>%
  mutate(group = str_detect(weathercan_name, "title"),
         group = cumsum(group)) %>%
  filter(weathercan_name != "probability") %>%
  group_by(group) %>% 
  mutate(title = ECCC_name[1]) %>%
  group_by(title) %>%
  filter(!str_detect(weathercan_name, "title")) %>%
  select(-group) %>%
  nest()
```

```{r, results = "asis", echo = FALSE}
for(t in 1:nrow(g)) {
  cat("<center><h4>", str_to_title(g$title[t]), "</h3></center>\n")
  print(knitr::kable(g$data[[t]], format = "html"))
  cat("\n")
}
```
---
title: "Climate Normals"
author: "Steffi LaZerte"
date: "2021-11-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Climate Normals}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



## Downloading Climate Normals

Climate Normals and Averages describe the average climate conditions specific to a particular location. These can be downloaded from Environment and Climate Change Canada using the `normals_dl()` function.

First we'll load the `weathercan` package for downloading the data and the `tidyr` package for unnesting the data (see below).


```r
library(weathercan)
library(tidyr)
library(dplyr)
library(naniar) # For exploring missing values
```

To download climate normals, we'll first find the stations we're interested in using the `stations_search()` function. We'll use the `normals_years = "current"` argument to filter to only stations with available climate normals for the `1981-2010` year range.


```r
stations_search("Winnipeg", normals_years = "current")
```

```
## # A tibble: 1 × 13
##   prov  station_name    station_id climate_id WMO_id TC_id   lat   lon  elev tz    normals
##   <chr> <chr>                <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr> <lgl>  
## 1 MB    WINNIPEG RICHA…       3698 5023222     71852 YWG    49.9 -97.2  239. Etc/… TRUE
```

Let's look at the climate normals from this station in Winnipeg, MB. Note that unlike the `weather_dl()` function, the `normals_dl()` function requires `climate_id` not `station_id`. By default the normals are downloaded for the years "1981-2010" (currently 1981-2010 and 1971-2000 are the only year ranges available)


```r
n <- normals_dl(climate_ids = "5023222")
n
```

```
## # A tibble: 1 × 7
##   prov  station_name                climate_id normals_years meets_wmo normals     frost  
##   <chr> <chr>                       <chr>      <chr>         <lgl>     <list>      <list> 
## 1 MB    WINNIPEG RICHARDSON INT'L A 5023222    1981-2010     TRUE      <tibble [1… <tibbl…
```

Because there are two different types of climate normals (weather measurements and first/last frost dates), the data are nested as two different datasets. We can see that the Airport (Richardson Int'l) has 197 average weather measurements/codes as well as first/last frost dates.

We can also see that this station has data quality sufficient to meet the WMO standards for temperature and precipitation (i.e. both these measurements have code >= A). See the [ECCC calculations document](https://climate.weather.gc.ca/doc/Canadian_Climate_Normals_1981_2010_Calculation_Information.pdf) for more details.

To extract either data set we can use the `unnest()` function from the `tidyr` package.


```r
normals <- unnest(n, normals)
frost <- unnest(n, frost)
```

Note that this extracts the measurements for all three stations (in the case of the `normals` data frame), but not all measurements are available for each station


```r
normals
```

```
## # A tibble: 13 × 203
##    prov  station_name          climate_id normals_years meets_wmo period temp_daily_avera…
##    <chr> <chr>                 <chr>      <chr>         <lgl>     <fct>              <dbl>
##  1 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Jan                -16.4
##  2 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Feb                -13.2
##  3 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Mar                 -5.8
##  4 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Apr                  4.4
##  5 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      May                 11.6
##  6 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Jun                 17  
##  7 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Jul                 19.7
##  8 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Aug                 18.8
##  9 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Sep                 12.7
## 10 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Oct                  5  
## 11 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Nov                 -4.9
## 12 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Dec                -13.2
## 13 MB    WINNIPEG RICHARDSON … 5023222    1981-2010     TRUE      Year                 3
```

To visualize missing data we can use the `gg_miss_var()` function from the `naniar` package.

```r
select(normals, -contains("_code")) %>%  # Remove '_code' columns
  gg_miss_var(facet = station_name)
```


```r
suppressWarnings({select(normals, -contains("_code")) %>%  # Remove '_code' columns
    gg_miss_var(facet = station_name)})
```

<img src="unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" width="100%" style="display: block; margin: auto;" />

Let's take a look at the frost data.


```r
if("normals" %in% names(frost)) frost <- select(frost, -normals) # tidyr v1
glimpse(frost)
```

```
## Rows: 7
## Columns: 13
## $ prov                                  <chr> "MB", "MB", "MB", "MB", "MB", "MB", "MB"
## $ station_name                          <chr> "WINNIPEG RICHARDSON INT'L A", "WINNIPEG R…
## $ climate_id                            <chr> "5023222", "5023222", "5023222", "5023222"…
## $ normals_years                         <chr> "1981-2010", "1981-2010", "1981-2010", "19…
## $ meets_wmo                             <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE
## $ frost_code                            <chr> "A", "A", "A", "A", "A", "A", "A"
## $ date_first_fall_frost                 <dbl> 265, 265, 265, 265, 265, 265, 265
## $ date_last_spring_frost                <dbl> 143, 143, 143, 143, 143, 143, 143
## $ length_frost_free                     <dbl> 121, 121, 121, 121, 121, 121, 121
## $ prob                                  <chr> "10%", "25%", "33%", "50%", "66%", "75%", …
## $ prob_first_fall_temp_below_0_on_date  <dbl> 255, 259, 261, 265, 268, 270, 276
## $ prob_length_frost_free                <dbl> 96, 109, 114, 119, 126, 129, 141
## $ prob_last_spring_temp_below_0_on_date <dbl> 158, 152, 148, 144, 140, 137, 129
```

### Finding stations with specific measurements

The include data frame, `normals_measurements` contains a list of stations with their corresponding measurements. Be aware that this data might be out of date!


```r
normals_measurements
```

```
## # A tibble: 307,891 × 5
##    prov  station_name climate_id normals   measurement            
##    <chr> <chr>        <chr>      <chr>     <chr>                  
##  1 AB    HORBURG      301C3D4    1981-2010 temp_daily_average     
##  2 AB    HORBURG      301C3D4    1981-2010 temp_daily_average_code
##  3 AB    HORBURG      301C3D4    1981-2010 temp_sd                
##  4 AB    HORBURG      301C3D4    1981-2010 temp_sd_code           
##  5 AB    HORBURG      301C3D4    1981-2010 temp_daily_max         
##  6 AB    HORBURG      301C3D4    1981-2010 temp_daily_max_code    
##  7 AB    HORBURG      301C3D4    1981-2010 temp_daily_min         
##  8 AB    HORBURG      301C3D4    1981-2010 temp_daily_min_code    
##  9 AB    HORBURG      301C3D4    1981-2010 temp_extreme_max       
## 10 AB    HORBURG      301C3D4    1981-2010 temp_extreme_max_code  
## # … with 307,881 more rows
```

For example, if you wanted all `climate_id`s for stations that have data on
soil temperature for 1981-2010 normals:


```r
normals_measurements %>%
  filter(stringr::str_detect(measurement, "soil"),
         normals == "1981-2010") %>%
  pull(climate_id) %>%
  unique()
```

```
##  [1] "3070560" "1100119" "112G8L1" "5021054" "5021848" "8102234" "8403600" "8501900"
##  [9] "8502800" "8202800" "8205990" "2403500" "6073960" "6104025" "6105976" "7040440"
## [17] "7042388" "4012400" "4019035" "4028060" "4043900" "4075518"
```

## Understanding Climate Normals

The measurements contained in the climate normals are very specific. To better understand how they are calculated please explore the following resources:

- ECCC Climate Normals Calculations ([1981-2010](https://climate.weather.gc.ca/doc/Canadian_Climate_Normals_1981_2010_Calculation_Information.pdf) | [1971-2000](https://climate.weather.gc.ca/doc/Canadian_Climate_Normals_1971_2000_Calculation_Information.pdf))
    - [`weathercan` Climate Normals Codes](flags.html)
- [ECCC Climate Normals Technical Documentation](https://www.canada.ca/en/environment-climate-change/services/climate-change/canadian-centre-climate-services/display-download/technical-documentation-climate-normals.html)
    - [`weathercan` Climate Normals Glossary](glossary_normals.html)



---
title: "Weather: Terms and Units"
author: "Steffi LaZerte"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Weather: Terms and Units}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
library(weathercan)
library(dplyr)
```

This table shows details regarding original column (measurement) names and units of all weather measurements. It further provides links back to the ECCC glossary for more details.

For details on climate normals measurements, see the `glossary_normals` vignette.

```{r, asis = TRUE, echo = FALSE}
temp <- glossary %>%
  mutate(http = stringr::str_detect(ECCC_ref, "http"),
         ECCC_ref = replace(ECCC_ref, http & !is.na(http), paste0("[ECCC glossary page](", ECCC_ref[http & !is.na(http)], ")")),
         ECCC_ref = replace(ECCC_ref, !http & !is.na(http), "[See the 'flags' vignette](flags.html)")) %>%
  select(Interval = interval, `ECCC Name` = ECCC_name, `Formatted weathercan name` = weathercan_name, units, Reference = ECCC_ref)
  
knitr::kable(temp)
```
---
title: "Getting Started"
author: "Steffi LaZerte"
date: "2021-11-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Getting Started}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---





```r
library(dplyr)
library(ggplot2)
library(weathercan)
```


## Stations

`weathercan` includes the function `stations()` which returns a list of stations and their details (including `station_id`).


```r
head(stations())
```

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
```

```r
glimpse(stations())
```

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
```

You can look through this data frame directly, or you can use the `stations_search` function:


```r
stations_search("Kamloops")
```

```
## # A tibble: 40 × 16
##    prov  station_name         station_id climate_id WMO_id TC_id   lat   lon  elev tz        interval start   end normals normals_1981_2010 normals_1971_2000
##    <chr> <chr>                     <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <chr>    <dbl> <dbl> <lgl>   <lgl>             <lgl>            
##  1 BC    KAMLOOPS                   1274 1163779        NA <NA>   50.7 -120.  379. Etc/GMT+8 day       1878  1982 FALSE   FALSE             FALSE            
##  2 BC    KAMLOOPS                   1274 1163779        NA <NA>   50.7 -120.  379. Etc/GMT+8 month     1878  1982 FALSE   FALSE             FALSE            
##  3 BC    KAMLOOPS A                 1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 day       1951  2013 TRUE    TRUE              TRUE             
##  4 BC    KAMLOOPS A                 1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      1953  2013 TRUE    TRUE              TRUE             
##  5 BC    KAMLOOPS A                 1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 month     1951  2013 TRUE    TRUE              TRUE             
##  6 BC    KAMLOOPS A                51423 1163781     71887 YKA    50.7 -120.  345. Etc/GMT+8 day       2013  2021 FALSE   FALSE             FALSE            
##  7 BC    KAMLOOPS A                51423 1163781     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      2013  2021 FALSE   FALSE             FALSE            
##  8 BC    KAMLOOPS AFTON MINES       1276 1163790        NA <NA>   50.7 -120.  701  Etc/GMT+8 day       1977  1993 FALSE   FALSE             TRUE             
##  9 BC    KAMLOOPS AFTON MINES       1276 1163790        NA <NA>   50.7 -120.  701  Etc/GMT+8 month     1977  1993 FALSE   FALSE             TRUE             
## 10 BC    KAMLOOPS AUT              42203 1163842     71741 ZKA    50.7 -120.  345  Etc/GMT+8 day       2006  2021 FALSE   FALSE             FALSE            
## # … with 30 more rows
```

You can narrow down your search by specifying time intervals (options are "hour", "day", or "month"):


```r
stations_search("Kamloops", interval = "hour")
```

```
## # A tibble: 3 × 16
##   prov  station_name station_id climate_id WMO_id TC_id   lat   lon  elev tz        interval start   end normals normals_1981_2010 normals_1971_2000
##   <chr> <chr>             <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <chr>    <dbl> <dbl> <lgl>   <lgl>             <lgl>            
## 1 BC    KAMLOOPS A         1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      1953  2013 TRUE    TRUE              TRUE             
## 2 BC    KAMLOOPS A        51423 1163781     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      2013  2021 FALSE   FALSE             FALSE            
## 3 BC    KAMLOOPS AUT      42203 1163842     71741 ZKA    50.7 -120.  345  Etc/GMT+8 hour      2006  2021 FALSE   FALSE             FALSE
```

You can specify more than one interval:


```r
stations_search("Kamloops", interval = c("hour", "month"))
```

```
## # A tibble: 21 × 16
##    prov  station_name            station_id climate_id WMO_id TC_id   lat   lon  elev tz        interval start   end normals normals_1981_2010 normals_1971_2000
##    <chr> <chr>                        <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <chr>    <dbl> <dbl> <lgl>   <lgl>             <lgl>            
##  1 BC    KAMLOOPS                      1274 1163779        NA <NA>   50.7 -120.  379. Etc/GMT+8 month     1878  1982 FALSE   FALSE             FALSE            
##  2 BC    KAMLOOPS A                    1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      1953  2013 TRUE    TRUE              TRUE             
##  3 BC    KAMLOOPS A                    1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 month     1951  2013 TRUE    TRUE              TRUE             
##  4 BC    KAMLOOPS A                   51423 1163781     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      2013  2021 FALSE   FALSE             FALSE            
##  5 BC    KAMLOOPS AFTON MINES          1276 1163790        NA <NA>   50.7 -120.  701  Etc/GMT+8 month     1977  1993 FALSE   FALSE             TRUE             
##  6 BC    KAMLOOPS AUT                 42203 1163842     71741 ZKA    50.7 -120.  345  Etc/GMT+8 hour      2006  2021 FALSE   FALSE             FALSE            
##  7 BC    KAMLOOPS AUT                 42203 1163842     71741 ZKA    50.7 -120.  345  Etc/GMT+8 month     2006  2006 FALSE   FALSE             FALSE            
##  8 BC    KAMLOOPS CDA                  1277 1163810        NA <NA>   50.7 -120.  345  Etc/GMT+8 month     1949  1977 FALSE   FALSE             FALSE            
##  9 BC    KAMLOOPS CHERRY CREEK         1278 1163814        NA <NA>   50.7 -121.  556. Etc/GMT+8 month     1970  1974 FALSE   FALSE             FALSE            
## 10 BC    KAMLOOPS CHERRY CREEK 2       1279 1163815        NA <NA>   50.6 -121.  701  Etc/GMT+8 month     1974  1977 FALSE   FALSE             FALSE            
## # … with 11 more rows
```


You can also search by proximity. These results include a new column `distance` specifying the distance in km from the coordinates:


```r
stations_search(coords = c(50.667492, -120.329049), dist = 20, interval = "hour")
```

```
## # A tibble: 3 × 17
##   prov  station_name station_id climate_id WMO_id TC_id   lat   lon  elev tz        interval start   end normals normals_1981_2010 normals_1971_2000 distance
##   <chr> <chr>             <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <chr>    <dbl> <dbl> <lgl>   <lgl>             <lgl>                <dbl>
## 1 BC    KAMLOOPS A         1275 1163780     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      1953  2013 TRUE    TRUE              TRUE                  8.64
## 2 BC    KAMLOOPS AUT      42203 1163842     71741 ZKA    50.7 -120.  345  Etc/GMT+8 hour      2006  2021 FALSE   FALSE             FALSE                 8.64
## 3 BC    KAMLOOPS A        51423 1163781     71887 YKA    50.7 -120.  345. Etc/GMT+8 hour      2013  2021 FALSE   FALSE             FALSE                 9.28
```

We can also perform more complex searches using `filter()` function from the `dplyr` package
direction on the data returned by stations():


```r
BCstations <- stations() %>%
  filter(prov %in% c("BC")) %>%
  filter(interval == "hour") %>%
  filter(lat > 49 & lat < 49.5) %>%
  filter(lon > -119 & lon < -116) %>%
  filter(start <= 2002) %>%
  filter(end >= 2016)
BCstations
```

```
## # A tibble: 3 × 16
##   prov  station_name                station_id climate_id WMO_id TC_id   lat   lon  elev tz      interval start   end normals normals_1981_2010 normals_1971_20…
##   <chr> <chr>                            <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>   <chr>    <dbl> <dbl> <lgl>   <lgl>             <lgl>           
## 1 BC    CRESTON CAMPBELL SCIENTIFIC       6838 114B1F0     71770 WJR    49.1 -116.  641. Etc/GM… hour      1994  2021 FALSE   FALSE             FALSE           
## 2 BC    NELSON CS                         6839 1145M29     71776 WNM    49.5 -117.  535. Etc/GM… hour      1994  2021 FALSE   FALSE             FALSE           
## 3 BC    WARFIELD RCS                     31067 1148705     71401 XWF    49.1 -118.  567. Etc/GM… hour      2001  2021 FALSE   FALSE             FALSE
```

```r
## weather_dl() accepts numbers so we can create a vector to input into weather:
stn_vector <- BCstations$station_id
stn_vector
```

```
## [1]  6838  6839 31067
```

You can update this list of stations with


```r
stations_dl()
```

And check when it was last updated with

```r
stations_meta()
```

```
## $ECCC_modified
## [1] "2021-10-31 23:34:00 UTC"
## 
## $weathercan_modified
## [1] "2021-11-30"
```


## Weather

Once you have your `station_id`(s) you can download weather data:


```r
kam <- weather_dl(station_ids = 51423, start = "2016-01-01", end = "2016-02-15")
```

```
## As of weathercan v0.3.0 time display is either local time or UTC
## See Details under ?weather_dl for more information.
## This message is shown once per session
```

```r
kam
```

```
## # A tibble: 1,104 × 37
##    station_name station_id station_operator prov    lat   lon  elev climate_id WMO_id TC_id date       time                year  month day   hour  weather  hmdx
##    <chr>             <dbl> <lgl>            <chr> <dbl> <dbl> <dbl> <chr>      <chr>  <chr> <date>     <dttm>              <chr> <chr> <chr> <chr> <chr>   <dbl>
##  1 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 00:00:00 2016  01    01    00:00 <NA>       NA
##  2 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 01:00:00 2016  01    01    01:00 Mostly…    NA
##  3 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 02:00:00 2016  01    01    02:00 <NA>       NA
##  4 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 03:00:00 2016  01    01    03:00 <NA>       NA
##  5 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 04:00:00 2016  01    01    04:00 Cloudy     NA
##  6 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 05:00:00 2016  01    01    05:00 <NA>       NA
##  7 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 06:00:00 2016  01    01    06:00 <NA>       NA
##  8 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 07:00:00 2016  01    01    07:00 Cloudy     NA
##  9 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 08:00:00 2016  01    01    08:00 <NA>       NA
## 10 KAMLOOPS A        51423 NA               BC     50.7 -120.  345. 1163781    71887  YKA   2016-01-01 2016-01-01 09:00:00 2016  01    01    09:00 Snow       NA
## # … with 1,094 more rows, and 19 more variables: hmdx_flag <chr>, precip_amt <dbl>, precip_amt_flag <chr>, pressure <dbl>, pressure_flag <chr>, rel_hum <dbl>,
## #   rel_hum_flag <chr>, temp <dbl>, temp_dew <dbl>, temp_dew_flag <chr>, temp_flag <chr>, visib <dbl>, visib_flag <chr>, wind_chill <dbl>,
## #   wind_chill_flag <chr>, wind_dir <dbl>, wind_dir_flag <chr>, wind_spd <dbl>, wind_spd_flag <chr>
```

You can also download data from multiple stations at once:


```r
kam.pg <- weather_dl(station_ids = c(48248, 51423), start = "2016-01-01", end = "2016-02-15")

kam.pg
```

```
## # A tibble: 2,208 × 37
##    station_name station_id station_operator prov    lat   lon  elev climate_id WMO_id TC_id date       time                year  month day   hour  weather  hmdx
##    <chr>             <dbl> <lgl>            <chr> <dbl> <dbl> <dbl> <chr>      <chr>  <chr> <date>     <dttm>              <chr> <chr> <chr> <chr> <chr>   <dbl>
##  1 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 00:00:00 2016  01    01    00:00 <NA>       NA
##  2 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 01:00:00 2016  01    01    01:00 <NA>       NA
##  3 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 02:00:00 2016  01    01    02:00 <NA>       NA
##  4 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 03:00:00 2016  01    01    03:00 <NA>       NA
##  5 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 04:00:00 2016  01    01    04:00 <NA>       NA
##  6 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 05:00:00 2016  01    01    05:00 <NA>       NA
##  7 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 06:00:00 2016  01    01    06:00 <NA>       NA
##  8 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 07:00:00 2016  01    01    07:00 <NA>       NA
##  9 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 08:00:00 2016  01    01    08:00 <NA>       NA
## 10 PRINCE GEOR…      48248 NA               BC     53.9 -123.   680 1096453    71302  VXS   2016-01-01 2016-01-01 09:00:00 2016  01    01    09:00 <NA>       NA
## # … with 2,198 more rows, and 19 more variables: hmdx_flag <chr>, precip_amt <dbl>, precip_amt_flag <chr>, pressure <dbl>, pressure_flag <chr>, rel_hum <dbl>,
## #   rel_hum_flag <chr>, temp <dbl>, temp_dew <dbl>, temp_dew_flag <chr>, temp_flag <chr>, visib <dbl>, visib_flag <chr>, wind_chill <dbl>,
## #   wind_chill_flag <chr>, wind_dir <dbl>, wind_dir_flag <chr>, wind_spd <dbl>, wind_spd_flag <chr>
```


And plot it:


```r
ggplot(data = kam.pg, aes(x = time, y = temp, group = station_name, colour = station_name)) +
  theme(legend.position = "top") +
  geom_line() +
  theme_minimal()
```

<img src="unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" style="display: block; margin: auto;" />

Or you can use the vector created above:


```r
stn_vec_df <- weather_dl(station_ids = stn_vector, start = "2016-01-01", end = "2016-02-15")

stn_vec_df
```

```
## # A tibble: 3,312 × 37
##    station_name station_id station_operator prov    lat   lon  elev climate_id WMO_id TC_id date       time                year  month day   hour  weather  hmdx
##    <chr>             <dbl> <lgl>            <chr> <dbl> <dbl> <dbl> <chr>      <chr>  <chr> <date>     <dttm>              <chr> <chr> <chr> <chr> <chr>   <dbl>
##  1 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 00:00:00 2016  01    01    00:00 <NA>       NA
##  2 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 01:00:00 2016  01    01    01:00 <NA>       NA
##  3 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 02:00:00 2016  01    01    02:00 <NA>       NA
##  4 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 03:00:00 2016  01    01    03:00 <NA>       NA
##  5 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 04:00:00 2016  01    01    04:00 <NA>       NA
##  6 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 05:00:00 2016  01    01    05:00 <NA>       NA
##  7 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 06:00:00 2016  01    01    06:00 <NA>       NA
##  8 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 07:00:00 2016  01    01    07:00 <NA>       NA
##  9 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 08:00:00 2016  01    01    08:00 <NA>       NA
## 10 CRESTON CAM…       6838 NA               BC     49.1 -116.  641. 114B1F0    71770  WJR   2016-01-01 2016-01-01 09:00:00 2016  01    01    09:00 <NA>       NA
## # … with 3,302 more rows, and 19 more variables: hmdx_flag <chr>, precip_amt <dbl>, precip_amt_flag <chr>, pressure <dbl>, pressure_flag <chr>, rel_hum <dbl>,
## #   rel_hum_flag <chr>, temp <dbl>, temp_dew <dbl>, temp_dew_flag <chr>, temp_flag <chr>, visib <dbl>, visib_flag <chr>, wind_chill <dbl>,
## #   wind_chill_flag <chr>, wind_dir <dbl>, wind_dir_flag <chr>, wind_spd <dbl>, wind_spd_flag <chr>
```

For more information on the data flags, see the [Flags vignette](flags.html), for more information on units and terms, see the [Terms and Units vignette](glossary.html).

## Climate Normals

To access climate normals, you first need to know the `climate_id` associated with the station you're interested in.


```r
stations_search("Winnipeg", normals_years = "current")
```

```
## # A tibble: 1 × 13
##   prov  station_name                station_id climate_id WMO_id TC_id   lat   lon  elev tz        normals normals_1981_2010 normals_1971_2000
##   <chr> <chr>                            <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <lgl>   <lgl>             <lgl>            
## 1 MB    WINNIPEG RICHARDSON INT'L A       3698 5023222     71852 YWG    49.9 -97.2  239. Etc/GMT+6 TRUE    TRUE              TRUE
```

The current year range is 1981-2010, but you can also search for stations in the
previous year range:


```r
stations_search("Winnipeg", normals_years = "1971-2000")
```

```
## # A tibble: 1 × 13
##   prov  station_name                station_id climate_id WMO_id TC_id   lat   lon  elev tz        normals normals_1981_2010 normals_1971_2000
##   <chr> <chr>                            <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <lgl>   <lgl>             <lgl>            
## 1 MB    WINNIPEG RICHARDSON INT'L A       3698 5023222     71852 YWG    49.9 -97.2  239. Etc/GMT+6 TRUE    TRUE              TRUE
```

Note that the Winnipeg station has normals for both year ranges.

Then you can download the climate normals with the `normals_dl()` function.


```r
n <- normals_dl("5023222")
```

There are two parts to the normals data, average weather measurements and average frost dates.


```r
library(tidyr)
unnest(n, normals)
```

```
## # A tibble: 13 × 203
##    prov  station_name       climate_id normals_years meets_wmo period temp_daily_avera… temp_daily_averag… temp_sd temp_sd_code temp_daily_max temp_daily_max_c…
##    <chr> <chr>              <chr>      <chr>         <lgl>     <fct>              <dbl> <chr>                <dbl> <chr>                 <dbl> <chr>            
##  1 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Jan                -16.4 A                      4.1 A                     -11.3 A                
##  2 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Feb                -13.2 A                      4.2 A                      -8.1 A                
##  3 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Mar                 -5.8 A                      3.1 A                      -0.8 A                
##  4 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Apr                  4.4 A                      2.7 A                      10.9 A                
##  5 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      May                 11.6 A                      2.1 A                      18.6 A                
##  6 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Jun                 17   A                      2   A                      23.2 A                
##  7 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Jul                 19.7 A                      1.4 A                      25.9 A                
##  8 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Aug                 18.8 A                      1.9 A                      25.4 A                
##  9 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Sep                 12.7 A                      1.3 A                      19   A                
## 10 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Oct                  5   A                      1.8 A                      10.5 A                
## 11 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Nov                 -4.9 A                      3.6 A                      -0.5 A                
## 12 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Dec                -13.2 A                      4.4 A                      -8.5 A                
## 13 MB    WINNIPEG RICHARDS… 5023222    1981-2010     TRUE      Year                 3   A                      1.2 A                       8.7 A                
## # … with 191 more variables: temp_daily_min <dbl>, temp_daily_min_code <chr>, temp_extreme_max <dbl>, temp_extreme_max_code <chr>,
## #   temp_extreme_max_date <date>, temp_extreme_max_date_code <chr>, temp_extreme_min <dbl>, temp_extreme_min_code <chr>, temp_extreme_min_date <date>,
## #   temp_extreme_min_date_code <chr>, rain <dbl>, rain_code <chr>, snow <dbl>, snow_code <chr>, precip <dbl>, precip_code <chr>, snow_mean_depth <dbl>,
## #   snow_mean_depth_code <chr>, snow_median_depth <dbl>, snow_median_depth_code <chr>, snow_depth_month_end <dbl>, snow_depth_month_end_code <chr>,
## #   rain_extreme_daily <dbl>, rain_extreme_daily_code <chr>, rain_extreme_daily_date <date>, rain_extreme_daily_date_code <chr>, snow_extreme_daily <dbl>,
## #   snow_extreme_daily_code <chr>, snow_extreme_daily_date <date>, snow_extreme_daily_date_code <chr>, precip_extreme_daily <dbl>,
## #   precip_extreme_daily_code <chr>, precip_extreme_daily_date <date>, precip_extreme_daily_date_code <chr>, snow_extreme_depth <dbl>, …
```

```r
unnest(n, frost)
```

```
## # A tibble: 7 × 14
##   prov  station_name     climate_id normals_years meets_wmo normals   frost_code date_first_fall_… date_last_spring… length_frost_fr… prob  prob_first_fall_tem…
##   <chr> <chr>            <chr>      <chr>         <lgl>     <list>    <chr>                  <dbl>             <dbl>            <dbl> <chr>                <dbl>
## 1 MB    WINNIPEG RICHAR… 5023222    1981-2010     TRUE      <tibble … A                        265               143              121 10%                    255
## 2 MB    WINNIPEG RICHAR… 5023222    1981-2010     TRUE      <tibble … A                        265               143              121 25%                    259
## 3 MB    WINNIPEG RICHAR… 5023222    1981-2010     TRUE      <tibble … A                        265               143              121 33%                    261
## 4 MB    WINNIPEG RICHAR… 5023222    1981-2010     TRUE      <tibble … A                        265               143              121 50%                    265
## 5 MB    WINNIPEG RICHAR… 5023222    1981-2010     TRUE      <tibble … A                        265               143              121 66%                    268
## 6 MB    WINNIPEG RICHAR… 5023222    1981-2010     TRUE      <tibble … A                        265               143              121 75%                    270
## 7 MB    WINNIPEG RICHAR… 5023222    1981-2010     TRUE      <tibble … A                        265               143              121 90%                    276
## # … with 2 more variables: prob_length_frost_free <dbl>, prob_last_spring_temp_below_0_on_date <dbl>
```


Alternatively, download the 1971-2000 normals:


```r
n <- normals_dl("5023222", normals_years = "1971-2000")
unnest(n, normals)
```

```
## # A tibble: 13 × 229
##    prov  station_name       climate_id normals_years meets_wmo period temp_daily_avera… temp_daily_averag… temp_sd temp_sd_code temp_daily_max temp_daily_max_c…
##    <chr> <chr>              <chr>      <chr>         <lgl>     <fct>              <dbl> <chr>                <dbl> <chr>                 <dbl> <chr>            
##  1 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Jan                -17.8 A                      3.9 A                     -12.7 A                
##  2 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Feb                -13.6 A                      4.2 A                      -8.5 A                
##  3 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Mar                 -6.1 A                      3.5 A                      -1.1 A                
##  4 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Apr                  4   A                      2.7 A                      10.3 A                
##  5 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      May                 12   A                      2.5 A                      19.2 A                
##  6 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Jun                 17   A                      1.8 A                      23.3 A                
##  7 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Jul                 19.5 A                      1.5 A                      25.8 A                
##  8 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Aug                 18.5 A                      1.8 A                      25   A                
##  9 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Sep                 12.3 A                      1.4 A                      18.6 A                
## 10 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Oct                  5.3 A                      1.6 A                      10.8 A                
## 11 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Nov                 -5.3 A                      3.3 A                      -0.9 A                
## 12 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Dec                -14.4 A                      4.2 A                      -9.7 A                
## 13 MB    WINNIPEG RICHARDS… 5023222    1971-2000     TRUE      Year                 2.6 A                      1.3 A                       8.3 A                
## # … with 217 more variables: temp_daily_min <dbl>, temp_daily_min_code <chr>, temp_extreme_max <dbl>, temp_extreme_max_code <chr>,
## #   temp_extreme_max_date <date>, temp_extreme_max_date_code <chr>, temp_extreme_min <dbl>, temp_extreme_min_code <chr>, temp_extreme_min_date <date>,
## #   temp_extreme_min_date_code <chr>, rain <dbl>, rain_code <chr>, snow <dbl>, snow_code <chr>, precip <dbl>, precip_code <chr>, snow_mean_depth <dbl>,
## #   snow_mean_depth_code <chr>, snow_median_depth <dbl>, snow_median_depth_code <chr>, snow_depth_month_end <dbl>, snow_depth_month_end_code <chr>,
## #   rain_extreme_daily <dbl>, rain_extreme_daily_code <chr>, rain_extreme_daily_date <date>, rain_extreme_daily_date_code <chr>, snow_extreme_daily <dbl>,
## #   snow_extreme_daily_code <chr>, snow_extreme_daily_date <date>, snow_extreme_daily_date_code <chr>, precip_extreme_daily <dbl>,
## #   precip_extreme_daily_code <chr>, precip_extreme_daily_date <date>, precip_extreme_daily_date_code <chr>, snow_extreme_depth <dbl>, …
```

```r
unnest(n, frost)
```

```
## # A tibble: 0 × 6
## # … with 6 variables: prov <chr>, station_name <chr>, climate_id <chr>, normals_years <chr>, meets_wmo <lgl>, normals <list>
```






---
title: "Interpolating"
author: "Steffi LaZerte"
date: "2021-11-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Interpolating}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




## Packages

You'll need several packages from the **tidyverse** in addition to **`weathercan`** to complete the following analysis.


```r
library(weathercan)
library(ggplot2)
library(dplyr)
```

## General usage
You can merge weather data with other data frames by linearly interpolating between points.

For example, here we have a dataset of weather data from Kamloops


```r
glimpse(kamloops)
```

```
## Rows: 4,368
## Columns: 37
## $ station_name     <chr> "KAMLOOPS A", "KAMLOOPS A", "KAMLOOPS A", "KAMLOOPS A", "KAMLOOPS A", "KAMLOOPS A", "KAMLO…
## $ station_id       <dbl> 51423, 51423, 51423, 51423, 51423, 51423, 51423, 51423, 51423, 51423, 51423, 51423, 51423,…
## $ station_operator <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ prov             <chr> "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", "BC", …
## $ lat              <dbl> 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, 50.7, …
## $ lon              <dbl> -120.45, -120.45, -120.45, -120.45, -120.45, -120.45, -120.45, -120.45, -120.45, -120.45, …
## $ elev             <dbl> 345.3, 345.3, 345.3, 345.3, 345.3, 345.3, 345.3, 345.3, 345.3, 345.3, 345.3, 345.3, 345.3,…
## $ climate_id       <chr> "1163781", "1163781", "1163781", "1163781", "1163781", "1163781", "1163781", "1163781", "1…
## $ WMO_id           <chr> "71887", "71887", "71887", "71887", "71887", "71887", "71887", "71887", "71887", "71887", …
## $ TC_id            <chr> "YKA", "YKA", "YKA", "YKA", "YKA", "YKA", "YKA", "YKA", "YKA", "YKA", "YKA", "YKA", "YKA",…
## $ date             <date> 2016-01-01, 2016-01-01, 2016-01-01, 2016-01-01, 2016-01-01, 2016-01-01, 2016-01-01, 2016-…
## $ time             <dttm> 2016-01-01 00:00:00, 2016-01-01 01:00:00, 2016-01-01 02:00:00, 2016-01-01 03:00:00, 2016-…
## $ year             <chr> "2016", "2016", "2016", "2016", "2016", "2016", "2016", "2016", "2016", "2016", "2016", "2…
## $ month            <chr> "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", …
## $ day              <chr> "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", "01", …
## $ hour             <chr> "00:00", "01:00", "02:00", "03:00", "04:00", "05:00", "06:00", "07:00", "08:00", "09:00", …
## $ weather          <chr> NA, "Mostly Cloudy", NA, NA, "Cloudy", NA, NA, "Cloudy", NA, "Snow", "Snow", "Snow", "Snow…
## $ hmdx             <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ hmdx_flag        <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ precip_amt       <dbl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ precip_amt_flag  <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ pressure         <dbl> 99.95, 99.93, 99.92, 99.90, 99.86, 99.82, 99.80, 99.78, 99.77, 99.78, 99.79, 99.74, 99.69,…
## $ pressure_flag    <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ rel_hum          <dbl> 74, 76, 74, 73, 70, 71, 69, 69, 71, 71, 71, 70, 69, 70, 68, 68, 70, 74, 73, 74, 74, 74, 77…
## $ rel_hum_flag     <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ temp             <dbl> -9.1, -9.6, -9.9, -9.5, -9.4, -9.8, -10.0, -10.2, -10.1, -9.7, -9.4, -9.0, -8.6, -8.2, -8.…
## $ temp_dew         <dbl> -12.9, -13.1, -13.7, -13.5, -13.9, -14.1, -14.7, -14.9, -14.4, -14.0, -13.7, -13.5, -13.3,…
## $ temp_dew_flag    <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ temp_flag        <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ visib            <dbl> 64.4, 64.4, 64.4, 64.4, 64.4, 64.4, 64.4, 64.4, 48.3, 48.3, 48.3, 48.3, 48.3, 48.3, 48.3, …
## $ visib_flag       <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ wind_chill       <dbl> -17, -17, -18, -17, -17, -17, -18, -17, -17, -16, -15, -14, -14, -13, -13, -13, -13, -14, …
## $ wind_chill_flag  <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ wind_dir         <dbl> 13, 11, 11, 11, 11, 10, 9, 7, 7, 10, 11, 10, 10, 13, 11, 10, 10, 9, 12, 10, 13, 12, 10, 12…
## $ wind_dir_flag    <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
## $ wind_spd         <dbl> 19, 20, 20, 18, 18, 16, 23, 15, 14, 15, 12, 11, 12, 9, 10, 12, 11, 12, 10, 11, 11, 6, 6, 4…
## $ wind_spd_flag    <chr> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
```

As well as a data set of finch visits to an RFID feeder

```r
glimpse(finches)
```

```
## Rows: 16,886
## Columns: 10
## $ animal_id <fct> 041868FF93, 041868FF93, 041868FF93, 06200003BB, 06200003BB, 06200003BB, 06200003BB, 06200003BB, 0…
## $ date      <date> 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, …
## $ time      <dttm> 2016-03-01 06:57:42, 2016-03-01 06:58:41, 2016-03-01 07:07:21, 2016-03-01 07:32:34, 2016-03-01 0…
## $ logger_id <fct> 2300, 2300, 2300, 2400, 2400, 2400, 2400, 2400, 2300, 2300, 2300, 2300, 2300, 2400, 2300, 2400, 2…
## $ species   <chr> "Mountain Chickadee", "Mountain Chickadee", "Mountain Chickadee", "House Finch", "House Finch", "…
## $ age       <chr> "AHY", "AHY", "AHY", "SY", "SY", "SY", "SY", "SY", "AHY", "AHY", "AHY", "AHY", "AHY", "SY", "AHY"…
## $ sex       <chr> "U", "U", "U", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "M", "F", "M", "M", "M", "M", "M…
## $ site_name <chr> "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "…
## $ lon       <dbl> -120.3622, -120.3622, -120.3622, -120.3635, -120.3635, -120.3635, -120.3635, -120.3635, -120.3622…
## $ lat       <dbl> 50.66967, 50.66967, 50.66967, 50.66938, 50.66938, 50.66938, 50.66938, 50.66938, 50.66967, 50.6696…
```

Although the times in the weather data do not exactly match those in the finch data, we can merge them together through linear [interpolation](https://en.wikipedia.org/wiki/Linear_interpolation). This function uses the `approx` function from the `stats` package under the hood.

Here we specify that we only want the temperature (`temp`) column:


```r
finches_temperature <- weather_interp(data = finches, weather = kamloops, cols = "temp")
```

```
## temp is missing 4 out of 4368 data, interpolation may be less accurate as a result.
```

```r
summary(finches_temperature)
```

```
##       animal_id         date                 time                     logger_id     species         
##  0620000513:7624   Min.   :2016-03-01   Min.   :2016-03-01 06:57:42   1500:6370   Length:16886      
##  041868D861:2767   1st Qu.:2016-03-05   1st Qu.:2016-03-05 13:54:13   2100: 968   Class :character  
##  0620000514:1844   Median :2016-03-09   Median :2016-03-09 16:54:47   2200:2266   Mode  :character  
##  06200004F8:1386   Mean   :2016-03-08   Mean   :2016-03-09 07:45:58   2300:3531                     
##  041868BED6: 944   3rd Qu.:2016-03-13   3rd Qu.:2016-03-13 08:24:58   2400:1477                     
##  06200003BB: 708   Max.   :2016-03-16   Max.   :2016-03-16 16:39:30   2700:2274                     
##  (Other)   :1613                                                                                    
##      age                sex             site_name              lon              lat             temp        
##  Length:16886       Length:16886       Length:16886       Min.   :-120.4   Min.   :50.67   Min.   :-0.2317  
##  Class :character   Class :character   Class :character   1st Qu.:-120.4   1st Qu.:50.67   1st Qu.: 5.0561  
##  Mode  :character   Mode  :character   Mode  :character   Median :-120.4   Median :50.67   Median : 7.1651  
##                                                           Mean   :-120.4   Mean   :50.67   Mean   : 7.4349  
##                                                           3rd Qu.:-120.4   3rd Qu.:50.67   3rd Qu.: 9.3319  
##                                                           Max.   :-120.4   Max.   :50.67   Max.   :16.3712  
## 
```

```r
glimpse(finches_temperature)
```

```
## Rows: 16,886
## Columns: 11
## $ animal_id <fct> 041868FF93, 041868FF93, 041868FF93, 06200003BB, 06200003BB, 06200003BB, 06200003BB, 06200003BB, 0…
## $ date      <date> 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, …
## $ time      <dttm> 2016-03-01 06:57:42, 2016-03-01 06:58:41, 2016-03-01 07:07:21, 2016-03-01 07:32:34, 2016-03-01 0…
## $ logger_id <fct> 2300, 2300, 2300, 2400, 2400, 2400, 2400, 2400, 2300, 2300, 2300, 2300, 2300, 2400, 2300, 2400, 2…
## $ species   <chr> "Mountain Chickadee", "Mountain Chickadee", "Mountain Chickadee", "House Finch", "House Finch", "…
## $ age       <chr> "AHY", "AHY", "AHY", "SY", "SY", "SY", "SY", "SY", "AHY", "AHY", "AHY", "AHY", "AHY", "SY", "AHY"…
## $ sex       <chr> "U", "U", "U", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "M", "F", "M", "M", "M", "M", "M…
## $ site_name <chr> "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "…
## $ lon       <dbl> -120.3622, -120.3622, -120.3622, -120.3635, -120.3635, -120.3635, -120.3635, -120.3635, -120.3622…
## $ lat       <dbl> 50.66967, 50.66967, 50.66967, 50.66938, 50.66938, 50.66938, 50.66938, 50.66938, 50.66967, 50.6696…
## $ temp      <dbl> 3.984667, 3.991222, 4.036750, 4.162833, 4.162917, 4.163000, 4.163083, 4.163167, 4.180417, 4.18075…
```

```r
ggplot(data = finches_temperature, aes(x = temp, fill = animal_id)) +
  theme_bw() +
  theme(legend.position = "none") +
  geom_histogram(binwidth = 1) +
  labs(x = "Temperature (C)", y = "Activity Count", fill = "Finch ID")
```

<img src="interp-unnamed-chunk-4-1.png" title="plot of chunk unnamed-chunk-4" alt="plot of chunk unnamed-chunk-4" width="100%" style="display: block; margin: auto;" />

Or summarized:

```r
finches_temperature <- finches_temperature %>%
  group_by(date) %>%
  summarize(n = length(time),
            temp = mean(temp))

ggplot(data = finches_temperature, aes(x = date, y = n)) +
  theme_bw() +
  theme(legend.position = "top") +
  geom_point(aes(shape = "Activity")) +
  geom_line(aes(y = temp * 100, colour = "Temperature")) +
  scale_colour_discrete(name = "") +
  scale_shape_discrete(name = "") +
  scale_y_continuous(name = "Activity", sec.axis = sec_axis(~. / 100, name = "Temperature (C)"))
```

<img src="interp-unnamed-chunk-5-1.png" title="plot of chunk unnamed-chunk-5" alt="plot of chunk unnamed-chunk-5" width="100%" style="display: block; margin: auto;" />

## Data gaps

By default, gaps of 2 hours (or 2 days, with a daily scale) will be interpolated over (i.e. they will be filled with values interpolated from either side of the gap), but longer gaps will be skipped and filled with `NA`s. You can adjust this behaviour with `na_gap`. Note that as Environment and Climate Change Canada data is downloaded on an hourly scale, it makes no sense to apply `na_gap` values of less than 1.

In this example, note the larger number of `NA`s in `temp` and how it corresponds to the missing variables in the weather dataset:


```r
finches_temperature <- weather_interp(data = finches, weather = kamloops,
                                      cols = "temp", na_gap = 1)
```

```
## temp is missing 4 out of 4368 data, interpolation may be less accurate as a result.
```

```r
summary(finches_temperature)
```

```
##       animal_id         date                 time                     logger_id     species         
##  0620000513:7624   Min.   :2016-03-01   Min.   :2016-03-01 06:57:42   1500:6370   Length:16886      
##  041868D861:2767   1st Qu.:2016-03-05   1st Qu.:2016-03-05 13:54:13   2100: 968   Class :character  
##  0620000514:1844   Median :2016-03-09   Median :2016-03-09 16:54:47   2200:2266   Mode  :character  
##  06200004F8:1386   Mean   :2016-03-08   Mean   :2016-03-09 07:45:58   2300:3531                     
##  041868BED6: 944   3rd Qu.:2016-03-13   3rd Qu.:2016-03-13 08:24:58   2400:1477                     
##  06200003BB: 708   Max.   :2016-03-16   Max.   :2016-03-16 16:39:30   2700:2274                     
##  (Other)   :1613                                                                                    
##      age                sex             site_name              lon              lat             temp        
##  Length:16886       Length:16886       Length:16886       Min.   :-120.4   Min.   :50.67   Min.   :-0.2317  
##  Class :character   Class :character   Class :character   1st Qu.:-120.4   1st Qu.:50.67   1st Qu.: 5.0746  
##  Mode  :character   Mode  :character   Mode  :character   Median :-120.4   Median :50.67   Median : 7.1668  
##                                                           Mean   :-120.4   Mean   :50.67   Mean   : 7.4433  
##                                                           3rd Qu.:-120.4   3rd Qu.:50.67   3rd Qu.: 9.3458  
##                                                           Max.   :-120.4   Max.   :50.67   Max.   :16.3712  
##                                                                                            NA's   :84
```

```r
finches_temperature %>%
  select(date, time, temp) %>%
  filter(is.na(temp))
```

```
## # A tibble: 84 × 3
##    date       time                 temp
##    <date>     <dttm>              <dbl>
##  1 2016-03-10 2016-03-10 16:00:12    NA
##  2 2016-03-10 2016-03-10 16:00:33    NA
##  3 2016-03-10 2016-03-10 16:00:36    NA
##  4 2016-03-10 2016-03-10 16:00:39    NA
##  5 2016-03-10 2016-03-10 16:00:42    NA
##  6 2016-03-10 2016-03-10 16:00:45    NA
##  7 2016-03-10 2016-03-10 16:00:48    NA
##  8 2016-03-10 2016-03-10 16:00:51    NA
##  9 2016-03-10 2016-03-10 16:00:54    NA
## 10 2016-03-10 2016-03-10 16:00:57    NA
## # … with 74 more rows
```

```r
kamloops %>%
  select(time, temp) %>%
  filter(is.na(temp))
```

```
## # A tibble: 4 × 2
##   time                 temp
##   <dttm>              <dbl>
## 1 2016-02-11 19:00:00    NA
## 2 2016-03-08 13:00:00    NA
## 3 2016-03-11 01:00:00    NA
## 4 2016-04-09 00:00:00    NA
```

## Multiple weather columns

We could also add in more than one column at a time:


```r
finches_weather <- weather_interp(data = finches, weather = kamloops,
                                  cols = c("temp", "wind_spd"))
```

```
## temp is missing 4 out of 4368 data, interpolation may be less accurate as a result.
```

```
## wind_spd is missing 4 out of 4368 data, interpolation may be less accurate as a result.
```

```r
summary(finches_weather)
```

```
##       animal_id         date                 time                     logger_id     species         
##  0620000513:7624   Min.   :2016-03-01   Min.   :2016-03-01 06:57:42   1500:6370   Length:16886      
##  041868D861:2767   1st Qu.:2016-03-05   1st Qu.:2016-03-05 13:54:13   2100: 968   Class :character  
##  0620000514:1844   Median :2016-03-09   Median :2016-03-09 16:54:47   2200:2266   Mode  :character  
##  06200004F8:1386   Mean   :2016-03-08   Mean   :2016-03-09 07:45:58   2300:3531                     
##  041868BED6: 944   3rd Qu.:2016-03-13   3rd Qu.:2016-03-13 08:24:58   2400:1477                     
##  06200003BB: 708   Max.   :2016-03-16   Max.   :2016-03-16 16:39:30   2700:2274                     
##  (Other)   :1613                                                                                    
##      age                sex             site_name              lon              lat             temp        
##  Length:16886       Length:16886       Length:16886       Min.   :-120.4   Min.   :50.67   Min.   :-0.2317  
##  Class :character   Class :character   Class :character   1st Qu.:-120.4   1st Qu.:50.67   1st Qu.: 5.0561  
##  Mode  :character   Mode  :character   Mode  :character   Median :-120.4   Median :50.67   Median : 7.1651  
##                                                           Mean   :-120.4   Mean   :50.67   Mean   : 7.4349  
##                                                           3rd Qu.:-120.4   3rd Qu.:50.67   3rd Qu.: 9.3319  
##                                                           Max.   :-120.4   Max.   :50.67   Max.   :16.3712  
##                                                                                                             
##     wind_spd     
##  Min.   : 1.000  
##  1st Qu.: 7.634  
##  Median :13.738  
##  Mean   :14.443  
##  3rd Qu.:19.907  
##  Max.   :44.939  
## 
```

```r
glimpse(finches_weather)
```

```
## Rows: 16,886
## Columns: 12
## $ animal_id <fct> 041868FF93, 041868FF93, 041868FF93, 06200003BB, 06200003BB, 06200003BB, 06200003BB, 06200003BB, 0…
## $ date      <date> 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, 2016-03-01, …
## $ time      <dttm> 2016-03-01 06:57:42, 2016-03-01 06:58:41, 2016-03-01 07:07:21, 2016-03-01 07:32:34, 2016-03-01 0…
## $ logger_id <fct> 2300, 2300, 2300, 2400, 2400, 2400, 2400, 2400, 2300, 2300, 2300, 2300, 2300, 2400, 2300, 2400, 2…
## $ species   <chr> "Mountain Chickadee", "Mountain Chickadee", "Mountain Chickadee", "House Finch", "House Finch", "…
## $ age       <chr> "AHY", "AHY", "AHY", "SY", "SY", "SY", "SY", "SY", "AHY", "AHY", "AHY", "AHY", "AHY", "SY", "AHY"…
## $ sex       <chr> "U", "U", "U", "M", "M", "M", "M", "M", "F", "F", "F", "F", "F", "M", "F", "M", "M", "M", "M", "M…
## $ site_name <chr> "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "Kamloops, BC", "…
## $ lon       <dbl> -120.3622, -120.3622, -120.3622, -120.3635, -120.3635, -120.3635, -120.3635, -120.3635, -120.3622…
## $ lat       <dbl> 50.66967, 50.66967, 50.66967, 50.66938, 50.66938, 50.66938, 50.66938, 50.66938, 50.66967, 50.6696…
## $ temp      <dbl> 3.984667, 3.991222, 4.036750, 4.162833, 4.162917, 4.163000, 4.163083, 4.163167, 4.180417, 4.18075…
## $ wind_spd  <dbl> 22.88500, 22.93417, 22.26500, 19.74333, 19.74167, 19.74000, 19.73833, 19.73667, 19.39167, 19.3850…
```

```r
finches_weather <- finches_weather %>%
  group_by(date) %>%
  summarize(n = length(time),
            temp = mean(temp),
            wind_spd = mean(wind_spd))

ggplot(data = finches_weather, aes(x = date, y = n)) +
  theme_bw() +
  theme(legend.position = "top") +
  geom_bar(stat = "identity") +
  geom_line(aes(y = temp * 50, colour = "Temperature"), size = 2) +
  geom_line(aes(y = wind_spd * 50, colour = "Wind Speed"), size = 2) +
  scale_colour_discrete(name = "") +
  scale_y_continuous(
    name = "Activity Counts",
    sec.axis = sec_axis(~. / 50, name = "Temperature (C) / Wind Speed (km/h)"))
```

<img src="interp-unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" width="100%" style="display: block; margin: auto;" />




---
title: "Flags and codes"
author: "Steffi LaZerte"
date: "2021-11-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Flags and codes}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




## What are flags/codes

The data output of the `weather_dl()` function include corresponding `_flag` columns for each data column. These columns are used by ECCC to add notes regarding measurements. 

Similarly, the data output of the `normals_dl()` function include corresponding `_code` columns. These columns are used by ECCC to add notes regarding the amount of data used to calculate the normals.


### Flags
In the `weather_dl()` function if `format = TRUE` (the default), data corresponding to flags `M`, `NA`, `[empty]` and `L` are all replaced with `NA`.

For example, a sample of unformatted data from Magog station in Quebec looks like:


```
## # A tibble: 12 × 6
##    station_name `Date/Time` `Total Precip (mm)` `Total Precip Fl… `Snow Grnd Last … `Snow Grnd Last…
##    <chr>        <chr>       <chr>               <chr>             <chr>             <chr>           
##  1 MAGOG        2017-03     30.4                ^                 <NA>              M               
##  2 MAGOG        2017-04     114.0               ^                 0                 <NA>            
##  3 MAGOG        2017-05     78.8                ^                 0                 <NA>            
##  4 MAGOG        2017-06     140.7               ^                 0                 <NA>            
##  5 MAGOG        2017-07     80.7                <NA>              0                 <NA>            
##  6 MAGOG        2017-08     135.8               <NA>              0                 <NA>            
##  7 MAGOG        2017-09     63.0                ^                 0                 <NA>            
##  8 MAGOG        2017-10     140.8               ^                 0                 <NA>            
##  9 MAGOG        2017-11     70.0                ^                 0                 <NA>            
## 10 MAGOG        2017-12     45.7                ^                 10                <NA>            
## 11 MAGOG        2018-01     34.6                ^                 2                 <NA>            
## 12 MAGOG        2018-02     77.2                ^                 0                 <NA>
```

In this output, you can see two flags: `^` in `Total Precip` and `M` in `Snow Grnd Last Day`

This same sample, formatted looks like:


```
## # A tibble: 12 × 5
##    date       total_precip total_precip_flag snow_grnd_last_day snow_grnd_last_day_flag
##    <date>            <dbl> <chr>                          <dbl> <chr>                  
##  1 2017-03-01         30.4 ^                                 NA M                      
##  2 2017-04-01        114   ^                                  0 <NA>                   
##  3 2017-05-01         78.8 ^                                  0 <NA>                   
##  4 2017-06-01        141.  ^                                  0 <NA>                   
##  5 2017-07-01         80.7 <NA>                               0 <NA>                   
##  6 2017-08-01        136.  <NA>                               0 <NA>                   
##  7 2017-09-01         63   ^                                  0 <NA>                   
##  8 2017-10-01        141.  ^                                  0 <NA>                   
##  9 2017-11-01         70   ^                                  0 <NA>                   
## 10 2017-12-01         45.7 ^                                 10 <NA>                   
## 11 2018-01-01         34.6 ^                                  2 <NA>                   
## 12 2018-02-01         77.2 ^                                  0 <NA>
```

As you can see, we still have the two flags, but the missing data flag (`M`) is now replaced with NA. The other flag `^` is not, as it indicates that "The value displayed is based on incomplete data" (see below).

### Flags - Weather Data

The flags index can be accessed through the built in data frame: `flags`


|code    |meaning                                                             |
|:-------|:-------------------------------------------------------------------|
|[empty] |Indicates an unobserved value                                       |
|†       |Data that is not subject to review by the National Climate Archives |
|^       |The value displayed is based on incomplete data                     |
|A       |Accumulated                                                         |
|B       |More than one occurrence and estimated                              |
|C       |Precipitation occurred, amount uncertain                            |
|E       |Estimated                                                           |
|F       |Accumulated and estimated                                           |
|L       |Precipitation may or may not have occurred                          |
|M       |Missing                                                             |
|N       |Temperature missing but known to be > 0                             |
|S       |More than one occurrence                                            |
|T       |Trace                                                               |
|Y       |Temperature missing but known to be < 0                             |
|NA      |Not Available                                                       |

### Codes
In the `normals_dl`() function, codes are associated with each variable:


```
## # A tibble: 13 × 7
##    period temp_daily_average temp_daily_average_code temp_daily_max temp_daily_max_c… temp_daily_min
##    <fct>               <dbl> <chr>                            <dbl> <chr>                      <dbl>
##  1 Jan                 -16.6 A                                -11.1 A                          -21.9
##  2 Feb                 -13.6 A                                 -8.1 A                          -19  
##  3 Mar                  -6.2 A                                 -1   A                          -11.4
##  4 Apr                   4   A                                 10.5 A                           -2.6
##  5 May                  10.6 A                                 17.8 A                            3.4
##  6 Jun                  15.9 A                                 22.4 A                            9.3
##  7 Jul                  18.5 A                                 25.2 A                           11.7
##  8 Aug                  17.7 A                                 24.9 A                           10.4
##  9 Sep                  11.8 A                                 18.9 A                            4.7
## 10 Oct                   4.1 A                                 10.4 A                           -2.2
## 11 Nov                  -5.6 A                                 -0.5 A                          -10.6
## 12 Dec                 -14   A                                 -9   A                          -19  
## 13 Year                  2.2 A                                  8.4 A                           -3.9
## # … with 1 more variable: temp_daily_min_code <chr>
```

For example, here, the code indicates that these temperature variables meet the WMO '3 and 5 rule' (no more than 3 consecutive and no more than 5 total missing for either temperature or precipitation).                                                                   


### Codes - Climate Normals

The codes index for climate normals can be accessed through the built-in data frame: `codes`


|code |meaning                                                                                                                       |
|:----|:-----------------------------------------------------------------------------------------------------------------------------|
|A    |WMO '3 and 5 rule' (i.e. no more than 3 consecutive and no more than 5 total missing for either temperature or precipitation) |
|B    |At least 25 years                                                                                                             |
|C    |At least 20 years                                                                                                             |
|D    |At least 15 years                                                                                                             |




---
title: "Reproducibility"
author: "Steffi LaZerte"
date: "2021-11-30"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Reproducibility}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---



When using data from external sources it's a good idea to take note of when data
was downloaded, which version (if possible) and with what. 

Reproducibility with `weathercan` can be achieved by taking note (or better yet, 
compiling reports) with the following information:

1. Your computer information (and date)
    - R version
2. Specific information on packages you're using
    - Citations if presenting in papers/reports
3. The stations list version

For example:


```r
# Work
library(weathercan)
s <- stations_search("Winnipeg", normals_years = "current")
w <- weather_dl(s, interval = "month", start = "2021-01-01")

# Reproducibility
stations_meta()
citation('weathercan')
devtools::session_info() # Install devtools if you don't have it
```



```
## $ECCC_modified
## [1] "2021-10-31 23:34:00 UTC"
## 
## $weathercan_modified
## [1] "2021-11-30"
```

```
## 
## To cite 'weathercan' in publications, please use:
## 
##   LaZerte, Stefanie E and Sam Albers (2018). weathercan: Download and format weather data from
##   Environment and Climate Change Canada. The Journal of Open Source Software 3(22):571.
##   doi:10.21105/joss.00571.
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
```

```
## ─ Session info  ───────────────────────────────────────────────────────────────────────────────────────────────────
##  hash: four-thirty, dolphin, water closet
## 
##  setting  value
##  version  R version 4.1.2 (2021-11-01)
##  os       Ubuntu 20.04.3 LTS
##  system   x86_64, linux-gnu
##  ui       RStudio
##  language en_CA:en
##  collate  en_CA.UTF-8
##  ctype    en_CA.UTF-8
##  tz       America/Winnipeg
##  date     2021-11-30
##  rstudio  1.4.1717 Juliet Rose (desktop)
##  pandoc   2.11.4 @ /usr/lib/rstudio/bin/pandoc/pandoc
## 
## ─ Packages ────────────────────────────────────────────────────────────────────────────────────────────────────────
##  package      * version  date (UTC) lib source
##  assertthat     0.2.1    2019-03-21 [1] CRAN (R 4.1.0)
##  bit            4.0.4    2020-08-04 [1] CRAN (R 4.1.0)
##  bit64          4.0.5    2020-08-30 [1] CRAN (R 4.1.0)
##  bitops         1.0-7    2021-04-24 [1] CRAN (R 4.1.0)
##  blob           1.2.2    2021-07-23 [1] CRAN (R 4.1.2)
##  cachem         1.0.6    2021-08-19 [1] CRAN (R 4.1.1)
##  callr          3.7.0    2021-04-20 [1] CRAN (R 4.1.0)
##  cli            3.1.0    2021-10-27 [1] CRAN (R 4.1.2)
##  codemetar      0.3.1    2021-06-02 [1] CRAN (R 4.1.0)
##  colorspace     2.0-2    2021-06-24 [1] CRAN (R 4.1.0)
##  cowplot        1.1.1    2020-12-30 [1] CRAN (R 4.1.2)
##  crayon         1.4.2    2021-10-29 [1] CRAN (R 4.1.2)
##  curl           4.3.2    2021-06-23 [1] CRAN (R 4.1.0)
##  DBI            1.1.1    2021-01-15 [1] CRAN (R 4.1.0)
##  dbplyr         2.1.1    2021-04-06 [1] CRAN (R 4.1.0)
##  desc           1.4.0    2021-09-28 [1] CRAN (R 4.1.1)
##  devtools       2.4.2    2021-06-07 [1] CRAN (R 4.1.2)
##  DiagrammeR     1.0.6.1  2020-05-08 [1] CRAN (R 4.1.0)
##  digest         0.6.28   2021-09-23 [1] CRAN (R 4.1.1)
##  dplyr        * 1.0.7    2021-06-18 [1] CRAN (R 4.1.0)
##  dygraphs       1.1.1.6  2018-07-11 [1] CRAN (R 4.1.2)
##  ellipsis       0.3.2    2021-04-29 [1] CRAN (R 4.1.0)
##  evaluate       0.14     2019-05-28 [1] CRAN (R 4.1.0)
##  fansi          0.5.0    2021-05-25 [1] CRAN (R 4.1.0)
##  farver         2.1.0    2021-02-28 [1] CRAN (R 4.1.0)
##  fastmap        1.1.0    2021-01-25 [1] CRAN (R 4.1.0)
##  fs             1.5.1    2021-11-30 [1] CRAN (R 4.1.2)
##  gdata          2.18.0   2017-06-06 [1] CRAN (R 4.1.0)
##  generics       0.1.1    2021-10-25 [1] CRAN (R 4.1.2)
##  ggplot2      * 3.3.5    2021-06-25 [1] CRAN (R 4.1.0)
##  glue         * 1.5.1    2021-11-30 [1] CRAN (R 4.1.2)
##  gtable         0.3.0    2019-03-25 [1] CRAN (R 4.1.0)
##  gtools         3.9.2    2021-06-06 [1] CRAN (R 4.1.2)
##  highr          0.9      2021-04-16 [1] CRAN (R 4.1.0)
##  hms            1.1.1    2021-09-26 [1] CRAN (R 4.1.1)
##  htmltools      0.5.2    2021-08-25 [1] CRAN (R 4.1.1)
##  htmlwidgets    1.5.4    2021-09-08 [1] CRAN (R 4.1.1)
##  httr           1.4.2    2020-07-20 [1] CRAN (R 4.1.0)
##  hunspell       3.0.1    2020-12-09 [1] CRAN (R 4.1.0)
##  igraph         1.2.8    2021-11-07 [1] CRAN (R 4.1.2)
##  jsonlite       1.7.2    2020-12-09 [1] CRAN (R 4.1.0)
##  knitr        * 1.36     2021-09-29 [1] CRAN (R 4.1.2)
##  labeling       0.4.2    2020-10-20 [1] CRAN (R 4.1.0)
##  lattice        0.20-45  2021-09-22 [4] CRAN (R 4.1.1)
##  lifecycle      1.0.1    2021-09-24 [1] CRAN (R 4.1.1)
##  lubridate    * 1.8.0    2021-10-07 [1] CRAN (R 4.1.1)
##  magrittr       2.0.1    2020-11-17 [1] CRAN (R 4.1.0)
##  memoise        2.0.1    2021-11-26 [1] CRAN (R 4.1.2)
##  munsell        0.5.0    2018-06-12 [1] CRAN (R 4.1.0)
##  naniar       * 0.6.1    2021-05-14 [1] CRAN (R 4.1.1)
##  parsedate      1.2.1    2021-04-20 [1] CRAN (R 4.1.0)
##  pillar         1.6.4    2021-10-18 [1] CRAN (R 4.1.1)
##  pingr          2.0.1    2020-06-22 [1] CRAN (R 4.1.0)
##  pkgbuild       1.2.1    2021-11-30 [1] CRAN (R 4.1.2)
##  pkgconfig      2.0.3    2019-09-22 [1] CRAN (R 4.1.0)
##  pkgdown        1.6.1    2020-09-12 [1] CRAN (R 4.1.1)
##  pkgload        1.2.3    2021-10-13 [1] CRAN (R 4.1.1)
##  prettyunits    1.1.1    2020-01-24 [1] CRAN (R 4.1.0)
##  processx       3.5.2    2021-04-30 [1] CRAN (R 4.1.0)
##  ps             1.6.0    2021-02-28 [1] CRAN (R 4.1.0)
##  purrr          0.3.4    2020-04-17 [1] CRAN (R 4.1.0)
##  R6             2.5.1    2021-08-19 [1] CRAN (R 4.1.1)
##  rappdirs       0.3.3    2021-01-31 [1] CRAN (R 4.1.0)
##  RavenR       * 2.1.4    2021-09-23 [1] CRAN (R 4.1.2)
##  rcmdcheck      1.4.0    2021-09-27 [1] CRAN (R 4.1.2)
##  RColorBrewer   1.1-2    2014-12-07 [1] CRAN (R 4.1.0)
##  Rcpp           1.0.7    2021-07-07 [1] CRAN (R 4.1.0)
##  RCurl          1.98-1.5 2021-09-17 [1] CRAN (R 4.1.2)
##  readr        * 2.1.0    2021-11-11 [1] CRAN (R 4.1.2)
##  rematch        1.0.1    2016-04-21 [1] CRAN (R 4.1.0)
##  remotes        2.4.2    2021-11-30 [1] CRAN (R 4.1.2)
##  rhub           1.1.1    2019-04-08 [1] CRAN (R 4.1.0)
##  rlang          0.4.12   2021-10-18 [1] CRAN (R 4.1.1)
##  rprojroot      2.0.2    2020-11-15 [1] CRAN (R 4.1.0)
##  RSQLite        2.2.8    2021-08-21 [1] CRAN (R 4.1.2)
##  rstudioapi     0.13     2020-11-12 [1] CRAN (R 4.1.0)
##  rvest        * 1.0.2    2021-10-16 [1] CRAN (R 4.1.1)
##  scales         1.1.1    2020-05-11 [1] CRAN (R 4.1.0)
##  sessioninfo    1.2.1    2021-11-02 [1] CRAN (R 4.1.2)
##  spelling       2.2      2020-10-18 [1] CRAN (R 4.1.0)
##  stringi        1.7.6    2021-11-29 [1] CRAN (R 4.1.2)
##  stringr      * 1.4.0    2019-02-10 [1] CRAN (R 4.1.0)
##  testthat     * 3.1.0    2021-10-04 [1] CRAN (R 4.1.1)
##  tibble         3.1.6    2021-11-07 [1] CRAN (R 4.1.2)
##  tidyhydat    * 0.5.4    2021-09-15 [1] CRAN (R 4.1.2)
##  tidyr        * 1.1.4    2021-09-27 [1] CRAN (R 4.1.1)
##  tidyselect     1.1.1    2021-04-30 [1] CRAN (R 4.1.0)
##  tzdb           0.2.0    2021-10-27 [1] CRAN (R 4.1.2)
##  urlchecker     1.0.0    2021-03-04 [1] CRAN (R 4.1.0)
##  usethis        2.1.3    2021-10-27 [1] CRAN (R 4.1.2)
##  utf8           1.2.2    2021-07-24 [1] CRAN (R 4.1.0)
##  uuid           1.0-3    2021-11-01 [1] CRAN (R 4.1.2)
##  vctrs          0.3.8    2021-04-29 [1] CRAN (R 4.1.0)
##  visdat         0.5.3    2019-02-15 [1] CRAN (R 4.1.0)
##  visNetwork     2.1.0    2021-09-29 [1] CRAN (R 4.1.2)
##  vroom          1.5.6    2021-11-10 [1] CRAN (R 4.1.2)
##  weathercan   * 0.6.2    2021-11-30 [1] local
##  whoami         1.3.0    2019-03-19 [1] CRAN (R 4.1.0)
##  withr          2.4.3    2021-11-30 [1] CRAN (R 4.1.2)
##  xfun           0.28     2021-11-04 [1] CRAN (R 4.1.2)
##  xml2           1.3.2    2020-04-23 [1] CRAN (R 4.1.0)
##  xopen          1.0.0    2018-09-17 [1] CRAN (R 4.1.0)
##  xts            0.12.1   2020-09-09 [1] CRAN (R 4.1.2)
##  zoo            1.8-9    2021-03-09 [1] CRAN (R 4.1.0)
## 
##  [1] /home/steffi/R/x86_64-pc-linux-gnu-library/4.1
##  [2] /usr/local/lib/R/site-library
##  [3] /usr/lib/R/site-library
##  [4] /usr/lib/R/library
## 
## ───────────────────────────────────────────────────────────────────────────────────────────────────────────────────
```
---
title: "weathercan and tidyhydat"
author: "Steffi LaZerte"
date: "2021-11-30"
output: rmarkdown::html_document
---



[`tidyhydat`](https://docs.ropensci.org/tidyhydat) is another R package for accessing
data from ECCC. In this case, `tidyhydat` gives access to the [National Water Data Archive
(HYDAT)](https://www.canada.ca/en/environment-climate-change/services/water-overview/quantity/monitoring/survey/data-products-services/national-archive-hydat.html).

HYDAT contains lots of data including stream flow (historical and real-time), water levels,
station metrics and information on data types and codes.

Here we'll go through a brief example of how hydrology data from `tidyhydat`
can compliment weather data from `weathercan` (and vice versa).

## Setup

### Loading packages


```r
library(weathercan)
library(tidyhydat)

library(dplyr)
library(ggplot2)
library(lubridate)
library(glue)
```

### Prep HYDAT data

`tidyhydat` data needs to be downloaded and cached locally in order to be used
(Note, this can take a while!)


```r
download_hydat()
```


## Exploring climate and hydrology

In the summer of 2020, my region (Brandon, Manitoba) experienced an incredibly
heavy rain fall event. Our downspout was ripped off the gutter and many people in the
area experienced flooding as rain poured into their basements.

Let's take a look at how this event was captured by weather and hydrometric stations
monitored by ECCC.

The event occurred in late June/early July, so let's give ourselves a two-month
range.


```r
dates <- c("2020-06-01", "2020-08-01")
```

We'll find a local Brandon weather station that has daily data for this range

```r
stations_search("brandon", interval = "day",
                starts_latest = 2020, ends_earliest = 2020)
```

```
## # A tibble: 2 × 16
##   prov  station_name station_id climate_id WMO_id TC_id   lat   lon  elev tz        interval start   end normals
##   <chr> <chr>             <dbl> <chr>       <dbl> <chr> <dbl> <dbl> <dbl> <chr>     <chr>    <dbl> <dbl> <lgl>  
## 1 MB    BRANDON A         50821 5010481     71140 YBR    49.9 -100.  409. Etc/GMT+6 day       2012  2021 FALSE  
## 2 MB    BRANDON RCS       49909 5010490     71136 PBO    49.9 -100.  409. Etc/GMT+6 day       2012  2021 FALSE  
## # … with 2 more variables: normals_1981_2010 <lgl>, normals_1971_2000 <lgl>
```

In this case "A" is for "Airport", let's go with that!


```r
rain <- weather_dl(station_ids = 50821, interval = "day", start = dates[1], end = dates[2])
```

Take a quick look:


```r
ggplot(data = rain, aes(x = date, y = total_rain)) +
  theme_bw() +
  geom_bar(stat = "identity") +
  scale_y_continuous(name = "Total Rain (mm)", expand = c(0,0))
```

<img src="unnamed-chunk-7-1.png" title="plot of chunk unnamed-chunk-7" alt="plot of chunk unnamed-chunk-7" width="80%" style="display: block; margin: auto;" />

Yikes! You can see why my downspout came off!

Now let's get some HYDAT data to compare. First we'll find a local station


```r
search_stn_name("brandon")
```

```
## # A tibble: 4 × 5
##   STATION_NUMBER STATION_NAME                                PROV_TERR_STATE_LOC LATITUDE LONGITUDE
##   <chr>          <chr>                                       <chr>                  <dbl>     <dbl>
## 1 05MH001        ASSINIBOINE RIVER AT BRANDON                MB                      49.9    -100. 
## 2 05MH006        LITTLE SOURIS RIVER NEAR BRANDON            MB                      49.7     -99.8
## 3 02OC007        MASKINONGE (LAC) A SAINT GABRIEL DE BRANDON QC                      46.3     -73.4
## 4 05MH013        ASSINIBOINE RIVER NEAR BRANDON              MB                      49.9    -100.
```

There are a couple of options, but whoops, one's from Quebec! Let's filter this
to only Manitoba and only stations with 2020 data with the `hy_stn_data_range()`
function.


```r
search_stn_name("brandon") %>%
  filter(PROV_TERR_STATE_LOC == "MB") %>%
  pull(STATION_NUMBER) %>%
  hy_stn_data_range() %>%
  filter(Year_from <= 2020, Year_to >= 2020)
```

```
##   Queried from version of HYDAT released on 2021-10-19
##    Observations:                      2
##    Station(s) returned:               1
##    Stations requested but not returned: 
##     All stations returned.
## # A tibble: 2 × 6
##   STATION_NUMBER DATA_TYPE SED_DATA_TYPE Year_from Year_to RECORD_LENGTH
##   <chr>          <chr>     <chr>             <int>   <int>         <int>
## 1 05MH001        H         <NA>               2014    2020             7
## 2 05MH001        Q         <NA>               1906    2020            75
```

Hmm, let's see what kind of data is available by looking at the included
`hy_data_types` data frame.


```r
filter(hy_data_types, DATA_TYPE %in% c("H", "Q"))
```

```
## # A tibble: 2 × 3
##   DATA_TYPE DATA_TYPE_EN DATA_TYPE_FR 
##   <chr>     <chr>        <chr>        
## 1 H         Water Level  Niveaux d'eau
## 2 Q         Flow         Debit
```

Great! We have both flow and water level data for a station number "05MH001",
"Assiniboine River at Brandon".


Let's grab the flow and water level data for this station.

```r
flow <- hy_daily_flows(station_number = "05MH001",
                       start_date = dates[1], end_date = dates[2])
level <- hy_daily_levels(station_number = "05MH001",
                         start_date = dates[1], end_date = dates[2])
```

### Ploting rain and flow


```r
g <- ggplot() +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  geom_bar(data = rain, aes(x = date, y = (total_rain * 2)), stat = "identity",
           alpha = 0.7, fill = "cornflowerblue") +
  geom_line(data = flow, aes(x = Date, y = Value)) +
  scale_y_continuous(name = bquote(Total~Flow~(m^3/s)), expand = c(0, 0),
                     limits = c(0, max(flow$Value * 1.1)),
                     sec.axis = sec_axis(trans = ~ . / 2, name = "Total Rain (mm)"))
g
```

<img src="unnamed-chunk-12-1.png" title="plot of chunk unnamed-chunk-12" alt="plot of chunk unnamed-chunk-12" width="80%" style="display: block; margin: auto;" />


Interesting, looks like there's a bit of a lag between the rain event and the
dramatic increase in water flow in the Assiniboine (unsurprisingly, this is
called "lag to peak").

Let's add a bit of information about this lag to peak.


```r
d <- data.frame(dates = c(rain$date[which.max(rain$total_precip)],
                          flow$Date[which.max(flow$Value)]),
                y = max(flow$Value) + 5)

g +
  geom_path(data = d, aes(x = dates, y = y),
              arrow = arrow(length = unit(0.25, "lines"), ends = "both", type = "closed")) +
  annotate(geom = "text",
           x = d$dates[1] + (d$dates[2] - d$dates[1])/2,
           y = d$y[1] + 10,
           label = glue("{d$dates[2] - d$dates[1]}-day delay"))
```

<img src="unnamed-chunk-13-1.png" title="plot of chunk unnamed-chunk-13" alt="plot of chunk unnamed-chunk-13" width="80%" style="display: block; margin: auto;" />

We can expect a lag like this because much of the flow being captured by the Brandon HYDAT
station is from precipitation in the upstream catchment area (not only from
local contributions), which takes time to travel.


### Ploting rain and water level


```r
g <- ggplot() +
  theme_bw() +
  theme(axis.title.x = element_blank()) +
  geom_bar(data = rain,
           aes(x = date, y = (total_rain/65) + min(level$Value)),
           stat = "identity", alpha = 0.7, fill = "cornflowerblue") +
  geom_line(data = level, aes(x = Date, y = Value)) +
  scale_y_continuous(name = "Water Level (m)", expand = c(0, 0),
                     sec.axis = sec_axis(trans = ~ (. - min(level$Value)) * 65,
                                         name = "Total Rain (mm)")) +
  coord_cartesian(ylim = c(min(level$Value), max(level$Value)*1.001))
g
```

<img src="unnamed-chunk-14-1.png" title="plot of chunk unnamed-chunk-14" alt="plot of chunk unnamed-chunk-14" width="80%" style="display: block; margin: auto;" />

Again there looks to be a lag, let's see if it's the same as before.


```r
d <- data.frame(dates = c(rain$date[which.max(rain$total_precip)],
                          level$Date[which.max(level$Value)]),
                y = max(level$Value)*1.00025)

g +
  geom_path(data = d, aes(x = dates, y = y),
              arrow = arrow(length = unit(0.25, "lines"), ends = "both", type = "closed")) +
  annotate(geom = "text",
           x = d$dates[1] + (d$dates[2] - d$dates[1])/2,
           y = d$y[1] * 1.0002,
           label = glue("{d$dates[2] - d$dates[1]}-day delay"))
```

<img src="unnamed-chunk-15-1.png" title="plot of chunk unnamed-chunk-15" alt="plot of chunk unnamed-chunk-15" width="80%" style="display: block; margin: auto;" />


### Ploting flow and water level

Looks like the flow and water level match up, perhaps we should take a closer look.


```r
ggplot() +
  theme_bw() +
  theme(legend.position = c(0.8, 0.8)) +
  geom_line(data = flow, aes(x = Date, y = Value, colour = "Flow"),
            size = 2) +
  geom_line(data = level, size = 1,
            aes(x = Date, y = (Value - min(Value) + 0.1) * 130, colour = "Level")) +
  scale_y_continuous(bquote(Total~Flow~(m^3/s)), expand = c(0, 0),
                     sec.axis = sec_axis(trans = ~ ./130 + min(level$Value) - 0.1,
                                         name = "Water Level (m)")) +
  scale_colour_manual(name = "Type",
                      values = c("Flow" = "cornflowerblue",
                                 "Level" = "grey30"))
```

<img src="unnamed-chunk-16-1.png" title="plot of chunk unnamed-chunk-16" alt="plot of chunk unnamed-chunk-16" width="80%" style="display: block; margin: auto;" />

Almost a perfect match between water level and flow (which makes sense).


> Hopefully this short article gives you a sense of how you might combine different
> types of ECCC data gathered via different R packages for a more comprehensive
> look at the world around.


---
title: "Mapping weather data"
author: "Steffi LaZerte"
date: "`r Sys.Date()`"
output: rmarkdown::html_document
---

<style>
div.leaflet-popup-content-wrapper {
  width: 700px;
  height: 100%;
}

div.leaflet-popup-tip-container {
  opacity: 0;
}

div.leaflet-popup-content {
  width: 90% !important;
}

img {
  max-width: 90% !important;
  min-width: 90% !important;
  border: none;
}

/* Fix the NA mis-alignment in the legend */
div.info.legend.leaflet-control br {
  clear: both;
}
</style>

This article is based on the blog post [Integrating data from weathercan](https://ropensci.org/blog/2018/03/06/weathercan/) written for [ROpenSci](https://ropensci.org) March 6th 2018.

In that article I demonstrated how we can incorporate data from `weathercan` into spatial visualizations. In this article I'd like to take that even further and show you how you can create interactive maps which highlight spatial variability in weather data.

Here, we'll take a look at annual temperatures throughout different Eco Regions in Manitoba, Canada.

```{r, include = FALSE}
knitr::opts_chunk$set(cache = FALSE)
```

## Setup

### Using extra CSS styles

The map that we create here may look different for you unless you include the tweaks to the CSS styles I have made. If you're using RMarkdown, you can supply these as a custom .css file, or inline with the `<style>` and `</style>` tags (like below).

```
<style>
div.leaflet-popup-content-wrapper {
  width: 700px;
  height: 100%;
}

div.leaflet-popup-tip-container {
  opacity: 0;
}

div.leaflet-popup-content {
  width: 90% !important;
}

img {
  max-width: 90% !important;
  min-width: 90% !important;
  border: none;
}

/* Fix the NA mis-alignment in the legend */
div.info.legend.leaflet-control br {
  clear: both;
}
</style>
```


### Loading packages

```{r, message = FALSE}
library(dplyr)
library(purrr)
library(stringr)
library(ggplot2)
library(tidyr)
library(weathercan)
library(leaflet)
library(sf)
library(htmltools)
```

### Download Manitoban Eco Regions shapefile

```{r, eval = FALSE}
download.file("http://mli2.gov.mb.ca/environment/shp_zip_files/env_ecological_areas_py_shp.zip",
              destfile = "./mapping_shp/ecological_shp.zip")
unzip("./mapping_shp/ecological_shp.zip")
file.remove("./mapping_shp/ecological_shp.zip")
```

### Download Manitoba weather data

We'll select all currently operating stations (`end >= 2018`) and download the daily weather data for 2017:

```{r, message = FALSE}
mb <- filter(stations(), prov == "MB", interval == "day", end >= 2018)
w <- weather_dl(mb$station_id, start = "2017-01-01", end = "2017-12-31", interval = "day")
```

## Calculating summaries
In this section, we'll summarize our data and create the means of showcasing it in our map. We'll calculate some basic summaries and we'll style some popups to contain this information (including the figures we'll create later).

We'll do a station-specific summaries/pop-ups for when users click on a station marker, and region-specific ones for when they click on the region polygon.

### Station summaries
```{r}
mb_stations <- w %>%
  group_by(station_id) %>%
  mutate(n = n(),
         n_missing = sum(is.na(mean_temp)),
         mean_temp = mean(mean_temp, na.rm = TRUE)) %>%
  filter((n - n_missing) > 0) %>% # Only keep stations with some temperature data
  select(station_name, station_id, lat, lon, mean_temp, n, n_missing) %>%
  distinct() %>%
  mutate(station_name = tools::toTitleCase(tolower(station_name)),
         info = paste0("<h3>Station: ", station_name, " (", station_id, ")</h2>",
                       "<hr>",
                       "<div>",
                       "<strong>Mean Temperature: </strong>", round(mean_temp, 1), "C<br>",
                       "<strong>No. days with data:  </strong>", n-n_missing, "<br>",
                       "<strong>No. days total:  </strong>", n, "</div>",
                       "<img src = './mapping_svg/", station_id, ".svg'>"),
         pretty_name = map(station_name, 
                           ~HTML(paste0("<strong>Station: </strong>", .x)))) %>%
  ungroup() %>%
  st_as_sf(coords = c("lon", "lat"), crs = "+proj=longlat")
```

### Eco Region summaries

Before we summarize this data, we'll filter the Eco Regions to just Manitoba and will join this to the stations data (after transforming the stations data to the same CRS), so we can figure out which stations belong to which regions.

```{r}
mb_ecoregions <- st_read("./mapping_shp/env_ecological_areas.shp") %>%
  filter(MANITOBA == "yes") %>%
  # Get the larger-scale regions and combine
  group_by(ECOREGION, REGION_NAM, REGION_NOM) %>%
  summarize()

mb_stations <- st_transform(mb_stations, crs = st_crs(mb_ecoregions))
  
mb_ecoregions <- st_join(mb_ecoregions, mb_stations) %>%
  group_by(ECOREGION, REGION_NAM, REGION_NOM) %>%
  summarize(n_stations = length(unique(station_id[!is.na(station_id)])),
            mean_temp = mean(mean_temp, na.rm = TRUE),
            mean_temp = replace(mean_temp, is.nan(mean_temp), NA)) %>%
  mutate(info = paste0("<h3>Region: ", REGION_NAM, "/", REGION_NOM, " (", ECOREGION, ")</h2>",
                       "<hr>",
                       "<div>",
                       "<strong>Mean Temperature: </strong>", round(mean_temp, 1), 
                       if_else(!is.na(mean_temp), "C<br>", "<br>"),
                       "<strong>No. stations:  </strong>", n_stations, "</div>"),
         info = if_else(n_stations > 0, 
                        paste0(info, "<img src = './mapping_svg/", ECOREGION, ".svg'>"),
                        info),
         pretty_name = map(REGION_NAM, 
                           ~HTML(paste0("<strong>Region: </strong>", .x)))) %>%
  ungroup()
```

```{r, include = FALSE}
# The above details work locally
# This chunk makes the figure paths work online for the weathercan website
mb_ecoregions <- mutate(mb_ecoregions, 
                        info = str_replace(info, "./mapping_svg/",
    "https://raw.githubusercontent.com/ropensci/weathercan/master/vignettes/articles/mapping_svg/"),
    info = str_replace(info, ".svg'>", ".svg?sanitize=true'>"))
                        
    
mb_stations <- mutate(mb_stations, info = str_replace(info, "./mapping_svg/",
    "https://raw.githubusercontent.com/ropensci/weathercan/master/vignettes/articles/mapping_svg/"),
    info = str_replace(info, ".svg'>", ".svg?sanitize=true'>"))
```



### Figures

We'll set up some functions to create and save figures for our map

```{r}
plot_station_fig <- function(d, station_id) {
  g <- ggplot(d, aes(x = date, y = mean_temp)) +
    theme_bw() +
    geom_line(na.rm = TRUE) +
    labs(x = "Date", y = "Mean Daily Temperature (C)")
  ggsave(paste0("./mapping_svg/", station_id, ".svg"), plot = g,
         width = 6, height = 3, dpi = 100)
}

plot_region_fig <- function(d, region) {
  g <- ggplot(d, aes(x = date, y = mean_temp, 
                group = station_name, colour = station_name)) +
    theme_bw() +
    geom_line(na.rm = TRUE) +
    scale_colour_viridis_d(end = 0.8) +
    labs(x = "Date", y = "Mean Daily Temperature (C)", colour = "Stations")
  ggsave(paste0("./mapping_svg/", region, ".svg"), plot = g,
         width = 6, height = 3, dpi = 100)
}
```

Now we'll apply these functions to our data. Note that this `figs` object isn't important, it's just used as a way to loop through the data and save the figures as svg files.
```{r}
figs <- st_join(mb_stations, mb_ecoregions) %>%
  st_set_geometry(NULL) %>%
  left_join(select(w, station_id, date, mean_temp), by = "station_id") %>%
  nest(-ECOREGION) %>%
  mutate(fig_region = map2(data, ECOREGION, ~plot_region_fig(.x, .y))) %>%
  unnest(data) %>%
  nest(-ECOREGION, -station_id) %>%
  mutate(fig_station = map2(data, station_id, ~plot_station_fig(.x, .y)))
```


## Mapping

Finally, we're ready to create our map!

We'll start by transforming our data into WGS84 for leaflet, then we'll create a palette and get some icons...
```{r}
# Required for Leaflet
mb_ecoregions <- st_transform(mb_ecoregions, crs = 4326)
mb_stations <- st_transform(mb_stations, crs = 4326)

# Setup Palette for polygons
pal_eco <- colorNumeric(palette = "viridis",
                          domain = mb_ecoregions$mean_temp)

# Get icons for stations (red if above average, blue if below)
station_icons <- awesomeIcons(iconColor = "black", library = "ion",
                              markerColor = ifelse(mb_stations$mean_temp > 
                                                     mean(mb_stations$mean_temp, 
                                                          na.rm = TRUE), 
                                                   "red", "blue"))
```

Now for the real magic!


```{r}
leaflet(width = "750px", height = "85vh") %>% 
  addTiles() %>%
  addPolygons(data = mb_ecoregions,
              color = "#444444", weight = 1, opacity = 1, fillOpacity = 0.5,
              fillColor = ~pal_eco(mean_temp),
              label = ~pretty_name, popup = ~info,
              popupOptions = popupOptions(keepInView = TRUE),
              highlightOptions = highlightOptions(bringToFront = TRUE, 
                                                  fillOpacity = 1)) %>%
  addAwesomeMarkers(data = mb_stations, group = "Stations",
                    icon = station_icons,
                    label = ~pretty_name, popup = ~info,
                    popupOptions = popupOptions(keepInView = TRUE)) %>%
  addLegend("bottomright", pal = pal_eco,
            values = mb_ecoregions$mean_temp,
            title = "Mean region temperature",
            labFormat = labelFormat(suffix = " C")) %>%
  addLegend("bottomright", 
            title = "Station temperature", 
            colors = c("#d63e2a", "#37a7da"), opacity = 1, 
            labels = c("Above the average", "Below the average")) %>%
  addLayersControl(
    overlayGroups = "Stations",
    options = layersControlOptions(collapsed = FALSE))


```

I think this is an interesting way of looking at the Eco Regions in Manitoba. 

First we see a clear (and expected) pattern of decreasing temperatures with latitude (South to North). Further, we also see a bit of a South West to North East pattern.

Looking at the Eco Regions like this gives a clear idea of some of the parameters that make these Eco Regions distinct. For example, the patch of lands in the lower left corner of the province are the Mid-Boreal Uplands and the Boreal Transition. See how much cooler they are (on average) than the rest of southern Manitoba. 

These southern Boreal regions represent a distinct ecology in southern Manitoba. If you zoom in on the map, you'll see that they also hold two parks: Duck Mountain Provincial Park and Riding Mountain National Park.

---
title: "Meteoland"
author: "Steffi LaZerte"
date: "`r Sys.Date()`"
output: rmarkdown::html_document
---

[`meteoland`](https://github.com/vegmod/meteoland) is a package for interpolating meteorological data over spatial scales. As of version v0.7.9 they support transforming and important `weathercan` output for use in spatial interpolations. In this article we will go over a hypothetical example.

## Setup
First we'll load the packages and find the stations we're interested in and 
```{r}
library(meteoland)
library(weathercan)
library(dplyr)
library(tidyr)
library(sf)
library(ggplot2)
library(rnaturalearth)

s <- stations() %>%
  filter(prov == "MB", interval != "month", start <= 2015, end >= 2015, !is.na(elev)) %>%
  group_by(station_id) %>%
  mutate(n = length(interval)) %>%
  filter(n == 2)

s_map <- s %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3347)

mb <- ne_states(country = "Canada", returnclass = "sf") %>%
  filter(name == "Manitoba") %>%
  st_transform(3347)

ggplot() +
  geom_sf(data = mb) +
  geom_sf(data = s_map)
```

Let's focus on northern Manitoba

```{r}
s <- filter(s, lat > 55)
s_map <- s %>%
  st_as_sf(coords = c("lon", "lat"), crs = 4326) %>%
  st_transform(3347)

ggplot() +
  geom_sf(data = mb) +
  geom_sf(data = s_map)
```


## Download data

```{r}
mb_hr <- weather_dl(station_id = unique(s$station_id), interval = "hour",
                    start = "2015-01-01", end = "2015-12-31", verbose = TRUE)
mb_day <- weather_dl(station_id = unique(s$station_id), interval = "day",
                     start = "2015-01-01", end = "2015-12-31", verbose = TRUE)
```

## Pass to `meteoland`

First we'll reshape our hourly and daily data into a `meteoland` interpolations data object.
```{r}
mb_north <- reshapeweathercan(mb_hr, mb_day, output = "MeteorologyInterpolationData")
```

We can get a sense of the data coverage (number of dates with data per station per variable).
```{r}
interpolation.coverage(mb_north, type = 'spatial') %>%
  head()
```

Or the number of stations with data per date per variable.
```{r}
interpolation.coverage(mb_north, type = 'temporal') %>%
  head()
```

Next we have to calibrate the variable we're interested in (here `Tmin`)
```{r}
tmin <- interpolation.calibration(mb_north, variable = "Tmin",
                                  N_seq = 20,
                                  alpha_seq = seq(5, 10, by = 1),
                                  verbose = TRUE)

mb_north@params$N_MinTemperature = tmin$N
mb_north@params$alpha_MinTemperature = tmin$alpha
```

Next we cross-validate the data
```{r}
cv <- interpolation.cv(mb_north, verbose = TRUE)
summary(cv)
```

We create a dummy `SpatialPointsTopography` object representing the points (in this case) that we wish to interpolate over. Note that here I'm using a mean elevation as a placeholder as I don't have actual elevation values. If you have slope and aspect, even better. Remember, interpolation is only as good as the data you give it!

```{r}
interp <- expand.grid(lat = seq(56, 58.5, 0.25),
                      lon = seq(-101, -95, 0.25), 
                      elev = mean(s_map$elev))

interp <- SpatialPointsTopography(as.matrix(interp[, c("lon", "lat")]),
                                  elevation = interp$elev,
                                  proj4string = CRS("+proj=longlat +ellps=WGS84"))
```

Now for the actual interpolation. Here we interpolate over all dates in the range of the original data.
```{r}
new_interp <- interpolationpoints(mb_north, interp)
```

We can summarize and plot our interpolations
```{r}
map <- summarypoints(new_interp, var = "MinTemperature")

mb_interp <- st_as_sf(map) %>%
  rename_at(.vars = vars(contains("matrix")), ~"min_temp") %>%
  mutate(station_id = NA, interp = TRUE)

mb_stations <- mb_day %>%
  group_by(station_id, lat, lon) %>%
  summarize(min_temp = mean(min_temp, na.rm = TRUE)) %>%
  st_as_sf(coords = c("lon", "lat"), crs = st_crs(mb_interp)) %>%
  mutate(interp = FALSE)

map_sf <- rbind(mb_interp, mb_stations) %>%
  mutate(interp = factor(interp, levels = c("TRUE", "FALSE")))

ggplot() +
  geom_sf(data = mb) +
  geom_sf(data = map_sf, aes(fill = min_temp, shape = interp, colour = interp), size = 4) +
  scale_shape_manual(values = c(21,23)) +
  coord_sf(ylim = c(2100000, 2550000))
```
---
title: "Using with the tidyverse"
author: "Sam Albers"
date: "`r Sys.Date()`"
output: rmarkdown::html_document
---

## Source material
This vignette is adapted very heavily from Hadley Wickham's incredible [*R for Data Science*](http://r4ds.had.co.nz/) book. You should support Hadley and the work he does by buying [it](https://www.amazon.com/Data-Science-Transform-Visualize-Model/dp/1491910399).

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warnings = FALSE, message = FALSE, fig.width = 10, fig.height = 6)
```

## Packages
In addition to *weathercan*, you'll need several packages from the *tidyverse*  to complete the following analysis.

```{R pck}
library(weathercan)
library(dplyr)
library(tidyr)
library(ggplot2)
library(broom)
library(modelr)
library(purrr)
```

## Using weathercan to load in data

Your first decision that you need to make when analyzing data from weather stations across canada is to determine for which stations you'd like to query from Environment and Climate Change Canada. In this example, to keep processing time low, we will query two stations with very long records that happen to be far apart. To make that choice we can use (tidyverse)[http://tidyverse.org/] tools and the included `stations()` function to access the data frame of stations in this package:

```{R stn_pick, echo = TRUE}
stations() %>%
  filter(station_id %in% c(707, 4859, 6693,5397, 2315),
         interval == "day") %>%
  select(prov, station_name, station_id, start, end)
```

These two weather stations will be our test data for this vignette. You can broaden or expand your analysis by choosing different or more station. Our next step is to use the `weather_dl()` function to load in the data. 

The following will take quite some time to download as it is downloading over 100 years of daily data for 5 stations.

```{R load_in, echo = TRUE}
pancan_df <- weather_dl(station_ids = c(707, 4859, 6693,5397, 2315), 
                        interval = "day") %>%
  filter(year >= 1920) %>%
  select(station_name, station_id, prov, lat, lon, elev, climate_id, WMO_id, TC_id, mean_temp, date)
```

## Plot the data
```{r raw_plt}
ggplot(pancan_df, aes(x = date, y = mean_temp, colour = station_name)) +
  geom_point() +
  geom_line()
```

This is quite a large dataset. 

## Creating list-columns
```{r nesting}
pancan_df_nest <- pancan_df %>%
  group_by(station_name, station_id, prov, lat, lon, elev, climate_id, WMO_id, TC_id) %>%
  nest()
pancan_df_nest
```

## Fit some models

Define the model

```{r mod_def}
clim_model <- function(df) {
  lm(mean_temp ~ date, data = df)
}
```

Run the model with the existing data

```{r add_lm}
pancan_df_nest <- pancan_df_nest %>% 
  mutate(model = map(data, clim_model))
pancan_df_nest
```

Then add the residuals to the model

```{r add_resid}
pancan_df_nest <- pancan_df_nest %>% 
  mutate(model = map(data, clim_model),
         resids = map2(data, model, add_residuals)) 
pancan_df_nest
```

## Working with list-columns
We can unnest the results then plot them

### `unnest()`
```{r resid}
resids <- unnest(pancan_df_nest, resids)
resids


ggplot(data = resids, aes(date, resid)) +
  geom_line(aes(group = station_name), alpha = 1 / 3) + 
  geom_point() +
  geom_hline(yintercept = 0) +
  facet_wrap(~ station_name, ncol = 1)
```


### Using broom
```{r broom}
glance_df <- pancan_df_nest %>% 
  mutate(glance = map(model, broom::glance)) %>% 
  unnest(glance, .drop = TRUE) %>%
  select(station_name, prov, r.squared, p.value, AIC)
```

```{r, echo = FALSE}
knitr::kable(glance_df)
```


## Looking at the predictions
```{r pred}
preds <- pancan_df_nest %>% 
  mutate(model = map(data, clim_model),
         preds = map2(data, model, add_predictions)) %>%
  unnest(preds)
preds

ggplot(data = preds, aes(x = date, y = mean_temp, colour = station_name)) +
  geom_point() +
  geom_line(aes(y = pred)) +
  facet_wrap(~ station_name, scales = "free_y", ncol = 1)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weather.R
\name{weather_dl}
\alias{weather_dl}
\alias{weather}
\title{Download weather data from Environment and Climate Change Canada}
\usage{
weather_dl(
  station_ids,
  start = NULL,
  end = NULL,
  interval = "hour",
  trim = TRUE,
  format = TRUE,
  string_as = NA,
  time_disp = "none",
  stn = NULL,
  encoding = "UTF-8",
  list_col = FALSE,
  verbose = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{station_ids}{Numeric/Character. A vector containing the ID(s) of the
station(s) you wish to download data from. See the \code{\link{stations}}
data frame or the \code{\link{stations_search}} function to find IDs.}

\item{start}{Date/Character. The start date of the data in YYYY-MM-DD format
(applies to all stations_ids). Defaults to start of range.}

\item{end}{Date/Character. The end date of the data in YYYY-MM-DD format
(applies to all station_ids). Defaults to end of range.}

\item{interval}{Character. Interval of the data, one of "hour", "day",
"month".}

\item{trim}{Logical. Trim missing values from the start and end of the
weather dataframe. Only applies if \code{format = TRUE}}

\item{format}{Logical. If TRUE, formats data for immediate use. If FALSE,
returns data exactly as downloaded from Environment and Climate Change
Canada. Useful for dealing with changes by Environment Canada to the format
of data downloads.}

\item{string_as}{Character. What value to replace character strings in a
numeric measurement with. See Details.}

\item{time_disp}{Character. Either "none" (default) or "UTC". See details.}

\item{stn}{DEFUNCT. Now use \code{stations_dl()} to update internal data and
\code{stations_meta()} to check the date it was last updated.}

\item{encoding}{Character. Text encoding for download.}

\item{list_col}{Logical. Return data as nested data set? Defaults to FALSE.
Only applies if \code{format = TRUE}}

\item{verbose}{Logical. Include progress messages}

\item{quiet}{Logical. Suppress all messages (including messages regarding
missing data, etc.)}
}
\value{
A tibble with station ID, name and weather data.
}
\description{
Downloads data from Environment and Climate Change Canada (ECCC) for one or
more stations. For details and units, see the glossary vignette
(\code{vignette("glossary", package = "weathercan")}) or the glossary online
\url{https://climate.weather.gc.ca/glossary_e.html}.
}
\details{
Data can be returned 'raw' (format = FALSE) or can be formatted.
Formatting transforms dates/times to date/time class, renames columns, and
converts data to numeric where possible. If character strings are contained
in traditionally numeric fields (e.g., weather speed may have values such
as "< 30"), they can be replaced with a character specified by \code{string_as}.
The default is NA. Formatting also replaces data associated with certain
flags with NA (M = Missing).

Start and end date can be specified, but if not, it will default to the
start and end date of the range (this could result in downloading a lot of
data!).

For hourly data, timezones are always "UTC", but the actual times are
either local time (default; \code{time_disp = "none"}), or UTC (\code{time_disp = "UTC"}). When \code{time_disp = "none"}, times reflect the local time without
daylight savings. This means that relative measures of time, such as
"nighttime", "daytime", "dawn", and "dusk" are comparable among stations in
different timezones. This is useful for comparing daily cycles. When
\code{time_disp = "UTC"} the times are transformed into UTC timezone. Thus
midnight in Kamloops would register as 08:00:00 (Pacific time is 8 hours
behind UTC). This is useful for tracking weather events through time, but
will result in odd 'daily' measures of weather (e.g., data collected in the
afternoon on Sept 1 in Kamloops will be recorded as being collected on Sept
2 in UTC).

Files are downloaded from the url stored in
\code{getOption("weathercan.urls.weather")}. To change this location use
\code{options(weathercan.urls.weather = "your_new_url")}.

Data is downloaded from ECCC as a series of files which are then bound
together. Each file corresponds to a different month, or year, depending on
the interval. Metadata (station name, lat, lon, elevation, etc.) is
extracted from the start of the most recent file (i.e. most recent dates)
for a given station. Note that important data (i.e. station name, lat, lon)
is unlikely to change between files (i.e. dates), but some data may or may
not be available depending on the date of the file (e.g., station operator
was added as of April 1st 2018, so will be in all data which includes dates
on or after April 2018).
}
\examples{
\dontshow{if (check_eccc()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

kam <- weather_dl(station_ids = 51423,
                  start = "2016-01-01", end = "2016-02-15")

stations_search("Kamloops A$", interval = "hour")
stations_search("Prince George Airport", interval = "hour")

kam.pg <- weather_dl(station_ids = c(48248, 51423),
                     start = "2016-01-01", end = "2016-02-15")

library(ggplot2)

ggplot(data = kam.pg, aes(x = time, y = temp,
                          group = station_name,
                          colour = station_name)) +
       geom_line()
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{normals_measurements}
\alias{normals_measurements}
\title{List of climate normals measurements for each station}
\format{
A data frame with 113,325 rows and 5 variables:
\describe{
\item{prov}{Province}
\item{station_name}{Station Name}
\item{climate_id}{Climate ID}
\item{normals}{Year range of climate normals}
\item{measurement}{Climate normals measurement available for this station}
}
}
\usage{
normals_measurements
}
\description{
A data frame listing the climate normals measurements available for each
station.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{pg}
\alias{pg}
\title{Hourly weather data for Prince George}
\format{
An example dataset of hourly weather data for Prince George:
\describe{
\item{station_name}{Station name}
\item{station_id}{Environment Canada's station ID number. Required for
downloading station data.}
\item{prov}{Province}
\item{lat}{Latitude of station location in degree decimal format}
\item{lon}{Longitude of station location in degree decimal format}
\item{date}{Date}
\item{time}{Time}
\item{year}{Year}
\item{month}{Month}
\item{day}{Day}
\item{hour}{Hour}
\item{qual}{Data quality}
\item{weather}{The state of the atmosphere at a specific time.}
\item{hmdx}{Humidex}
\item{hmdx_flag}{Humidex data flag}
\item{pressure}{Pressure (kPa)}
\item{pressure_flag}{Pressure data flag}
\item{rel_hum}{Relative humidity}
\item{rel_hum_flag}{Relative humidity data flag}
\item{temp}{Temperature}
\item{temp_dew}{Dew Point Temperature}
\item{temp_dew_flag}{Dew Point Temperatureflag}
\item{visib}{Visibility (km)}
\item{visib_flag}{Visibility data flag}
\item{wind_chill}{Wind Chill}
\item{wind_chill_flag}{Wind Chill flag}
\item{wind_dir}{Wind Direction (10's of degrees)}
\item{wind_dir_flag}{wind Direction Flag}
\item{wind_spd}{Wind speed km/hr}
\item{wind_spd_flag}{Wind speed flag}
\item{elev}{Elevation (m)}
\item{climate_id}{Climate identifier}
\item{WMO_id}{World Meteorological Organization Identifier}
\item{TC_id}{Transport Canada Identifier}
}
}
\source{
\url{https://climate.weather.gc.ca/index_e.html}
}
\usage{
pg
}
\description{
Downloaded with \code{\link{weather}()}. Terms are more thoroughly defined
here \url{https://climate.weather.gc.ca/glossary_e.html}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{kamloops}
\alias{kamloops}
\title{Hourly weather data for Kamloops}
\format{
An example dataset of hourly weather data for Kamloops:
\describe{
\item{station_name}{Station name}
\item{station_id}{Environment Canada's station ID number. Required for
downloading station data.}
\item{prov}{Province}
\item{lat}{Latitude of station location in degree decimal format}
\item{lon}{Longitude of station location in degree decimal format}
\item{date}{Date}
\item{time}{Time}
\item{year}{Year}
\item{month}{Month}
\item{day}{Day}
\item{hour}{Hour}
\item{qual}{Data quality}
\item{weather}{The state of the atmosphere at a specific time.}
\item{hmdx}{Humidex}
\item{hmdx_flag}{Humidex data flag}
\item{pressure}{Pressure (kPa)}
\item{pressure_flag}{Pressure data flag}
\item{rel_hum}{Relative humidity}
\item{rel_hum_flag}{Relative humidity data flag}
\item{temp}{Temperature}
\item{temp_dew}{Dew Point Temperature}
\item{temp_dew_flag}{Dew Point Temperature flag}
\item{visib}{Visibility (km)}
\item{visib_flag}{Visibility data flag}
\item{wind_chill}{Wind Chill}
\item{wind_chill_flag}{Wind Chill flag}
\item{wind_dir}{Wind Direction (10's of degrees)}
\item{wind_dir_flag}{wind Direction Flag}
\item{wind_spd}{Wind speed km/hr}
\item{wind_spd_flag}{Wind speed flag}
\item{elev}{Elevation (m)}
\item{climate_id}{Climate identifier}
\item{WMO_id}{World Meteorological Organization Identifier}
\item{TC_id}{Transport Canada Identifier}
}
}
\source{
\url{https://climate.weather.gc.ca/index_e.html}
}
\usage{
kamloops
}
\description{
Downloaded with \code{\link{weather}()}. Terms are more thoroughly defined
here \url{https://climate.weather.gc.ca/glossary_e.html}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{glossary}
\alias{glossary}
\title{Glossary of units and terms}
\format{
A data frame with 77 rows and 5 variables:
\describe{
\item{interval}{Data interval type, 'hour', 'day', or 'month'.}
\item{ECCC_name}{Original column name when downloaded directly from ECCC}
\item{weathercan_name}{R-compatible name given when downloaded with the
\code{weather_dl()} function using the default argument \code{format =
  TRUE}.}
\item{units}{Units of the measurement.}
\item{ECCC_ref}{Link to the glossary or reference page on the ECCC
website.}
}
}
\usage{
glossary
}
\description{
A reference dataset matching information on columns in data downloaded using
the \code{weather_dl()} function. Indicates the units of the data, and
contains a link to the ECCC glossary page explaining the measurement.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{codes}
\alias{codes}
\title{Meaning of climate normal 'codes'}
\format{
A data frame with 4 rows and 2 variables:
\describe{
\item{code}{Code}
\item{meaning}{Explanation of the code}
}
}
\usage{
codes
}
\description{
A reference dataset containing \code{codes} matched to their meaning. Data
downloaded using the \code{normals_dl()} function contains columns indicating
\code{code}. These are presented here for interpretation.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/normals.R
\name{normals_dl}
\alias{normals_dl}
\title{Download climate normals from Environment and Climate Change Canada}
\usage{
normals_dl(
  climate_ids,
  normals_years = "1981-2010",
  format = TRUE,
  stn = NULL,
  verbose = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{climate_ids}{Character. A vector containing the Climate ID(s) of the
station(s) you wish to download data from. See the \code{\link{stations}}
data frame or the \code{\link{stations_search}} function to find Climate
IDs.}

\item{normals_years}{Character. The year range for which you want climate
normals. Default "1981-2010".}

\item{format}{Logical. If TRUE (default) formats measurements to numeric and
date accordingly. Unlike \code{weather_dl()}, \code{normals_dl()} will always format
column headings as normals data from ECCC cannot be directly made into a
data frame without doing so.}

\item{stn}{DEFUNCT. Now use \code{stations_dl()} to update internal data and
\code{stations_meta()} to check the date it was last updated.}

\item{verbose}{Logical. Include progress messages}

\item{quiet}{Logical. Suppress all messages (including messages regarding
missing data, etc.)}
}
\value{
tibble with nested normals and first/last frost data
}
\description{
Downloads climate normals from Environment and Climate Change Canada (ECCC)
for one or more stations (defined by \code{climate_id}s). For details and units,
see the \code{\link{glossary_normals}} data frame or the \code{glossary_normals} vignette:
\code{vignette("glossary_normals", package = "weathercan")}
}
\details{
Climate normals from ECCC include two types of data, averages by
month for a variety of measurements as well as data relating to the
frost-free period. Because these two data sources are quite different, we
return them as nested data so the user can extract them as they wish. See
examples for how to use the \code{unnest()} function from the
\href{https://tidyr.tidyverse.org/}{\code{tidyr}}
package to extract the two different datasets.

The data also returns a column called \code{meets_wmo} this reflects whether or
not the climate normals for this station met the WMO standards for
temperature and precipitation (i.e. both have code >= A). Each measurement
column has a corresponding \verb{_code} column which reflects the data quality
of that measurement (see the \href{https://climate.weather.gc.ca/doc/Canadian_Climate_Normals_1981_2010_Calculation_Information.pdf}{1981-2010 ECCC calculations document}
or the \href{https://climate.weather.gc.ca/doc/Canadian_Climate_Normals_1971_2000_Calculation_Information.pdf}{1971-2000 ECCC calculations document}
for more details)

Climate normals are downloaded from the url stored in option
\code{weathercan.urls.normals}. To change this location use:
\code{options(weathercan.urls.normals = "your_new_url")}.
}
\examples{
\dontshow{if (check_eccc()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

# Find the climate_id
stations_search("Brandon A", normals_years = "current")

# Download climate normals 1981-2010
n <- normals_dl(climate_ids = "5010480")
n

# Pull out last frost data
library(tidyr)
f <- unnest(n, frost)
f

# Pull out normals
nm <- unnest(n, normals)
nm

# Download climate normals 1971-2000
n <- normals_dl(climate_ids = "5010480", normals_years = "1971-2000")
n

# Note that some do not have last frost dates
n$frost

# Download multiple stations for 1981-2010,
n <- normals_dl(climate_ids = c("301C3D4", "301FFNJ", "301N49A"))
n

# Note, putting both into the same data set can be done but makes for
# a very unweildly dataset (there is lots of repetition)
nm <- unnest(n, normals)
f <- unnest(n, frost)
both <- dplyr::full_join(nm, f)
both
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{check_eccc}
\alias{check_eccc}
\title{Check access to ECCC}
\usage{
check_eccc()
}
\value{
FALSE if not, TRUE if so
}
\description{
Checks if whether there is internet access, weather data, normals data,
and eccc sites are available and accessible, and whether we're NOT running
on cran
}
\examples{

check_eccc()

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stations.R
\name{stations_dl}
\alias{stations_dl}
\title{Get available stations}
\usage{
stations_dl(skip = NULL, verbose = FALSE, quiet = FALSE)
}
\arguments{
\item{skip}{Numeric. Number of lines to skip at the beginning of the csv. If
NULL, automatically derived.}

\item{verbose}{Logical. Include progress messages}

\item{quiet}{Logical. Suppress all messages (including messages regarding
missing data, etc.)}
}
\description{
This function can be used to download a Station Inventory CSV file from
Environment and Climate Change Canada. This is only necessary if the station
you're interested was only recently added. The 'stations' data set included
in this package contains station data downloaded when the package was last
compiled. This function may take a few minutes to run.
}
\details{
The stations list is downloaded from the url stored in the option
\code{weathercan.urls.stations}. To change this location use
\code{options(weathercan.urls.stations = "your_new_url")}.

The list of which stations have climate normals is downloaded from the url
stored in the option \code{weathercan.urls.stations.normals}. To change this
location use \code{options(weathercan.urls.normals = "your_new_url")}.

Currently there are two sets of climate normals available: 1981-2010 and
1971-2000. Whether a station has climate normals for a given year range is
specified in \code{normals_1981_2010} and \code{normals_1971_2000}, respectively.

The column \code{normals} represents the most current year range of climate
normals (i.e. currently 1981-2010)
}
\examples{
\dontshow{if (check_eccc()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

# Update stations data frame
stations_dl()

# Updated stations data frame is now automatically used
stations_search("Winnipeg")
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/weathercan-pkg.R
\docType{package}
\name{weathercan-package}
\alias{weathercan-package}
\alias{weathercan}
\title{Easy downloading of weather data from Environment and Climate Change Canada}
\description{
\code{weathercan} is an R package for simplifying the downloading of
Historical Climate Data from the Environment and Climate Change Canada (ECCC)
website (\url{https://climate.weather.gc.ca})
}
\details{
Bear in mind that these downloads can be fairly large and performing
repeated, large downloads may use up Environment Canada's bandwidth
unnecessarily. Try to stick to what you need.

There are four main aspects of this package:
\enumerate{
\item Access \strong{stations} lists
\itemize{
\item \code{\link{stations}} (a data frame listing stations)
\item \code{\link[=stations_search]{stations_search()}} identify stations by name or proximity to a
location
\item \code{\link[=stations_dl]{stations_dl()}} re-download/update stations data
}
\item Download \strong{weather} data
}
\itemize{
\item \code{\link[=weather_dl]{weather_dl()}}
}
\enumerate{
\item Merge \strong{weather} data into other data sets through interpolation
over time
}
\itemize{
\item \code{\link[=weather_interp]{weather_interp()}}
}
\enumerate{
\item Download \strong{climate normals} data
}
\itemize{
\item \code{\link[=normals_dl]{normals_dl()}}
}

We also include several practice data sets:
\itemize{
\item \code{\link{finches}}
\item \code{\link{kamloops}}
\item \code{\link{kamloops_day}}
\item \code{\link{pg}}
}

As well as several vignettes:
\itemize{
\item General Usage: \code{vignette("usage")}
\item Merging and Interpolating: \code{vignette("interpolation")}
\item Flags and Codes: \code{vignette("flags")}
\item Weather Data Glossary: \code{vignette("glossary")}
\item Climate Normals Glossary: \code{vignette("glossary_normals")}
}

\href{https://docs.ropensci.org/weathercan/}{Online} we also have some
advanced articles:
\itemize{
\item Using \code{weathercan} with \href{https://www.tidyverse.org/}{tidyverse}
(\href{https://docs.ropensci.org/weathercan/articles/articles/use_with_tidyverse.html}{here})
\item Mapping weather data
(\href{https://docs.ropensci.org/weathercan/articles/articles/mapping.html}{here})
}
}
\references{
Environment and Climate Change Canada: \url{https://www.canada.ca/en/environment-climate-change.html}

Glossary of terms \url{https://climate.weather.gc.ca/glossary_e.html}

ECCC Historical Climate Data: \url{https://climate.weather.gc.ca/}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{glossary_normals}
\alias{glossary_normals}
\title{Glossary of terms for Climate Normals}
\format{
A data frame with 18 rows and 3 variables:
\describe{
\item{ECCC_name}{Original measurement type from ECCC}
\item{weathercan_name}{R-compatible name given when downloaded with the
\code{normals_dl()} function}
\item{description}{Description of the measurement type from ECCC}
}
}
\usage{
glossary_normals
}
\description{
A reference dataset matching information on columns in climate normals data
downloaded using the \code{normals_dl()} function. Indicates the names and
descriptions of different data measurements.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{finches}
\alias{finches}
\title{RFID Data on finch visits to feeders}
\format{
An example dataset of finch RFID data for interpolation:
\describe{
\item{bird_id}{Bird ID number}
\item{time}{Time}
\item{feeder_id}{feeder ID}
\item{species}{Species}
\item{lat}{Latitude of station location in degree decimal format}
\item{lon}{Longitude of station location in degree decimal format}
}
}
\usage{
finches
}
\description{
RFID Data on finch visits to feeders
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{flags}
\alias{flags}
\title{Meaning of coded 'flags'}
\format{
A data frame with 16 rows and 2 variables:
\describe{
\item{code}{Flag code}
\item{meaning}{Explanation of the code}
}
}
\usage{
flags
}
\description{
A reference dataset containing 'flags' matched to their meaning. Data
downloaded using the \code{weather_dl()} function contains columns indicating
'flags' these codes are presented here for interpretation.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stations.R
\name{stations}
\alias{stations}
\title{Access Station data downloaded from Environment and Climate Change Canada}
\format{
A data frame:
\describe{
\item{prov}{Province}
\item{station_name}{Station name}
\item{station_id}{Environment Canada's station ID number. Required for
downloading station data.}
\item{climate_id}{Climate ID number}
\item{WMO_id}{Climate ID number}
\item{TC_id}{Climate ID number}
\item{lat}{Latitude of station location in degree decimal format}
\item{lon}{Longitude of station location in degree decimal format}
\item{elev}{Elevation of station location in metres}
\item{tz}{Local timezone excluding any Daylight Savings}
\item{interval}{Interval of the data measurements ('hour', 'day', 'month')}
\item{start}{Starting year of data record}
\item{end}{Ending year of data record}
\item{normals}{Whether current climate normals are available for that station}
\item{normals_1981_2010}{Whether 1981-2010 climate normals are available for that station}
\item{normals_1971_2000}{Whether 1981-2010 climate normals are available for that station}
}
}
\source{
\url{https://climate.weather.gc.ca/index_e.html}
}
\usage{
stations()
}
\description{
This function access the built-in stations data frame. You can update this
data frame with \code{stations_dl()} which will update the locally stored data.
}
\details{
You can check when this was last updated with \code{stations_meta()}.

A dataset containing station information downloaded from Environment and
Climate Change Canada. Note that a station may have several station IDs,
depending on how the data collection has changed over the years. Station
information can be updated by running \code{stations_dl()}.
}
\examples{
\dontshow{if (check_eccc()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

stations()
stations_meta()

library(dplyr)
filter(stations(), interval == "hour", normals == TRUE, prov == "MB")

\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/interpolate.R
\name{weather_interp}
\alias{weather_interp}
\alias{add_weather}
\title{Interpolate and add weather data to a dataframe}
\usage{
weather_interp(
  data,
  weather,
  cols = "all",
  interval = "hour",
  na_gap = 2,
  quiet = FALSE
)
}
\arguments{
\item{data}{Dataframe. Data with dates or times to which weather data should
be added.}

\item{weather}{Dataframe. Weather data downloaded with \code{\link{weather}}
which should be interpolated and added to \code{data}.}

\item{cols}{Character. Vector containing the weather columns to add or 'all'
for all relevant columns. Note that some measure are omitted because they
cannot be linearly interpolated (e.g., wind direction).}

\item{interval}{What interval is the weather data recorded at? "hour" or
"day".}

\item{na_gap}{How many hours or days (depending on the interval) is it
acceptable to skip over when interpolating over NAs (see details).}

\item{quiet}{Logical. Suppress all messages (including messages regarding
missing data, etc.)}
}
\description{
When data and the weather measurements do not perfectly line up, perform a
linear interpolation between two weather measurements and merge the results
into the provided dataset. Only applies to numerical weather columns (see
\code{weather} for more details).
}
\details{
\strong{Dealing with NA values} If there are NAs in the weather data,
\code{na_gap} can be used to specify a tolerance. For example, a tolerance of
2 with an interval of "hour", means that a two hour gap in data can be
interpolated over (i.e. if you have data for 9AM and 11AM, but not 10AM, the
data between 9AM and 11AM will be interpolated. If, however, you have 9AM and
12PM, but not 10AM or 11AM, no interpolation will happen and data between 9AM
and 12PM will be returned as NA.)
}
\examples{
\dontshow{if (check_eccc()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

# Weather data only
head(kamloops)

# Data about finch observations at RFID feeders in Kamloops, BC
head(finches)

# Match weather to finches
finch_weather <- weather_interp(data = finches, weather = kamloops)
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{kamloops_day}
\alias{kamloops_day}
\title{Daily weather data for Kamloops}
\format{
An example dataset of daily weather data for Kamloops:
\describe{
\item{station_name}{Station name}
\item{station_id}{Environment Canada's station ID number. Required for
downloading station data.}
\item{prov}{Province}
\item{lat}{Latitude of station location in degree decimal format}
\item{lon}{Longitude of station location in degree decimal format}
\item{date}{Date}
\item{year}{Year}
\item{month}{Month}
\item{day}{Day}
\item{cool_deg_days}{Cool degree days}
\item{cool_deg_days_flag}{Cool degree days flag}
\item{dir_max_gust}{Direction of max wind gust}
\item{dir_max_gust_flag}{Direction of max wind gust flag}
\item{heat_deg_days}{Heat degree days}
\item{heat_deg_days_flag}{Heat degree days flag}
\item{max_temp}{Maximum temperature}
\item{max_temp_flag}{Maximum temperature flag}
\item{mean_temp}{Mean temperature}
\item{mean_temp_flag}{Mean temperature flag}
\item{min_temp}{Minimum temperature}
\item{min_temp_flag}{Minimum temperature flag}
\item{snow_grnd}{Snow on the ground (cm)}
\item{snow_grnd_flag}{Snow on the ground flag}
\item{spd_max_gust}{Speed of the max gust km/h}
\item{spd_max_gust_flag}{Speed of the max gust flag}
\item{total_precip}{Total precipitation (any form)}
\item{total_precip_flag}{Total precipitation flag}
\item{total_rain}{Total rain (any form)}
\item{total_rain_flag}{Total rain flag}
\item{total_snow}{Total snow (any form)}
\item{total_snow_flag}{Total snow flag}
\item{elev}{Elevation (m)}
\item{climate_id}{Climate identifier}
\item{WMO_id}{World Meteorological Organization Identifier}
\item{TC_id}{Transport Canada Identifier}
}
}
\source{
\url{https://climate.weather.gc.ca/index_e.html}
}
\usage{
kamloops_day
}
\description{
Downloaded with \code{\link{weather}()}. Terms are more thoroughly defined
here \url{https://climate.weather.gc.ca/glossary_e.html}
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stations.R
\name{stations_search}
\alias{stations_search}
\title{Search for stations by name or location}
\usage{
stations_search(
  name = NULL,
  coords = NULL,
  dist = 10,
  interval = c("hour", "day", "month"),
  normals_years = NULL,
  normals_only = NULL,
  stn = NULL,
  starts_latest = NULL,
  ends_earliest = NULL,
  verbose = FALSE,
  quiet = FALSE
)
}
\arguments{
\item{name}{Character. A vector of length 1 or more with text against which
to match. Will match station names that contain all components of
\code{name}, but they can be in different orders and separated by other
text.}

\item{coords}{Numeric. A vector of length 2 with latitude and longitude of a
place to match against. Overrides \code{lat} and \code{lon} if also
provided.}

\item{dist}{Numeric. Match all stations within this many kilometres of the
\code{coords}.}

\item{interval}{Character. Return only stations with data at these intervals.
Must be any of "hour", "day", "month".}

\item{normals_years}{Character. One of \code{NULL} (default), \code{current},
\code{1981-2010}, or \code{1971-2000}. \code{current} returns only stations from most
recent normals year range. Default \code{NULL} does not filter by climate
normals. Specific year ranges return stations with normals in that period.
See Details for more specifics.}

\item{normals_only}{DEPRECATED. Logical. Return only stations with climate
normals?}

\item{stn}{DEFUNCT. Now use \code{stations_dl()} to update internal data and
\code{stations_meta()} to check the date it was last updated.}

\item{starts_latest}{Numeric. Restrict results to stations with data
collection beginning in or before the specified year.}

\item{ends_earliest}{Numeric. Restrict results to stations with data
collection ending in or after the specified year.}

\item{verbose}{Logical. Include progress messages}

\item{quiet}{Logical. Suppress all messages (including messages regarding
missing data, etc.)}
}
\value{
Returns a subset of the stations data frame which match the search
parameters. If the search was by location, an extra column 'distance' shows
the distance in kilometres from the location to the station. If no stations
are found withing \code{dist}, the closest 10 stations are returned.
}
\description{
Returns stations that match the name provided OR which are within \code{dist}
km of the location provided. This is designed to provide the user with
information with which to decide which station to then get weather data from.
}
\details{
To search by coordinates, users must make sure they have the
\href{https://cran.r-project.org/package=sp}{sp} package installed.

The \code{current}, most recent, climate normals year range is \code{1981-2010}.
}
\examples{
\dontshow{if (check_eccc()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}

stations_search(name = "Kamloops")
stations_search(name = "Kamloops", interval = "hour")

stations_search(name = "Ottawa", starts_latest = 1950, ends_earliest = 2010)

stations_search(name = "Ottawa", normals_years = "current")   # 1981-2010
stations_search(name = "Ottawa", normals_years = "1981-2010") # Same as above
stations_search(name = "Ottawa", normals_years = "1971-2000") # 1971-2010

if(requireNamespace("sp")) {
  stations_search(coords = c(53.915495, -122.739379))
}
\dontshow{\}) # examplesIf}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/stations.R
\name{stations_meta}
\alias{stations_meta}
\title{Show stations list meta data}
\usage{
stations_meta()
}
\description{
Date of ECCC update and date downloaded via weathercan.
}
\examples{
\dontshow{if (check_eccc()) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
stations_meta()
\dontshow{\}) # examplesIf}
}
