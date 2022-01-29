# tidyhydat use cases

```
@article{beaton2019identifying,
  title={Identifying historic river ice breakup timing using MODIS and Google Earth Engine in support of operational flood monitoring in Northern Ontario},
  author={Beaton, A and Whaley, R and Corston, K and Kenny, F},
  journal={Remote Sensing of Environment},
  volume={224},
  pages={352--364},
  year={2019},
  publisher={Elsevier}
}
```

```
@article{moore2020detecting,
  title={Detecting the Effects of Sustained Glacier Wastage on Streamflow in Variably Glacierized Catchments},
  author={Moore, Robert Daniel and Pelto, Ben and Menounos, Brian and Hutchinson, David},
  journal={Frontiers in Earth Science},
  volume={8},
  pages={136},
  year={2020},
  publisher={Frontiers}
}
```
<!-- README.md is generated from README.Rmd. Please edit that file -->

# tidyhydat <img src="man/figures/logo.png" align="right" />

<!-- badges: start -->

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![Coverage
status](https://codecov.io/gh/ropensci/tidyhydat/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tidyhydat?branch=master)
[![R build
status](https://github.com/ropensci/tidyhydat/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/tidyhydat/actions)

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/tidyhydat)](https://cran.r-project.org/package=tidyhydat)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/tidyhydat?color=brightgreen)](https://CRAN.R-project.org/package=tidyhydat)
[![cran
checks](https://cranchecks.info/badges/worst/tidyhydat)](https://cran.r-project.org/web/checks/check_results_tidyhydat.html)
[![r-universe](https://ropensci.r-universe.dev/badges/tidyhydat)](https://ropensci.r-universe.dev/ui#builds)

[![](http://badges.ropensci.org/152_status.svg)](https://github.com/ropensci/software-review/issues/152)
[![DOI](http://joss.theoj.org/papers/10.21105/joss.00511/status.svg)](https://doi.org/10.21105/joss.00511)
[![DOI](https://zenodo.org/badge/100978874.svg)](https://zenodo.org/badge/latestdoi/100978874)
<!-- badges: end -->

## Project Status

This package is maintained by the Data Science Partnerships Program in
the [British Columbia Ministry of Citizens’
Services](https://www2.gov.bc.ca/gov/content/governments/organizational-structure/ministries-organizations/ministries/citizens-services).

## What does `tidyhydat` do?

-   Provides functions (`hy_*`) that access hydrometric data from the
    HYDAT database, a national archive of Canadian hydrometric data and
    return tidy data.
-   Provides functions (`realtime_*`) that access Environment and
    Climate Change Canada’s real-time hydrometric data source.
-   Provides functions (`search_*`) that can search through the
    approximately 7000 stations in the database and aid in generating
    station vectors
-   Keep functions as simple as possible. For example, for daily flows,
    the `hy_daily_flows()` function queries the database, *tidies* the
    data and returns a [tibble](https://tibble.tidyverse.org/) of daily
    flows.

## Installation

You can install `tidyhydat` from CRAN:

    install.packages("tidyhydat")

To install the development version of the `tidyhydat` package, you can
install directly from the rOpenSci development server:

    install.packages("tidyhydat", repos = "https://dev.ropensci.org")

## Usage

More documentation on `tidyhydat` can found at the rOpenSci doc page:
<https://docs.ropensci.org/tidyhydat/>

When you install `tidyhydat`, several other packages will be installed
as well. One of those packages, `dplyr`, is useful for data
manipulations and is used regularly here. To use actually use `dplyr` in
a session you must explicitly load it. A helpful `dplyr` tutorial can be
found
[here](https://cran.r-project.org/package=dplyr/vignettes/dplyr.html).

    library(tidyhydat)
    library(dplyr)

### HYDAT download

To use many of the functions in the `tidyhydat` package you will need to
download a version of the HYDAT database, Environment and Climate Change
Canada’s database of historical hydrometric data then tell R where to
find the database. Conveniently `tidyhydat` does all this for you via:

    download_hydat()

This downloads (with your permission) the most recent version of HYDAT
and then saves it in a location on your computer where `tidyhydat`’s
function will look for it. Do be patient though as this can take a long
time! To see where HYDAT was saved you can run `hy_default_db()`. Now
that you have HYDAT downloaded and ready to go, you are all set to begin
looking at Canadian hydrometric data.

### Real-time

To download real-time data using the datamart we can use approximately
the same conventions discussed above. Using `realtime_dd()` we can
easily select specific stations by supplying a station of interest:

    realtime_dd(station_number = "08LG006")
    #>   Queried on: 2021-12-16 18:51:26 (UTC)
    #>   Date range: 2021-11-15 to 2021-11-15 
    #> # A tibble: 248 x 8
    #>    STATION_NUMBER PROV_TERR_STATE_LOC Date                Parameter Value Grade
    #>    <chr>          <chr>               <dttm>              <chr>     <dbl> <chr>
    #>  1 08LG006        BC                  2021-11-15 08:00:00 Flow        162 <NA> 
    #>  2 08LG006        BC                  2021-11-15 08:05:00 Flow        164 <NA> 
    #>  3 08LG006        BC                  2021-11-15 08:10:00 Flow        165 <NA> 
    #>  4 08LG006        BC                  2021-11-15 08:15:00 Flow        167 <NA> 
    #>  5 08LG006        BC                  2021-11-15 08:20:00 Flow        169 <NA> 
    #>  6 08LG006        BC                  2021-11-15 08:25:00 Flow        171 <NA> 
    #>  7 08LG006        BC                  2021-11-15 08:30:00 Flow        172 <NA> 
    #>  8 08LG006        BC                  2021-11-15 08:35:00 Flow        174 <NA> 
    #>  9 08LG006        BC                  2021-11-15 08:40:00 Flow        174 <NA> 
    #> 10 08LG006        BC                  2021-11-15 08:45:00 Flow        176 <NA> 
    #> # ... with 238 more rows, and 2 more variables: Symbol <chr>, Code <chr>

### Plotting

Plot methods are also provided to quickly visualize realtime data:

    realtime_ex <- realtime_dd(station_number = "08LG006")

    plot(realtime_ex)

![](man/figures/README-unnamed-chunk-7-1.png)

and also historical data:

    hy_ex <- hy_daily_flows(station_number = "08LA001", start_date = "2013-01-01")

    plot(hy_ex)

![](man/figures/README-unnamed-chunk-8-1.png)

## Getting Help or Reporting an Issue

To report bugs/issues/feature requests, please file an
[issue](https://github.com/ropensci/tidyhydat/issues/).

These are very welcome!

## How to Contribute

If you would like to contribute to the package, please see our
[CONTRIBUTING](https://github.com/ropensci/tidyhydat/blob/master/CONTRIBUTING.md)
guidelines.

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/tidyhydat/blob/master/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.

## Citation

Get citation information for `tidyhydat` in R by running:


    Albers S (2017). "tidyhydat: Extract and Tidy Canadian Hydrometric
    Data." _The Journal of Open Source Software_, *2*(20). doi:
    10.21105/joss.00511 (URL: https://doi.org/10.21105/joss.00511), <URL:
    http://dx.doi.org/10.21105/joss.00511>.

    A BibTeX entry for LaTeX users is

      @Article{,
        title = {tidyhydat: Extract and Tidy Canadian Hydrometric Data},
        author = {Sam Albers},
        doi = {10.21105/joss.00511},
        url = {http://dx.doi.org/10.21105/joss.00511},
        year = {2017},
        publisher = {The Open Journal},
        volume = {2},
        number = {20},
        journal = {The Journal of Open Source Software},
      }

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

## License

Copyright 2017 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the “License”); you may
not use this file except in compliance with the License. You may obtain
a copy of the License at

<https://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an “AS IS” BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.
---
title: "Release Steps"
output: github_document
---

## Update `allstations` data and documentation
```
source("data-raw/HYDAT_internal_data/process_internal_data.R")
```

## Check if version is appropriate
http://shiny.andyteucher.ca/shinyapps/rver-deps/

## Build and check within `R/devtools`
```
devtools::check_win_devel()
devtools::check_win_release()
devtools::check() ## build locally
```

## Build and check on rhub
```
library(rhub)

check_with_rdevel()
check_for_cran()
check_on_windows()
```

## Run this in the console
```
R CMD build tidyhydat
R CMD check tidyhydat_0.5.1.tar.gz --as-cran ## or whatever the package name is
```

## Documentation
- Update NEWS
- Update cran-comments

## Actually release it
```
devtools::release()
```

## Once it is release create signed release on github
```
git tag -s [version] -m "[version]"
git push --tags
```

```
# Copyright 2018 Province of British Columbia
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and limitations under the License.
```
# tidyhydat (development version)

### MINOR IMPROVEMENT
* `download_hydat()` now has an `ask` parameter that can be used to bypass the keypress confirmation when downloading the HYDAT database (@rchlumsk, #165). 

# tidyhydat 0.5.4
- When add a local timezone column, use the most common timezone in the data rather than the first one. This just seems more likely to be useful to users
- Add more documentation to `realtime_add_local_datetime` to make how timezones are dealt with clearer (#157)
- Expose the query time for realtime functions as an attribute (#160)
- Add Government of Canada as data contributor (#156)

# tidyhydat 0.5.3
- Allow pkg to loaded without internet and rather just issue an message when it is not. (#149)
- Added `add = TRUE` to all `on.exit` call so add not to overwrite previous call (#151)
- Remove redundant and ill-advised  `closeAllConnections` (#153)
- Update internal data
- Convert most of the docs to markdown (#121)

# tidyhydat 0.5.2
- add internal function `hy_check` to verify that HYDAT contains all the right tables and that those tables contain data. 

# tidyhydat 0.5.1
- Replace `class(x) ==` with `inherits`
- Fix bug and added corresponding tests where a request for multiple stations to `realtime_dd` would fail if any data was missing
- Update internal data
- Fix all non-secure or borken links

# tidyhydat 0.5.0

### MINOR FIXES
- Revise multi prov test to realtime because of network load and prone to intermittent failure
- Adding rOpenSci doc site to DESCRIPTION
- Fix character NA's in `hy_stations` (#125)
- Allow downloading HYDAT to alternative locations (#129)
- Provide better documentation of how change default db location in `hy_set_default_db()`

# tidyhydat 0.4.0


### IMPROVEMENTS
* All functions now return either "hy" or "realtime" class with associated print and plot methods (#119)
* prov_terr_state_loc now accepts a "CA" value to specify only stations located in Canada (#112)
* functions that access internet resources now fail with an informative error message (#116)
* tests that require internet resources are skipped when internet is down 
* Add small join example to calculate runoff to introduction vignette (#120)

### BUG FIXES
* `pull_station_number` now only returns unique values (#109)
* Adding a offset column that reflects OlsonNames() and is thus DST independent (#110)
* Caught all `R_CHECK_LENGTH_1_CONDITION` instances

# tidyhydat 0.3.5

### IMPROVEMENTS
* New function: `realtime_add_local_datetime()` adds a local datetime column to `realtime_dd()` tibble (#64)
* New function: `pull_station_number()` wraps `pull(STATION_NUMBER)` for convenience

### MINOR BREAKING CHANGES
* In effort to standardize, the case of column names for some rarely used function outputs were changed to reflect more commonly used function outputs. This may impact some workflows where columns are referenced by names (#99).   

### BUG FIXES
* Functions that have a `start_date` and `end_date` actually work with said argument (#98)
* `hy_annual_instant_peaks()` now parses the date correctly into UTC and includes a datetime and time zone column.  (#64)
* `hy_stn_data_range()` now returns actual `NA`'s rather than string NA's (#97)

### MINOR IMPROVEMENT
* `download_hydat()` now returns an informative error if the download fails due to proxy-related connection issues (@rywhale, #101). 

# tidyhydat 0.3.4

### IMPROVEMENT
* Added rlang as a dependency and applied tidyeval idiom to more safety control variable environments
* 15% speed improvement in `realtime_dd` by eliminating loop (#91)
* 40% speed improvement when querying full provinces (#89)
* reorganized file naming so that helper functions are placed in utils-* files

### BUG FIXES
* Fixed `hy_monthly_flows` and `hy_monthly_levels` date issue (#24)

### MINOR IMPROVEMENT
* realtime tidying now not duplicated and is handled by a function
* simplified `tidyhydat:::station_choice` and added more unit testing
* no longer outputting a message when `station_number = "ALL"`.
* Exporting pipe (`%>%`)

# tidyhydat 0.3.3

### NEW FEATURES
  * Open a connection to the HYDAT database directly using `hy_src()` for advanced functionality (PR#77).
  * New vignette outlining `hy_src()` (PR#77)
  * Add some tools to improve the usability of the test database (PR#77).
  * `download_hydat()` now uses `httr::GET()`

### MINOR IMPROVEMENTS
  * Better downloading messages
  
### BUG FIXES
  * Fixed package startup message so it can be supressed. (#79)
  * Fixed bug that resulted in `download_hydat` choice wasn't respected.
  * `onAttach()` now checks 115 days after last HYDAT release to prevent slow package load times if HYDAT is longer than 3 months between RELEASES.
  * Fixed margin error in `hy_plot()`
  * Fixed a bug in `realtime_plot()` that prevented a lake level station from being called
  * Fixed a bug in `hy_daily()` that threw an error when only a level station was called
  * Added new tests for `hy_daily()` and `realtime_plot()`
  * Added `HYD_STATUS` and `REAL_TIME` columns to `allstations`. 


# tidyhydat 0.3.2

### NEW FEATURES
  * New `hy_daily()` function which combines all daily data into one dataframe.
  * Add a quick base R plotting feature for quick visualization of realtime and historical data.
  * Add `realtime_daily_mean` function that quickly converts higher resolution data into daily means.
  * New vignette outlining some example usage.
  
### BUG FIXES
  * Fixed bug in `download_hydat()` that create a path that wasn't OS-independent.
  * Fixed a bug on `download_hydat()` where by sometimes R had trouble overwriting an existing version of existing database. Now the old database is simply deleted before the new one is downloaded.
  * `hy_annual_instant_peaks()` now returns a date object with HOUR, MINUTE and TIME_ZONE returned as separed columns. (#10)
  * All variable values of LEVEL and FLOW have been changed to Level and Flow to match the output of `hy_data_types`. (#60)
  * Tidier and coloured error messages throughout.
  * Review field incorrectly specified the rOpenSci review page. Removed the link from the DESCRIPTION.
  


# tidyhydat 0.3.1

### NEW FEATURES

  * When package is loaded, tidyhydat checks to see if HYDAT is even present
  * When package is loaded, it now tests to see if their a new version of HYDAT if the current date is greater than 3 months after the last release date of HYDAT. 
  * Prep for CRAN release
  * Starting to use raw SQL for table queries
  * Removing 2nd vignette from build. Still available on github

# tidyhydat 0.3.0 

### NEW FEATURES

  * New NEWS template!
  * Moved `station_number` to first argument to facilitate piped analysis (#54)
  * `search_stn_name` and `search_stn_number` now query both realtime and historical data sources and have tests for a more complete list (#56)
  * With credential stored in .Renviron file, `ws_token` can successfully be called by `ws_token()`.
  * `.onAttach()` checks if HYDAT is downloaded on package load.

### MINOR IMPROVEMENTS
  * Significant function and argument name changes (see below)
  * Adding `rappdirs` to imports and using to generate download path for `download_hydat()` (#44)
  * Adding `rappdirs` so that all the hy_* functions access hydat from `rappdirs::user_data_dir()` via `hy_dir()` (#44)
  * Revised and cleaned up documentation including two vignettes (#48)
  * `FULL MONTH` evaluate to a logic (#51)
  * All download tests are skipped on cran (#53)  
  * Removed time limit for `download_realtime_ws()` with some documentation on actual limits. [(3234c22)](https://github.com/ropensci/tidyhydat/commit/3234c2246c97fed5860e8dfb9adc3d6f0aa503fe)


### BUG FIXES

  * Add informative error message for a single missing station input (#38)
  * No longer trying to build .Rd file for `.onload` (#47)
  * Fixed `SED_MONTHLY_LOADS` (#51)
  

### FUNCTION NAME CHANGES (#45)
  * hy_agency_list <- AGENCY_LIST
  * hy_annual_instant_peaks <- ANNUAL_INSTANT_PEAKS
  * hy_annual_stats <- ANNUAL_STATISTICS
  * hy_daily_flows <- DLY_FLOWS
  * hy_daily_levels <- DLY_LEVELS
  * hy_monthly_flows <- MONTHLY_FLOWS
  * hy_monthly_levels <- MONTHLY_LEVELS
  * hy_sed_daily_loads <- SED_DLY_LOADS
  * hy_sed_daily_suscon <- SED_DLY_SUSCON
  * hy_sed_monthly_loads <- SED_MONTHLY_LOADS
  * hy_sed_monthly_suscon <- SED_MONTHLY_SUSCON
  * hy_sed_samples <- SED_SAMPLES
  * hy_sed_samples_psd <- SED_SAMPLES_PSD
  * hy_stations <- STATIONS
  * hy_stn_remarks <- STN_REMARKS
  * hy_stn_datum_conv <- STN_DATUM_CONVERSION
  * hy_stn_datum_unrelated <- STN_DATUM_UNRELATED
  * hy_stn_data_range <- STN_DATA_RANGE
  * hy_stn_data_coll <- STN_DATA_COLLECTION
  * hy_stn_op_schedule <- STN_OPERATION_SCHEDULE
  * hy_stn_regulation <- STN_REGULATION
  * hy_agency_list <- AGENCY_LIST
  * hy_reg_office_list <- REGIONAL_OFFICE_LIST
  * hy_datum_list <- DATUM_LIST
  * hy_version <- VERSION
  * realtime_dd <- download_realtime_dd
  * realtime_stations <- realtime_network_meta
  * search_stn_name <- search_name
  * search_stn_number <- search_number
  
### ARGUMENT NAME CHANGES (#45)
  * station_number <- STATION_NUMBER
  * prov_terr_state_loc <- PROV_TERR_STATE_LOC



# tidyhydat 0.2.9

* Explicitly state in docs that time is in UTC (#32)
* Added test for realtime_network_meta and moved to httr to download.
* download functions all use httr now
* removed need for almost all @import statement by referencing them all directly (#34)
* Fixed error message when directly calling some tidyhydat function using :: (#31)
* To reduce overhead, `output_symbol` has been added as an argument so code can be produced if desired (#33)

# tidyhydat 0.2.8

* Added examples to every function
* Completed test suite including `download_realtime_ws` (#27)
* Fixed bugs in several `STN_*` functions
* Added `STN_DATUM_RELATED`
* Updated documentation

# tidyhydat 0.2.7

* Updated documentation
* Updated README
* Created a small database so that unit testing occurs remotely (#1)
* Fixed `STN_DATA_RANGE` bug (#26)

# tidyhydat 0.2.6

* using `styler` package to format code to tidyverse style guide
* added `PROV_TERR_STATE_LOC` to `allstations`
* added `search_number` function
* added `MONTHLY` functions
* created function families
* added `on.exit()` to internal code; a better way to disconnect
* Updated documentation

# tidyhydat 0.2.5

* fixed minor bug in download_realtime_ws so that better error message is outputted when no data is returned

# tidyhydat 0.2.4

* download_realtime_dd can now accept stations from multiple provinces or simply select multiple provinces
* better error messages for get_ws_token and download_realtime_ws
* All functions that previously accepted STATION_NUMBER == "ALL" now throw an error. 
* Added function to download hydat

# tidyhydat 0.2.3

* Remove significant redundancy in station selecting mechanism
* Added package startup message when HYDAT is out of date  
* Add internal allstations data
* Added all the tables as functions or data from HYDAT
* Made missing station ouput truncated at 10 missing stations

# tidyhydat 0.2.2

* Adding several new tables
* removed need for both prov and stn args
* reduced some repetition in code

# tidyhydat 0.2.1

* added STN_REGULATION
* tidied ANNUAL_STATISTICS
* added a series of lookup tables (DATUM_LIST, AGENCY_LIST, REGIONAL_OFFICE_LIST)
* cleared up output of STATIONS

# tidyhydat 0.2.0

* standardize hydat outputs to consistent tibble structure
* Adding search_name function
* final names for download functions
* functions output an information message about stations retrieved

# tidyhydat 0.1.1

*Renamed real-time function as download_realtime and download_realtime2
*Added more units tests
*Wrote vignette for package utilization
*Brought all data closer to a "tidy" state

# tidyhydat 0.1.0

*Added ability for STATIONS to retrieve ALL stations in the HYDAT database
*Added ability for STATIONS to retrieve ALL stations in the HYDAT database
*Standardize documentation; remove hydat_path default
*Better error handling for download_realtime
*Update documentation
*Adding param_id data, data-raw and documentation
*Dates filter to ANNUAL_STATISTICS and DLY_FLOWS; func and docs
*DLY_LEveLS function and docs
*download_ws and get_ws_token function and docs
*UPDATE README

# tidyhydat 0.0.4

*Added ability for STATIONS to retrieve ALL stations in the HYDAT database
*Added ability for STATIONS to retrieve ALL stations in the HYDAT database
*Standardize documentation; remove hydat_path default
*Better error handling for download_realtime
*Update documentation
*Adding param_id data, data-raw and documentation
*Dates filter to ANNUAL_STATISTICS and DLY_FLOWS; func and docs
*DLY_LEveLS function and docs
*download_ws and get_ws_token function and docs
*UPDATE README

# tidyhydat 0.0.3

*fixed db connection problem; more clear documentation
*better error handling; more complete realtime documentation
*harmonized README with standardized arguments

# tidyhydat 0.0.2

*Added example analysis to README
*Added devex badge; license to all header; import whole readr package
*Able to take other protidyhydat inces than BC now
*Update documentation; README

# tidyhydat 0.0.1

*Initial package commit
*Add license and include bcgotidyhydat  files in RBuildIgnore
*Two base working function; package level R file and associated documentation
*Only importing functions used in the function
*Update README with example
*Added download_ functions
*Added ANNUAL_STATISTICS query/table and docs
*Updated docs and made DLY_FLOWS more rigorous
# Contributor Code of Conduct

As contributors and maintainers of this project, and in the interest of
fostering an open and welcoming community, we pledge to respect all people who
contribute through reporting issues, posting feature requests, updating
documentation, submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free
experience for everyone, regardless of level of experience, gender, gender
identity and expression, sexual orientation, disability, personal appearance,
body size, race, ethnicity, age, religion, or nationality.

Examples of unacceptable behavior by participants include:

* The use of sexualized language or imagery
* Personal attacks
* Trolling or insulting/derogatory comments
* Public or private harassment
* Publishing other's private information, such as physical or electronic
  addresses, without explicit permission
* Other unethical or unprofessional conduct

Project maintainers have the right and responsibility to remove, edit, or
reject comments, commits, code, wiki edits, issues, and other contributions
that are not aligned to this Code of Conduct, or to ban temporarily or
permanently any contributor for other behaviors that they deem inappropriate,
threatening, offensive, or harmful.

By adopting this Code of Conduct, project maintainers commit themselves to
fairly and consistently applying these principles to every aspect of managing
this project. Project maintainers who do not follow or enforce the Code of
Conduct may be permanently removed from the project team.

This Code of Conduct applies both within project spaces and in public spaces
when an individual is representing the project or its community.

Instances of abusive, harassing, or otherwise unacceptable behavior may be
reported by contacting a project maintainer at sam.albers@gov.bc.ca. All complaints will be reviewed and investigated 
and will result in a response that is deemed necessary and appropriate to the 
circumstances. Maintainers are obligated to maintain confidentiality with regard 
to the reporter of an incident.


This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.3.0, available at
[http://contributor-covenant.org/version/1/3/0/][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/3/0/
## How to contribute
Government employees, public and members of the private sector are encouraged to contribute to the repository by **forking and submitting a pull request**. 

(If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git) and  check out a more detailed guide to [pull requests](https://help.github.com/articles/using-pull-requests/).)

Pull requests will be evaluated by the repository guardians on a schedule and if deemed beneficial will be committed to the master.

All contributors retain the original copyright to their stuff, but by contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users **under the terms of the license under which this project is distributed.**

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/tidyhydat/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/ropensci/tidyhydat.git`
* Make sure to track progress upstream (i.e., on our version of `tidyhydat` at `ropensci/tidyhydat`) by doing `git remote add upstream https://github.com/ropensci/tidyhydat.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/tidyhydat`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Prefer to Email? Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

### Thanks for contributing!
tidyhydat 0.5.4
=========================
- When add a local timezone column, use the most common timezone in the data rather than the first one. This just seems more likely to be useful to users
- Add more documentation to `realtime_add_local_datetime` to make how timezones are dealt with clearer (#157)
- Expose the query time for realtime functions as an attribute (#160)

## Test environments
* win-builder (via `devtools::check_win_devel()` and `devtools::check_win_release()`)
* local Windows 10, R 4.1.0 (via R CMD check --as-cran)
* ubuntu-20.04, r: 'release' (github actions)
* ubuntu-20.04, r: 'devel' (github actions)
* macOS,        r: 'release' (github actions)
* windows,      r: 'release' (github actions)
* Fedora Linux, R-devel, clang, gfortran - r-hub
* Debian Linux, R-release, GCC (debian-gcc-release) - r-hub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit - r-hub


tidyhydat 0.5.3
=========================
## NEWS
- Allow pkg to loaded without internet and rather just issue an message when it is not. (#149)
- Added `add = TRUE` to all `on.exit` call so add not to overwrite previous call (#151)
- Remove redundant and ill-advised  `closeAllConnections` (#153)
- Update internal data
- Convert most of the docs to markdown (#121)

## Test environments
* win-builder (via `devtools::check_win_devel()` and `devtools::check_win_release()`)
* local Windows 10, R 4.0.5 (via R CMD check --as-cran)
* ubuntu-20.04, r: 'release' (github actions)
* ubuntu-20.04, r: 'devel' (github actions)
* macOS,        r: 'release' (github actions)
* windows,      r: 'release' (github actions)
* Fedora Linux, R-devel, clang, gfortran - r-hub
* Debian Linux, R-release, GCC (debian-gcc-release) - r-hub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit - r-hub



tidyhydat 0.5.2
=========================
## NEWS
- add internal function `hy_check` to verify that HYDAT contains all the right tables and that those tables contain data. 

## Test environments
* win-builder (via `devtools::check_win_devel()` and `devtools::check_win_release()`)
* local Windows 10, R 4.0.3 (via R CMD check --as-cran)
* ubuntu-20.04, r: 'release' (github actions)
* ubuntu-20.04, r: 'devel' (github actions)
* macOS,        r: 'release' (github actions)
* windows,      r: 'release' (github actions)
* Fedora Linux, R-devel, clang, gfortran - r-hub
* Debian Linux, R-release, GCC (debian-gcc-release) - r-hub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit - r-hub



tidyhydat 0.5.1
=========================
## Re-submission note
- Fix all non-secure or broken links

## NEWS
- Replace `class(x) ==` with `inherits`
- Fix bug and added corresponding tests where a request for multiple stations to `realtime_dd` would fail if any data was missing
- Update internal data

## Test environments
* win-builder (via `devtools::check_win_devel()` and `devtools::check_win_release()`)
* local Windows 10, R 4.0.2 (via R CMD check --as-cran)
* ubuntu-20.04, r: 'release' (github actions)
* ubuntu-20.04, r: 'devel' (github actions)
* macOS,        r: 'release' (github actions)
* windows,      r: 'release' (github actions)
* Fedora Linux, R-devel, clang, gfortran - r-hub
* Debian Linux, R-release, GCC (debian-gcc-release) - r-hub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit - r-hub


tidyhydat 0.5.0
=========================
## Re-submission note
* Vignette linked changed to canonical form. 
* GitHub actions badge link now points only to GitHub repo.
* `CONTRIBUTING.md` and `CODE_OF_CONDUCT.md` now link directly to GitHub repo.
* CRAN checks link changed to `CRAN.R-project.org` link
* Add title to README and fixed bug in realtime `plot` method and reduce size of README as per Uwe comments

## Test environments
* win-builder (via `devtools::check_win_devel()` and `devtools::check_win_release()`)
* local Windows 10, R 3.6.1 (via R CMD check --as-cran)
* ubuntu, R 3.6.1 (travis-ci) (release)
* ubuntu, (travis-ci) (devel)
* ubuntu-16.04, r: '3.3' (github actions)
* ubuntu-16.04, r: '3.4' (github actions)
* ubuntu-16.04, r: '3.5' (github actions)
* ubuntu-16.04, r: '3.6' (github actions)
* Fedora Linux, R-devel, clang, gfortran - r-hub
* Debian Linux, R-release, GCC (debian-gcc-release) - r-hub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit - r-hub
* macOS 10.11 El Capitan, R-release (experimental) - r-hub


tidhydat 0.4.0
=========================
## Test environments
* win-builder (via `devtools::check_win_devel()` and `devtools::check_win_release()`)
* local Windows 10, R 3.5.3 (via R CMD check --as-cran)
* ubuntu, R 3.5.3 (travis-ci) (release)
* ubuntu, R 3.5.3 (travis-ci) (devel)
* Fedora Linux, R-devel, clang, gfortran - r-hub
* Debian Linux, R-release, GCC (debian-gcc-release) - r-hub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit - r-hub
* macOS 10.11 El Capitan, R-release (experimental) - r-hub


tidyhydat 0.3.5
=========================
### IMPROVEMENTS
* New function: `realtime_add_local_datetime()` adds a local datetime column to `realtime_dd()` tibble (#64)
* New function: `pull_station_number()` wraps `pull(STATION_NUMBER)` for convenience

### MINOR BREAKING CHANGES
* In effort to standardize, the case of column names for some rarely used function outputs were changed to reflect more commonly used function outputs. This may impact some workflows where columns are referenced by names (#99).   

### BUG FIXES
* Functions that have a `start_date` and `end_date` actually work with said argument (#98)
* `hy_annual_instant_peaks()` now parses the date correctly into UTC and includes a datetime and time zone column.  (#64)
* `hy_stn_data_range()` now returns actual `NA`'s rather than string NA's (#97)

### MINOR IMPROVEMENT
* `download_hydat()` now returns an informative error if the download fails due to proxy-related connection issues (@rywhale, #101). 

## Test environments
* win-builder (via `devtools::build_win()`)
* local Windows 10, R 3.4.3 (via R CMD check --as-cran)
* ubuntu, R 3.4.3 (travis-ci) (release)
* Debian Linux, R-release, GCC (debian-gcc-release) - r-hub
* Windows Server 2008 R2 SP1, R-devel, 32/64 bit - r-hub
* macOS 10.11 El Capitan, R-release (experimental) - r-hub
* macOS 10.9 Mavericks, R-oldrel (experimental) (macos-mavericks-oldrel) - r-hub
 
## R CMD check results

* No warnings
* No notes
* No errors



## Downstream dependencies

There are currently no downstream dependencies.
## How to contribute
Government employees, public and members of the private sector are encouraged to contribute to the repository by **forking and submitting a pull request**. 

(If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git) and  check out a more detailed guide to [pull requests](https://help.github.com/articles/using-pull-requests/).)

Pull requests will be evaluated by the repository guardians on a schedule and if deemed beneficial will be committed to the master.

All contributors retain the original copyright to their stuff, but by contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users **under the terms of the license under which this project is distributed.**

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/tidyhydat/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/ropensci/tidyhydat.git`
* Make sure to track progress upstream (i.e., on our version of `tidyhydat` at `ropensci/tidyhydat`) by doing `git remote add upstream https://github.com/ropensci/tidyhydat.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/tidyhydat`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Prefer to Email? Email the person listed as maintainer in the `DESCRIPTION` file of this repo.

### Thanks for contributing!
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

<!--- Did you remember to include tests? Unless you're just changing
grammar, please include new tests for your change -->
---
name: Bug report
about: Create a report to help us improve
title: ''
labels: ''
assignees: ''

---

**Describe the bug**
A clear and concise description of what the bug is.

**To Reproduce**
Steps to reproduce the behavior:
1. Go to '...'
2. Click on '....'
3. Scroll down to '....'
4. See error

**Expected behavior**
A clear and concise description of what you expected to happen.

**Screenshots**
If applicable, add screenshots to help explain your problem.

**Desktop (please complete the following information):**
 - OS: [e.g. iOS]
 - Browser [e.g. chrome, safari]
 - Version [e.g. 22]

**Smartphone (please complete the following information):**
 - Device: [e.g. iPhone6]
 - OS: [e.g. iOS8.1]
 - Browser [e.g. stock browser, safari]
 - Version [e.g. 22]

**Additional context**
Add any other context about the problem here.
---
name: Feature request
about: Suggest an idea for this project
title: ''
labels: ''
assignees: ''

---

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want to happen.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
## tidyhydat internal data

- The file `inst/test_db/tinyhydat.sqlite3` is the first 100 rows (or less) of each table in the hydat database. This is used exclusively for testing purposes and is generated by `tinyhydat_proc.R`

- All other internal data (`allstations`,`hy_data_types`,`hy_data_symbols`) are drawn from HYDAT and are provided under the ([Open Government License - Canada](http://open.canada.ca/en/open-government-licence-canada)). 
---
title: "tidyhydat: Extract and Tidy Canadian Hydrometric Data"
authors:
- affiliation: 1
  name: Sam J. Albers
  orcid: 0000-0002-9270-7884
date: "2017-12-14"
output:
  html_document:
    keep_md: yes
bibliography: paper.bib
tags:
- R
- tidy data
- hydrology
- Canada
affiliations:
- index: 1
  name: Hydrology and Hydrometric Programs, Ministry of Environment and Climate Change Strategy, British Columbia Provincial Government
---

> Tidy datasets are all alike but every messy dataset is messy in its own way - @wickham2014tidy

# Introduction
Environment and Climate Change Canada (ECCC) through the Water Survey of Canada (WSC) maintains several national hydrometric data sources. These data are partially funded by provincial partners and constitute the main data products of a national integrated hydrometric network. Historical data are stored in the [HYDAT database](http://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/). HYDAT is the Canadian national Water Data Archive, published quarterly by the Government of Canada's Department of Environment and Climate Change. It is a relational database that contains daily, monthly and annual data on water flow, water levels and sediment.

Real-time data are provided by ECCC over the web. Files are updated to a [datamart](http://dd.weather.gc.ca/hydrometric/) on an hourly basis though the lag between actual hydrometric measurement and the availability of hydrometric data is approximately 2.5 hours. The objective of this document is to outline the usage of `tidyhydat` [@alberstidyhydat], an R package that accesses these hydrometric data sources and *tidies* them. `tidyhydat` is part of the [rOpenSci](https://ropensci.org/packages/) suite of packages and resides at  https://github.com/ropensci/tidyhydat. The objective of `tidyhydat` is to provide a standard method of accessing ECCC data sources using a consistent and easy to use interface that employs tidy data principles developed by @wickham2014tidy within the R project [@RCore]. 

## Why use R in hydrology?
There are many statistical computing projects that offer great functionality for users. For `tidyhydat` I have chosen to use R. R is a mature open-source project that provides significant potential for advanced modelling, visualization and data manipulation. For hydrologists considering data analysis tools there are several commonly cited reasons to use R:

- R is and always will be free to use and modify.
- R is easily extensible and comprehensive. It is complimented by a rich suite of packages that implement a vast array of classical and modern statistical methods, exceptionally high-quality graphing capabilities and powerful data manipulation tools to handle a wide variety of data formats.
- R facilitates the scientific method by allowing for a fully reproducible data workflow that can be repeated by others when code is shared.  
- R has a friendly community which is an important infrastructure element of any open source project. 

There have been recent calls to use R more broadly in the field of hydrology [@moore2017watershed]. The `tidyhydat` package is an effort to push this call forward by being a standard package by which hydrologists and other users interact with WSC data in R. Conducting hydrological analysis in a programming environment like R allows hydrologists the ability to create fully reproducible workflows, automate repetitive tasks and provide the same rigour to the data analysis process that hydrologists apply to field equipment and experimental design [@wilson2014best].

## Why use tidy data?
Embedded within `tidyhydat` is the principle of *tidy data*. @wickham2014tidy defines tidy data by three principles:

- Each variable forms a column
- Each observation forms a row
- Each type of observational unit forms a table

It is illustrative here to provide an example of the types of data *tidying* processes that `tidyhydat` does for you automatically. The raw `DLY_FLOWS` table in the HYDAT database returns data that looks like this:

```
## # Source:   table<DLY_FLOWS> [?? x 73]
## # Database: sqlite 3.19.3
## #   [C:\Users\salbers\R\win-library\3.4\tidyhydat\test_db\tinyhydat.sqlite3]
##    STATION_NUMBER  YEAR MONTH FULL_MONTH NO_DAYS MONTHLY_MEAN
##             <chr> <int> <int>      <int>   <int>        <dbl>
##  1        05AA008  1910     7          0      31           NA
##  2        05AA008  1910     8          1      31         3.08
##  3        05AA008  1910     9          1      30         3.18
##  4        05AA008  1910    10          1      31         5.95
##  5        05AA008  1911     1          1      31         1.42
##  6        05AA008  1911     2          1      28         1.31
##  7        05AA008  1911     3          1      31         1.65
##  8        05AA008  1911     4          1      30         6.33
##  9        05AA008  1911     5          1      31        18.20
## 10        05AA008  1911     6          1      30        24.20
## # ... with more rows, and 67 more variables: MONTHLY_TOTAL <dbl>,
## #   FIRST_DAY_MIN <int>, MIN <dbl>, FIRST_DAY_MAX <int>, MAX <dbl>,
## #   FLOW1 <dbl>, FLOW_SYMBOL1 <chr>, FLOW2 <dbl>, FLOW_SYMBOL2 <chr>,
## #   FLOW3 <dbl>, FLOW_SYMBOL3 <chr>, FLOW4 <dbl>, FLOW_SYMBOL4 <chr>,
## #   FLOW5 <dbl>, FLOW_SYMBOL5 <chr>, FLOW6 <dbl>, FLOW_SYMBOL6 <chr>,
## #   FLOW7 <dbl>, FLOW_SYMBOL7 <chr>, FLOW8 <dbl>, FLOW_SYMBOL8 <chr>,
## #   FLOW9 <dbl>, FLOW_SYMBOL9 <chr>, FLOW10 <dbl>, FLOW_SYMBOL10 <chr>,
## #   FLOW11 <dbl>, FLOW_SYMBOL11 <chr>, FLOW12 <dbl>, FLOW_SYMBOL12 <chr>,
## #   FLOW13 <dbl>, FLOW_SYMBOL13 <chr>, FLOW14 <dbl>, FLOW_SYMBOL14 <chr>,
## #   FLOW15 <dbl>, FLOW_SYMBOL15 <chr>, FLOW16 <dbl>, FLOW_SYMBOL16 <chr>,
## #   FLOW17 <dbl>, FLOW_SYMBOL17 <chr>, FLOW18 <dbl>, FLOW_SYMBOL18 <chr>,
## #   FLOW19 <dbl>, FLOW_SYMBOL19 <chr>, FLOW20 <dbl>, FLOW_SYMBOL20 <chr>,
## #   FLOW21 <dbl>, FLOW_SYMBOL21 <chr>, FLOW22 <dbl>, FLOW_SYMBOL22 <chr>,
## #   FLOW23 <dbl>, FLOW_SYMBOL23 <chr>, FLOW24 <dbl>, FLOW_SYMBOL24 <chr>,
## #   FLOW25 <dbl>, FLOW_SYMBOL25 <chr>, FLOW26 <dbl>, FLOW_SYMBOL26 <chr>,
## #   FLOW27 <dbl>, FLOW_SYMBOL27 <chr>, FLOW28 <dbl>, FLOW_SYMBOL28 <chr>,
## #   FLOW29 <dbl>, FLOW_SYMBOL29 <chr>, FLOW30 <dbl>, FLOW_SYMBOL30 <chr>,
## #   FLOW31 <dbl>, FLOW_SYMBOL31 <chr>
```

This data structure clearly violates the principles of tidy data - this is messy data. For example, column headers (e.g. `FLOW1`) contain the day number - a value. HYDAT is structured like this for very reasonable historical reasons. It does, however, significantly limit a hydrologists ability to efficiently use hydrometric data. 

`tidyhydat` aims to make interacting with WSC data sources simpler. I have applied tidy data principles so that users can avoid thinking about the basic data process of importing and tidying and focus on the iterative process of visualizing and modelling their data [@wickham2016r]. After loading `tidyhydat` itself, we simply need to supply a `station_number` argument to the `hy_daily_flows()` function:


```r
library(tidyhydat)
hy_daily_flows(station_number = "08MF005")
```

```
## # A tibble: 37,561 x 5
##    STATION_NUMBER       Date Parameter Value Symbol
##             <chr>     <date>     <chr> <dbl>  <chr>
##  1        08MF005 1912-03-01      FLOW   538   <NA>
##  2        08MF005 1912-03-02      FLOW   538   <NA>
##  3        08MF005 1912-03-03      FLOW   538   <NA>
##  4        08MF005 1912-03-04      FLOW   538   <NA>
##  5        08MF005 1912-03-05      FLOW   538   <NA>
##  6        08MF005 1912-03-06      FLOW   538   <NA>
##  7        08MF005 1912-03-07      FLOW   479   <NA>
##  8        08MF005 1912-03-08      FLOW   479   <NA>
##  9        08MF005 1912-03-09      FLOW   459   <NA>
## 10        08MF005 1912-03-10      FLOW   459   <NA>
## # ... with 37,551 more rows
```

As you can see, this is much tidier data and is much easier to work with. In addition to these tidy principles, specific to `tidyhydat`, we can also define that *for a common data source, variables should be referred to by a common name*. For example, hydrometric stations are given a unique 7 digit identifier that contains important watershed information. This identifier is variously referred to as `STATION_NUMBER` or `ID` depending on the exact ECCC data source. To tidy this hydrometric data, we have renamed, where necessary, each instance of the unique identifier as `STATION_NUMBER`. This consistency to data formats, and in particular tidy data, situates `tidyhydat` well to interact seamlessly with the powerful tools being developed in the `tidyverse` [@wickham2017tidyverse] and provides a path in R to realize some of the goals outlined by @moore2017watershed.

# References
## revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.0.5 (2021-03-31) |
|os       |Windows 10 x64               |
|system   |x86_64, mingw32              |
|ui       |RStudio                      |
|language |(EN)                         |
|collate  |English_Canada.1252          |
|ctype    |English_Canada.1252          |
|tz       |America/Los_Angeles          |
|date     |2021-05-17                   |

# Dependencies

|package    |old      |new      |<U+0394>  |
|:----------|:--------|:--------|:--|
|tidyhydat  |0.5.2    |0.5.3    |*  |
|askpass    |1.1      |1.1      |   |
|assertthat |0.2.1    |0.2.1    |   |
|BH         |1.75.0-0 |1.75.0-0 |   |
|bit        |4.0.4    |4.0.4    |   |
|bit64      |4.0.5    |4.0.5    |   |
|blob       |1.2.1    |1.2.1    |   |
|cachem     |1.0.4    |1.0.4    |   |
|cli        |2.5.0    |2.5.0    |   |
|clipr      |0.7.1    |0.7.1    |   |
|cpp11      |0.2.7    |0.2.7    |   |
|crayon     |1.4.1    |1.4.1    |   |
|curl       |4.3.1    |4.3.1    |   |
|DBI        |1.1.1    |1.1.1    |   |
|dbplyr     |2.1.1    |2.1.1    |   |
|dplyr      |1.0.6    |1.0.6    |   |
|ellipsis   |0.3.2    |0.3.2    |   |
|fansi      |0.4.2    |0.4.2    |   |
|fastmap    |1.1.0    |1.1.0    |   |
|generics   |0.1.0    |0.1.0    |   |
|glue       |1.4.2    |1.4.2    |   |
|hms        |1.1.0    |1.1.0    |   |
|httr       |1.4.2    |1.4.2    |   |
|jsonlite   |1.7.2    |1.7.2    |   |
|lifecycle  |1.0.0    |1.0.0    |   |
|lubridate  |1.7.10   |1.7.10   |   |
|magrittr   |2.0.1    |2.0.1    |   |
|memoise    |2.0.0    |2.0.0    |   |
|mime       |0.10     |0.10     |   |
|openssl    |1.4.4    |1.4.4    |   |
|pillar     |1.6.1    |1.6.1    |   |
|pkgconfig  |2.0.3    |2.0.3    |   |
|plogr      |0.2.0    |0.2.0    |   |
|purrr      |0.3.4    |0.3.4    |   |
|R6         |2.5.0    |2.5.0    |   |
|rappdirs   |0.3.3    |0.3.3    |   |
|Rcpp       |1.0.6    |1.0.6    |   |
|readr      |1.4.0    |1.4.0    |   |
|rlang      |0.4.11   |0.4.11   |   |
|RSQLite    |2.2.7    |2.2.7    |   |
|sys        |3.4      |3.4      |   |
|tibble     |3.1.1    |3.1.1    |   |
|tidyr      |1.1.3    |1.1.3    |   |
|tidyselect |1.1.1    |1.1.1    |   |
|utf8       |1.2.1    |1.2.1    |   |
|vctrs      |0.3.8    |0.3.8    |   |
|withr      |2.4.2    |2.4.2    |   |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*