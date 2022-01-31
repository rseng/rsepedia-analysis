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

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
title: README
output: md_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->
    
```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-"
)
```

# tidyhydat <img src="man/figures/logo.png" align="right" />


<!-- badges: start -->
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0) [![Coverage status](https://codecov.io/gh/ropensci/tidyhydat/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/tidyhydat?branch=master) 
[![R build status](https://github.com/ropensci/tidyhydat/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/tidyhydat/actions)

[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/tidyhydat)](https://cran.r-project.org/package=tidyhydat) [![CRAN Downloads](https://cranlogs.r-pkg.org/badges/tidyhydat?color=brightgreen)](https://CRAN.R-project.org/package=tidyhydat) [![cran checks](https://cranchecks.info/badges/worst/tidyhydat)](https://cran.r-project.org/web/checks/check_results_tidyhydat.html)  [![r-universe](https://ropensci.r-universe.dev/badges/tidyhydat)](https://ropensci.r-universe.dev/ui#builds)


[![](http://badges.ropensci.org/152_status.svg)](https://github.com/ropensci/software-review/issues/152)  [![DOI](http://joss.theoj.org/papers/10.21105/joss.00511/status.svg)](https://doi.org/10.21105/joss.00511) [![DOI](https://zenodo.org/badge/100978874.svg)](https://zenodo.org/badge/latestdoi/100978874) 
<!-- badges: end -->


## Project Status

This package is maintained by the Data Science Partnerships Program in the [British Columbia Ministry of Citizens' Services](https://www2.gov.bc.ca/gov/content/governments/organizational-structure/ministries-organizations/ministries/citizens-services).


## What does `tidyhydat` do?

- Provides functions (`hy_*`) that access hydrometric data from the HYDAT database, a national archive of Canadian hydrometric data and return tidy data.
- Provides functions (`realtime_*`) that access Environment and Climate Change Canada's real-time hydrometric data source.
- Provides functions (`search_*`) that can search through the approximately 7000 stations in the database and aid in generating station vectors
- Keep functions as simple as possible. For example, for daily flows, the `hy_daily_flows()` function queries the database, *tidies* the data and returns a [tibble](https://tibble.tidyverse.org/) of daily flows.

## Installation
You can install `tidyhydat` from CRAN:
```{r, echo=TRUE, eval=FALSE}
install.packages("tidyhydat")
```


To install the development version of the `tidyhydat` package, you can install directly from the rOpenSci development server:
```{r, echo=TRUE, eval=FALSE}
install.packages("tidyhydat", repos = "https://dev.ropensci.org")
```

## Usage
More documentation on `tidyhydat` can found at the rOpenSci doc page: https://docs.ropensci.org/tidyhydat/

When you install `tidyhydat`, several other packages will be installed as well. One of those packages, `dplyr`, is useful for data manipulations and is used regularly here. To use actually use `dplyr` in a session you must explicitly load it. A helpful `dplyr` tutorial can be found [here](https://cran.r-project.org/package=dplyr/vignettes/dplyr.html).
  
```{r, eval = TRUE, echo=TRUE, message=FALSE, warning=FALSE}
library(tidyhydat)
library(dplyr)
```
  
### HYDAT download
To use many of the functions in the `tidyhydat` package you will need to download a version of the HYDAT database, Environment and Climate Change Canada's database of historical hydrometric data then tell R where to find the database. Conveniently `tidyhydat` does all this for you via:
```{r, eval=FALSE}
download_hydat()
```
This downloads (with your permission) the most recent version of HYDAT and then saves it in a location on your computer where `tidyhydat`'s function will look for it. Do be patient though as this can take a long time! To see where HYDAT was saved you can run `hy_default_db()`. Now that you have HYDAT downloaded and ready to go, you are all set to begin looking at Canadian hydrometric data.

### Real-time
To download real-time data using the datamart we can use approximately the same conventions discussed above. Using `realtime_dd()` we can easily select specific stations by supplying a station of interest:
```{r}
realtime_dd(station_number = "08LG006")
```

### Plotting

Plot methods are also provided to quickly visualize realtime data:
```{r}
realtime_ex <- realtime_dd(station_number = "08LG006")

plot(realtime_ex)
```

and also historical data:
```{r, fig.height=7, fig.width=12}
hy_ex <- hy_daily_flows(station_number = "08LA001", start_date = "2013-01-01")

plot(hy_ex)
```

## Getting Help or Reporting an Issue

To report bugs/issues/feature requests, please file an [issue](https://github.com/ropensci/tidyhydat/issues/).

These are very welcome!

## How to Contribute

If you would like to contribute to the package, please see our 
[CONTRIBUTING](https://github.com/ropensci/tidyhydat/blob/master/CONTRIBUTING.md) guidelines.
  
Please note that this project is released with a [Contributor Code of Conduct](https://github.com/ropensci/tidyhydat/blob/master/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.
  
## Citation
Get citation information for `tidyhydat` in R by running:
```{r, echo=FALSE, comment=""}
citation("tidyhydat")
```



[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

## License
  
  Copyright 2017 Province of British Columbia
  
  Licensed under the Apache License, Version 2.0 (the "License");
  you may not use this file except in compliance with the License.
  You may obtain a copy of the License at 
  
  https://www.apache.org/licenses/LICENSE-2.0
  
  Unless required by applicable law or agreed to in writing, software
  distributed under the License is distributed on an "AS IS" BASIS,
  WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
  See the License for the specific language governing permissions and
  limitations under the License.
  
---
title: "tidyhydat: Extract and Tidy Canadian Hydrometric Data"
authors:
- affiliation: 1
  name: Sam Albers
  orcid: 0000-0002-9270-7884
date: "`r Sys.Date()`"
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
  index: 1
  name: Hydrology and Hydrometric Programs, Ministry of Environment and Climate Change Strategy, British Columbia Provincial Government
---

> "Tidy datasets are all alike but every messy dataset is messy in its own way - "
@wickham2014tidy

```{r options, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, warning = FALSE, messages = FALSE, fig.width = 8, fig.height = 12)
```

```{r packages, warning=FALSE, message=FALSE, echo = FALSE}
library(tidyhydat)
library(dplyr)
library(dbplyr)
```

# Introduction
Environment and Climate Change Canada (ECCC) through the Water Survey of Canada (WSC) maintains several national hydrometric data sources. These data are partially funded by provincial partners and constitute the main data products of a national integrated hydrometric network. Historical data are stored in the [HYDAT database](http://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/). HYDAT is the Canadian national Water Data Archive, published quarterly by the Government of Canada's Department of Environment and Climate Change. It is relational database that contains daily, monthly and annual data on water flow, water levels and sediment.

Real-time data are provided by ECCC over the web. Files are updated to a [datamart](http://dd.weather.gc.ca/hydrometric/) on an hourly basis though the lag between actual hydrometric measurement and the availability of hydrometric data is approximately 2.5 hours. The objective of this document is the outline the usage of `tidyhydat` [@alberstidyhydat], an R package that accesses these hydrometric data sources and *tidies* them. `tidyhydat` is part of the [rOpenSci](https://ropensci.org/packages/) suite of packages and resides at  https://github.com/ropensci/tidyhydat. The objective of `tidyhydat` is to provide a standard method of accessing ECCC data sources using a consistent and easy to use interface that employs tidy data principles developed by @wickham2014tidy within the R project [@RCore]. 

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
```{r, echo = FALSE}
hydat_con <- DBI::dbConnect(RSQLite::SQLite(), system.file("test_db/tinyhydat.sqlite3", package = "tidyhydat"))
tbl(hydat_con, "DLY_FLOWS") 
```

This data structure clearly violates the principles of tidy data - this is messy data. For example, column headers (e.g. `FLOW1`) contain the day number - a value. HYDAT is structured like this for very reasonable historical reasons. It does, however, significantly limit a hydrologists ability to efficiently use hydrometric data. 

`tidyhydat` aims to make interacting with WSC data sources simpler. I have applied tidy data principles so that users can avoid thinking about the basic data process of importing and tidying and focus on the iterative process of visualizing and modelling their data [@wickham2016r]. After loading `tidyhydat` itself, we simply need to supply a `station_number` argument to the `hy_daily_flows()` function:

```{r, echo = TRUE, message=FALSE}
library(tidyhydat)
hy_daily_flows(station_number = "08MF005")
```

As you can see, this is much tidier data and is much easier to work with. In addition to these tidy principles, specific to `tidyhydat`, we can also define that *for a common data source, variables should be referred to by a common name*. For example, hydrometric stations are given a unique 7 digit identifier that contains important watershed information. This identifier is variously referred to as `STATION_NUMBER` or `ID` depending on the exact ECCC data source. To tidy this hydrometric data, we have renamed, where necessary, each instance of the unique identifier as `STATION_NUMBER`. This consistency to data formats, and in particular tidy data, situates `tidyhydat` well to interact seamlessly with the powerful tools being developed in the `tidyverse` [@wickham2017tidyverse] and provides a path in R to realize some of the goals outlined by @moore2017watershed.

# References
---
title: "tidyhydat"
subtitle: "Making the case for reproducible workflows in hydrology"
author: "Sam Albers"
date: 2018-03-06
output:
  xaringan::moon_reader:
    keep_md: true
    lib_dir: libs
    css: ["default", "hygge_sam.css"]
    nature:
      highlightStyle: github
      highlightLines: true
      countIncrementalSlides: false
---

```{r setup, include=FALSE}
options(htmltools.dir.version = FALSE)


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

knitr::opts_chunk$set(
  collapse = TRUE,
  echo = FALSE,
  comment = "#>",
  fig.path = "graphics/prod/figs"
)

options(scipen = 10)
```

```{r, pck-load, warning=FALSE, message=FALSE}
library(tidyhydat)
library(knitr)
library(tidyverse)
library(lubridate)
library(corrr)
library(leaflet)
library(sf)
library(mapview)

```

class: inverse
background-image: url(https://upload.wikimedia.org/wikipedia/commons/3/3e/Clearwater_River_Wells_Gray_Park.jpg)
background-size: cover

## Outline

.VeryLarge[
- Common Analysis Problems
- What is R and why use it?
- What is tidyhydat?
- Some R basics
- An example of how R can help
- Leveraging R and what I'm not showing you
- Where to get help
- Questions
]

---
class: inverse, center, middle
# Common Analysis Problems
---
class: center, basic
## Accessing Hydrometric Data
```{r data-explorer, out.width = "85%"}
include_graphics("graphics/ec_data_explorer2.gif")
```

### 11 clicks!

---

class: basic, center

### Stakeholder/Manager: "Hey, this is a really cool analysis but we need to add five stations. Can you run it again?"
<img src="https://media.giphy.com/media/l4FGuE8LtZg2wKmZi/giphy.gif"/>

--
### Make it reproducible!

---
class: basic, center

### Get off the factory line
### How much time do you spend copying and pasting?

<img src="https://media.giphy.com/media/3ohhwmbdlEf9CtCgBG/giphy.gif"/>


--
### Automate!
--

### But how...
---
class: inverse, left, middle
## ...Use R!
.pull-left[
(or more generally any programmatic code based analysis approach...)
]
<center><img src="https://www.r-project.org/logo/Rlogo.png" alt="Drawing" style="width: 450px;" /></center>
---
.pull-left[
### What is R?
.large[
- Free and open source
- Statistical programming language
- Publication quality graphics
- But definitely not intimidating...
]
]

--
.pull-right[
### Why use R?

.large[
- Efficient
- Reproducible
- Scalable
]
]

--
### Not guaranteed to help with this...
<center><img src="https://media.giphy.com/media/HteV6g0QTNxp6/giphy.gif" style="width: 450px;"/> </center>

---

## Questions worth asking...
.large[
- Are your methods <span style="color:#309688">reproducible</span>?
- What is your analysis recipe? 
- Can you share it?
]

<center><img src="graphics/data_recipe.png" alt="Drawing" style="width: 750px;" /></center>
---

## Excuse me,  do you have a  moment to talk about Excel?

<img src="http://ichef.bbci.co.uk/news/976/cpsprodpb/CED0/production/_89744925_clippy.jpg"/>

---
class:basic

| R                                         | Excel                                                  |
|-------------------------------------------|--------------------------------------------------------|
| Data and analysis are separate            | Data and analysis are usually stored in the same place |

<br>

<center><img src="https://i0.wp.com/civilengineerspk.com/wp-content/uploads/2016/04/excel2010_06.jpg" alt="Drawing" style="width: 450px;" /></center>

.footnote[
From: http://blog.yhat.com/posts/R-for-excel-users.html.
]

---
class:basic

| R                                         | Excel                                                  |
|-------------------------------------------|--------------------------------------------------------|
| Data and analysis are separate            | Data and analysis are usually stored in the same place |
| Data structure is strict                  | Data structure is flexible                             |

<br>

<center><img src="https://i.stack.imgur.com/d9cMS.jpg" alt="Drawing" style="width: 450px;" /></center>


.footnote[
From: http://blog.yhat.com/posts/R-for-excel-users.html.
]

---
class:basic

| R                                         | Excel                                                  |
|-------------------------------------------|--------------------------------------------------------|
| Data and analysis are separate            | Data and analysis are usually stored in the same place |
| Data structure is strict                  | Data structure is flexible                             |
| Operations are achieved through scripting | Operations are achieved through pointing and clicking  |

<br>

<center><img src="https://media.giphy.com/media/UMbqRmfi03SQE/giphy.gif" alt="Drawing" style="width: 400px;" /></center>



.footnote[
From: http://blog.yhat.com/posts/R-for-excel-users.html.
]

---
class:basic

| R                                         | Excel                                                  |
|-------------------------------------------|--------------------------------------------------------|
| Data and analysis are separate            | Data and analysis are usually stored in the same place |
| Data structure is strict                  | Data structure is flexible                             |
| Operations are achieved through scripting | Operations are achieved through pointing and clicking  |
| Iteration is automated                    | Iteration is usually done by hand                      |

### R provides a clear pathway for <span style="color:#309688">efficiency</span> and <span style="color:#309688">reproducibility</span> through automation and code

.footnote[
From: http://blog.yhat.com/posts/R-for-excel-users.html.
]

---
class:basic

> The objective of tidyhydat is to provide a standard method of accessing ECCC hydrometric data sources (historical and real time) using a consistent and easy to use interface that employs tidy data principles within the R project.

<center><img src="https://raw.githubusercontent.com/ropensci/tidyhydat/master/tools/readme/tidyhydat_large.png" alt="Drawing" style="width: 300px;" /></center>

--

# <center><span style="color:#309688">tidy</span><span style="color:red">|</span><span style="color:#234075">hydat</span></center>


---

## hydat::Water Survey of Canada Network

```{r out.width='100%', fig.height=6, eval=TRUE, message=FALSE}
stns <- hy_stations() %>% 
  filter(HYD_STATUS == "ACTIVE")

st_as_sf(stns, coords = c("LONGITUDE","LATITUDE"),
             crs = 4326,
             agr= "constant") %>%
mapview(zcol = "STATION_NAME", legend = FALSE, map.types = "Esri.WorldImagery", cex = 4,
        popup = popupTable(., zcol = c("STATION_NUMBER", "STATION_NAME", "PROV_TERR_STATE_LOC")))

#leaflet(data = stns) %>% 
#  addTiles() %>% 
#  addMarkers(~LONGITUDE, ~LATITUDE, label=~as.character(STATION_NAME), clusterOptions = markerClusterOptions()) %>% #
#  setView(-96, 63, zoom = 3)
```


---

## <span style="color:#309688">tidy</span>::tidy data
> Tidy datasets are all alike but every messy dataset is messy in its own way<sup>1</sup>

--

### Each variable forms a column

### Each observation forms a row

<center><img src="https://media.giphy.com/media/kXBVtKjLxINji/giphy.gif" style="width: 450px;"/> </center>

.footnote[
[1] [Wickham, Hadley. 2014. Tidy Data. Journal of Statistical Software 59 (10). Foundation for Open Access Statistics: 1–23.](https://www.jstatsoft.org/article/view/v059i10)
]
---

## <span style="color:#309688">tidy</span>::<span style="color:red">un</span>tidy data
```{r, echo = FALSE}
src <- hy_src()
tbl(src, "DLY_FLOWS") %>%
  filter(STATION_NUMBER == "08MF005") %>%
  select(-contains("_SYMBOL"), )
```

---
## <span style="color:#309688">tidy</span>::tidy data
```{r, echo = FALSE, message=FALSE}
hy_daily_flows(station_number = "08MF005")
```

<center><img src="https://media.giphy.com/media/kXBVtKjLxINji/giphy.gif" style="width: 450px;"/> </center>
---

## <span style="color:#309688">tidy</span>::tidyhydat
<center><img src="graphics/data-science-wheel.png" alt="Drawing" style="width: 800px;" /></center>

--

<center><img src="https://cdn.onlinewebfonts.com/svg/img_411268.png" style="width: 200px;"/> </center>

---
class: inverse, center, middle
# An Example 
---

class: basiclh

## tidyhydat & some basic R
```{r excel_ex, echo = TRUE, eval = FALSE, message=FALSE}
=SUM(A1:A23)
=AVERAGE(A1:A23)
```


---

class: basiclh

## tidyhydat & some basic R
```{r flow_ex, echo = TRUE, eval = TRUE, message=FALSE}
flows_data <- hy_daily_flows(station_number = c("08MF005","09CD001","05KJ001","02KF005"))
```

- `<-`: <span style="color:#309688">assignment operator</span>
- `flows_data`: <span style="color:#309688">object</span>
- `hy_daily_flows`: <span style="color:#309688">function</span>
- `station_number`: <span style="color:#309688">argument</span>
--

```{r, echo = TRUE, eval = TRUE, message=FALSE}
flows_data
```


---
class: basiclh, center

### Analyze the correlation between:
```{r, eval=TRUE, echo=FALSE, message=FALSE}
stns_tbl <- hy_stations(c("08MF005","09CD001","05KJ001","02KF005"))[,c("STATION_NUMBER", "STATION_NAME")]

x <- stns_tbl %>%
  rename(`Station Name`=STATION_NAME, `Station Number`=STATION_NUMBER) %>%
  knitr::kable(format = 'html') %>%
  kableExtra::kable_styling(bootstrap_options = c("striped", "hover"))
gsub("<thead>.*</thead>", "", x)
```

```{r, message=FALSE, fig.width = 9, fig.height=5}
stns <- hy_stations(stns_tbl$STATION_NUMBER) 


st_as_sf(stns, coords = c("LONGITUDE","LATITUDE"),
             crs = 4326,
             agr= "constant") %>%
mapview(zcol = "STATION_NAME", legend = FALSE, cex = 6,
        popup = popupTable(., zcol = c("STATION_NUMBER", "STATION_NAME", "PROV_TERR_STATE_LOC")))
```


---

## Build the analysis
```{r, echo = TRUE}
flows_data
```
### `flows_data`: <span style="color:#309688">object</span>

---

## Build the analysis
```{r, echo = TRUE}
flows_data %>%
  spread(key = STATION_NUMBER, value = Value) #<<
```

### `%>%`: <span style="color:#309688">"then"</span>
### `spread`: <span style="color:#309688">function</span>
---

## Build the analysis
```{r, echo = TRUE}
flows_data %>%
  spread(key = STATION_NUMBER, value = Value) %>%
  select(-Date, -Symbol, -Parameter) #<<
```

### `select`: <span style="color:#309688">function</span>
---

## Build the analysis
```{r, echo = TRUE}
flows_data %>%
  spread(key = STATION_NUMBER, value = Value) %>%
  select(-Date, -Symbol, -Parameter) %>%
  correlate() #<<
```

### `correlation`: <span style="color:#309688">function</span>
---

## Build the analysis
```{r, echo = TRUE}
flows_data %>%
  spread(key = STATION_NUMBER, value = Value) %>%
  select(-Date, -Symbol, -Parameter) %>%
  correlate() %>%
  stretch() #<<
```

### `stretch`: <span style="color:#309688">function</span>

---

## Scalable
```{r scalable1, echo = FALSE, message=FALSE, fig.width = 10, fig.height=6}
stns <- hy_stations(prov_terr_state_loc = "NU") %>%
  filter(HYD_STATUS == "ACTIVE")

st_as_sf(stns, coords = c("LONGITUDE","LATITUDE"),
             crs = 4326,
             agr= "constant") %>%
mapview(zcol = "STATION_NAME", legend = FALSE, 
        popup = popupTable(., zcol = c("STATION_NUMBER", "STATION_NAME", "PROV_TERR_STATE_LOC", "HYD_STATUS")))
```

---

## Scalable
```{r scalable2, echo = TRUE, message=FALSE}
stns <- hy_stations(prov_terr_state_loc = "NU") %>%
  filter(HYD_STATUS == "ACTIVE")

nu_flows <- hy_daily_flows(station_number = stns$STATION_NUMBER)
nu_flows
```

---

## Scalable
```{r scalable3, echo = TRUE, message=FALSE}
nu_flows %>% #<<
  spread(STATION_NUMBER, Value) %>%
  select(-Date, -Symbol, -Parameter) %>%
  correlate() %>% 
  stretch() 
```

---

## Efficient, Reproducible and Scalable
<center><img src="graphics/data-science-wheel.png" alt="Drawing" style="width: 800px;" /></center>
---

## What else is available in `tidyhydat`?

### All tables in HYDAT 
.Large[
- Instantaneous peaks
- Daily, monthly and yearly temporal summaries
- Discharge, level, sediment, particle size
- Data ranges
- Station metadata
]

---

## What else is available in `tidyhydat`?
```{r, echo = TRUE, eval = TRUE}
search_stn_name("fraser")
```

---

## Pointing and clicking
```{r water_office, out.width = "95%"}
include_graphics("graphics/wateroffice.gif")
```
---
## What else is available in `tidyhydat`?
```{r, message=FALSE, warning=FALSE, fig.width=11, fig.height=5, echo = TRUE, eval = TRUE}
realtime_plot("08MF005", Parameter = "Flow")
```


---
## What else is available in R?
```{r, warning=FALSE, message=FALSE, echo=TRUE}
raw_stns <- hy_stations() %>%
  select(STATION_NUMBER:PROV_TERR_STATE_LOC, DRAINAGE_AREA_GROSS)
  
mad_long_avg <- hy_annual_stats(raw_stns$STATION_NUMBER) %>%
  filter(Sum_stat == "MEAN", Parameter == "Flow") %>%
  group_by(STATION_NUMBER) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  right_join(raw_stns)
mad_long_avg #<<
```

---
## What else is available in R?

```{r, warning=FALSE, message=FALSE, echo=TRUE, fig.height=5, fig.width=11}
library(ggplot2)
ggplot(mad_long_avg,aes(x = Value, y = DRAINAGE_AREA_GROSS, colour = PROV_TERR_STATE_LOC)) +
  geom_point() +
  scale_y_continuous(trans = "log10") +
  scale_x_continuous(trans = "log10") +
  scale_colour_viridis_d(name = "Jurisdiction") +
  labs(x = "Mean long term annual discharge (m^3)", y = "Gross drainage area (km^2)") +
  theme_minimal()
```

---

## It can be daunting!

<center><img src="https://media.giphy.com/media/iKhphMBOk5ScE/giphy.gif" style="width: 650px;"/></center>

---

## Resources for R

<img src="https://cdn.sstatic.net/Sites/stackoverflow/company/img/logos/so/so-logo.svg?v=2bb144720a66" alt="Drawing" style="width: 400px;" />


<img src="https://www.rstudio.com/wp-content/uploads/2017/11/logoRStudioCommunity.svg" alt="Drawing" style="width: 400px;" />

<center><img src="https://www.r-project.org/logo/Rlogo.png" alt="Drawing" style="width: 300px;" /></center>





---
## Contribute to `tidyhydat`

### Openly developed on GitHub <img src = "https://assets-cdn.github.com/images/modules/logos_page/GitHub-Mark.png" style = "width: 35px; vertical-align: middle; margin-bottom: 3px;">
<https://github.com/ropensci/tidyhydat>

Any contribution helps. You don't have to be an R programmer!

.pull-left[
- Questions
- Ideas / Feature-requests
- Bugs
- Bug-fixes
- Development
]
.pull-right[
<center><img src="https://raw.githubusercontent.com/ropensci/tidyhydat/master/tools/readme/tidyhydat_large.png" alt="Drawing" style="width: 150px;" /></center>
]
--

### For example...
```{r, eval=FALSE, echo = TRUE}
Authors@R: c(person("Sam", "Albers",email = "sam.albers@gov.bc.ca", role = c("aut", "cre")),
    person("David", "Hutchinson", email = "david.hutchinson@canada.ca", role = "ctb"), #<<
    person("Dewey", "Dunnington", email = "dewey@fishandwhistle.net", role = "ctb"), #<<
    person("Province of British Columbia", role = "cph"))
```
---
class: inverse, center
## Some Helpful Links

Installing R & RStudio with local package libraries

    -https://github.com/bcgov/bcgov-data-science-resources/wiki/Installing-R-&-RStudio

Installing tidyhydat

    -https://cran.rstudio.com/web/packages/tidyhydat/README.html

Getting started with `tidyhydat`

    -https://cran.rstudio.com/web/packages/tidyhydat/vignettes/tidyhydat_an_introduction.html
    -https://cran.rstudio.com/web/packages/tidyhydat/vignettes/tidyhydat_example_analysis.html

BC Gov data science resource wiki

    -https://github.com/bcgov/bcgov-data-science-resources/wiki
    
---

class: basic
background-image: url(https://media.giphy.com/media/TnDoEoXfT7YoE/giphy.gif)
background-size: cover

## Questions?
.content-box-blue[
Slides available from 

    -https://github.com/ropensci/tidyhydat/blob/master/presentations/tidyhydat_intro.pdf
    -https://github.com/ropensci/tidyhydat/blob/master/presentations/tidyhydat_intro.Rmd
 
Contact <sam.albers@gov.bc.ca>
]


---
title: "Introduction to working with Canadian Water Data in R"
subtitle: "Using tidyhydat and weathercan"
author: "Sam Albers <br> Digital Platforms and Data Division <br> Office of the Chief Information Officer <br> Ministry of Citizens' Services <br> Province of BC <br><br> CWRA Webinar <br>"
date: 2019-09-25
output:
  xaringan::moon_reader:
    keep_md: true
    lib_dir: libs
    css: ["default", "default-fonts", "hygge"]
    nature:
      highlightStyle: github
      highlightLines: true
      countIncrementalSlides: false
      beforeInit: "https://platform.twitter.com/widgets.js"
      ratio: '16:9'
---

layout: true

---

```{r, include=FALSE}
# Copyright 2019 Province of British Columbia
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


```{r setup, include=FALSE}
options(htmltools.dir.version = FALSE)
options(width = 90)
options(max_print = 5)

knitr::opts_chunk$set(
  collapse = TRUE,
  echo = FALSE,
  comment = "#>",
  fig.path = "graphics/prod/figs"
)

options(scipen = 10)
```

```{r, pck-load, warning=FALSE, message=FALSE}
library(tidyhydat)
library(weathercan)
library(knitr)
library(tidyverse)
library(lubridate)
library(corrr)
library(sf)
library(rnaturalearth)
library(fontawesome)
```


```{r, theme, warning=FALSE, echo=FALSE}
bg_black <- "#272822"

theme_set(theme_void() %+replace%
            theme(legend.text = element_text(colour = "white", size = 18),
                  legend.title = element_text(colour = "white", size = 18),
                  plot.background = element_rect(fill = bg_black, color = bg_black),
                  axis.text = element_text(colour = "white", size = 16),
                  axis.title = element_text(colour = "white", size = 18),
                  axis.title.y = element_text(angle = 90, vjust = 1),
                  plot.title = element_text(colour = "white", size = 22, hjust = 0)))


scale_colour_continuous <- scale_colour_viridis_c
scale_fill_continuous <- scale_fill_viridis_c
scale_colour_discrete <- scale_colour_viridis_d
scale_fill_discrete <- scale_fill_viridis_d
```



## Outline

.VeryLarge[
- Who am I?
- Learning Outcomes
- Review R and RStudio and rationale behind using them
- Introduce packages:
  - `dplyr`
  - `tidyhydat`
  - `weathercan`
- Provide an example of using them together
- `tidyhydat` and `weathercan` development
- Where and how to get help in R
- Questions
]

---
## Sam Albers



.pull-left[
- Data Scientist with BC government
- Environmental Scientist by training
- Been using R for 10 years
- Maintainer for `tidyhydat`, `rsoi`
- Contributor on many other packages including `weathercan`
- Maintainer of the Hydrology task view
]
.pull-right[
<center><img src="graphics/intro_me.jpg" alt="Drawing" style="width: 350px;" /></center>

[`r fa(name = "twitter")` @big_bad_sam](https://twitter.com/big_bad_sam)  
[`r fa(name = "github")` @boshek](http://github.com/boshek)  
[`r fa(name = "paper-plane")` sam.albers@gov.bc.ca](sam.albers@gov.bc.ca)
]

---
## What are we hoping to learn?
.pull-left[
.VeryLarge[
- Describe visual elements of RStudio
- Define and assign data to variable
- Manage your workspaces and projects
- Call a function
- Understand the six main `dplyr` verbs
- Overview of tidyhydat and weathercan functions
- Describe usage of `tidyhydat` and `weathercan`
- How to ask for help in R
]
]

.pull-right[
<img src="https://i.imgflip.com/3b0y51.jpg"/>

]

---
class: inverse, center, middle
# Common Analysis Problems
---
class: center, basic
## Accessing Environment and Climate Change Canada Data
```{r data-explorer, out.width = "85%"}
include_graphics("graphics/ec_data_explorer2.gif")
```

### 11 clicks!

---

class: basic, center

### Stakeholder/Manager: "Hey, this is a really cool analysis but we need to add five stations. Can you run it again?"
<img src="https://media.giphy.com/media/l4FGuE8LtZg2wKmZi/giphy.gif"/>

--
### Make it reproducible!

---

## Questions worth asking...
.large[
- Are your methods <span style="color:#309688">reproducible</span>?
- What is your analysis recipe? 
- Can you share it?
]

<center><img src="graphics/data_recipe.png" alt="Drawing" style="width: 750px;" /></center>


---
class: inverse, left, middle
## ...Use R!
.pull-left[
(or more generally any programmatic code based analysis approach...)
]
<center><img src="https://www.r-project.org/logo/Rlogo.png" alt="Drawing" style="width: 450px;" /></center>
---
.pull-left[
### What is R?
.large[
- Free and open source
- Statistical programming language
- Publication quality graphics
- Much of the innovation occurs in contributed packages
- But definitely not intimidating...
]

### Some example code
```{r echo = TRUE, message=FALSE}
all_time_greats <- c(99, 66, 4, 9)
```

- `<-`: <span style="color:#309688">assignment operator</span>
- `all_time_greats`: <span style="color:#309688">object</span>
- `c`: <span style="color:#309688">function</span>

]

--
.pull-right[

### What is RStudio?

.large[
- Provides a place to write and run code
- A means to organize projects
- Referred to as an IDE
]
--
### Not guaranteed to help with this...
<center><img src="https://media.giphy.com/media/HteV6g0QTNxp6/giphy.gif" style="width: 450px;"/> </center>

]

---
class: inverse, left, middle
# R and RStudio

---

.pull-left[
## The Problem
- Many tasks when analyzing environmental data are repetitive yet interactive
- Typically hydrologists/water professionals aren't computer scientists
- Helpful to abstract away unneeded complexity when possible
- A clean and easy to remember syntax reduces your cognitive load when doing analysis



<center><img src="https://www.herocollector.com/Content/ArticleImages/7a716739-72cb-40d5-acfc-dfc35783d8a5.jpg" style="width: 450px;"/></center>



]

--

.pull-right[
## Enter `dplyr`
> a consistent set of verbs that help you solve the most common data manipulation challenges

- Independent of the data source
- Designed for data science

<center><img src="https://raw.githubusercontent.com/rstudio/hex-stickers/master/PNG/dplyr.png" style="width: 300px;"/></center>

]


---

##`dplyr` verbs

Functions with English meanings that map directly to the action being taken when that function is called

Installation: `install.packages("dplyr")`


.pull-left[
- `%>%` a special symbol to chain operations. Read it as "then"
- `select()` picks variables based on their names.
- `filter()` picks cases based on their values.
- `summarise()` reduces multiple values down to a single summary.
- `arrange()` changes the ordering of the rows.
- `mutate()` adds new variables that are functions of existing variables

For a offline tutorial: http://swcarpentry.github.io/r-novice-gapminder/13-dplyr/index.html
]


.pull-right[
<center><img src="https://raw.githubusercontent.com/allisonhorst/stats-illustrations/master/rstats-artwork/dplyr_wrangling.png" style="width: 450px;"/></center>

Artwork by [@allison_horst](https://twitter.com/allison_horst)
] 

---
class: inverse, left, middle
# dplyr code break

---
class:basic

> The objective of tidyhydat is to provide a standard method of accessing ECCC hydrometric data sources (historical and real time) using a consistent and easy to use interface that employs tidy data principles within the R project.

<center><img src="https://github.com/ropensci/tidyhydat/raw/master/man/figures/tidyhydat_large.png" alt="Drawing" style="width: 300px;" /></center>


Installation: `install.packages("tidyhydat")`

---

## hydat::Water Survey of Canada Network

.pull-left[
```{r eval=TRUE, message=FALSE, cache=TRUE, fig.width=6.9, fig.height=5.9, fig.align='center'}
stns <- hy_stations() %>% 
  filter(HYD_STATUS == "ACTIVE")

stns_sf <- st_as_sf(stns, coords = c("LONGITUDE","LATITUDE"),
             crs = 4326,
             agr= "constant") 

can <- ne_countries(country = "Canada", returnclass = "sf")

ggplot() +
  geom_sf(data = can, fill = NA) +
  geom_sf(data = stns_sf, size = 1, aes(colour = PROV_TERR_STATE_LOC)) +
  guides(colour = FALSE) +
  coord_sf(crs = 102009, datum = NA) +
  theme_void()
```
]

.pull-right[

## `r round(fs::file_size(file.path(hy_dir(), "Hydat.sqlite3"))/1E9,2)` GB
## `r nrow(hy_stations())` stations in database
## SQLite database
## Self contained
]

---
class:basic

> The objective of weathercan is to provide a standard method of accessing ECCC climate data sources using a consistent and easy to use interface that employs tidy data principles within the R project.

<center><img src="https://raw.githubusercontent.com/ropensci/weathercan/master/inst/assets/weathercan_logo.png" alt="Drawing" style="width: 300px;" /></center>


Installation: `install.packages("weathercan")`
---

## weathercan::Climate Data

.pull-left[
```{r eval=TRUE, message=FALSE, cache=TRUE, fig.width=6.9, fig.height=5.9, fig.align='center'}
stations_unique <- stations %>% 
  select(prov, station_name, lat, lon) %>% 
  unique()



weather_stns_sf <- stations_unique %>%
  filter(!station_name == "POINT LEPREAU") %>% 
  filter(!is.na(lat), !is.na(lon)) %>% 
  st_as_sf(coords = c("lon","lat"),
             crs = 4326,
             agr= "constant") 

can <- ne_countries(country = "Canada", returnclass = "sf")

ggplot() +
  geom_sf(data = can, fill = NA) +
  geom_sf(data = weather_stns_sf, size = 1, aes(colour = prov)) +
  guides(colour = FALSE) +
  coord_sf(crs = 102009, datum = NA) +
  theme_void()
```
]

.pull-right[


## `r length(unique(stations$station_name))` stations 
## Available online
]
---

class: inverse, center, middle
# Looking closer at `tidyhydat` and `weathercan`

---

class: basiclh

## tidyhydat

Download the database:
```{r download, echo = TRUE, eval = FALSE, message=FALSE}
download_hydat()
```

Access some flow data
```{r flow_ex, echo = TRUE, eval = TRUE, message=FALSE}
flows_data <- hy_daily_flows(station_number = c("08MF005","09CD001","05KJ001","02KF005"))
```

- `<-`: <span style="color:#309688">assignment operator</span>
- `flows_data`: <span style="color:#309688">object</span>
- `hy_daily_flows`: <span style="color:#309688">function</span>
- `station_number`: <span style="color:#309688">argument</span>

---

## What else is available in `tidyhydat`?

### All tables in HYDAT 
.Large[
- See `help(package = "tidyhydat")`
- Realtime data
- Instantaneous peaks
- Daily, monthly and yearly temporal summaries
- Discharge, level, sediment, particle size
- Data ranges
- Station metadata
]

---
## What else is available in `tidyhydat`?
```{r, message=FALSE, warning=FALSE, fig.width=11, fig.height=5, echo = TRUE, eval = TRUE}
plot(flows_data)
```
---

## What else is available in `tidyhydat`?
```{r, echo = TRUE, eval = TRUE}
search_stn_name("fraser")
```

---
class: inverse, left, middle
# tidyhydat code break

---
## weathercan
```{r, echo = TRUE, eval = TRUE}
vic_gonzales <- weather_dl(station_ids = "114", interval = "day", start = "2019-01-01", end = "2019-01-31")
vic_gonzales
```
---
## What else is available in `weathercan`?
.Large[
- See `help(package = "weathercan")`
- Normals
- Climate normals measurements
- Station metadata
]

---
class: inverse, left, middle
# weathercan code break

---
## Path of Hurricane Dorian
```{r, cache=TRUE, echo=FALSE}
## Download the data
zip_path <- tempfile()

download.file("https://www.nhc.noaa.gov/gis/best_track/al052019_best_track.zip",
              destfile = zip_path)
unzip(zip_path, exdir = "data/dorian")

## Read in
dorian <- read_sf("data/dorian/AL052019_pts.shp") %>% 
  st_transform(3347)
north_america = ne_countries(continent = "North America", returnclass = "sf") %>% 
  st_transform(3347)
```

```{r, fig.width=13}
ggplot() +
  geom_sf(data = north_america, aes(fill = sovereignt), alpha = 0.3) +
  geom_sf(data = dorian, colour = "black") +
  guides(fill = FALSE) +
  theme_void()
```

---
## Point where Dorian is over Canadian land
```{r, warning=FALSE, fig.width=13}
canada = ne_states(country = "Canada", returnclass = "sf") %>% 
  st_transform(3347)

dorian_canada <- st_intersection(dorian, canada)


ggplot() +
  geom_sf(data = canada, aes(fill = name), alpha = 0.3) +
  geom_sf(data = dorian_canada, colour = "purple", size = 3) +
  guides(fill = FALSE) +
  theme_void()
```

---
## Nova Scotia with buffer

```{r, fig.width=13}

dorian_buffer <- st_buffer(dorian_canada, dist = 7E4)

maritimes <- canada %>% 
  filter(name %in% c("New Brunswick", "Nova Scotia", "Prince Edward Island"))


ggplot() +
  geom_sf(data = maritimes, aes(fill = name), alpha = 0.3) +
  geom_sf(data = dorian_canada, colour = "purple", size = 3) +
  geom_sf(data = dorian_buffer, colour = "orange", alpha = 0.5) +
  guides(fill = FALSE) +
  theme_void()
```

---

## Hydrometric Stations
```{r, warning=FALSE, fig.width=13}
hydro_stations <- realtime_stations() %>% 
  st_as_sf(coords = c("LONGITUDE", "LATITUDE"),
           crs = 4326,
           agr = "constant") %>% 
  st_transform(3347)

climate_stations <- stations %>% 
  filter(end == 2019) %>% 
  filter(interval == "hour") %>% 
  filter(!is.na(lon), !is.na(lat)) %>% 
  st_as_sf(coords = c("lon", "lat"),
           crs = 4326,
           agr = "constant") %>% 
  st_transform(3347)

hydro_dorian <- st_intersection(hydro_stations, dorian_buffer)

climate_dorian <- st_intersection(climate_stations, dorian_buffer)

ggplot() +
  geom_sf(data = maritimes, aes(fill = name), alpha = 0.3) +
  geom_sf(data = dorian_canada, colour = "purple", size = 3) +
  geom_sf(data = dorian_buffer, colour = "orange", alpha = 0.5) +
  geom_sf(data = hydro_dorian, aes(colour = STATION_NUMBER)) +
  guides(fill = FALSE) +
  theme_void()
```

---
## Hydro Data
```{r, fig.width=13, echo = TRUE}
hydro_dorian$STATION_NUMBER
hydro_data <- realtime_dd(station_number = hydro_dorian$STATION_NUMBER) %>% 
  filter(Parameter == "Level") 

hydro_data
```

---
## Hydro Data
```{r, fig.width=15}
dorian_date <- as.POSIXct(paste0(dorian_canada$YEAR, dorian_canada$MONTH, dorian_canada$DAY, " 00:00:01"), "%Y%m%d %H:%M:%S", tz = "GMT")

hydro_data %>% 
  ggplot(aes(x = Date, y = Value, colour = STATION_NUMBER)) +
  geom_line() +
  geom_vline(xintercept = dorian_date, colour = "black", linetype = 2) +
  facet_wrap(~STATION_NUMBER, scales = "free_y") +
  labs(y = "Level (m)") +
  theme_minimal()
```

---
## Climate Stations
```{r, fig.width=13}
ggplot() +
  geom_sf(data = maritimes, aes(fill = name), alpha = 0.3) +
  geom_sf(data = dorian_canada, colour = "purple", size = 3) +
  geom_sf(data = dorian_buffer, colour = "orange", alpha = 0.5) +
  geom_sf(data = climate_dorian, aes(colour = station_name)) +
  guides(fill = FALSE) +
  theme_void()
```

---
## Climate Data
```{r, echo=TRUE, message=FALSE}
climate_dorian$station_id
climate_data <- weather_dl(station_ids = climate_dorian$station_id, 
                           start = "2019-09-01", interval = "hour", quiet = TRUE)

climate_data
```
---

## Climate Data
```{r, fig.width=15, warning=FALSE}
ggplot(climate_data) +
  geom_line(aes(x = time, y = wind_spd, colour = station_name)) +
  geom_vline(xintercept = dorian_date, colour = "black", linetype = 2) +
  facet_wrap(~station_name) +
  labs(y = "Wind Speed (km/h)") +
  theme_minimal()
```


---
## What else is available in R - ggplot2

```{r, warning=FALSE, message=FALSE, echo=TRUE, eval = FALSE, fig.height=5, fig.width=11}
library(ggplot2)

canada_stations <- hy_stations(prov_terr_state_loc = "CA") %>% 
  filter(DRAINAGE_AREA_GROSS < 10000)

ggplot(canada_stations, aes(x = DRAINAGE_AREA_GROSS, fill = HYD_STATUS)) +
  geom_density(alpha = 0.5) +
  labs(x = "Mean long term annual discharge (m^3)", y = "Gross drainage area (km^2)") +
  theme_minimal() +
  facet_wrap(~PROV_TERR_STATE_LOC, scales = "free_y")
```
---

## What else is available in R - ggplot2
```{r, warning=FALSE, message=FALSE, echo=FALSE, eval = TRUE, fig.width=13}
library(ggplot2)

canada_stations <- hy_stations(prov_terr_state_loc = "CA") %>% 
  filter(DRAINAGE_AREA_GROSS < 10000)

ggplot(canada_stations, aes(x = DRAINAGE_AREA_GROSS, fill = HYD_STATUS)) +
  geom_density(alpha = 0.5) +
  labs(x = "Gross drainage area (km^2)") +
  facet_wrap(~PROV_TERR_STATE_LOC, scales = "free_y") +
  theme_minimal()
```


---

## It can be hard!

<center><img src="https://media.giphy.com/media/iKhphMBOk5ScE/giphy.gif" style="width: 650px;"/></center>

---

## Resources for R

<a href = "https://stackoverflow.com/"><img src="https://cdn.sstatic.net/Sites/stackoverflow/company/img/logos/so/so-logo.svg?v=2bb144720a66" alt="Drawing" style="width: 400px;" />


<img src="https://www.rstudio.com/wp-content/uploads/2017/11/logoRStudioCommunity.svg" alt="Drawing" style="width: 400px;" />

<img src="https://www.r-project.org/logo/Rlogo.png" alt="Drawing" style="width: 300px;" />

---
class: inverse, left, middle
## Reprex

> Prepare Reproducible Example Code via the Clipboard

<center><img src="https://raw.githubusercontent.com/tidyverse/reprex/master/man/figures/logo.png" alt="Drawing" style="width: 400px;" /></center>



---
## Contribute to `tidyhydat` and `weathercan`

### Openly developed on GitHub <img src = "https://github.githubassets.com/images/modules/logos_page/GitHub-Mark.png" style = "width: 35px; vertical-align: middle; margin-bottom: 3px;">

<https://github.com/ropensci/tidyhydat>

<https://github.com/ropensci/weathercan>


.pull-left[
Any contribution helps. You don't have to be an R programmer!
- Questions
- Ideas / Feature requests
- Bugs
- Bug-fixes
- Development
]
.pull-right[
<center><img src="https://github.com/ropensci/tidyhydat/raw/master/man/figures/tidyhydat_large.png" alt="Drawing" style="width: 150px;" /></center>
<center><img src="https://raw.githubusercontent.com/ropensci/weathercan/master/inst/assets/weathercan_logo.png" alt="Drawing" style="width: 150px;" /></center>
]

---
## Ways to contribute

- Cite as you would with a paper
- Documentation - write a vignette!
- Use the package - find bugs

### tidyhydat
- SQL code embedded to efficiently do analysis - leverage the database

### weathercan
- Print and plot methods

---

## Ways to cite
<p>📝Albers S (2017).
&ldquo;tidyhydat: Extract and Tidy Canadian Hydrometric Data.&rdquo;
<em>The Journal of Open Source Software</em>, <b>2</b>(20).
doi: <a href="https://doi.org/10.21105/joss.00511">10.21105/joss.00511</a>, <a href="http://dx.doi.org/10.21105/joss.00511">http://dx.doi.org/10.21105/joss.00511</a>. 
</p>


<p>📝LaZerte S, Albers S (2018).
&ldquo;weathercan: Download and format weather data from Environment and Climate Change Canada.&rdquo;
<em>The Journal of Open Source Software</em>, <b>3</b>(22), 571.
<a href="http://joss.theoj.org/papers/10.21105/joss.00571">http://joss.theoj.org/papers/10.21105/joss.00571</a>. 
</p>


```{r, eval=FALSE}
th <- citation("tidyhydat")
print(th, style = "html")
```

```{r, eval=FALSE}
wc <- citation("weathercan")
print(wc, style = "html")
```


<center><img src="https://camo.githubusercontent.com/3e3b4c621878afddfe80f1e22d718ef947292f29/68747470733a2f2f7261776769742e636f6d2f726f70656e7363692f6c6f676f732f6d61737465722f69636f6e5f6c6574746572696e675f636f6c6f722e737667" alt="Drawing" style="width: 600px;" /></center>
---

class: center
## Some Helpful Links

Intro R & RStudio: <https://r4ds.had.co.nz>

Getting started with `tidyhydat`: <https://docs.ropensci.org/tidyhydat>
    
Getting started with `weathercan`: <https://ropensci.github.io/weathercan>

Hydrology CRAN task view: <https://CRAN.R-project.org/view=Hydrology>

rOpenSci: <https://ropensci.org>


    
📝 But we all have to work in excel so read this: 
<https://www.tandfonline.com/doi/full/10.1080/00031305.2017.1375989>

<center><img src="https://cache.desktopnexus.com/cropped-wallpapers/589/589090-1536x864-[DesktopNexus.com].jpg?st=OwEAPaoek3cQGSZM3J0L0w&e=1569384966" style="width: 450px;"/></center>
    
---

class: basic
background-image: url(https://media.giphy.com/media/TnDoEoXfT7YoE/giphy.gif)
background-size: cover

## Questions?
.content-box-blue[
Slides available from 

.small[
https://github.com/ropensci/tidyhydat/blob/master/presentations/tidyhydat_weathercan/tidyhydat_weather.pdf
https://github.com/ropensci/tidyhydat/blob/master/presentations/tidyhydat_weathercan/tidyhydat_weather.Rmd
 ]
Contact <sam.albers@gov.bc.ca>
]


---
title: "tidyhydat: An Introduction"
author: "Sam Albers"
date: "`r Sys.Date()`"
output:
  html_vignette:
     keep_md: true
vignette: >
  %\VignetteIndexEntry{tidyhydat: An Introduction}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r options, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE, 
                      message = FALSE, 
                      eval = nzchar(Sys.getenv("hydat_eval")),
                      fig.width=7, fig.height=7)
```
## Package loading
In addition to tidyhydat, this vignette makes use of the [dplyr](https://dplyr.tidyverse.org/) package for data manipulations and [ggplot2](https://ggplot2.tidyverse.org/) for plotting.
```{r packages, warning=FALSE, message=FALSE, echo = TRUE}
library(tidyhydat)
library(dplyr)
library(ggplot2)
```

# `tidyhydat` package
This vignette will outline a few key options that will hopefully make `tidyhydat` useful.  

## HYDAT download
To use many of the functions in the `tidyhydat` package you will need to download a version of the HYDAT database, Environment and Climate Change Canada's database of historical hydrometric data then tell R where to find the database. Conveniently `tidyhydat` does all this for you via:
```{r, eval=FALSE}
download_hydat()
```
This downloads the most recent version of HYDAT and then saves it in a location on your computer where `tidyhydat`'s function will look for it. Do be patient though as this takes a long time! To see where HYDAT was saved you can run `hy_dir()`. Now that you have HYDAT downloaded and ready to go, you are all set to begin some hydrologic analysis.
  
## Usage
Most functions in `tidyhydat` follow a common argument structure. We will use the `hy_daily_flows()` function for the following examples though the same approach applies to most functions in the package (See `ls("package:tidyhydat")` for a list of exported objects). Much of the functionality of `tidyhydat` originates with the choice of hydrometric stations that you are interested in. A user will often find themselves creating vectors of station numbers. There are several ways to do this. 

The simplest case is if you would like to extract only station. You can supply this directly to the `station_number` argument:
```{r example1, warning=FALSE}
hy_daily_flows(station_number = "08LA001")
```

Another method is to use `hy_stations()` to generate your vector which is then given the `station_number` argument. For example, we could take a subset for only those active stations within Prince Edward Island (Province code:PE) and then create vector for `hy_daily_flows()`:

```{r example2, warning=FALSE}
PEI_stns <- hy_stations() %>%
  filter(HYD_STATUS == "ACTIVE") %>%
  filter(PROV_TERR_STATE_LOC == "PE") %>%
  pull_station_number()

PEI_stns

hy_daily_flows(station_number = PEI_stns)
```

We can also merge our station choice and data extraction into one unified pipe which accomplishes a single goal. For example if for some reason we wanted all the stations in Canada that had the name "Canada" in them we unify that selection and data extraction process into a single pipe:
```{r, example3}
search_stn_name("canada") %>%
  pull_station_number() %>%
  hy_daily_flows()
```

We saw above that if we were only interested in a subset of dates we could use the `start_date` and `end_date` arguments. A date must be supplied to both these arguments in the form of YYYY-MM-DD. If you were interested in all daily flow data from station number "08LA001" for 1981, you would specify all days in 1981 :
```{r warning=FALSE, warning=FALSE, message=FALSE, eval=FALSE}
hy_daily_flows(station_number = "08LA001", 
               start_date = "1981-01-01", end_date = "1981-12-31")
```

This generally outlines the usage of the HYDAT functions within `tidyhydat`. 

## Real-time functions
In addition to the approved and vetted data in the HYDAT database ECCC also offers unapproved data that is subject to revision. `tidyhydat` provides three functions to access these data sources. Remember these are **unapproved** data and should treated as such:

- `realtime_stations()`
- `realtime_dd()`

Not every stations is currently part of the real-time network. Therefore `realtime_stations()` points to a (hopefully) updated ECCC data file of active real-time stations. We can use the `realtime_stations()` functionality to get a vector of stations by jurisdiction. For example, we can choose all the stations in Prince Edward Island using the following:
```{r, eval=FALSE}
realtime_stations(prov_terr_state_loc = "PE")
```

`hy_stations()` and `realtime_stations()` perform similar tasks albeit on different data sources. `hy_stations()` extracts directly from HYDAT. In addition to real-time stations, `hy_stations()` outputs discontinued and non-real-time stations:
```{r stations, eval=FALSE}
hy_stations(prov_terr_state_loc = "PE")
```

This is contrast to `realtime_stations()` which downloads all real-time stations. Though this is not always the case, it is best to use `realtime_stations()` when dealing with real-time data and `hy_stations()` when interacting with HYDAT. It is also appropriate to filter the output of `hy_stations()` by the `REAL_TIME` column.   

### Meterological Service of Canada datamart - `realtime_dd()`
To download real-time data using the datamart we can use approximately the same conventions discussed above. Using `realtime_dd()` we can easily select specific stations by supplying a station of interest:
```{r, eval=FALSE}
realtime_dd(station_number = "08LG006")
```
Another option is to provide simply the province as an argument and download all stations from that province:
```{r, eval=FALSE}
realtime_dd(prov_terr_state_loc = "PE")
```

## Search functions
You can also make use of auxiliary functions in `tidyhydat` called `search_stn_name()` and `search_stn_number()` to look for matches when you know part of a name of a station. For example:
```{r, echo=TRUE}
search_stn_name("liard")
```
Similarly, `search_stn_number()` can be useful if you are interested in all stations from the *08MF* sub-sub-drainage:
```{r, echo=TRUE}
search_stn_number("08MF")
```

## Using joins 
Sometimes it is required to make use of information from two tables from HYDAT. In some cases, we need to combine the information into one table using a common column. Here we will illustrate calculating runoff by combining the `hy_stations` tables with the `hy_daily_flows` table by the `STATION_NUMBER` column:
```{r}
stns <- c("08NH130", "08NH005")
runoff_data <- hy_daily_flows(station_number = stns, start_date = "2000-01-01") %>%
  left_join(
    hy_stations(station_number = stns) %>%
      select(STATION_NUMBER, STATION_NAME, DRAINAGE_AREA_GROSS),
    by = "STATION_NUMBER") %>%
  ## conversion to mm/d
  mutate(runoff = Value / DRAINAGE_AREA_GROSS * 86400 / 1e6 * 1e3) 


ggplot(runoff_data) + 
  geom_line(aes(x = Date, y = runoff, colour = STATION_NAME)) +
  labs(y = "Mean daily runoff [mm/d]") +
  theme_minimal() +
  theme(legend.position = "bottom")
```


# License

    Copyright 2017 Province of British Columbia

    Licensed under the Apache License, Version 2.0 (the "License");
    you may not use this file except in compliance with the License.
    You may obtain a copy of the License at 

       http://www.apache.org/licenses/LICENSE-2.0

    Unless required by applicable law or agreed to in writing, software
    distributed under the License is distributed on an "AS IS" BASIS,
    WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
    See the License for the specific language governing permissions and
    limitations under the License.
---
title: "Two examples of using tidyhydat"
author: "Sam Albers"
date: "`r Sys.Date()`"
output:
  html_vignette:
     keep_md: true
vignette: >
  %\VignetteIndexEntry{Two examples of using tidyhydat}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(echo = TRUE, 
                      warning = FALSE, 
                      message = FALSE, 
                      eval = nzchar(Sys.getenv("hydat_eval")),
                      fig.width=7, fig.height=7)
```

# Introduction

The real power of `tidyhydat` comes with its ability to efficiently import hydrometric data into R which then can be used in conjunction with R's vast array of packages; `tidyhydat` gets hydrometric data into R and allows the user to start conducting data analysis as quickly as possible. This vignette is designed to take a user through two short examples of importing data using `tidyhydat` for use in hydrological visualization and analysis and will therefore use additional tools beyond simply the functions in `tidyhydat`. 

## Familiarity with the tidyverse
`tidyhydat` is designed to exist as a "tidy tool" that is complementary to the [tidyverse](https://CRAN.R-project.org/package=tidyverse) suite of packages. This vignette assumes some familiarity with basic tidy tools particularly `dplyr`, `tidyr`, `ggplot2` and `lubridate`. The [R for Data Science](https://r4ds.had.co.nz/) book by Garrett Grolemund and Hadley Wickham is an enormously useful book both for general R use and for `tidyhydat`. The following is a list of useful links to begin learning particular aspect of each `tidyverse` package:

 - `dplyr`: https://dplyr.tidyverse.org/
 - `tidyr`: https://tidyr.tidyverse.org/
 - `ggplot2`: https://ggplot2.tidyverse.org/
 - `lubridate`: https://lubridate.tidyverse.org/

# Example 1: Basic Data Extraction and Plotting

In this section, we will be using `dplyr`, a data manipulation package, `ggplot2`, an extremely simple yet powerful data visualization tool, `lubridate` a package that aids significantly in working with date and times in R and of course `tidyhydat`. Take a deep breathe. Though this might seem like an overwhelming onslaught of information on packages, each of these packages will save you considerable time with relatively minimal learning investment.
```{r pkg_load_1}
library(tidyhydat)
library(dplyr)
library(ggplot2)
library(lubridate)
```

## Objective
Use `tidyhydat` to find the longest flow record in the Canadian hydrometric network then illustrate the use and ease of `ggplot2` to create some visualizations of the discharge data. 

Your first step is to download the HYDAT database which `tidyhydat` facilitates for you:
```{r dl_hy, eval=FALSE}
download_hydat()
```
You can see where the database was downloaded by running `hy_dir()`. This should be the only instance where you will need to interact directly with HYDAT. Each `tidyhydat` function prefixed with `hy` will automatically know where to look for HYDAT saving you the trouble. When using `tidyhydat`, often your first task is to find the station(s) that you are interested in. Because we are interested in the longest record we can extract that information with `hy_stn_data_range()` and then feed that information to `hy_daily_flows()` like this: 

```{r, eval= FALSE, warning=FALSE, message=FALSE}
longest_record_data <- hy_stn_data_range() %>%
  filter(DATA_TYPE == "Q", RECORD_LENGTH == max(RECORD_LENGTH)) %>%
  pull_station_number() %>%
  hy_daily_flows()
```
Let's break this down line by line to understand how `tidyhydat` uses tidy tools. First we are interested in getting data on record length:
```{r data_range}
hy_stn_data_range()
```
Our objective here is to filter from this data for the station that has the longest record of flow (`DATA_TYPE == "Q"`). You'll also notice this symbol `%>%` which in R is called a [pipe](https://magrittr.tidyverse.org/reference/pipe.html). In code, read it as the word *then*. So for the data_range data we want to grab the data *then* filter it by flow ("Q") in `DATA_TYPE` and then by the maximum value of `RECORD_LENGTH`:
```{r filter}
hy_stn_data_range() %>%
  filter(DATA_TYPE == "Q", RECORD_LENGTH == max(RECORD_LENGTH))
```
*then* pull the `STATION_NUMBER` that has the longest record:
```{r pull}
hy_stn_data_range() %>%
  filter(DATA_TYPE == "Q", RECORD_LENGTH == max(RECORD_LENGTH)) %>%
  pull_station_number()
```
*then* feed that number to `hy_daily_flows()`
```{r, full1}
longest_record_data <- hy_stn_data_range() %>%
  filter(DATA_TYPE == "Q", RECORD_LENGTH == max(RECORD_LENGTH)) %>%
  pull_station_number() %>%
  hy_daily_flows()
```

The result of this collection of simple functions is that we've extracted data the entire daily flow dataset (from 1860-2016) from station `02HA003`.If we wanted more information on this station we could query other tables in HYDAT for more information. The `hy_stations()` function is very useful and outputs considerable metadata on a given station (converted here to a list for viewing purposes):
```{r hy_stns}
hy_stations(station_number = unique(longest_record_data$STATION_NUMBER)) %>%
  as.list()
```


We now know that this station is actually *NIAGARA RIVER AT QUEENSTON* in Ontario and has an upstream drainage basin area of 686,000 square kilometres. As a first step toward visualization, let's simply plot the time series for the entire record with a smoother added:
```{r old_rec}
longest_record_data %>%
  ggplot(aes(x = Date, y = Value)) +
  geom_line() +
  geom_point() +
  geom_smooth() +
  labs(y = "Discharge (m)") +
  theme_minimal()
```
You can see very clearly where continuous monitoring was established. However, this type of plot obscures much of the data that we are interested in. We could plot all the years on the same axis and separate year by line colour. That requires a little bit of manipulation ahead of plotting using both `dplyr` and `lubridate`. Remember you can run a pipe adding one line at a time and comparing the outputs to get a sense of what each step is doing:
```{r old_rec_yr}
longest_record_data %>%
  mutate(dayofyear = yday(Date), Year = year(Date)) %>%
  mutate(dayofyear_formatted = as.Date(dayofyear - 1, origin = "2016-01-01")) %>% ## leap year as placeholder
  ggplot(aes(x = dayofyear_formatted, y = Value, colour = Year)) +
  geom_line() +
  scale_x_date(date_labels = "%b %d") +
  labs(y = "Discharge (m)") +
  theme_minimal()
```

This still is not a very useful plot mostly because our colour range still doesn't resolve year very well and any intra-annual is not visible with this plot. Consider, for example, if rather than a line plot we use a tile plot and modified our colour scale to include more colours:
```{r tile_plt}
longest_record_data %>%
  mutate(dayofyear = yday(Date), Year = year(Date)) %>%
  mutate(dayofyear_formatted = as.Date(dayofyear - 1, origin = "2016-01-01")) %>% 
  ggplot(aes(x = dayofyear_formatted, y = Year, fill = Value)) +
  geom_tile() +
  scale_x_date(date_labels = "%b") +
  scale_y_reverse(expand = c(0, 0)) +
  scale_fill_viridis_c(name = "Discharge (m^3/s) ") +
  labs(y = "Year", x = "Date") +
  theme_minimal() +
  theme(legend.position="bottom")
```
This type of plot provides a clear indicator of discharge patterns over time. Our analysis ends here but one can imagine drawing in climate indices and variation to help explain periods of extreme discharge. 

# Example 2: Distribution of Historical Data
## Objective
Evaluate current conditions in the station from Nunavut with longest record that is also in the Reference Hydrometric Basin Network relative to historical values.

For our second example, we will illustrate how to calculate a common metric used in hydrology: percentiles. For this approach we will use another piped series of functions to zero in on the station that we want. The packages we need for this example are loaded first:
```{r pkg_load_2}
library(tidyhydat)
library(dplyr)
library(ggplot2)
library(lubridate)
```

The pipe below finds the active station that has realtime data in the [Reference Hydrometric Basin Network](https://www.canada.ca/en/environment-climate-change/services/water-overview/quantity/monitoring/survey/data-products-services/reference-hydrometric-basin-network.html) (RHBN) in the territory of Nunavut that has the longest record then grabs all the daily flow information:
  
```{r}
nunavut_stn_flows <- hy_stations() %>%
  filter(HYD_STATUS == "ACTIVE") %>%
  filter(REAL_TIME == TRUE) %>%
  filter(RHBN == TRUE) %>%
  filter(PROV_TERR_STATE_LOC == "NU") %>%
  pull_station_number() %>%
  hy_stn_data_range() %>%
  filter(DATA_TYPE == "Q") %>% 
  filter(RECORD_LENGTH == max(RECORD_LENGTH)) %>%
  pull_station_number() %>%
  hy_daily_flows()
``` 

Armed with our data, we can now evaluate the historical distribution of data then to put some real time data into historical context. First we need to calculate where an individual observation is distributed against all other observations on that day for every day. The `prctile` column tells us what percentage of values on that same day over the entire data record fall above and below the observations on that row. We are restricting our analysis to the last thirty days from when this vignette was last compiled. The `ecdf()` function creates an equation to calculate percentiles based on the `Value` column (i.e discharge) then takes each individual observation of `Value` and calculates the percentile. 
```{r}
pct_flow <- nunavut_stn_flows %>%
  mutate(dayofyear = yday(Date), Year = year(Date)) %>%
  filter(dayofyear %in% yday(seq.Date(from = (Sys.Date()-30), 
                                      to = Sys.Date(), by = "day"))) %>%
  group_by(dayofyear) %>%
  mutate(prctile = ecdf(Value)(Value)) %>%
  mutate(Date_no_year = dmy(paste0(day(Date),"-",month(Date),"-",year(Sys.Date())))) %>%
  ungroup()
```

To collect real time data, we can use the `realtime_dd()` function in `tidyhydat`. Because real time data is collected in hourly (or higher) intervals and we are operating on a daily basis, we first need to take a daily mean of discharge. Again this is accomplished in a pipe:
```{r realtime}
nunavut_realtime <- realtime_dd(unique(nunavut_stn_flows$STATION_NUMBER)) %>%
  mutate(Date_day = as.Date(Date)) %>%
  group_by(Date_day) %>%
  summarise(Value = mean(Value, na.rm = TRUE)) %>%
  ungroup()
```

Finally we can plot all of this on to one figure. This type of plot is useful to combine historical data sources with real time data and provide a visual assessment of how typical flows in a given system are. 

**CAUTION**: All real time data are presented AS IS and represents unapproved data.
```{r, pcrtile_plt}
ggplot(pct_flow, aes(x = Date_no_year, y = Value)) +
  geom_point(aes(colour = prctile)) +
  geom_line(data = nunavut_realtime, aes(x = Date_day), colour = "black") +
  geom_point(data = nunavut_realtime, aes(x = Date_day, shape = factor(year(Date_day))), colour = "black") +
  scale_colour_viridis_c(name = "Discharge Percentile") +
  scale_shape_discrete(name = "Year") +
  theme_minimal() +
  labs(title = "Historical flow relative to current year",
       subtitle = "Current year flows are displayed in black",
       caption = "Real time data is presents AS IS and represents unapproved data",
       x = "Date", y = "Discharge (m^3/s)")
```
---
title: "Stepping into the HYDAT Database"
author: "Dewey Dunnington"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Stepping into the HYDAT Database}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
library(tidyhydat)
library(dplyr)
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

The HYDAT database is a massive hydrologic data resource. The functions in this package are designed to get the most out of the HYDAT database as quickly as possible, however in the process of reformatting the data to be more useful, the package modifies the original tables. This vignette is intended to demonstrate how to access tables within the database for your own custom HYDAT analysis, should additional information be needed. This should not be necessary for the vast majority of users, and is intended only for advanced R users.

## Downloading HYDAT

Before loading the HYDAT database, the latest version of the database must be downloaded using `hydat_download()`. This is a fairly lengthy operation (the download is around 1 GB) and may require several cups of coffee worth of your time.

```{r, eval = FALSE}
hydat_download()
```

```{r, include=FALSE}
# we are actually going to use the test database
# so the vignette can be reproducibly rebuilt without
# needing to call hydat_download
prev_default <- hy_set_default_db(hy_test_db())
```

## Working with HYDAT tables

The HYDAT database is a SQLite database, which can be accessed in R using the [dplyr](https://dplyr.tidyverse.org/) and [dbplyr](https://dbplyr.tidyverse.org/) packages. This package has simplified the connection process, so all you have to do to connect to the database is use `hy_src()`.

```{r}
src <- hy_src()
```

To list the tables, use `src_tbls()` from the **dplyr** package.

```{r}
src_tbls(src)
```

To inspect any particular table, use the `tbl()` function with the `src` and the table name.

```{r}
tbl(src, "STN_OPERATION_SCHEDULE")
```

Working with SQL tables in dplyr is much like working with regular data frames, except no data is actually read from the database until necessary. Because some of these tables are large (particularly those containing the actual data), you will want to `filter()` the tables before you `collect()` them (the `collect()` operation loads them into memory as a `data.frame`).

```{r}
tbl(src, "STN_OPERATION_SCHEDULE") %>%
  filter(STATION_NUMBER == "05AA008") %>%
  collect()
```

When you are finished with the database (i.e., the end of the script), it is good practice to close the connection (you may get a loud red warning if you don't!).

```{r}
hy_src_disconnect(src)
```

```{r, include=FALSE}
# set the default location back to whatever it was before
tidyhydat:::hy_set_default_db(NULL)
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_monthly_flows.R
\name{hy_monthly_flows}
\alias{hy_monthly_flows}
\title{Extract monthly flows information from the HYDAT database}
\format{
A tibble with 8 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Year}{Year of record.}
\item{Month}{Numeric month value}
\item{Full_Month}{Logical value is there is full record from Month}
\item{No_days}{Number of days in that month}
\item{Sum_stat}{Summary statistic being used.}
\item{Value}{Value of the measurement in m^3/s.}
\item{Date_occurred}{Observation date. Formatted as a Date class. MEAN is a annual summary
and therefore has an NA value for Date.}
}
}
\source{
HYDAT
}
\usage{
hy_monthly_flows(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}
}
\value{
A tibble of monthly flows.
}
\description{
Tidy data of monthly flows information from the monthly_flows HYDAT table. \code{station_number} and
\code{prov_terr_state_loc} can both be supplied. If both are omitted all values from the \code{hy_stations} table are returned.
That is a large vector for \code{hy_monthly_flows}.
}
\examples{
\dontrun{
hy_monthly_flows(station_number = c("02JE013","08MF005"),
  start_date = "1996-01-01", end_date = "2000-01-01")

hy_monthly_flows(prov_terr_state_loc = "PE")
          }

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_agency_list}
\alias{hy_agency_list}
\title{hy_agency_list function}
\source{
HYDAT
}
\usage{
hy_agency_list(hydat_path = NULL)
}
\arguments{
\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}
}
\value{
A tibble of agencies
}
\description{
AGENCY look-up Table
}
\examples{
\dontrun{
hy_agency_list()
}

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_monthly_levels.R
\name{hy_monthly_levels}
\alias{hy_monthly_levels}
\title{Extract monthly levels information from the HYDAT database}
\format{
A tibble with 8 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Year}{Year of record.}
\item{Month}{Numeric month value}
\item{Full_month}{Logical value is there is full record from Month}
\item{No_days}{Number of days in that month}
\item{Sum_stat}{Summary statistic being used.}
\item{Value}{Value of the measurement in metres.}
\item{Date_occurred}{Observation date. Formatted as a Date class. MEAN is a annual summary
and therefore has an NA value for Date.}
}
}
\source{
HYDAT
}
\usage{
hy_monthly_levels(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}
}
\value{
A tibble of monthly levels.
}
\description{
Tidy data of monthly river or lake levels information from the DLY_LEVELS HYDAT table. \code{station_number} and
\code{prov_terr_state_loc} can both be supplied. If both are omitted all values from the \code{hy_stations} table are returned.
That is a large vector for \code{hy_monthly_levels}.
}
\examples{
\dontrun{
hy_monthly_levels(station_number = c("02JE013","08MF005"), 
  start_date = "1996-01-01", end_date = "2000-01-01")

hy_monthly_levels(prov_terr_state_loc = "PE")
          }
}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_stn_datum_unrelated}
\alias{hy_stn_datum_unrelated}
\title{Extract station datum unrelated from HYDAT database}
\format{
A tibble with 4 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{DATUM_ID}{Unique code identifying a datum}
\item{Year_from}{First year of use}
\item{Year_to}{Last year of use}
}
}
\usage{
hy_stn_datum_unrelated(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of hy_stn_datum_unrelated
}
\description{
hy_stn_datum_unrelated look-up Table
}
\examples{
\dontrun{
hy_stn_datum_unrelated()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realtime.R
\name{realtime_stations}
\alias{realtime_stations}
\title{Download a tibble of active realtime stations}
\format{
A tibble with 6 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{STATION_NAME}{Official name for station identification}
\item{LATITUDE}{North-South Coordinates of the gauging station in decimal degrees}
\item{LONGITUDE}{East-West Coordinates of the gauging station in decimal degrees}
\item{PROV_TERR_STATE_LOC}{The province, territory or state in which the station is located}
\item{TIMEZONE}{Timezone of the station}
}
}
\usage{
realtime_stations(prov_terr_state_loc = NULL)
}
\arguments{
\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\description{
An up to date dataframe of all stations in the Realtime Water Survey of Canada
hydrometric network operated by Environment and Climate Change Canada
}
\examples{
\dontrun{
## Available inputs for prov_terr_state_loc argument:
unique(realtime_stations()$prov_terr_state_loc)

realtime_stations(prov_terr_state_loc = "BC")
realtime_stations(prov_terr_state_loc = c("QC","PE"))
}
}
\seealso{
Other realtime functions: 
\code{\link{realtime_dd}()}
}
\concept{realtime functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download.R
\name{download_hydat}
\alias{download_hydat}
\title{Download and set the path to HYDAT}
\usage{
download_hydat(dl_hydat_here = NULL, ask = TRUE)
}
\arguments{
\item{dl_hydat_here}{Directory to the HYDAT database. The path is chosen by the \code{rappdirs} package and is OS specific and can be view by \code{\link[=hy_dir]{hy_dir()}}.
This path is also supplied automatically to any function that uses the HYDAT database. A user specified path can be set though this is not the advised approach.
It also downloads the database to a directory specified by \code{\link[=hy_dir]{hy_dir()}}.}

\item{ask}{Whether to ask (as \code{TRUE}/\code{FALSE}) if HYDAT should be downloaded. If \code{FALSE} the keypress question is skipped.}
}
\description{
Download the HYDAT sqlite database. This database contains all the historical hydrometric data for Canada's integrated hydrometric network.
The function will check for a existing sqlite file and won't download the file if the same version is already present.
}
\examples{
\dontrun{
download_hydat()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_annual_stats.R
\name{hy_annual_stats}
\alias{hy_annual_stats}
\title{Extract annual statistics information from the HYDAT database}
\format{
A tibble with 8 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Parameter}{Parameter being measured. Only possible values are FLOW and LEVEL}
\item{Year}{Year of record.}
\item{Sum_stat}{Summary statistic being used.}
\item{Value}{Value of the measurement. If Parameter equals FLOW the units are m^3/s. If Parameter equals LEVEL the
units are metres.}
\item{Date}{Observation date. Formatted as a Date class. MEAN is a annual summary
and therefore has an NA value for Date.}
\item{Symbol}{Measurement/river conditions}
}
}
\source{
HYDAT
}
\usage{
hy_annual_stats(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_year = "ALL",
  end_year = "ALL"
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_year}{First year of the returned record}

\item{end_year}{Last year of the returned record}
}
\value{
A tibble of hy_annual_stats.
}
\description{
Provides wrapper to turn the ANNUAL_STATISTICS table in HYDAT into a tidy data frame of annual statistics.
Statistics provided include MEAN, MAX and MIN on an annual basis.
}
\examples{
\dontrun{
  ## Multiple stations province not specified
  hy_annual_stats(station_number = c("08NM083","05AE027"))

  ## Multiple province, station number not specified
  hy_annual_stats(prov_terr_state_loc = c("AB","SK"))
}

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_sed_samples.R
\name{hy_sed_samples}
\alias{hy_sed_samples}
\title{Extract instantaneous sediment sample information from the HYDAT database}
\format{
A tibble with 19 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{SED_DATA_TYPE}{Contains the type of sampling method used in collecting sediment for a station}
\item{Date}{Contains the time to the nearest minute of when the sample was taken}
\item{SAMPLE_REMARK_CODE}{Descriptive Sediment Sample Remark in English}
\item{TIME_SYMBOL}{An "E" symbol means the time is an estimate only}
\item{FLOW}{Contains the instantaneous discharge in cubic metres per second at the time the sample was taken}
\item{SYMBOL_EN}{Indicates a condition where the daily mean has a larger than expected error}
\item{SAMPLER_TYPE}{Contains the type of measurement device used to take the sample}
\item{SAMPLING_VERTICAL_LOCATION}{The location on the cross-section of the river
at which the single sediment samples are collected. If one of the standard
locations is not used the distance in meters will be shown}
\item{SAMPLING_VERTICAL_EN}{Indicates sample location relative to the
regular measurement cross-section or the regular sampling site}
\item{TEMPERATURE}{Contains the instantaneous water temperature
in Celsius at the time the sample was taken}
\item{CONCENTRATION_EN}{Contains the instantaneous concentration sampled in milligrams per litre}
\item{SV_DEPTH2}{Depth 2 for split vertical depth integrating (m)}
}
}
\source{
HYDAT
}
\usage{
hy_sed_samples(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}
}
\value{
A tibble of instantaneous sediment samples data
}
\description{
Provides wrapper to turn the hy_sed_samples table in HYDAT into a tidy data frame of instantaneous sediment sample information.
\code{station_number} and \code{prov_terr_state_loc} can both be supplied. If both are omitted all values from the \code{hy_stations}
table are returned. That is a large vector for \code{hy_sed_samples}.
}
\examples{
\dontrun{
hy_sed_samples(station_number = "01CA004")
          }

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_sed_daily_suscon.R
\name{hy_sed_daily_suscon}
\alias{hy_sed_daily_suscon}
\title{Extract daily suspended sediment concentration information from the HYDAT database}
\format{
A tibble with 5 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Date}{Observation date. Formatted as a Date class.}
\item{Parameter}{Parameter being measured. Only possible value is Suscon}
\item{Value}{Discharge value. The units are mg/l.}
\item{Symbol}{Measurement/river conditions}
}
}
\source{
HYDAT
}
\usage{
hy_sed_daily_suscon(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL,
  symbol_output = "code"
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{symbol_output}{Set whether the raw code, or the \code{english} or the \code{french} translations are outputted. Default
value is \code{code}.}
}
\value{
A tibble of daily suspended sediment concentration
}
\description{
Provides wrapper to turn the SED_DLY_SUSCON table in HYDAT into a tidy data frame of daily suspended sediment concentration information.
\code{station_number} and \code{prov_terr_state_loc} can both be supplied. If both are omitted all values from the \code{hy_stations}
table are returned. That is a large vector for \code{hy_sed_daily_suscon}.
}
\examples{
\dontrun{
hy_sed_daily_suscon(station_number = "01CE003")
          }

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_sed_monthly_suscon.R
\name{hy_sed_monthly_suscon}
\alias{hy_sed_monthly_suscon}
\title{Extract monthly flows information from the HYDAT database}
\format{
A tibble with 8 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Year}{Year of record.}
\item{Month}{Numeric month value}
\item{Full_Month}{Logical value is there is full record from Month}
\item{No_days}{Number of days in that month}
\item{Sum_stat}{Summary statistic being used.}
\item{Value}{Value of the measurement in mg/l.}
\item{Date_occurred}{Observation date. Formatted as a Date class. MEAN is a annual summary
and therefore has an NA value for Date.}
}
}
\source{
HYDAT
}
\usage{
hy_sed_monthly_suscon(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}
}
\value{
A tibble of monthly suspended sediment concentrations.
}
\description{
Tidy data of monthly suspended sediment concentration information from the SED_DLY_SUSCON HYDAT table.  \code{station_number} and
\code{prov_terr_state_loc} can both be supplied. If both are omitted all values from the \code{hy_stations} table are returned.
That is a large vector for \code{hy_sed_monthly_suscon}.
}
\examples{
\dontrun{
hy_sed_monthly_suscon(station_number = "08MF005")
          }
}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-search.R
\name{search_stn_name}
\alias{search_stn_name}
\alias{search_stn_number}
\title{A search function for hydrometric station name or number}
\usage{
search_stn_name(search_term, hydat_path = NULL)

search_stn_number(search_term, hydat_path = NULL)
}
\arguments{
\item{search_term}{Only accepts one word.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}
}
\value{
A tibble of stations that match the \code{search_term}
}
\description{
Use this search function when you only know the partial station name or want to search.
}
\examples{
\dontrun{
search_stn_name("Cowichan")

search_stn_number("08HF")
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{allstations}
\alias{allstations}
\title{All Canadian stations}
\format{
A tibble with 5 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{STATION_NAME}{Official name for station identification}
\item{PROV_TERR_STATE_LOC}{The province, territory or state in which the station is located}
\item{HYD_STATUS}{Current status of discharge or level monitoring in the hydrometric network}
\item{REAL_TIME}{Logical. Indicates if a station has the capacity to deliver data in
real-time or near real-time}
\item{LATITUDE}{North-South Coordinates of the gauging station in decimal degrees}
\item{LONGITUDE}{East-West Coordinates of the gauging station in decimal degrees}
\item{station_tz}{Timezone of station calculated using the lutz package based on LAT/LONG of stations}
\item{standard_offset}{Offset from UTC of local standard time}
}
}
\source{
HYDAT, Meteorological Service of Canada datamart
}
\usage{
allstations
}
\description{
A shorthand to avoid having always call \code{hy_stations} or \code{realtime_stations}.
Populated by both realtime and historical data from HYDAT.
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realtime.R
\name{realtime_dd}
\alias{realtime_dd}
\title{Download a tibble of realtime river data from the last 30 days from the Meteorological Service of Canada datamart}
\format{
A tibble with 8 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{PROV_TERR_STATE_LOC}{The province, territory or state in which the station is located}
\item{Date}{Observation date and time for last thirty days. Formatted as a POSIXct class in UTC for consistency.}
\item{Parameter}{Parameter being measured. Only possible values are Flow and Level}
\item{Value}{Value of the measurement. If Parameter equals Flow the units are m^3/s.
If Parameter equals Level the units are metres.}
\item{Grade}{reserved for future use}
\item{Symbol}{reserved for future use}
\item{Code}{quality assurance/quality control flag for the discharge}
\item{station_tz}{Station timezone based on tidyhydat::allstations$station_tz}
}
}
\usage{
realtime_dd(station_number = NULL, prov_terr_state_loc = NULL)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of water flow and level values. The date and time of the query (in UTC) is also
stored as an attribute.
}
\description{
Download realtime river data from the last 30 days from the Meteorological Service of Canada (MSC) datamart.
The function will prioritize downloading data collected at the highest resolution. In instances where data is
not available at high (hourly or higher) resolution daily averages are used. Currently, if a station does not
exist or is not found, no data is returned.
}
\examples{
\dontrun{
## Download from multiple provinces
realtime_dd(station_number=c("01CD005","08MF005"))

## To download all stations in Prince Edward Island:
pei <- realtime_dd(prov_terr_state_loc = "PE")

## Access the time of query
attributes(pei)$query_time
}

}
\seealso{
Other realtime functions: 
\code{\link{realtime_stations}()}
}
\concept{realtime functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_sed_samples_psd.R
\name{hy_sed_samples_psd}
\alias{hy_sed_samples_psd}
\title{Extract instantaneous sediment sample particle size distribution information from the HYDAT database}
\format{
A tibble with 5 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{SED_DATA_TYPE}{Contains the type of sampling method used in collecting sediment for a station}
\item{Date}{Contains the time to the nearest minute of when the sample was taken}
\item{PARTICLE_SIZE}{Particle size (mm)}
\item{PERCENT}{Contains the percentage values for indicated particle sizes for samples collected}
}
}
\source{
HYDAT
}
\usage{
hy_sed_samples_psd(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}
}
\value{
A tibble of sediment sample particle size data
}
\description{
Provides wrapper to turn the hy_sed_samples_psd table in HYDAT into a tidy data frame of instantaneous sediment sample
particle size distribution.  \code{station_number} and \code{prov_terr_state_loc} can both be supplied. If both
are omitted all values from the \code{\link[=hy_stations]{hy_stations()}} table are returned. That is a large vector for \code{hy_sed_samples_psd}.
}
\examples{
\dontrun{
hy_sed_samples_psd(station_number = "01CA004")
          }

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_reg_office_list}
\alias{hy_reg_office_list}
\title{Extract regional office list from HYDAT database}
\source{
HYDAT
}
\usage{
hy_reg_office_list(hydat_path = NULL)
}
\arguments{
\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}
}
\value{
A tibble of offices
}
\description{
OFFICE look-up Table
}
\examples{
\dontrun{
hy_reg_office_list()
}


}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_stn_data_range}
\alias{hy_stn_data_range}
\title{Extract station data range from HYDAT database}
\format{
A tibble with 6 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{DATA_TYPE}{Code for the type of data}
\item{SED_DATA_TYPE}{Code for the type of instantaneous sediment data}
\item{Year_from}{First year of use}
\item{Year_to}{Last year of use}
\item{RECORD_LENGTH}{Number of years of data available in the HYDAT database}
}
}
\source{
HYDAT
}
\usage{
hy_stn_data_range(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of hy_stn_data_range
}
\description{
hy_stn_data_range look-up Table
}
\examples{
\dontrun{
hy_stn_data_range(station_number = c("02JE013","08MF005"))
}

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_sed_monthly_loads.R
\name{hy_sed_monthly_loads}
\alias{hy_sed_monthly_loads}
\title{Extract monthly flows information from the HYDAT database}
\format{
A tibble with 8 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Year}{Year of record.}
\item{Month}{Numeric month value}
\item{Full_Month}{Logical value is there is full record from Month}
\item{No_days}{Number of days in that month}
\item{Sum_stat}{Summary statistic being used.}
\item{Value}{Value of the measurement in tonnes.}
\item{Date_occurred}{Observation date. Formatted as a Date class. MEAN is a annual summary
and therefore has an NA value for Date.}
}
}
\source{
HYDAT
}
\usage{
hy_sed_monthly_loads(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}
}
\value{
A tibble of monthly sediment loads.
}
\description{
Tidy data of monthly loads information from the SED_DLY_LOADS HYDAT table. \code{station_number} and
\code{prov_terr_state_loc} can both be supplied. If both are omitted all values from the \code{hy_stations} table are returned.
That is a large vector for \code{hy_sed_monthly_loads}.
}
\examples{
\dontrun{
hy_sed_monthly_loads(station_number = "01CE003")
          }
          
}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{pull_station_number}
\alias{pull_station_number}
\title{Convenience function to pull station number from tidyhydat functions}
\usage{
pull_station_number(.data)
}
\arguments{
\item{.data}{A table of data}
}
\value{
A vector of station_numbers
}
\description{
This function mimics \code{dplyr::pull} to avoid having to always type
dplyr::pull(STATION_NUMBER). Instead we can now take advantage of autocomplete.
This can be used with \code{realtime_} and \code{hy_} functions.
}
\examples{
\dontrun{

hy_stations(prov_terr_state_loc = "PE") \%>\%
 pull_station_number() \%>\%
 hy_annual_instant_peaks()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realtime_plot.R
\name{realtime_plot}
\alias{realtime_plot}
\title{Convenience function to plot realtime data}
\usage{
realtime_plot(station_number = NULL, Parameter = c("Flow", "Level"))
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. Can only be one value.}

\item{Parameter}{Parameter of interest. Either "Flow" or "Level". Defaults to "Flow".}
}
\value{
A plot of recent realtime values
}
\description{
This is an easy way to visualize a single station using base R graphics.
More complicated plotting needs should consider using \code{ggplot2}. Inputting more
5 stations will result in very busy plots and longer load time. Legend position will
sometimes overlap plotted points.
}
\examples{
\dontrun{
## One station
realtime_plot("08MF005")

## Multiple stations
realtime_plot(c("07EC002","01AD003"))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_db.R
\name{hy_set_default_db}
\alias{hy_set_default_db}
\title{Set the default database path}
\usage{
hy_set_default_db(hydat_path = NULL)
}
\arguments{
\item{hydat_path}{The path to the a HYDAT sqlite3 database file
(e.g., \link{hy_test_db})}
}
\value{
returns the previous value of \link{hy_default_db}.
}
\description{
For many reasons, it may be convenient to set the default
database location to somewhere other than the global default. Users
may wish to use a previously downloaded version of the database for
reproducibility purposes, store hydat somewhere other than hy_dir().
}
\examples{
\dontrun{
# set default to the test database
hy_set_default_db(hy_test_db())

# get the default value
hy_default_db()

# set back to the default db location
hy_set_default_db(NULL)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-tidy-eval.R
\name{tidyeval}
\alias{tidyeval}
\alias{quo}
\alias{quos}
\alias{enquo}
\alias{sym}
\alias{syms}
\alias{ensym}
\alias{expr}
\alias{exprs}
\alias{enexpr}
\alias{quo_name}
\title{Tidy eval helpers}
\description{
These functions provide tidy eval-compatible ways to capture
symbols (\code{sym()}, \code{syms()}, \code{ensym()}), expressions (\code{expr()},
\code{exprs()}, \code{enexpr()}), and quosures (\code{quo()}, \code{quos()}, \code{enquo()}).
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_stations.R
\name{hy_stations}
\alias{hy_stations}
\title{Extract station information from the HYDAT database}
\format{
A tibble with 15 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{STATION_NAME}{Official name for station identification}
\item{PROV_TERR_STATE_LOC}{The province, territory or state in which the station is located}
\item{REGIONAL_OFFICE_ID}{The identifier of the regional office responsible for the station.
Links to \link[tidyhydat]{hy_reg_office_list}}
\item{HYD_STATUS}{Current status of discharge or level monitoring in the hydrometric network}
\item{SED_STATUS}{Current status of sediment monitoring in the hydrometric network}
\item{LATITUDE}{North-South Coordinates of the gauging station in decimal degrees}
\item{LONGITUDE}{East-West Coordinates of the gauging station in decimal degrees}
\item{DRAINAGE_AREA_GROSS}{The total surface area that drains to the gauge site (km^2)}
\item{DRAINAGE_AREA_EFFECT}{The portion of the drainage basin that contributes runoff to
the gauge site, calculated by subtracting any noncontributing portion from the
gross drainage area (km^2)}
\item{RHBN}{Logical. Reference Hydrometric Basin Network station. The Reference Hydrometric
Basin Network (RHBN) is a sub-set of the national network that has been identified
for use in the detection, monitoring, and assessment of climate change.}
\item{REAL_TIME}{Logical. Indicates if a station has the capacity to deliver data in
real-time or near real-time}
\item{CONTRIBUTOR_ID}{Unique ID of an agency that contributes data to the
HYDAT database. The agency is non-WSC and non WSC funded}
\item{OPERATOR_ID}{Unique ID of an agency that operates a hydrometric station}
\item{DATUM_ID}{Unique ID for a datum}
}
}
\source{
HYDAT
}
\usage{
hy_stations(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of stations and associated metadata
}
\description{
Provides wrapper to turn the hy_stations table in HYDAT into a tidy data frame of station information. \code{station_number} and
\code{prov_terr_state_loc} can both be supplied. If both are omitted all values from the \code{hy_stations} table are returned. This
is the entry point for most analyses is tidyhydat as establish the stations for consideration is likely the first step in many
instances.
}
\examples{
\dontrun{
## Multiple stations province not specified
hy_stations(station_number = c("08NM083","08NE102"))

## Multiple province, station number not specified
hy_stations(prov_terr_state_loc = c("AB","YT"))
}


}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{hy_data_types}
\alias{hy_data_types}
\title{DATA TYPES look-up table}
\format{
A tibble with 5 rows and 3 variables:
\describe{
\item{DATA_TYPE}{Data type code}
\item{DATA_TYPE_EN}{Descriptive data type (English)}
\item{DATA_TYPE_FR}{Descriptive data type (French)}
}
}
\source{
HYDAT
}
\usage{
hy_data_types
}
\description{
A look table for data types
}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tidyhydat-package.R
\docType{package}
\name{tidyhydat-package}
\alias{tidyhydat}
\alias{tidyhydat-package}
\title{tidyhydat: Extract and Tidy Canadian 'Hydrometric' Data}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Provides functions to access historical and real-time national 'hydrometric' data from Water Survey of Canada data sources (<https://dd.weather.gc.ca/hydrometric/csv/> and <https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/>) and then applies tidy data principles.
}
\references{
To download the latest version of hydat please:
\itemize{
\item use the \code{\link[=download_hydat]{download_hydat()}} function.
\item If that fails you can download directly from this link: \url{https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/}
}

For more information on tidy data please see
\itemize{
\item Wickham, Hadley. 2014. Tidy Data. The Journal of Statistical Software. 59. \doi{10.18637/jss.v059.i10}
\item tidy data vignette: \url{https://CRAN.R-project.org/package=tidyr/vignettes/tidy-data.html}
}

For more information on HYDAT and ECCC data sources
\itemize{
\item Please see this description of the \href{https://www.canada.ca/en/environment-climate-change/services/water-overview/quantity/monitoring/survey/data-products-services/national-archive-hydat.html}{database}
\item This page is landing page for technical description of \href{https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/}{HYDAT}
\item This pdf links to a document that outlines database table \href{https://collaboration.cmc.ec.gc.ca/cmc/hydrometrics/www/HYDAT_Definition_EN.pdf}{definitions}
\item This FAQ provides a helpful list of ECCC data source \href{https://wateroffice.ec.gc.ca/contactus/faq_e.html}{questions}
}
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/tidyhydat/}
  \item \url{https://github.com/ropensci/tidyhydat/}
  \item Report bugs at \url{https://github.com/ropensci/tidyhydat/issues/}
}

}
\author{
\strong{Maintainer}: Sam Albers \email{sam.albers@gov.bc.ca} (\href{https://orcid.org/0000-0002-9270-7884}{ORCID})

Other contributors:
\itemize{
  \item David Hutchinson \email{david.hutchinson@canada.ca} [contributor]
  \item Dewey Dunnington \email{dewey@fishandwhistle.net} [contributor]
  \item Ryan Whaley \email{rdgwhaley@gmail.com} [contributor]
  \item Province of British Columbia [copyright holder]
  \item Government of Canada [data contributor]
  \item Luke Winslow (Reviewed for rOpenSci) [reviewer]
  \item Laura DeCicco (Reviewed for rOpenSci) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_stn_regulation.R
\name{hy_stn_regulation}
\alias{hy_stn_regulation}
\title{Extract station regulation from the HYDAT database}
\format{
A tibble with 4 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Year_from}{First year of use}
\item{Year_to}{Last year of use}
\item{REGULATED}{logical}
}
}
\source{
HYDAT
}
\usage{
hy_stn_regulation(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of stations, years of regulation and the regulation status
}
\description{
Provides wrapper to turn the hy_stn_regulation table in HYDAT into a tidy data frame of station regulation.
\code{station_number} and \code{prov_terr_state_loc} can both be supplied. If both are omitted all values
from the \code{hy_stations} table are returned.
}
\examples{
\dontrun{
## Multiple stations province not specified
hy_stn_regulation(station_number = c("08NM083","08NE102"))

## Multiple province, station number not specified
hy_stn_regulation(prov_terr_state_loc = c("AB","YT"))
}

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_daily_levels.R
\name{hy_daily_levels}
\alias{hy_daily_levels}
\title{Extract daily levels information from the HYDAT database}
\format{
A tibble with 5 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Date}{Observation date. Formatted as a Date class.}
\item{Parameter}{Parameter being measured. Only possible value is Level}
\item{Value}{Level value. The units are metres.}
\item{Symbol}{Measurement/river conditions}
}
}
\source{
HYDAT
}
\usage{
hy_daily_levels(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL,
  symbol_output = "code"
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{symbol_output}{Set whether the raw code, or the \code{english} or the \code{french} translations are outputted. Default
value is \code{code}.}
}
\value{
A tibble of daily levels
}
\description{
Provides wrapper to turn the DLY_LEVELS table in HYDAT into a tidy data frame.  The primary value returned by this
function is discharge. \code{station_number} and \code{prov_terr_state_loc} can both be supplied. If both are omitted all
values from the \code{hy_stations} table are returned. That is a large vector for \code{hy_daily_levels}.
}
\examples{
\dontrun{
hy_daily_levels(station_number = c("02JE013","08MF005"), 
  start_date = "1996-01-01", end_date = "2000-01-01")

hy_daily_levels(prov_terr_state_loc = "PE")
          }

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_daily.R
\name{hy_daily}
\alias{hy_daily}
\title{Extract all daily water level and flow measurements}
\format{
A tibble with 5 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Date}{Observation date. Formatted as a Date class.}
\item{Parameter}{Parameter being measured.}
\item{Value}{Discharge value. The units are m^3/s.}
\item{Symbol}{Measurement/river conditions}
}
}
\source{
HYDAT
}
\usage{
hy_daily(
  station_number = NULL,
  prov_terr_state_loc = NULL,
  hydat_path = NULL,
  ...
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{...}{See \code{\link[=hy_daily_flows]{hy_daily_flows()}} arguments}
}
\value{
A tibble of daily flows and levels
}
\description{
A thin wrapper around \code{hy_daily_flows} and `hy_daily_levels`` that returns a data frames that
contains both parameters. All arguments are passed directly to these functions.
}
\examples{
\dontrun{
hy_daily(station_number = c("02JE013","08MF005"))
}
}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_stn_data_coll}
\alias{hy_stn_data_coll}
\title{Extract station data collection from HYDAT database}
\format{
A tibble with 6 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{DATA_TYPE}{The type of data}
\item{Year_from}{First year of use}
\item{Year_to}{Last year of use}
\item{MEASUREMENT}{The sampling method used in the collection of
sediment data or the type of the gauge used in the collection of the hydrometric data}
\item{OPERATION}{The schedule of station operation
for the collection of sediment or hydrometric data}
}
}
\source{
HYDAT
}
\usage{
hy_stn_data_coll(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of hy_stn_data_coll
}
\description{
hy_stn_data_coll look-up Table
}
\examples{
\dontrun{
hy_stn_data_coll(station_number = c("02JE013","08MF005"))
}

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_daily_flows.R
\name{hy_daily_flows}
\alias{hy_daily_flows}
\title{Extract daily flows information from the HYDAT database}
\format{
A tibble with 5 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Date}{Observation date. Formatted as a Date class.}
\item{Parameter}{Parameter being measured. Only possible value is Flow}
\item{Value}{Discharge value. The units are m^3/s.}
\item{Symbol}{Measurement/river conditions}
}
}
\source{
HYDAT
}
\usage{
hy_daily_flows(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL,
  symbol_output = "code"
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{symbol_output}{Set whether the raw code, or the \code{english} or the \code{french} translations are outputted. Default
value is \code{code}.}
}
\value{
A tibble of daily flows
}
\description{
Provides wrapper to turn the DLY_FLOWS table in HYDAT into a tidy data frame of daily flows.
\code{station_number} and \code{prov_terr_state_loc} can both be supplied. If both are omitted all
values from the \code{hy_stations} table are returned. That is a large tibble for \code{hy_daily_flows}.
}
\examples{
\dontrun{
#download_hydat()
hy_daily_flows(station_number = c("08MF005"), 
  start_date = "1996-01-01", end_date = "2000-01-01")

hy_daily_flows(prov_terr_state_loc = "PE")
          }

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_db.R
\name{hy_src}
\alias{hy_src}
\alias{hy_src_disconnect}
\title{Open a connection to the HYDAT database}
\usage{
hy_src(hydat_path = NULL)

hy_src_disconnect(src)
}
\arguments{
\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{src}{A  as returned by \code{\link[=hy_src]{hy_src()}}.}
}
\value{
A SQLite DBIConnection
}
\description{
This function gives low-level access to the underlying HYDAT database used by
other functions. Many of these tables are too large to load into memory,
so it is best to use dplyr to \code{\link[dplyr:filter]{dplyr::filter()}} them before using
\code{\link[dplyr:compute]{dplyr::collect()}} to read them into memory.
}
\examples{
\dontrun{
library(dplyr)

# src is a src_sqlite
src <- hy_src(hydat_path = hy_test_db())
src_tbls(src)

# to get a table, use dplyr::tbl()
tbl(src, "STATIONS")

# one you're sure the results are what you want
# get a data.frame using collect()
tbl(src, "STATIONS") \%>\%
  filter(PROV_TERR_STATE_LOC == "BC") \%>\%
  collect()
  
# close the connection to the database
hy_src_disconnect(src)
}
}
\seealso{
\code{\link[=download_hydat]{download_hydat()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_stn_op_schedule}
\alias{hy_stn_op_schedule}
\title{Extract station operation schedule from HYDAT database}
\format{
A tibble with 6 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{DATA_TYPE}{The type of data}
\item{Year}{Year of operation schedule}
\item{Month_from}{First month of use}
\item{Month_to}{Last month of use}
}
}
\source{
HYDAT
}
\usage{
hy_stn_op_schedule(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of hy_stn_op_schedule
}
\description{
hy_stn_op_schedule look-up Table
}
\examples{
\dontrun{
hy_stn_op_schedule(station_number = c("02JE013"))
}

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{hy_data_symbols}
\alias{hy_data_symbols}
\title{DATA SYMBOLS look-up table}
\format{
A tibble with 5 rows and 3 variables:
\describe{
\item{SYMBOL_ID}{Symbol code}
\item{SYMBOL_EN}{Description of Symbol (English)}
\item{SYMBOL_FR}{Description of Symbol (French)}
}
}
\source{
HYDAT
}
\usage{
hy_data_symbols
}
\description{
A look table for data symbols
}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_plot.R
\name{hy_plot}
\alias{hy_plot}
\title{This function is deprecated in favour of generic plot methods}
\usage{
hy_plot(
  station_number = NULL,
  Parameter = c("Flow", "Level", "Suscon", "Load")
)
}
\arguments{
\item{station_number}{A (or several) seven digit Water Survey of Canada station number.}

\item{Parameter}{Parameter of interest. Either "Flow" or "Level".}
}
\description{
This is an easy way to visualize a single station using base R graphics.
More complicated plotting needs should consider using \code{ggplot2}. Inputting more
5 stations will result in very busy plots and longer load time. Legend position will
sometimes overlap plotted points.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{\%>\%}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{dplyr}{\code{\link[dplyr:reexports]{\%>\%}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{station_choice}
\alias{station_choice}
\title{Function to chose a station based on consistent arguments for hydat functions.}
\usage{
station_choice(hydat_con, station_number, prov_terr_state_loc)
}
\arguments{
\item{hydat_con}{A database connection}

\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\description{
A function to avoid duplication in HYDAT functions.  This function is not intended for external use.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_dir}
\alias{hy_dir}
\title{Output OS-independent path to the HYDAT sqlite database}
\usage{
hy_dir(...)
}
\arguments{
\item{...}{arguments potentially passed to \code{rappdirs::user_data_dir}}
}
\description{
Provides the download location for \link{download_hydat} in an OS independent manner.
}
\examples{
\dontrun{
hy_dir()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_sed_daily_loads.R
\name{hy_sed_daily_loads}
\alias{hy_sed_daily_loads}
\title{Extract daily sediment load information from the HYDAT database}
\format{
A tibble with 4 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{Date}{Observation date. Formatted as a Date class.}
\item{Parameter}{Parameter being measured. Only possible value is Load}
\item{Value}{Discharge value. The units are tonnes.}
}
}
\source{
HYDAT
}
\usage{
hy_sed_daily_loads(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_date = NULL,
  end_date = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}

\item{end_date}{Leave blank if all dates are required. Date format needs to be in YYYY-MM-DD. Date is inclusive.}
}
\value{
A tibble of daily suspended sediment loads
}
\description{
Provides wrapper to turn the SED_DLY_LOADS table in HYDAT into a tidy data frame of daily sediment load information.
\code{station_number} and \code{prov_terr_state_loc} can both be supplied. If both are omitted all values from the
\code{hy_stations} table are returned. That is a large vector for \code{hy_sed_daily_loads}.
}
\examples{
\dontrun{
hy_sed_daily_loads(prov_terr_state_loc = "PE")
          }

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_plot.R, R/realtime_plot.R
\name{plot}
\alias{plot}
\alias{plot.hy}
\alias{plot.realtime}
\title{Plot historical and realtime data}
\usage{
\method{plot}{hy}(x = NULL, ...)

\method{plot}{realtime}(x = NULL, Parameter = c("Flow", "Level"), ...)
}
\arguments{
\item{x}{Object created by either a hy_daily_* or realtime_dd data retrieval function}

\item{...}{passed to \code{\link[=plot]{plot()}}}

\item{Parameter}{Parameter of interest. Either "Flow" or "Level". Defaults to "Flow".}
}
\description{
This method plots either daily time series data from HYDAT or realtime data from
the datamart. These plots are intended to be convenient and quick methods to
visualize hydrometric data.
}
\section{Methods (by class)}{
\itemize{
\item \code{realtime}: plot.realtime
}}

\examples{
\dontrun{
# One station
fraser <- hy_daily_flows("08MF005")
plot(fraser)
}

\dontrun{
# One station
fraser_realtime <- realtime_dd("08MF005")
plot(fraser_realtime)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_stn_datum_conv}
\alias{hy_stn_datum_conv}
\title{Extract station datum conversions from HYDAT database}
\format{
A tibble with 4 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{DATUM_FROM}{Identifying a datum from which water level is being converted}
\item{DATUM_TO}{Identifying a datum to which water level is being converted}
\item{CONVERSTION_FACTOR}{The conversion factor applied to water levels referred to
one datum to obtain water levels referred to another datum}
}
}
\usage{
hy_stn_datum_conv(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of hy_stn_datum_conv
}
\description{
hy_stn_datum_conv look-up Table
}
\examples{
\dontrun{
hy_stn_datum_conv(station_number = c("02JE013","08MF005"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_stn_remarks}
\alias{hy_stn_remarks}
\title{Extract station remarks from HYDAT database}
\format{
A tibble with 4 variables:
\describe{
\item{STATION_NUMBER}{Unique 7 digit Water Survey of Canada station number}
\item{REMARK_TYPE}{Type of Remark}
\item{Year}{Year of the remark}
\item{REMARK}{Remark}
}
}
\usage{
hy_stn_remarks(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}
}
\value{
A tibble of hy_stn_remarks
}
\description{
hy_stn_remarks look-up Table
}
\examples{
\dontrun{
hy_stn_remarks(station_number = c("02JE013","08MF005"))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_annual_instant_peaks.R
\name{hy_annual_instant_peaks}
\alias{hy_annual_instant_peaks}
\title{Extract annual max/min instantaneous flows and water levels from HYDAT database}
\source{
HYDAT
}
\usage{
hy_annual_instant_peaks(
  station_number = NULL,
  hydat_path = NULL,
  prov_terr_state_loc = NULL,
  start_year = NULL,
  end_year = NULL
)
}
\arguments{
\item{station_number}{A seven digit Water Survey of Canada station number. If this argument is omitted, the value of \code{prov_terr_state_loc}
is returned.}

\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}

\item{prov_terr_state_loc}{Province, state or territory. If this argument is omitted, the value of \code{station_number}
is returned. See \code{unique(allstations$prov_terr_state_loc)}. Will also accept \code{CA} to return only Canadian stations.}

\item{start_year}{First year of the returned record}

\item{end_year}{Last year of the returned record}
}
\value{
A tibble of hy_annual_instant_peaks.
}
\description{
Provides wrapper to turn the ANNUAL_INSTANT_PEAKS table in HYDAT into a tidy data frame of instantaneous flows and water levels.
\code{station_number} and \code{prov_terr_state_loc} can both be supplied.
}
\examples{
\dontrun{
## Multiple stations province not specified
hy_annual_instant_peaks(station_number = c("08NM083","08NE102"))

## Multiple province, station number not specified
hy_annual_instant_peaks(prov_terr_state_loc = c("AB","YT"))
}

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realtime.R
\name{realtime_add_local_datetime}
\alias{realtime_add_local_datetime}
\title{Add local datetime column to realtime tibble}
\usage{
realtime_add_local_datetime(.data, set_tz = NULL)
}
\arguments{
\item{.data}{Tibble created by \code{realtime_dd}}

\item{set_tz}{A timezone string in the format of \code{OlsonNames()}}
}
\description{
Adds \code{local_datetime} and \code{tz_used} columns based on either the most common timezone in the original data or
a user supplied timezone. This function is meant to used in a pipe with the \code{realtime_dd()} function.
}
\details{
\code{Date} from \code{realtime_dd} is supplied in UTC which is the easiest format to work with across timezones. This function
does not change \code{Date} from UTC. Rather \code{station_tz} specifies the local timezone name and is useful in instances where
\code{realtime_add_local_datetime} adjusts local_datetime to a common timezone that is not the \code{station_tz}. This function is most
useful when all stations exist within the same timezone.
}
\examples{
\dontrun{

realtime_dd(c("08MF005","02LA004")) \%>\%
 realtime_add_local_datetime()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy_db.R
\name{hy_test_db}
\alias{hy_test_db}
\alias{hy_downloaded_db}
\alias{hy_default_db}
\title{Get the location of the HYDAT database}
\usage{
hy_test_db()

hy_downloaded_db()

hy_default_db()
}
\value{
The file location of a HYDAT database.
}
\description{
The full HYDAT database needs to be downloaded from \link{download_hydat}, but for testing
purposes, a small test database is included in this package. Use
\code{hydat_path = hy_test_db()} in hy_* functions to explicitly use the test database;
use \code{hydat_path = hy_downloaded_db()} to explicitly use the full, most recent
downloaded database (this is also the path returned by \code{hy_default_db()}).
}
\examples{
\dontrun{
hy_test_db()
hy_downloaded_db()
hy_default_db()
}

}
\seealso{
\link{hy_src}, \link{hy_set_default_db}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/realtime.R
\name{realtime_daily_mean}
\alias{realtime_daily_mean}
\title{Calculate daily means from higher resolution realtime data}
\usage{
realtime_daily_mean(.data, na.rm = FALSE)
}
\arguments{
\item{.data}{A data argument that is designed to take only the output of realtime_dd}

\item{na.rm}{a logical value indicating whether NA values should be stripped before the computation proceeds.}
}
\description{
This function is meant to be used within a pipe as a means of easily moving from higher resolution
data to daily means.
}
\examples{
\dontrun{
realtime_dd("08MF005") \%>\% realtime_daily_mean()
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_version}
\alias{hy_version}
\title{Extract version number from HYDAT database}
\source{
HYDAT
}
\usage{
hy_version(hydat_path = NULL)
}
\arguments{
\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}
}
\value{
version number and release date
}
\description{
A function to get version number of hydat
}
\examples{
\dontrun{
hy_version()
}


}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_datum_list}()},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()}
}
\concept{HYDAT functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/hy.R
\name{hy_datum_list}
\alias{hy_datum_list}
\title{Extract datum list from HYDAT database}
\source{
HYDAT
}
\usage{
hy_datum_list(hydat_path = NULL)
}
\arguments{
\item{hydat_path}{The path to the hydat database or NULL to use the default location
used by \link{download_hydat}. It is also possible to pass in an existing
\link[dplyr]{src_sqlite} such that the database only needs to be opened once per
user-level call.}
}
\value{
A tibble of DATUMS
}
\description{
DATUM look-up Table
}
\examples{
\dontrun{
hy_datum_list()
}

}
\seealso{
Other HYDAT functions: 
\code{\link{hy_agency_list}()},
\code{\link{hy_annual_instant_peaks}()},
\code{\link{hy_annual_stats}()},
\code{\link{hy_daily_flows}()},
\code{\link{hy_daily_levels}()},
\code{\link{hy_daily}()},
\code{\link{hy_data_symbols}},
\code{\link{hy_data_types}},
\code{\link{hy_monthly_flows}()},
\code{\link{hy_monthly_levels}()},
\code{\link{hy_reg_office_list}()},
\code{\link{hy_sed_daily_loads}()},
\code{\link{hy_sed_daily_suscon}()},
\code{\link{hy_sed_monthly_loads}()},
\code{\link{hy_sed_monthly_suscon}()},
\code{\link{hy_sed_samples_psd}()},
\code{\link{hy_sed_samples}()},
\code{\link{hy_stations}()},
\code{\link{hy_stn_data_coll}()},
\code{\link{hy_stn_data_range}()},
\code{\link{hy_stn_op_schedule}()},
\code{\link{hy_stn_regulation}()},
\code{\link{hy_version}()}
}
\concept{HYDAT functions}
