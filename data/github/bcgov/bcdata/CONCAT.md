
<!--
Copyright 2018 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.
-->

# bcdata <a href='https://bcgov.github.io/bcdata/'><img src='man/figures/logo.png' align="right" height="139" /></a>

<!-- badges: start -->

[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![R build
status](https://github.com/bcgov/bcdata/workflows/R-CMD-check/badge.svg)](https://github.com/bcgov/bcdata)
[![Codecov test
coverage](https://codecov.io/gh/bcgov/bcdata/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcgov/bcdata?branch=master)
[![CRAN_Status_Badge](https://www.r-pkg.org/badges/version/bcdata)](https://cran.r-project.org/package=bcdata)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/bcdata?color=brightgreen)](https://CRAN.R-project.org/package=bcdata)
[![cran
checks](https://cranchecks.info/badges/worst/bcdata)](https://CRAN.R-project.org/web/checks/check_results_bcdata.html)
[![img](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)](https://github.com/bcgov/repomountie/blob/master/doc/lifecycle-badges.md)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02927/status.svg)](https://doi.org/10.21105/joss.02927)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737824.svg)](https://doi.org/10.5281/zenodo.4737824)
<!-- badges: end -->

An R package 📦 for searching & retrieving data from the [B.C. Data
Catalogue](https://catalogue.data.gov.bc.ca).

-   `bcdc_browse()` - Open the catalogue in your default browser
-   `bcdc_search()` - Search records in the catalogue
-   `bcdc_search_facets()` - List catalogue facet search options
-   `bcdc_get_record()` - Print a catalogue record
-   `bcdc_tidy_resources()` - Get a data frame of resources for a record
-   `bcdc_get_data()` - Get catalogue data
-   `bcdc_query_geodata()` - Get & query catalogue geospatial data
    available through a [Web Feature
    Service](https://en.wikipedia.org/wiki/Web_Feature_Service)

**Note:** The `bcdata` package supports downloading *most* file types,
including zip archives. It will do its best to identify and read data
from zip files, however if there are multiple data files in the zip, or
data files that `bcdata` doesn’t know how to import, it will fail. If
you encounter a file type in the B.C. Data Catalogue not currently
supported by `bcdata` please file an
[issue](https://github.com/bcgov/bcdata/issues/).

### Reference

[bcdata package 📦 home page and reference
guide](https://bcgov.github.io/bcdata/)

### Installation

You can install `bcdata` directly from
[CRAN](https://cran.r-project.org/package=bcdata):

``` r
install.packages("bcdata")
```

To install the development version from GitHub, use the
[remotes](https://cran.r-project.org/package=remotes) package:

``` r
install.packages("remotes")

remotes::install_github("bcgov/bcdata")
library(bcdata)
```

### Vignettes

-   [Get Started with
    bcdata](https://bcgov.github.io/bcdata/articles/bcdata.html)
-   [Querying Spatial Data with
    bcdata](https://bcgov.github.io/bcdata/articles/efficiently-query-spatial-data-in-the-bc-data-catalogue.html)
-   [Exploring Silviculture Data with
    bcdata](https://bcgov.github.io/bcdata/articles/explore-silviculture-data-using-bcdata.html)

### Methods for `bcdc_promise`

The `bcdc_query_geodata()` returns an object of the class
`bcdc_promise`. We have written an ever growing list methods for this
class. You can use these methods directly on a object returned by
`bcdc_query_geodata()`. Here are all the methods for the `bcdc_promise`
class:

-   `as_tibble`
-   `collect`
-   `filter`
-   `head`
-   `mutate`
-   `names`
-   `print`
-   `select`
-   `show_query`
-   `tail`

### BCDC Authentication

If you are an authorized editor of the B.C. Data Catalogue you may want
to access records that are not publicly available (e.g., in DRAFT,
waiting to be published). This can be done by authenticating with the
catalogue with an API key.

***Important Note:*** *Your API key is like a password and you must take
care to keep it private. Do not share it, and be careful to not include
it in any scripts or accidentally commit it to GitHub.*

You can log in to the catalogue to obtain your API key, then store it as
an environment variable in your [`.Renviron`
file](https://rstats.wtf/r-startup.html#renviron). The environment
variable must be called `BCDC_KEY`, set like this:

    BCDC_KEY=your-api-key

This way, the relevant bcdata functions will read that key and use it to
authorize your calls to the catalogue, allowing you to access additional
records that you are authorized to see if you were logged into the
catalogue web interface. Functions that benefit from this are:

-   `bcdc_search()`
-   `bcdc_list()`
-   `bcdc_get_record()`
-   `bcdc_get_data()`

### Getting Help or Reporting an Issue

To report bugs/issues/feature requests, please file an
[issue](https://github.com/bcgov/bcdata/issues/).

### How to Contribute

If you would like to contribute to the package, please see our
[CONTRIBUTING](https://github.com/bcgov/bcdata/blob/master/CONTRIBUTING.md)
guidelines.

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/bcgov/bcdata/blob/master/CODE_OF_CONDUCT.md).
By participating in this project you agree to abide by its terms.

### Citation


    To cite bcdata in publications please use:

      Teucher et al., (2021). bcdata: An R package for searching and
      retrieving data from the B.C. Data Catalogue. Journal of Open Source
      Software, 6(61), 2927, https://doi.org/10.21105/joss.02927

    A BibTeX entry for LaTeX users is

      @Article{,
        doi = {10.21105/joss.02927},
        year = {2021},
        publisher = {The Open Journal},
        volume = {6},
        number = {61},
        pages = {2927},
        author = {Andrew C. Teucher and Sam J. Albers and Stephanie L. Hazlitt},
        title = {bcdata: An R package for searching and retrieving data from the B.C. Data Catalogue},
        journal = {Journal of Open Source Software},
      }

### License

Copyright 2018 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the “License”); you may
not use this file except in compliance with the License. You may obtain
a copy of the License at

<https://www.apache.org/licenses/LICENSE-2.0>

Unless required by applicable law or agreed to in writing, software
distributed under the License is distributed on an “AS IS” BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and
limitations under the License.

------------------------------------------------------------------------

*This project was created using the
[bcgovr](https://github.com/bcgov/bcgovr) package.*
# bcdata (development version)

* Added `bcdc_get_citation` to generate bibliographic entries (via `utils::bibentry`) for individuals records. #273
* Results from `bcdc_search()` (objects of class `"bcdc_recordlist"`) now print 50 records by default, instead of 10. In addition, there is a new `[` method for `"bcdc_recordlist"` objects, allowing you to subset these lists and still have a nice printout (#288).
* Ensure that `bcdc_get_data()` fails informatively when a given resource doesn't exist in a record (#290)


# bcdata 0.3.0

* The [BC Data Catalogue and its API have been updated](https://www2.gov.bc.ca/gov/content?id=8A553CABCCDD434D8614D1CA92B03400), requiring changes to the `bcdata` package, most of which are internal only (#283). These should be mostly invisible to the user, except for the removal of the `type` search facet in `bcdc_search()` and `bcdc_search_facets()`. If you use an API key (authorized catalogue editors only), you will need to login to the new catalogue and get your updated key and set the value your `BCDC_KEY` environment variable to the new key.

### IMPROVEMENTS
* Add `names` method for `bcdc.promise` objects. You can now call `names` on an object produced by `bcdc_query_geodata`. This is handy when trying to figure out exact column spelling etc. #278

### BUG FIXES
* Fix bug where sticky column were incorrectly identified in `bcdc_describe_feature` (#279)

# bcdata 0.2.4

* Code in `.onLoad()` that sent a request to the wfs getCapabilities endpoint could cause the package to fail to load. This was moved into an internal function `bcdc_get_capabilities()` that makes the request the first time it's required, and stores the result for the remainder of the session (#271)
* testthat is now used conditionally to only run tests if the testthat package is installed.

# bcdata 0.2.3

### IMPROVEMENTS
- Setting the `bcdata.single_download_limit` limit dynamically from the getCapabilities endpoint. #256
- `bcdc_describe_feature` now joins an object description column to the returned object to provide more information about a field directly in R. #241, #259
- Better documentation and information surrounding the `bcdata.max_geom_pred_size` option. #243, #258 
- Add new function `bcdc_check_geom_size` so users can check prior to submitting a WFS request with `filter` to see 
if the request will require a bounding box conversion. #243, #258
- Better documentation and messaging about when and why paginated requests are required by `bcdc_query_geodata()`. #240, #264
- Better documentation and print method for what records are suitable for use with `bcdc_query_geodata()`. #265, #267

# bcdata 0.2.2
### IMPROVEMENTS
* Added `bcdc_list_groups` and `bcdc_list_group_records` to provide the ability to query on the group endpoint of the catalogue API. #234
* Added new option `bcdata.single_download_limit` to enable setting the maximum number of records an object can be before forcing a paginated download (#252)

### BUG FIXES
* Fixed bug in `collect.bcdc_promise` where the wrong parameter name in `crul::Paginator$new()` resulted in an error in paginated wfs requests (#250, thanks @meztez)
* Fixed a bug where the name of `bcdata.chunk_limit` option had a typo, so that it was not actually used properly (#252)

# bcdata 0.2.1

### BUG FIXES
* Remove link for pipe documentation for simplicity.
* Fixed bug where using many `as.` functions (e.g., `as.Date()`, `as.character()`, `as.numeric()`) in a filter statement would fail. (#218, #219)

### MAINTENANCE
* Updated internal SQL translation to use `DBI` S4 generics (`DBI::dbQuoteIdentifier()` is now used instead of 
  `dbplyr::sql_escape_ident()` and `DBI::dbQuoteString()` instead of `dbplyr::sql_escape_string()`), to comply 
  with upcoming `dbplyr` 2.0 release (#225, #225; https://github.com/tidyverse/dbplyr/issues/385)
* Wrapped all examples that call web resources in `try()` to avoid spurious check failures (#229).

# bcdata 0.2.0

### BREAKING CHANGES
* Rename `selectable` column from `bcdc_describe_feature` to `sticky` and modify corresponding docs and tests (#180).

### IMPROVEMENTS
* Add explore-silviculture-data-using-bcdata vignette/article. h/t @hgriesbauer 
* Add `head` and `tail` methods for `bcdc.promise` objects. Thanks to @hgriesbauer for the suggestion! (#182, #186)
* Provide `as_tibble` as an alias for `collect` in line with `dbplyr` behaviour (#166)
* Geometry predicates can now take a `bbox` object as well as an `sf*` object (#176)
* When reading in excel files, `bcdc_get_data` now outputs a messages indicating the presence and names of any sheets (#190)
* `bcdc_get_data()` & `bcdc_query_geodata()` will now work with full B.C. data catalogue url including resource (#125, #196)
* `bcdc_sf` objects now have an `time_downloaded` attribute
* Authorized B.C. Data Catalogue editors can now authenticate with the catalogue by setting 
a `BCDC_KEY` environment variable with their catalogue API token (https://github.com/bcgov/bcdata#bcdc-authentication; #208).

### BUG FIXES
* Fix `select`, `filter` and `mutate` roxygen so that bcdata specific documentation to these methods is available
* Add tests for attributes

# bcdata 0.1.2

### IMPROVEMENTS
* Add `bcdc_tidy_resources` for retrieving a data frame containing the metadata for all resources from a single B.C. Data Catalogue record (PR#149, #147)
* Add a more decorative record print method  (#73)
* More reliable detection of layer name for a wfs call in `bcdc_query_geodata()` (#129, #138, #139)
* Add `mutate` method for bcdc_promise that only fails and suggest an alternative approach. (PR#134)
* Add back in querying vignette
* Using `tidyselect` so that `select.bcdc_promise` behaviour is typical of `dplyr::select` ($140, #138)
* Using GitHub actions for CI. 

### MINOR BREAKING CHANGES
* Remove `BEYOND()` and `RELATE()` geometry predicates as they are currently not fully supported by geoserver
                                                        
### BUG FIXES
* Now precompiling vignettes so that queries are submitted locally and no actually requests are made from CRAN (#151)
* Fix `NOTE: Namespace in Imports field not imported from: ‘methods’` error on CRAN (#145)
* Fixed a bug where functions nested inside geometry predicates were not evaluated (#146, #154)
* Fixed a bug where `DWITHIN` wasn't working because `units` needed to be unquoted (#154)
* Fixed a bug where `BBOX()` used in a `filter()` statement combined with `bcdc_query_geodata()` did not work (#135, #137, #131)
* Fixed a bug where layer names with a number in them would not work in `bcdc_query_geodata()` (#126, #127)


# bcdata 0.1.1

* Expand and standardize checking w[ms]f features to make package more resistant to slight warehouse API changes. 
* Data retrieval functions now work with BCGW name (#106)
* Add CITATION file (#104)
* Increased test coverage (#112)
* Skipping all tests on CRAN that require a web connection
* Better and more informative error message when experiencing http failures occur (#121)
* Added print methods for `show_query`
* Change examples to donttest
* Added verbose argument to `bcdc_get_record` to enable suppressing console writing

# bcdata 0.1.0

* Initial Release
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
reported by contacting a project maintainer at andy.tecuher@gov.bc.ca, sam.albers@gov.bc.ca, or stephanie.hazlitt@gov.bc.ca. All complaints will be reviewed and investigated 
and will result in a response that is deemed necessary and appropriate to the 
circumstances. Maintainers are obligated to maintain confidentiality with regard 
to the reporter of an incident.


This Code of Conduct is adapted from the [Contributor Covenant][homepage],
version 1.3.0, available at
[http://contributor-covenant.org/version/1/3/0/][version]

[homepage]: http://contributor-covenant.org
[version]: http://contributor-covenant.org/version/1/3/0/

---
*This project was created using the [bcgovr](https://github.com/bcgov/bcgovr) package.*
## How to contribute
Government employees, public and members of the private sector are encouraged to contribute to the repository by **forking and submitting a pull request**. 

(If you are new to GitHub, you might start with a [basic tutorial](https://help.github.com/articles/set-up-git) and  check out a more detailed guide to [pull requests](https://help.github.com/articles/using-pull-requests/).)

Pull requests will be evaluated by the repository guardians on a schedule and if deemed beneficial will be committed to the master.

All contributors retain the original copyright to their stuff, but by contributing to this project, you grant a world-wide, royalty-free, perpetual, irrevocable, non-exclusive, transferable license to all users **under the terms of the license under which this project is distributed.**

---
*This project was created using the [bcgovr](https://github.com/bcgov/bcgovr) package.*
### CRAN check issues

There were no CRAN check issues to be fixed.

## Test environments

* local macOS install (macOS Big Sur 11.6), R 4.1.1
* local Windows 10 install, R 4.1.1
* ubuntu 18.04.6 (on github actions), R 4.1.1, R 4.0.5, R 3.5
* ubuntu 18.04.6 (on github actions), R-devel (2021-10-25 r81105)
* Windows Server 2019 (on github actions), R 4.1.1
* macOS Catalina 10.15.7 (on github actions), R 4.1.1
* win-builder (R-devel: 2021-10-25 r81104)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

We checked 1 reverse dependency, comparing R CMD check results across CRAN and dev versions of this package.

We saw 0 new problems
We failed to check 0 packages
---
name: Bug report or feature request
about: Describe a bug you've seen or make a case for a new feature
---

Please briefly describe your problem and what output you expect. Please include a link to the B.C. data catalogue record that you are trying to access if applicable. 

Please include a minimal reproducible example (AKA a reprex). If you've never heard of a [reprex](http://reprex.tidyverse.org/) before, start by reading <https://www.tidyverse.org/help/#reprex>.

Brief description of the problem

```r
# insert reprex here
```
## revdepcheck results

We checked 1 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                        |
|:--------|:----------------------------|
|version  |R version 4.1.1 (2021-08-10) |
|os       |macOS Big Sur 11.6           |
|system   |x86_64, darwin17.0           |
|ui       |X11                          |
|language |(EN)                         |
|collate  |en_CA.UTF-8                  |
|ctype    |en_CA.UTF-8                  |
|tz       |America/Vancouver            |
|date     |2021-10-27                   |

# Dependencies

|package           |old     |new     |Δ  |
|:-----------------|:-------|:-------|:--|
|bcdata            |0.2.4   |0.3.0   |*  |
|assertthat        |0.2.1   |0.2.1   |   |
|base64enc         |0.1-3   |0.1-3   |   |
|bit               |4.0.4   |4.0.4   |   |
|bit64             |4.0.5   |4.0.5   |   |
|blob              |1.2.2   |1.2.2   |   |
|cellranger        |1.1.0   |1.1.0   |   |
|classInt          |0.4-3   |0.4-3   |   |
|cli               |3.1.0   |3.1.0   |   |
|clipr             |0.7.1   |0.7.1   |   |
|colorspace        |2.0-2   |2.0-2   |   |
|cpp11             |0.4.0   |0.4.0   |   |
|crayon            |1.4.1   |1.4.1   |   |
|crosstalk         |1.1.1   |1.1.1   |   |
|crul              |1.1.0   |1.1.0   |   |
|curl              |4.3.2   |4.3.2   |   |
|DBI               |1.1.1   |1.1.1   |   |
|dbplyr            |2.1.1   |2.1.1   |   |
|digest            |0.6.28  |0.6.28  |   |
|dplyr             |1.0.7   |1.0.7   |   |
|e1071             |1.7-9   |1.7-9   |   |
|ellipsis          |0.3.2   |0.3.2   |   |
|fansi             |0.5.0   |0.5.0   |   |
|farver            |2.1.0   |2.1.0   |   |
|fastmap           |1.1.0   |1.1.0   |   |
|generics          |0.1.1   |0.1.1   |   |
|ggplot2           |3.3.5   |3.3.5   |   |
|glue              |1.4.2   |1.4.2   |   |
|gridExtra         |2.3     |2.3     |   |
|gtable            |0.3.0   |0.3.0   |   |
|hms               |1.1.1   |1.1.1   |   |
|htmltools         |0.5.2   |0.5.2   |   |
|htmlwidgets       |1.5.4   |1.5.4   |   |
|httpcode          |0.3.0   |0.3.0   |   |
|isoband           |0.2.5   |0.2.5   |   |
|jsonlite          |1.7.2   |1.7.2   |   |
|labeling          |0.4.2   |0.4.2   |   |
|lazyeval          |0.2.2   |0.2.2   |   |
|leaflet           |2.0.4.1 |2.0.4.1 |   |
|leaflet.extras    |1.0.0   |1.0.0   |   |
|leaflet.providers |1.9.0   |1.9.0   |   |
|lifecycle         |1.0.1   |1.0.1   |   |
|magrittr          |2.0.1   |2.0.1   |   |
|markdown          |1.1     |1.1     |   |
|mime              |0.12    |0.12    |   |
|munsell           |0.5.0   |0.5.0   |   |
|pillar            |1.6.4   |1.6.4   |   |
|pkgconfig         |2.0.3   |2.0.3   |   |
|png               |0.1-7   |0.1-7   |   |
|prettyunits       |1.1.1   |1.1.1   |   |
|progress          |1.2.2   |1.2.2   |   |
|proxy             |0.4-26  |0.4-26  |   |
|purrr             |0.3.4   |0.3.4   |   |
|R6                |2.5.1   |2.5.1   |   |
|raster            |3.5-2   |3.5-2   |   |
|RColorBrewer      |1.1-2   |1.1-2   |   |
|Rcpp              |1.0.7   |1.0.7   |   |
|readr             |2.0.2   |2.0.2   |   |
|readxl            |1.3.1   |1.3.1   |   |
|rematch           |1.0.1   |1.0.1   |   |
|rlang             |0.4.12  |0.4.12  |   |
|s2                |1.0.7   |1.0.7   |   |
|scales            |1.1.1   |1.1.1   |   |
|sf                |1.0-3   |1.0-3   |   |
|sp                |1.4-5   |1.4-5   |   |
|stringi           |1.7.5   |1.7.5   |   |
|stringr           |1.4.0   |1.4.0   |   |
|terra             |1.4-11  |1.4-11  |   |
|tibble            |3.1.5   |3.1.5   |   |
|tidyselect        |1.1.1   |1.1.1   |   |
|triebeard         |0.3.0   |0.3.0   |   |
|tzdb              |0.2.0   |0.2.0   |   |
|units             |0.7-2   |0.7-2   |   |
|urltools          |1.7.3   |1.7.3   |   |
|utf8              |1.2.2   |1.2.2   |   |
|vctrs             |0.3.8   |0.3.8   |   |
|viridis           |0.6.2   |0.6.2   |   |
|viridisLite       |0.4.0   |0.4.0   |   |
|vroom             |1.5.5   |1.5.5   |   |
|withr             |2.4.2   |2.4.2   |   |
|wk                |0.5.0   |0.5.0   |   |
|xfun              |0.27    |0.27    |   |
|xml2              |1.3.2   |1.3.2   |   |
|yaml              |2.2.1   |2.2.1   |   |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*