
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

*Wow, no problems at all. :)**Wow, no problems at all. :)*---
output:
  github_document:
html_preview: true
---

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



```{r setup, echo = FALSE, warning = FALSE, message = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/"
)
library(bcdata)
```


# bcdata <a href='https://bcgov.github.io/bcdata/'><img src='man/figures/logo.png' align="right" height="139" /></a>



<!-- badges: start -->
[![License](https://img.shields.io/badge/License-Apache%202.0-blue.svg)](https://opensource.org/licenses/Apache-2.0)
[![R build status](https://github.com/bcgov/bcdata/workflows/R-CMD-check/badge.svg)](https://github.com/bcgov/bcdata)
[![Codecov test coverage](https://codecov.io/gh/bcgov/bcdata/branch/master/graph/badge.svg)](https://app.codecov.io/gh/bcgov/bcdata?branch=master)
[![CRAN\_Status\_Badge](https://www.r-pkg.org/badges/version/bcdata)](https://cran.r-project.org/package=bcdata) [![CRAN Downloads](https://cranlogs.r-pkg.org/badges/bcdata?color=brightgreen)](https://CRAN.R-project.org/package=bcdata) [![cran checks](https://cranchecks.info/badges/worst/bcdata)](https://CRAN.R-project.org/web/checks/check_results_bcdata.html)
[![img](https://img.shields.io/badge/Lifecycle-Maturing-007EC6)](https://github.com/bcgov/repomountie/blob/master/doc/lifecycle-badges.md)
[![DOI](https://joss.theoj.org/papers/10.21105/joss.02927/status.svg)](https://doi.org/10.21105/joss.02927)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.4737824.svg)](https://doi.org/10.5281/zenodo.4737824)
<!-- badges: end -->

An R package 📦 for searching & retrieving data from the [B.C. Data Catalogue]( https://catalogue.data.gov.bc.ca).

- `bcdc_browse()` - Open the catalogue in your default browser
- `bcdc_search()` - Search records in the catalogue
- `bcdc_search_facets()` - List catalogue facet search options
- `bcdc_get_record()` - Print a catalogue record
- `bcdc_tidy_resources()` - Get a data frame of resources for a record
- `bcdc_get_data()` - Get catalogue data
- `bcdc_query_geodata()` - Get & query catalogue geospatial data available through a [Web Feature Service](https://en.wikipedia.org/wiki/Web_Feature_Service)

**Note:** The `bcdata` package supports downloading _most_ file types, including zip archives. It will do its best to identify and read data from
zip files, however if there are multiple data files in the zip, or data files that `bcdata` doesn't know how to import, it will fail. 
If you encounter a file type in the B.C. Data Catalogue not currently supported by `bcdata` please file an [issue](https://github.com/bcgov/bcdata/issues/). 

### Reference
[bcdata package 📦 home page and reference guide](https://bcgov.github.io/bcdata/)

### Installation
You can install `bcdata` directly from [CRAN](https://cran.r-project.org/package=bcdata): 

```{r eval=FALSE}
install.packages("bcdata")
```

To install the development version from GitHub, use the [remotes](https://cran.r-project.org/package=remotes) package:

```{r eval=FALSE}
install.packages("remotes")

remotes::install_github("bcgov/bcdata")
library(bcdata)
```


### Vignettes

- [Get Started with bcdata](https://bcgov.github.io/bcdata/articles/bcdata.html)
- [Querying Spatial Data with bcdata](https://bcgov.github.io/bcdata/articles/efficiently-query-spatial-data-in-the-bc-data-catalogue.html)
- [Exploring Silviculture Data with bcdata](https://bcgov.github.io/bcdata/articles/explore-silviculture-data-using-bcdata.html)

### Methods for `bcdc_promise`

The `bcdc_query_geodata()` returns an object of the class `bcdc_promise`. We have written an ever growing list methods for this class. You can use these methods directly on a object returned by `bcdc_query_geodata()`. Here are all the methods for the `bcdc_promise` class:

```{r echo=FALSE, results='asis'}
bcdc_methods <- methods(class = "bcdc_promise")
bcdc_methods <- sort(attributes(bcdc_methods)$info[,c("generic"), ])

cat(paste0("- `", bcdc_methods, "`", collapse = "\n"))
```


### BCDC Authentication

If you are an authorized editor of the B.C. Data Catalogue you may want to
access records that are not publicly available (e.g., in DRAFT, waiting to be
published). This can be done by authenticating with the catalogue with an API
key.

_**Important Note:**_ *Your API key is like a password and you must take care to
keep it private. Do not share it, and be careful to not include it in any
scripts or accidentally commit it to GitHub.*

You can log in to the catalogue to obtain your API key, then store it as an
environment variable in your [`.Renviron` file](https://rstats.wtf/r-startup.html#renviron). 
The environment variable must be called `BCDC_KEY`, set like this:

```
BCDC_KEY=your-api-key
```

This way, the relevant bcdata functions will read that key and use it to
authorize your calls to the catalogue, allowing you to access additional records
that you are authorized to see if you were logged into the catalogue web
interface. Functions that benefit from this are:

- `bcdc_search()`
- `bcdc_list()`
- `bcdc_get_record()`
- `bcdc_get_data()`

### Getting Help or Reporting an Issue

To report bugs/issues/feature requests, please file an [issue](https://github.com/bcgov/bcdata/issues/).

### How to Contribute

If you would like to contribute to the package, please see our 
[CONTRIBUTING](https://github.com/bcgov/bcdata/blob/master/CONTRIBUTING.md) guidelines.

Please note that this project is released with a [Contributor Code of Conduct](https://github.com/bcgov/bcdata/blob/master/CODE_OF_CONDUCT.md). By participating in this project you agree to abide by its terms.

### Citation
```{r, echo=FALSE, comment=""}
citation("bcdata")
```


### License

Copyright 2018 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

https://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.

---
*This project was created using the [bcgovr](https://github.com/bcgov/bcgovr) package.* 
---
title: "Querying Spatial Data with bcdata"
date: "2021-10-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Querying Spatial Data with bcdata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<!--
Copyright 2019 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.
-->




This vignette illustrates how to use `bcdata::bcdc_query_geodata` to request and query geospatial data that has an associated [Web Feature Service](https://en.wikipedia.org/wiki/Web_Feature_Service) from the [B.C. Data Catalogue](https://catalogue.data.gov.bc.ca/dataset). To illustrate, we will request and merge two spatial data sets from the catalogue---school district and greenspaces spatial data---and then examine the amount of park space contained within the boundaries of the Greater Victoria, Prince George and Kamloops/Thompson British Columbia school districts.

## Getting Started
First you need to load the package. We will also load the `sf` and `dplyr` packages to help us work with spatial data. You can learn more about the `sf` package [here](https://r-spatial.github.io/sf/) and `dplyr` [here](https://dplyr.tidyverse.org/):


```r
library(bcdata)
library(sf)
library(dplyr)
```


## Geospatial Data in the B.C. Data Catalogue
The [B.C. Data Catalogue](https://catalogue.data.gov.bc.ca/dataset) provides many data sets with spatial information through a [Web Feature Service (WFS)](https://en.wikipedia.org/wiki/Web_Feature_Service). Technically speaking, this means if we have an internet connection we can issue [HTTP](https://en.wikipedia.org/wiki/Hypertext_Transfer_Protocol) requests to the catalogue and seamlessly import the response data into R as an `sf` objects. The `bcdata` package provides a means to a) choose which layer you want and b) use `dplyr` verbs to specifically tailor your request.  A `dbplyr` backend is implemented so that requests are executed lazily meaning results are not transferred over the web until the user specifically calls the `collect` function. This approach mimics the `dplyr` verb translation to `SQL` seen for many database types. A good introduction to principles of `dbplyr` is available [here](https://cran.r-project.org/package=dbplyr/vignettes/dbplyr.html).

## School District Data
Our first step is to extract the [school district polygons](https://catalog.data.gov.bc.ca/dataset/78ec5279-4534-49a1-97e8-9d315936f08b) from the B.C. Data Catalogue. This layer is described using this command:


```r
bcdc_get_record("78ec5279-4534-49a1-97e8-9d315936f08b")
#> B.C. Data Catalogue Record: School Districts of BC
#> Name: school-districts-of-bc (ID: 78ec5279-4534-49a1-97e8-9d315936f08b)
#> Permalink: https://catalogue.data.gov.bc.ca/dataset/78ec5279-4534-49a1-97e8-9d315936f08b
#> Licence: Open Government Licence - British Columbia
#> Description: The School Districts dataset contains the spatial representation (polygon)
#>  of the current extent of the administrative areas defined under section 176(1) of the
#>  School Act for the purposes of preservation and promotion of the fundamental principle
#>  of local autonomy and control of public education at the public and governmental levels
#>  through locally elected school boards.
#> Available Resources (1):
#>  1. WMS getCapabilities request (wms)
#> Access the full 'Resources' data frame using:
#>  bcdc_tidy_resources('78ec5279-4534-49a1-97e8-9d315936f08b')
#> Query and filter this data using:
#>  bcdc_query_geodata('78ec5279-4534-49a1-97e8-9d315936f08b')
```

This data is the boundary of each school district. The key information in this metadata is that the layer has a resource in `"wms"` format ---which means it is available through a Web Feature Service. From this we know we can make use of `bcdc_query_geodata`.


```r
bcdc_query_geodata("78ec5279-4534-49a1-97e8-9d315936f08b")
#> Querying 'school-districts-of-bc' record
#> • Using collect() on this object will return 59 features and 9 fields
#> • At most six rows of the record are printed here
#> ───────────────────────────────────────────────────────────────────────────────────────────────────────
#> Simple feature collection with 6 features and 9 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 956376 ymin: 475108.4 xmax: 1635228 ymax: 901924.4
#> Projected CRS: NAD83 / BC Albers
#> # A tibble: 6 × 10
#>   id                 ADMIN_AREA_SID SCHOOL_DISTRICT… SCHOOL_DISTRICT_… FEATURE_CODE FEATURE_AREA_SQM
#>   <chr>                       <int> <chr>                        <int> <chr>                   <dbl>
#> 1 WHSE_TANTALIS.TA_…            300 Arrow Lakes                     10 FH10000300        7392472526.
#> 2 WHSE_TANTALIS.TA_…            301 Revelstoke                      19 FH10000300        9416076465.
#> 3 WHSE_TANTALIS.TA_…            302 Kootenay-Columb…                20 FH10000300        3072672101.
#> 4 WHSE_TANTALIS.TA_…            303 Vernon                          22 FH10000300        5588468673.
#> 5 WHSE_TANTALIS.TA_…            304 Central Okanagan                23 FH10000300        2916757936.
#> 6 WHSE_TANTALIS.TA_…            305 Cariboo-Chilcot…                27 FH10000300       61213520885.
#> # … with 4 more variables: FEATURE_LENGTH_M <dbl>, OBJECTID <int>, SE_ANNO_CAD_DATA <chr>,
#> #   geometry <POLYGON [m]>
```

This is the initial query to the data in the catalogue. What has been returned is *not* the actual data but rather a subset to help you tune your query. The printed output of this query offers several useful pieces of information. Because we have queried with a unique ID, we are shown the name of the record. We also receive instruction that using `collect()` will retrieve a given number of features and fields present for this query. Lastly, there is a reminder that what is printed is only the first 6 rows of the record. Since we are limiting the scope of analysis to the Greater Victoria, Prince George and Kamloops/Thompson school districts, we want to ask the catalogue for only those polygons just like we would in a typical `dplyr` workflow:


```r
bcdc_query_geodata("78ec5279-4534-49a1-97e8-9d315936f08b") %>%
  filter(SCHOOL_DISTRICT_NAME %in% c("Greater Victoria", "Prince George","Kamloops/Thompson"))
#> Querying 'school-districts-of-bc' record
#> • Using collect() on this object will return 1 features and 9 fields
#> • At most six rows of the record are printed here
#> ───────────────────────────────────────────────────────────────────────────────────────────────────────
#> Simple feature collection with 1 feature and 9 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 1126789 ymin: 821142.1 xmax: 1528155 ymax: 1224202
#> Projected CRS: NAD83 / BC Albers
#> # A tibble: 1 × 10
#>   id                 ADMIN_AREA_SID SCHOOL_DISTRICT… SCHOOL_DISTRICT_… FEATURE_CODE FEATURE_AREA_SQM
#>   <chr>                       <int> <chr>                        <int> <chr>                   <dbl>
#> 1 WHSE_TANTALIS.TA_…            328 Prince George                   57 FH10000300       51888780641.
#> # … with 4 more variables: FEATURE_LENGTH_M <dbl>, OBJECTID <int>, SE_ANNO_CAD_DATA <chr>,
#> #   geometry <POLYGON [m]>
```

To further tune our query, we can also request only the columns we want. Really we only want the school district column and the spatial information. During an actual analysis, it is possible that you may need to initially collect more data than you want to determine value to subset by. For example, there is currently no way to ask the catalogue for all possible unique values of `SCHOOL_DISTRICT_NAME`. Is that case the data will need to be brought into R and unique values will need to be determined there.


```r
bcdc_query_geodata("78ec5279-4534-49a1-97e8-9d315936f08b") %>%
  filter(SCHOOL_DISTRICT_NAME %in% c("Greater Victoria", "Prince George","Kamloops/Thompson")) %>%
  select(SCHOOL_DISTRICT_NAME)
#> Querying 'school-districts-of-bc' record
#> • Using collect() on this object will return 1 features and 5 fields
#> • At most six rows of the record are printed here
#> ───────────────────────────────────────────────────────────────────────────────────────────────────────
#> Simple feature collection with 1 feature and 5 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 1126789 ymin: 821142.1 xmax: 1528155 ymax: 1224202
#> Projected CRS: NAD83 / BC Albers
#> # A tibble: 1 × 6
#>   id        ADMIN_AREA_SID SCHOOL_DISTRICT_… SCHOOL_DISTRICT… OBJECTID                      geometry
#>   <chr>              <int> <chr>                        <int>    <int>                 <POLYGON [m]>
#> 1 WHSE_TAN…            328 Prince George                   57      562 ((1137478 1221549, 1137399 1…
```

Note that in the `select` statement, we did not explicitly ask for the spatial data and also that there are several columns present that we didn't select. This is because within each data set in the data catalogue, there are several columns that will always be returned regardless of what is selected. If you really don't want those columns after you `collect` the data, which we will take care of right now, you can drop them:


```r
districts <- bcdc_query_geodata("78ec5279-4534-49a1-97e8-9d315936f08b") %>%
  filter(SCHOOL_DISTRICT_NAME %in% c("Greater Victoria", "Prince George","Kamloops/Thompson")) %>%
  select(SCHOOL_DISTRICT_NAME) %>%
  collect()
```

Again note here that we have assigned the object a name and added the `collect` statement. This step happens when you have selected the data you want and wish to begin working with it in R like a normal `sf` object. For example, we can now plot these three school districts:


```r
plot(st_geometry(districts))
```

<img src="vignette-fig-districts-1.png" title="plot of chunk districts" alt="plot of chunk districts" width="100%" />

Now that we have the spatial boundaries narrowed by specific school districts we can perform some spatial operations to determine parks in the school districts.

## Greenspaces Data
For the purposes of this example, let's consider [this greenspace](https://catalogue.data.gov.bc.ca/dataset/6a2fea1b-0cc4-4fc2-8017-eaf755d516da) layer in the catalogue. This layer is described here:


```r
bcdc_get_record("6a2fea1b-0cc4-4fc2-8017-eaf755d516da")
#> B.C. Data Catalogue Record: Local and Regional Greenspaces
#> Name: local-and-regional-greenspaces (ID: 6a2fea1b-0cc4-4fc2-8017-eaf755d516da)
#> Permalink: https://catalogue.data.gov.bc.ca/dataset/6a2fea1b-0cc4-4fc2-8017-eaf755d516da
#> Licence: Open Government Licence - British Columbia
#> Description: This dataset contains spatial and attribute information for local and
#>  regional greenspaces in British Columbia. Local and regional greenspaces are municipal
#>  or regional district lands designated by local government agencies and managed for
#>  public enjoyment, ecosystem or wildlife values. Spatial boundaries were sourced from
#>  municipal and regional district web sites, which in some cases provide datasets under
#>  Open Government Licence, and in other cases, publicize parks and greenspaces on web maps
#>  or pdf maps. Boundaries were edge-matched to the ParcelMap BC cadastre.  This spatial
#>  layer contains multipart polygons.
#> Available Resources (2):
#>  1. WMS getCapabilities request (wms)
#>  2. LocalRegionalParksAttributeValues (xlsx)
#> Access the full 'Resources' data frame using:
#>  bcdc_tidy_resources('6a2fea1b-0cc4-4fc2-8017-eaf755d516da')
#> Query and filter this data using:
#>  bcdc_query_geodata('6a2fea1b-0cc4-4fc2-8017-eaf755d516da')
```

Again we recognize this is [WFS-enabled](https://en.wikipedia.org/wiki/Web_Feature_Service) geospatial data, which means we can make use of `bcdc_query_geodata`.


```r
bcdc_query_geodata("6a2fea1b-0cc4-4fc2-8017-eaf755d516da")
#> Querying 'local-and-regional-greenspaces' record
#> • Using collect() on this object will return 8708 features and 19 fields
#> • At most six rows of the record are printed here
#> ───────────────────────────────────────────────────────────────────────────────────────────────────────
#> Simple feature collection with 6 features and 19 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 1228153 ymin: 453419.9 xmax: 1240644 ymax: 467184.1
#> Projected CRS: NAD83 / BC Albers
#> # A tibble: 6 × 20
#>   id          LOCAL_REG_GREENS… PARK_NAME   PARK_TYPE PARK_PRIMARY_USE REGIONAL_DISTRI… MUNICIPALITY
#>   <chr>                   <int> <chr>       <chr>     <chr>            <chr>            <chr>       
#> 1 WHSE_BASEM…                59 Hazelnut M… Local     Park             Metro Vancouver  Surrey      
#> 2 WHSE_BASEM…                60 Hazelwood … Local     Park             Metro Vancouver  Surrey      
#> 3 WHSE_BASEM…                61 Hemlock Pa… Local     Park             Metro Vancouver  Surrey      
#> 4 WHSE_BASEM…                62 Heritage W… Local     Park             Metro Vancouver  Surrey      
#> 5 WHSE_BASEM…                63 Heron Park  Local     Park             Metro Vancouver  Surrey      
#> 6 WHSE_BASEM…                64 Hillcrest … Local     Park             Metro Vancouver  Surrey      
#> # … with 13 more variables: CIVIC_NUMBER <int>, CIVIC_NUMBER_SUFFIX <chr>, STREET_NAME <chr>,
#> #   LATITUDE <dbl>, LONGITUDE <dbl>, WHEN_UPDATED <date>, WEBSITE_URL <chr>,
#> #   LICENCE_COMMENTS <chr>, FEATURE_AREA_SQM <dbl>, FEATURE_LENGTH_M <dbl>, OBJECTID <int>,
#> #   SE_ANNO_CAD_DATA <chr>, geometry <POLYGON [m]>
```

Since we are interested in only "Park" data we can subset our query:


```r
bcdc_query_geodata("6a2fea1b-0cc4-4fc2-8017-eaf755d516da") %>%
  filter(PARK_PRIMARY_USE == "Park")
#> Querying 'local-and-regional-greenspaces' record
#> • Using collect() on this object will return 4373 features and 19 fields
#> • At most six rows of the record are printed here
#> ───────────────────────────────────────────────────────────────────────────────────────────────────────
#> Simple feature collection with 6 features and 19 fields
#> Geometry type: POLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 1228153 ymin: 453419.9 xmax: 1240644 ymax: 467184.1
#> Projected CRS: NAD83 / BC Albers
#> # A tibble: 6 × 20
#>   id          LOCAL_REG_GREENS… PARK_NAME   PARK_TYPE PARK_PRIMARY_USE REGIONAL_DISTRI… MUNICIPALITY
#>   <chr>                   <int> <chr>       <chr>     <chr>            <chr>            <chr>       
#> 1 WHSE_BASEM…                59 Hazelnut M… Local     Park             Metro Vancouver  Surrey      
#> 2 WHSE_BASEM…                60 Hazelwood … Local     Park             Metro Vancouver  Surrey      
#> 3 WHSE_BASEM…                61 Hemlock Pa… Local     Park             Metro Vancouver  Surrey      
#> 4 WHSE_BASEM…                62 Heritage W… Local     Park             Metro Vancouver  Surrey      
#> 5 WHSE_BASEM…                63 Heron Park  Local     Park             Metro Vancouver  Surrey      
#> 6 WHSE_BASEM…                64 Hillcrest … Local     Park             Metro Vancouver  Surrey      
#> # … with 13 more variables: CIVIC_NUMBER <int>, CIVIC_NUMBER_SUFFIX <chr>, STREET_NAME <chr>,
#> #   LATITUDE <dbl>, LONGITUDE <dbl>, WHEN_UPDATED <date>, WEBSITE_URL <chr>,
#> #   LICENCE_COMMENTS <chr>, FEATURE_AREA_SQM <dbl>, FEATURE_LENGTH_M <dbl>, OBJECTID <int>,
#> #   SE_ANNO_CAD_DATA <chr>, geometry <POLYGON [m]>
```

Here we see that this greatly reduces the number of features that we are dealing with (and correspondingly the amount of data that needs to be transferred over the web). Remember also that we still have not actually requested the full data set. This is just still a preview. Also this query still includes all municipal parks in BC while we only want the ones in the three school districts - the polygons defined by the `districts` object. To find that subset of parks we can make use of the built-in geometric operators which allow us to perform spatial operations remotely fine tuning our query even further. Here using the `INTERSECTS` function is appropriate and since this is a last tuning step, we can call `collect` and assign a name to this object. These requests can sometimes take quite a long:


```r
districts_parks <- bcdc_query_geodata("6a2fea1b-0cc4-4fc2-8017-eaf755d516da") %>%
  filter(PARK_PRIMARY_USE == "Park") %>%
  filter(INTERSECTS(districts)) %>%
  collect()
#> The object is too large to perform exact spatial operations using bcdata.
#> Object size: 948576 bytes
#> BC Data Threshold: 500000 bytes
#> Exceedance: 448576 bytes
#> See ?bcdc_check_geom_size for more details
#> A bounding box was drawn around the object passed to INTERSECTS and all features within the box will be returned.
```

Plotting both the filtered parks data and the district polygons reveals an important consideration when using `bcdata`:

<img src="vignette-fig-district_parks-1.png" title="plot of chunk district_parks" alt="plot of chunk district_parks" width="100%" />

In this example, many parks not contained within the school districts are included in the `districts_parks` object. This is because rather than a full intersection, `bcdata` draws a bounding box around all the polygons that are doing the intersection (in this case `district`) and does the intersection based on that bounding box. This behaviour is imposed by the Web Feature Server but controlled via the `bcdata.max_geom_pred_size` option (See `?bcdc_options` for default values). Using this example, you can check to see if the size of the `districts` object exceeded the current value of `bcdata.max_geom_pred_size`:


```r
bcdc_check_geom_size(districts)
#> The object is too large to perform exact spatial operations using bcdata.
#> Object size: 948576 bytes
#> BC Data Threshold: 500000 bytes
#> Exceedance: 448576 bytes
#> See ?bcdc_check_geom_size for more details
```

Drawing the bounding box illustrates this point:

<img src="vignette-fig-bbox-1.png" title="plot of chunk bbox" alt="plot of chunk bbox" width="100%" />

We are left with two options to get around this problem. First, we can simply do some additional processing with the `sf` package. Specifically we can use a spatial join to assign parks into their respective district:


```r
districts_parks_join <- districts_parks %>%
  st_join(districts, left = FALSE)
```

<img src="vignette-fig-dp_join-1.png" title="plot of chunk dp_join" alt="plot of chunk dp_join" width="100%" />

A second approach is to set an internal option (`bcdata.max_geom_pred_size`) and increase the threshold of when a bounding box is drawn. Options are set in R like this:

```r
options("bcdata.max_geom_pred_size" = {object size in bytes})
```

The value of `bcdata.max_geom_pred_size` is set conservatively so that requests to the Web Feature Service are more consistently successful. Increasing this value may result in invalid requests.

Finally, to address our original question of which school district has the most municipal park space we can calculate the area of each park polygon and then sum those areas by school district:


```r
districts_parks_join %>%
  mutate(area = st_area(geometry)) %>%
  st_set_geometry(NULL) %>%
  group_by(SCHOOL_DISTRICT_NAME) %>%
  summarise(total_area = sum(area)) %>%
  arrange(total_area)
#> # A tibble: 1 × 2
#>   SCHOOL_DISTRICT_NAME total_area
#>   <chr>                     [m^2]
#> 1 Prince George         12254346.
```

## Additional Useful Functions
There are a couple of other functions in `bcdata` that are useful to know when working with spatial data from the catalogue. `bcdc_describe_feature` gives the column names, whether the column is selectable, and the column types in both R and on the remote server:


```r
bcdc_describe_feature("6a2fea1b-0cc4-4fc2-8017-eaf755d516da")
#> # A tibble: 20 × 5
#>    col_name                sticky remote_col_type          local_col_type column_comments           
#>    <chr>                   <lgl>  <chr>                    <chr>          <chr>                     
#>  1 id                      TRUE   xsd:string               character      <NA>                      
#>  2 LOCAL_REG_GREENSPACE_ID TRUE   xsd:decimal              numeric        LOCAL_REG_GREENSPACE_ID i…
#>  3 PARK_NAME               FALSE  xsd:string               character      PARK NAME is the name of …
#>  4 PARK_TYPE               FALSE  xsd:string               character      PARK_TYPE is the type of …
#>  5 PARK_PRIMARY_USE        FALSE  xsd:string               character      PARK PRIMARY USE defines …
#>  6 REGIONAL_DISTRICT       FALSE  xsd:string               character      REGIONAL_DISTRICT is the …
#>  7 MUNICIPALITY            FALSE  xsd:string               character      MUNICIPALITY is the name …
#>  8 CIVIC_NUMBER            FALSE  xsd:decimal              numeric        CIVIC_NUMBER is the stree…
#>  9 CIVIC_NUMBER_SUFFIX     FALSE  xsd:string               character      CIVIC_NUMBER_SUFFIX is th…
#> 10 STREET_NAME             FALSE  xsd:string               character      STREET_NAME is the name o…
#> 11 LATITUDE                FALSE  xsd:decimal              numeric        LATITUDE is the geographi…
#> 12 LONGITUDE               FALSE  xsd:decimal              numeric        LONGITUDE is the geograph…
#> 13 WHEN_UPDATED            FALSE  xsd:date                 date           WHEN_UPDATED is the date …
#> 14 WEBSITE_URL             FALSE  xsd:string               character      WEBSITE_URL contains a li…
#> 15 LICENCE_COMMENTS        FALSE  xsd:string               character      LICENCE_COMMENTS describe…
#> 16 FEATURE_AREA_SQM        FALSE  xsd:decimal              numeric        FEATURE_AREA_SQM is the s…
#> 17 FEATURE_LENGTH_M        FALSE  xsd:decimal              numeric        FEATURE_LENGTH_M is the s…
#> 18 SHAPE                   FALSE  gml:GeometryPropertyType sfc geometry   SHAPE is the column used …
#> 19 OBJECTID                TRUE   xsd:decimal              numeric        OBJECTID is a column requ…
#> 20 SE_ANNO_CAD_DATA        FALSE  xsd:hexBinary            numeric        SE_ANNO_CAD_DATA is a bin…
```

This is a helpful initial step to learn column names and types when you construct your query.

Another useful function is `show_query()` which provides information on the request issued to the remote server:


```r
bcdc_query_geodata("6a2fea1b-0cc4-4fc2-8017-eaf755d516da") %>%
  filter(PARK_PRIMARY_USE == "Park") %>%
  filter(INTERSECTS(districts)) %>%
  show_query()
#> <url>
#> <body>
#> SERVICE: WFS VERSION: 2.0.0 REQUEST: GetFeature outputFormat: application/json typeNames:
#>  WHSE_BASEMAPPING.GBA_LOCAL_REG_GREENSPACES_SP SRSNAME: EPSG:3005 CQL_FILTER:
#>  (("PARK_PRIMARY_USE" = 'Park') AND (INTERSECTS(SHAPE, POLYGON ((1126789 821142.1,
#>  1528155 821142.1, 1528155 1224202, 1126789 1224202, 1126789 821142.1)))))
#> 
#> <full query url>
#> https://openmaps.gov.bc.ca/geo/pub/wfs?SERVICE=WFS&VERSION=2.0.0&REQUEST=GetFeature&outputFormat=application%2Fjson&typeNames=WHSE_BASEMAPPING.GBA_LOCAL_REG_GREENSPACES_SP&SRSNAME=EPSG%3A3005&CQL_FILTER=%28%28%22PARK_PRIMARY_USE%22%20%3D%20%27Park%27%29%20AND%20%28INTERSECTS%28SHAPE%2C%20POLYGON%20%28%281126789%20821142.1%2C%201528155%20821142.1%2C%201528155%201224202%2C%201126789%201224202%2C%201126789%20821142.1%29%29%29%29%29
```

This output is what being created by the dplyr code outlined above.

## Using B.C. Geographic Warehouse (BCGW) layer names

If you are familiar with the [B.C. Geographic Warehouse (BCGW)](https://www2.gov.bc.ca/gov/content/data/geographic-data-services/bc-spatial-data-infrastructure/bc-geographic-warehouse),
you may already know the name of a layer that you want from the BCGW.
`bcdc_query_geodata()` (as well as all other related functions)
supports supplying that name directly. For example, the
[record for the B.C. airports layer](https://catalogue.data.gov.bc.ca/dataset/bc-airports#object-description)
shows that the object name is `WHSE_IMAGERY_AND_BASE_MAPS.GSR_AIRPORTS_SVW`, and
we can use that in `bcdc_query_geodata()`:


```r
# Look at the columns available:
bcdc_describe_feature("WHSE_IMAGERY_AND_BASE_MAPS.GSR_AIRPORTS_SVW")
#> # A tibble: 42 × 5
#>    col_name                      sticky remote_col_type local_col_type column_comments              
#>    <chr>                         <lgl>  <chr>           <chr>          <chr>                        
#>  1 id                            TRUE   xsd:string      character      <NA>                         
#>  2 CUSTODIAN_ORG_DESCRIPTION     TRUE   xsd:string      character      CUSTODIAN_ORG_DESCRIPTION co…
#>  3 BUSINESS_CATEGORY_CLASS       TRUE   xsd:string      character      BUSINESS_CATEGORY_CLASS desi…
#>  4 BUSINESS_CATEGORY_DESCRIPTION TRUE   xsd:string      character      BUSINESS_CATEGORY_DESCRIPTIO…
#>  5 OCCUPANT_TYPE_DESCRIPTION     TRUE   xsd:string      character      OCCUPANT_TYPE_DESCRIPTION co…
#>  6 SOURCE_DATA_ID                TRUE   xsd:string      character      SOURCE_DATA_ID is a unique o…
#>  7 SUPPLIED_SOURCE_ID_IND        TRUE   xsd:string      character      SUPPLIED_SOURCE_ID_IND is an…
#>  8 AIRPORT_NAME                  TRUE   xsd:string      character      AIRPORT_NAME is a business n…
#>  9 DESCRIPTION                   FALSE  xsd:string      character      DESCRIPTION describes the Oc…
#> 10 PHYSICAL_ADDRESS              FALSE  xsd:string      character      PHYSICAL_ADDRESS contains th…
#> # … with 32 more rows

# Query the data with bcdc_query_geodata and filter + select:
bcdc_query_geodata("WHSE_IMAGERY_AND_BASE_MAPS.GSR_AIRPORTS_SVW") %>%
  filter(DESCRIPTION == "airport") %>%
  select(AIRPORT_NAME, LOCALITY, NUMBER_OF_RUNWAYS)
#> Querying 'WHSE_IMAGERY_AND_BASE_MAPS.GSR_AIRPORTS_SVW' record
#> • Using collect() on this object will return 37 features and 11 fields
#> • At most six rows of the record are printed here
#> ───────────────────────────────────────────────────────────────────────────────────────────────────────
#> Simple feature collection with 6 features and 11 fields
#> Geometry type: POINT
#> Dimension:     XY
#> Bounding box:  xmin: 833323.9 ymin: 406886.6 xmax: 1266385 ymax: 1054950
#> Projected CRS: NAD83 / BC Albers
#> # A tibble: 6 × 12
#>   id         CUSTODIAN_ORG_DESCRI… BUSINESS_CATEGO… BUSINESS_CATEGO… OCCUPANT_TYPE_D… SOURCE_DATA_ID
#>   <chr>      <chr>                 <chr>            <chr>            <chr>            <chr>         
#> 1 WHSE_IMAG… "Ministry of Forest,… airTransportati… Air Transportat… BC Airports      455           
#> 2 WHSE_IMAG… "Ministry of Forest,… airTransportati… Air Transportat… BC Airports      464           
#> 3 WHSE_IMAG… "Ministry of Forest,… airTransportati… Air Transportat… BC Airports      482           
#> 4 WHSE_IMAG… "Ministry of Forest,… airTransportati… Air Transportat… BC Airports      483           
#> 5 WHSE_IMAG… "Ministry of Forest,… airTransportati… Air Transportat… BC Airports      484           
#> 6 WHSE_IMAG… "Ministry of Forest,… airTransportati… Air Transportat… BC Airports      487           
#> # … with 6 more variables: SUPPLIED_SOURCE_ID_IND <chr>, AIRPORT_NAME <chr>, LOCALITY <chr>,
#> #   NUMBER_OF_RUNWAYS <int>, SEQUENCE_ID <int>, geometry <POINT [m]>
```
---
title: "Exploring Silviculture Data with bcdata"
date: "2021-10-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Exploring Silviculture Data with bcdata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<!--
Copyright 2019 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.
-->




This vignette will demonstrate how you can use `bcdata` to access and explore British Columbia silviculture data.


## Silviculture Data

Silviculture data for British Columbia are largely collected and stored through the [Reporting Silviculture Updates and Land Status Tracking System (RESULTS) database](https://www2.gov.bc.ca/gov/content?id=F62992DA5B324446AD5E1B5BFFA404CC)---including information on the management of openings, disturbances, silviculture activities and obligation declarations as required by the [Forest and Range Practices Act](https://www2.gov.bc.ca/gov/content?id=D65006C0326B4E0F838194044C10CA86). [RESULTS data sets are publically available through the B.C. Data Catalogue](https://catalogue.data.gov.bc.ca/dataset?q=RESULTS+silviculture) under the [Open Government Licence - British Columbia](https://www2.gov.bc.ca/gov/content?id=A519A56BC2BF44E4A008B33FCF527F61).


## Reforestation with Climate Considerations

[Western larch](https://en.wikipedia.org/wiki/Larix_occidentalis) has been identified as a tree species that will be well-adapted to projected future climates in northern British Columbia.  In 2013, forest policy was developed that allowed foresters to plant western larch outside of its natural range as a climate change adaptation silviculture measure.

#### How much western larch has been planted in the Prince George Natural Resource District?

Let's use the `bcdata` package to search, retrieve and explore the RESULTS silviculture data and answer this question.


## Getting Started

To start, let's load the `bcdata` package. We will also load the `dplyr` and `ggplot2` packages to help us work with the data. You can learn more about the `dplyr` package [here](https://dplyr.tidyverse.org/) and the `ggplot2` package [here](https://ggplot2.tidyverse.org/):


```r
library(bcdata)
library(dplyr)
library(ggplot2)
```


## Getting the Data

We can gather the data we need to answer our question using the [RESULTS - silviculture forest cover dataset](https://catalogue.data.gov.bc.ca/dataset/258bb088-4113-47b1-b568-ce20bd64e3e3). First, let's take a look at the metadata record using `bcdc_get_record()`:


```r
# Get the metadata using the human-readable record name
bcdc_get_record("results-forest-cover-silviculture")
#> B.C. Data Catalogue Record: RESULTS - Forest Cover Silviculture
#> Name: results-forest-cover-silviculture (ID: 258bb088-4113-47b1-b568-ce20bd64e3e3)
#> Permalink: https://catalogue.data.gov.bc.ca/dataset/258bb088-4113-47b1-b568-ce20bd64e3e3
#> Licence: Open Government Licence - British Columbia
#> Description: RESULTS opening's forest cover poylgons with silviculture component
#>  provided.  Current forest cover subimssion into RESULTS must contain attribute and map
#>  information.  However, there are historical forest cover polygon infomration where maps
#>  are not available.  Forest Cover is provided at three critical milestones of at
#>  harvesting, at regeneration, and at free growing.  This is a part o fthe Silviculture
#>  and Land Status Tracking dataset, which includes tracking achievement of silviculture
#>  obligations on Crown Land
#> Available Resources (1):
#>  1. WMS getCapabilities request (wms)
#> Access the full 'Resources' data frame using:
#>  bcdc_tidy_resources('258bb088-4113-47b1-b568-ce20bd64e3e3')
#> Query and filter this data using:
#>  bcdc_query_geodata('258bb088-4113-47b1-b568-ce20bd64e3e3')
```

We see that this is a [Web Feature Service-enabled](https://en.wikipedia.org/wiki/Web_Feature_Service) geospatial data set--the list of data resources includes `WMS getCapabilities request`--so we can query and retrieve this geospatial data set using `bcdc_query_geodata()`:


```r
# Query the data using the permanent ID of the record to guard against name changes
bcdc_query_geodata("258bb088-4113-47b1-b568-ce20bd64e3e3")
#> Querying 'results-forest-cover-silviculture' record
#> • Using collect() on this object will return 898559 features and 159 fields
#> • Accessing this record requires pagination and will make 899 separate requests to the WFS.
#> • See ?bcdc_options
#> • At most six rows of the record are printed here
#> ───────────────────────────────────────────────────────────────────────────────────────────────────────
#> Simple feature collection with 6 features and 159 fields
#> Geometry type: MULTIPOLYGON
#> Dimension:     XY
#> Bounding box:  xmin: 1184167 ymin: 526455.5 xmax: 1801754 ymax: 1083545
#> Projected CRS: NAD83 / BC Albers
#> # A tibble: 6 × 160
#>   id                 FOREST_COVER_ID STOCKING_STANDARD… OPENING_ID STANDARDS_UNIT_… SILV_POLYGON_NU…
#>   <chr>                        <int>              <int>      <int> <chr>            <chr>           
#> 1 WHSE_FOREST_VEGET…         4177991            2243439    1724262 1                C               
#> 2 WHSE_FOREST_VEGET…         3994007                 NA    1248495 <NA>             91              
#> 3 WHSE_FOREST_VEGET…         3994067                 NA    1120935 <NA>             86              
#> 4 WHSE_FOREST_VEGET…         3994009            1404558    1248495 A                AA              
#> 5 WHSE_FOREST_VEGET…         3994057            1211894    1120935 B                2B              
#> 6 WHSE_FOREST_VEGET…         3994063                 NA    1120935 <NA>             82              
#> # … with 154 more variables: SILV_POLYGON_AREA <dbl>, SILV_POLYGON_NET_AREA <dbl>,
#> #   SILV_NON_MAPPED_AREA <int>, STOCKING_STATUS_CODE <chr>, STOCKING_TYPE_CODE <chr>,
#> #   STOCKING_CLASS_CODE <chr>, SILV_RESERVE_CODE <chr>, SILV_RESERVE_OBJECTIVE_CODE <chr>,
#> #   TREE_COVER_PATTERN_CODE <chr>, REENTRY_YEAR <chr>, REFERENCE_YEAR <int>, SITE_INDEX <int>,
#> #   SITE_INDEX_SOURCE_CODE <chr>, BGC_ZONE_CODE <chr>, BGC_SUBZONE_CODE <chr>, BGC_VARIANT <chr>,
#> #   BGC_PHASE <chr>, BEC_SITE_SERIES <chr>, BEC_SITE_TYPE <chr>, BEC_SERAL <chr>,
#> #   IS_SILV_IMPLIED_IND <chr>, FOREST_COVER_SILV_TYPE <chr>, S_FOREST_COVER_LAYER_ID <int>, …
```

This query shows that this data set has many features and over 150 fields.  Each feature is a treatment unit within a harvested opening, and contains information on the leading five tree species that are present in each treatment unit, including stems per hectare, age, and height.

Note that we have only _queried_ the data set so far---the data set would be too large (~1GB) to download efficiently.  So, let's use `filter()` with `bcdc_query_geodata()` to refine our query _before_ we collect the data and import it into R.


## Refining a Geospatial Data Query

To address our question, we need the treatment data (1) from the Prince George Natural Resource District _and_ (2) that contain western larch.

First, we can use the `bcdata` package to download the spatial boundary for the Prince George Natural Resource District&mdash;`DPG` is the `ORG_UNIT` for Prince George Natural Resource District:


```r
## Create a spatial feature object named dpg
dpg <- bcdc_query_geodata("natural-resource-nr-district") %>%
  filter(ORG_UNIT=="DPG") %>% # filter for Prince George Natural Resource District
  collect() # and collect the data
```

Let's plot this spatial object and double check we have we what we need:


```r
dpg %>%
  ggplot() +
  geom_sf() +
  theme_minimal()
```

<img src="vignette-fig-plot-dpg-1.png" title="plot of chunk plot-dpg" alt="plot of chunk plot-dpg" width="100%" />

Now we have a spatial object that we can use as a bounding box to filter and download records in the RESULTS - silviculture layer from the Prince George Natural Resource District.

We only need to download the treatments that have western larch planted. We can use the `bcdc_describe_feature()` helper function to examine the column names and types of the layer. In this case, we want to keep rows where the five `S_SPECIES_CODE_*` columns contain `"LW"`, the code for western larch.


```r
# Make a vector of tree species we are interested in
# (in this case only LW for western larch)
spp_list = c("LW")

# Query and filter the data layer
trees_dpg <-
  bcdc_query_geodata("258bb088-4113-47b1-b568-ce20bd64e3e3") %>%
  filter(INTERSECTS(dpg)) %>% #filter for records that are within the DPG boundary
  filter(
      S_SPECIES_CODE_1 %in% spp_list |  #filter for LW records
      S_SPECIES_CODE_2 %in% spp_list |
      S_SPECIES_CODE_3 %in% spp_list |
      S_SPECIES_CODE_4 %in% spp_list |
      S_SPECIES_CODE_5 %in% spp_list
  ) %>%
  collect() #collect/download the data
```


## Exploring the Data

Let's look at the dimensions of this now much more manageable data object we have downloaded from the B.C. Data Catalogue:


```r
dim(trees_dpg)
#> [1] 159 160
```

We can see there are several treatment units planted with western larch, and we can make a quick map of these harvested openings for the Prince George Natural Resource District:


```r
trees_dpg %>%
  ggplot() +
  geom_sf() +
  geom_sf(data = dpg, fill = NA) + #add the DPG spatial boundary
  theme_minimal()
```

<img src="vignette-fig-map-larch-plantations-dpg-1.png" title="plot of chunk map-larch-plantations-dpg" alt="plot of chunk map-larch-plantations-dpg" width="100%" />


We can also create some quick descriptive summaries of the data, treating the geospatial attribute table as a data frame in R, and answer our original question---how much western larch has been planted in the Prince George Natural Resource District?


#### What is the size and age distribution of larch plantations in the Prince George Natural Resource District in the year 2020?


```r
trees_dpg %>%
  mutate(age = 2020 - REFERENCE_YEAR + S_SPECIES_AGE_1) %>% #create a plantation age column
  ggplot() +  #start a plot
  aes(x = age, y = FEATURE_AREA_SQM/10000) + #convert feature area to hectares
  geom_bar(stat = "sum") +
  scale_x_continuous(name = "Stand Age (years)",
                     limits = c(2, 30),
                     breaks = seq(from = 5, to = 30, by = 5)) +
  labs(y = "Sum Treatment Unit Area (ha)",
       title = "Area Planted with Western Larch by Stand Age\nin Prince George Natural Resource District",
       caption = paste0("Data sourced from the B.C. Data Catalogue\n on ",
                        Sys.Date(),
                        " using the bcdata R package")) +
  theme_minimal() +
  theme(legend.position = "none")
```

<img src="vignette-fig-unnamed-chunk-1-1.png" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" width="100%" />

#### What is the Biogeoclimatic Ecosystem Classification (BEC) distribution of western larch plantations in the Prince George Natural Resource District?

We can download [British Columbia biogeoclimatic (BEC) data](https://catalogue.data.gov.bc.ca/dataset/bec-map) from the B.C. Data Catalogue using `bcdata` and join to our existing `trees_dpg` geospatial dataset using the `st_join` function from the `sf` package. You can learn more about `sf` [here](https://r-spatial.github.io/sf/).



```r
library(sf) #load the sf package

# Load the BEC data for Prince George Natural Resource District
bgc_dpg <- bcdc_query_geodata("WHSE_FOREST_VEGETATION.BEC_BIOGEOCLIMATIC_POLY") %>%
  filter(BBOX(dpg, crs = `EPSG:3005`)) %>% #note filtering with a BoundingBox
  collect()

# Join the BEC data with our tress_dpg geospatial data set
trees_bec_dpg <- trees_dpg %>%
  st_join(bgc_dpg[, "MAP_LABEL"]) #join BEC data for each polygon
```

Now, we can summarize the area planted with western larch by biogeoclimatic unit:


```r
trees_bec_dpg %>%
  group_by(MAP_LABEL) %>% # group polygons by biogeoclimatic unit
  summarise(Area = sum(FEATURE_AREA_SQM)/10000) %>%
  ggplot() +
  aes(x = MAP_LABEL, y = Area) +
  geom_col() +
  labs(y = "Sum Treatment Unit Area (ha)",
       x = "Biogeoclimatic Unit",
       title = "Area Planted with Western Larch by Biogeoclimatic Unit\nin Prince George Natural Resource District",
       caption = paste0("Data sourced from the B.C. Data Catalogue\n on ",
                        Sys.Date(),
                        " using the bcdata R package")) +
  theme_minimal() +
  theme(legend.position = "none")
```

<img src="vignette-fig-unnamed-chunk-3-1.png" title="plot of chunk unnamed-chunk-3" alt="plot of chunk unnamed-chunk-3" width="100%" />
---
title: "Get Started with bcdata"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Get Started with bcdata}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

<!--
Copyright 2019 Province of British Columbia

Licensed under the Apache License, Version 2.0 (the "License");
you may not use this file except in compliance with the License.
You may obtain a copy of the License at

http://www.apache.org/licenses/LICENSE-2.0

Unless required by applicable law or agreed to in writing, software distributed under the License is distributed on an "AS IS" BASIS,
WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
See the License for the specific language governing permissions and limitations under the License.
-->




The `bcdata` [R](https://www.r-project.org/) package contains functions for searching & retrieving data from the [B.C. Data Catalogue]( https://catalogue.data.gov.bc.ca).

The [B.C. Data Catalogue](https://www2.gov.bc.ca/gov/content?id=79B5224167334667A44C9E8B5143D0C5) is the place to find British Columbia Government data, applications and web services. Much of the data are released under the [Open Government Licence --- British Columbia](https://www2.gov.bc.ca/gov/content/data/open-data/open-government-licence-bc), as well as numerous other [licences](https://catalogue.data.gov.bc.ca/dataset?download_audience=Public).



You can install `bcdata` directly from CRAN:


```r
install.packages("bcdata")

library(bcdata)
```

### `bcdc_browse()`

`bcdata::bcdc_browse()` let's you access the [B.C. Data Catalogue web interface](https://catalogue.data.gov.bc.ca) directly from R---opening the catalogue search page in your default browser:


```r
## Take me to the B.C. Data Catalogue home page
bcdc_browse()
```

If you know the catalogue "human-readable" record name or permanent ID you can open directly to the record web page:


```r
## Take me to the B.C. Winery Locations catalogue record using the record name
bcdc_browse("bc-winery-locations")

## Take me to the B.C. Winery Locations catalogue record using the record permanent ID
bcdc_browse("1d21922b-ec4f-42e5-8f6b-bf320a286157")
```

### `bcdc_search()`

`bcdc_search()` let's you search records in the B.C. Data Catalogue, returning the search results in your R session.

Let's search the catalogue for records that contain the word "recycling":


```r
## Give me the catalogue search results for 'recycling'
bcdc_search("recycling")
#> List of B.C. Data Catalogue Records
#> Number of records: 3
#> Titles:
#> 1: BC FIRST Tire Recycling Data 1991-2006 (csv)
#>  ID: a29ad492-29a2-44b9-8693-d27a8cc8e686
#>  Name: bc-first-tire-recycling-data-1991-2006
#> 2: Tire Stewardship BC Tire Recycling Data (csv)
#>  ID: f791329b-c2dc-4f82-9993-209780f2a1c6
#>  Name: tire-stewardship-bc-tire-recycling-data
#> 3: Environmental Protection Information Resources e-Library (other)
#>  ID: dae0f2c3-b4f4-4d16-a96d-d7fe7c1581f3
#>  Name: environmental-protection-information-resources-e-library
#> 
#> Access a single record by calling `bcdc_get_record(ID)` with the ID from the desired
#>  record.
```

You can set the number of records to be returned from the search and/or you can customize your search using the catalogue search _facets_ `license_id`, `download_audience`, `res_format`, `sector`, and `organization`:


```r
## Give me the first catalogue search result for 'recycling'
bcdc_search("recycling", n = 1)
#> List of B.C. Data Catalogue Records
#> Number of records: 1
#> Titles:
#> 1: BC FIRST Tire Recycling Data 1991-2006 (csv)
#>  ID: a29ad492-29a2-44b9-8693-d27a8cc8e686
#>  Name: bc-first-tire-recycling-data-1991-2006
#> 
#> Access a single record by calling `bcdc_get_record(ID)` with the ID from the desired
#>  record.

## Give me the catalogue search results for 'recycling' where the
## data is licenced under Open Government Licence – British Columbia
bcdc_search("recycling", license_id = "2")
#> List of B.C. Data Catalogue Records
#> Number of records: 1
#> Titles:
#> 1: BC FIRST Tire Recycling Data 1991-2006 (csv)
#>  ID: a29ad492-29a2-44b9-8693-d27a8cc8e686
#>  Name: bc-first-tire-recycling-data-1991-2006
#> 
#> Access a single record by calling `bcdc_get_record(ID)` with the ID from the desired
#>  record.
```

You can see all valid values for the catalogue search facets using `bcdata::bcdc_search_facets()`:


```r
## Valid values for search facet 'license_id'
bcdc_search_facets(facet = "license_id")
#>         facet count                                             display_name name
#> 1  license_id    62                           Statistics Canada Open Licence   21
#> 2  license_id     1               Queen's Printer Licence - British Columbia   25
#> 3  license_id    12                      Open Government Licence – TransLink   48
#> 4  license_id    12 Open Government Licence – Municipality of North Cowichan   44
#> 5  license_id     4                 Open Government Licence - Destination BC   43
#> 6  license_id    59                         Open Government Licence - Canada   24
#> 7  license_id  1566               Open Government Licence - British Columbia    2
#> 8  license_id     2                  Open Government Licence - BC Assessment   47
#> 9  license_id     4                   Open Data Licence for ICBC Information   49
#> 10 license_id     2 Open Data Commons - Public Domain Dedication and Licence   45
#> 11 license_id    15                           Elections BC Open Data Licence   42
#> 12 license_id  1395                                              Access Only   22
```

Finally, you can retrieve the _metadata_ for a single catalogue record by using the record name or permanent ID with `bcdc_get_record()`. It is advised to use the permanent ID rather than the human-readable name in non-interactive situations---like scripts---to guard against future name changes of a record:


```r
## Give me the catalogue record metadata for `bc-first-tire-recycling-data-1991-2006`
bcdc_get_record("a29ad492-29a2-44b9-8693-d27a8cc8e686")
#> B.C. Data Catalogue Record: BC FIRST Tire Recycling Data 1991-2006
#> Name: bc-first-tire-recycling-data-1991-2006 (ID: a29ad492-29a2-44b9-8693-d27a8cc8e686)
#> Permalink: https://catalogue.data.gov.bc.ca/dataset/a29ad492-29a2-44b9-8693-d27a8cc8e686
#> Licence: Open Government Licence - British Columbia
#> Description: Financial Incentives for Recycling Scrap Tires (FIRST) collection and
#>  recycling data (tonnes) from 1991 to 2006. In 2007 [Tire Stewardship
#>  BC](http://www.tsbc.ca/), a not for profit society, launched the new scrap tire
#>  recycling program replacing the government-run program that had been in place since
#>  1991. Tire Stewardship BC collection and recycling data is available
#>  [here](https://catalogue.data.gov.bc.ca/dataset/f791329b-c2dc-4f82-9993-209780f2a1c6).
#> Available Resources (1):
#>  1. BC FIRST Tire Recycling Data 1991-2006 (csv)
#> Access the full 'Resources' data frame using:
#>  bcdc_tidy_resources('a29ad492-29a2-44b9-8693-d27a8cc8e686')
```

### `bcdc_get_data()`

Once you have located the B.C. Data Catalogue record with the data you want, you can use `bcdata::bcdc_get_data()` to download and read the data from the record.  You can use the record name, permanent ID or the result from `bcdc_get_record()`. Let's look at the B.C. Highway Web Cameras data:


```r
## Get the data resource for the `bc-highway-cams` catalogue record
bcdc_get_data("bc-highway-cams")
#> # A tibble: 917 × 19
#>    links_bchighwaycam  links_imageDisplay links_imageThumbn… links_replayTheDay    id highway_number
#>    <chr>               <chr>              <chr>              <chr>              <dbl> <chr>         
#>  1 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     2 5             
#>  2 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     5 3             
#>  3 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     6 16            
#>  4 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     7 1             
#>  5 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     8 1             
#>  6 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     9 19            
#>  7 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    10 97            
#>  8 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    11 1             
#>  9 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    12 1             
#> 10 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    13 99            
#> # … with 907 more rows, and 13 more variables: highway_locationDescription <chr>, camName <chr>,
#> #   caption <chr>, credit <chr>, orientation <chr>, latitude <dbl>, longitude <dbl>,
#> #   imageStats_updatePeriodMean <chr>, imageStats_updatePeriodStdDev <dbl>, markedDelayed <dbl>,
#> #   updatePeriodMean <dbl>, updatePeriodStdDev <dbl>, fetchMean <dbl>

## OR use the permanent ID, which is better for scripts or non-interactive use
bcdc_get_data("6b39a910-6c77-476f-ac96-7b4f18849b1c")
#> # A tibble: 917 × 19
#>    links_bchighwaycam  links_imageDisplay links_imageThumbn… links_replayTheDay    id highway_number
#>    <chr>               <chr>              <chr>              <chr>              <dbl> <chr>         
#>  1 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     2 5             
#>  2 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     5 3             
#>  3 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     6 16            
#>  4 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     7 1             
#>  5 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     8 1             
#>  6 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     9 19            
#>  7 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    10 97            
#>  8 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    11 1             
#>  9 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    12 1             
#> 10 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    13 99            
#> # … with 907 more rows, and 13 more variables: highway_locationDescription <chr>, camName <chr>,
#> #   caption <chr>, credit <chr>, orientation <chr>, latitude <dbl>, longitude <dbl>,
#> #   imageStats_updatePeriodMean <chr>, imageStats_updatePeriodStdDev <dbl>, markedDelayed <dbl>,
#> #   updatePeriodMean <dbl>, updatePeriodStdDev <dbl>, fetchMean <dbl>

## OR use the result from bcdc_get_record()
my_record <- bcdc_get_record("6b39a910-6c77-476f-ac96-7b4f18849b1c")
bcdc_get_data(my_record)
#> # A tibble: 917 × 19
#>    links_bchighwaycam  links_imageDisplay links_imageThumbn… links_replayTheDay    id highway_number
#>    <chr>               <chr>              <chr>              <chr>              <dbl> <chr>         
#>  1 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     2 5             
#>  2 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     5 3             
#>  3 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     6 16            
#>  4 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     7 1             
#>  5 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     8 1             
#>  6 http://images.driv… http://images.dri… http://images.dri… http://images.dri…     9 19            
#>  7 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    10 97            
#>  8 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    11 1             
#>  9 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    12 1             
#> 10 http://images.driv… http://images.dri… http://images.dri… http://images.dri…    13 99            
#> # … with 907 more rows, and 13 more variables: highway_locationDescription <chr>, camName <chr>,
#> #   caption <chr>, credit <chr>, orientation <chr>, latitude <dbl>, longitude <dbl>,
#> #   imageStats_updatePeriodMean <chr>, imageStats_updatePeriodStdDev <dbl>, markedDelayed <dbl>,
#> #   updatePeriodMean <dbl>, updatePeriodStdDev <dbl>, fetchMean <dbl>
```

A catalogue record can have one or multiple data files---or "resources". If there is only one resource, `bcdc_get_data()` will return that resource by default, as in the above `bc-highway-cams` example. If there are multiple data resources you will need to specify which resource you want. Let's look at a catalogue record that contains multiple data resources---BC Schools - Programs Offered in Schools:


```r
## Get the record ID for the `bc-schools-programs-offered-in-schools` catalogue record
bcdc_search("school programs", n = 1)
#> List of B.C. Data Catalogue Records
#> Number of records: 1
#> Titles:
#> 1: BC Schools - Programs Offered in Schools (txt, xlsx)
#>  ID: b1f27d1c-244a-410e-a361-931fac62a524
#>  Name: bc-schools-programs-offered-in-schools
#> 
#> Access a single record by calling `bcdc_get_record(ID)` with the ID from the desired
#>  record.

## Get the metadata for the `bc-schools-programs-offered-in-schools` catalogue record
bcdc_get_record("b1f27d1c-244a-410e-a361-931fac62a524")
#> B.C. Data Catalogue Record: BC Schools - Programs Offered in Schools
#> Name: bc-schools-programs-offered-in-schools (ID: b1f27d1c-244a-410e-a361-931fac62a524)
#> Permalink: https://catalogue.data.gov.bc.ca/dataset/b1f27d1c-244a-410e-a361-931fac62a524
#> Licence: Open Government Licence - British Columbia
#> Description: BC Schools English Language Learners, French Immersion, Francophone, Career
#>  Preparation, Aboriginal Support Services, Aboriginal Language and Culture, Continuing
#>  Education and Career Technical Programs offered in BC schools up to 2013/2014.
#> Available Resources (2):
#>  1. ProgramsOfferedinSchools.txt (txt)
#>  2. ProgramsOfferedinSchools.xlsx (xlsx)
#> Access the full 'Resources' data frame using:
#>  bcdc_tidy_resources('b1f27d1c-244a-410e-a361-931fac62a524')
```

We see there are two data files or resources available in this record, so we need to tell `bcdc_get_data()` which one we want. When used interactively, `bcdc_get_data()` will prompt you with the list of available resources through `bcdata` and ask you to select the resource you want. The resource ID for each data set is available _in_ the metadata record ☝️:


```r
## Get the txt data resource from the `bc-schools-programs-offered-in-schools`
## catalogue record
bcdc_get_data("b1f27d1c-244a-410e-a361-931fac62a524", resource = 'a393f8cf-51ec-42c6-8449-4cea4c75385c')
#> # A tibble: 16,152 × 24
#>    `Data Level` `School Year` `Facility Type` `Public Or Independ… `District Numbe… `District Name` 
#>    <chr>        <chr>         <chr>           <chr>                <chr>            <chr>           
#>  1 SCHOOL LEVEL 2005/2006     STANDARD        BC Public School     005              Southeast Koote…
#>  2 SCHOOL LEVEL 2006/2007     STANDARD        BC Public School     005              Southeast Koote…
#>  3 SCHOOL LEVEL 2007/2008     STANDARD        BC Public School     005              Southeast Koote…
#>  4 SCHOOL LEVEL 2005/2006     STANDARD        BC Public School     005              Southeast Koote…
#>  5 SCHOOL LEVEL 2006/2007     STANDARD        BC Public School     005              Southeast Koote…
#>  6 SCHOOL LEVEL 2007/2008     STANDARD        BC Public School     005              Southeast Koote…
#>  7 SCHOOL LEVEL 2008/2009     STANDARD        BC Public School     005              Southeast Koote…
#>  8 SCHOOL LEVEL 2009/2010     STANDARD        BC Public School     005              Southeast Koote…
#>  9 SCHOOL LEVEL 2010/2011     STANDARD        BC Public School     005              Southeast Koote…
#> 10 SCHOOL LEVEL 2011/2012     STANDARD        BC Public School     005              Southeast Koote…
#> # … with 16,142 more rows, and 18 more variables: School Number <chr>, School Name <chr>,
#> #   Has Eng Lang Learner Prog <lgl>, Has Core French <lgl>, Has Early French Immersion <lgl>,
#> #   Has Late French Immersion <lgl>, Has Prog Francophone <lgl>,
#> #   Has Any French Immersion Prog <lgl>, Has Any French Prog <lgl>, Has Aborig Supp Services <lgl>,
#> #   Has Other Appr Aborig Prog <lgl>, Has Aborig Lang And Cult <lgl>, Has Continuing Ed Prog <lgl>,
#> #   Has Distributed Learn Prog <lgl>, Has Career Prep Prog <lgl>, Has Coop Prog <lgl>,
#> #   Has Apprenticeship Prog <lgl>, Has Career Technical Prog <lgl>
```

Alternatively, you can retrieve the full details of the available resources for a given record as a data frame using `bcdc_tidy_resources()`:


```r
## Get a data frame of data resources for the `bc-schools-programs-offered-in-schools`
## catalogue record
bcdc_tidy_resources("b1f27d1c-244a-410e-a361-931fac62a524")
#> # A tibble: 2 × 9
#>   name     url           id       format ext   package_id   location  wfs_available bcdata_available
#>   <chr>    <chr>         <chr>    <chr>  <chr> <chr>        <chr>     <lgl>         <lgl>           
#> 1 Program… http://www.b… a393f8c… txt    txt   b1f27d1c-24… catalogu… FALSE         TRUE            
#> 2 Program… http://www.b… 1e34098… xlsx   xlsx  b1f27d1c-24… catalogu… FALSE         TRUE
```

`bcdc_get_data()` will also detect if the data resource is a geospatial file, and automatically reads and returns it as an [`sf` object](https://r-spatial.github.io/sf/) in your R session.

Let's get the air zones for British Columbia:


```r
## Find the B.C. Air Zones catalogue record
bcdc_search("air zones", res_format = "geojson")
#> List of B.C. Data Catalogue Records
#> Number of records: 1
#> Titles:
#> 1: British Columbia Air Zones (shp, kml, geojson)
#>  ID: e8eeefc4-2826-47bc-8430-85703d328516
#>  Name: british-columbia-air-zones
#> 
#> Access a single record by calling `bcdc_get_record(ID)` with the ID from the desired
#>  record.

## Get the metadata for the B.C. Air Zones catalogue record
bc_az_metadata <- bcdc_get_record("e8eeefc4-2826-47bc-8430-85703d328516")

## Get the B.C. Air Zone geospatial data
bc_az <- bcdc_get_data(bc_az_metadata, resource = "c495d082-b586-4df0-9e06-bd6b66a8acd9")

## Plot the B.C. Air Zone geospatial data with ggplot()
bc_az %>%
  ggplot() +
  geom_sf() +
  theme_minimal()
```

<img src="vignette-fig-air_zones-1.png" title="plot of chunk air_zones" alt="plot of chunk air_zones" width="100%" />


**Note:** The `bcdata` package supports downloading _most_ file types, including zip archives. It will do its best to identify and read data from
zip files, however if there are multiple data files in the zip, or data files that `bcdata` doesn't know how to import, it will fail.


### `bcdc_query_geodata()`

Many geospatial data sets in the B.C. Data Catalogue are available through a [Web Feature Service](https://en.wikipedia.org/wiki/Web_Feature_Service). While `bcdc_get_data()` will retrieve the geospatial data for you, sometimes the geospatial file is very large---and slow to download---and/or you may only want _some_ of the data. `bcdc_query_geodata()` let's you query catalogue geospatial data available as a Web Feature Service using `select` and `filter` functions (just like in [`dplyr`](https://dplyr.tidyverse.org/). The `bcdc::collect()` function returns the `bcdc_query_geodata()` query results as an [`sf` object](https://r-spatial.github.io/sf/) in your R session.

Let's get the Capital Regional District boundary from the [B.C. Regional Districts geospatial data](https://catalogue.data.gov.bc.ca/dataset/d1aff64e-dbfe-45a6-af97-582b7f6418b9)---the whole file takes 30-60 seconds to download and I only need the one polygon, so why not save some time:


```r
## Find the B.C. Regional Districts catalogue record
bcdc_search("regional districts administrative areas", res_format = "wms", n = 1)
#> List of B.C. Data Catalogue Records
#> Number of records: 1
#> Titles:
#> 1: Regional Districts - Legally Defined Administrative Areas of BC (other, xlsx, wms,
#>  kml)
#>  ID: d1aff64e-dbfe-45a6-af97-582b7f6418b9
#>  Name: regional-districts-legally-defined-administrative-areas-of-bc
#> 
#> Access a single record by calling `bcdc_get_record(ID)` with the ID from the desired
#>  record.

## Get the metadata for the B.C. Regional Districts catalogue record
bc_regional_districts_metadata <- bcdc_get_record("d1aff64e-dbfe-45a6-af97-582b7f6418b9")

## We can see in the search results, and in the metadata record, that this record has a `"wms"`
## resource format, indicating that it is available as a Web Feature Service and thus
## we can query it using `bcdc_query_geodata()`

## Have a quick look at the geospatial columns to help with filter or select
bcdc_describe_feature(bc_regional_districts_metadata)
#> # A tibble: 21 × 5
#>    col_name                 sticky remote_col_type local_col_type column_comments                   
#>    <chr>                    <lgl>  <chr>           <chr>          <chr>                             
#>  1 id                       TRUE   xsd:string      character       <NA>                             
#>  2 LGL_ADMIN_AREA_ID        TRUE   xsd:decimal     numeric        "An operationally-generated uniqu…
#>  3 ADMIN_AREA_NAME          FALSE  xsd:string      character      "The authoritative, officially ap…
#>  4 ADMIN_AREA_ABBREVIATION  FALSE  xsd:string      character      "A short form or commonly-known a…
#>  5 ADMIN_AREA_BOUNDARY_TYPE FALSE  xsd:string      character      "BOUNDARY TYPE is a high-level gr…
#>  6 ADMIN_AREA_GROUP_NAME    FALSE  xsd:string      character      "The name given to the larger adm…
#>  7 CHANGE_REQUESTED_ORG     FALSE  xsd:string      character      "The government acronym of the Mi…
#>  8 UPDATE_TYPE              FALSE  xsd:string      character      "A short description of the lates…
#>  9 WHEN_UPDATED             FALSE  xsd:date        date           "The date and time the record was…
#> 10 MAP_STATUS               FALSE  xsd:string      character      "That the digital map has been ap…
#> # … with 11 more rows

## Get the Capital Regional District polygon from the B.C. Regional
## Districts geospatial data
my_regional_district <- bcdc_query_geodata(bc_regional_districts_metadata) %>%
  filter(ADMIN_AREA_NAME == "Capital Regional District") %>%
  collect()

## Plot the Capital Regional District polygon with ggplot()
my_regional_district  %>%
  ggplot() +
  geom_sf() +
  theme_minimal()
```

<img src="vignette-fig-regional_districts-1.png" title="plot of chunk regional_districts" alt="plot of chunk regional_districts" width="100%" />

The vignette [Querying Spatial Data with bcdata](https://bcgov.github.io/bcdata/articles/efficiently-query-spatial-data-in-the-bc-data-catalogue.html) provides a full demonstration on how to use `bcdata::bcdc_query_geodata()` to fine tune a [Web Feature Service](https://www2.gov.bc.ca/gov/content?id=95D78D544B244F34B89223EF069DF74E) request for geospatial data from the B.C. Data Catalogue.
---
title: "bcdata Service Documentation"
date: "2021-10-27"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{bcdata Service Documentation}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---


# DataBC Services used by `bcdata`

This document is an attempt at a comprehensive list of services and API endpoints accessed by the bcdata R package, as well which return values we rely on from those endpoints.

## BC Data Catalogue

### API version: 
  - PROD: `https://catalogue.data.gov.bc.ca/api/3`
  - BETA: `https://beta-catalogue.data.gov.bc.ca/api/3`

### Endpoints:
  - `/action/package_show`
  - `/action/package_search`
    - `license_id`
    - `download_audience`
    - `res_format`
    - `sector`
    - `organization`
  - `/action/package_list`
  - `/action/group_show`

### Response values used:
  - package:
    - `title`
    - `name`
    - `id`
    - `license_title`
    - `type`
    - `notes`
    - `layer_name`
    - `resources` (see below)
    
  - resource:
    - `id`
    - `package_id`
    - `object_name`
        (This is not always the same as the `typeNames` parameter in `resource.url`, as that is sometimes a simplified view - eg., `WHSE_ADMIN_BOUNDARIES.ADM_NR_DISTRICTS_SPG` vs 
`WHSE_ADMIN_BOUNDARIES.ADM_NR_DISTRICTS_SP`)
    - `details`
      - `column_comments`
      - `column_name`
    - `bcdc_type` (not actually using yet but [may be useful](https://github.com/bcgov/bcdata/pull/283#issuecomment-924442166))
    - `format`
    - `resource_storage_location`
    - `name`
    - `url`
    
  - group: 
    - `description`
    - `packages`

## Web Services

### API Version:

  - TEST: https://test.openmaps.gov.bc.ca
  - DELIVERY: https://delivery.openmaps.gov.bc.ca
  - PROD: https://openmaps.gov.bc.ca
  
  Endpoints:
  
  - wfs: `geo/pub/wfs`
  - wms: `geo/pub/wms`

  Query Parameters for `geo/pub/wfs`:
  
  - query is sent in the body of a `POST` request (with `encode = "form"`). If a dataset has > n records (default n = 1000), pagination is used to send sequential requests. Pagination is executed using `count`, `sortBY`, and `startIndex`.
    - SERVICE = "WFS"
    - VERSION = "2.0.0"
    - REQUEST = "GetCapabilities"
    - REQUEST = "GetFeature"
      - outputFormat = "application/json"
      - typeNames (extracted from `resource.url` and compared against `resource.object_name`)
      - SRSNAME (default `EPSG:3005`)
      - CQL_FILTER
      - count
      - propertyName
      - sortBy
      - startIndex
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/describe-feature.R
\name{bcdc_describe_feature}
\alias{bcdc_describe_feature}
\title{Describe the attributes of a Web Feature Service}
\usage{
bcdc_describe_feature(record)
}
\arguments{
\item{record}{either a \code{bcdc_record} object (from the result of \code{bcdc_get_record()}),
a character string denoting the name or ID of a resource (or the URL) or a BC Geographic
Warehouse (BCGW) name.

It is advised to use the permanent ID for a record or the BCGW name rather than the
human-readable name to guard against future name changes of the record.
If you use the human-readable name a warning will be issued once per
session. You can silence these warnings altogether by setting an option:
\code{options("silence_named_get_data_warning" = TRUE)} - which you can set
in your .Rprofile file so the option persists across sessions.}
}
\value{
\code{bcdc_describe_feature} returns a tibble describing the attributes of a B.C. Data Catalogue record.
The tibble returns the following columns:
\itemize{
\item col_name: attributes of the feature
\item sticky: whether a column can be separated from the record in a Web Feature Service call via the \code{dplyr::select} method
\item remote_col_type: class of what is return by the web feature service
\item local_col_type: the column class in R
\item column_comments: additional metadata specific to that column
}
}
\description{
Describe the attributes of column of a record accessed through the Web Feature Service.
This can be a useful tool to examine a layer before issuing a query with \code{bcdc_query_geodata}.
}
\examples{
\donttest{
try(
  bcdc_describe_feature("bc-airports")
)

try(
  bcdc_describe_feature("WHSE_IMAGERY_AND_BASE_MAPS.GSR_AIRPORTS_SVW")
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc_options.R
\name{bcdc_options}
\alias{bcdc_options}
\title{Retrieve options used in bcdata, their value if set and the default value.}
\usage{
bcdc_options()
}
\description{
This function retrieves bcdata specific options that can be set. These options can be set
using \verb{option(\{name of the option\} = \{value of the option\})}. The default options are purposefully
set conservatively to hopefully ensure successful requests. Resetting these options may result in
failed calls to the data catalogue. Options in R are reset every time R is re-started. See examples for
addition ways to restore your initial state.
}
\details{
\code{bcdata.max_geom_pred_size} is the maximum size in bytes of an object used for a geometric operation. Objects
that are bigger than this value will have a bounding box drawn and apply the geometric operation
on that simpler polygon. The \link{bcdc_check_geom_size} function can be used to assess whether a given spatial object
exceed the value of this option. Users can iteratively try to increase the maximum geometric predicate size and see
if the bcdata catalogue accepts the request.

\code{bcdata.chunk_limit} is an option useful when dealing with very large data sets. When requesting large objects
from the catalogue, the request is broken up into smaller chunks which are then recombined after they've
been downloaded. This is called "pagination". bcdata does this all for you but using this option you can set the size of the chunk
requested. On faster internet connections, a bigger chunk limit could be useful while on slower connections,
it is advisable to lower the chunk limit. Chunks must be less than 10000.

\code{bcdata.single_download_limit} is the maximum number of records an object can be before forcing a paginated download
(see entry for \code{bcdata.chunk_limit} for details on pagination).
Tweaking this option in conjunction with \code{bcdata.chunk_limit} can often resolve failures in large and complex downloads.
The default is 10000 records.
}
\examples{
\donttest{
## Save initial conditions
try(
  original_options <- options()
)

## See initial options
try(
  bcdc_options()
)

try(
  options(bcdata.max_geom_pred_size = 1E6)
)

## See updated options
try(
  bcdc_options()
)

## Reset initial conditions
try(
 options(original_options)
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc_search.R
\name{bcdc_search}
\alias{bcdc_search}
\title{Search the B.C. Data Catalogue}
\usage{
bcdc_search(
  ...,
  license_id = NULL,
  download_audience = "Public",
  res_format = NULL,
  sector = NULL,
  organization = NULL,
  n = 100
)
}
\arguments{
\item{...}{search terms}

\item{license_id}{the type of license (see \code{bcdc_search_facets("license_id")}).}

\item{download_audience}{download audience
(see \code{bcdc_search_facets("download_audience")}). Default \code{"Public"}}

\item{res_format}{format of resource (see \code{bcdc_search_facets("res_format")})}

\item{sector}{sector of government from which the data comes
(see \code{bcdc_search_facets("sector")})}

\item{organization}{government organization that manages the data
(see \code{bcdc_search_facets("organization")})}

\item{n}{number of results to return. Default \code{100}}
}
\value{
A list containing the records that match the search
}
\description{
Search the B.C. Data Catalogue
}
\examples{
\donttest{
try(
  bcdc_search("forest")
)

try(
  bcdc_search("regional district", res_format = "fgdb")
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc_search.R
\name{bcdc_tidy_resources}
\alias{bcdc_tidy_resources}
\title{Provide a data frame containing the metadata for all resources from a single B.C. Data Catalogue record}
\usage{
bcdc_tidy_resources(record)
}
\arguments{
\item{record}{either a \code{bcdc_record} object (from the result of \code{bcdc_get_record()}),
a character string denoting the name or ID of a resource (or the URL) or a BC Geographic
Warehouse (BCGW) name.

It is advised to use the permanent ID for a record or the BCGW name rather than the
human-readable name to guard against future name changes of the record.
If you use the human-readable name a warning will be issued once per
session. You can silence these warnings altogether by setting an option:
\code{options("silence_named_get_data_warning" = TRUE)} - which you can set
in your .Rprofile file so the option persists across sessions.}
}
\value{
A data frame containing the metadata for all the resources for a record
}
\description{
Returns a rectangular data frame of all resources contained within a record. This is particularly useful
if you are trying to construct a vector of multiple resources in a record. The data frame also provides
useful information on the formats, availability and types of data available.
}
\examples{
\donttest{
try(
  airports <- bcdc_get_record("bc-airports")
)

try(
  bcdc_tidy_resources(airports)
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc-web-services.R
\name{bcdc_preview}
\alias{bcdc_preview}
\title{Get preview map from the B.C. Web Map Service}
\usage{
bcdc_preview(record)
}
\arguments{
\item{record}{either a \code{bcdc_record} object (from the result of \code{bcdc_get_record()}),
a character string denoting the name or ID of a resource (or the URL) or a BC Geographic
Warehouse (BCGW) name.

It is advised to use the permanent ID for a record or the BCGW name rather than the
human-readable name to guard against future name changes of the record.
If you use the human-readable name a warning will be issued once per
session. You can silence these warnings altogether by setting an option:
\code{options("silence_named_get_data_warning" = TRUE)} - which you can set
in your .Rprofile file so the option persists across sessions.}
}
\description{
Note this does not return the actual map features, rather
opens an image preview of the layer in a
\href{https://rstudio.github.io/leaflet/}{Leaflet} map window
}
\examples{
\donttest{
try(
  bcdc_preview("regional-districts-legally-defined-administrative-areas-of-bc")
)

try(
  bcdc_preview("points-of-well-diversion-applications")
)

# Using BCGW name
try(
  bcdc_preview("WHSE_LEGAL_ADMIN_BOUNDARIES.ABMS_REGIONAL_DISTRICTS_SP")
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc_search.R
\name{bcdc_list_groups}
\alias{bcdc_list_groups}
\alias{bcdc_list_group_records}
\title{Retrieve group information for B.C. Data Catalogue}
\usage{
bcdc_list_groups()

bcdc_list_group_records(group)
}
\arguments{
\item{group}{Name of the group}
}
\description{
Returns a tibble of groups or records. Groups can be viewed here:
https://catalogue.data.gov.bc.ca/group or accessed directly from R using \code{bcdc_list_groups}
}
\section{Functions}{
\itemize{
\item \code{bcdc_list_groups}: 
}}

\examples{
\donttest{
try(
  bcdc_list_group_records('environmental-reporting-bc')
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-classes.R, R/utils-mutate.R
\name{mutate.bcdc_promise}
\alias{mutate.bcdc_promise}
\alias{mutate}
\title{Throw an informative error when attempting mutate on a \code{bcdc_promise} object}
\usage{
\method{mutate}{bcdc_promise}(.data, ...)
}
\arguments{
\item{.data}{object of class \code{bcdc_promise} (likely passed from \code{\link[=bcdc_query_geodata]{bcdc_query_geodata()}})}

\item{...}{One or more unquoted expressions separated by commas. See details.}
}
\description{
The CQL syntax to generate WFS calls does not current allow arithmetic operations. Therefore
this function exists solely to generate an informative error that suggests an alternative
approach to use mutate with bcdata

See \code{dplyr::\link[dplyr]{mutate}} for details.
}
\section{Methods (by class)}{
\itemize{
\item \code{bcdc_promise}: mutate.bcdc_promise
}}

\examples{
\donttest{

## Mutate columns
try(
  bcdc_query_geodata("bc-airports") \%>\%
    mutate(LATITUDE * 100)
)
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc_search.R
\name{bcdc_get_record}
\alias{bcdc_get_record}
\title{Show a single B.C. Data Catalogue record}
\usage{
bcdc_get_record(id)
}
\arguments{
\item{id}{the human-readable name, permalink ID, or
URL of the record.

It is advised to use the permanent ID for a record rather than the
human-readable name to guard against future name changes of the record.
If you use the human-readable name a warning will be issued once per
session. You can silence these warnings altogether by setting an option:
\code{options("silence_named_get_record_warning" = TRUE)} - which you can put
in your .Rprofile file so the option persists across sessions.}
}
\value{
A list containing the metadata for the record
}
\description{
Show a single B.C. Data Catalogue record
}
\examples{
\donttest{
try(
  bcdc_get_record("https://catalogue.data.gov.bc.ca/dataset/bc-airports")
)

try(
  bcdc_get_record("bc-airports")
)

try(
  bcdc_get_record("https://catalogue.data.gov.bc.ca/dataset/76b1b7a3-2112-4444-857a-afccf7b20da8")
)

try(
  bcdc_get_record("76b1b7a3-2112-4444-857a-afccf7b20da8")
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-classes.R, R/utils-show-query.R
\name{show_query.bcdc_promise}
\alias{show_query.bcdc_promise}
\alias{show_query.bcdc_sf}
\alias{show_query}
\title{Show SQL and URL used for Web Feature Service request from B.C. Data Catalogue}
\usage{
\method{show_query}{bcdc_promise}(x, ...)

\method{show_query}{bcdc_sf}(x, ...)
}
\arguments{
\item{x}{object of class bcdc_promise or bcdc_sf}
}
\description{
Display Web Feature Service query CQL

See \code{dplyr::\link[dplyr:explain]{show_query}} for details.
}
\section{Methods (by class)}{
\itemize{
\item \code{bcdc_promise}: show_query.bcdc_promise

\item \code{bcdc_sf}: show_query.bcdc_promise
}}

\examples{
\donttest{
try(
  bcdc_query_geodata("bc-environmental-monitoring-locations") \%>\%
    filter(PERMIT_RELATIONSHIP == "DISCHARGE") \%>\%
    show_query()
)
  }

\donttest{
try(
  air <- bcdc_query_geodata("bc-airports") \%>\%
    collect()
)

try(
  show_query(air)
)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{bcdc_get_data}
\alias{bcdc_get_data}
\title{Download and read a resource from a B.C. Data Catalogue record}
\usage{
bcdc_get_data(record, resource = NULL, verbose = TRUE, ...)
}
\arguments{
\item{record}{either a \code{bcdc_record} object (from the result of \code{bcdc_get_record()}),
a character string denoting the name or ID of a resource (or the URL) or a BC Geographic
Warehouse (BCGW) name.

It is advised to use the permanent ID for a record or the BCGW name rather than the
human-readable name to guard against future name changes of the record.
If you use the human-readable name a warning will be issued once per
session. You can silence these warnings altogether by setting an option:
\code{options("silence_named_get_data_warning" = TRUE)} - which you can set
in your .Rprofile file so the option persists across sessions.}

\item{resource}{optional argument used when there are multiple data files
within the same record. See examples.}

\item{verbose}{When more than one resource is available for a record,
should extra information about those resources be printed to the console?
Default \code{TRUE}}

\item{...}{arguments passed to other functions. Tabular data is passed to a function to handle
the import based on the file extension. \code{\link[=bcdc_read_functions]{bcdc_read_functions()}} provides details on which functions
handle the data import. You can then use this information to look at the help pages of those functions.
See the examples for a workflow that illustrates this process.
For spatial Web Feature Service data the \code{...} arguments are passed to \code{\link[=bcdc_query_geodata]{bcdc_query_geodata()}}.}
}
\value{
An object of a type relevant to the resource (usually a tibble or an sf object)
}
\description{
Download and read a resource from a B.C. Data Catalogue record
}
\examples{
\donttest{
# Using the record and resource ID:
try(
  bcdc_get_data(record = '76b1b7a3-2112-4444-857a-afccf7b20da8',
                resource = '4d0377d9-e8a1-429b-824f-0ce8f363512c')
)

try(
  bcdc_get_data('1d21922b-ec4f-42e5-8f6b-bf320a286157')
)

# Using a `bcdc_record` object obtained from `bcdc_get_record`:
try(
  record <- bcdc_get_record('1d21922b-ec4f-42e5-8f6b-bf320a286157')
)

try(
  bcdc_get_data(record)
)

# Using a BCGW name
try(
  bcdc_get_data("WHSE_IMAGERY_AND_BASE_MAPS.GSR_AIRPORTS_SVW")
)

# Using sf's sql querying ability
try(
  bcdc_get_data(
    record = '30aeb5c1-4285-46c8-b60b-15b1a6f4258b',
    resource = '3d72cf36-ab53-4a2a-9988-a883d7488384',
    layer = 'BC_Boundary_Terrestrial_Line',
    query = "SELECT SHAPE_Length, geom FROM BC_Boundary_Terrestrial_Line WHERE SHAPE_Length < 100"
  )
)

## Example of correcting import problems

## Some initial problems reading in the data
try(
  bcdc_get_data('d7e6c8c7-052f-4f06-b178-74c02c243ea4')
)

## From bcdc_get_record we realize that the data is in xlsx format
try(
 bcdc_get_record('8620ce82-4943-43c4-9932-40730a0255d6')
)

## bcdc_read_functions let's us know that bcdata
## uses readxl::read_excel to import xlsx files
try(
 bcdc_read_functions()
)

## bcdata let's you know that this resource has
## multiple worksheets
try(
 bcdc_get_data('8620ce82-4943-43c4-9932-40730a0255d6')
)

## we can control what is read in from an excel file
## using arguments from readxl::read_excel
try(
  bcdc_get_data('8620ce82-4943-43c4-9932-40730a0255d6', sheet = 'Regional Districts')
)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc_search.R
\name{bcdc_list}
\alias{bcdc_list}
\title{Return a full list of the names of B.C. Data Catalogue records}
\usage{
bcdc_list()
}
\value{
A character vector of the names of B.C. Data Catalogue records
}
\description{
Return a full list of the names of B.C. Data Catalogue records
}
\examples{
\donttest{
try(
  bcdc_list()
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cql-geom-predicates.R
\name{CQL}
\alias{CQL}
\title{CQL escaping}
\usage{
CQL(...)
}
\arguments{
\item{...}{Character vectors that will be combined into a single CQL statement.}
}
\value{
An object of class \code{c("CQL", "SQL")}
}
\description{
Write a CQL expression to escape its inputs, and return a CQL/SQL object.
Used when writing filter expressions in \code{\link[=bcdc_query_geodata]{bcdc_query_geodata()}}.
}
\details{
See \href{https://docs.geoserver.org/stable/en/user/tutorials/cql/cql_tutorial.html}{the CQL/ECQL for Geoserver website}.
}
\examples{
CQL("FOO > 12 & NAME LIKE 'A&'")
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
Pipe operator
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cql-geom-predicates.R
\name{cql_geom_predicates}
\alias{cql_geom_predicates}
\alias{EQUALS}
\alias{DISJOINT}
\alias{INTERSECTS}
\alias{TOUCHES}
\alias{CROSSES}
\alias{WITHIN}
\alias{CONTAINS}
\alias{OVERLAPS}
\alias{BBOX}
\alias{DWITHIN}
\title{CQL Geometry Predicates}
\usage{
EQUALS(geom)

DISJOINT(geom)

INTERSECTS(geom)

TOUCHES(geom)

CROSSES(geom)

WITHIN(geom)

CONTAINS(geom)

OVERLAPS(geom)

BBOX(coords, crs = NULL)

DWITHIN(
  geom,
  distance,
  units = c("meters", "feet", "statute miles", "nautical miles", "kilometers")
)
}
\arguments{
\item{geom}{an \code{sf}/\code{sfc}/\code{sfg} or \code{bbox} object (from the \code{sf} package)}

\item{coords}{the coordinates of the bounding box as four-element numeric
vector \code{c(xmin, ymin, xmax, ymax)}, a \code{bbox} object from the \code{sf}
package (the result of running \code{sf::st_bbox()} on an \code{sf} object), or
an \code{sf} object which then gets converted to a bounding box on the fly.}

\item{crs}{(Optional) A numeric value or string containing an SRS code. If
\code{coords} is a \code{bbox} object with non-empty crs, it is taken from that.
(For example, \code{'EPSG:3005'} or just \code{3005}. The default is to use the CRS of
the queried layer)}

\item{distance}{numeric value for distance tolerance}

\item{units}{units that distance is specified in. One of
\code{"feet"}, \code{"meters"}, \code{"statute miles"}, \code{"nautical miles"}, \code{"kilometers"}}
}
\value{
a CQL expression to be passed on to the WFS call
}
\description{
Functions to construct a CQL expression to be used
to filter results from \code{\link[=bcdc_query_geodata]{bcdc_query_geodata()}}.
See \href{https://docs.geoserver.org/stable/en/user/filter/ecql_reference.html#spatial-predicate}{the geoserver CQL documentation for details}.
The sf object is automatically converted in a
bounding box to reduce the complexity of the Web Feature Service call. Subsequent in-memory
filtering may be needed to achieve exact results.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cql-translator.R
\docType{class}
\name{wfsConnection-class}
\alias{wfsConnection-class}
\alias{dbQuoteIdentifier,wfsConnection,ANY-method}
\alias{dbQuoteString,wfsConnection,ANY-method}
\title{wfsConnection class}
\usage{
\S4method{dbQuoteIdentifier}{wfsConnection,ANY}(conn, x)

\S4method{dbQuoteString}{wfsConnection,ANY}(conn, x)
}
\description{
wfsConnection class
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-classes.R, R/utils-select.R
\name{select.bcdc_promise}
\alias{select.bcdc_promise}
\alias{select}
\title{Select columns from bcdc_query_geodata() call}
\usage{
\method{select}{bcdc_promise}(.data, ...)
}
\arguments{
\item{.data}{object of class \code{bcdc_promise} (likely passed from \code{\link[=bcdc_query_geodata]{bcdc_query_geodata()}})}

\item{...}{One or more unquoted expressions separated by commas. See details.}
}
\description{
Similar to a \code{dplyr::select} call, this allows you to select which columns you want the Web Feature Service to return.
A key difference between \code{dplyr::select} and \code{bcdata::select} is the presence of "sticky" columns that are
returned regardless of what columns are selected. If any of these "sticky" columns are selected
only "sticky" columns are return. \code{bcdc_describe_feature} is one way to tell if columns are sticky in advance
of issuing the Web Feature Service call.

See \code{dplyr::\link[dplyr]{select}} for details.
}
\section{Methods (by class)}{
\itemize{
\item \code{bcdc_promise}: select.bcdc_promise
}}

\examples{
\donttest{
try(
  feature_spec <- bcdc_describe_feature("bc-airports")
)

try(
  ## Columns that can selected:
  feature_spec[feature_spec$sticky == TRUE,]
)

## Select columns
try(
  bcdc_query_geodata("bc-airports") \%>\%
    select(DESCRIPTION, PHYSICAL_ADDRESS)
)

## Select "sticky" columns
try(
  bcdc_query_geodata("bc-airports") \%>\%
    select(LOCALITY)
)
}


}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-as_tibble.R, R/utils-classes.R,
%   R/utils-collect.R
\name{as_tibble}
\alias{as_tibble}
\alias{collect.bcdc_promise}
\alias{as_tibble.bcdc_promise}
\alias{collect}
\title{as_tibble}
\usage{
\method{collect}{bcdc_promise}(x, ...)

\method{as_tibble}{bcdc_promise}(x, ...)
}
\arguments{
\item{x}{object of class \code{bcdc_promise}}
}
\description{
See \code{tibble::\link[tibble]{as_tibble}} for details.

After tuning a query, \code{collect()} is used to actually bring the data into memory.
This will retrieve an sf object into R. The \code{as_tibble()} function can be used
interchangeably with \code{collect} which matches \code{dbplyr} behaviour.

See \code{dplyr::\link[dplyr:compute]{collect}} for details.
}
\examples{
\donttest{
try(
  bcdc_query_geodata("bc-airports") \%>\%
    collect()
)

try(
  bcdc_query_geodata("bc-airports") \%>\%
    as_tibble()
)
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc_search.R
\name{bcdc_search_facets}
\alias{bcdc_search_facets}
\title{Get the valid values for a facet (that you can use in \code{\link[=bcdc_search]{bcdc_search()}})}
\usage{
bcdc_search_facets(
  facet = c("license_id", "download_audience", "res_format", "sector", "organization",
    "groups")
)
}
\arguments{
\item{facet}{the facet(s) for which to retrieve valid values. Can be one or
more of:
\verb{"license_id", "download_audience", "res_format", "sector", "organization", "groups"}}
}
\value{
A data frame of values for the selected facet
}
\description{
Get the valid values for a facet (that you can use in \code{\link[=bcdc_search]{bcdc_search()}})
}
\examples{
\donttest{
try(
  bcdc_search_facets("res_format")
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc_browse.R
\name{bcdc_browse}
\alias{bcdc_browse}
\title{Load the B.C. Data Catalogue URL into an HTML browser}
\usage{
bcdc_browse(
  query = NULL,
  browser = getOption("browser"),
  encodeIfNeeded = FALSE
)
}
\arguments{
\item{query}{Default (NULL) opens a browser to \code{https://catalogue.data.gov.bc.ca}.
This argument will also accept a B.C. Data Catalogue record ID or name to take you
directly to that page. If the provided ID or name doesn't lead to a valid webpage,
bcdc_browse will search the data catalogue for that string.}

\item{browser}{a non-empty character string giving the name of the
    program to be used as the HTML browser.  It should be in the PATH,
    or a full path specified.  Alternatively, an \R function to be
    called to invoke the browser.

    Under Windows \code{NULL} is also allowed (and is the default), and
    implies that the file association mechanism will be used.
  }

\item{encodeIfNeeded}{Should the URL be encoded by
    \code{\link[utils]{URLencode}} before passing to the browser?  This is not
    needed (and might be harmful) if the \code{browser} program/function
    itself does encoding, and can be harmful for \samp{file://} URLs on some
    systems and for \samp{http://} URLs passed to some CGI applications.
    Fortunately, most URLs do not need encoding.}
}
\value{
A browser is opened with the B.C. Data Catalogue URL loaded if the
session is interactive. The URL used is returned as a character string.
}
\description{
This is a wrapper around utils::browseURL with the URL for the B.C. Data Catalogue as
the default
}
\examples{
\donttest{
## Take me to the B.C. Data Catalogue home page
try(
  bcdc_browse()
)

## Take me to the B.C. airports catalogue record
try(
 bcdc_browse("bc-airports")
)

## Take me to the B.C. airports catalogue record
try(
  bcdc_browse("76b1b7a3-2112-4444-857a-afccf7b20da8")
)
}
}
\seealso{
\code{\link[utils]{browseURL}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc-get-citation.R
\name{bcdc_get_citation}
\alias{bcdc_get_citation}
\title{Generate a bibentry from a Data Catalogue Record}
\usage{
bcdc_get_citation(record)
}
\arguments{
\item{record}{either a \code{bcdc_record} object (from the result of \code{bcdc_get_record()}),
a character string denoting the name or ID of a resource (or the URL)

It is advised to use the permanent ID for a record rather than the
human-readable name to guard against future name changes of the record.
If you use the human-readable name a warning will be issued once per
session. You can silence these warnings altogether by setting an option:
\code{options("silence_named_get_data_warning" = TRUE)} - which you can set
in your .Rprofile file so the option persists across sessions.}
}
\description{
Generate a "TechReport" bibentry object directly from a catalogue record.
The primary use of this function is as a helper to create a \code{.bib} file for use
in reference management software to cite data from the B.C. Data Catalogue.
This function is likely to be starting place for this process and manual
adjustment will often be needed. The bibentries are not designed to be
authoritative and may not reflect all fields required for individual
citation requirements.
}
\examples{

}
\seealso{
\code{\link[utils:bibentry]{utils::bibentry()}}

try(
bcdc_get_citation("76b1b7a3-2112-4444-857a-afccf7b20da8")
)
\subsection{Or directly on a record object}{

try(
rec <- bcdc_get_record("76b1b7a3-2112-4444-857a-afccf7b20da8")
bcdc_get_citation(rec)
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cql-geom-predicates.R
\name{bcdc_check_geom_size}
\alias{bcdc_check_geom_size}
\title{Check spatial objects for WFS spatial operations}
\usage{
bcdc_check_geom_size(x)
}
\arguments{
\item{x}{object of class sf, sfc or sfg}
}
\value{
invisibly return logical indicating whether the check pass. If the return
value is TRUE, the object will not need a bounding box drawn. If the return value is
FALSE, the check will fails and a bounding box will be drawn.
}
\description{
Check a spatial object to see if it exceeds the current set value of
'bcdata.max_geom_pred_size' option, which controls how the object is treated when used inside a spatial predicate function in \code{\link[=filter.bcdc_promise]{filter.bcdc_promise()}}. If the object does exceed the size
threshold a bounding box is drawn around it and all features
within the box will be returned. Further options include:
\itemize{
\item Try adjusting the value of the 'bcdata.max_geom_pred_size' option
\item Simplify the spatial object to reduce its size
\item Further processing on the returned object
}
}
\details{
See the \href{https://bcgov.github.io/bcdata/articles/efficiently-query-spatial-data-in-the-bc-data-catalogue.html}{Querying Spatial Data with bcdata}
for more details.
}
\examples{
\donttest{
try({
  airports <- bcdc_query_geodata("bc-airports") \%>\% collect()
  bcdc_check_geom_size(airports)
})
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdc-web-services.R
\name{bcdc_query_geodata}
\alias{bcdc_query_geodata}
\title{Query data from the B.C. Web Feature Service}
\usage{
bcdc_query_geodata(record, crs = 3005)
}
\arguments{
\item{record}{either a \code{bcdc_record} object (from the result of \code{bcdc_get_record()}),
a character string denoting the name or ID of a resource (or the URL) or a BC Geographic
Warehouse (BCGW) name.

It is advised to use the permanent ID for a record or the BCGW name rather than the
human-readable name to guard against future name changes of the record.
If you use the human-readable name a warning will be issued once per
session. You can silence these warnings altogether by setting an option:
\code{options("silence_named_get_data_warning" = TRUE)} - which you can set
in your .Rprofile file so the option persists across sessions.}

\item{crs}{the epsg code for the coordinate reference system. Defaults to
\code{3005} (B.C. Albers). See https://epsg.io.}
}
\value{
A \code{bcdc_promise} object. This object includes all of the information
required to retrieve the requested data. In order to get the actual data as
an \code{sf} object, you need to run \code{\link[=collect]{collect()}} on the \code{bcdc_promise}.
}
\description{
Queries features from the B.C. Web Feature Service. See
\code{\link[=bcdc_tidy_resources]{bcdc_tidy_resources()}} - if a resource has a value of
\code{"wms"} in the \code{format} column it is available as a Web
Feature Service, and you can query and download it
using \code{bcdc_query_geodata()}. The response will be
paginated if the number of features is above the number
set by the \code{bcdata.single_download_limit} option.
Please see \code{\link[=bcdc_options]{bcdc_options()}} for defaults and more
information.
}
\details{
Note that this function doesn't actually return the data, but rather an
object of class \code{bcdc_promise}, which includes all of the information
required to retrieve the requested data. In order to get the actual data as
an \code{sf} object, you need to run \code{\link[=collect]{collect()}} on the \code{bcdc_promise}. This
allows further refining the call to \code{bcdc_query_geodata()} with \code{\link[=filter]{filter()}}
and/or \code{\link[=select]{select()}} statements before pulling down the actual data as an \code{sf}
object with \code{\link[=collect]{collect()}}. See examples.
}
\examples{

\donttest{
# Returns a bcdc_promise, which can be further refined using filter/select:
try(
  bcdc_query_geodata("bc-airports", crs = 3857)
)

# To obtain the actual data as an sf object, collect() must be called:
try(
  bcdc_query_geodata("bc-airports", crs = 3857) \%>\%
    filter(PHYSICAL_ADDRESS == 'Victoria, BC') \%>\%
    collect()
)

try(
  bcdc_query_geodata("groundwater-wells") \%>\%
    filter(OBSERVATION_WELL_NUMBER == "108") \%>\%
    select(WELL_TAG_NUMBER, INTENDED_WATER_USE) \%>\%
    collect()
)

## A moderately large layer
try(
  bcdc_query_geodata("bc-environmental-monitoring-locations")
)

try(
  bcdc_query_geodata("bc-environmental-monitoring-locations") \%>\%
    filter(PERMIT_RELATIONSHIP == "DISCHARGE")
)


## A very large layer
try(
  bcdc_query_geodata("terrestrial-protected-areas-representation-by-biogeoclimatic-unit")
)

## Using a BCGW name
try(
  bcdc_query_geodata("WHSE_IMAGERY_AND_BASE_MAPS.GSR_AIRPORTS_SVW")
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bcdata-package.R
\docType{package}
\name{bcdata-package}
\alias{bcdata}
\alias{bcdata-package}
\title{bcdata: Search and Retrieve Data from the BC Data Catalogue}
\description{
\if{html}{\figure{logo.png}{options: align='right' alt='logo' width='120'}}

Search, query, and download tabular and 'geospatial' data from the British Columbia Data Catalogue (<https://catalogue.data.gov.bc.ca/>). Search catalogue data records based on keywords, data licence, sector, data format, and B.C. government organization. View metadata directly in R, download many data formats, and query 'geospatial' data available via the B.C. government Web Feature Service ('WFS') using 'dplyr' syntax.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://bcgov.github.io/bcdata/}
  \item \url{https://catalogue.data.gov.bc.ca/}
  \item \url{https://github.com/bcgov/bcdata/}
  \item Report bugs at \url{https://github.com/bcgov/bcdata/issues}
}

}
\author{
\strong{Maintainer}: Andy Teucher \email{andy.teucher@gov.bc.ca} (\href{https://orcid.org/0000-0002-7840-692X}{ORCID})

Authors:
\itemize{
  \item Sam Albers \email{sam.albers@gov.bc.ca} (\href{https://orcid.org/0000-0002-9270-7884}{ORCID}) [contributor]
  \item Stephanie Hazlitt \email{stephanie.hazlitt@gov.bc.ca} (\href{https://orcid.org/0000-0002-3161-2304}{ORCID}) [contributor]
}

Other contributors:
\itemize{
  \item Province of British Columbia [copyright holder]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-classes.R, R/utils-filter.R
\name{filter.bcdc_promise}
\alias{filter.bcdc_promise}
\alias{filter}
\title{Filter a query from bcdc_query_geodata()}
\usage{
\method{filter}{bcdc_promise}(.data, ...)
}
\arguments{
\item{.data}{object of class \code{bcdc_promise} (likely passed from \code{\link[=bcdc_query_geodata]{bcdc_query_geodata()}})}

\item{...}{Logical predicates with which to filter the results. Multiple
conditions are combined with \code{&}. Only rows where the condition evaluates to
\code{TRUE} are kept. Accepts normal R expressions as well as any of the special
\link[=cql_geom_predicates]{CQL geometry functions} such as \code{WITHIN()} or \code{INTERSECTS()}.
If you know \code{CQL} and want to write a \code{CQL} query directly, write it enclosed
in quotes, wrapped in the \code{\link[=CQL]{CQL()}} function. e.g., \code{CQL("ID = '42'")}}
}
\description{
Filter a query from Web Feature Service using dplyr
methods. This filtering is accomplished lazily so that
the full sf object is not read into memory until
\code{collect()} has been called.

See \code{dplyr::\link[dplyr]{filter}} for details.
}
\section{Methods (by class)}{
\itemize{
\item \code{bcdc_promise}: filter.bcdc_promise
}}

\examples{
\donttest{
try(
  crd <- bcdc_query_geodata("regional-districts-legally-defined-administrative-areas-of-bc") \%>\%
    filter(ADMIN_AREA_NAME == "Cariboo Regional District") \%>\%
    collect()
)

try(
  ret1 <- bcdc_query_geodata("fire-perimeters-historical") \%>\%
    filter(FIRE_YEAR == 2000, FIRE_CAUSE == "Person", INTERSECTS(crd)) \%>\%
    collect()
)
  }
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_data.R
\name{bcdc_read_functions}
\alias{bcdc_read_functions}
\title{Formats supported and loading functions}
\usage{
bcdc_read_functions()
}
\description{
Provides a tibble of formats supported by bcdata and the associated function that
reads that data into R. This function is meant as a resource to determine which parameters
can be passed through the \code{bcdc_get_data} function to the reading function. This is
particularly important to know if the data requires using arguments from the read in function.
}
