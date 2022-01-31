
<!-- README.md is generated from README.Rmd. Please edit that file -->

# `ruODK`: An R Client for the ODK Central API <img src="man/figures/ruODK2.png" align="right" alt="Especially in these trying times, it is important to ask: ruODK?" width="120" />

<!-- badges: start -->

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3953158.svg)](https://doi.org/10.5281/zenodo.3953158)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Last-changedate](https://img.shields.io/github/last-commit/ropensci/ruODK.svg)](https://github.com/ropensci/ruODK/commits/main)
[![GitHub
issues](https://img.shields.io/github/issues/ropensci/ruodk.svg?style=popout)](https://github.com/ropensci/ruODK/issues/)
[![CI - GitHub
Actions](https://github.com/ropensci/ruODK/workflows/tic/badge.svg)](https://github.com/ropensci/ruODK/actions)
[![CI -
Appveyor](https://ci.appveyor.com/api/projects/status/1cs19xx0t64bmd2q/branch/master?svg=true)](https://ci.appveyor.com/project/florianm/ruodk/branch/main)
[![Test
coverage](https://codecov.io/gh/ropensci/ruODK/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/ruODK)
[![CodeFactor](https://www.codefactor.io/repository/github/ropensci/ruodk/badge)](https://www.codefactor.io/repository/github/ropensci/ruodk)
[![Hosted RStudio with
ruODK](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dbca-wa/urODK/master?urlpath=rstudio)
<!-- badges: end -->

Especially in these trying times, it is important to ask “r u ODK?”.

`ruODK` is an R client to access and parse data from ODK Central.

[OpenDataKit](https://getodk.org/) (ODK) is [free-and open-source
software](https://getodk.org/software/) that helps millions of people
collect data quickly, accurately, offline, and at scale. The software is
in active use in every country in the world and is supported by a large
and helpful community.

`ruODK` is a community contribution to the ODK ecosystem, but not
directly affiliated with ODK.

`ruODK` assumes some familiarity of its users with the ODK ecosystem and
workflows. For a detailed overview, read the extensive [ODK
documentation](https://docs.getodk.org/) and visit the friendly [ODK
forum](https://forum.getodk.org/).

[ODK Central](https://docs.getodk.org/central-intro/) is a cloud-based
data clearinghouse for digitally captured data, replacing the older
software [ODK Aggregate](https://docs.getodk.org/aggregate-intro/). ODK
Central manages user accounts and permissions, stores form definitions,
and allows data collection clients like [ODK
Collect](https://docs.getodk.org/collect-intro/) to connect to it for
form download and submission upload.

![An ODK setup with ODK Build, Central, Collect, and
ruODK](https://www.lucidchart.com/publicSegments/view/952c1350-3003-48c1-a2c8-94bad74cdb46/image.png)

A typical [ODK workflow](https://docs.getodk.org/#how-is-odk-used): An
XForm is designed e.g. in [ODK Build](https://build.getodk.org/),
[published to ODK Central](https://docs.getodk.org/central-forms/), and
downloaded onto an Android device running ODK Collect. After data have
been captured digitally using [ODK
Collect](https://docs.getodk.org/collect-intro/), the data are uploaded
and stored in ODK Central. The next step from there is to extract the
data, optionally upload it into another data warehouse, and then to
analyse and generate insight from it.

While data can be retrieved in bulk through the GUI, ODK Central’s API
provides access to its data and functionality through both an OData and
a RESTful API with a comprehensive and interactive
[documentation](https://odkcentral.docs.apiary.io/#reference/odata-endpoints).

`ruODK` is aimed at the technically minded researcher who wishes to
access and process data from ODK Central using the programming language
R.

Benefits of using the R ecosystem in combination with ODK:

-   Scalability: Both R and ODK are free and open source software.
    Scaling to many users does not incur license fees.
-   Ubiquity: R is known to many scientists and is widely taught at
    universities.
-   Automation: The entire data access and analysis workflow can be
    automated through R scripts.
-   Reproducible reporting (e.g. 
    [Sweave](https://support.rstudio.com/hc/en-us/articles/200552056-Using-Sweave-and-knitr),
    [RMarkdown](https://rmarkdown.rstudio.com/)), interactive web apps
    ([Shiny](https://shiny.rstudio.com/)), workflow scaling
    ([drake](https://docs.ropensci.org/drake/)).
-   Rstudio-as-a-Service (RaaS) at
    [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dbca-wa/urODK/master?urlpath=rstudio)

`ruODK`’s scope:

-   To wrap all ODK Central API endpoints with a focus on **data
    access**.
-   To provide working examples of interacting with the ODK Central API.
-   To provide convenience helpers for the day to day tasks when working
    with ODK Central data in R: **data munging** the ODK Central API
    output into tidy R formats.

<!-- TODO: vignette "workflows" -->

`ruODK`’s use cases:

-   Smaller projects: Example
    [rOzCBI](https://dbca-wa.github.io/rOzCBI/)
    1.  Data collection: ODK Collect
    2.  Data clearinghouse: ODK Central
    3.  Data analysis and reporting: `Rmd` (ruODK)
    4.  Publishing and dissemination:
        [`ckanr`](https://docs.ropensci.org/ckanr/),
        [`CKAN`](https://ckan.org/)
-   Larger projects:
    1.  Data collection: ODK Collect
    2.  Data clearinghouse: ODK Central
    3.  ETL pipeline into data warehouses: `Rmd` (ruODK)
    4.  QA: in data warehouse
    5.  Reporting: `Rmd`
    6.  Publishing and dissemination:
        [`ckanr`](https://docs.ropensci.org/ckanr/),
        [`CKAN`](https://ckan.org/)

Out of scope:

-   To wrap “management” API endpoints. ODK Central is a [VueJS/NodeJS
    application](https://github.com/opendatakit/central-frontend/) which
    provides a comprehensive graphical user interface for the management
    of users, roles, permissions, projects, and forms.
-   To provide extensive data visualisation. We show only minimal
    examples of data visualisation and presentation, mainly to
    illustrate the example data. Once the data is in your hands as tidy
    tibbles… urODK!

## A quick preview

<img src="man/figures/odata.svg" alt="ruODK screencast" width="100%" />

## Install

You can install the latest release of `ruODK` from the [rOpenSci
R-Universe](https://ropensci.r-universe.dev):

``` r
# Enable the rOpenSci universe
options(repos = c(ropensci = 'https://ropensci.r-universe.dev',
                  CRAN = 'https://cloud.r-project.org'))
install.packages('ruODK')
```

Alternatively, you can install the development version from the `main`
branch.

``` r
if (!requireNamespace("remotes")) install.packages("remotes")
# Full install
remotes::install_github(
   "ropensci/ruODK@main", 
   dependencies = TRUE, 
   upgrade = "always",
   build_vignettes = TRUE)

# Minimal install without vignettes
remotes::install_github(
   "ropensci/ruODK@main", 
   dependencies = TRUE, 
   upgrade = "ask",
   build_vignettes = FALSE)
```

If the install fails, read the error messages carefully and install any
unmet dependencies (system libraries or R packages).

If the install fails on building the vignettes, you can set
`build_vignettes=FALSE` and read the vignettes from the online docs
instead.

If the installation still fails, or the above does not make any sense,
feel free to submit a [bug
report](https://github.com/ropensci/ruODK/issues/new/choose).

## ODK Central

### Access to an ODK Central instance

First, we need an ODK Central instance and some data to play with!

Either [request a free trial](https://getodk.org/#odk-cloud) or follow
the [setup instructions](https://docs.getodk.org/central-intro/) to
build and deploy your very own ODK Central instance.

### ODK Central setup

The ODK Central [user manual](https://docs.getodk.org/central-using/)
provides up-to-date descriptions of the steps below.

-   [Create a web user
    account](https://docs.getodk.org/central-users/#creating-a-web-user)
    on an ODK Central instance. Your username will be an email address.
-   [Create a project](https://docs.getodk.org/central-projects/) and
    give the web user at least [read
    permissions](https://docs.getodk.org/central-projects/#managing-project-managers).
-   Create an XForm, e.g. using ODK Build, or use the [example
    forms](https://github.com/ropensci/ruODK/tree/master/inst/extdata)
    provided by `ruODK`. The `.odkbuild` versions can be loaded into
    [ODK Build](https://build.getodk.org/), while the `.xml` versions
    can be directly imported into ODK Central.
-   [Publish the form](https://docs.getodk.org/central-forms/) to ODK
    Central.
-   Collect some data for this form on ODK Collect and let ODK Collect
    submit the finalised forms to ODK Central.

## Configure `ruODK`

For all available detailed options to configure authentication for
`ruODK`, read
[`vignette("setup", package = "ruODK")`](https://docs.ropensci.org/ruODK/articles/setup.html).

## Use ruODK

A detailed walk-through with some data visualisation examples is
available in the
[`vignette("odata-api", package="ruODK")`](https://docs.ropensci.org/ruODK/articles/odata-api.html).

See also
[`vignette("restful-api", package="ruODK")`](https://docs.ropensci.org/ruODK/articles/restful-api.html)
for examples using the alternative RESTful API.

## Try ruODK

[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dbca-wa/urODK/master?urlpath=rstudio)
will launch a disposable, hosted RStudio instance with `ruODK` installed
and the companion package [`urODK`](https://github.com/dbca-wa/urODK)
opened as starting point for a hands-on workshop or instant demo of
`ruODK` usage.

Create a new RMarkdown workbook from `ruODK` template “ODK Central via
OData” and follow the instructions within.

## Contribute

Contributions through [issues](https://github.com/ropensci/ruODK/issues)
and PRs are welcome!

See the [contributing
guide](https://docs.ropensci.org/ruODK/CONTRIBUTING.html) on best
practices and further readings for code contributions.

## Attribution

`ruODK` was developed, and is maintained, by Florian Mayer for the
Western Australian [Department of Biodiversity, Conservation and
Attractions (DBCA)](https://www.dbca.wa.gov.au/). The development was
funded both by DBCA core funding and external funds from the [North West
Shelf Flatback Turtle Conservation
Program](https://flatbacks.dbca.wa.gov.au/).

To cite package `ruODK` in publications use:

``` r
citation("ruODK")
#> 
#> To cite ruODK in publications use (with the respective version number:
#> 
#>   Mayer, Florian Wendelin. (2020, Nov 19).  ruODK: An R Client for the
#>   ODK Central API (Version X.X.X).  Zenodo.
#>   https://doi.org/10.5281/zenodo.3953158
#> 
#> A BibTeX entry for LaTeX users is
#> 
#>   @Misc{,
#>     title = {ruODK: Client for the ODK Central API},
#>     author = {Florian W. Mayer},
#>     note = {R package version X.X.X},
#>     year = {2020},
#>     url = {https://github.com/ropensci/ruODK},
#>   }
```

## Acknowledgements

The Department of Biodiversity, Conservation and Attractions (DBCA)
acknowledges the traditional owners of country throughout Western
Australia and their continuing connection to the land, waters and
community. We pay our respects to them, their culture and to their
Elders past and present.

This software was created on Whadjuk boodja (ground) both as a
contribution to the ODK ecosystem and for the conservation of the
biodiversity of Western Australia, and in doing so, caring for country.

## Package functionality

See
[`vignette("comparison", package="ruODK")`](https://docs.ropensci.org/ruODK/articles/comparison.html)
for a comprehensive comparison of ruODK to other software packages from
both an ODK and an OData angle.
# `ruODK` 1.3.0
This release supports ODK Central 1.3.0 and represents an over-due version
bump to reflect the supported ODK Central version.
The test server is now updated to ODK Central 1.3.0, and all tests pass.

There are still some newer and as yet unsupported API endpoints in ODK Central, 
which serve administrative purposes of the front-end. Contributions are welcome,
get started on [these issues](https://github.com/ropensci/ruODK/milestone/2)
and the contributing guide. As ruODK focuses on data retrieval, these
administrative endpoints are non-critical to ruODK's purpose.

## Major fixes
## Minor fixes
## Documentation
## Data
* Packaged data has been re-created to represent the latest server outputs.
## Maintenance
* All tests pass, GitHub Actions is as per usual brittle at the installation 
  step.
  

# `ruODK` 1.2.0.0000
We are shaping up to a release targetting the ODK Central 1.2 release.
ODK Central is undergoing some bug fixes and patches, while ruODK's test server
will be migrated to another instance. The latter is required to enable tests
which create/update/delete entities in ODK Central.

## Major fixes
## Minor fixes
## Documentation
## Data
## Maintenance
* All DEPENDS and SUGGESTS bumped to latest available under R release.
* The minimum supported version is `rversions::r_oldrel()`: 4.0.5 (2021-03-31).
* ruODK is developed under `rversions::r_release()`: 4.1.1 (2021-08-10).

# `ruODK` 1.2.0
## Major fixes
* ODK Central returns Geoshapes as Multipolygon. `split_geoshape` adjusted for
  `odkce_version` >= 1.2. (#131)
* `readr::parse_datetime` stopped supporting timezone "Z". (#130)
## Minor fixes
## Documentation
## Data
* All data refreshed from test server running ODK Central 1.2.
## Maintenance


# `ruODK` 0.10.2
## Major fixes
## Minor fixes
## Documentation
* Update installation instructions with source install from rOpenSci R-Universe
  and troubleshooting (#128, thanks @yanokwa)
## Data
## Maintenance

# `ruODK` 0.10.2
## Major fixes
* Fix ODK Central v1.2 time out on NULL query parameters `skip` and `top`. 
  ruODK now only supplies non-NULL query parameters and has an additional
  seat-belt to drop any query parameter that is an empty string.
  (#126, thanks @yanokwa, @mtyszler, @thaliehln)
## Minor fixes
## Documentation
## Data
## Maintenance

# `ruODK` 0.10.1
## Major fixes
* `submission_export` now terminates immediately if an incorrect passphrase is 
  given. ODK Central can exceed memory limits if `submission_export` is run 
  repeatedly with an incorrect passphrase. (#30, thanks @Thaliehln)
## Minor fixes
* Add `encryption_key_list`
## Documentation
## Data
## Maintenance

# `ruODK` 0.10.0
## Major fixes
## Minor fixes
* Make `ru_msg_*` conditional to `get_ru_verbose()`.
## Documentation
* Reference re-ordered into topics.
* Re-worded the example preamble on setup.
* Added metadata to pkgdown config.
* Enable Markdown docs. Embedded lifecycle badges should work. R CMD Check does
  not complain which is nice.
## Data
## Maintenance
* Prepare to cover remaining API endpoints.

# `ruODK` 0.9.11
## Major fixes
## Minor fixes
* Add `published_at` to `form_list` and `form_detail`, update examples, tests, 
  test fixtures to show that draft forms can be detected by a NA `published_at`
  in ODK Central versions having form drafts, and by NA `hash` and `version`
  in ODK Central versions before form drafts.
## Documentation
## Data
## Maintenance

# `ruODK` 0.9.10
This is a "everything so far works" release. 
There are a few ODK Central API endpoints waiting to be implemented still.

## Major fixes
## Minor fixes
* Default ODK Central version bumped to current release (1.1)
## Documentation
* PDF manual updated
* Welcoming new contributors Marcelo (@mtyszler) and Hélène (@Thaliehln)
## Data
## Maintenance
* Updated Zenodo archive at <https://zenodo.org/record/4609910>
* Updated Docker image at <https://hub.docker.com/u/dbcawa/ruodk> 
  (RStudio server with ruODK)

# `ruODK` 0.9.9
## Major fixes
* `submission_export` now supports ODK Central v1.1 features to omit media
  attachments (`media = FALSE`), and to omit repeat data (`include_repeats=FALSE`).
  Calling `submission_export` on an older version of ODK Central 
  (as determined through `get_default_odkc_version()`) with these new parameters
  will emit a (verbose only) "noop" message and not act further on them.
## Minor fixes
* Bugfix to `submission_export` on encrypted forms with multiple encryption 
  keys. (Thanks @Thaliehln #117 #30)

# `ruODK` 0.9.8
## Minor fixes
* Add support for passphrase to the `ru_setup` family (#116)


# `ruODK` 0.9.7
## Major fixes
* `odata_submission_get()` bugfix: `handle_ru_attachments()` 
  now finds and downloads media attachments from both main submissions and 
  nested form groups. (#114)
* `odata_submission_get()` bugfix: missing media attachments (due to upload 
  error from ODK Collect to ODK Central) are tolerated without interrupting the
  batch download. A diagnostic warning message will be emitted for each failed
  download. (#114)
## Minor fixes
## Documentation
## Data
## Maintenance

# `ruODK` 0.9.6
## Major fixes
* Support encryption (#30 #110, @Thaliehln).
  * Note that `ruODK` only supports one passphrase. When switching between
    multiple encrypted forms, it would make sense to store the different 
    passphrases in separate environment variables, and refer to these environment
    variables explicitly in function calls.
  * The updated ruODK::submission_export should now export data 
    from both encrypted projects and non-encrypted projects.
    The HTTP method is changed from GET to POST and encryption key ID / 
    passphrase are provided via POST body using a JSON format. 
    Encrypted forms can be extracted and inspected like non-encrypted forms:
    
```{r, eval=FALSE}
se <- submission_export()
t <- tempdir()
f <- unzip(se, exdir = t)
fs::dir_ls(t)
fid <- get_test_fid()
sub <- fs::path(t, glue::glue("Locations.csv")) %>% readr::read_csv()
sub %>% knitr::kable(.)
```

# `ruODK` 0.9.5
## Major fixes
* `form_schema_ext` retrieves choice lists when choice filters are present
   (#105,  @mtyszler).
## Minor fixes
## Documentation
## Data
## Maintenance

# `ruODK` 0.9.4
## Major fixes
* `form_schema_ext` performance enhancement (#106, thanks @mtyszler).
## Maintenance
* Tests use vcr to cache server response (#104).
  Delete the local cache `tests/fixtures` to re-generate the vcr cache, or 
  enjoy much faster running tests using cached server response.

# `ruODK` 0.9.3
This is a point release to create a new RStudio Server image.
## Minor fixes
* Form schema now also works on draft forms (#103, HT @dmenne).
## Maintenance
* Automated code reviews by <codefactor.io>.
* GitHub Actions welcomes Ubuntu 20.04 LTS and MacOS X back to the passing 
  environments.

# `ruODK` 0.9.2
## Major fixes
* `form_schema_ext()` Shows the extended schema of one form, including 
  (multi-language) labels and choice lists. (#77, thanks @mtyszler for the PR)
* Development continues in the default branch `main`.
## Minor fixes
* `form_list` now handles draft forms with `NA` hash and version (#86, 
  thanks @dmenne for the PR).
* Migrate package tests to a new ODK Central instance and update contributing 
  guidelines with new settings. (#14)
* Drop Import of `tidyselect` in favour of using `dplyr::all_of()`.
* All calls to `httr::RETRY(times=)` are configurable via setting `retries`. (#94)
## Documentation
* Mapview popups in vignette "Spatial" are back after an upstream
  bugfix is in progress (https://github.com/r-spatial/mapview/issues/312).

# `ruODK` 0.9.1
## Major fixes
ODK Central versions 0.7 to 0.9 export geotraces and geoshapes via OData with
a trailing empty coordinate. `ruODK` removes any trailing empty coordinates from
both GeoJSON and WKT formats. (#88, HT @TimonWeitkamp for the bug report)

## Documentation
A new vignette "Spatial" demonstrates how to parse spatial data into native 
formats, such as `sf`, and gives pointers on what to do next with them.

 
# `ruODK` 0.9.0
This is the release on passing 
[rOpenSci peer review](https://github.com/ropensci/software-review/issues/335).

Thanks to the rOpenSci editors and reviewers @maelle, @karissawhiting and 
@jmt2080ad, as well as to @OdiljonQurbon, @dickoa, @arestrom and @dmenne.

A DOI was minted at <https://doi.org/10.5281/zenodo.3953159>.

# `ruODK` 0.6.6
This version addresses ROpenSci reviewer comments from @karissawhiting and 
@jmt2080ad, contributions from @dickoa, as well as ideas and suggestions from 
@OdiljonQurbon, @arestrom and @dmenne.

This version supports ODK Central 0.9 while providing backwards compatibility
for ODK Central <= 0.7.

## Major fixes
* New environment variables `ODKC_(TEST_)VERSION` allow `ruODK` to toggle
  between deprecated/removed and new/added API endpoints, e.g. `form_schema`. (#61)
* Split and rename WKT POINT (ODK geopoint) fields with 
  `odata_submission_get(wkt=T)`. (#31 #7 HT @OdiljonQurbon)
* `submission_get` now accepts a vector of (all or selected) submission instance 
  IDs (`iid`), similar to `odata_submission_get()`. (#38)
* All `httr::GET()` are now replaced with `httr::RETRY("GET", ...)` (#48)
* Refactor `odata_submission_parse()` into `odata_submission_rectangle()`,
  `handle_ru_{geopoints, geotraces, geoshapes, datetimes, attachments}`. (#54 #69)
* Maintain backwards compatibility to ODK Central v7 for changed spatial output
  for geotrace and geoshape (#70)

## Minor fixes
* Drop `. <- NULL` in favour of `utils::globalVariables(".")`. (#35)
* Print settings now hides passwords but instructs how to show them. (#37)
* `ru_setup()` now prints settings. (#37)
* `parse_datetime()` renamed to `ru_datetime()` to avoid naming conflict with 
  `readr::parse_datetime()`. (#43)
* Add a global default for verbosity. (#51 HT @arestrom)
* Add a global default for time zone. (#53 HT @arestrom)
* Use `httr::modify_url` to build URLs rather than `glue::glue` (#66)
* Silenced spurious messages from `tibble::as_tibble()` which is called from 
  `odata_submission_rectangle()`. Use `ru_verbose` to toggle useful diagnostic
  messages. (#79 HT @dmenne)
* Renamed `master` branch to `main`, updated docs (HT @arestrom #81)

## Dependencies
* Moved `rlist` to Imports, as it is now used in `odata_submission_get()`. (#6)
* Dropped `rlist` dependency in favour of `tidyr`/`dplyr`. (#85)
* Use development versions of `rlang`, `purrr` and `lifecycle` to get their
  latest bug fixes.

## Data
* Use canned data in all vignettes, so they can build without authentication. 
  (#33)
* Update canned data (and `make-data.R`) to work with CI-built pkgdown site.
  (#52)
* Remove nested list `form_schema_raw` which is only generated from legacy 
  ODK Central (< 0.8) (#61)
* Added ODK Central < v0.7 form schema for tests.

## Documentation
* Updated workshop companion package [urODK](https://github.com/dbca-wa/urODK).
* Rename vignettes to `odata-api` and `restful-api`. (#34)
* Warn against using plain text credentials in vignette `setup`. (#34)
* More documentation improvements at 
  [#34](https://github.com/ropensci/ruODK/issues/34).
* Add screencast to the README. HT to asciicast! (#45)
* Improve logo - more turtles, but questionable photoshopping.
* Add examples where missing. (#32)
* Build pkgdown site via GH actions. (#52)
* Minor typographic changes: end every function title with a full stop.
* Broken links and other inconsistencies fixed after contributions from the 
  ODK Forum, @dickoa, @arestrom, @dmenne. 
  Thanks for the first community feedback! (#76 #74 #79 #80 #81)

# Docker
* Added Dockerfile to build an RStudio Server image based on `rocker/geospatial`
  with `ruODK` and dependencies installed.
* Build a separate Dockerfile for [`urODK`](https://github.com/dbca-wa/urODK)
  to launch a hosted RStudio instance in Binderhub. (#83)

# `ruODK` 0.6.6
* The big one has landed: `odata_submission_get()` now defaults to parse 
  submissions from nested lists into a tibble, parse dates and datetimes,
  download and link attachments. (#6)

# `ruODK` 0.6.5
## Documentation
* Use lifecycle badges on functions. Add lifecycle to dependencies, version bump
  `usethis`. (#29)
  
## Code
* Refactor list wrangling code to use `map_*(.default=NA)`, removing some 
  internal helpers (thanks to @jennybc)
* Use dummy imports to silence R CMD check NOTES as per [googledrive](https://github.com/tidyverse/googledrive/blob/050a982cba630503702bdde05a77d727baa36d48/R/googledrive-package.R)'s 
  example (thanks to @jennybc)
* Drop unused internal helper functions

# `ruODK` 0.6.4
## Data
*  Use three new test forms to make package smaller and tests faster.
   Use the main test form for example data, including data from two nested tables.
   Use the main test form in all vignettes and README.
   Use a small form without attachments for tests repeatedly exporting to ZIP.
   Use another small form with only one submission and two attachments for tests
   downloading attachments.
   The test credentials are unchanged (#23)

# `ruODK` 0.6.3
## Dependencies
*  [tidyr 1.0.0](https://www.tidyverse.org/articles/2019/09/tidyr-1-0-0/) is out!
  Move `{tidyr}` dependency from GitHub master to CRAN version (#27)
*  Add `{usethis}` to Suggests, it is used in the setup step
  
## Documentation
*  Add [David Henry](https://github.com/schemetrica)'s 
  [Pentaho Kettle tutorial](https://forum.opendatakit.org/t/automating-data-delivery-using-the-odata-endpoint-in-odk-central/22010) 
  to the software review in the README (#28)
  
## DIY for workshops
*  Add inaugural RMarkdown template "odata" (#26)
*  Add [binder](https://mybinder.org/) launch button (one click start for #26)

# `ruODK` 0.6.2
## Code
*  Simplifiy `ru_setup()` to use OData Service URL.
*  Change all functions to default to `get_default{pid,fid,url,un,pw}()`, partly
  moving project ID (pid) and form ID (fid) to kwargs. This changes all examples,
  tests, vignettes, READMEs.
*  Reduce installed package size by sharing attachment files. Add new parameter
  `separate=FALSE` to `attachment_get` to prevent separating attachment files 
  into subfolders named after their submission `uuid` (#22)
  
## Documentation
*  Add a high level overview diagram to README and `inst/joss/paper.md` to
  illustrate `ruODK`'s intended purpose in the ODK ecosystem (#19)
*  Added link to explain 
  [environment variables and R startup](https://whattheyforgot.org/r-startup.html) 
  to vignette "setup". @maelle
*  Add comparison of similar software to README (#25)

# `ruODK` 0.6.1
*  ROpenSci submission review [milestone](https://github.com/ropensci/ruODK/milestone/3), 
  [discussion](https://github.com/ropensci/software-review/issues/335).
*  Updates to documentation (#13 #15 #17)
*  Group function docs (#18)
*  Update contribution guidelines and add account request issue template:
  How to run `ruODK` tests and build the vignettes (#15 #20)
*  Add dedicated `ru_setup()` and `ru_settings()`. 
  Pat down functions for missing credentials and yell loudly but clearly about
  httr errors. (#16)
*  Drop `@importFrom` to reduce duplication. All external functions are prefixed
  with their package name already.
*  Add convenience helpers `attachment_link()` and `parse_datetime()`.
*  Use `janitor::clean_names()` on column names over home-grown helpers.
*  Since `submission_detail` is now metadata only, add `submission_get` to download
  submission data.

# `ruODK` 0.6.0
*  Version bump and lifecycle bump to indicate that `ruODK` is ready to be used
  against ODK Central 0.6.

# `ruODK` 0.6.0.9000
*  Version bump to match ODK Central's version.
*  Updated to match ODK Central's API at 0.6 (RESTful and OData) (c.f.).
*  Add commented out crosslinks to source code and related tests for convenience.
*  Encryption endpoints (new in 0.6) are not yet supported.
*  Audit logs supported, as they are a read-only export.

# `ruODK` 0.3.1
## Preparation for ROpenSci pre-submission
*  Check name with [`available::available("ruODK")`](https://devguide.ropensci.org/building.html#naming-your-package):
  *  Name valid: ✔
  *  Available on CRAN: ✔ 
  *  Available on Bioconductor: ✔
  *  Available on GitHub:  ✔ 
  *  Rude misinterpretations: none
  *  In summary: Package name seems to be OK. Well, ODK. OK, ruODK.
*  Added metadata via 
  [`codemetar::write_codemeta("ruODK")`](https://devguide.ropensci.org/building.html#creating-metadata-for-your-package).
*  [Cross-platform](https://devguide.ropensci.org/building.html#platforms): 
  Runs on GNU/Linux (TravisCI) and on Windows (AppVeyor)
*  [Function naming](https://devguide.ropensci.org/building.html#function-and-argument-naming)
  follows `object_verb`. 
  *  `ruODK` uses verb singulars (e.g. `submission_list` to 
  list multiple submissions), while ODK Central's API URLs use verb plurals.
  *  `ruODK` uses `snake_case`.
  *  Exception to `object_verb`: 
    Functions operating on the OData endpoints are named `odata_object_verb`.
    Helper functions not related to API endpoints are named `verb_object`. 
*  [Code style](https://devguide.ropensci.org/building.html#code-style) done
  by `styler::style_package()`, see section "Release" in `README.md`.
*  `ruODK` [has a `README.Rmd`](https://devguide.ropensci.org/building.html#readme) 
  and has a 
  [website generated by `pkgdown`](https://devguide.ropensci.org/building.html#website).
*  `ruODK` has documentation 
  [generated by `roxygen2`](https://devguide.ropensci.org/building.html#documentation).
*  [Test coverage](https://devguide.ropensci.org/building.html#testing) 100%, but
  could use more edge cases.
*  Tests use a live ODK Central instance, which is kept up to date by ODK.
  This means that ruODK's test always run against the latest master of ODK Central.
  `ruODK` does not (but maybe should?) use `webmockr` and `vcr`.
*  Spellchecks with `spelling::spell_check_package()`: added technical terms and
  function names to `inst/WORDLIST`.
*  Added citation and section in `README`.
*  Added `inst/joss/paper.md` for submission to JOSS.
*  Added [examples](https://devguide.ropensci.org/building.html#examples) to 
  function docs.

## TODO
*  Implement remaining missing functions ([ticket](https://github.com/ropensci/ruODK/issues/9)).

# `ruODK` 0.3.0
*  Use tidyverse issue template
*  Follow `goodpractice`
*  Created vignette "Setup"
*  Add AppVeyor
*  Refactor storage path of attachments to not contain "uuid:" (for Windows compat)
*  Started on REST API: `project_list`, `project_detail`, `form_list`, 
  `form_detail`. Naming scheme is `object_verb`. 
*  For now, functions related to the OData endpoint
  are named `verb_object`, maybe we should rename them to `odata_object_verb`.
*  Refactor URLs - build from project and form IDs
  
# `ruODK` 0.2.4
*  Cleaned up logo, thanks to @issa-tseng for suggestions

# `ruODK` 0.2.3
*  Added a new form, Flora Quadrat 0.3 to inst/extdata.

# `ruODK` 0.2.2
*  Added more test coverage.

# `ruODK` 0.2.1
*  Added various `usethis` goodies, test stubs, badges.

# `ruODK` 0.2.0
*  Recursively unnests raw data into wide format. Manual post-processing is still
  required to rename (anonymous in ODK and auto-named through R) coordinates,
  and to handle attachments.

# `ruODK` 0.1.0
*  Parses metadata, submissions, 
  and handles attachments (retaining already downloaded attachments).
*  Includes example forms as both `.xml` and `.odbbuild` in `inst/extdata`.
*  Includes example data as retrieved from ODK Central.
*  Includes vignette demonstrating tidying of retrieved data.
*  Handles repeating groups.

Roadmap:

*  Handle paginated submissions.
*  Retain already downloaded submissions.
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behaviour by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behaviour may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(<http://contributor-covenant.org>), version 1.0.0, available at 
<http://contributor-covenant.org/version/1/0/0/>
# Contributing 
This contributing guide has been derived from the `tidyverse` boilerplate.
Where it seems over the top, common sense is appreciated, and every contribution
is appreciated.

## Non-technical contributions to ruODK
Feel free to [report issues](https://github.com/ropensci/ruODK/issues):

* Bug reports are for unplanned malfunctions.
* Feature requests are for ideas and new features.
* Account requests are for getting access to the ODK Central instances run by DBCA
  (DBCA campaigns only) or the CI server (contributors, to run tests).

## Technical contributions to `ruODK`

If you would like to contribute to the code base, follow the process below.

*  [Prerequisites](#prerequisites)
*  [PR Process](#pr-process)
  *  [Fork, clone, branch](#fork-clone-branch)
  *  [Check](#check)
  *  [Style](#style)
  *  [Document](#document)
  *  [Test](#test)
  *  [NEWS](#news)
  *  [Re-check](#re-check)
  *  [Commit](#commit)
  *  [Push and pull](#push-and-pull)
  *  [Review, revise, repeat](#review-revise-repeat)
*   [Resources](#resources)
*   [Code of Conduct](#code-of-conduct)

This explains how to propose a change to `ruODK` via a pull request using
Git and GitHub. 

For more general info about contributing to `ruODK`, see the 
[Resources](#resources) at the end of this document.

### Prerequisites
To test the package, you will need valid credentials for the ODK Central instance 
used as a test server.
Create an [account request issue](https://github.com/ropensci/ruODK/issues/new/choose).

Before you do a pull request, you should always file an issue and make sure
the maintainers agree that it is a problem, and is happy with your basic proposal 
for fixing it. 
If you have found a bug, follow the issue template to create a minimal
[reprex](https://www.tidyverse.org/help/#reprex).

### Checklists
Some changes have intricate internal and external dependencies, which are easy
to miss and break. These checklists aim to avoid these pitfalls.

Test and update reverse dependencies (wastdr, urODK, etlTurtleNesting, etc.).

#### Adding a dependency
* Update DESCRIPTION
* Update GH Actions install workflows - do R package deps have system deps? Can GHA install them in all environments?
* Update Dockerfile
* Update urODK binder install.R
* Update installation instructions

#### Renaming a vignette
* Search-replace all links to the vignette throughout 
  * ruODK, 
  * urODK, 
  * ODK Central "OData" modal
  * ODK Central docs

#### Adding or updating a test form
* Update tests
* Update examples
* Update packaged data if test form submissions are included
* Add new cassette to vcr cache for each test using the test form

#### Adding or updating package data
* Update tests using the package data
* Update examples
* Update README if showing package data

#### Adding a settings variable
* Update ru_setup, ru_settings, update and add to settings tests
* Update .Renviron
* Update GitHub secrets
* Update tic.yml (add new env vars)
* Update vignette "Setup"

### PR process

#### Fork, clone, branch

The first thing you'll need to do is to [fork](https://help.github.com/articles/fork-a-repo/) 
the [`ruODK` GitHub repo](https://github.com/ropensci/ruODK), and 
then clone it locally. We recommend that you create a branch for each PR.

#### Check

Before changing anything, make sure the package still passes the below listed
flavours of `R CMD check` locally for you. 

```r
goodpractice::goodpractice(quiet = FALSE, )
devtools::check(cran = TRUE, remote = TRUE, incoming = TRUE)
chk <- rcmdcheck::rcmdcheck(args = c("--as-cran"))
```

#### Style

Match the existing code style. This means you should follow the tidyverse 
[style guide](http://style.tidyverse.org). Use the 
[styler](https://CRAN.R-project.org/package=styler) package to apply the style 
guide automatically.

Be careful to only make style changes to the code you are contributing. If you
find that there is a lot of code that doesn't meet the style guide, it would be
better to file an issue or a separate PR to fix that first.

```r
styler::style_pkg()
lintr:::addin_lint_package()
devtools::document(roclets = c("rd", "collate", "namespace"))
spelling::spell_check_package()
spelling::spell_check_files("README.Rmd", lang = "en_AU")
spelling::update_wordlist()
```

#### Document

We use [roxygen2](https://cran.r-project.org/package=roxygen2), specifically with the 
[Markdown syntax](https://cran.r-project.org/web/packages/roxygen2/vignettes/markdown.html),
to create `NAMESPACE` and all `.Rd` files. All edits to documentation
should be done in roxygen comments above the associated function or
object. Then, run `devtools::document()` to rebuild the `NAMESPACE` and `.Rd` 
files.

See the `RoxygenNote` in [DESCRIPTION](DESCRIPTION) for the version of
roxygen2 being used. 

```r
spelling::spell_check_package()
spelling::spell_check_files("README.Rmd", lang = "en_AU")
spelling::update_wordlist()
codemetar::write_codemeta("ruODK")
if (fs::file_info("README.md")$modification_time <
  fs::file_info("README.Rmd")$modification_time) {
  rmarkdown::render("README.Rmd", encoding = "UTF-8", clean = TRUE)
  if (fs::file_exists("README.html")) fs::file_delete("README.html")
}
```

#### Test

We use [testthat](https://cran.r-project.org/package=testthat). Contributions
with test cases are easier to review and verify. 

To run tests and build the vignettes, you'll need access to the 
[ruODK test server](https://odkc.dbca.wa.gov.au/).
If you haven't got an account yet, create an [accont request issue](https://github.com/ropensci/ruODK/issues/new/choose)
to request access to this ODK Central instance.

The tests require the following additions to your `.Renviron`:

```r
# Required for testing
ODKC_TEST_SVC="https://odkc.dbca.wa.gov.au/v1/projects/2/forms/Flora-Quadrat-04.svc"
ODKC_TEST_URL="https://odkc.dbca.wa.gov.au"
ODKC_TEST_PID=2
ODKC_TEST_PID_ENC=3
ODKC_TEST_PP="ThePassphrase"
ODKC_TEST_FID="Flora-Quadrat-04"
ODKC_TEST_FID_ZIP="Spotlighting-06"
ODKC_TEST_FID_ATT="Flora-Quadrat-04-att"
ODKC_TEST_FID_GAP="Flora-Quadrat-04-gap"
ODKC_TEST_FID_WKT="Locations"
ODKC_TEST_FID_I8N0="I8n_no_lang"
ODKC_TEST_FID_I8N1="I8n_label_lng"
ODKC_TEST_FID_I8N2="I8n_label_choices"
ODKC_TEST_FID_I8N3="I8n_no_lang_choicefilter"
ODKC_TEST_FID_I8N4="I8n_lang_choicefilter"
ODKC_TEST_FID_ENC="Locations"
ODKC_TEST_VERSION=1.0
RU_VERBOSE=TRUE
RU_TIMEZONE="Australia/Perth"
RU_RETRIES=3
ODKC_TEST_UN="..."
ODKC_TEST_PW="..."

# Your ruODK default settings for everyday use
ODKC_URL="..."
ODKC_PID=1
ODKC_FID="..."
ODKC_UN="..."
ODKC_PW="..."
```

Keep in mind that `ruODK` defaults to use `ODKC_{URL,UN,PW}`, so for everyday 
use outside of contributing, you will want to use your own `ODKC_{URL,UN,PW}`
account credentials.

```r
devtools::test()
devtools::test_coverage()
```

#### NEWS

For user-facing changes, add a bullet to `NEWS.md` that concisely describes
the change. Small tweaks to the documentation do not need a bullet. The format
should include your GitHub username, and links to relevant issue(s)/PR(s), as
seen below.

```md
* `function_name()` followed by brief description of change (#issue-num, @your-github-user-name)
```

#### Re-check

Before submitting your changes, make sure that the package either still
passes `R CMD check`, or that the warnings and/or notes have not _changed_
as a result of your edits.

```r
devtools::check()
goodpractice::goodpractice(quiet = FALSE)
```

#### Commit

When you've made your changes, write a clear commit message describing what
you've done. If you've fixed or closed an issue, make sure to include keywords
(e.g. `fixes #101`) at the end of your commit message (not in its
title) to automatically close the issue when the PR is merged.

#### Push and pull

Once you've pushed your commit(s) to a branch in _your_ fork, you're ready to
make the pull request. Pull requests should have descriptive titles to remind
reviewers/maintainers what the PR is about. You can easily view what exact
changes you are proposing using either the [Git diff](http://r-pkgs.had.co.nz/git.html#git-status) 
view in RStudio, or the [branch comparison view](https://help.github.com/articles/creating-a-pull-request/) 
you'll be taken to when you go to create a new PR. If the PR is related to an 
issue, provide the issue number and slug in the _description_ using 
auto-linking syntax (e.g. `#15`).

#### Check the docs
Double check the output of the 
[rOpenSci documentation CI](https://dev.ropensci.org/job/ruODK/lastBuild/console) 
for any breakages or error messages.

#### Review, revise, repeat

The latency period between submitting your PR and its review may vary. 
When a maintainer does review your contribution, be sure to use the same 
conventions described here with any revision commits.

### Resources

*  [Happy Git and GitHub for the useR](http://happygitwithr.com/) by Jenny Bryan.
*  [Contribute to the tidyverse](https://www.tidyverse.org/contribute/) covers 
   several ways to contribute that _don't_ involve writing code.
*  [Contributing Code to the Tidyverse](http://www.jimhester.com/2017/08/08/contributing/) by Jim Hester.
*  [R packages](http://r-pkgs.had.co.nz/) by Hadley Wickham.
   *  [Git and GitHub](http://r-pkgs.had.co.nz/git.html)
   *  [Automated checking](http://r-pkgs.had.co.nz/check.html)
   *  [Object documentation](http://r-pkgs.had.co.nz/man.html)
   *  [Testing](http://r-pkgs.had.co.nz/tests.html)
*  [dplyr’s `NEWS.md`](https://github.com/tidyverse/dplyr/blob/master/NEWS.md) 
   is a good source of examples for both content and styling.
*  [Closing issues using keywords](https://help.github.com/articles/closing-issues-using-keywords/) 
   on GitHub.
*  [Autolinked references and URLs](https://help.github.com/articles/autolinked-references-and-urls/) 
   on GitHub.
*  [GitHub Guides: Forking Projects](https://guides.github.com/activities/forking/).

### Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](CODE_OF_CONDUCT.md). By participating in this project you agree to
abide by its terms.

## Maintaining `ruODK`
The steps to prepare a new `ruODK` release are in `data-raw/make_release.R`.
It is not necessary to run them as a contributor, but immensely convenient for
the maintainer to have them there in one place.

## Package maintenance
The code steps run by the package maintainer to prepare a release live at
`data-raw/make_release.R`. Being an R file, rather than a Markdown file like
this document, makes it easier to execute individual lines.

Pushing the Docker image requires privileged access to the Docker repository.
## Test environments
* Local machine
  * Running under: Ubuntu 20.04.2 LTS
  * Platform: x86_64-pc-linux-gnu (64-bit)
  * R version 4.0.4 (2021-02-15)
* AppVeyor CI
  * Running under: Windows Server 2012 R2 x64 (build 9600)
  * Platform: x86_64-w64-mingw32/x64 (64-bit)
  * R version 4.0.3 Patched (2020-11-08 r79411)
* GitHub Actions: R devel, release, oldrel for each of
  * Windows-latest (Windows Server 2019)
  * Windows Server 2016
  * MacOS-lastest (MacOS X Catalina 10.05)
  * Ubuntu 20.04
  * Ubuntu 18.04

## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new release.

Resolved NOTE comments:
* Possibly invalid URLs:
  * The package comparison section in the README contains Markdown badges with 
    CRAN links to packages that are not yet or not any more on CRAN. These
    links are correct, and while they currently do not resolve, they will do so 
    once the packages are (re-)submitted to CRAN. Currently removed.
  * The README contains an ODK Central form OData service URL to illustrate 
    setting up ruODK. The URL redirects to a login screen if followed directly.
    This is expected behaviour. Currently not appearing as warning.
* The PDF version of the manual is now included.
* The example data contains UTF-8 strings. This is a realistic scenario. 
  The note has disappeared after the R version 4 release.
* Test coverage: All functionality supporting the current ODK Central release is 
  covered by tests. 
  The only exception is `form_schema{_parse}`, which supports a breaking 
  change between ODK Central 0.7 and 0.8. The test server runs ODK Central 0.8,
  a production server (used by the package author, but not accessible to other 
  maintainers) runs 0.7 successfully. The tests for v 0.7 use packaged data.
---
name: Bug report
about: ruODK is doing something wrong
title: ''
labels: bug
assignees: florianm

---

## Problem
<!-- Please briefly describe your problem and what output you expect. -->

## Reproducible example
<!-- 
     If this issue refers to a bug or unexpected behaviour of ruODK, 
     please include a minimal reproducible example 
     ([reprex](https://reprex.tidyverse.org/), see also 
     <https://www.tidyverse.org/help/#reprex>). 
     
     Note: If this issue involves an authenticated web request, 
     do not include your credentials, server address, or any identifying
     information. The failing function and the (sanitised) error output alone
     will be a welcome help to narrow down the problem.
-->

```{r}
# insert reprex here
```

<details>
<summary><strong>Session Info</strong></summary>
  
```{r}
# utils::sessionInfo()
```
</details>
---
name: Account request
about: I'd like to run the ruODK package tests
title: 'Account request for ODK Central test instance for ruODK package tests'
labels: account-request
assignees: florianm

---
I would like to run the tests.
Please provide access to the ODK Central test instance for ruODK <https://odkc.dbca.wa.gov.au/>.

[ ] I have sent an email to the `ruODK` package maintainer requesting access.
[ ] I will not use the server to run production campaigns.

<!-- 
  This account request is to interact with ruODK's development.
  
  If you'd like to evaluate ODK Central for your own purposes, 
  visit <https://getodk.org/> to request a free cloud trial. 
-->---
name: Feature request
about: ruODK is missing something
title: ''
labels: feature
assignees: florianm

---

## Feature
<!-- Please briefly describe the new or missing feature. -->
---
title: 'ruODK: An R Client for the ODK Central API'
authors:
- affiliation: 1
  name: Florian Wendelin Mayer
  orcid: 0000-0003-4269-4242
date: "20 July 2020"
output: pdf_document
bibliography: paper.bib
tags:
- database
- open-data
- opendatakit
- odk
- api
- data
- dataset
affiliations:
- index: 1
  name: Department of Biodiversity, Conservation and Attractions, Western Australia
---
![ruODK logo](../../man/figures/ruODK.png)

# Summary
ruODK [@ruodk] is an R Client for the ODK Central API.

# Background

Open Data Kit (ODK, [@odk], [@hartung]) is a suite of open source tools that 
help organizations collect and manage data.

The core ODK tools are [@odkdocs]:

-  ODK Collect, an Android app that replaces paper-based forms.
-  ODK Aggregate, a proven server for data storage and analysis tool.
-  ODK Central, a modern server with a RESTful API.
-  ODK Build, a drag-and-drop form designer.
-  ODK XLSForm, an Excel-based form designer.
-  ODK Briefcase, a desktop tool that pulls and exports data from Aggregate and Collect.

The core workflow of ODK is:

-  A form for digital data capture is designed, either by hand, or using form 
  builders like ODK Build.
-  A data clearinghouse like ODK Aggregate or ODK Central disseminates these
  form templates to authorised data collection devices (Android devices running
  ODK Collect).
-  Once data has been collected, ODK Collect submits the data back to the data
  clearinghouse, e.g. ODK Central.
-  From there, it is the responsibility of the maintainers to export and further
  process and use the data.

At the time of writing, [ODK Central](https://docs.getodk.org/central-intro/)
is being introduced as a replacement for ODK Aggregate.

While the use of ODK Aggregate is well established, the use of ODK Central is new
to most users. As ODK Central provides all of its content through a well-documented
API [@odkapi], there is an opportunity to create API wrappers to assist with the
data retrieval from ODK Central.

ruODK [@ruodk], currently available from GitHub [@github], is the first dedicated
R Client for the ODK Central API.

A typical data flow involves a form builder (ODK Build), a data clearinghouse
(ODK Central), an electronic data capture app (ODK Collect), and a way to
automate data access and analysis (ruODK).

![A typical ODK setup with ODK Build, Central, Collect, and ruODK](https://www.lucidchart.com/publicSegments/view/cd47b81f-04cf-49d7-af3f-eda5f8755203/image.png)

# Scope

`ruODK` aims:

-  To wrap all ODK Central API endpoints with a focus on **data access**. 
  While this is mostly not a hard task, there is still a small barrier to novice
  R users, and some duplication of code.
-  To provide working examples of interacting with the ODK Central API.
-  To provide convenience helpers for the day to day tasks when working with 
  ODK Central data in R: transforming the ODK Central API output into tidy 
  R formats.
  
## Out of scope

-  To wrap "management" API endpoints. The ODK Central GUI already provides a 
  highly capable interface for the management of users, roles, permissions, 
  projects, and forms.
  ODK Central is a [VueJS application](https://github.com/getodk/central-frontend/) 
  working on the "management" API endpoints of the ODK Central back-end.
-  To provide extensive data visualisation capability. 
  We show only minimal examples of data visualisation and presentation, mainly 
  to illustrate the example data.
  
# Typical use cases

## Smaller projects
Smaller projects, such as ephemeral research projects or proof of concept studies, 
may generate a manageable number of form submissions over a limited amount of time.

Here, it may be sufficient to export the entire set of submissions, analyse
and visualise it, and generate some products such as a CSV export of the data,
vector or raster figures and maps for publications, and a rendered report.

The vignettes "odata-api" and "restful-api" show examples of this use case. 
They differ in that the vignette "odata-api" retrieves submissions through the 
OData API, while the vignette "restful-api" retrieves submissions through the 
RESTful API.

## Larger projects
Long-term projects, such as environmental monitoring programs, may generate many
thousand submissions over a long, ongoing period of time.

Such projects typically store their data in dedicated databases. There, the data
can be value-added through QA/QC/review, through integration with other data
sources, and through further data processing. This requires the data to be
regularly and incrementally exported from ODK Central, transformed into the target
data formats (which may differ from ODK Central's output), and loaded into the 
target database.

In this case, ruODK assists with the retrieval of data, inspection of new vs. 
existing submissions, and export into an intermediary format, such as CSV 
snapshots.  

# References---
output: github_document
---
<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/",
  out.width = "100%")
```
# `ruODK`: An R Client for the ODK Central API <img src="man/figures/ruODK2.png" align="right" alt="Especially in these trying times, it is important to ask: ruODK?" width="120" />

<!-- badges: start -->
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3953158.svg)](https://doi.org/10.5281/zenodo.3953158)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![Last-changedate](https://img.shields.io/github/last-commit/ropensci/ruODK.svg)](https://github.com/ropensci/ruODK/commits/main)
[![GitHub issues](https://img.shields.io/github/issues/ropensci/ruodk.svg?style=popout)](https://github.com/ropensci/ruODK/issues/)
[![CI - GitHub Actions](https://github.com/ropensci/ruODK/workflows/tic/badge.svg)](https://github.com/ropensci/ruODK/actions)
[![CI - Appveyor](https://ci.appveyor.com/api/projects/status/1cs19xx0t64bmd2q/branch/master?svg=true)](https://ci.appveyor.com/project/florianm/ruodk/branch/main)
[![Test coverage](https://codecov.io/gh/ropensci/ruODK/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/ruODK)
[![CodeFactor](https://www.codefactor.io/repository/github/ropensci/ruodk/badge)](https://www.codefactor.io/repository/github/ropensci/ruodk)
[![Hosted RStudio with ruODK](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dbca-wa/urODK/master?urlpath=rstudio)
<!-- badges: end -->

Especially in these trying times, it is important to ask "r u ODK?".

`ruODK` is an R client to access and parse data from ODK Central.

[OpenDataKit](https://getodk.org/) (ODK) is 
[free-and open-source software](https://getodk.org/software/) 
that helps millions of people collect data quickly, accurately, offline, 
and at scale. The software is in active use in every country in the world and is 
supported by a large and helpful community.

`ruODK` is a community contribution to the ODK ecosystem, but not directly 
affiliated with ODK.

`ruODK` assumes some familiarity of its users with the ODK ecosystem and workflows.
For a detailed overview, read the extensive [ODK documentation](https://docs.getodk.org/)
and visit the friendly [ODK forum](https://forum.getodk.org/).

[ODK Central](https://docs.getodk.org/central-intro/) is a cloud-based
data clearinghouse for digitally captured data, replacing the older software 
[ODK Aggregate](https://docs.getodk.org/aggregate-intro/). 
ODK Central manages user accounts and permissions, stores form definitions, 
and allows data collection clients like 
[ODK Collect](https://docs.getodk.org/collect-intro/) to connect to it for 
form download and submission upload.

![An ODK setup with ODK Build, Central, Collect, and ruODK](
https://www.lucidchart.com/publicSegments/view/952c1350-3003-48c1-a2c8-94bad74cdb46/image.png)

A typical [ODK workflow](https://docs.getodk.org/#how-is-odk-used):
An XForm is designed e.g. in [ODK Build](https://build.getodk.org/),
[published to ODK Central](https://docs.getodk.org/central-forms/), 
and downloaded onto an Android device running ODK Collect.
After data have been captured digitally using 
[ODK Collect](https://docs.getodk.org/collect-intro/), the data are uploaded
and stored in ODK Central. The next step from there is to extract the data, 
optionally upload it into another data warehouse, and then to analyse and 
generate insight from it. 

While data can be retrieved in bulk through the GUI, ODK Central's API provides 
access to its data and functionality through both an OData and a RESTful API 
with a comprehensive and interactive 
[documentation](https://odkcentral.docs.apiary.io/#reference/odata-endpoints).

`ruODK` is aimed at the technically minded researcher who wishes to access and
process data from ODK Central using the programming language R.

Benefits of using the R ecosystem in combination with ODK:

*  Scalability: Both R and ODK are free and open source software. Scaling to
   many users does not incur license fees.
*  Ubiquity: R is known to many scientists and is widely taught at universities.
*  Automation: The entire data access and analysis workflow can be automated 
   through R scripts.
*  Reproducible reporting (e.g. 
   [Sweave](https://support.rstudio.com/hc/en-us/articles/200552056-Using-Sweave-and-knitr), 
   [RMarkdown](https://rmarkdown.rstudio.com/)), interactive web apps 
   ([Shiny](https://shiny.rstudio.com/)), 
   workflow scaling ([drake](https://docs.ropensci.org/drake/)).
* Rstudio-as-a-Service (RaaS) at 
  [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dbca-wa/urODK/master?urlpath=rstudio)

`ruODK`'s scope:

*  To wrap all ODK Central API endpoints with a focus on **data access**.
*  To provide working examples of interacting with the ODK Central API.
*  To provide convenience helpers for the day to day tasks when working with 
   ODK Central data in R: **data munging** the ODK Central API output into tidy 
   R formats.
  
<!-- TODO: vignette "workflows" -->
`ruODK`'s use cases:

*  Smaller projects: Example [rOzCBI](https://dbca-wa.github.io/rOzCBI/)
   1. Data collection: ODK Collect
   2. Data clearinghouse: ODK Central
   3. Data analysis and reporting: `Rmd` (ruODK)
   4. Publishing and dissemination: [`ckanr`](https://docs.ropensci.org/ckanr/), 
      [`CKAN`](https://ckan.org/)
     
*  Larger projects:
   1. Data collection: ODK Collect
   2. Data clearinghouse: ODK Central
   3. ETL pipeline into data warehouses: `Rmd` (ruODK)
   4. QA: in data warehouse
   5. Reporting: `Rmd`
   6. Publishing and dissemination: [`ckanr`](https://docs.ropensci.org/ckanr/),
      [`CKAN`](https://ckan.org/)
  
Out of scope:

*  To wrap "management" API endpoints. ODK Central is a
   [VueJS/NodeJS application](https://github.com/opendatakit/central-frontend/) 
   which provides a comprehensive graphical user interface for the management of 
   users, roles, permissions, projects, and forms.
*  To provide extensive data visualisation. We show only minimal examples of data 
   visualisation and presentation, mainly to illustrate the example data.
   Once the data is in your hands as tidy tibbles... urODK!

## A quick preview
<img src="man/figures/odata.svg" alt="ruODK screencast" width="100%" />

## Install
You can install the latest release of `ruODK` from the 
[rOpenSci R-Universe](https://ropensci.r-universe.dev):

```{r r-universe, eval = FALSE}
# Enable the rOpenSci universe
options(repos = c(ropensci = 'https://ropensci.r-universe.dev',
                  CRAN = 'https://cloud.r-project.org'))
install.packages('ruODK')
```

Alternatively, you can install the development version from the `main` branch.

```{r gh-installation, eval = FALSE}
if (!requireNamespace("remotes")) install.packages("remotes")
# Full install
remotes::install_github(
   "ropensci/ruODK@main", 
   dependencies = TRUE, 
   upgrade = "always",
   build_vignettes = TRUE)

# Minimal install without vignettes
remotes::install_github(
   "ropensci/ruODK@main", 
   dependencies = TRUE, 
   upgrade = "ask",
   build_vignettes = FALSE)
```

If the install fails, read the error messages carefully and install any unmet
dependencies (system libraries or R packages).

If the install fails on building the vignettes, you can set 
`build_vignettes=FALSE` and read the vignettes from the online docs instead.

If the installation still fails, or the above does not make any sense, 
feel free to submit a [bug report](https://github.com/ropensci/ruODK/issues/new/choose).

## ODK Central
### Access to an ODK Central instance
First, we need an ODK Central instance and some data to play with!

Either [request a free trial](https://getodk.org/#odk-cloud)
or follow the [setup instructions](https://docs.getodk.org/central-intro/)
to build and deploy your very own ODK Central instance.

### ODK Central setup
The ODK Central [user manual](https://docs.getodk.org/central-using/) 
provides up-to-date descriptions of the steps below.

*  [Create a web user account](https://docs.getodk.org/central-users/#creating-a-web-user) 
   on an ODK Central instance. Your username will be an email address.
*  [Create a project](https://docs.getodk.org/central-projects/) and give 
   the web user at least
   [read permissions](https://docs.getodk.org/central-projects/#managing-project-managers).
*  Create an XForm, e.g. using ODK Build, or use the 
   [example forms](https://github.com/ropensci/ruODK/tree/master/inst/extdata) 
   provided by `ruODK`. The `.odkbuild` versions can be loaded into 
   [ODK Build](https://build.getodk.org/), while the `.xml` versions can be 
   directly imported into ODK Central.
*  [Publish the form](https://docs.getodk.org/central-forms/)
   to ODK Central.
*  Collect some data for this form on ODK Collect and let ODK Collect submit the
   finalised forms to ODK Central.

## Configure `ruODK`
For all available detailed options to configure authentication for `ruODK`, read 
[`vignette("setup", package = "ruODK")`](https://docs.ropensci.org/ruODK/articles/setup.html).

## Use ruODK
A detailed walk-through with some data visualisation examples is available in the 
[`vignette("odata-api", package="ruODK")`](https://docs.ropensci.org/ruODK/articles/odata-api.html).

See also [`vignette("restful-api", package="ruODK")`](https://docs.ropensci.org/ruODK/articles/restful-api.html)
for examples using the alternative RESTful API.

## Try ruODK
[![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/dbca-wa/urODK/master?urlpath=rstudio) 
will launch a disposable, hosted RStudio instance with `ruODK` installed and
the companion package [`urODK`](https://github.com/dbca-wa/urODK) opened as 
starting point for a hands-on workshop or instant demo of `ruODK` usage. 

Create a new RMarkdown workbook from `ruODK` template "ODK Central via OData" 
and follow the instructions within.

## Contribute
Contributions through [issues](https://github.com/ropensci/ruODK/issues) and PRs 
are welcome!

See the [contributing guide](https://docs.ropensci.org/ruODK/CONTRIBUTING.html)
on best practices and further readings for code contributions.

## Attribution
`ruODK` was developed, and is maintained, by Florian Mayer for the Western Australian 
[Department of Biodiversity, Conservation and Attractions (DBCA)](https://www.dbca.wa.gov.au/).
The development was funded both by DBCA core funding and external funds from
the [North West Shelf Flatback Turtle Conservation Program](https://flatbacks.dbca.wa.gov.au/).

To cite package `ruODK` in publications use:

```{r citation}
citation("ruODK")
```

## Acknowledgements
The Department of Biodiversity, Conservation and Attractions (DBCA) acknowledges 
the traditional owners of country throughout Western Australia and their continuing 
connection to the land, waters and community. We pay our respects to them, their 
culture and to their Elders past and present.

This software was created on Whadjuk boodja (ground) both as a contribution to 
the ODK ecosystem and for the conservation of the biodiversity of 
Western Australia, and in doing so, caring for country.

## Package functionality
See [`vignette("comparison", package="ruODK")`](https://docs.ropensci.org/ruODK/articles/comparison.html)
for a comprehensive comparison of ruODK to other software packages from
both an ODK and an OData angle.
---
title: "Report Title"
author: "Report Author"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_depth: 3
    toc_float: true
    fig_width: 10
    fig_height: 6
    code_folding: show
    theme: lumen
  pdf_document:
    latex_engine: xelatex
  word_document: default
---
<!--  
This Rmd template provides an example workflow to access submissions and 
attachments of a given ODK Central form using its OData API and disseminate
data, attachments, and this report via CKAN or Google Drive.

You will need:
* An ODK Central form's OData Service URL (ending in .svc), accessed from the 
"Analyse using OData" button on the form submission tab.
* Username and password of an ODK Central webuser authorised to view this form.
* Optional: base URL and write-permitted API key for a CKAN instance.
* Optional: a Google account and the means to authenticate (e.g. phone for 2FA).
-->

<!--
Step 1: Setup ruODK, ckanr
Secrets are set as environment variables to keep them outside of this document. 
Read vignette("setup") for setup and authentication options.
-->
<!-- 1a. Open the .Renviron file:  -->
```{r edit_renviron, echo=FALSE, eval=FALSE}
usethis::edit_r_environ()
```

<!-- 1b. Paste the following block using your own, working credentials:  -->
```{r paste_env_vars, echo=FALSE, eval=FALSE}
# ODK Central web user with read-permitted role on project
ODKC_UN="my@email.com"
ODKC_PW="my-odkc-password"

# CKAN user with permissions to create a dataset
CKAN_URL="https://demo.ckan.org"
CKAN_KEY="my-ckan-api-key"
```
<!-- 1c: Restart the R session to load these environment variables. -->

<!-- 
* Configure ODK Central base URL (url), project ID (pid), and form ID (fid) 
in one step using the OData Service URL.
* Paste your own Service URL over `svc="..."`.
* Paste your local timezone over `Australia/Perth`.
* Adjust `loc` to a local subfolder for downloaded attachments.
-->
```{r setup, message=FALSE, warning=FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
# Data wrangling
# library(stringr)
library(dplyr)
library(fs)
library(ruODK)

# Visualisation
library(skimr)
library(knitr)
library(DT)

# Dissemination
library(readr)
# library(ckanr)
# library(googledrive)

# Spatial: see also vignette "spatial"
library(leaflet)
# library(rgeos)
# library(sf)

# Set ruODK defaults to an ODK Central form's OData Service URL
ruODK::ru_setup(
  svc=paste0("https://my-odkc-instance.com/v1/projects/1/forms/form_id.svc"),
  tz = "Australia/Perth", # your local timezone
  odkc_version = 1.0,     # your ODKC version, only needed if <=0.7
  verbose = TRUE)

loc <- fs::path("media")
fs::dir_create(loc)
```

# Download data
<!-- Get table names. `ft` will contain one row per separate table. -->
```{r svc_get}
ft <- ruODK::odata_service_get()
ft %>% knitr::kable(.)
```

<!-- 
Get and parse submissions.
Repeat for each form table name in ft$url. 
The `submissions_id` of a nested table links to the main table's `id`.
Once happy with the parsing results, verbose can be set to FALSE.
See vignette "spatial" for more operations on geo fields.
-->
```{r dl_submissions}
data <- ruODK::odata_submission_get(
  table = ft$url[1], 
  local_dir = loc, 
  wkt = TRUE
)

data_sub1 <- ruODK::odata_submission_get(
  table = ft$url[2], 
  local_dir = loc
) %>%
  dplyr::left_join(
    data, 
    by = c("submissions_id" = "id")
  )

data_sub2 <- ruODK::odata_submission_get(
  table = ft$url[3], 
  local_dir = loc
) %>%
  dplyr::left_join(
    data, 
    by = c("submissions_id" = "id")
)

# repeat for any remaining nested tables
```

# Analyse and visualise data
<!-- 
  Summarise, analyse, visualise the tibble(s) data (, data_sub1, data_sub2) 
-->
```{r data_vis}
skimr::skim(data)
dplyr::glimpse(data)
DT::datatable(head(data))
```

```{r data_map}
leaflet::leaflet(width = 800, height = 600) %>%
  leaflet::addProviderTiles("OpenStreetMap.Mapnik", group = "Place names") %>%
  leaflet::addProviderTiles("Esri.WorldImagery", group = "Aerial") %>%
  leaflet::clearBounds() %>%
  leaflet::addAwesomeMarkers(
    data = data,
    # 
    # # Adjust to your coordinate field names
    # 
    lng = ~location_longitude, 
    lat = ~location_latitude, 
    icon = leaflet::makeAwesomeIcon(text = "Q", markerColor = "red"),
    # 
    # # With your own field names
    # 
    # label = ~ glue::glue("{location_area_name} {encounter_start_datetime}"),
    # popup = ~ glue::glue(
    #   "<h3>{location_area_name}</h3>",
    #   "Survey start {encounter_start_datetime}</br>",
    #   "Reporter {reporter}</br>",
    #   "Device {device_id}</br>",
    #   "<h5>Site</h5>",
    #   '<div><img src="{habitat_morphological_type_photo}"',
    #   ' height="150px" alt="My photo"></img></div>'
    # ),
    # 
    # # If there are many submissions, cluster markers for performance:
    clusterOptions = leaflet::markerClusterOptions()
  ) %>%
  leaflet::addLayersControl(
    baseGroups = c("Place names", "Aerial"),
    options = leaflet::layersControlOptions(collapsed = FALSE)
  )
```

# Export
<!--
The form submissions are now extracted and visualised. What's next:

* Save data to local files (e.g. CSV).
* Compile report (e.g. to HTML).
* Compress all outputs as ZIP.
* Upload these artifacts to a CKAN data catalogue.
* Upload same artifacts to Google Drive.

Notes: 

* Generate the HTML report once off without the next chunk (`eval=F`), 
as the chunk refers to the rendered output file (HTML) before the file is 
created initially.
* Run report always twice to generate (run 1) and upload (run 2) the latest HTML.
-->
```{r data_export, eval=FALSE}
#------------------------------------------------------------------------------#
# Prepare report and products as local files
#
rep_fn <- "my_report.html" # The file name you save this template under
data_fn <- here::here(loc, "data.csv") %>% as.character()           # Main data
data_sub1_fn <- here::here(loc, "data_sub1.csv") %>% as.character() # Sub table 1
data_sub2_fn <- here::here(loc, "data_sub2.csv") %>% as.character() # Sub table 2
zip_fn <- "products.zip" # Attachments as one zip file (top level)

# Write data tbls to CSV files
readr::write_csv(data, path = data_fn)
readr::write_csv(data_sub1, path = data_sub1_fn)
readr::write_csv(data_sub2, path = data_sub2_fn)

# Compress everything into `zip_fn`, retain relative path to `loc`
zip(zipfile = zip_fn, files = fs::dir_ls(loc))

#------------------------------------------------------------------------------#
# CKAN
#
# Upload to a CKAN data catalogue (needs url and API key of a write permitted user)
# See ROpenSci package ckanr
ckanr::ckanr_setup(url = Sys.getenv("CKAN_URL"), key = Sys.getenv("CKAN_KEY"))
ckan_ds_name <- "my-ckan-dataset-slug"

# Run once to create resources on an existing dataset, then comment out
d <- ckanr::package_show(ckan_ds_name)
res_data_main <- ckanr::resource_create(
  package_id = d$id, name="Main data", upload = data_fn)
res_data_sub1 <- ckanr::resource_create(
  package_id = d$id, name="Nested data table 1", upload = data_sub1_fn)
res_data_sub2 <- ckanr::resource_create(
  package_id = d$id, name="Nested data table 2", upload = data_sub2_fn)
# add remaining tables
if (fs::file_exists(rep_fn)){
res_report <- ckanr::resource_create(
  package_id = d$id, name="Data report", upload = rep_fn)
}

if (fs::file_exists(zip_fn)){
res_zip <- ckanr::resource_create(
  package_id = d$id, name="All data and attachments", upload = zip_fn)
}
# Paste res_data_main$id over RID and keep here, repeat for each resource
r <- ckanr::resource_update(res_data_main$id, path = data_fn)
r <- ckanr::resource_update(res_data_sub1$id, path = data_sub1_fn)
r <- ckanr::resource_update(res_data_sub2$id, path = data_sub2_fn)
if (fs::file_exists(rep_fn)){
  r <- ckanr::resource_update(res_report$id, path = rep_fn)
}
r <- ckanr::resource_update(res_zip$id, path = zip_fn)

#------------------------------------------------------------------------------#
# Google Drive
#
# Run once per machine, then comment out:
googledrive::drive_auth(use_oob = TRUE)

# Upload to Google Drive
gd_fn <- "My Google Drive folder name"
googledrive::drive_ls(gd_fn) %>% googledrive::drive_rm(.)  # Wipe older outputs
if (fs::file_exists(rep_fn)){
  googledrive::drive_upload(rep_fn, path=rep_fn)           # Report as HTML
}
googledrive::drive_upload(data_fn, path=data_fn)           # Main data as CSV
googledrive::drive_upload(data_sub1_fn, path=data_sub1_fn) # Sub table 1 as CSV
googledrive::drive_upload(data_sub2_fn, path=data_sub2_fn) # Sub table 2 as CSV
googledrive::drive_upload(zip_fn, path=zip_fn)             # All outputs as ZIP
```
---
title: "Comparison of related software packages"
description: >
  An overview of related software packages in the ODK / OData space.
output: 
  rmarkdown::html_vignette:
    self_contained: false
vignette: >
  %\VignetteIndexEntry{Comparison of related software packages}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

There are several other R packages interacting with the ODK ecosystem, and/or
[OData](https://www.odata.org/).

## Comparison of ODK related software packages (non-ODK core)

| Package                         | [`ruODK`](https://docs.ropensci.org/ruODK/) | [`odkr`](https://rapidsurveys.io/odkr/) | [`odk`](https://cran.r-project.org/package=odk) | [`odkmeta`](https://github.com/nap2000/odkmeta) | [`koboloadeR`](https://unhcr.github.io/koboloadeR/docs/index.html) | [Pentaho Kettle tutorial](https://github.com/schemetrica/automating-data-delivery-odk-central)
|------------------------------|---------------|---------------|---------------|---------------|---------------|---------------|
| Elevator pitch               | "[ckanr](https://github.com/ropensci/ckanr) for ODK Central"  | "Drive ODK Briefcase through R" | "Export ODK Aggregate to SPSS" | "Export ODK Aggregate to STATA" | "Metapackage for the extended ODK ecosystem" | "What ruODK does, but as GUI" |
| Last commit                  | [![Last-changedate](https://img.shields.io/github/last-commit/ropensci/ruODK.svg)](https://github.com/ropensci/ruODK/commits/main) | [![Last-changedate](https://img.shields.io/github/last-commit/rapidsurveys/odkr.svg)](https://github.com/rapidsurveys/odkr/commits/master) |    Nov 2017    |  [![Last-changedate](https://img.shields.io/github/last-commit/nap2000/odkmeta.svg)](https://github.com/nap2000/odkmeta/commits/master) | [![Last-changedate](https://img.shields.io/github/last-commit/unhcr/koboloadeR.svg)](https://github.com/unhcr/koboloadeR/commits/master) | [![Last-changedate](https://img.shields.io/github/last-commit/schemetrica/automating-data-delivery-odk-central.svg)](https://github.com/schemetrica/automating-data-delivery-odk-central/commits/master) |
| Website                      | [![](https://img.shields.io/static/v1?label=code&message=GitHub&color=brightgreen)](https://github.com/ropensci/ruODK) [![docs](https://img.shields.io/static/v1?label=docs&message=pkgdown&color=brightgreen)](https://docs.ropensci.org/ruODK/) | [![](https://img.shields.io/static/v1?label=code&message=GitHub&color=brightgreen)](https://github.com/rapidsurveys/odkr) [![docs](https://img.shields.io/static/v1?label=docs&message=pkgdown&color=brightgreen)](https://rapidsurveys.io/odkr/)   | [![docs](https://img.shields.io/static/v1?label=docs&message=rdrr.io&color=brightgreen)](https://rdrr.io/cran/odk/) |    [![](https://img.shields.io/static/v1?label=code&message=GitHub&color=brightgreen)](https://github.com/nap2000/odkmeta) | [![](https://img.shields.io/static/v1?label=code&message=GitHub&color=brightgreen)](https://github.com/unhcr/koboloadeR) [![docs](https://img.shields.io/static/v1?label=docs&message=pkgdown&color=brightgreen)](https://unhcr.github.io/koboloadeR/docs/index.html) | [![](https://img.shields.io/static/v1?label=code&message=GitHub&color=brightgreen)](https://github.com/schemetrica/automating-data-delivery-odk-central) |
| Test coverage             | ![tic](https://github.com/ropensci/ruODK/workflows/tic/badge.svg)  [![codecov](https://codecov.io/gh/ropensci/ruODK/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/ruODK) [![Build status](https://ci.appveyor.com/api/projects/status/1cs19xx0t64bmd2q/branch/main?svg=true)](https://ci.appveyor.com/project/florianm/ruodk/branch/main) | [![codecov](https://codecov.io/gh/rapidsurveys/odkr/branch/master/graph/badge.svg)](https://codecov.io/gh/rapidsurveys/odkr) | ❌ | In repo | [![Travis build status](https://travis-ci.org/unhcr/koboloadeR.svg?branch=gh-pages)](https://travis-ci.org/unhcr/koboloadeR) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/unhcr/koboloadeR?branch=gh-pages&svg=true)](https://ci.appveyor.com/project/unhcr/koboloadeR) [![codecov](https://codecov.io/gh/unhcr/koboloadeR/branch/gh-pages/graph/badge.svg)](https://codecov.io/gh/unhcr/koboloadeR) | NA |
| Working examples             | README, 3 vignettes, pkgdown, Rmd templates | README, pkgdown | CRAN PDF | README | README, 9 vignettes, shiny apps, pkgdown | Tutorial with screenshots |
| Available on CRAN            | [![CRAN status](https://www.r-pkg.org/badges/version/ruODK)](https://cran.r-project.org/package=ruODK) | [![CRAN status](https://www.r-pkg.org/badges/version/odkr)](https://cran.r-project.org/package=odkr) |  [![version](http://www.r-pkg.org/badges/version/odk)](https://CRAN.R-project.org/package=odk)  | NA | NA | NA |
| Technologies            | Tidyverse R, XForms | Base R | Base R | Stata | R metapackage, XlsForms | Pentaho Kettle GUI |
| External dependencies        | None | Java, ODK Briefcase | SPSS | Stata | Java, ODK Briefcase, wraps `odkr` | [Pentaho Kettle](http://www.ibridge.be/), Java |
| Affiliation                  |  [ROpenSci](https://ropensci.org/) | [RapidSurveys](https://rapidsurveys.io/) |  [Muntashir-Al-Arefin](https://stackoverflow.com/users/8875690/muntashir-al-arefin) | [ODK Central developer Matt White](https://github.com/matthew-white) | [UNHCR](https://github.com/unhcr) | [Schemetrica](https://github.com/schemetrica) |
| Covers ODK Central OData API        | ✅  |  ❌ |  ❌  | ❌  | ❌ | ✅  |
| Covers ODK Central REST API         | ✅  |  ❌ |  ❌  | ❌  | ❌ | ❌  |
| Covers ODK Central bulk export      | ✅  |  ❌ |  ❌  | ✅  | ❌ | ✅  |  
| Covers ODK Central OpenRosa API     | ❌  no need, gets all data through OData/REST API |  ✅  via ODK Briefcase |  ❌  | ❌  | ✅  via ODK Briefcase  | ✅  |  
| Data post-processing                | ✅  |  ✅ |  ❌ | ✅  | ✅ | ✅  |  
| Data visualisation examples         | ✅  |  ❌ | ❌  | ❌  | ✅  | ❌ |  
 
In summary:

`ruODK` provides a dependency-free interface to ODK Central.

`koboloadeR` is a metapackage containing lots of ancillary packages, with some
heavy dependencies on Java and ODK Briefcase (which in turn can access ODK Central).
Although built around the XlsForm standard and paradigm, `koboloadeR` is well worth 
exploring as a larger context to data wrangling in the ODK ecosystem.

Schemetrica's tutorial illustrates data ETL from ODK Central and deserves a special
mention, as it is both very recent and aimed specifically against ODK Central.
The GUI paradigm of Pentaho Kettle addresses a different audience to the scripting
paradigm of `ruODK`. It should be mentioned that Kettle's composable data 
manipulation steps can be used for many other use cases apart from ODK Central.

## Comparison of OData related R packages

| Package                         | [`ruODK`](https://docs.ropensci.org/ruODK/) | [`odataR`](https://github.com/HanOostdijk/odataR) | [`cbsodataR`](https://github.com/edwindj/cbsodataR) | [`OData`](https://cran.r-project.org/package=OData) | [OData JDBC R tutorial](https://www.cdata.com/kb/tech/odata-jdbc-r.rst)
|------------------------------|---------------|---------------|---------------|---------------|--------------|
| Elevator pitch               | "[ckanr](https://github.com/ropensci/ckanr) for ODK Central"  | "OData client for https://opendata.cbs.nl (and similar)" | "OData client for https://www.cbs.nl" | "Minimal OData example" | "Minimal RJDBC example" | 
| Last commit                  | [![Last-changedate](https://img.shields.io/github/last-commit/ropensci/ruODK.svg)](https://github.com/ropensci/ruODK/commits/main) | [![Last-changedate](https://img.shields.io/github/last-commit/HanOostdijk/odataR.svg)](https://github.com/HanOostdijk/odataR/commits/master) | [![Last-changedate](https://img.shields.io/github/last-commit/edwindj/cbsodataR.svg)](https://github.com/edwindj/cbsodataR/commits/master) | Dec 2016 | ❓ |
| Website | [![](https://img.shields.io/static/v1?label=code&message=GitHub&color=brightgreen)](https://github.com/ropensci/ruODK) [![docs](https://img.shields.io/static/v1?label=docs&message=pkgdown&color=brightgreen)](https://docs.ropensci.org/ruODK/) | [![](https://img.shields.io/static/v1?label=code&message=GitHub&color=brightgreen)](https://github.com/HanOostdijk/odataR) | [![docs](https://img.shields.io/static/v1?label=code&message=CRAN&color=green)](https://cran.r-project.org/package=cbsodataR) [![](https://img.shields.io/static/v1?label=code&message=GitHub&color=brightgreen)](https://github.com/edwindj/cbsodataR) [![docs](https://img.shields.io/static/v1?label=docs&message=pkgdown&color=brightgreen)](https://edwindj.github.io/cbsodataR/) | [![](https://img.shields.io/static/v1?label=code&message=CRAN&color=green)](https://cran.r-project.org/package=OData) | [![docs](https://img.shields.io/static/v1?label=docs&message=html&color=brightgreen)](https://www.cdata.com/kb/tech/odata-jdbc-r.rst) |
| Test coverage            | ![tic](https://github.com/ropensci/ruODK/workflows/tic/badge.svg)  [![Build status](https://ci.appveyor.com/api/projects/status/1cs19xx0t64bmd2q/branch/main?svg=true)](https://ci.appveyor.com/project/florianm/ruodk/branch/main) [![codecov](https://codecov.io/gh/ropensci/ruODK/branch/main/graph/badge.svg)](https://codecov.io/gh/ropensci/ruODK) | ❌ | [![Travis-CI Build Status](https://travis-ci.org/edwindj/cbsodataR.png?branch=master)](https://travis-ci.org/edwindj/cbsodataR) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/edwindj/cbsodatar?branch=master)](https://ci.appveyor.com/project/edwindj/cbsodatar) | ❌ | ❌ |
| Targets ODK Central | ✅ | ❌ | ❌ | ❌ | ❌ |
| Works with ODK Central | ✅ | ❓ |  ❓  | ❌ | ❌ |
| Data wrangling helpers for post-processing | ✅ | some |  some  | ❌ | ❌ |
| Actively maintained to work against ODK Central | ✅ | ❌ | ❌ | ❌ | ❌ |
| Technologies | R, httr, xml2, tidyr, purrr | R, jsonlite, tidyverse | R, tidyverse | R, XML, RJSONIO | R, RJDBC, Java |
| External dependencies | ✅  None  |  ✅  None   | ✅  None | ✅  None | ❌ JDBC, Java | 
| Available on CRAN     | [![CRAN status](https://www.r-pkg.org/badges/version/ruODK)](https://cran.r-project.org/package=ruODK) | NA |  [![version](http://www.r-pkg.org/badges/version/cbsodataR)](https://CRAN.R-project.org/package=cbsodataR) | [![version](http://www.r-pkg.org/badges/version/OData)](https://CRAN.R-project.org/package=OData) |  NA |

In summary: 

`ruODK` is the only R package explicitly aimed at ODK Central's OData and 
RESTful API, as well as providing context and helpers around specific recurring 
data wrangling tasks.

The value of OData lies in its self-descriptive nature, which allows tools to 
introspect the data structures and types. Both GUI-driven tools like MS PowerBI
and `ruODK` use this introspection to assist users in wrangling their own data.

The script-based approach of `ruODK` allows to automate the data extraction, 
transformation, and reporting pipeline, and therefore provide reproducible 
reporting.---
title: "Spatial data"
description: >
  Examples and conversation starters for spatial data.
output: 
  rmarkdown::html_vignette:
    self_contained: false
vignette: >
  %\VignetteIndexEntry{Spatial data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include=FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

For forms with spatial types, such as geopoint, geotrace, or geoshape, 
`ruODK` gives two options to access the captured spatial data.

Firstly, to make spatial data as simple and accessible as possible, 
`ruODK` extracts the lat/lon/alt/acc from geopoints, as well as from 
the first coordinate of geotraces and geoshapes into separate columns.
This works for both GeoJSON and WKT.
The extracted columns are named as the original geofield appended with 
`_latitude`, `_longitude`, `_altitude`, and `_accuracy`, respectively.

Secondly, this vignette demonstrates how to turn the spatial data types returned 
from ODK Central into spatially enabled objects. 
To do so, we have to address two challenges.

The first challenge is to select which of the potentially many spatial fields
which an ODK form can capture shall be used as the primary geometry of a native
spatial object, such as an `sf` SimpleFeature class.
If several spatial fields are captured, it is up to the user to choose which 
field to use as primary geometry.

The second challenge is that the parsed data from ODK Central is a plain table 
(`tbl_df`) in which some columns contain spatial data. Well Know Text (WKT) is 
parsed as text columns, whereas GeoJSON (nested JSON) is parsed as list columns.

Most spatial packages require either atomic coordinates in separate columns,
which works well for points (latitude, longitude, altitude), or the data to be
spatially enabled.
This vignette shows how to transform a `tbl_df` with a column containing (point, 
line, or polygon) WKT into a spatially enabled `sf` object.

```{r, message=FALSE, warning=FALSE}
library(ruODK)

# Visualisation
library(leaflet)
library(ggplot2)
library(lattice)
# library(tmap) # Suggested but not included here yet

# Spatial
can_run <- require(sf) && require(leafem) && require(mapview)
# Fix https://github.com/r-spatial/mapview/issues/313
# See also https://github.com/r-spatial/mapview/issues/312
# option 'fgb' requires GDAL >= 3.1.0
if(require(mapview)) mapview::mapviewOptions(
  fgb = FALSE,
  basemaps = c("Esri.WorldImagery", 
               "Esri.WorldShadedRelief", 
               "OpenTopoMap",
               "OpenStreetMap"),
  layers.control.pos = "topright")
```

## Data

The original data shown in this vignette are hosted on a ODK Central server which 
is used for the `ruODK` package tests. The form we show here contains every 
spatial widget supported by ODK Build for every supported spatial field type.

With working credentials to the ODK Central test server we could download the data directly.

```{r dl_submissions, eval=FALSE}
# Set ruODK defaults to an ODK Central form, choose tz and verbosity
ruODK::ru_setup(
  url = get_test_url(),
  pid = get_test_pid(),
  fid = get_test_fid_wkt(),
  un = get_test_un(), 
  pw = get_test_pw(),
  odkc_version = 0.8,
  tz = "Australia/Perth",
  verbose = TRUE)
data_wkt <- ruODK::odata_submission_get(wkt = TRUE)
data_gj <- ruODK::odata_submission_get(wkt = FALSE)
```

To allow users to build this vignette without credentials to the ODK Central 
test server, `ruODK` provides above form data also as package data. 

```{r}
data("geo_wkt", package = "ruODK")
data("geo_gj", package = "ruODK")
```


## Map geopoints
We can turn data with a text column containing WKT into an `sf` (SimpleFeatures) 
object.

In addition, we can leave the `tbl_df` as non-spatial object, and instead 
use the separately extracted latitude, longitude, altitude, and accuracy 
individually e.g. to plot a Leaflet map.

```{r, fig.width=9, fig.height=7, eval=can_run}
geo_sf_point <- geo_wkt %>% sf::st_as_sf(wkt="point_location_point_gps")
# ODK Collect captures WGS84 (EPSG:4326)
sf::st_crs(geo_sf_point) <- 4326
```

### Mapview using sf
See also further information under heading "Outlook".
```{r}
mapview::mapview(geo_sf_point, col.regions = sf::sf.colors(10))
```

### GGplot using sf
```{r}
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = geo_sf_point, ggplot2::aes(fill = device_id))
```

### Leaflet using sf
```{r}
leaflet::leaflet(data = geo_sf_point) %>% 
  leaflet::addTiles() %>% 
  leaflet::addMarkers(label = ~ device_id, popup = ~ device_id)
```

### Leaflet using extracted coordinate components in tbl_df
```{r}
leaflet::leaflet(data = geo_wkt) %>% 
  leaflet::addTiles() %>% 
  leaflet::addMarkers(
    lng = ~ point_location_point_gps_longitude,
    lat = ~ point_location_point_gps_latitude,
    label = ~ device_id, 
    popup = ~ device_id)
```


## Map geotraces (lines)
We use `sf::st_as_sf` on a text column containing a WKT geotrace.

```{r, fig.width=9, fig.height=7, eval=can_run}
geo_sf_line <- geo_wkt %>% sf::st_as_sf(wkt="path_location_path_gps")
# ODK Collect captures WGS84 (EPSG:4326)
sf::st_crs(geo_sf_line) <- 4326
```

### Mapview using sf
```{r}
mapview::mapview(geo_sf_line, col.regions = sf::sf.colors(10))
```

### GGplot using sf
```{r}
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = geo_sf_line, ggplot2::aes(fill = device_id))
```

### Leaflet using sf and extracted coordinates
You can show either first extracted coordinate components from plain `tbl_df`
or show the full polygons using `{leafem}`.
See the mapview article on [extra functionality](https://r-spatial.github.io/mapview/articles/articles/mapview_06-add.html).

```{r}
leaflet::leaflet(data = geo_wkt) %>% 
  leaflet::addTiles() %>% 
  leaflet::addMarkers(
    lng = ~ path_location_path_gps_longitude,
    lat = ~ path_location_path_gps_latitude,
    label = ~ device_id, 
    popup = ~ device_id)%>% 
  leafem::addFeatures(geo_sf_line, label = ~ device_id, popup = ~ device_id)
```

## Map geoshapes (polygons)
Again, we'll use `sf::st_as_sf` but select a WKT geoshape column.

```{r, fig.width=9, fig.height=7, eval=can_run}
geo_sf_poly <- geo_wkt %>% sf::st_as_sf(wkt="shape_location_shape_gps")
# ODK Collect captures WGS84 (EPSG:4326)
sf::st_crs(geo_sf_poly) <- 4326
```

### Mapview using sf
```{r}
mapview::mapview(geo_sf_poly, col.regions = sf::sf.colors(10))
```

### GGplot using sf
```{r}
ggplot2::ggplot() + 
  ggplot2::geom_sf(data = geo_sf_poly, ggplot2::aes(fill = device_id))
```

### Leaflet using sf and extracted coordinates
You can show either first extracted coordinate components from plain `tbl_df`
or show the full polygons using `{leafem}`.
See the mapview article on [extra functionality](https://r-spatial.github.io/mapview/articles/articles/mapview_06-add.html).
```{r}
leaflet::leaflet(data = geo_wkt) %>% 
  leaflet::addTiles() %>% 
  leaflet::addMarkers(
    lng = ~ shape_location_shape_gps_longitude,
    lat = ~ shape_location_shape_gps_latitude,
    label = ~ device_id, 
    popup = ~ device_id) %>% 
  leafem::addFeatures(geo_sf_poly, label = ~ device_id, popup = ~ device_id)
```

# Outlook
The above examples show how to turn spatial data into an `sf` object, and give
very rudimentary visualisation examples to bridge the gap between
spatial data coming from ODK and creating maps and further spatial analyses in R.

See the [sf homepage](https://r-spatial.github.io/sf/) for more context and examples. 
The [sf cheatsheet](https://github.com/rstudio/cheatsheets/blob/master/sf.pdf) 
deserves a spatial mention.

Review the options for 
[mapview popups](https://r-spatial.github.io/mapview/articles/articles/mapview_04-popups.html)
and the whole [mapview](https://r-spatial.github.io/mapview/index.html) homepage for a 
comprehensive overview of mapview.

The powerful visualisation package [`tmap`](https://cran.r-project.org/package=tmap)
supports `sf` objects and produces both printable and static maps as well as interactive `leaflet` maps.
See the [vignette "Get started"](https://cran.r-project.org/web/packages/tmap/vignettes/tmap-getstarted.html).

There are several other good entry points for all things R and spatial, 
including but not limited to:

* The [R Spatial CRAN Task View](https://cran.r-project.org/view=Spatial)
* The [RSpatial](https://rspatial.org/) website
* [Geospatial data in R and beyond](https://www.maths.lancs.ac.uk/~rowlings/Teaching/UseR2012/) 
  by [Barry Rowlingson](http://barry.rowlingson.com/)
* [GIS with R](https://www.jessesadler.com/post/gis-with-r-intro/) 
  by [Jesse Sadler](https://www.jessesadler.com/page/cv/)
* GIS and mapping by [Olivier Gimenez](https://oliviergimenez.github.io/): 
  [Slides](https://oliviergimenez.github.io/intro_spatialR/) and
  [code](https://github.com/oliviergimenez/intro_spatialR)
  
The above list of examples and resources is far from comprehensive.
Feel free to [contribute or suggest](https://github.com/ropensci/ruODK/issues) 
other working examples for turning data from `ruODK` into spatial formats.
---
title: "Accessing the OData API"
description: >
  Accessing submission data via the OData pathway.
output: 
  rmarkdown::html_vignette:
    self_contained: false
vignette: >
  %\VignetteIndexEntry{Accessing the OData API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
resource_files:
- media/1604290006239.jpg
- media/1604290049411.jpg
- media/1604290379613.jpg
- media/1604290400891.jpg
- media/1604290499008.jpg
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

```{r pkgs}
# This vignette requires the following packages:
library(DT)
library(leaflet)
# library(listviewer)
library(magrittr)
library(tibble)
library(tidyr)
library(lubridate)
library(knitr)
library(dplyr)
library(ggplot2)
library(ruODK)
```
<!-- TODO split out vignettes -->
<!-- * This vignette: the typical workflow -->
<!-- * New vignette: spatial is not special -->
<!-- * New vignette: DIY rectangling -->
<!-- * New vignette: data visualisation (refer to that from odata and restful api) -->
This vignette demonstrates `ruODK`'s workflow to extract data from ODK Central's
OData service endpoint, and to prepare the data and the media attachments for
further analysis and visualisation.

The demonstrated workflow is roughly equivalent to ODK Central's "Export all data",
which downloads all submissions and all repeating subgroups as CSV spreadsheets,
and all media attachments in a local subfolder "attachments".

An alternative pathway to getting data out of ODK Central is to use the REST API
as documented (with live examples in multiple programming languages) at the
[ODK Central API docs](https://odkcentral.docs.apiary.io/).

# Configure ruODK
The OData service URL is shown in the form's "Submissions" tab >
"Analyse via OData" on ODK Central. It contains base URL, project ID, and 
form ID and is used by `ruODK::ru_setup()`.

```{r ru_setup_demo, eval=FALSE}
# ODK Central's OData URL contains base URL, project ID, and form ID
# ODK Central credentials can live in .Renviron
# See vignette("setup") for setup and authentication options.
ruODK::ru_setup(
  svc = Sys.getenv("ODKC_TEST_SVC"),
  un = Sys.getenv("ODKC_TEST_UN"),
  pw = Sys.getenv("ODKC_TEST_PW"),
  tz = "Australia/Perth",
  verbose = TRUE
)

# File attachment download location
loc <- fs::path("media")
```

This vignette shows how to access data, but under the bonnet uses the 
included package data. This allows to rebuild the vignette offline and without
credentials to the originating ODK Central server.

```{r use_pkg_data}
# Canned data
data("fq_svc")
data("fq_form_schema")
data("fq_meta")
data("fq_raw")
data("fq_raw_strata")
data("fq_raw_taxa")
data("fq_data")
data("fq_data_strata")
data("fq_data_taxa")
```

To extract data from the OData API endpoints, we have to:

* discover data endpoints from the OData service document,
* inspect the metadata schema to infer data types,
* download submissions from the data endpoints,
* download media attachments and adjust their file paths to the downloaded files.

## OData service document
Let's start off with the service document.

```{r load_odata_service, eval=F}
fq_svc <- ruODK::odata_service_get()
```

The same data is included as example data `fq_svc`.

```{r view_service}
fq_svc %>% knitr::kable(.)
```

`ruODK` provides the names and urls of the service endpoints as `tibble`.
We see the main data available under the url `Submissions`, and repeating 
form groups called `taxon_encounter` and `vegetation_stratum` in the ODK form under 
the url `Submissions.taxon_encounter` and `Submissions.vegetation_stratum`, 
respectively.

The main value we get out of the service document are these names of the 
form groups, which can differ between forms.

## OData metadata document
Let's inspect the form metadata to review our data schema.
While we can download the submission data without it, the metadata document 
contains information about field data types and attachment names.

```{r load_metadata, eval=FALSE}
fq_meta <- ruODK::odata_metadata_get()
```

```{r view_metadata, fig.width=7}
if (requireNamespace("listviewer")) {
  listviewer::jsonedit(fq_meta)
} else {
  ru_msg_info("Please install package listviewer!")
}
```

As an alternative to the OData metadata document, ODK Central also offers
form metadata as a much cleaner JSON document, which `ruODK` can read and parse
into a clean `tibble` of field type, name, and path. 

`ruODK` uses this introspection to parse submission data.

```{r form_schema, eval=FALSE}
fq_form_schema <- ruODK::form_schema()
```

```{r form_schema_show}
fq_form_schema %>% knitr::kable(.)
```

## OData submission data documents
Now let's download the form submissions and, separately, repeating form groups.
`ruODK::odata_submission_get()` defaults to download the submission data,
parse it into a tidy tibble, parses dates and datetimes, downloads
and links file attachments, and handles spatial datatypes.

This vignette is built with canned data, so the verbose messages are not shown.

With `wkt=TRUE`, we'll receive spatial types as Well Known Text, which `ruODK` parses
as plain text.
With `wkt=FALSE` (the default), we'll receive spatial types as GeoJSON, which
`ruODK` parses into a nested list.
`ruODK` retains the original spatial field, and annotates the data with 
extracted longitude, latitude, altitude, and (where given) accuracy.
These additional fields are prefixed with the original field name to prevent 
name collisions between possibly multiple location fields.

```{r load_odata, eval=FALSE}
fq_data <- ruODK::odata_submission_get(
  table = fq_svc$name[1], 
  local_dir = loc,
  wkt=TRUE)

fq_data_strata <- ruODK::odata_submission_get(
  table = fq_svc$name[2], 
  local_dir = loc
)

fq_data_taxa <- ruODK::odata_submission_get(
  table = fq_svc$name[3], 
  local_dir = loc,
  wkt=TRUE)
```

# Detour: Data rectangling
<!-- TODO: rewrite -->
<!-- * odata_submission_get(parse=F) -->
<!-- * odata_submission_rectangle(form_schema=fs) -->
<!-- * handle_ru_{datetimes, attachments, geopoints, geotraces, geoshapes} -->
The function `ruODK::odata_submission_get()` received the original XML response 
as a nested list of lists.
To analyse and visualise the data, this nested list of lists must be transformed 
into a rectangular shape.
The function `ruODK::odata_submission_rectangle()` is used internally to 
recursively un-nest list columns using `tidyr::unnest_wider()`.
Unnamed columns, notably the anonymous lat/lon/alt coordinates, are named 
automatically to become unique (a feature of `tidyr::unnest_*()`), and then 
sanitised using the helper `janitor::clean_names()`.

By default, form group names are used as prefix to the field names. 
This behaviour can be disabled by handing the argument `names_sep=NULL` to
`tidyr::unnest_wider()` through running 
`ruODK::odata_submission_get() %>% ruODK::odata_submission_rectangle(names_sep = NULL)`.

The vectorised function `ruODK::attachment_get()` is then used internally 
to download and link attachments like photos and other media to a local, 
relative path. This will take some time during the first run. 
Once the files exist locally, the download will be skipped.

When used through `ruODK::odata_submission_get()`, `ruODK` will introspect the
form schema to detect and then parse media attachment fields automatically.
When used manually, field names of media attachment fields can be (partially or 
fully) specified, see `??ruODK::attachment_get()`.

The date formats are parsed from ISO8601 timestamps into POSIXct objects with
`ruODK::handle_ru_datetimes()`. We use our local timezone (GMT+08) in this example.
`ruODK` introspects the form schema to detect and then parse date and datetime 
fields automatically.

The repeated subgroup `taxon_encounter` is left joined to the main submission data
to receive a (repeated) copy of the main submission data 
(such as location, time and habitat description). 
We will do the same to the other repeated subgroup `vegetation_stratum`.

For clarity, we enable verbose messages from `ruODK::odata_submission_get()` 
and preserve the message output in the code chunk options with `message=TRUE`.
In real-world use cases, messages can be disabled through the chunk option
`message=FALSE`.

We use a custom local path for attachments (`loc`). This results in a smaller
installed package size for `ruODK`, as it shares the attachment files with the
other vignettes. The default is a local folder `media`.

The raw and unparsed example data is provided as data objects 
`fq_raw` (main submissions of form Flora Quadrat 0.4), 
`fq_raw_taxa` (repeated group "Taxon Encounter" within a Flora Quadrat), and
`fq_raw_strata` (repeated group "Vegetation Stratum" within a Flora Quadrat).

The parsed versions are included as data objects `fq_data`, `fq_data_strata`, 
and `fq_data_taxa`. To enable users without ODK Central credentials to build
this vignette (e.g. on package installation with `build_vignettes=TRUE`),
we show the real functions (such as `ruODK::odata_submission_get()`), but do not
evaluate them. Instead, we use "canned data". The `ruODK` test suite ensures 
that canned data are equivalent to live data. 

The result of this code chunk should be exactly the same as the compact version
with `odata_submission_get(parse=TRUE)`.

```{r}
# Candidates for ruODK::handle_ru_datetimes()
fq_form_schema %>% 
  dplyr::filter(type %in% c("dateTime", "date")) %>% 
  knitr::kable(.)

# Candidates for ruODK::handle_ru_attachments()
fq_form_schema %>% 
  dplyr::filter(type == "binary") %>% 
  knitr::kable(.)

# Candidates for ruODK::handle_ru_geopoints()
fq_form_schema %>% 
  dplyr::filter(type == "geopoint") %>% 
  knitr::kable(.)
```


```{r, eval=FALSE}
# The raw submission data
fq_raw <- ruODK::odata_submission_get(table = fq_svc$name[1], parse = FALSE)
fq_strata <- ruODK::odata_submission_get(table = fq_svc$name[2], parse = FALSE)
fq_taxa <- ruODK::odata_submission_get(table = fq_svc$name[3], parse = FALSE)

# Parse main data
fq_data <- fq_raw %>%
  ruODK::odata_submission_rectangle() %>%
  ruODK::handle_ru_datetimes(fq_form_schema) %>%
  ruODK::handle_ru_geopoints(fq_form_schema) %>%
  ruODK::handle_ru_geotraces(fq_form_schema) %>%
  ruODK::handle_ru_geoshapes(fq_form_schema) %>% 
  ruODK::handle_ru_attachments(fq_form_schema, local_dir = t)

# Parse nested group "taxa"
fq_data_taxa <- fq_taxa %>% 
  ruODK::odata_submission_rectangle() %>%
  ruODK::handle_ru_datetimes(fq_form_schema) %>%
  ruODK::handle_ru_geopoints(fq_form_schema) %>%
  ruODK::handle_ru_geotraces(fq_form_schema) %>%
  ruODK::handle_ru_geoshapes(fq_form_schema) %>% 
  ruODK::handle_ru_attachments(fq_form_schema, local_dir = t) %>%
  dplyr::left_join(fq_data, by = c("submissions_id" = "id"))

# Parse nested group "strata"
fq_data_strata <- fq_strata %>% 
  ruODK::odata_submission_rectangle() %>%
  ruODK::handle_ru_datetimes(fq_form_schema) %>%
  ruODK::handle_ru_geopoints(fq_form_schema) %>%
  ruODK::handle_ru_geotraces(fq_form_schema) %>%
  ruODK::handle_ru_geoshapes(fq_form_schema) %>% 
  ruODK::handle_ru_attachments(fq_form_schema, local_dir = t) %>%
  dplyr::left_join(fq_data, by = c("submissions_id" = "id"))
```

Note: A manually resized version of the original photos in this example live in
the package source under [`articles/media`
](https://github.com/ropensci/ruODK/tree/main/vignettes/media). 
To minimise package size, they were resized with imagemagick:

```{sh, eval=FALSE}
find vignettes/media -type f -exec mogrify -resize 200x150 {} \\;
```

## DIY rectangling
For those wishing to go one step further, this section demonstrates the inner
workings of `ruODK`, the recursive use of `tidyr::unnest_wider()`.

The unnesting could also be done manually by building up a pipeline, which
stepwise unnests each list column.
This requires knowledge of the data structure, which can either be looked up
from the metadata, or by inspecting the raw data, `fq_raw`.

The following command has been built by stepwise adding `tidyr::unnest_wider()` 
expressions to the pipe until all list columns were eliminated.

The trailing `invisible()` allows us to toggle parts of the pipe by catching the
dangling `%>%`.

```{r rectangle_diy}
fq_data_diy <- tibble::tibble(value = fq_raw$value) %>%
  tidyr::unnest_wider(value) %>%
  # 1. find list columns:
  tidyr::unnest_wider(`__system`) %>%
  tidyr::unnest_wider(meta) %>%
  # add more lines here to unnest other form groups
  #
  # 2. rename column names
  dplyr::rename(
    uuid = `__id`
    # add more columns, e.g.
    # longitude=`...1`, latitude=`...2`, altitude=`...3`
  ) %>%
  # 3. handle media attachments
  # dplyr::mutate(photo_1 = attachment_get(data_url, uuid, photo_1)) %>%
  invisible()
```

# Visualise data
<!-- TODO: split out and show working examples for each vis (table, map, pivot, export) -->
<!-- for each data type (date, point, trace, shape) -->
This section provides some examples of standard data visualisations.

## Datatable
The package `DT` provides an interactive (and searchable) datatable.

```{r vis_data}
DT::datatable(fq_data)
DT::datatable(fq_data_taxa)
DT::datatable(fq_data_strata)
# DT::datatable(head(fq_data_diy))
```

## Map
The R package `leaflet` provides interactive maps.

Constructing label and popup requires knowledge of the dataset structure.

```{r map_data}
leaflet::leaflet(width = 800, height = 600) %>%
  leaflet::addProviderTiles("OpenStreetMap.Mapnik", group = "Place names") %>%
  leaflet::addProviderTiles("Esri.WorldImagery", group = "Aerial") %>%
  leaflet::clearBounds() %>%
  leaflet::addAwesomeMarkers(
    data = fq_data,
    lng = ~location_corner1_longitude, 
    lat = ~location_corner1_latitude,
    icon = leaflet::makeAwesomeIcon(text = "Q", markerColor = "red"),
    label = ~ glue::glue("{location_area_name} {encounter_start_datetime}"),
    popup = ~ glue::glue(
      "<h3>{location_area_name}</h3>",
      "Survey start {encounter_start_datetime}</br>",
      "Device {device_id}</br>",
      "<h5>Site</h5>",
      '<div><img src="{location_quadrat_photo}"',
      ' height="150px" alt="Quadrat photo"></img></div>',
      "<h5>Mudmap</h5>",
      '<div><img src="{perimeter_mudmap_photo}',
      ' height="150px" alt="Mudmap"></img></div>',
      "<h5>Habitat</h5>",
      "Morphological type: {habitat_morphological_type}</br>",
      '<div><img src="{habitat_morphological_type_photo}"',
      'height="150px" alt="Morphological type"></img></div>'
    ),
    clusterOptions = leaflet::markerClusterOptions()
  ) %>%
  leaflet::addLayersControl(
    baseGroups = c("Place names", "Aerial"),
    options = leaflet::layersControlOptions(collapsed = FALSE)
  )
```

```{r map_data_tae}
leaflet::leaflet(width = 800, height = 600) %>%
  leaflet::addProviderTiles("OpenStreetMap.Mapnik", group = "Place names") %>%
  leaflet::addProviderTiles("Esri.WorldImagery", group = "Aerial") %>%
  leaflet::clearBounds() %>%
  leaflet::addAwesomeMarkers(
    data = fq_data_taxa,
    lng = ~ location_corner1_longitude, 
    lat = ~ location_corner1_latitude,
    icon = leaflet::makeAwesomeIcon(text = "T", markerColor = "green"),
    label = ~ glue::glue("{field_name} {encounter_start_datetime}"),
    popup = ~ glue::glue(
      "<h3>{field_name}</h3>",
      "Survey start {encounter_start_datetime}</br>",
      "Device {device_id}</br>",
      "<h5>Taxon</h5>",
      '<div><img src="media/{photo_in_situ}"',
      ' height="150px" alt="Taxon in situ"></img></div>',
      "Specimen barcode: {voucher_specimen_barcode}</br>",
      "Life form: {life_form}</br>"
    ),
    clusterOptions = leaflet::markerClusterOptions()
  ) %>%
  leaflet::addLayersControl(
    baseGroups = c("Place names", "Aerial"),
    options = leaflet::layersControlOptions(collapsed = FALSE)
  )
```

## Summarising data
See Hadley Wickam's 
[R for Data Science](https://r4ds.had.co.nz/exploratory-data-analysis.html)
for more ideas on data exploration.

```{r eda, fig.height=5, fig.width=7}
# How many submissions per device?
fq_data %>% 
  dplyr::group_by(meta_instance_id) %>% 
  dplyr::tally() %>% 
  knitr::kable()

# How many species sightings per life form?
fq_data_taxa %>% 
  dplyr::group_by(life_form) %>% 
  dplyr::tally() %>% 
  knitr::kable()

# GGplot of a pivot table
fq_data_taxa  %>%
  dplyr::group_by(life_form) %>%
  dplyr::tally() %>%
  ggplot2::ggplot(ggplot2::aes(x = life_form, y = n)) +
  ggplot2::labs(
    title = "Title",
    subtitle = "Subtitle",
    x = "Life form",
    y = "Abundance"
  ) +
  ggplot2::geom_point() +
  ggplot2::theme_classic()

# GGplot with groups
fq_data_taxa %>%
  ggplot2::ggplot(
    ggplot2::aes(
      x=encounter_start_datetime, 
      y=field_name, 
      colour=life_form,
      shape=meta_instance_id)) + 
  ggplot2::labs(
    title="Title", 
    subtitle="Subtitle", 
    x="Observation date", 
    y="Species", 
    colour="Life form",
    shape="Data collection device") +
  ggplot2::geom_point() + 
  ggplot2::theme_classic() + 
  ggplot2::theme(axis.text.x = element_text(angle = 45, hjust = 1))
```

# Export
The rectangled data can now be exported. e.g. to CSV. Note that all list columns 
must be either unnested or dropped before exporting to CSV.

```{r export, eval=F}
fq_data %>% readr::write_csv("flora_quadrats.csv")
fq_data_tae %>% readr::write_csv("flora_quadrats_taxon_encounters.csv")
fq_data_veg %>% readr::write_csv("flora_quadrats_vegetation_strata.csv")
```

# ruReady to ODK?
In this vignette, we took a scenic tour through the general workflow of accessing
and wrangling ODK Central data using `ruODK`.

For your convenience, `ruODK` includes a template RMarkdown workbook with the
essential steps of the above workflow and colour-by-numbers instructions, which
can be used as a starting point for projects using data from ODK Central.

To create a new RMarkdown document from the `ruODK` template, 
run `rmarkdown::draft("test.Rmd", "odata", package="ruODK")`.

Users of RStudio can alternatively "Create a new RMarkdown document" 
"From template" and select `ruODK`'s template "ODK Central via OData".

Make sure to install a fresh version of `ruODK` to get the latest and greatest
template.
---
title: "Accessing the RESTful API"
description: >
  Accessing submission data via the REST pathway.
output: 
  rmarkdown::html_vignette:
    self_contained: false
vignette: >
  %\VignetteIndexEntry{Accessing the RESTful API}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
```

# Scope
This vignette provides a guided walk-through of the "getting data out" 
functions of the RESTful API endpoints which list and view details.

`ruODK` users would mix and match parts of the demonstrated workflows to build
their own data pipelines, e.g.:

* to build a quick analysis from all data, freshly downloaded from a smaller 
  project, or
* to build an interactive ETL pipeline to selectively download only new submissions
  for further processing and upload into downstream data warehouses.
  
A typical and more stream-lined workflow is provided in the RMarkdown template 
"ODK Central via OData" which is supplied by `ruODK`.

## Three ways to happiness

ODK Central offers no less than three different ways to access data:

* viewing ODK Central data in MS PowerBI, MS Excel, Tableau, or `ruODK` 
  through the OData service endpoints, or
* downloading all submissions including attachments as one (possibly gigantic) 
  zip archive either through the "Export Submissions" button in the ODK Central
  form submissions page or through `ruODK`, or
* viewing ODK Central data through `ruODK`'s RESTful API functions.

While the `vignette("odata", package="ruODK")` 
(online [here](https://docs.ropensci.org/ruODK/articles/odata-api.html)) 
illustrates the first option, this vignette demonstrates the remaining two.

Not implemented (yet) are the "managing ODK Central" functions which create, 
update, and delete projects, forms, users, roles, and permissions. 
We haven't yet found a strong use case to automate those functions - 
ODK Central (driven by humans) does those jobs beautifully on an expected scale.

# Setup ruODK
See [`vignette("Setup", package = "ruODK")`](https://docs.ropensci.org/ruODK/articles/setup.html) 
for detailed options to configure `ruODK`.

Here, we'll grab the OData service URL from the form used in this vignette,
plus username and password of a web user of ODK Central with access to that form.

`ruODK::ru_setup()` will populate the default url, project ID, and form ID which
are used by `ruODK`'s other functions (unless specified otherwise).

```{r ru_setup}
library(ruODK)
ruODK::ru_setup(
  svc = "https://sandbox.central.getodk.org/v1/projects/14/forms/build_Flora-Quadrat-0-4_1564384341.svc",
  un = Sys.getenv("ODKC_TEST_UN"),
  pw = Sys.getenv("ODKC_TEST_PW"),
  tz = "Australia/Perth",
  verbose = TRUE
)
t <- fs::dir_create("media")
```

```{r load_canned_data, echo=FALSE}
# We load canned data, so end users can build vignettes without authenticated 
# calls to ODK Central
data("fq_project_list")
data("fq_project_detail")
data("fq_form_list")
data("fq_form_xml")
data("fq_form_schema")
data("fq_zip_data")
data("fq_zip_strata")
data("fq_zip_taxa")
data("fq_submission_list")
data("fq_submissions")
data("fq_attachments")
```

# Projects

List projects. We see the project ID, a name, the number of forms and app users,
dates of last form submissions plus project management timestamps (created, 
updated).

The important bit here is the project ID.

```{r project_list, eval=F}
fq_project_list <- ruODK::project_list()
```

```{r project_list_show}
fq_project_list %>% knitr::kable(.)
```

Inspect a project using its ID. We receive a tibble with exactly one row,
all the columns of `ruODK::project_list` plus a column `verbs`, which contains all
available API actions for a project.

```{r project_details, eval=FALSE}
fq_project_detail <- ruODK::project_detail()
```

```{r project_detail_show}
# Project details (without verbs)
fq_project_detail %>% dplyr::select(-"verbs") %>%  knitr::kable(.)

# Available verbs
fq_project_detail$verbs[[1]] %>% unlist(.)
```

Nothing apart from the verbs is new compared to the data returned by 
`ruODK::project_list`.

To learn more about the functionality behind the verbs, refer to the interactive 
[ODK Central API documentation](https://odkcentral.docs.apiary.io/#reference/project-management).

To retrieve data from ODK Central, the functions shown below will suffice.

# Forms
## List forms for a project

To download form submissions, we need to know project ID and form ID.

There are several ways of retrieving the form ID:

* Browsing forms in the ODK Central's project overviews,
* Stealing the form ID from the OData service endpoint URL as shown on 
  ODK Central's form submission page,
* Listing form metadata for a given project ID with `ruODK::form_list()`.

```{r form_list, eval=FALSE}
fq_form_list <- ruODK::form_list()
```

```{r form_list_show}
fq_form_list %>% knitr::kable(.)
```

Further to the metadata shown here, a column `xml` contains the entire XForms
definition (originally XML) as nested list.

If the original XML is needed rather than the R equivalent (nested list), 
we can use `ruODK::form_xml` with parameter `parse=FALSE`:

```{r form_xml, eval=F}
fq_form_xml <- ruODK::form_xml(parse=FALSE)
```

```{r form_xml_show}
if (require(listviewer)){
  listviewer::jsonedit(fq_form_xml)
} else {
 ru_msg_warn("Install package listviewer to browse the form XML.")  
}
```

## Inspect form schema

The `form_schema` represents all form fields of the XForms definition.

See the 
[ODK Central API docs](https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-individual-form/getting-form-schema-fields) 
and the examples of `??ruODK::form_schema()` for more detail.


```{r form_schema, eval=FALSE}
fq_form_schema <- ruODK::form_schema()
```

```{r form_schema_view}
fq_form_schema %>% knitr::kable(.)
```

## Show details of one form

The details of a form are exactly the same as the output of `ruODK::form_list()`.

```{r form_detail, eval=F}
fq_form_detail <- ruODK::form_detail()
```

```{r form_detail_view}
fq_form_detail %>% knitr::kable(.)
```

# Submissions
We are getting closer to the actual data! This section shows two of the options
for data access: dump all submissions, or extract a subset.

## Get all submissions for one form
Smaller datasets lend themselves to be exported in one go.
ODK Central offers one giant zip file containing all submissions, any
repeating groups, and any attachments both on the form submission page, and as
API endpoint which is provided as `ruODK::submission_export()`.

The default behaviour of `ruODK::submission_export()` is to write the zip file
to the project root (`here::here()`), and to overwrite existing previous downloads.
See `?ruODK::submission_export()` for alternative download and retention options.

In the following chuck, we illustrate common tasks:

* Download the zip file.
* Unpack the zip file.
* Join repeating form group data `data_taxon` to main data `data_quadrat` to
  annotate `data_taxon` with data from `data_quadrat`.
* Sanitise the column names.
* Prepend all attachment filenames (e.g. `data_quadrat$location_quadrat_photo`,
  `data_taxon$photo_in_situ`) with `media/`.


```{r submission_export, eval=F}
# Predict filenames (with knowledge of form)
fid <- ruODK::get_test_fid()
fid_csv <- fs::path(t, glue::glue("{fid}.csv"))
fid_csv_veg <- fs::path(t, glue::glue("{fid}-vegetation_stratum.csv"))
fid_csv_tae <- fs::path(t, glue::glue("{fid}-taxon_encounter.csv"))

# Download the zip file
se <- ruODK::submission_export(local_dir = t, overwrite = FALSE, verbose = TRUE)

# Unpack the zip file
f <- unzip(se, exdir = t)
fs::dir_ls(t)

# Prepend attachments with media/ to turn into relative file paths
fq_zip_data <- fid_csv %>% 
  readr::read_csv(na = c("", "NA", "na")) %>% # form uses "na" for NA
  janitor::clean_names(.) %>% 
  dplyr::mutate(id = meta_instance_id) %>% 
  ruODK::handle_ru_datetimes(fq_form_schema) %>% 
  ruODK::handle_ru_geopoints(fq_form_schema) %>% 
  ruODK::attachment_link(fq_form_schema)

fq_zip_strata <- fid_csv_veg %>% 
  readr::read_csv(na = c("", "NA", "na")) %>%
  janitor::clean_names(.) %>% 
  dplyr::mutate(id = parent_key) %>% 
  # ruODK::handle_ru_datetimes(fq_form_schema) parent_key%>% # no dates
  # ruODK::handle_ru_geopoints(fq_form_schema) %>%  # no geopoints
  # ruODK::ruODK::attachment_link(fq_form_schema) %>% # no att.
  dplyr::left_join(fq_zip_data, by = c("parent_key" = "meta_instance_id"))

fq_zip_taxa <- fid_csv_tae %>%
  readr::read_csv(na = c("", "NA", "na")) %>%
  janitor::clean_names(.) %>% 
  dplyr::mutate(id = parent_key) %>% 
  # ruODK::handle_ru_datetimes(fq_form_schema) %>% 
  # ruODK::handle_ru_geopoints(fq_form_schema) %>% 
  # ruODK::ruODK::attachment_link(fq_form_schema) %>%
  dplyr::left_join(fq_zip_data, by = c("parent_key" = "meta_instance_id"))
```

```{r zip_view}
head(fq_zip_data)
head(fq_zip_strata)
head(fq_zip_taxa)
# Further: create map with popups, see vignette "odata"
```

## List submissions for one form
Not always is it appropriate to download all submissions and all attachments
at once. 

If forms feed into downstream data warehouses, the typical ETL workflow is to 

* List all submissions from ODK Central
* Select the subset of new submissions to download, e.g.
  * Submissions younger than the oldest submission date in the data warehouse.
  * Submissions whose `instance_id` is not already present in the data warehouse.
* Download only the selected submissions.
* Download attachments of only the selected submissions.

```{r submission_list, eval=F}
fq_submission_list <- ruODK::submission_list()
```

```{r submission_list_view}
fq_submission_list %>% knitr::kable(.)
```

The list of submissions critically contains each submission's unique ID in
`instance_id`. If the submissions shall be downloaded and uploaded into another
data warehouse, the `instance_id` can be used to determine whether a record
already exists in the downstream warehouse or not.
This workflow is preferable where the majority of submissions is already 
imported into another downstream data warehouse, and we only want to add new 
submissions, as in submissions which are not already imported into the data 
warehouse.

Furthermore, the `instance_id`s can now be used to retrieve the actual 
submissions.

## Get submission data

In order to import each submission, we need to retrieve the data by 
`instance_id`.

```{r submission_data, eval=F}
# One submission
fq_one_submission <- ruODK::get_one_submission(fq_submission_list$instance_id[[1]])

# Multiple submissions
fq_submissions <- ruODK::submission_get(fq_submission_list$instance_id)
```

## Parse submissions
The data in `sub` is one row of the bulk downloaded submissions in `data_quadrat`.
The data in `submissions` represents all (or let's pretend, the selected) 
submissions in `data_quadrat`.
The field `xml` contains the actual submission data including repeating groups.

The structure is different to the output of `ruODK::odata_submission_get`,
therefore `ruODK::odata_submission_rectangle` does not work for those, as
here we might have repeating groups included in a submission.

This structure could be used for upload into data warehouses accepting nested data
as e.g. JSON.

```{r view_submission_data, fig.width=7}
if (requireNamespace("listviewer")) {
  listviewer::jsonedit(fq_submissions, mode = "code")
} else {
  ru_msg_info("Please install package listviewer!")
}
```

# Outlook
The approach shown here yields nested and stand-alone records, which is useful
if the subsequent use requires records in nested JSON or XML format. 
Complex forms with repeating sub-groups will result in highly nested lists, 
whose structure heavily depends on the completeness of the submissions.

The other approach shown in 
[`vignette("odata-api", package="ruODK")`](https://docs.ropensci.org/ruODK/articles/odata-api.html) 
yields rectangled data in several normalised tables, which is useful for 
analysis and visualisation.
---
title: "Setup"
description: >
  Provide sensitive credentials to ruODK through different pathways.
output: 
  rmarkdown::html_vignette:
    self_contained: false
vignette: >
  %\VignetteIndexEntry{Setup}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}  
---

```{r, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
options(rmarkdown.html_vignette.check_title = FALSE)
```

```{r setup}
library(ruODK)
```

# Configure ruODK

`ruODK` functions work on a given ODK Central instance using a web user's
credentials (username and password). Some functions also require project and
form IDs.

`ruODK`'s functions accept these parameters either as explicit keyword
arguments, or fall back to defaults.

Note: Always consider code to be public. Never use plain text credentials in
code.

`ruODK` suggests as best practice to set the defaults (or parts thereof) using
`ruODK::ru_setup()` or through permanently setting environment variables.

## Best practice: `ru_setup`

`ruODK` provides helpers for settings, `ru_setup()` to set and `ru_settings()`
to get settings.

While normal users of `ruODK` will only need a default pid, fid, url, username,
and password, contributors to `ruODK` can include optional test server settings
which are required to run the test suite and build the vignettes.

`ruODK` infers the base URL, project and form ID from the OData Service URL
which is shown in ODK Central on the Form Submissions tab.

Unless specified as function parameter, `ruODK` converts dates and times to the
default timezone, in this example set to "Australia/Perth".

Furthermore, some functions offer verbose messages, which can assist to debug
unexpected behaviour. Unless specified in the settings, or in the respective
function calls, `ruODK` defaults to hide verbose messages.

```{r ru_setup, eval=FALSE}
# ruODK user using OData service URL, username (an email), and password
# Never use plaintext username and password, use Sys.getenv() instead
ruODK::ru_setup(
  svc = Sys.getenv("ODKC_SVC"), 
  un = Sys.getenv("ODKC_UN"), 
  pw = Sys.getenv("ODKC_PW"),
  tz = "Australia/Perth",
  verbose = TRUE
)

# ruODK contributors: see contributing guidelines for .Renviron variables

# Review settings
ruODK::ru_settings()
```

Now we can call `ruODK` functions without specifying `url`, `un`, and `pw`, and
let `ruODK` fall back to the defaults:

```{r, eval=FALSE}
ruODK::project_list()
ruODK::submission_list()
```

## Permanent defaults: Environment variables

Read a great overview of R's startup process, and how environment variables are
sourced at the beginning of a new session at
[whattheyforgot.org](https://whattheyforgot.org/r-startup.html).

`ruODK`'s functions default to the getters
`ruODK::get_default_{pid,fid,url,un.pw}()`. These getters in turn look up their
values from environment variables.

The getters and setters are documented in the "Settings" family of the 
[ruODK function reference](https://docs.ropensci.org/ruODK/reference/index.html).

A convenient way to have often used environment variables available is to add
them to `~/.Renviron` using `usethis::edit_r_environ(scope = "user")`. This
loads them freshly into a new session, eliminating the need to run `ru_setup()`.
Note that the environment variables can be cleared or overwritten through
calling `ru_setup()` or `Sys.setenv()` with respective arguments.

`ru_setup()` will not change any omitted arguments.

```{r open_renviron, eval=FALSE}
usethis::edit_r_environ(scope = "user")
```

```{r renviron, eval=FALSE}
ODKC_PID=1
ODKC_FID="build_Flora-Quadrat-0-2_1558575936"
ODKC_URL="https://odkcentral.dbca.wa.gov.au"
ODKC_UN="me@email.com"
ODKC_PW="..."
ODKC_PP="..."
ODKC_VERSION=0.8

# Test settings
ODKC_TEST_SVC="https://odkc.dbca.wa.gov.au/v1/projects/2/forms/Flora-Quadrat-04.svc"
ODKC_TEST_URL="https://odkc.dbca.wa.gov.au"
ODKC_TEST_PID=2
ODKC_TEST_PID_ENC=3
ODKC_TEST_PP="ThePassphrase"
ODKC_TEST_FID="Flora-Quadrat-04"
ODKC_TEST_FID_ZIP="Spotlighting-06"
ODKC_TEST_FID_ATT="Flora-Quadrat-04-att"
ODKC_TEST_FID_GAP="Flora-Quadrat-04-gap"
ODKC_TEST_FID_WKT="Locations"
ODKC_TEST_FID_I8N0="I8n_no_lang"
ODKC_TEST_FID_I8N1="I8n_label_lng"
ODKC_TEST_FID_I8N2="I8n_label_choices"
ODKC_TEST_FID_I8N3="I8n_no_lang_choicefilter"
ODKC_TEST_FID_I8N4="I8n_lang_choicefilter"
ODKC_TEST_FID_ENC="Locations"
ODKC_TEST_VERSION=1.0
RU_VERBOSE=TRUE
RU_TIMEZONE="Australia/Perth"
RU_RETRIES=3
ODKC_TEST_UN="..."
ODKC_TEST_PW="..."
```

As an alternative to setting environment variables through `~/.Renviron`, you
can set them through `Sys.setenv()`:

```{r setenv, eval=FALSE}
Sys.setenv(ODKC_URL="https://odkc.dbca.wa.gov.au")
Sys.setenv(ODKC_UN="me@mail.com")
Sys.setenv(ODKC_PW="...")
Sys.setenv(ODKC_PP="...")
Sys.setenv(ODKC_TEST_URL="...")
Sys.setenv(ODKC_TEST_UN="...")
Sys.setenv(ODKC_TEST_PW="...")
Sys.setenv(ODKC_TEST_PID=2)
# continue for the remaining test variables shown above
```

## The hard way: Per function call

We can supply those credentials to each `ruODK` function call.

Note, including sensitive credentials in plain text into code is bad practice.
This option is shown only for completeness.

```{r, eval=FALSE}
# Use ruODK without ru_setup
ruODK::project_list(
  url="https://my-odkc.com", 
  un = Sys.getenv("ODKC_UN"), 
  pw = Sys.getenv("ODKC_PW")
)

# Tests use default test settings explicitly
ruODK::project_list(
  url=ruODK::get_test_url(),
  un=ruODK::get_test_un(), 
  pw=ruODK::get_test_pw()
)
```

An example use case are the `ruODK` tests, which explicitly set `url`, `un`,
`pw`, `pp`, `pid` and `fid` from the test variables
`ruODK::get_test_{url, un, pw, pp, pid, fid}()`. Note that this uses functions
instead of plain text versions of sensitive credentials. Alternatively,
variables could also be used to set credentials per function call.

## Moving across forms

`ruODK`'s functions default to the default values for project ID (`pid`), form
ID (`fid`), base URL (`url`), username (`un`), and password (`pw`).

A typical workflow is to run several functions of `ruODK` against one form (and
the overarching project). By running `ru_setup()` with the form's OData service
URL and a web user's username and password, subsequent functions can omit `pid`,
`fid`, `url`, `un`, and `pw`.

```{r, eval=FALSE}
# Server 1, Project 1, Form 1
ruODK::ru_setup(
  svc = "https://central1.org/v1/projects/1/forms/form1.svc",
  un = Sys.getenv("ODKC_UN"), 
  pw = Sys.getenv("ODKC_PW")
)

ruODK::project_detail()
ruODK::form_detail()
ruODK::submission_list()

# Server 1, Project 1, Form 3
ruODK::ru_setup(svc = "https://central1.org/v1/projects/1/forms/form3.svc")

ruODK::project_detail()
ruODK::form_detail()
ruODK::submission_list()

# Server 1, Project 5, Form 4
ruODK::ru_setup(svc = "https://central1.org/v1/projects/5/forms/form4.svc")

ruODK::project_detail()
ruODK::form_detail()
ruODK::submission_list()

# Server 2, Project 11, Form 12
ruODK::ru_setup(
  svc = "https://central2.org/v1/projects/11/forms/form12.svc",
  un = Sys.getenv("ODKC_UN2"), 
  pw = Sys.getenv("ODKC_PW2")
)

ruODK::project_detail()
ruODK::form_detail()
ruODK::submission_list()

# Server 2, Project 11, Form 13
ruODK::ru_setup(svc = "https://central2.org/v1/projects/11/forms/form13.svc")

ruODK::project_detail()
ruODK::form_detail()
ruODK::submission_list()
```

## Legacy support for deprecated ODK Central features

Users of the latest ODK Central version can safely ignore this section. Users of
older ODK Central versions will be exposed to breaking changes coming from
upstream ODK Central.

Very occasionally, ODK Central will deprecate API endpoints. `ruODK` aims to
keep up with the latest ODK Central API, and will adapt its behaviour to support
both new and old, deprecated behaviour.

As there is no corresponding API endpoint to determine the ODK Central version,
`ruODK` has introduced the environment variables `ODKC_VERSION` and
`ODKC_TEST_VERSION`. If unset, `ruODK` defaults to the latest ODK Central
version, which is the version deployed to the ODK Central sandbox (on which
`ruODK`'s package tests run).

The numeric format of `major.minor` (e.g. `0.7` or `0.8`) future-proofs `ruODK`
to be able to handle any other future breaking changes gracefully through the
same variable.

A recent example is the change of the API endpoint for `form_schema`, which
replaced the `.schema.json` API in favour of the `/fields` API. This required
`ruODK` to toggle its behaviour between the relatively involved un-nesting of
the JSON schema (versions before 0.8) with parsing a clean, flat list of field
names, types, and XForms paths (version 0.8 onward).
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_ru_geoshapes.R
\name{handle_ru_geoshapes}
\alias{handle_ru_geoshapes}
\title{Split all geoshapes of a submission tibble into their components.}
\usage{
handle_ru_geoshapes(
  data,
  form_schema,
  wkt = FALSE,
  odkc_version = get_default_odkc_version(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{data}{Submissions rectangled into a tibble. E.g. the output of\preformatted{ruODK::odata_submission_get(parse = FALSE) \%>\%
ruODK::odata_submission_rectangle(form_schema = ...)
}}

\item{form_schema}{The \code{form_schema} for the submissions.
E.g. the output of \code{ruODK::form_schema()}.}

\item{wkt}{Whether geofields are GeoJSON (if FALSE) or WKT
strings (if TRUE), default: FALSE.}

\item{odkc_version}{The ODK Central version as decimal number (major.minor).
\code{ruODK} uses this parameter to adjust for breaking changes in ODK Central.

Default: \code{\link{get_default_odkc_version}} or 1.1 if unset.

Set default \code{get_default_odkc_version} through
\code{ru_setup(odkc_version=1.1)}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The submissions tibble with all geoshapes retained in their original
format, plus columns of their first point's coordinate components as
provided by \code{\link{split_geoshape}}.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
For a given tibble of submissions, find all columns which are listed
in the form schema as type \code{geoshape}, and extract their components.
Extracted components are longitude (X), latitude (Y), altitude (Z, where
given), and accuracy (M, where given) of the first point of the geoshape.

The original column is retained to allow parsing into other spatially
enabled formats.
}
\examples{
\dontrun{
library(magrittr)
data("geo_fs")
data("geo_wkt_raw")
data("geo_gj_raw")

# GeoJSON
geo_gj_parsed <- geo_gj_raw \%>\%
  ruODK::odata_submission_rectangle(form_schema = geo_fs) \%>\%
  ruODK::handle_ru_geoshapes(form_schema = geo_fs, wkt = FALSE)

dplyr::glimpse(geo_gj_parsed)

# WKT
geo_wkt_parsed <- geo_wkt_raw \%>\%
  ruODK::odata_submission_rectangle(form_schema = geo_fs) \%>\%
  ruODK::handle_ru_geoshapes(form_schema = geo_fs, wkt = TRUE)

dplyr::glimpse(geo_wkt_parsed)
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_meta}
\alias{fq_meta}
\title{OData metadata document for an ODK Central form.}
\format{
A list of lists
}
\source{
See \code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
}
\usage{
fq_meta
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The OData response for the metadata of an ODK Central form.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_geoshape.R
\name{split_geoshape}
\alias{split_geoshape}
\title{Annotate a dataframe containing a geoshape column with lon, lat, alt of the
geotrace's first point.}
\usage{
split_geoshape(
  data,
  colname,
  wkt = FALSE,
  odkc_version = get_default_odkc_version()
)
}
\arguments{
\item{data}{(dataframe) A dataframe with a geoshape column.}

\item{colname}{(chr) The name of the geoshape column.
This column will be retained.}

\item{wkt}{Whether geofields are GeoJSON (if FALSE) or WKT
strings (if TRUE), default: FALSE.}

\item{odkc_version}{The ODK Central version as decimal number (major.minor).
\code{ruODK} uses this parameter to adjust for breaking changes in ODK Central.

Default: \code{\link{get_default_odkc_version}} or 1.1 if unset.

Set default \code{get_default_odkc_version} through
\code{ru_setup(odkc_version=1.1)}.

See \code{vignette("Setup", package = "ruODK")}.}
}
\value{
The given dataframe with the geoshape column colname, plus
three new columns, \code{colname_longitude}, \code{colname_latitude},
\code{colname_altitude}.
The three new columns are prefixed with the original \code{colname} to
avoid naming conflicts with any other geoshape columns.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This function is used by \code{\link{handle_ru_geopoints}}
on all \code{geopoint} fields as per \code{\link{form_schema}}.
}
\examples{
\dontrun{
library(magrittr)
data("geo_fs")
data("geo_wkt_raw")
data("geo_gj_raw")

# Find variable names of geoshapes
geo_fields <- geo_fs \%>\%
  dplyr::filter(type == "geoshape") \%>\%
  magrittr::extract2("ruodk_name")
geo_fields[1] # First geoshape in data: shape_location_shape_gps

# Rectangle but don't parse submission data (GeoJSON and WKT)
geo_gj_rt <- geo_gj_raw \%>\%
  odata_submission_rectangle(form_schema = geo_fs)
geo_wkt_rt <- geo_wkt_raw \%>\%
  odata_submission_rectangle(form_schema = geo_fs)

# Data with first geoshape split
gj_first_gt <- split_geoshape(geo_gj_rt, geo_fields[1], wkt = FALSE)
cn_gj <- names(gj_first_gt)
testthat::expect_true("shape_location_shape_gps_longitude" \%in\% cn_gj)
testthat::expect_true("shape_location_shape_gps_latitude" \%in\% cn_gj)
testthat::expect_true("shape_location_shape_gps_altitude" \%in\% cn_gj)
testthat::expect_true(
  is.numeric(gj_first_gt$shape_location_shape_gps_longitude)
)
testthat::expect_true(
  is.numeric(gj_first_gt$shape_location_shape_gps_latitude)
)
testthat::expect_true(
  is.numeric(gj_first_gt$shape_location_shape_gps_altitude)
)

wkt_first_gt <- split_geoshape(geo_wkt_rt, geo_fields[1], wkt = TRUE)
cn_wkt <- names(wkt_first_gt)
testthat::expect_true("shape_location_shape_gps_longitude" \%in\% cn_wkt)
testthat::expect_true("shape_location_shape_gps_latitude" \%in\% cn_wkt)
testthat::expect_true("shape_location_shape_gps_altitude" \%in\% cn_wkt)
testthat::expect_true(
  is.numeric(wkt_first_gt$shape_location_shape_gps_longitude)
)
testthat::expect_true(
  is.numeric(wkt_first_gt$shape_location_shape_gps_latitude)
)
testthat::expect_true(
  is.numeric(wkt_first_gt$shape_location_shape_gps_altitude)
)
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_msg.R
\name{ru_msg_info}
\alias{ru_msg_info}
\title{Print a blue info message with an info symbol.}
\usage{
ru_msg_info(message, verbose = get_ru_verbose())
}
\arguments{
\item{message}{(chr) A message to print}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
ru_msg_info("This is an info message.")
}
\seealso{
Other messaging: 
\code{\link{ru_msg_abort}()},
\code{\link{ru_msg_noop}()},
\code{\link{ru_msg_success}()},
\code{\link{ru_msg_warn}()}
}
\concept{messaging}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_geotrace.R
\name{split_geotrace}
\alias{split_geotrace}
\title{Annotate a dataframe containing a geotrace column with lon, lat, alt of the
geotrace's first point.}
\usage{
split_geotrace(
  data,
  colname,
  wkt = FALSE,
  odkc_version = get_default_odkc_version()
)
}
\arguments{
\item{data}{(dataframe) A dataframe with a geotrace column.}

\item{colname}{(chr) The name of the geotrace column.
This column will be retained.}

\item{wkt}{Whether geofields are GeoJSON (if FALSE) or WKT
strings (if TRUE), default: FALSE.}

\item{odkc_version}{The ODK Central version as decimal number (major.minor).
\code{ruODK} uses this parameter to adjust for breaking changes in ODK Central.

Default: \code{\link{get_default_odkc_version}} or 1.1 if unset.

Set default \code{get_default_odkc_version} through
\code{ru_setup(odkc_version=1.1)}.

See \code{vignette("Setup", package = "ruODK")}.}
}
\value{
The given dataframe with the geotrace column colname, plus
three new columns, \code{colname_longitude}, \code{colname_latitude},
\code{colname_altitude}.
The three new columns are prefixed with the original \code{colname} to
avoid naming conflicts with any other geotrace columns.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This function is used by \code{\link{handle_ru_geopoints}}
on all \code{geopoint} fields as per \code{\link{form_schema}}.

The format of the geotrace (GeoJSON, WKT, ODK Linestring) is determined via
parameters \code{wkt} and \code{odkc_version}, rather than inferred from the class of
the column. ODK Linestrings are character vectors without a leading
"LINESTRING (", WKT are character vectors with a leading "LINESTRING (",
and GeoJSON are list columns.
}
\examples{
\dontrun{
library(magrittr)
data("geo_fs")
data("geo_wkt_raw")
data("geo_gj_raw")

# Find variable names of geotraces
geo_fields <- geo_fs \%>\%
  dplyr::filter(type == "geotrace") \%>\%
  magrittr::extract2("ruodk_name")
geo_fields[1] # First geotrace in data: path_location_path_gps

# Rectangle but don't parse submission data (GeoJSON and WKT)
geo_gj_rt <- geo_gj_raw \%>\%
  odata_submission_rectangle(form_schema = geo_fs)
geo_wkt_rt <- geo_wkt_raw \%>\%
  odata_submission_rectangle(form_schema = geo_fs)

# Data with first geotrace split
gj_first_gt <- split_geotrace(geo_gj_rt, geo_fields[1], wkt = FALSE)
testthat::expect_true(
  "path_location_path_gps_longitude" \%in\% names(gj_first_gt)
)
testthat::expect_true(
  "path_location_path_gps_latitude" \%in\% names(gj_first_gt)
)
testthat::expect_true(
  "path_location_path_gps_altitude" \%in\% names(gj_first_gt)
)
testthat::expect_true(
  is.numeric(gj_first_gt$path_location_path_gps_longitude)
)
testthat::expect_true(
  is.numeric(gj_first_gt$path_location_path_gps_latitude)
)
testthat::expect_true(
  is.numeric(gj_first_gt$path_location_path_gps_altitude)
)

wkt_first_gt <- split_geotrace(geo_wkt_rt, geo_fields[1], wkt = TRUE)
testthat::expect_true(
  "path_location_path_gps_longitude" \%in\% names(wkt_first_gt)
)
testthat::expect_true(
  "path_location_path_gps_latitude" \%in\% names(wkt_first_gt)
)
testthat::expect_true(
  "path_location_path_gps_altitude" \%in\% names(wkt_first_gt)
)
testthat::expect_true(
  is.numeric(wkt_first_gt$path_location_path_gps_longitude)
)
testthat::expect_true(
  is.numeric(wkt_first_gt$path_location_path_gps_latitude)
)
testthat::expect_true(
  is.numeric(wkt_first_gt$path_location_path_gps_altitude)
)
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/form_xml.R
\name{form_xml}
\alias{form_xml}
\title{Show the XML representation of one form as list.}
\usage{
form_xml(
  parse = TRUE,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{parse}{Whether to parse the XML into a nested list, default: TRUE}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
The form XML as a nested list.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# With explicit pid and fid
fxml_defaults <- form_xml(1, "build_xformsId")

# With defaults
fxml <- form_xml()
listviewer::jsonedit(fxml)

# form_xml returns a nested list
class(fxml)
# > "list"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-individual-form/retrieving-form-xml}

Other form-management: 
\code{\link{form_detail}()},
\code{\link{form_list}()},
\code{\link{form_schema_ext}()},
\code{\link{form_schema}()}
}
\concept{form-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_msg.R
\name{ru_msg_abort}
\alias{ru_msg_abort}
\title{rlang::abort() with a red error message with a cross symbol.}
\usage{
ru_msg_abort(message)
}
\arguments{
\item{message}{(chr) A message to print}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
ru_msg_abort("This is an error, abort.")
}
}
\seealso{
Other messaging: 
\code{\link{ru_msg_info}()},
\code{\link{ru_msg_noop}()},
\code{\link{ru_msg_success}()},
\code{\link{ru_msg_warn}()}
}
\concept{messaging}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_zip_data}
\alias{fq_zip_data}
\title{A tibble of the main data table of records from a test form.}
\format{
A tibble of main records from a test form.
}
\source{
\code{\link{submission_export}}
run on the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}.
}
\usage{
fq_zip_data
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_setup.R
\name{odata_svc_parse}
\alias{odata_svc_parse}
\title{Retrieve URL, project ID, and form ID from an ODK Central OData service URL.}
\usage{
odata_svc_parse(svc)
}
\arguments{
\item{svc}{(character) The OData service URL of a form as provided by the
ODK Central form submissions tab.
Example: "https://sandbox.central.getodk.org/v1/projects/14/forms/build_Flora-Quadrat-0-2_1558575936.svc"}
}
\value{
A named list with three components (all of type character):
\itemize{
\item \code{url} The ODK Central base URL.
\item \code{pid} The project ID.
\item\code{fid} The form ID.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other ru_settings: 
\code{\link{ru_settings}()},
\code{\link{ru_setup}()},
\code{\link{yell_if_error}()},
\code{\link{yell_if_missing}()}
}
\concept{ru_settings}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{geo_gj}
\alias{geo_gj}
\title{The parsed submissions of a form containing geofields in GeoJSON.}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 1 rows and 52 columns.
}
\source{
\code{\link{odata_submission_get}(wkt=FALSE, parse=TRUE)}
run on the test form
\code{system.file("extdata", "Locations.xml", package = "ruODK")}.
}
\usage{
geo_gj
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_raw_taxa}
\alias{fq_raw_taxa}
\title{OData submission data for a subgroup of an ODK Central form.}
\format{
A list of lists
}
\source{
See \code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
}
\usage{
fq_raw_taxa
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The OData response for a subgroup of an ODK Central form.

This subgroup represents an individual plant taxon which is encountered by
the enumerators. Typically, one voucher specimen is taken for each distinct
encountered plant taxon. A field name is allocated by the enumerators, which
can be the proper canonical name (if known) or any other moniker.
The voucher specimens are later determined by taxonomic experts, who then
provide the real, terminal taxonomic name for a given voucher specimen.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_setup.R
\name{yell_if_missing}
\alias{yell_if_missing}
\title{Abort on missing ODK Central credentials (url, username, password).}
\usage{
yell_if_missing(url, un, pw, pid = NULL, fid = NULL, iid = NULL)
}
\arguments{
\item{url}{A URL (character)}

\item{un}{A username (character)}

\item{pw}{A password (character)}

\item{pid}{A project ID (numeric, optional)}

\item{fid}{A form ID (character, optional)}

\item{iid}{A submission instance ID (character, optional)}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This is a helper function to pat down \code{\link{ruODK}} functions
for missing credentials and stop with a loud but informative yell.
}
\examples{
testthat::expect_error(yell_if_missing("", "username", "password"))
testthat::expect_error(yell_if_missing("url", "", "password"))
testthat::expect_error(yell_if_missing("url", "username", ""))
testthat::expect_error(yell_if_missing(NULL, "", ""))
testthat::expect_error(yell_if_missing("", "", ""))
testthat::expect_error(yell_if_missing("", "", "", ""))
testthat::expect_error(yell_if_missing("", "", "", "", ""))
testthat::expect_error(yell_if_missing("", "", "", "", "", ""))
}
\seealso{
Other ru_settings: 
\code{\link{odata_svc_parse}()},
\code{\link{ru_settings}()},
\code{\link{ru_setup}()},
\code{\link{yell_if_error}()}
}
\concept{ru_settings}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_form_schema}
\alias{fq_form_schema}
\title{JSON form schema for an ODK Central form.}
\format{
The output of \code{ruODK::form_schema()}, a tibble with columns "type",
"name" and "path" and one row per form field.
}
\source{
See \code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
and \code{ruODK::form_schema()}.
}
\usage{
fq_form_schema
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The parsed form schema of an ODK Central form.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.

This data is used to build vignettes offline and without the need for
credentials to an ODK Central server. The test suite ensures that the
"canned" data is identical to the "live" data.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachment_list.R
\name{attachment_list}
\alias{attachment_list}
\title{List all attachments for a list of submission instances.}
\usage{
attachment_list(
  iid,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{iid}{A list of submission instance IDs, e.g. from
\code{\link{submission_list}$instance_id}.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A tibble containing some high-level details of the submission
attachments.
One row per submission attachment, columns are submission attributes:\preformatted{    * name: The attachment filename, e.g. 12345.jpg
    * exists: Whether the attachment for that submission exists on the
      server.
}
}
\description{
List all attachments for a list of submission instances.
}
\examples{
\dontrun{
# Step 1: Setup ruODK with OData Service URL (has url, pid, fid)
ruODK::ru_setup(svc = "...")

# Step 2: List all submissions of form
sl <- submission_list()

# Step 3a: Get attachment list for first submission
al <- get_one_submission_attachment_list(sl$instance_id[[1]])

# Ste 3b: Get all attachments for all submissions
all <- attachment_list(sl$instance_id)
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/attachments/listing-expected-submission-attachments}

\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-form-attachments/listing-expected-form-attachments}

Other submission-management: 
\code{\link{encryption_key_list}()},
\code{\link{submission_detail}()},
\code{\link{submission_export}()},
\code{\link{submission_get}()},
\code{\link{submission_list}()}
}
\concept{submission-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/form_schema_ext.R
\name{form_schema_ext}
\alias{form_schema_ext}
\title{Show the extended schema of one form.}
\usage{
form_schema_ext(
  flatten = FALSE,
  odata = FALSE,
  parse = TRUE,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  odkc_version = get_default_odkc_version(),
  retries = get_retries(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{flatten}{Whether to flatten the resulting list of lists (\code{TRUE})
or not (\code{FALSE}, default). Only applies to ODK Central version < 0.8.}

\item{odata}{Whether to sanitise the field names to match the way they will
be outputted for OData. While the original field names as given in the
XForms definition may be used as-is for CSV output, OData has some
restrictions related to the domain-qualified identifier syntax it uses.
Only applies to ODK Central version < 0.8.
Default: \code{FALSE}.}

\item{parse}{Whether to parse the form schema into a tibble of form field
type and name. This uses \code{\link{form_schema_parse}} internally.
If used together with \code{flatten=TRUE}, \code{\link{form_schema}} will raise
a warning and return the unparsed, flattened form schema.
Only applies to ODK Central version < 0.8.
Default: TRUE.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{odkc_version}{The ODK Central version as decimal number (major.minor).
\code{ruODK} uses this parameter to adjust for breaking changes in ODK Central.

Default: \code{\link{get_default_odkc_version}} or 1.1 if unset.

Set default \code{get_default_odkc_version} through
\code{ru_setup(odkc_version=1.1)}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
A tibble containing the form definition.
For ODK Central 0.8, and with default parameters
(\code{parse=TRUE}) for ODK Central 0.7, \code{\link{form_schema}} returns
a tibble with the columns:

\itemize{
\item \code{name} The field name as given in the form schema.
\item \code{type} The field type, e.g. "string", "select1", etc.
\item \code{path} The XForms path of the field,
\item \code{ruodk_name} The predicted field name as generated by
\code{\link{odata_submission_get}}, prefixed by the path, additionally
cleaned with \code{\link[janitor]{make_clean_names}} to match the
cleaned column names from \code{\link{odata_submission_rectangle}}.
\item \code{label} The field label as given in the form schema.
If specific languages are available,
this column will return the \code{default} language or it will be empty
if this is not specified.
\item \code{label_lang} The field label in languange \emph{_lang} as
given in the form schema.
\item \code{choices} A list of lists containing at least \code{values} and,
if available, \code{labels} of the choices as given in the form schema.
If specific languages are available, this column will return the
\code{default} language or it will be empty if this is not specified.
Please notice that whenever choice filters are applied, this will return
the unfiltered choice list.
\item \code{choices_lang} A list of lists containing at least
\code{values} and, if available, \code{labels} of the choices in language
\emph{_lang} as given in the form schema.
Please notice that whenever choice filters are applied, this will return
the unfiltered choice list.

}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
}
\details{
ODK Central has introduced a new API endpoint in version 0.8 which
returns a parsed and flattened list of fields. This replaces the nested
form schema which is challenging to parse. This list is returned
by \code{\link{form_schema}}.

However this still misses important elements, in particular \code{labels} and
\code{choice_lists}.

\code{\link{form_schema_ext}} returns the same object as
\code{\link{form_schema}}
adding \code{labels} and \code{choice lists} in all languages available.
This is done by using the return object from \code{\link{form_xml}}.

It has the exact function signature as \code{\link{form_schema}}.
In that sense, any call to \code{\link{form_schema}} can be replaced
by \code{\link{form_schema_ext}}

This function, however, has been prepared with ODK Central version 0.8 or
higher. If you use it with an earlier version, a warning will be given.
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# With current ODK Central (>0.7)
# get extended schema:
fsx <- form_schema_ext()

# print choice list in english:
fsx[fsx$name == "test_yn", "choices_english_(en)"][[1]]

# view the extended schema:
fsx
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-individual-form/getting-form-schema-fields}

\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-individual-form/retrieving-form-schema-json}

\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-individual-form/retrieving-form-xml}

Other form-management: 
\code{\link{form_detail}()},
\code{\link{form_list}()},
\code{\link{form_schema}()},
\code{\link{form_xml}()}
}
\concept{form-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_svc}
\alias{fq_svc}
\title{OData service document for an ODK Central form.}
\format{
A tibble with one row per submission data endpoint.
}
\source{
OData service document for
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
}
\usage{
fq_svc
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The OData response for the metadata of an ODK Central form.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{geo_gj_raw}
\alias{geo_gj_raw}
\title{The unparsed submissions of a form containing geofields in GeoJSON.}
\format{
An object of class \code{list} of length 2.
}
\source{
\code{\link{odata_submission_get}(wkt=FALSE, parse=FALSE)}
run on the test form
\code{system.file("extdata", "Locations.xml", package = "ruODK")}.
}
\usage{
geo_gj_raw
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odata_submission_get.R
\name{odata_submission_get}
\alias{odata_submission_get}
\title{Retrieve and rectangle form submissions, parse dates, geopoints, download and
link attachments.}
\usage{
odata_submission_get(
  table = "Submissions",
  skip = NULL,
  top = NULL,
  count = FALSE,
  wkt = FALSE,
  filter = NULL,
  parse = TRUE,
  download = TRUE,
  orders = c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd"),
  local_dir = "media",
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  odkc_version = get_default_odkc_version(),
  tz = get_default_tz(),
  retries = get_retries(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{table}{The submission EntityType, or in plain words, the table name.
Default: \code{Submissions} (the main table).
Change to \code{Submissions.GROUP_NAME} for repeating form groups.
The group name can be found through \code{\link{odata_service_get}}.}

\item{skip}{The number of rows to be omitted from the results.
Example: 10, default: \code{NA} (none skipped).}

\item{top}{The number of rows to return.
Example: 100, default: \code{NA} (all returned).}

\item{count}{If TRUE, an \code{@odata.count} property will be returned in the
response from ODK Central. Default: \code{FALSE}.}

\item{wkt}{If TRUE, geospatial data will be returned as WKT (Well Known Text)
strings. Default: \code{FALSE}, returns GeoJSON structures.
Note that accuracy is only returned through GeoJSON.}

\item{filter}{If provided, will filter responses to those matching the query.
For an \code{odkc_version} below 1.1, this parameter will be discarded.
In ODK Central v1.1, only the fields \code{system/submitterId} and
\code{system/submissionDate} are available to reference.
In ODK Central v1.2, other fields may become available.
The operators \code{lt}, \code{lte}, \code{eq}, \code{neq}, \code{gte}, \code{gt}, \code{not}, \code{and}, and \code{or}
are supported, and the built-in functions
\code{now}, \code{year}, \code{month}, \code{day}, \code{hour}, \code{minute}, \code{second.}
\code{ruODK} does not validate the query string given to \code{filter}.
It is highly recommended to refer to the ODK Central API documentation
as well as the
\href{https://docs.oasis-open.org/odata/odata/v4.01/odata-v4.01-part1-protocol.html#_Toc31358948}{OData spec on filters}.
for filter options and capabilities}

\item{parse}{Whether to parse submission data based on form schema.
Dates and datetimes will be parsed into local time.
Attachments will be downloaded, and the field updated to the local file
path.
Point locations will be split into components; GeoJSON (\code{wkt=FALSE})
will be split into latitude, longitude, altitude and accuracy
(with anonymous field names), while WKT will be split into
longitude, latitude,and altitude (missing accuracy) prefixed by
the original field name.
See details for the handling of geotraces and geoshapes.
Default: TRUE.}

\item{download}{Whether to download attachments to \code{local_dir} or not.
If in the future ODK Central supports hot-linking attachments,
this parameter will replace attachment file names with their fully
qualified attachment URL.
Default: TRUE.}

\item{orders}{(vector of character) Orders of datetime elements for
lubridate.
Default:
\code{c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd")}.}

\item{local_dir}{The local folder to save the downloaded files to,
default: \code{"media"}.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{odkc_version}{The ODK Central version as decimal number (major.minor).
\code{ruODK} uses this parameter to adjust for breaking changes in ODK Central.

Default: \code{\link{get_default_odkc_version}} or 1.1 if unset.

Set default \code{get_default_odkc_version} through
\code{ru_setup(odkc_version=1.1)}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{tz}{A timezone to convert dates and times to.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
timezone can be set globally or per function.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
A list of lists.
\itemize{
\item \code{value} contains the submissions as list of lists.
\item \code{@odata.context} is the URL of the metadata.
\item \code{@odata.count} is the total number of rows in the table.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
\code{\link{odata_submission_get}} downloads submissions from
(default) the main form group (submission table) including any non-repeating
form groups, or from any other table as specified by parameter \code{table}.

With parameter \code{parse=TRUE} (default), submission data is parsed into a
tibble. Any fields of type \code{dateTime} or \code{date} are parsed into
dates, with an optional parameter \code{tz} to specify the local timezone.

A parameter \code{local_dir} (default: \code{media}) specifies a local
directory for downloaded attachment files.
Already existing, previously downloaded attachments will be retained.

With parameter \code{wkt=TRUE}, spatial fields will be returned as WKT, rather
than GeoJSON. In addition, fields of type \code{geopoint} will be split into
latitude, longitude, and altitude, prefixed with the original field name.
E.g. a field \code{start_location} of type \code{geopoint} will be split into
\code{start_location_latitude}, \code{start_location_longitude}, and
\code{start_location_altitude}. The field name prefix will allow multiple fields
of type \code{geopoint} to be split into their components without naming
conflicts.

Geotraces (lines) and gepshapes (polygons) will be retained in their original
format, plus columns of their first point's coordinate components as
provided by \code{\link{split_geotrace}} and \code{\link{split_geoshape}},
respectively.

Entirely unpopulated form fields, as well as notes and form groups, will be
excluded from the resulting tibble.
Submitting at least one complete form instance will prevent the accidental
exclusion of an otherwise mostly empty form field.

The only remaining manual step is to optionally join any sub-tables to the
master table.

The parameter \code{verbose} enables diagnostic messages along the download
and parsing process.

With parameter \code{parse=FALSE}, submission data is presented as nested
list, which is the R equivalent of the JSON structure returned from the API.
From there, \code{\link{odata_submission_rectangle}} can rectangle the data
into a tibble, and subsequent lines of \code{\link{handle_ru_datetimes}},
\code{\link{handle_ru_attachments}}, \code{\link{handle_ru_geopoints}},
\code{\link{handle_ru_geotraces}}, and \code{\link{handle_ru_geoshapes}}
parse dates, download and link file attachments, and extract coordinates from
geofields.
\code{ruODK} offers this manual and explicit pathway as an option to
investigate and narrow down unexpected or unwanted behaviour.
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

form_tables <- ruODK::odata_service_get()
data <- odata_submission_get() # default: main data table
data <- odata_submission_get(table = form_tables$url[1]) # same, explicitly
data_sub1 <- odata_submission_get(table = form_tables$url[2]) # sub-table 1
data_sub2 <- odata_submission_get(table = form_tables$url[3]) # sub-table 2

# Skip one row, return the next 1 rows (top), include total row count
data <- odata_submission_get(
  table = form_tables$url[1],
  skip = 1,
  top = 1,
  count = TRUE
)

# Filter submissions
data <- odata_submission_get(
  table = form_tables$url[1],
  filter = "year(__system/submissionDate) lt year(now())"
)
data <- odata_submission_get(
  table = form_tables$url[1],
  filter = "year(__system/submissionDate) lt 2020"
)
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/odata-endpoints/odata-form-service}

\url{https://odkcentral.docs.apiary.io/#reference/odata-endpoints/odata-form-service/data-document}

Other odata-api: 
\code{\link{odata_metadata_get}()},
\code{\link{odata_service_get}()}
}
\concept{odata-api}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_setup.R
\name{ru_settings}
\alias{ru_settings}
\alias{get_default_pid}
\alias{get_default_fid}
\alias{get_default_url}
\alias{get_default_un}
\alias{get_default_pw}
\alias{get_default_pp}
\alias{get_default_tz}
\alias{get_test_url}
\alias{get_test_un}
\alias{get_test_pw}
\alias{get_test_pid}
\alias{get_test_fid}
\alias{get_test_fid_zip}
\alias{get_test_fid_att}
\alias{get_test_fid_gap}
\alias{get_test_fid_wkt}
\alias{get_test_pp}
\alias{get_ru_verbose}
\alias{get_default_odkc_version}
\alias{get_test_odkc_version}
\alias{get_retries}
\title{Get or set \code{ruODK} settings.}
\usage{
ru_settings()

get_default_pid()

get_default_fid()

get_default_url()

get_default_un()

get_default_pw()

get_default_pp()

get_default_tz()

get_test_url()

get_test_un()

get_test_pw()

get_test_pid()

get_test_fid()

get_test_fid_zip()

get_test_fid_att()

get_test_fid_gap()

get_test_fid_wkt()

get_test_pp()

get_ru_verbose()

get_default_odkc_version()

get_test_odkc_version()

get_retries()
}
\value{
\code{\link{ru_settings}} prints your default ODK Central project ID,
form ID, url, username, and password, corresponding optional test
server as well as verbosity and HTTP request settings.
\code{\link{ru_setup}} sets your production and test settings, while
\code{get_(default/test)_*} get each of those respective settings.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
ru_settings()
}
\seealso{
\code{\link{ru_setup}},
\code{\link{get_default_pid}},
\code{\link{get_default_fid}},
\code{\link{get_default_url}},
\code{\link{get_default_un}},
\code{\link{get_default_pw}},
\code{\link{get_default_pp}},
\code{\link{get_default_tz}},
\code{\link{get_default_odkc_version}},
\code{\link{get_retries}},
\code{\link{get_test_pid}},
\code{\link{get_test_fid}},
\code{\link{get_test_fid_zip}},
\code{\link{get_test_fid_att}},
\code{\link{get_test_fid_gap}},
\code{\link{get_test_fid_wkt}},
\code{\link{get_test_url}},
\code{\link{get_test_un}},
\code{\link{get_test_pw}},
\code{\link{get_test_pp}},
\code{\link{get_test_odkc_version}},
\code{\link{get_ru_verbose}}.

Other ru_settings: 
\code{\link{odata_svc_parse}()},
\code{\link{ru_setup}()},
\code{\link{yell_if_error}()},
\code{\link{yell_if_missing}()}
}
\concept{ru_settings}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachment_get.R
\name{attachment_get}
\alias{attachment_get}
\title{Download attachments and return the local path.}
\usage{
attachment_get(
  sid,
  fn,
  local_dir = "media",
  separate = FALSE,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{sid}{One or many ODK submission UUIDs, an MD5 hash.}

\item{fn}{One or many ODK form attachment filenames,
e.g. "1558330537199.jpg".}

\item{local_dir}{The local folder to save the downloaded files to,
default: "media".}

\item{separate}{(logical) Whether to separate locally downloaded files into
a subfolder named after the submission uuid within \code{local_dir},
default: FALSE.
The defaults mirror the behaviour of \code{\link{submission_export}}, which
keeps all attachment files together in a folder \code{media}.
Enable this option if downloaded files collide on idential names. This can
happen if two data collection devices by chance generate the same filename
for two respective media files, e.g. \code{DCIM0001.jpg}.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The relative file path for the downloaded attachment(s)
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This function is the workhorse for
\code{\link{handle_ru_attachments}}.
This function is vectorised and can handle either one or many records.
Parameters submission_uuid and attachment_filename accept single or exactly
the same number of multiple values.
The other parameters are automatically repeated.

The media attachments are downloaded into a folder given by \code{local_dir}:

workdir/media/filename1.jpg

workdir/media/filename2.jpg

workdir/media/filename3.jpg
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

a_local_dir <- here::here()

# Step 2: Get unparsed submissions
fresh_raw <- odata_submission_get(parse = FALSE)

# Step 3: Get attachment field "my_photo"
fresh_parsed <- fresh_raw \%>\%
  odata_submission_rectangle() \%>\%
  dplyr::mutate(
    my_photo = attachment_get(id,
      my_photo,
      local_dir = a_local_dir,
      verbose = TRUE
    )
    # Repeat for all other attachment fields
  )
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-form-attachments/downloading-a-form-attachment}

\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/attachments/downloading-an-attachment}

Other utilities: 
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/form_schema.R
\name{form_schema}
\alias{form_schema}
\title{Show the schema of one form.}
\usage{
form_schema(
  flatten = FALSE,
  odata = FALSE,
  parse = TRUE,
  draft = FALSE,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  odkc_version = get_default_odkc_version(),
  retries = get_retries(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{flatten}{Whether to flatten the resulting list of lists (TRUE) or not
(FALSE, default). Only applies to ODK Central version < 0.8.}

\item{odata}{Whether to sanitise the field names to match the way they will
be outputted for OData. While the original field names as given in the
XForms definition may be used as-is for CSV output, OData has some
restrictions related to the domain-qualified identifier syntax it uses.
Only applies to ODK Central version < 0.8.
Default: FALSE.}

\item{parse}{Whether to parse the form schema into a tibble of form field
type and name. This uses \code{\link{form_schema_parse}} internally.
If used together with \code{flatten=TRUE}, \code{\link{form_schema}} will raise
a warning and return the unparsed, flattened form schema.
Only applies to ODK Central version < 0.8.
Default: TRUE.}

\item{draft}{Whether the form is published (FALSE) or a draft (TRUE).
Default: TRUE.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{odkc_version}{The ODK Central version as decimal number (major.minor).
\code{ruODK} uses this parameter to adjust for breaking changes in ODK Central.

Default: \code{\link{get_default_odkc_version}} or 1.1 if unset.

Set default \code{get_default_odkc_version} through
\code{ru_setup(odkc_version=1.1)}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
A tibble or nested list (v0.7) containing the form definition.
At the lowest nesting level, each form field consists of a list of two
nodes, \code{name} (the underlying field name) and \code{type} (the XForms field
type, as in "string", "select1", "geopoint", "binary" and so on).
These fields are nested in lists of tuples \code{name} (the XForms screen name),
\code{children} (the fields as described above), \code{type} ("structure" for non-
repeating screens, "repeat" for repeating screens).
A list with \code{name} "meta" may precede the structure, if several metadata
fields are captured (e.g. "instanceId", form start datetimes etc.).
In all cases for ODK Central 0.8, and with default parameters
(\code{parse=TRUE}) for ODK Central 0.7, \code{\link{form_schema}} returns
a tibble with the columns:

\itemize{
\item \code{name} The field name as given in the form schema.
\item \code{type} The field type, e.g. "string", "select1", etc.
\item \code{path} The XForms path of the field,
\item \code{ruodk_name} The predicted field name as generated by
\code{\link{odata_submission_get}}, prefixed by the path, additionally
cleaned with \code{\link[janitor]{make_clean_names}} to match the
cleaned column names from \code{\link{odata_submission_rectangle}}.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
ODK Central has introduced a new API endpoint in version 0.8 which
returns a parsed and flattened list of fields. This replaces the nested
form schema which is challenging to parse.

While users of newest ODK Central versions ( > 0.8) can ignore the legacy
support for ODK Central's earlier form schema API, users of ODK Central
version < 0.8 can set an environment variable \code{ODKC_VERSION} to their
ODKC's version in format \code{<major>.<minor>} e.g. \code{0.7}.
This variable caters for future breaking changes.

Either way, \code{\link{form_schema}} will always return a tibble with
columns \code{name}, \code{type}, \code{path} and \code{ruodk_name}.
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# With explicit pid and fid
fs_defaults <- form_schema(pid = 1, fid = "build_xformsId")

# With current ODK Central (v0.8)
fs <- form_schema()

# With defaults, ODK Central v0.7
fs_nested <- form_schema(
  flatten = FALSE,
  odata = FALSE,
  parse = FALSE,
  odkc_version = 0.7
)
listviewer::jsonedit(fs_nested)

fs_flattened <- form_schema(
  flatten = TRUE,
  odata = FALSE,
  parse = FALSE,
  odkc_version = 0.7
)
listviewer::jsonedit(fs_flattened)

# form_schema returns a nested list. There's nothing to change about that.
class(fs_nested)
# > "list"

class(fs_flattened)
# > "list"

# This assumes knowledge of that exact form being tested.
# First node: type "structure" (a field group) named "meta".
fs_nested[[1]]$type
# > "structure"

fs_nested[[1]]$name
# > "meta"

# The first node contains children, which means it's an XForms field group.
names(fs_nested[[1]])
# > "name" "children" "type"

# Next node: a "meta" field of type "string" capturing the  "instanceId".
# First child node of "meta": type "string", name "instanceId".
fs_nested[[1]]$children[[1]]$type
# > "string"
fs_nested[[1]]$children[[1]]$name
# > "instanceID"

# In the flattened version, the field's and it's ancestors' names are the
# components of "path".
fs_flattened[[1]]$path
# > "meta". "instanceId"

fs_flattened[[1]]$type
# > "string"

# Last node: a "meta" field capturing the datetime of form completion
fs_flattened[[length(fs_flattened)]]$type
# > "dateTime"
fs_nested[[length(fs_nested)]]$type
# > "dateTime"

# Parsed into a tibble of form field type/name:
# Useful to inform further parsing of submission data (attachments, dates)
fs <- form_schema(parse = TRUE, odkc_version = 0.7)
fs <- form_schema(odkc_version = 0.8)

# Attachments: used by handle_ru_attachments
fs \%>\% dplyr::filter(type == "binary")

# dateTime: used by handle_ru_datetimes
fs \%>\% dplyr::filter(type == "dateTime")

# Point location: used by handle_ru_geopoints
fs \%>\% dplyr::filter(type == "geopoint")
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-individual-form/getting-form-schema-fields}

\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-individual-form/retrieving-form-schema-json}

Other form-management: 
\code{\link{form_detail}()},
\code{\link{form_list}()},
\code{\link{form_schema_ext}()},
\code{\link{form_xml}()}
}
\concept{form-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odata_metadata_get.R
\name{odata_metadata_get}
\alias{odata_metadata_get}
\title{Retrieve metadata from an OData URL ending in .svc as list of lists.}
\usage{
odata_metadata_get(
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A nested list containing Edmx (dataset schema definition) and
.attrs (Version).
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

meta <- odata_metadata_get()
listviewer::jsonedit(meta)
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/odata-endpoints/odata-form-service/metadata-document}

Other odata-api: 
\code{\link{odata_service_get}()},
\code{\link{odata_submission_get}()}
}
\concept{odata-api}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/form_schema_parse.R
\name{form_schema_parse}
\alias{form_schema_parse}
\title{Parse a form_schema into a tibble of fields with name, type, and path.}
\usage{
form_schema_parse(fs, path = "Submissions", verbose = get_ru_verbose())
}
\arguments{
\item{fs}{The output of form_schema as nested list}

\item{path}{The base path for form fields. Default: "Submissions".
\code{\link{form_schema_parse}} recursively steps into deeper nesting
levels, which are reflected as separate OData tables.
The returned value in \code{path} reflects the XForms group name, which
translates to separate screens in ODK Collect.
Non-repeating form groups will be flattened out into the main Submissions
table. Repeating groups are available as separate OData tables.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This function is used by \code{\link{form_schema}} for older versions of
ODK Central (pre 0.8). These return the form schema as XML, requiring the
quite involved code of \code{\link{form_schema_parse}}, while newer ODK
Central versions return JSON, which is parsed directly in
\code{\link{form_schema}}.

The \code{form_schema} returned from ODK Central versions < 0.8 is a nested list
of lists containing the form definition.
The form definition consists of fields (with a type and name), and form
groups, which are rendered as separate ODK Collect screens.
Form groups in turn can also contain form fields.

\code{\link{form_schema_parse}} recursively unpacks the form and extracts the
name and type of each field. This information then informs
\code{\link{handle_ru_attachments}}, \code{\link{handle_ru_datetimes}},
\code{\link{handle_ru_geopoints}}, \code{\link{handle_ru_geotraces}}, and
\code{\link{handle_ru_geoshapes}}.
}
\examples{
\dontrun{
# Option 1: in two steps, ODKC Version 0.7
fs <- form_schema(flatten = FALSE, parse = FALSE, odkc_version = 0.7)
fsp <- form_schema_parse(fs)

# Option 2: in one go
fsp <- form_schema(parse = TRUE)

fsp
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_zip_strata}
\alias{fq_zip_strata}
\title{A tibble of a repeated sub-group of records from a test form.}
\format{
A tibble of repeated sub-group of records from a test form.
}
\source{
\code{\link{submission_export}}
run on the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}.
}
\usage{
fq_zip_strata
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/user_list.R
\name{user_list}
\alias{user_list}
\title{List all users.}
\usage{
user_list(
  qry = NULL,
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries(),
  orders = c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd"),
  tz = get_default_tz(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{qry}{A query string to filter users by. The query string is not case
sensitive and can contain special characters, such as \code{@}.

The query string must be at least 5 alphabetic characters long to
return good enough matches.
E.g. \code{janet} will match a user with display name \verb{Janette Doe}.
E.g., \verb{@dbca.wa} will match users with an email from \code{dbca.wa.gov.au},
whereas \verb{@dbca.w} or \verb{@dbca} will return no matches.

Default: \code{NULL}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{orders}{(vector of character) Orders of datetime elements for
lubridate.

Default:
\code{c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd")}.}

\item{tz}{A timezone to convert dates and times to.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
timezone can be set globally or per function.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
A tibble with user details as per the ODK Central API docs.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}
}
\details{
Currently, there are no paging or filtering options, so listing Users will
get you every User in the system, every time.
Optionally, a q query string parameter may be provided to filter the returned
users by any given string.
The search is performed via a trigram similarity index over both the Email
and Display Name fields, and results are ordered by match score, best matches
first.
If a q parameter is given, and it exactly matches an email address that
exists in the system, that user's details will always be returned,
even for actors who cannot user.list.
The request must still authenticate as a valid Actor.
This allows non-Administrators to choose a user for an action
(e.g. grant rights) without allowing full search.

Actors who cannot \code{user.list} will always receive [] with a 200 OK response.
ruODK does not (yet) warn if this is the case, and you (the requesting Actor)
have no permission to \code{user.list}.
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# All users
ul <- user_list()

# Search users
# Given a user with display name "Janette Doe" and email "@org.com.au"
user_list(qry = "jan") # no results, query string too short
user_list(qry = "jane") # no results, query string too short
user_list(qry = "janet") # returns Janette
user_list(qry = "@org") # no results, query string too short
user_list(qry = "@org.c") # no results, query string too short
user_list(qry = "@org.co") # returns all users matching "@org.co"

# Actor not allowed to user.list
user_list() # If this is empty, you might not have permissions to list users
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/accounts-and-users/users/listing-all-users}

\url{https://www.postgresql.org/docs/9.6/pgtrgm.html}
}
\concept{user-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_submission_list}
\alias{fq_submission_list}
\title{A tibble of submission metadata.}
\format{
A tibble of submission metadata.
}
\source{
The output of \code{\link{submission_list}}
run on the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}.
}
\usage{
fq_submission_list
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_form_list}
\alias{fq_form_list}
\title{A tibble of forms.}
\format{
A tibble of forms
}
\source{
The output of \code{\link{form_list}}.
run on the project.
}
\usage{
fq_form_list
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_setup.R
\name{ru_setup}
\alias{ru_setup}
\title{Configure default \code{\link{ruODK}} settings.}
\usage{
ru_setup(
  svc = NULL,
  pid = NULL,
  fid = NULL,
  url = NULL,
  un = NULL,
  pw = NULL,
  pp = NULL,
  tz = NULL,
  odkc_version = NULL,
  retries = NULL,
  verbose = NULL,
  test_svc = NULL,
  test_pid = NULL,
  test_fid = NULL,
  test_fid_zip = NULL,
  test_fid_att = NULL,
  test_fid_gap = NULL,
  test_fid_wkt = NULL,
  test_url = NULL,
  test_un = NULL,
  test_pw = NULL,
  test_pp = NULL,
  test_odkc_version = NULL
)
}
\arguments{
\item{svc}{(optional, character) The OData service URL of a form.
This parameter will set \code{pid}, \code{fid}, and \code{url}.
It is sufficient to supply \code{svc}, \code{un}, and \code{pw}.}

\item{pid}{(optional, character) The ID of an existing project on \code{url}.
This will override the project ID from \code{svc}.
A numeric value for \code{pid} will be converted to character.}

\item{fid}{(optional, character) The alphanumeric ID of an existing form
in \code{pid}. This will override the form ID from \code{svc}.}

\item{url}{An ODK Central URL,
e.g. "https://sandbox.central.getodk.org".
This will override the ODK Central base URL from \code{svc}.}

\item{un}{An ODK Central username which is the email of a "web user" in the
specified ODK Central instance \code{url} (optional, character).}

\item{pw}{The password for user \code{un} (optional, character).}

\item{pp}{The passphrase (optional, character) for an encrypted form.}

\item{tz}{Global default time zone.
\code{ruODK}'s time zone is determined in order of precedence:
\itemize{
\item Function parameter:
e.g. \code{\link{odata_submission_get}(tz = "Australia/Perth")}
\item \code{ruODK} setting: \code{\link{ru_setup}(tz = "Australia/Perth")}
\item Environment variable \code{RU_TIMEZONE} (e.g. set in \code{.Renviron})
\item UTC (GMT+00)
}}

\item{odkc_version}{The ODK Central version as major/minor version, e.g. 1.1.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{verbose}{Global default for \code{ruODK} verbosity.
\code{ruODK} verbosity is determined in order of precedence:
\itemize{
\item Function parameter:
e.g. \code{\link{odata_submission_get}(verbose = TRUE)}
\item \code{ruODK} setting: \code{\link{ru_setup}(verbose = TRUE)}
\item Environment variable \code{RU_VERBOSE} (e.g. set in \code{.Renviron})
\item \code{FALSE}.
}}

\item{test_svc}{(optional, character) The OData service URL of a test form.
This parameter will set \code{test_pid}, \code{test_fid}, and
\code{test_url}. It is sufficient to supply \code{test_svc},
\code{test_un}, and \code{test_pw} to configure testing.}

\item{test_pid}{(optional, character) The numeric ID of an existing project
on \code{test_url}. This will override the project ID from \code{test_svc}.
A numeric value for \code{test_pid} will be converted to character.}

\item{test_fid}{(optional, character) The alphanumeric ID of an existing form
in \code{test_pid}. This will override the form ID from \code{test_svc}.
This form is used as default form in all tests, examples, vignettes, data,
and Rmd templates.}

\item{test_fid_zip}{(optional, character) The alphanumeric ID of an existing
form in \code{test_pid}. This will override the form ID from
\code{test_svc}.
Provide the form ID of a form with few submissions and without attachments.
This form is used to test the repeated download of all form submissions.}

\item{test_fid_att}{(optional, character) The alphanumeric ID of an existing
form in \code{test_pid}. This will override the form ID from
\code{test_svc}.
Provide the form ID of a form with few submissions and few attachments.
This form is used to test downloading and linking attachments.}

\item{test_fid_gap}{(optional, character) The alphanumeric ID of an existing
form in \code{test_pid}. This will override the form ID from
\code{test_svc}.
Provide the form ID of a form with gaps in the first submission.
This form is used to test parsing incomplete submissions.}

\item{test_fid_wkt}{(optional, character) The alphanumeric ID of an existing
form in \code{test_pid}. This will override the form ID from
\code{test_svc}.
Provide the form ID of a form with geopoints, geotraces, and geoshapes.}

\item{test_url}{(optional, character) A valid ODK Central URL for testing.
This will override the ODK Central base URL from \code{svc}.}

\item{test_un}{(optional, character) A valid ODK Central username (email)
privileged to view the test project(s) at \code{test_url}.}

\item{test_pw}{(optional, character) The valid ODK Central password for
\code{test_un}.}

\item{test_pp}{(optional, character) The valid passphrase to decrypt the
data of encrypted project \code{test_pid} for download. Only used for tests.}

\item{test_odkc_version}{The ODK Central test server's version as major/minor
version, e.g. 1.1.}
}
\description{
Settings are returned invisibly and additionally printed depending on
\code{\link{get_ru_verbose}}.
}
\details{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}

\code{\link{ru_setup}} sets ODK Central connection details.
\code{\link{ruODK}}'s functions default to use the default project ID,
form ID, URL, username, and password unless specified explicitly.

Any parameters not specified will remain unchanged. It is therefore possible
to set up username and password initially with
\code{ru_setup(un="XXX", pw="XXX")}, and switch between forms with
\code{ru_setup(svc="XXX")}, supplying the form's OData service URL.
ODK Central conveniently provides the OData service URL in the form
submission tab, which in turn contains base URL, project ID, and form ID.

\code{\link{ruODK}}'s automated tests require a valid ODK Central URL, and a
privileged username and password of a "web user" on that ODK Central
instance, as well as an existing project and form.
}
\examples{
# `ruODK` users only need default settings to their ODK Central:
ru_setup(url = "https://my-odkc.com", un = "me@email.com", pw = "...")

# `ruODK` contributors and maintainers need specific ODK Central
# instances to run tests and build vignettes, see contributing guide:
ru_setup(
  url = "https://odkcentral.dbca.wa.gov.au",
  un = "me@email.com",
  pw = "...",
  pp = "...",
  test_url = "https://sandbox.central.getodk.org",
  test_un = "me@email.com",
  test_pw = "...",
  test_pp = "...",
  test_pid = 14,
  test_fid = "build_Flora-Quadrat-0-2_1558575936",
  test_fid_zip = "build_Spotlighting-0-6_1558333698",
  test_fid_att = "build_Flora-Quadrat-0-1_1558330379",
  test_fid_gap = "build_Turtle-Track-or-Nest-1-0_1569907666",
  test_fid_wkt = "build_Locations_1589344221",
  retries = 3,
  verbose = TRUE
)
}
\seealso{
Other ru_settings: 
\code{\link{odata_svc_parse}()},
\code{\link{ru_settings}()},
\code{\link{yell_if_error}()},
\code{\link{yell_if_missing}()}
}
\concept{ru_settings}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/submission_list.R
\name{submission_list}
\alias{submission_list}
\title{List all submissions of one form.}
\usage{
submission_list(
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries(),
  orders = c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd"),
  tz = get_default_tz()
)
}
\arguments{
\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{orders}{(vector of character) Orders of datetime elements for
lubridate.

Default:
\code{c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd")}.}

\item{tz}{A timezone to convert dates and times to.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
timezone can be set globally or per function.}
}
\value{
A tibble containing some high-level details of the form submissions.
One row per submission, columns are submission attributes:\preformatted{    * instance_id: uuid, string. The unique ID for each submission.
    * submitter_id: user ID, integer.
    * created_at: time of submission upload, dttm
    * updated_at: time of submission update on server, dttm or NA
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# Set default credentials, see vignette("setup")
ruODK::ru_setup(
  svc = paste0(
    "https://sandbox.central.getodk.org/v1/projects/14/",
    "forms/build_Flora-Quadrat-0-2_1558575936.svc"
  ),
  un = "me@email.com",
  pw = "..."
)

sl <- submission_list()
sl \%>\% knitr::kable(.)

fl <- form_list()

# submission_list returns a tibble
class(sl)
# > c("tbl_df", "tbl", "data.frame")

# Submission attributes are the tibble's columns
names(sl)
# > "instance_id" "submitter_id" "device_id" "created_at" "updated_at"

# Number of submissions (rows) is same as advertised in form_list
form_list_nsub <- fl \%>\%
  filter(fid == get_test_fid()) \%>\%
  magrittr::extract2("submissions") \%>\%
  as.numeric()
nrow(sl) == form_list_nsub
# > TRUE
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/submissions/listing-all-submissions-on-a-form}

Other submission-management: 
\code{\link{attachment_list}()},
\code{\link{encryption_key_list}()},
\code{\link{submission_detail}()},
\code{\link{submission_export}()},
\code{\link{submission_get}()}
}
\concept{submission-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_data}
\alias{fq_data}
\title{Parsed submission data for an ODK Central form.}
\format{
The output of \code{\link{odata_submission_get}} for a set of example
data. A tidy tibble referencing the attachments included in the vignettes
and documentation at a relative path \verb{attachments/media/<filename>.<ext>}.
}
\source{
See \code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
and \code{\link{odata_submission_get}}.
}
\usage{
fq_data
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The parsed OData response for the submissions of an ODK Central form.
This form represents a Flora Quadrat, which is a ca 50 by 50 m quadrat of
a uniform plant community.

The XML and .odkbuild versions for this form are available as
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
and \code{system.file("extdata", "FloraQuadrat04.odkbuild", package = "ruODK")},
respectively.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odata_submission_rectangle.R
\name{odata_submission_rectangle}
\alias{odata_submission_rectangle}
\title{Rectangle the output of \code{\link{odata_submission_get}(parse=FALSE)}
into a tidy tibble and unnest all levels.}
\usage{
odata_submission_rectangle(
  data,
  names_repair = "universal",
  names_sep = "_",
  form_schema = NULL,
  verbose = get_ru_verbose()
)
}
\arguments{
\item{data}{A nested list of lists as given by
\code{\link{odata_submission_get}}.}

\item{names_repair}{The argument \code{names_repair} for
\code{tidyr::unnest_wider}, default: "universal".}

\item{names_sep}{The argument \code{names_sep} for
\code{tidyr::unnest_wider}, default: "_".
Un-nested variables inside a list column will be prefixed by the list
column name, separated by \code{names_sep}.
This avoids unsightly repaired names such as \code{latitude...1}.}

\item{form_schema}{An optional form_schema,
like the output of \code{\link{form_schema}}. If a form schema is supplied,
location fields will not be unnested. While WKT location fields contain
plain text and will never be unnested, GeoJSON location fields would cause
errors during unnesting.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The submissions as un-nested tibble
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# Using canned data
data_parsed <- odata_submission_rectangle(fq_raw, verbose = TRUE)
# Field "device_id" is known part of fq_raw
testthat::expect_equal(
  data_parsed$device_id[[1]],
  fq_raw$value[[1]]$device_id
)

# fq_raw has two submissions
testthat::expect_equal(length(fq_raw$value), nrow(data_parsed))
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
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
To learn more about tidy eval and how to use these tools, read
\url{http://rlang.tidyverse.org/articles/tidy-evaluation.html}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{unnest_all}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachment_list.R
\name{get_one_submission_attachment_list}
\alias{get_one_submission_attachment_list}
\title{List all attachments of one submission.}
\usage{
get_one_submission_attachment_list(
  iid,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{iid}{The \code{instance_id}, a UUID, as returned by
\code{\link{submission_list}}.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A tibble containing some high-level details of the submission
attachments.
One row per submission attachment, columns are submission attributes:\preformatted{    * name: The attachment filename, e.g. 12345.jpg
    * exists: Whether the attachment for that submission exists on the
      server.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
When a Submission is created, either over the OpenRosa or the REST interface,
its XML data is analysed to determine which file attachments it references:
these may be photos or video taken as part of the survey, or an audit/timing
log, among other things. Each reference is an expected attachment, and these
expectations are recorded permanently alongside the Submission. With this
subresource, you can list the expected attachments, see whether the server
actually has a copy or not, and download, upload, re-upload, or clear binary
data for any particular attachment.

You can retrieve the list of expected Submission attachments at this route,
along with a boolean flag indicating whether the server actually has a copy
of the expected file or not. If the server has a file, you can then append
its filename to the request URL to download only that file.
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

sl <- submission_list()

al <- get_one_submission_attachment_list(sl$instance_id[[1]])
al \%>\% knitr::kable(.)

# attachment_list returns a tibble
class(al)
# > c("tbl_df", "tbl", "data.frame")

# Submission attributes are the tibble's columns
names(al)
# > "name" "exists"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/attachments/listing-expected-submission-attachments}

\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-form-attachments/listing-expected-form-attachments}

Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_ru_attachments.R
\name{handle_ru_attachments}
\alias{handle_ru_attachments}
\title{Download and link submission attachments according to a form schema.}
\usage{
handle_ru_attachments(
  data,
  form_schema,
  local_dir = "media",
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{data}{Submissions rectangled into a tibble. E.g. the output of\preformatted{ruODK::odata_submission_get(parse = FALSE) \%>\%
ruODK::odata_submission_rectangle()
}}

\item{form_schema}{The \code{form_schema} for the submissions.
E.g. the output of \code{ruODK::form_schema()}.}

\item{local_dir}{The local folder to save the downloaded files to,
default: "media".}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The submissions tibble with all attachments downloaded and linked to
a \code{local_dir}.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
For a given tibble of submissions, download and link attachments
for all columns which are marked in the form schema as type "binary".
}
\examples{
\dontrun{
library(magrittr)
data("fq_raw")
data("fq_form_schema")
t <- tempdir()
fs::dir_ls(t) \%>\% fs::file_delete()
fq_with_att <- fq_raw \%>\%
  ruODK::odata_submission_rectangle() \%>\%
  ruODK::handle_ru_attachments(
    form_schema = fq_form_schema,
    local_dir = t,
    pid = ruODK::get_test_pid(),
    fid = ruODK::get_test_fid(),
    url = ruODK::get_test_url(),
    un = ruODK::get_test_un(),
    pw = ruODK::get_test_pw(),
    verbose <- ruODK::get_ru_verbose()
  )
# There should be files in local_dir
testthat::expect_true(fs::dir_ls(t) \%>\% length() > 0)
}

}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_list.R
\name{project_list}
\alias{project_list}
\title{List all projects.}
\usage{
project_list(
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries(),
  orders = c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd"),
  tz = get_default_tz()
)
}
\arguments{
\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{orders}{(vector of character) Orders of datetime elements for
lubridate.

Default:
\code{c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd")}.}

\item{tz}{A timezone to convert dates and times to.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
timezone can be set globally or per function.}
}
\value{
A tibble with one row per project and all project metadata
as columns as per ODK Central API docs.
}
\description{
While the API endpoint will return all projects,
\code{\link{project_list}} will fail with incorrect or missing
authentication.
}
\details{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

pl <- project_list()
knitr::kable(pl)

# project_list returns a tibble
class(pl)
# > "tbl_df" "tbl" "data.frame"

# columns are project metadata
names(pl)
# > "id" "name" "forms" "app_users" "created_at" "updated_at"
# > "last_submission" "archived"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/project-management/projects/listing-projects}

Other project-management: 
\code{\link{project_create}()},
\code{\link{project_detail}()}
}
\concept{project-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_create.R
\name{project_create}
\alias{project_create}
\title{Create a new project.}
\usage{
project_create(
  name,
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw()
)
}
\arguments{
\item{name}{The desired name of the project. Can contain whitespace.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}
}
\value{
A tibble with one row per project and all project metadata
as columns as per ODK Central API docs.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#experimental}{\figure{lifecycle-experimental.svg}{options: alt='[Experimental]'}}}{\strong{[Experimental]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

p <- project_create("Test Project")
knitr::kable(p)

# project_create returns a tibble
class(p)
# > "tbl_df" "tbl" "data.frame"

# columns are project metadata
names(p)
# > "id" "name" "archived"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/project-management/projects/creating-a-project}

Other project-management: 
\code{\link{project_detail}()},
\code{\link{project_list}()}
}
\concept{project-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_data_taxa}
\alias{fq_data_taxa}
\title{Parsed submission data for a subgroup of an ODK Central form.}
\format{
The output of \code{\link{odata_submission_get}} for a set of example
data. A tidy tibble referencing the attachments included in the vignettes
and documentation at a relative path \verb{attachments/media/<filename>.<ext>}.
}
\source{
See \code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
and \code{\link{odata_submission_get}}.
}
\usage{
fq_data_taxa
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The parsed OData response for a subgroup of an ODK Central form.

This subgroup represents an individual plant taxon which is encountered by
the enumerators. Typically, one voucher specimen is taken for each distinct
encountered plant taxon. A field name is allocated by the enumerators, which
can be the proper canonical name (if known) or any other moniker.
The voucher specimens are later determined by taxonomic experts, who then
provide the real, terminal taxonomic name for a given voucher specimen.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{geo_gj88}
\alias{geo_gj88}
\title{The parsed submissions of a form containing geofields in GeoJSON
with trailing empty coordinates present.}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 1 rows and 51 columns.
}
\source{
\code{\link{odata_submission_get}(wkt=FALSE, parse=TRUE)}
run on the test form
\code{system.file("extdata", "Locations.xml", package = "ruODK")}.
}
\usage{
geo_gj88
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This issue was fixed in #88.
ODK Central versions 0.7 - 0.9 export geotraces and geoshapes with trailing
empty coordinates. ruODK has a patch to drop trailing empty coordinates.
This dataset is used to test the patch in ruODK.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/split_geopoint.R
\name{split_geopoint}
\alias{split_geopoint}
\title{Annotate a dataframe containing a geopoint column with lon, lat, alt.}
\usage{
split_geopoint(data, colname, wkt = FALSE)
}
\arguments{
\item{data}{(dataframe) A dataframe with a geopoint column.}

\item{colname}{(chr) The name of the geopoint column.
This column will be retained.}

\item{wkt}{Whether geofields are GeoJSON (if FALSE) or WKT
strings (if TRUE), default: FALSE.}
}
\value{
The given dataframe with the WKT POINT column colname, plus
three new columns, \code{colname_longitude}, \code{colname_latitude},
\code{colname_altitude}.
The three new columns are prefixed with the original \code{colname} to
avoid naming conflicts with any other geopoint columns.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This function is used by \code{\link{handle_ru_geopoints}}
on all \code{geopoint} fields as per \code{\link{form_schema}}.
}
\examples{
\dontrun{
df_wkt <- tibble::tibble(
  stuff = c("asd", "sdf", "sdf"),
  loc = c(
    "POINT (115.99 -32.12 20.01)",
    "POINT (116.12 -33.34 15.23)",
    "POINT (114.01 -31.56 23.56)"
  )
)
df_wkt_split <- df \%>\% split_geopoint("loc", wkt = TRUE)
testthat::expect_equal(
  names(df_wkt_split),
  c("stuff", "loc", "loc_longitude", "loc_latitude", "loc_altitude")
)

# With package data
data("geo_fs")
data("geo_wkt_raw")
data("geo_gj_raw")

# Find variable names of geopoints
geo_fields <- geo_fs \%>\%
  dplyr::filter(type == "geopoint") \%>\%
  magrittr::extract2("ruodk_name")
geo_fields[1] # First geotrace in data: point_location_point_gps

# Rectangle but don't parse submission data (GeoJSON and WKT)
geo_gj_rt <- geo_gj_raw \%>\%
  odata_submission_rectangle(form_schema = geo_fs)
geo_wkt_rt <- geo_wkt_raw \%>\%
  odata_submission_rectangle(form_schema = geo_fs)

# Data with first geopoint split
gj_first_gt <- split_geopoint(geo_gj_rt, geo_fields[1], wkt = FALSE)
gj_first_gt$point_location_point_gps_longitude

wkt_first_gt <- split_geopoint(geo_wkt_rt, geo_fields[1], wkt = TRUE)
wkt_first_gt$point_location_point_gps_longitude
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachment_get.R
\name{prepend_uuid}
\alias{prepend_uuid}
\title{Prepend a leading "uuid:" to any string, e.g. an md5 hash.}
\usage{
prepend_uuid(md5hash)
}
\arguments{
\item{md5hash}{A string, e.g. an md5 hash.}
}
\value{
The string with a prepended "uuid:"
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This is the inverse of the helper function \code{\link{strip_uuid}}.
}
\examples{
\dontrun{
prepend_uuid("1234")
prepend_uuid("d3bcefea-32a8-4dbc-80ca-4ecb0678e2b0")
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{geo_wkt88}
\alias{geo_wkt88}
\title{The parsed submissions of a form containing geofields in WKT
with trailing empty coordinates present.}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 1 rows and 48 columns.
}
\source{
\code{\link{odata_submission_get}(wkt=TRUE, parse=TRUE)}
run on the test form
\code{system.file("extdata", "Locations.xml", package = "ruODK")}.
}
\usage{
geo_wkt88
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This issue was fixed in #88.
ODK Central versions 0.7 - 0.9 export geotraces and geoshapes with trailing
empty coordinates. ruODK has a patch to drop trailing empty coordinates.
This dataset is used to test the patch in ruODK.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_form_detail}
\alias{fq_form_detail}
\title{A tibble of form metadata.}
\format{
A tibble of form metadata.
}
\source{
The output of \code{\link{form_detail}}
run on submissions of the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}.
}
\usage{
fq_form_detail
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/submission_get.R
\name{submission_get}
\alias{submission_get}
\title{Get submissions for a list of submission instance IDs.}
\usage{
submission_get(
  iid,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{iid}{A list of submission instance IDs, e.g. from
\code{\link{submission_list}$instance_id}.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A nested list of submission data.
}
\description{
Uses \code{\link{get_one_submission}} on a list of submission instance IDs
(\code{iid}) as returned from \code{\link{submission_list}$instance_id}.
By giving the list of \code{iid} to download explicitly, that list can be
modified using information not accessible to \code{ruODK},
e.g. \code{iid} can be restricted to "only not already downloaded submissions".
}
\details{
Forms with submission audit enabled will also receive the submission audit
as \code{audit.csv}. This will overwrite all previous \code{audit.csv} files.
To get the combined submission audit logs as one single, concatenated
\code{audit.csv} file, use \code{submission_export}.
}
\examples{
\dontrun{
# Step 1: Setup ruODK with OData Service URL (has url, pid, fid)
ruODK::ru_setup(svc = "...")

# Step 2: List all submissions of form
sl <- submission_list()

# Step 3: Get submissions
subs <- submission_get(sl$instance_id)
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/submissions/retrieving-submission-xml}

Other submission-management: 
\code{\link{attachment_list}()},
\code{\link{encryption_key_list}()},
\code{\link{submission_detail}()},
\code{\link{submission_export}()},
\code{\link{submission_list}()}
}
\concept{submission-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/audit_get.R
\name{audit_get}
\alias{audit_get}
\title{Get server audit log entries.}
\usage{
audit_get(
  action = NULL,
  start = NULL,
  end = NULL,
  limit = NULL,
  offset = NULL,
  url = Sys.getenv("ODKC_URL"),
  un = Sys.getenv("ODKC_UN"),
  pw = Sys.getenv("ODKC_PW"),
  retries = get_retries()
)
}
\arguments{
\item{action}{string. The action to filter the logs, e.g. "user.create".
See \url{https://odkcentral.docs.apiary.io/#reference/system-endpoints/server-audit-logs/}
for the full list of available actions.}

\item{start}{string. The ISO8601 timestamp of the earliest log entry to
return.
E.g. \verb{2000-01-01z} or \verb{2000-12-31T23:59.999z}, \verb{2000-01-01T12:12:12+08} or
\code{2000-01-01+08}.}

\item{end}{string. The ISO8601 timestamp of the last log entry to return.}

\item{limit}{integer. The max number of log entries to return.}

\item{offset}{integer. The number of log entries to skip.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A tibble containing server audit logs.
One row per audited action, columns are submission attributes:
\itemize{
\item actor_id: integer. The ID of the actor, if any, that initiated the
action.
\item action: string. The action that was taken.
\item actee_id: uuid, string. The ID of the permissioning object against
which the action was taken.
\item details: list. Additional details about the action that vary
according to the type of action.
\item logged_at: dttm. Time of action on server.
}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
Parameters to filter the audit logs:
\verb{action=form.create&start=2000-01-01z&end=2000-12-31T23\%3A59.999z}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

logs <- audit_get()

# With search parameters
logs <- audit_get(
  action = "project.update",
  start = "2019-08-01Z",
  end = "2019-08-31Z",
  limit = 100,
  offset = 0
)

# With partial search parameters
logs <- audit_get(
  limit = 100,
  offset = 0
)

logs \%>\% knitr::kable(.)

# audit_get returns a tibble
class(logs)
# > c("tbl_df", "tbl", "data.frame")

# Audit details
names(logs)
# > "actor_id" "action" "actee_id" "details" "logged_at"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/system-endpoints/server-audit-logs/getting-audit-log-entries}
}
\concept{server-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_ru_datetimes.R
\name{handle_ru_datetimes}
\alias{handle_ru_datetimes}
\title{Parse datetimes of submission data according to a form schema.}
\usage{
handle_ru_datetimes(
  data,
  form_schema,
  orders = c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd"),
  tz = get_default_tz(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{data}{Submissions rectangled into a tibble. E.g. the output of\preformatted{ruODK::odata_submission_get(parse = FALSE) \%>\%
ruODK::odata_submission_rectangle()
}}

\item{form_schema}{The \code{form_schema} for the submissions.
E.g. the output of \code{ruODK::form_schema()}.}

\item{orders}{(vector of character) Orders of datetime elements for
lubridate.

Default:
\code{c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd")}.}

\item{tz}{A timezone to convert dates and times to.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
timezone can be set globally or per function.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The submissions tibble with all date/dateTime columns mutated as
lubridate datetimes.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
For a given tibble of submissions, parse all columns which are
marked in the form schema as type "date" or "dateTime" using a set of
lubridate orders and a given timezone.
}
\examples{
\dontrun{
library(magrittr)
data("fq_raw")
data("fq_form_schema")

fq_with_dates <- fq_raw \%>\%
  ruODK::odata_submission_rectangle() \%>\%
  ruODK::handle_ru_datetimes(form_schema = fq_form_schema)

dplyr::glimpse(fq_with_dates)
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_submissions}
\alias{fq_submissions}
\title{A nested list of submission data.}
\format{
A nested list of submission data.
}
\source{
The output of \code{\link{submission_get}}
run on the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
using submission instance IDs from \code{\link{submission_list}}.
}
\usage{
fq_submissions
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_setup.R
\name{yell_if_error}
\alias{yell_if_error}
\title{Warn about failed web requests and give helpful troubleshooting tips.}
\usage{
yell_if_error(response, url, un, pw, pid = NULL, fid = NULL)
}
\arguments{
\item{response}{A httr response object}

\item{url}{A URL (character)}

\item{un}{A username (character)}

\item{pw}{A password (character)}

\item{pid}{A project ID (numeric, optional)}

\item{fid}{A form ID (character, optional)}
}
\value{
The response object
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
A wrapper around \code{httr::stop_for_status} with a more helpful error
message.
Examples: see tests for \code{\link{project_list}}.
This function is used internally but may be useful for debugging and
\code{\link{ruODK}} development.
}
\seealso{
Other ru_settings: 
\code{\link{odata_svc_parse}()},
\code{\link{ru_settings}()},
\code{\link{ru_setup}()},
\code{\link{yell_if_missing}()}
}
\concept{ru_settings}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-lifecycle.R
\name{lifecycle_shim}
\alias{lifecycle_shim}
\title{Shim to allow Import of lifecycle, required for building docs.}
\usage{
lifecycle_shim()
}
\description{
HT Jim Hester, Lionel Henry, Jenny Bryan for advice
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachment_link.R
\name{attachment_link}
\alias{attachment_link}
\title{Prefix attachment columns from CSV export with a local attachment file path.}
\usage{
attachment_link(data_tbl, form_schema, att_path = "media")
}
\arguments{
\item{data_tbl}{The downloaded submissions from
\code{\link{submission_export}} read into a \code{tibble} by
\code{readr::read_csv}.}

\item{form_schema}{The \code{form_schema} for the submissions.
E.g. the output of \code{ruODK::form_schema()}.}

\item{att_path}{A local path, default: "media" (as per .csv.zip export).
Selected columns of the dataframe (containing attchment filenames) are
prefixed with \code{att_path}, thus turning them into relative paths.}
}
\value{
The dataframe with attachment columns modified to contain relative
paths to the downloaded attachment files.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
t <- tempdir()
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# Predict filenames (with knowledge of form)
fid <- get_default_fid()
fid_csv <- fs::path(t, glue::glue("{fid}.csv"))
fid_csv_tae <- fs::path(t, glue::glue("{fid}-taxon_encounter.csv"))
fs <- form_schema()

# Download the zip file
se <- ruODK::submission_export(
  local_dir = t,
  overwrite = FALSE,
  verbose = TRUE
)

# Unpack the zip file
f <- unzip(se, exdir = t)
fs::dir_ls(t)

# Prepend attachments with media/ to turn into relative file paths
data_quadrat <- fid_csv \%>\%
  readr::read_csv(na = c("", "NA", "na")) \%>\%
  janitor::clean_names() \%>\%
  handle_ru_datetimes(fs) \%>\%
  attachment_link(fs)
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils-pipe.R
\name{\%>\%}
\alias{\%>\%}
\title{Pipe operator}
\usage{
lhs \%>\% rhs
}
\description{
See \code{magrittr::\link[magrittr:pipe]{\%>\%}} for details.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{geo_wkt}
\alias{geo_wkt}
\title{The parsed submissions of a form containing geofields in WKT.}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 1 rows and 49 columns.
}
\source{
\code{\link{odata_submission_get}(wkt=TRUE, parse=TRUE)}
run on the test form
\code{system.file("extdata", "Locations.xml", package = "ruODK")}.
}
\usage{
geo_wkt
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_zip_taxa}
\alias{fq_zip_taxa}
\title{A tibble of a repeated sub-group of records from a test form.}
\format{
A tibble of repeated sub-group of records from a test form.
}
\source{
\code{\link{submission_export}}
run on the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}.
}
\usage{
fq_zip_taxa
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/form_detail.R
\name{form_detail}
\alias{form_detail}
\title{Show details for one form.}
\usage{
form_detail(
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A tibble with one row and all form metadata as columns.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# With explicit credentials, see tests
fl <- form_list()

# The first form in the test project
f <- form_detail(fid = fl$fid[[1]])

# form_detail returns exactly one row
nrow(f)
# > 1

# form_detail returns all form metadata as columns: name, xmlFormId, etc.
names(f)

# > "name" "fid" "version" "state" "submissions" "created_at"
# > "created_by_id" "created_by" "updated_at" "published_at"
# > "last_submission" "hash"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/'-individual-form}

Other form-management: 
\code{\link{form_list}()},
\code{\link{form_schema_ext}()},
\code{\link{form_schema}()},
\code{\link{form_xml}()}
}
\concept{form-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{geo_wkt_raw}
\alias{geo_wkt_raw}
\title{The unparsed submissions of a form containing geofields in WKT.}
\format{
An object of class \code{list} of length 2.
}
\source{
\code{\link{odata_submission_get}(wkt=TRUE, parse=FALSE)}
run on the test form
\code{system.file("extdata", "Locations.xml", package = "ruODK")}.
}
\usage{
geo_wkt_raw
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_raw_strata}
\alias{fq_raw_strata}
\title{OData submission data for a subgroup of an ODK Central form.}
\format{
A list of lists
}
\source{
See \code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
}
\usage{
fq_raw_strata
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The OData response for the subgroup of an ODK Central form.

This subgroup represents vegetation strata as per the NVIS classification.
A vegetation stratum is a layer of plants with the same height, and dominated
by one or few plant taxa. Plant communities can be made of up to five strata,
with two to three being most common.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odata_submission_rectangle.R
\name{listcol_names}
\alias{listcol_names}
\title{A functional to extract names of list columns from a tibble.}
\usage{
listcol_names(tbl)
}
\arguments{
\item{tbl}{A tibble, possibly with list columns}
}
\value{
A vector of list column names
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/encryption_key_list.R
\name{encryption_key_list}
\alias{encryption_key_list}
\title{List all encryption keys for a form.}
\usage{
encryption_key_list(
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries(),
  orders = c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd"),
  tz = get_default_tz()
)
}
\arguments{
\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{orders}{(vector of character) Orders of datetime elements for
lubridate.

Default:
\code{c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz", "Ymd", "ymd")}.}

\item{tz}{A timezone to convert dates and times to.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
timezone can be set globally or per function.}
}
\value{
A tibble of encryption keys.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}
}
\details{
This endpoint provides a listing of all known encryption keys needed to
decrypt all Submissions for a given Form. It will return at least the
\code{base64RsaPublicKey} property (as column \code{public}) of all known versions
of the form that have submissions against them.
If managed keys are being used and a hint was provided, that will be returned
as well.
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

x <- encryption_key_list(
  pid = Sys.getenv("ODKC_TEST_PID_ENC"),
  fid = Sys.getenv("ODKC_TEST_FID_ENC"),
  url = get_test_url(),
  un = get_test_un(),
  pw = get_test_pw()
)

names(x)
# > [1] "id" "public" "managed" "hint" "created_at"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/submissions/listing-encryption-keys}

Other submission-management: 
\code{\link{attachment_list}()},
\code{\link{submission_detail}()},
\code{\link{submission_export}()},
\code{\link{submission_get}()},
\code{\link{submission_list}()}
}
\concept{submission-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ruODK.R
\docType{package}
\name{ruODK-package}
\alias{ruODK}
\alias{ruODK-package}
\title{ruODK: An R Client for the ODK Central API}
\description{
\code{\link{ruODK}} is an R Client for the ODK Central API.

Please see the \code{ruODK} website for full documentation:
\url{https://docs.ropensci.org/ruODK/}

\code{ruODK} is "pipe-friendly" and re-exports \verb{\\\%>\\\%} and \verb{\\\%||\\\%}, but does not
require their use.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/ruODK}
  \item \url{https://github.com/ropensci/ruODK}
  \item Report bugs at \url{https://github.com/ropensci/ruODK/issues}
}

}
\author{
\strong{Maintainer}: Florian W. Mayer \email{Florian.Mayer@dbca.wa.gov.au} (\href{https://orcid.org/0000-0003-4269-4242}{ORCID})

Other contributors:
\itemize{
  \item Maëlle Salmon \email{maelle.salmon@yahoo.se} (\href{https://orcid.org/0000-0002-2815-0399}{ORCID}) [reviewer]
  \item Karissa Whiting (\href{https://orcid.org/0000-0002-4683-1868}{ORCID}) [reviewer]
  \item Jason Taylor [reviewer]
  \item Marcelo Tyszler (\href{https://orcid.org/0000-0002-4573-0002}{ORCID}) [contributor]
  \item Hélène Langet (\href{https://orcid.org/0000-0002-6758-2397}{ORCID}) [contributor]
  \item DBCA [copyright holder, funder]
  \item NWSFTCP [funder]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_msg.R
\name{ru_msg_success}
\alias{ru_msg_success}
\title{Print a green success message with a tick symbol.}
\usage{
ru_msg_success(message, verbose = get_ru_verbose())
}
\arguments{
\item{message}{(chr) A message to print}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
ru_msg_success("This is a success message.")
}
\seealso{
Other messaging: 
\code{\link{ru_msg_abort}()},
\code{\link{ru_msg_info}()},
\code{\link{ru_msg_noop}()},
\code{\link{ru_msg_warn}()}
}
\concept{messaging}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachment_get.R
\name{strip_uuid}
\alias{strip_uuid}
\title{Strip the leading "uuid:" from a UUID hash.}
\usage{
strip_uuid(uuid)
}
\arguments{
\item{uuid}{A string which may contain any number of "uuid:"}
}
\value{
The string with every occurrence of "uuid:" deleted.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This is a helper function used by \code{\link{attachment_get}}.
}
\examples{
\dontrun{
strip_uuid("uuid:1234")
strip_uuid("uuid:d3bcefea-32a8-4dbc-80ca-4ecb0678e2b0")
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_ru_geopoints.R
\name{handle_ru_geopoints}
\alias{handle_ru_geopoints}
\title{Split all geopoints of a submission tibble into their components.}
\usage{
handle_ru_geopoints(data, form_schema, wkt = FALSE, verbose = get_ru_verbose())
}
\arguments{
\item{data}{Submissions rectangled into a tibble. E.g. the output of\preformatted{ruODK::odata_submission_get(parse = FALSE) \%>\%
ruODK::odata_submission_rectangle()
}}

\item{form_schema}{The \code{form_schema} for the submissions.
E.g. the output of \code{ruODK::form_schema()}.}

\item{wkt}{Whether geofields are GeoJSON (if FALSE) or WKT
strings (if TRUE), default: FALSE.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The submissions tibble with all geopoints retained in their original
format, plus columns of their coordinate components as provided by
\code{\link{split_geopoint}}.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
For a given tibble of submissions, find all columns which are listed
in the form schema as type \code{geopoint}, and extract their components.
Extracted components are longitude (X), latitude (Y), altitude (Z, where
given), and accuracy (M, where given).

The original column is retained to allow parsing into other spatially
enabled formats.
}
\examples{
library(magrittr)
data("geo_fs")
data("geo_gj_raw")
data("geo_wkt_raw")

# GeoJSON
geo_gj_parsed <- geo_gj_raw \%>\%
  ruODK::odata_submission_rectangle(form_schema = geo_fs) \%>\%
  ruODK::handle_ru_geopoints(form_schema = geo_fs, wkt = FALSE)

dplyr::glimpse(geo_gj_parsed)

# WKT
geo_wkt_parsed <- geo_wkt_raw \%>\%
  ruODK::odata_submission_rectangle(form_schema = geo_fs) \%>\%
  ruODK::handle_ru_geopoints(form_schema = geo_fs, wkt = TRUE)

dplyr::glimpse(geo_wkt_parsed)
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_form_xml}
\alias{fq_form_xml}
\title{A nested list of a form definition.}
\format{
A nested list of a form definition.
}
\source{
The output of \code{\link{form_xml}}
run on the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}.
}
\usage{
fq_form_xml
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/isodt_to_local.R
\name{isodt_to_local}
\alias{isodt_to_local}
\title{Parse an ISO8601 datetime string to a timezone aware datetime.}
\usage{
isodt_to_local(
  datetime_string,
  orders = c("YmdHMS", "YmdHMSz"),
  tz = get_default_tz()
)
}
\arguments{
\item{datetime_string}{(character) An ISO8601 datetime string as produced by
XForms exported from ODK Central.}

\item{orders}{(vector of character) Orders of datetime elements for
lubridate.
Default: \code{c("YmdHMS", "YmdHMSz", "Ymd HMS", "Ymd HMSz")}.}

\item{tz}{A timezone to convert dates and times to.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
timezone can be set globally or per function.}
}
\value{
A lubridate PosixCT datetime in the given timezone.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This function is used internally by \code{ruODK} to parse ISO timestamps
to timezone-aware local times.
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachment_get.R
\name{get_one_attachment}
\alias{get_one_attachment}
\title{Download one media attachment.}
\usage{
get_one_attachment(
  pth,
  fn,
  src,
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{pth}{A local file path to save the attachment to.}

\item{fn}{The attachment filename, as per ODK form submission. If NA, no file
will be downloaded, but NA will be returned.
If the file does not exist in ODK Central, a warning will be emitted.}

\item{src}{The attachment's download URL, generated by
\code{\link{attachment_url}}. The src must contain the \verb{uuid:} prefix.
Note that the main Submissions table contains the submission id in the
field \code{id}, whereas nested sub-tables contain the submission id in the
field \code{submissions_id}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The relative local path to the downloaded attachment or NA.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This is a helper function used by \code{\link{attachment_get}}.
This function is not vectorised, but mapped by \code{\link{attachment_get}}
to a tibble of input parameters.
}
\examples{
\dontrun{
# Step 1: Setup ruODK with OData Service URL (has url, pid, fid)
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# Step 2: Construct attachment_url
att_url <- ruODK:::attachment_url(
  "uuid:d3bcefea-32a8-4dbc-80ca-4ecb0678e2b0",
  "filename.jpg"
)

# Step 3: Get one attachment
local_fn <- get_one_attachment("media/filename.jpg", "filename.jpg", att_url)

# In real life: done in bulk behind the scenes during odata_submission_get()
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attachment_get.R
\name{attachment_url}
\alias{attachment_url}
\title{Build the download URL for one or many submission UUIDs and filenames.}
\usage{
attachment_url(
  uuid,
  fn,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url()
)
}
\arguments{
\item{uuid}{The UUID of one form submission, or a vector of UUIDs.}

\item{fn}{The attachment filename, as per ODK form submission, or a vector of
attachment filenames.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}
}
\value{
The inferred download URL.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
This is a helper function used by \code{\link{attachment_get}}.
This function is vectorised and accepts single values or vectors for uuid and
fn.
}
\examples{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

ruODK:::attachment_url(
  "uuid:d3bcefea",
  "filename.jpg",
  pid = 1,
  fid = "form1",
  url = "https://my.odkcentral.org"
)
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{geo_fs}
\alias{geo_fs}
\title{The form_schema of a form containing geofields in GeoJSON.}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 19 rows and 4 columns.
}
\source{
\code{\link{form_schema}}
run on the test form
\code{system.file("extdata", "Locations.xml", package = "ruODK")}.
}
\usage{
geo_fs
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/handle_ru_geotraces.R
\name{handle_ru_geotraces}
\alias{handle_ru_geotraces}
\title{Split all geotraces of a submission tibble into their components.}
\usage{
handle_ru_geotraces(
  data,
  form_schema,
  wkt = FALSE,
  odkc_version = get_default_odkc_version(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{data}{Submissions rectangled into a tibble. E.g. the output of\preformatted{ruODK::odata_submission_get(parse = FALSE) \%>\%
ruODK::odata_submission_rectangle(form_schema = ...)
}}

\item{form_schema}{The \code{form_schema} for the submissions.
E.g. the output of \code{ruODK::form_schema()}.}

\item{wkt}{Whether geofields are GeoJSON (if FALSE) or WKT
strings (if TRUE), default: FALSE.}

\item{odkc_version}{The ODK Central version as decimal number (major.minor).
\code{ruODK} uses this parameter to adjust for breaking changes in ODK Central.

Default: \code{\link{get_default_odkc_version}} or 1.1 if unset.

Set default \code{get_default_odkc_version} through
\code{ru_setup(odkc_version=1.1)}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The submissions tibble with all geotraces retained in their original
format, plus columns of their first point's coordinate components as
provided by \code{\link{split_geotrace}}.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
For a given tibble of submissions, find all columns which are listed
in the form schema as type \code{geotrace}, and extract their components.
Extracted components are longitude (X), latitude (Y), altitude (Z, where
given), and accuracy (M, where given) of the first point of the geotrace.

The original column is retained to allow parsing into other spatially
enabled formats.
}
\examples{
\dontrun{
library(magrittr)
data("geo_fs")
data("geo_wkt_raw")
data("geo_gj_raw")

# GeoJSON
geo_gj_parsed <- geo_gj_raw \%>\%
  ruODK::odata_submission_rectangle(form_schema = geo_fs) \%>\%
  ruODK::handle_ru_geotraces(form_schema = geo_fs, wkt = FALSE)

dplyr::glimpse(geo_gj_parsed)

# WKT
geo_wkt_parsed <- geo_wkt_raw \%>\%
  ruODK::odata_submission_rectangle(form_schema = geo_fs) \%>\%
  ruODK::handle_ru_geotraces(form_schema = geo_fs, wkt = TRUE)

dplyr::glimpse(geo_wkt_parsed)
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_attachments}
\alias{fq_attachments}
\title{A tibble of submission attachments.}
\format{
A tibble of submission attachments.
}
\source{
The output of \code{\link{attachment_list}}
run on submissions of the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}.
}
\usage{
fq_attachments
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_msg.R
\name{ru_msg_noop}
\alias{ru_msg_noop}
\title{Print a green noop message with a filled circle symbol.}
\usage{
ru_msg_noop(message, verbose = get_ru_verbose())
}
\arguments{
\item{message}{(chr) A message to print}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
ru_msg_noop("This is a noop message.")
}
\seealso{
Other messaging: 
\code{\link{ru_msg_abort}()},
\code{\link{ru_msg_info}()},
\code{\link{ru_msg_success}()},
\code{\link{ru_msg_warn}()}
}
\concept{messaging}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odata_service_get.R
\name{odata_service_get}
\alias{odata_service_get}
\title{Retrieve service metadata from an OData URL ending in .svc as tibble.}
\usage{
odata_service_get(
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A tibble with one row per submission data endpoint.
Columns: name, kind, url.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

svc <- odata_service_get()
svc
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/odata-endpoints/odata-form-service/service-document}

Other odata-api: 
\code{\link{odata_metadata_get}()},
\code{\link{odata_submission_get}()}
}
\concept{odata-api}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/form_list.R
\name{form_list}
\alias{form_list}
\title{List all forms.}
\usage{
form_list(
  pid = get_default_pid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A tibble with one row per form and all form metadata as columns.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# With default pid
fl <- form_list()

# With explicit pid
fl <- form_list(pid = 1)

class(fl)
# > c("tbl_df", "tbl", "data.frame")

# Filter out draft forms (published_at=NA)
only_published_forms <- fl \%>\% dplyr::filter(is.na(published_at))

# Note: older ODK Central versions < 1.1 have published_at = NA for both
# published and draft forms. Drafts have NA for version and hash.
only_published_forms <- fl \%>\% dplyr::filter(is.na(version) & is.na(hash))
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/forms}

Other form-management: 
\code{\link{form_detail}()},
\code{\link{form_schema_ext}()},
\code{\link{form_schema}()},
\code{\link{form_xml}()}
}
\concept{form-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_raw}
\alias{fq_raw}
\title{OData submission data for an ODK Central form.}
\format{
A list of lists
}
\source{
See \code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
}
\usage{
fq_raw
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The OData response for the submissions of an ODK Central form.
This form represents a Flora Quadrat, which is a ca 50 by 50 m quadrat of
a uniform plant community.

The XML and .odkbuild versions for this form are available as
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
and \code{system.file("extdata", "FloraQuadrat04.odkbuild", package = "ruODK")},
respectively.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/project_detail.R
\name{project_detail}
\alias{project_detail}
\title{List all details of one project.}
\usage{
project_detail(
  pid = get_default_pid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A tibble with exactly one row for the project and all project
metadata as columns as per ODK Central API docs.
Column names are renamed from ODK's \code{camelCase} to \code{snake_case}.
Values differ to values returned by ODK Central API:
\itemize{
\item archived: FALSE (if NULL) else TRUE
\item dates: NA if NULL
}
}
\description{
While the API endpoint will return all details for one project,
\code{\link{project_detail}} will fail with incorrect or missing
authentication.
}
\details{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

pd <- project_detail()

pd \%>\%
  dplyr::select(-"verbs") \%>\%
  knitr::kable(.)
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/project-management/projects/getting-project-details}

Other project-management: 
\code{\link{project_create}()},
\code{\link{project_list}()}
}
\concept{project-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/submission_get.R
\name{get_one_submission}
\alias{get_one_submission}
\title{Download one submission.}
\usage{
get_one_submission(
  iid,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{iid}{The \code{instance_id}, a UUID, as returned by
\code{\link{submission_list}}.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A nested list of submission data.
}
\description{
This function is the workhorse for the vectorised function submission_get,
which gets all submissions for a list of submission IDs.
}
\details{
Note this function returns a nested list containing any repeating subgroups.
As the presence and length of repeating subgroups is non-deterministic and
entirely depends on the completeness of the submission data, we cannot
rectangle them any further here. Rectangling requires knowledge of the form
schema and the completeness of submission data.

\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

# With explicit credentials, see tests
sl <- submission_list()

sub <- get_one_submission(sl$instance_id[[1]])
listviewer::jsonedit(sub)

# The details for one submission depend on the form fields
length(sub)
# > 11

# The items are the field names. Repeated groups have the same name.
names(sub)
# > "meta"                     "encounter_start_datetime" "reporter"
# > "device_id"                "location"                 "habitat"
# > "vegetation_structure"     "perimeter"                "taxon_encounter"
# > "taxon_encounter"          "encounter_end_datetime"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/submissions/retrieving-submission-xml}

Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_data_strata}
\alias{fq_data_strata}
\title{Parsed submission data for a subgroup of an ODK Central form.}
\format{
The output of \code{\link{odata_submission_get}} for a set of example
data. A tidy tibble referencing the attachments included in the vignettes
and documentation at a relative path \verb{attachments/media/<filename>.<ext>}.
}
\source{
See \code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}
and \code{\link{odata_submission_get}}.
}
\usage{
fq_data_strata
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
The parsed OData response for the subgroup of an ODK Central form.

This subgroup represents vegetation strata as per the NVIS classification.
A vegetation stratum is a layer of plants with the same height, and dominated
by one or few plant taxa. Plant communities can be made of up to five strata,
with two to three being most common.

This data is kept up to date with the data used in vignettes and package
tests. The data is comprised of test records with nonsensical data.
The forms used to capture this data are development versions of real-world
forms.
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/submission_export.R
\name{submission_export}
\alias{submission_export}
\title{Export all form submissions including repeats and attachments to CSV.}
\usage{
submission_export(
  local_dir = here::here(),
  overwrite = TRUE,
  media = TRUE,
  repeats = TRUE,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  pp = get_default_pp(),
  retries = get_retries(),
  odkc_version = get_default_odkc_version(),
  verbose = get_ru_verbose()
)
}
\arguments{
\item{local_dir}{The local folder to save the downloaded files to,
default: \code{here::here}.}

\item{overwrite}{Whether to overwrite previously downloaded zip files,
default: FALSE}

\item{media}{Whether to include media attachments, default: TRUE.
This feature only has effect on ODK Central v1.1 and higher.
Setting this feature to FALSE with an odkc_version < 1.1 and will display a
verbose noop message, but still return all media attachments.}

\item{repeats}{Whether to include repeat data (if TRUE), or whether
to return the root table only (FALSE). Default: TRUE.
Requesting \code{repeats=FALSE} will also omit any media, and override the
parameter \code{media}.
Setting this feature to FALSE with an odkc_version < 1.1 and will display a
verbose noop message, but still include all repeat data.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pp}{The passphrase for an encrypted form.

Default: NULL.

Passphrases can be stored e.g. as environment variables.

See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}

\item{odkc_version}{The ODK Central version as decimal number (major.minor).
\code{ruODK} uses this parameter to adjust for breaking changes in ODK Central.

Default: \code{\link{get_default_odkc_version}} or 1.1 if unset.

Set default \code{get_default_odkc_version} through
\code{ru_setup(odkc_version=1.1)}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The absolute path to the exported ZIP file named after the form ID.
The exported ZIP file will have the extension \code{.zip} unless only the
root table was requested (with \code{repeats=FALSE}), in which case the
exported file will have the extension \code{.csv}.
In contrast to ODK Central, which exports to \code{submissions.csv(.zip)},
the exported ZIP file is named after
the form to avoid accidentally overwriting the ZIP export from
another form.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#maturing}{\figure{lifecycle-maturing.svg}{options: alt='[Maturing]'}}}{\strong{[Maturing]}}
}
\details{
This function exports all the Submission data associated with a Form as one
zip file containing one or more CSV files, as well as all multimedia
attachments associated with the included Submissions.

For an incremental download of a subset of submissions, use
\code{\link{submission_list}} or \code{\link{odata_submission_get}} with
filter queries.
\subsection{Contents}{

The inclusion of subtables (from repeating form groups) can be toggled
through \code{repeats}, whereas the inclusion of media attachments can be toggled
through \code{media}.
}

\subsection{Download location}{

The file will be downloaded to the project root unless specified otherwise
(via \code{local_dir}). Subsequently, the zip file can be extracted.
Attachment filenames (e.g. "12345.jpg") should be prepended with \code{media}
(resulting in e.g. \verb{media/12345.jpg}) in order to represent the relative
path to the actual attachment file (as extracted from the zip file).
}

\subsection{Encryption}{

ODK Central supports two modes of encryption - learn about them
\href{https://odkcentral.docs.apiary.io/#reference/encryption}{here}.
\code{ruODK} supports project managed encryption, however the support is limited
to exactly one encryption key. The supplied passphrase will be used against
the first returned encryption key. Remaining encryption keys are ignored by
\code{ruODK}.

If an incorrect passphrase is given, the request is terminated immediately.
It has been reported that multiple requests with incorrect passphrases
can crash ODK Central.
}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

se <- submission_export()

# Unzip and inspect the loot
t <- tempdir()
f <- unzip(se, exdir = t)
fs::dir_ls(t)
fid <- get_test_fid()
sub <- fs::path(t, glue::glue("{fid}.csv")) \%>\% readr::read_csv()
sub \%>\% knitr::kable(.)
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/submissions/exporting-form-submissions-to-csv}

Other submission-management: 
\code{\link{attachment_list}()},
\code{\link{encryption_key_list}()},
\code{\link{submission_detail}()},
\code{\link{submission_get}()},
\code{\link{submission_list}()}
}
\concept{submission-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/odata_submission_rectangle.R
\name{unnest_all}
\alias{unnest_all}
\title{Recursively unnest_wide all list columns in a tibble.}
\usage{
unnest_all(
  nested_tbl,
  names_repair = "universal",
  names_sep = "_",
  form_schema = NULL,
  verbose = get_ru_verbose()
)
}
\arguments{
\item{nested_tbl}{A nested tibble}

\item{names_repair}{The argument \code{names_repair} for
\code{tidyr::unnest_wider}, default: "universal".}

\item{names_sep}{The argument \code{names_sep} for
\code{tidyr::unnest_wider}, default: "_".
Un-nested variables inside a list column will be prefixed by the list
column name, separated by \code{names_sep}. This avoids unsightly repaired
names such as \code{latitude...1}. Set to \code{NULL} to disable prefixing.}

\item{form_schema}{An optional form_schema,
like the output of \code{\link{form_schema}}. If a form schema is supplied,
location fields will not be unnested. While WKT location fields contain
plain text and will never be unnested, GeoJSON location fields would cause
errors during unnesting.}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\value{
The un-nested tibble in wide format
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\details{
\code{\link{odata_submission_rectangle}} uses this function
internally.
Interested users can use this function to break down \code{ruODK}'s automated
steps into smaller components.

The quite verbose output of \code{tidyr::unnest_wider} is captured
and hidden from the user.
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}}
}
\concept{utilities}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/submission_detail.R
\name{submission_detail}
\alias{submission_detail}
\title{Show metadata for one submission.}
\usage{
submission_detail(
  iid,
  pid = get_default_pid(),
  fid = get_default_fid(),
  url = get_default_url(),
  un = get_default_un(),
  pw = get_default_pw(),
  retries = get_retries()
)
}
\arguments{
\item{iid}{The \code{instance_id}, a UUID, as returned by
\code{\link{submission_list}}.}

\item{pid}{The numeric ID of the project, e.g.: 2.

Default: \code{\link{get_default_pid}}.

Set default \code{pid} through \code{ru_setup(pid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{fid}{The alphanumeric form ID, e.g. "build_Spotlighting-0-8_1559885147".

Default: \code{\link{get_default_fid}}.

Set default \code{fid} through \code{ru_setup(fid="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{url}{The ODK Central base URL without trailing slash.

Default: \code{\link{get_default_url}}.

Set default \code{url} through \code{ru_setup(url="...")}.

See \code{vignette("Setup", package = "ruODK")}.}

\item{un}{The ODK Central username (an email address).
Default: \code{\link{get_default_un}}.
Set default \code{un} through \code{ru_setup(un="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{pw}{The ODK Central password.
Default: \code{\link{get_default_pw}}.
Set default \code{pw} through \code{ru_setup(pw="...")}.
See \code{vignette("Setup", package = "ruODK")}.}

\item{retries}{The number of attempts to retrieve a web resource.

This parameter is given to \code{\link[httr]{RETRY}(times = retries)}.

Default: 3.}
}
\value{
A nested list of submission metadata.
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
# See vignette("setup") for setup and authentication options
# ruODK::ru_setup(svc = "....svc", un = "me@email.com", pw = "...")

sl <- submission_list()

sub <- submission_detail(sl$instance_id[[1]])

# The details for one submission return exactly one row
nrow(sub)
# > 1

# The columns are metadata about the submission and the submitter
names(sub)
# > "instance_id" "submitter_id" "submitter_type" "submitter_display_name"
# > "submitter_created_at" "device_id" "created_at"
}
}
\seealso{
\url{https://odkcentral.docs.apiary.io/#reference/forms-and-submissions/submissions/getting-submission-details}

Other submission-management: 
\code{\link{attachment_list}()},
\code{\link{encryption_key_list}()},
\code{\link{submission_export}()},
\code{\link{submission_get}()},
\code{\link{submission_list}()}
}
\concept{submission-management}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fs_v7}
\alias{fs_v7}
\title{The parsed XML form_schema of a form from ODK Central v0.6.}
\format{
An object of class \code{tbl_df} (inherits from \code{tbl}, \code{data.frame}) with 12 rows and 3 columns.
}
\source{
\code{\link{form_schema_parse}(fs_v7_raw)}
}
\usage{
fs_v7
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{fq_project_detail}
\alias{fq_project_detail}
\title{A tibble of project metadata.}
\format{
A tibble of project metadata.
}
\source{
The output of \code{\link{project_detail}}
run on the project containing the test form
\code{system.file("extdata", "FloraQuadrat04.xml", package = "ruODK")}.
}
\usage{
fq_project_detail
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fq_project_list}
\alias{fq_project_list}
\title{A tibble of project metadata.}
\format{
A tibble of project metadata.
}
\source{
The output of \code{\link{project_list}}
run on all projects on the configured ODK Central server.
}
\usage{
fq_project_list
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7_raw}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/drop_null_coords.R
\name{drop_null_coords}
\alias{drop_null_coords}
\title{Drop any NULL coordinates from a GeoJSON geometry.}
\usage{
drop_null_coords(x)
}
\arguments{
\item{x}{A GeoJSON geometry parsed as nested list.
E.g. \code{geo_gj$path_location_path_gps}.}
}
\value{
The nested list minus the last element (if NULL).
}
\description{
This helper patches a bug/feature in ODK Central (versions 0.7-0.9), where
geotrace / geoshape GeoJSON contains a last coordinate pair with NULL
lat/lon (no alt/acc), and WKT ends in \verb{, undefined NaN}.
}
\details{
While \code{\link{split_geotrace}} and \code{\link{split_geoshape}} modify
the WKT inline, it is more maintainable to separate the GeoJSON cleaner
into this function.

This helper drops the last element of a GeoJSON coordinate list if it is
\code{list(NULL, NULL)}.
}
\examples{
# A snapshot of geo data with trailing empty coordinates.
data("geo_gj88")

len_coords <- length(geo_gj88$path_location_path_gps[[1]]$coordinates)

length(geo_gj88$path_location_path_gps[[1]]$coordinates[[len_coords]]) \%>\%
  testthat::expect_equal(2)

geo_gj88$path_location_path_gps[[1]]$coordinates[[len_coords]][[1]] \%>\%
  testthat::expect_null()

geo_gj88$path_location_path_gps[[1]]$coordinates[[len_coords]][[2]] \%>\%
  testthat::expect_null()

# The last coordinate pair is a list(NULL, NULL).
# Invalid coordinates like these are a choking hazard for geospatial
# packages. We should remove them before we can convert ODK data into native
# spatial formats, such as sf.
str(geo_gj88$path_location_path_gps[[1]]$coordinates[[len_coords]])

geo_gj_repaired <- geo_gj88 \%>\%
  dplyr::mutate(
    path_location_path_gps = path_location_path_gps \%>\%
      purrr::map(drop_null_coords)
  )

len_coords_repaired <- length(
  geo_gj_repaired$path_location_path_gps[[1]]$coordinates
)
testthat::expect_equal(len_coords_repaired + 1, len_coords)
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{predict_ruodk_name}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ru_msg.R
\name{ru_msg_warn}
\alias{ru_msg_warn}
\title{rlang::warn() with a yellow warning message with a warning symbol.}
\usage{
ru_msg_warn(message, verbose = get_ru_verbose())
}
\arguments{
\item{message}{(chr) A message to print}

\item{verbose}{Whether to display debug messages or not.

Read \code{vignette("setup", package = "ruODK")} to learn how \code{ruODK}'s
verbosity can be set globally or per function.}
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\examples{
\dontrun{
ru_msg_warn("This is a warning.")
}
}
\seealso{
Other messaging: 
\code{\link{ru_msg_abort}()},
\code{\link{ru_msg_info}()},
\code{\link{ru_msg_noop}()},
\code{\link{ru_msg_success}()}
}
\concept{messaging}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\encoding{UTF-8}
\name{fs_v7_raw}
\alias{fs_v7_raw}
\title{The unparsed XML form_schema of a form from ODK Central v0.6 as nested list.}
\format{
An object of class \code{list} of length 6.
}
\source{
\code{\link{form_schema}(odkc_version = 0.7, parse = FALSE)}
}
\usage{
fs_v7_raw
}
\description{
\ifelse{html}{\href{https://lifecycle.r-lib.org/articles/stages.html#stable}{\figure{lifecycle-stable.svg}{options: alt='[Stable]'}}}{\strong{[Stable]}}
}
\seealso{
Other included: 
\code{\link{fq_attachments}},
\code{\link{fq_data_strata}},
\code{\link{fq_data_taxa}},
\code{\link{fq_data}},
\code{\link{fq_form_detail}},
\code{\link{fq_form_list}},
\code{\link{fq_form_schema}},
\code{\link{fq_form_xml}},
\code{\link{fq_meta}},
\code{\link{fq_project_detail}},
\code{\link{fq_project_list}},
\code{\link{fq_raw_strata}},
\code{\link{fq_raw_taxa}},
\code{\link{fq_raw}},
\code{\link{fq_submission_list}},
\code{\link{fq_submissions}},
\code{\link{fq_svc}},
\code{\link{fq_zip_data}},
\code{\link{fq_zip_strata}},
\code{\link{fq_zip_taxa}},
\code{\link{fs_v7}},
\code{\link{geo_fs}},
\code{\link{geo_gj88}},
\code{\link{geo_gj_raw}},
\code{\link{geo_gj}},
\code{\link{geo_wkt88}},
\code{\link{geo_wkt_raw}},
\code{\link{geo_wkt}}
}
\concept{included}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/form_schema_parse.R
\name{predict_ruodk_name}
\alias{predict_ruodk_name}
\title{Predict a field name after \code{tidyr::unnest_wider(names_sep="_")} prefixes
the form path.}
\usage{
predict_ruodk_name(name_str, path_str)
}
\arguments{
\item{name_str}{An Xforms field name string.}

\item{path_str}{A path string,
e.g. "Submissions" or "Submissions.group_name".}
}
\value{
The name as built by \code{tidyr::unnest_wider(names_sep="_")}.
}
\description{
Predict a field name after \code{tidyr::unnest_wider(names_sep="_")} prefixes
the form path.
}
\examples{
\dontrun{
predict_ruodk_name("bar", "Submissions.foo")
# > "foo_bar"
predict_ruodk_name("bar", "Submissions")
# > "bar"
predict_ruodk_name("rock", "Submissions.foo_fighters")
# > "foo_fighters_rock"
}
}
\seealso{
Other utilities: 
\code{\link{attachment_get}()},
\code{\link{attachment_link}()},
\code{\link{attachment_url}()},
\code{\link{drop_null_coords}()},
\code{\link{form_schema_parse}()},
\code{\link{get_one_attachment}()},
\code{\link{get_one_submission_attachment_list}()},
\code{\link{get_one_submission}()},
\code{\link{handle_ru_attachments}()},
\code{\link{handle_ru_datetimes}()},
\code{\link{handle_ru_geopoints}()},
\code{\link{handle_ru_geoshapes}()},
\code{\link{handle_ru_geotraces}()},
\code{\link{isodt_to_local}()},
\code{\link{odata_submission_rectangle}()},
\code{\link{prepend_uuid}()},
\code{\link{split_geopoint}()},
\code{\link{split_geoshape}()},
\code{\link{split_geotrace}()},
\code{\link{strip_uuid}()},
\code{\link{tidyeval}},
\code{\link{unnest_all}()}
}
\concept{utilities}
\keyword{internal}
