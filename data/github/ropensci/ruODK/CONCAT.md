
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

# References