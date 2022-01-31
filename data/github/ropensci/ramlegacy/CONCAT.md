
<!-- README.md is generated from README.Rmd. Please edit that file -->

# ramlegacy

[![CRAN
status](https://www.r-pkg.org/badges/version/ramlegacy)](https://cran.r-project.org/package=ramlegacy)
[![Travis Build
Status](https://travis-ci.com/ropensci/ramlegacy.svg?branch=master)](https://travis-ci.com/ropensci/ramlegacy)
[![AppVeyor Build
Status](https://ci.appveyor.com/api/projects/status/github/ropensci/ramlegacy?branch=master&svg=true)](https://ci.appveyor.com/project/kshtzgupta1/ramlegacy)
[![Coverage
status](https://codecov.io/gh/ropensci/ramlegacy/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/ramlegacy)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![](https://badges.ropensci.org/264_status.svg)](https://github.com/ropensci/software-review/issues/264)

  - **Authors**: Kshitiz Gupta, [Carl
    Boettiger](https://www.carlboettiger.info/)
  - **License**: [MIT](https://opensource.org/licenses/MIT)
  - [Package source code on
    Github](https://github.com/ropensci/ramlegacy)
  - [**Submit Bugs and feature
    requests**](https://github.com/ropensci/ramlegacy/issues)

`ramlegacy` is an R package that supports caching and reading in
different versions of the RAM Legacy Stock Assessment Data Base, an
online compilation of stock assessment results for commercially
exploited marine populations from around the world. More information
about the database can be found [here.](https://www.ramlegacy.org/)

## What does `ramlegacy` do?

  - Provides a function `download_ramlegacy()`, to download all the
    available versions of the RAM Legacy Stock Assessment Excel Database
    and cache them on the user’s computer as serialized RDS objects.
    This way once a version has been downloaded it doesn’t need to be
    re-downloaded for subsequent analysis.
  - Supports reading in specified tables or all tables from a cached
    version of the database through a function `load_ramlegacy()`
  - Provides a function `ram_dir()` to view the path of the location
    where the downloaded database was cached.

## Installation

You can install the development version from
[Github](https://github.com/ropensci/ramlegacy) with:

``` r
install.packages("devtools")
library(devtools)
install_github("ropensci/ramlegacy")
```

To ensure that the vignette is installed along with the package make
sure to remove `--no-build-vignettes` from the `build_opts` in
`install_github`

## Usage

Please see the ramlegacy
[vignette](https://docs.ropensci.org/ramlegacy/articles/ramlegacy.html)
for more detailed examples and additional package functionality.

Start by loading the package using `library`.

``` r
library(ramlegacy)
```

### download\_ramlegacy

`download_ramlegacy()` downloads the specified version of **RAM Legacy
Stock Assessment Excel Database** and then saves it as an RDS object in
user’s application data directory as detected by the
[rappdirs](https://CRAN.R-project.org/package=rappdirs) package. This
location is also where `load_ramlegacy()` by default will look for the
downloaded database.

``` r
# downloads version 4.44
download_ramlegacy(version = "4.44")
```

If version is not specified then `download_ramlegacy` defaults to
downloading current latest version (4.44) :

``` r
# downloads current latest version 4.44
download_ramlegacy()
```

The latest versions of the RAM Legacy Database are [archived in
Zenodo](https://zenodo.org/communities/rlsadb/) but the older versions
(v4.3, v3.0, v2.5, v2.0, v1.0) are not. To ensure access to these older
versions of the database `download_ramlegacy` supports downloading them
from this [Github
repository](https://www.github.com/kshtzgupta1/ramlegacy-assets/):

``` r
# downloads older version 4.3
download_ramlegacy(version = "4.3")
```

### load\_ramlegacy

After the specified version of the database has been downloaded and
cached on your local machine through `download_ramlegacy` you can call
`load_ramlegacy` to obtain a list of specific tables/all the tables from
that version of the database. If version is not specified but tables is
then `load_ramlegacy` defaults to returning a list containing the
specified dataframes from the latest version (currently 4.44). If both
version and tables are not specified then `load_ramlegacy` defaults to
returning a list containing all the dataframes in the latest version
(currently 4.44)

``` r
# get a list containing area and bioparams tables from
# version 4.3 of the database
load_ramlegacy(version = "4.3", tables = c("area", "bioparams"))

# get a list containing area and bioparams tables from version 4.44
# of the database
load_ramlegacy(version = "4.44", tables = c("area", "bioparams"))

# if tables is specified but version is not then the function defaults
# to returning a list containing the specified tables from the current
# latest version 4.44
load_ramlegacy(tables = c("area", "bioparams"))

# since both tables and version are not specified the function returns
# a list containing all the tables from the current latest version 4.44
load_ramlegacy()
```

To learn more about the different tables present in the database, what
the various acronyms mean and the different stock summaries accompanying
the databases please see this
[page.](https://docs.ropensci.org/ramlegacy/articles/tables_description.html)

### ram\_dir

To view the exact path where a certain version of the database was
downloaded and cached by `download_ramlegacy` you can run `ram_dir(vers
= 'version')`, specifying the version number inside the function call:

``` r
# download version 4.44
download_ramlegacy(version = "4.44")

# view the location where version 4.44 of the database was
# downloaded and cached
ram_dir(vers = "4.44")
```

## Similar Projects

1.  [`ramlegacy`](https://github.com/seananderson/ramlegacy) Sean
    Anderson has a namesake package that appears to be a stalled project
    on Github (last updated 9 months ago). However, unlike this package
    which supports downloading and reading in the Excel version of the
    database, Sean Anderson’s project downloads the Microsoft Access
    version and converts it to a local sqlite3 database.

2.  [`RAMlegacyr`](https://github.com/ashander/RAMlegacyr) `RAMlegacyr`
    is an older package last updated in 2015. Similar to Sean Anderson’s
    project, the package seems to be an R interface for the Microsoft
    Access version of the RAM Legacy Stock Assessment Database and
    provides a set of functions using RPostgreSQL to connect to the
    database.

## Citation

Current and older versions of the RAM Legacy Database are [archived in
Zenodo](https://zenodo.org/communities/rlsadb/), each version with its
own unique DOI. The suggested format for citing data is:

RAM Legacy Stock Assessment Database. 2018. Version
4.44-assessment-only. Released 2018-12-22. Accessed \[Date accessed
YYYY-MM-DD\]. Retrieved from
[DOI:10.5281/zenodo.2542919.](https://zenodo.org/record/2542919#.XE-rFs9KjBI)

The primary publication describing the RAM Legacy Stock Assessment
Database, and suggested citation for general use is:

Ricard, D., Minto, C., Jensen, O.P. and Baum, J.K. (2012) Evaluating the
knowledge base and status of commercially exploited marine species with
the RAM Legacy Stock Assessment Database. Fish and Fisheries 13 (4)
380-398.
[DOI: 10.1111/j.1467-2979.2011.00435.x](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-2979.2011.00435.x)

Several [publications](https://sites.uw.edu/ramlegac/publications/) have
relied on the RAM Legacy Stock Assessment
Database.

[![ropensci\_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# ramlegacy 0.1.0

* Initial Release

# ramlegacy 0.2.0

* Removed the loading behavior of `library(ramlegacy)` leaving the loading solely to `load_ramlegacy`
* Instead of assigning all the tables to global environment `load_ramlegacy` now just returns a list of tables that can be specified through `tables` argument.
* Modified `download_ramlegacy` and `load_ramlegacy` so that that while the default is still to download and read in the dataframes from the rappdirs directory the functions now also support downloading to a location chosen by the user and reading from that location. 

## For more fine-grained list of changes or to report a bug, consult 

* [The issues log](https://github.com/ropensci/ramlegacy/issues)
* [The commit log](https://github.com/ropensci/ramlegacy/commits/master)

## Versioning Guidelines

Releases will be numbered with the following semantic versioning format:

`<major>.<minor>.<patch>`

And constructed with the following guidelines:

* Breaking backward compatibility bumps the major (and resets the minor 
  and patch)
* New additions without breaking backward compatibility bumps the minor 
  (and resets the patch)
* Bug fixes and misc changes bumps the patch
* Following the RStudio convention, a .99 is appended after the patch
  number to indicate the development version on Github.  Any version
  Coming from Github will now use the .99 extension, which will never
  appear in a version number for the package on CRAN. 
# Contributor Code of Conduct

As contributors and maintainers of this project, we pledge to respect all people who 
contribute through reporting issues, posting feature requests, updating documentation,
submitting pull requests or patches, and other activities.

We are committed to making participation in this project a harassment-free experience for
everyone, regardless of level of experience, gender, gender identity and expression,
sexual orientation, disability, personal appearance, body size, race, ethnicity, age, or religion.

Examples of unacceptable behavior by participants include the use of sexual language or
imagery, derogatory comments or personal attacks, trolling, public or private harassment,
insults, or other unprofessional conduct.

Project maintainers have the right and responsibility to remove, edit, or reject comments,
commits, code, wiki edits, issues, and other contributions that are not aligned to this 
Code of Conduct. Project maintainers who do not follow the Code of Conduct may be removed 
from the project team.

Instances of abusive, harassing, or otherwise unacceptable behavior may be reported by 
opening an issue or contacting one or more of the project maintainers.

This Code of Conduct is adapted from the Contributor Covenant 
(https:contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
# Contributing to `ramlegacy`

Thank you for any and all contributions! Following these guidelines will help streamline the process of contributing and make sure that we're all on the same page. While we ask that you read this guide and follow it to the best of your abilities, we welcome contributions from all, regardless of your level of experience.

By participating in this project, you agree to abide by the [code of conduct](https://github.com/ropensci/ramlegacy/blob/master/CONDUCT.md).

# Types of contributions 

Following types of contributions are always welcome:

- Identify areas for future development ([open an Issue](https://github.com//ramlegacy/issues))
- Identify issues/bugs ([open an Issue](https://github.com/ropensci/ramlegacy/issues))
- Write tutorials/vignettes ([open a Pull Request](https://github.com/ropensci/ramlegacy/pulls) to contribute to the ones here, or make your own elsewhere and send us a link)
- Add functionality ([open a Pull Request](https://github.com/ropensci/ramlegacy/pulls))
- Fix bugs ([open a Pull Request](https://github.com/ropensci/ramlegacy/pulls))


# How to contribute code

- Fork the repository
- Clone the repository from GitHub to your computer e.g,. `git clone https://github.com/ropensci/ramlegacy.git`
- Make sure to track progress upstream (i.e., on our version of `ramlegacy` at `ropensci/ramlegacy`)
  - `git remote add upstream https://github.com/ropensci/ramlegacy.git`
  - Before making changes make sure to pull changes in from upstream with `git pull upstream`
- Make your changes
  - For changes beyond minor typos, add an item to NEWS.md describing the changes and add yourself to the DESCRIPTION file as a contributor
- Push to your GitHub account
- Submit a pull request to home base (likely master branch, but check to make sure) at `ropensci/ramlegacy`

# Code formatting

- In general follow the convention of <https://r-pkgs.had.co.nz/r.html#style> (snake_case functions and argument names, etc.)
- Use explicit package imports (i.e. package_name::package_function) and avoid @import if at all possible
## Resubmission
This is a resubmission. In this version I have:

* Corrected the database reference in DESCRIPTION to be in the form: authors (year) <doi:....>
* Put Excel within quotes in DESCRIPTION
* Explained the acronym RAM in DESCRIPTION
* Ensured that functions do not write by default in the user's home filespace.
* Ensured that functions in examples/vignettes/tests write to tempdir(). Please note that the code
  where we explicitly set the destination file path to `tempdir()` is within `\dontshow{}` for the 
  sake of having clear unambiguous examples.

## Test environments
* Local - ubuntu 18.04 (R 3.6.0)
* Travis CI - ubuntu 14.04 (R 3.6.0)
* Appveyor - Windows Server 2012 R2 x64 (R 3.6.0 (patched))
* win-builder (oldrelease, release, and devel)

## Note about examples in download_ramlegacy and load_ramlegacy
Because of the large size of the RAM Legacy database the examples in `download_ramlegacy` and `load_ramlegacy` take longer than 5 seconds to run and that's why there were placed within `donttest{}`. The runtime for examples in `download_ramlegacy` is 286.102 seconds and the runtime for examples in `load_ramlegacy` is 94.094 seconds.

## R CMD check results

There were no ERRORs, no WARNINGs, no NOTEs.

## Downstream dependencies

There are no downstream dependencies for this package.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "man/figures/README-",
  out.width = "100%"
)
```
ramlegacy
=============
[![CRAN status](https://www.r-pkg.org/badges/version/ramlegacy)](https://cran.r-project.org/package=ramlegacy)
[![Travis Build Status](https://travis-ci.com/ropensci/ramlegacy.svg?branch=master)](https://travis-ci.com/ropensci/ramlegacy) [![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/ropensci/ramlegacy?branch=master&svg=true)](https://ci.appveyor.com/project/kshtzgupta1/ramlegacy) [![Coverage status](https://codecov.io/gh/ropensci/ramlegacy/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/ramlegacy) 
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
 [![](https://badges.ropensci.org/264_status.svg)](https://github.com/ropensci/software-review/issues/264)



- **Authors**: Kshitiz Gupta, [Carl Boettiger](https://www.carlboettiger.info/)
- **License**: [MIT](https://opensource.org/licenses/MIT)
- [Package source code on Github](https://github.com/ropensci/ramlegacy)
- [**Submit Bugs and feature requests**](https://github.com/ropensci/ramlegacy/issues)

`ramlegacy` is an R package that supports caching and reading in different versions of the RAM Legacy Stock Assessment Data Base, an online compilation of stock assessment results for commercially exploited marine populations from around the world. More information about the database can be found [here.](https://www.ramlegacy.org/)


## What does `ramlegacy` do?

  - Provides a function `download_ramlegacy()`, to download all the available versions
    of the RAM Legacy Stock Assessment Excel Database and cache them on the user's computer
    as serialized RDS objects. This way once a version has been downloaded it doesn't need
    to be re-downloaded for subsequent analysis.
  - Supports reading in specified tables or all tables from a cached version of the database through a 
    function `load_ramlegacy()`
  - Provides a function `ram_dir()` to view the path of the location where the downloaded database was cached.

## Installation
You can install the development version from [Github](https://github.com/ropensci/ramlegacy) with:

```{r, eval = FALSE, echo = TRUE}
install.packages("devtools")
library(devtools)
install_github("ropensci/ramlegacy")
```

To ensure that the vignette is installed along with the package make sure to remove `--no-build-vignettes` from the `build_opts` in `install_github`

## Usage
Please see the ramlegacy [vignette](https://docs.ropensci.org/ramlegacy/articles/ramlegacy.html) for more detailed examples and additional package functionality. 

Start by loading the package using `library`.

```{r load, echo = TRUE, eval = FALSE}
library(ramlegacy)
```

### download_ramlegacy

`download_ramlegacy()` downloads the specified version of **RAM Legacy Stock Assessment Excel Database** and then saves it as an RDS object in user’s application data directory as detected by the [rappdirs](https://CRAN.R-project.org/package=rappdirs) package. This location is also where `load_ramlegacy()` by default will look for the downloaded database. 

```{r, download_ramlegacy example1, eval = F, echo = T}
# downloads version 4.44
download_ramlegacy(version = "4.44")
```

If version is not specified then `download_ramlegacy` defaults to downloading current latest version (4.44) :

```{r, download_ramlegacy_example2, eval = F, echo = T}
# downloads current latest version 4.44
download_ramlegacy()
```

The latest versions of the RAM Legacy Database are [archived in Zenodo](https://zenodo.org/communities/rlsadb/) but the older versions (v4.3, v3.0, v2.5, v2.0, v1.0) are not. To ensure access to these older versions of the database `download_ramlegacy` supports downloading them from this [Github repository](https://www.github.com/kshtzgupta1/ramlegacy-assets/):

```{r, download_ramlegacy_example3, eval = F, echo = T}
# downloads older version 4.3
download_ramlegacy(version = "4.3")
```

### load_ramlegacy
After the specified version of the database has been downloaded and cached on your local machine through `download_ramlegacy` you can call `load_ramlegacy` to obtain a list of specific tables/all the tables from that version of the database. If version is not specified but tables is then `load_ramlegacy` defaults to returning a list containing the specified dataframes from the latest version (currently 4.44). If both version and tables are not specified then `load_ramlegacy` defaults to returning a list containing all the dataframes in the latest version (currently 4.44)

```{r, load_ramlegacy example1, eval = F, echo = T}
# get a list containing area and bioparams tables from
# version 4.3 of the database
load_ramlegacy(version = "4.3", tables = c("area", "bioparams"))

# get a list containing area and bioparams tables from version 4.44
# of the database
load_ramlegacy(version = "4.44", tables = c("area", "bioparams"))

# if tables is specified but version is not then the function defaults
# to returning a list containing the specified tables from the current
# latest version 4.44
load_ramlegacy(tables = c("area", "bioparams"))

# since both tables and version are not specified the function returns
# a list containing all the tables from the current latest version 4.44
load_ramlegacy()
```

To learn more about the different tables present in the database, what the various acronyms mean and the different stock summaries accompanying the databases please see this [page.](https://docs.ropensci.org/ramlegacy/articles/tables_description.html)


### ram_dir
To view the exact path where a certain version of the database was downloaded and cached by `download_ramlegacy` you can run `ram_dir(vers = 'version')`, specifying the version number inside the function call:
```{r, ram_dir_example1, eval = F, echo = T}
# download version 4.44
download_ramlegacy(version = "4.44")

# view the location where version 4.44 of the database was
# downloaded and cached
ram_dir(vers = "4.44")
```

## Similar Projects

1. [`ramlegacy`](https://github.com/seananderson/ramlegacy) 
Sean Anderson has a namesake package that appears to be a stalled project on Github (last updated 9 months ago). However, unlike this package which supports downloading and reading in the Excel version of the database, Sean Anderson's project downloads the Microsoft Access version and converts it to a local sqlite3 database. 
 

2. [`RAMlegacyr`](https://github.com/ashander/RAMlegacyr)
`RAMlegacyr` is an older package last updated in 2015. Similar to Sean Anderson's project, the package seems to be an R interface for the Microsoft Access version of the RAM Legacy Stock Assessment Database and provides a set of functions using RPostgreSQL to connect to the database.

## Citation
Current and older versions of the RAM Legacy Database are [archived in Zenodo](https://zenodo.org/communities/rlsadb/), each version with its own unique DOI. The suggested format for citing data is:

RAM Legacy Stock Assessment Database. 2018. Version 4.44-assessment-only. Released 2018-12-22. Accessed [Date accessed YYYY-MM-DD]. Retrieved from [DOI:10.5281/zenodo.2542919.](https://zenodo.org/record/2542919#.XE-rFs9KjBI)

The primary publication describing the RAM Legacy Stock Assessment Database, and suggested citation for general use is:

Ricard, D., Minto, C., Jensen, O.P. and Baum, J.K. (2012) Evaluating the knowledge base and status of commercially exploited marine species with the RAM Legacy Stock Assessment Database. Fish and Fisheries 13 (4) 380-398. [DOI: 10.1111/j.1467-2979.2011.00435.x](https://onlinelibrary.wiley.com/doi/abs/10.1111/j.1467-2979.2011.00435.x)

Several [publications](https://sites.uw.edu/ramlegac/publications/) have relied on the RAM Legacy Stock Assessment Database.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
---
title: 'Tables and Stock Summary Description'
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: |
  %\VignetteIndexEntry{tables_description}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## Description of the dataframes present in the database:

* metadata: Table with summarized metadata (only available in newer versions starting from v4.40)
* stock: This stores the stock database table
* assessment: This stores the assessment database table
* taxonomy: This stores the taxonomy database table
* management: This stores the management database table
* assessor: This stores the assessor database table
* assessmetod: This stores the assessmetod database table
* area: This stores the area database table
* biometrics: This stores the biometrics database table
* tsmetrics: This stores the tsmetrics database table
* timeseries: The time series data is a matrix object with the following headers/columns: (1) assessid (2) stockid (3) stocklong (4) tsid (5) tsyear (6) tsvalue
* bioparams: The parameter data is a matrix object with the following headers/columns: (1) assessid (2) stockid (3) stocklong (4) bioid (5) biovalue (6) bioyear (7) bionotes
* timeseries_values_views: This stores the timeseries values with timeseries type along the columns and stocks along the rows
* timeseries_units_views: This stores the timeseries values with timeseries type along the columns and stocks along the rows
* timeseries_ids_views: This stores the timeseries IDs with timeseries type along the columns and stocks along the rows
* timeseries_assessments_views: This stores the timeseries assessments with timeseries type along the columns and stocks along the rows
* timeseries_notes_views: This stores the timeseries notes with timeseries type along the columns and stocks along the rows
* timeseries_sources_views: This stores the timeseries sources with timeseries type along the columns and stocks along the rows
* timeseries_years_views: This stores the timeseries years with timeseries type along the columns and stocks along the rows
* bioparams_values_views: This stores the reference point values, with reference point type along the columns and stocks along the rows
* bioparams_units_views: This stores the reference point units, with reference point type along the columns and stocks along the rows
* bioparams_ids_views: This stores the reference point IDs, with reference point type along the columns and stocks along the rows
* bioparams_assessments_views: This stores the reference point assessments, with reference point type along the columns and stocks along the rows
* bioparams_sources_views: This stores the reference point sources, with reference point type along the columns and stocks along the rows
* bioparams_notes_views: This stores the reference point notes, with reference point type along the columns and stocks along the rows

## Newer Versions (v4.40 onwards) also contain tables of individual most-used time series:

* tb.data: Total Biomass
* ssb.data: Spawning Stock Biomass
* tn.data: Total Abundance
* r.data: Recruits
* tc.data: Total Catch
* tl.data: Total Landings
* recc.data: Recreational Catch
* f.data: Fishing Mortality
* er.data: Exploitation Rate
* divtb.data: TB/TBmsy
* divssb.data: SSB/SSBmsy
* ivf.data: F/Fmsy
* diver.data: ER/ERmsy
* divbpref.data: B/Bmsypref
* divupref.data: U/Umsypref
* tbbest.data: TBbest
* tcbest.data: TCbest
* erbest.data: ERbest
* divtb.mgt.data: TB/TBmgt
* divssb.mgt.data: SSB/SSBmgt
* divf.mgt.data: F/Fmgt
* diver.mgt.data: ER/ERmgt
* divbpref.mgt.data: B/Bmgtpref
* divupref.mgt.data: U/Umgtpref
* cpair.data: Cpair
* tac.data: TAC
* cadv.data: Cadvised
* survb.data: survB
* cpue.data: CPUE
* effort.data: EFFORT


## Stock Summary Description

### Information Contained in the Stock Summary Documents

* A list of all stocks present in the region
    + specifies whether the stock includes B/Bmsy, U/Umsy, B/Bmgt, and U/Umgt time series data
    + specifies whether TB or SSB is preferred, and whether ER or F is preferred, if both types are present

### Information available for each stock:

* Metadata (scientific name, area, management authority, assessor, assessment years)
* Reference points available for stock
* Time series available for stock
* Plots – based on available data, may include:
    +  Kobe plot (MSYpref) – Plots U/Umsy vs. B/Bmsy if available, filling in gaps with U/Umgt and B/Bmgt if         available
    + Kobe plot (MGTpref) – Plots U/Umgt vs. B/Bmgt if available, filling in gaps with U/Umsy and B/Bmsy if          available
    + Spawner Recruitment – Plots recruits vs. spawning-stock biomass
    + Surplus Production – Plots annual surplus production vs. total biomass
    + Total Biomass – Preferably plots the most recent TB time series from an assessment with a MSY,                management target, or limit reference point available, otherwise plots most recent TB time series.
    + Spawning stock biomass – Preferably plots the most recent SSB time series from an assessment with a MSY,       management target, or limit reference point available, otherwise plots most recent SSB time series.
    + Total abundance (in numbers) – Preferably plots the most recent TN time series from an assessment with a       MSY, management target, or limit reference point available, otherwise plots most recent TN time series.
    + Fishing mortality (instantaneous rate) – Preferably plots the most recent F time series from an               assessment with a MSY, management target, or limit reference point available, otherwise plots most            recent F time series.
    + Exploitation rate (annual proportion) – Preferably plots the most recent ER time series from an               assessment with a MSY, management target, or limit reference point available, otherwise plots most            recent ER time series.
    + Recruits – Plots the most recent recruitment time series. If known, the year range during which               recruitment deviates were estimated in the assessment model is shaded grey.
    + TC, TL, RecC – Plots most recent catch and/or landings as well as recreational catch  and MSY if they         are also available from that assessment. (The year range of moratoriums are shaded grey if available.)
    + TAC, Cpair, Cadv – Plots most recent total allowable catch or landings (or other measures of total            quota), the corresponding catch or landings, and the corresponding scientifically advised catch if            available from the same assessment.
    + survB – Plots the most recent fishery-independent survey biomass time series.
    + CPUE – Plots the most recent fishery-independent catch-per-unit-effort time series.
    + EFFORT – Plots the most recent fishing effort time series.
    + CdivMSY – Plots the ratio of the most recent catch time series to MSY. MSY is drawn from the assessment       directly if available, or else is calculated from Bmsy and Umsy reference points (including proxies) if       available.

### Acronym and Abbreviation Definitions

* B – Biomass; may be either total biomass or spawning stock biomass
* Bmgt – Biomass at management target reference point
* Bmsy – Biomass at MSY
* Cadv – Scientifically advised catch
* CdivMSY – Catch divided by MSY
* Cpair – Catch corresponding to TAC and Cadv
* CPUE – Catch per unit effort
* EFFORT – Measure of fishing effort (depends on fishery)
* ER – Exploitation rate (annual proportion)
* F – Fishing mortality (instantaneous rate)
* MSY – Maximum Sustainable Yield
* RecC – Recreational catch
* SSB – Spawning stock biomass
* survB – Survey biomass
* TAC – Total allowable catch
* TB – Total biomass
* TC – Total catch 
* TL – Total landings
* TN – Total abundance; in numbers
* U – Harvest rate; may be either exploitation rate or fishing mortality
* Umsy – Harvest rate at MSY
* Umgt – Harvest rate at management target reference point
---
title: 'Introduction to ramlegacy'
author: "Carl Boettiger & Kshitiz Gupta" 
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: |
  %\VignetteIndexEntry{ramlegacy}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

## RAM Legacy Stock Assessment Database
From the database's [website](www.ramlegacy.org):

> The RAM Legacy Stock Assessment Database is a compilation of stock assessment results for commercially exploited marine populations from around the world. It is inspired by Dr. Ransom A. Myers' original [stock-recruitment database](http://wayback.archive-it.org/org-990/20170328023224/http://www.mscs.dal.ca/~myers/welcome.html), which is no longer being updated."

As of April 2019, the available versions are 1.0, 2.0, 2.5, 3.0, 4.3, 4.40, 4.41, 4.44.

The latest versions of the database (4.44, 4.41, 4.40) are available at [Zenodo.](https://zenodo.org/communities/rlsadb/) This is where `download_ramlegacy` downloads the database from.

The older versions of the database (4.3, 3.0, 2.5, 2.0, 1.0) are available in **Excel** (Microsoft Excel Open XML Format Spreadsheet) format at this [github repo](https://github.com/kshtzgupta1/ramlegacy-assets).

## About the package
The goal of `ramlegacy` is to cache a specified version(s) of the RAM Legacy Stock Assessment Excel Database and this way save time and effort spent in re-downloading different versions of the database as part of any future analysis involving the data. The package also supports reading in any specified dataframe(s) of the database through a function called `load_ramlegacy()`. The user can always view the location to which the database is cached and read from by running `ram_dir()` and specifying its version argument. This vignette will outline a few more details about these functions and the package that will hopefully make `ramlegacy` easy to use.

## Package Installation
You can install the development version from [github](https://github.com/ropensci/ramlegacy) with the `install_github` function available in `devtools` package:

```{r, eval = FALSE, echo = TRUE}
install.packages("devtools")
library(devtools)
install_github("ropensci/ramlegacy")
```

To ensure that the vignette is installed along with the package make sure to remove `--no-build-vignettes` from the `build_opts` in `install_github`

## download_ramlegacy
You can download any available version of the RAM Legacy Stock Assessment Excel Database ("1.0","2.0", "2.5", "3.0", "4.3", "4.40", "4.41", and "4.44") through `download_ramlegacy()`:

```{r, download_ramlegacy_example1, eval=FALSE, echo = T}
# download version 1.0
download_ramlegacy(version = "1.0")

# download version 2.0
download_ramlegacy(version = "2.0")

# download version 2.5
download_ramlegacy(version = "2.5")

# download version 3.0
download_ramlegacy(version = "3.0")

# download version 4.3
download_ramlegacy(version = "4.3")

# download version 4.44
download_ramlegacy(version "4.44")
```

If version is not specified then `download_ramlegacy` will default to downloading the
latest version of the database:

```{r, download_ramlegacy_example2, echo = T, eval = F}
# downloads latest version (currently 4.44)
download_ramlegacy()
```
If you want to download multiple versions please download them one at a time as passing them all at once will throw an error.

To enable caching behavior, before downloading the specified version `download_ramlegacy` checks whether that version has already been downloaded. If it is then `download_ramlegacy` gives the user the option to either download again and overwrite the existing version or not download and just exit the function. This behavior is modified if `overwrite`, a function argument to `download_ramlegacy`, is set to true. Then even if the version has already been downloaded calling `download_ramlegacy` will download it again overwriting the already present version. After downloading the database download_ramlegacy caches it by converting it to an RDS object which can then be read in through calling `load_ramlegacy()`. 

`download_ramlegacy` also has three other arguments: `ram_path`, `ram_url`, and `quiet`. By default `ram_path` has been set to a location provided by the rappdirs package and can be viewed by calling `ram_dir(vers = "version")`. Although this is **not the recommended approach** download_ramlegacy supports downloading to a user-specified path. `ram_url` is fixed to the zenodo URL of the database. Please do not **pass in any other url** to ram_url. By default, `download_ramlegacy` after it is called keeps the user updated of the download progress through status messages. Setting `quiet` to TRUE enables the user to suppress all status messages.

### File Structure of Downloaded Database
This section illustrates the contents and file structure of the directory to which the database is saved once you call `download_ramlegacy`

If you have downloaded multiple versions through `download_ramlegacy()` the contents of ramlegacy directory (the path to this directory can be viewed by calling `ram_dir()`) on your local computer will look something like this with each version getting its own subdirectory within the directory:


![ramlegacy directory structure](subdir.png)


For older versions (1.0, 2.0, 2.5, 3.0, 4.3) the contents of the version subdirectory look like this:


![older version subdirectory structure](v4.3.png)


For these older versions, the version subdirectory would contain the excel database and serialized database saved under the names `RLSADB v <version no> (assessment data only).xlsx` and `v <version number> .rds` respectively.

For newer versions (4.40, 4.41, 4.44) along with the excel file of the database file other documents like stock summary files are also downloaded:

![newer version subdirectory structure](v4.44.png)


The actual excel file of the database and the serialized rds database will be present as `RLSADB v <version no> (assessment data only).xlsx` and `v <version no> .rds` respectively inside the `DB Files With Assessment Data` folder.

## Description of Tables and Stock Summary Files within the Database
For a description of the dataframes present in the database as well as information on the stock summary files 
present in the database please see this [page.](https://ropensci.github.io/ramlegacy/articles/tables_description.html)

## load_ramlegacy
`load_ramlegacy` can be used to obtain a list of particular dataframes and passed in to the `tables` argument. To get a list of all the dataframes within a specific version of the database set `tables = NULL`. Make sure that the version you want to load was downloaded and is already present locally at the path you specify using `ram_path` otherwise `load_ramlegacy` will throw an error. 

```{r, load_ramlegacy_example1, echo = T, eval = F}
# returns a list containing area and bioparams tables from version 4.44 database
load_ramlegacy(version = "4.44", tables = c("area", "bioparams"))
```

In case in you want to load multiple versions please load them one at a time as passing them all at once will throw an error. Similar, to `download_ramlegacy`, `load_ramlegacy` has a `ram_path` argument which is set to local directory where the specified version of the RAM Legacy Stock Excel Assessment Database was downloaded. You can view that path through `ram_dir(vers = "version")`. Although this is **not the recommended approach** load_ramlegacy supports reading in the database's dataframes from a user-specified path provided that the database is present at the specified path as an rds file.
  
## ram_dir
To see the path of the directory where the specified version of database was downloaded on the user's computer, the user can call `ram_dir`. Note that this is also the location from which `load_ramlegacy` reads in the database. 

For example, to see where version 4.3 of the database was downloaded you can call:

```{r, ram_dir_example1, echo = T, eval = F}
# view the path where version 4.3 of the database was cached
ram_dir(vers = "4.3")
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ramlegacy-package.R
\docType{package}
\name{ramlegacy-package}
\alias{ramlegacy}
\alias{ramlegacy-package}
\title{ramlegacy: download, cache and read RAM Legacy Stock Assessment Database}
\description{
ramlegacy contains functions to download, cache and read in specified tables
from the excel version of the RAM Legacy Stock Assessment Database,
an online compilation of stock assessment results for commercially exploited
marine populations from around the world. More information about the database
can be found at <https://ramlegacy.org/>.
}
\section{Description of the dataframes present in the database}{


\itemize{
\item metadata: Table with summarized metadata (only available in newer
 versions starting from v4.40)
\item stock: This stores the stock database table
\item assessment: This stores the assessment database table
\item taxonomy: This stores the taxonomy database table
\item management: This stores the management database table
\item assessor: This stores the assessor database table
\item assessmetod: This stores the assessmetod database table
\item area: This stores the area database table
\item biometrics: This stores the biometrics database table
\item tsmetrics: This stores the tsmetrics database table
\item timeseries: The time series data is a matrix object with the following
headers/columns: (1) assessid (2) stockid (3) stocklong (4) tsid (5) tsyear
(6) tsvalue
\item bioparams: The parameter data is a matrix object with the following
headers/columns: (1) assessid (2) stockid (3) stocklong (4) bioid (5) biovalue
(6) bioyear (7) bionotes
\item timeseries_values_views: This stores the timeseries values with timeseries
type along the columns and stocks along the rows
\item timeseries_units_views: This stores the timeseries values with timeseries
type along the columns and stocks along the rows
\item timeseries_ids_views: This stores the timeseries IDs with timeseries type
 along the columns and stocks along the rows
\item timeseries_assessments_views: This stores the timeseries assessments with
timeseries type along the columns and stocks along the rows
\item timeseries_notes_views: This stores the timeseries notes with timeseries
type along the columns and stocks along the rows
\item timeseries_sources_views: This stores the timeseries sources with timeseries
type along the columns and stocks along the rows
\item timeseries_years_views: This stores the timeseries years with timeseries
type along the columns and stocks along the rows
\item bioparams_values_views: This stores the reference point values, with
reference point type along the columns and stocks along the rows
\item bioparams_units_views: This stores the reference point units, with
reference point type along the columns and stocks along the rows
\item bioparams_ids_views: This stores the reference point IDs, with reference
point type along the columns and stocks along the rows
\item bioparams_assessments_views: This stores the reference point assessments,
with reference point type along the columns and stocks along the rows
\item bioparams_sources_views: This stores the reference point sources, with
reference point type along the columns and stocks along the rows
\item bioparams_notes_views: This stores the reference point notes, with
reference point type along the columns and stocks along the rows
}
}

\section{Newer versions (v4.40 onwards) also contains tables of individual most-used time series}{

\itemize{
\item tb.data: Total Biomass
\item ssb.data: Spawning Stock Biomass
\item tn.data: Total Abundance
\item r.data: Recruits
\item tc.data: Total Catch
\item tl.data: Total Landings
\item recc.data: Recreational Catch
\item f.data: Fishing Mortality
\item er.data: Exploitation Rate
\item divtb.data: TB/TBmsy
\item divssb.data: SSB/SSBmsy
\item ivf.data: F/Fmsy
\item diver.data: ER/ERmsy
\item divbpref.data: B/Bmsypref
\item divupref.data: U/Umsypref
\item tbbest.data: TBbest
\item tcbest.data: TCbest
\item erbest.data: ERbest
\item divtb.mgt.data: TB/TBmgt
\item divssb.mgt.data: SSB/SSBmgt
\item divf.mgt.data: F/Fmgt
\item diver.mgt.data: ER/ERmgt
\item divbpref.mgt.data: B/Bmgtpref
\item divupref.mgt.data: U/Umgtpref
\item cpair.data: Cpair
\item tac.data: TAC
\item cadv.data: Cadvised
\item survb.data: survB
\item cpue.data: CPUE
\item effort.data: EFFORT
}
}

\references{
Ricard, D., Minto, C., Jensen, O.P. and Baum, J.K. (2012)
Evaluating the knowledge base and status of commercially exploited marine
species with the RAM Legacy Stock Assessment Database.
Fish and Fisheries 13 (4) 380-398. <doi:10.1111/j.1467-2979.2011.00435.x>
}
\seealso{
\url{www.ramlegacy.org}

\url{www.github.com/ropensci/ramlegacy}

\url{www.github.com/ropensci/ramlegacy/issues}
}
\author{
\strong{Maintainer}: Kshitiz Gupta \email{kshtzgupta1@berkeley.edu} [copyright holder]

Authors:
\itemize{
  \item Carl Boettiger \email{cboettig@gmail.com} (http://orcid.org/0000-0002-1642-628X) [copyright holder]
}

Other contributors:
\itemize{
  \item Sam Albers \email{sam.albers@gmail.com} (https://orcid.org/0000-0002-9270-7884) [reviewer]
  \item Jamie Afflerbach \email{afflerbach@nceas.ucsb.edu} (https://orcid.org/0000-0002-5215-9342) [reviewer]
  \item RAM Legacy Stock Assessment Database [data contributor]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/load_ramlegacy.R
\name{load_ramlegacy}
\alias{load_ramlegacy}
\title{Read-in dataframes from downloaded RAM Legacy Database}
\usage{
load_ramlegacy(version = NULL, tables = NULL, ram_path = NULL)
}
\arguments{
\item{version}{A character vector of length 1 specifying the version number
of the database. As of writing this package, the available versions are
"1.0", "2.0", "2.5", "3.0", "4.3", "4.40", "4.41" and "4.44". If version
argument is not specified then it defaults to newest version (v4.44).}

\item{tables}{A character vector specifying the names of particular
dataframes to load from a specified version. If not specified then
a list containing all the dataframes within the requested version
is returned. For a description of the different tables present
in the database please see below.}

\item{ram_path}{path to the local directory where the specified version of
the RAM Legacy Stock Excel Assessment Database was downloaded. By default
this path is set to within the rappdirs directory and can be viewed using
calling the function \code{\link{ram_dir}} and specifying the version number
inside the function call. Although this is not the \strong{recommended}
approach \code{load_ramlegacy} supports reading in the database's
dataframes from a user-specified path provided that the database is present
at the specified path as an rds object.}
}
\description{
Returns a list of specific dataframes or a list of all the
 dataframes present in the requested version of the database.
}
\section{Description of the dataframes present in the database}{


\itemize{
\item metadata: Table with summarized metadata (only available in newer
 versions starting from v4.40)
\item stock: This stores the stock database table
\item assessment: This stores the assessment database table
\item taxonomy: This stores the taxonomy database table
\item management: This stores the management database table
\item assessor: This stores the assessor database table
\item assessmetod: This stores the assessmetod database table
\item area: This stores the area database table
\item biometrics: This stores the biometrics database table
\item tsmetrics: This stores the tsmetrics database table
\item timeseries: The time series data is a matrix object with the following
headers/columns: (1) assessid (2) stockid (3) stocklong (4) tsid (5) tsyear
(6) tsvalue
\item bioparams: The parameter data is a matrix object with the following
headers/columns: (1) assessid (2) stockid (3) stocklong (4) bioid (5) biovalue
(6) bioyear (7) bionotes
\item timeseries_values_views: This stores the timeseries values with timeseries
type along the columns and stocks along the rows
\item timeseries_units_views: This stores the timeseries values with timeseries
type along the columns and stocks along the rows
\item timeseries_ids_views: This stores the timeseries IDs with timeseries type
 along the columns and stocks along the rows
\item timeseries_assessments_views: This stores the timeseries assessments with
timeseries type along the columns and stocks along the rows
\item timeseries_notes_views: This stores the timeseries notes with timeseries
type along the columns and stocks along the rows
\item timeseries_sources_views: This stores the timeseries sources with timeseries
type along the columns and stocks along the rows
\item timeseries_years_views: This stores the timeseries years with timeseries
type along the columns and stocks along the rows
\item bioparams_values_views: This stores the reference point values, with
reference point type along the columns and stocks along the rows
\item bioparams_units_views: This stores the reference point units, with
reference point type along the columns and stocks along the rows
\item bioparams_ids_views: This stores the reference point IDs, with reference
point type along the columns and stocks along the rows
\item bioparams_assessments_views: This stores the reference point assessments,
with reference point type along the columns and stocks along the rows
\item bioparams_sources_views: This stores the reference point sources, with
reference point type along the columns and stocks along the rows
\item bioparams_notes_views: This stores the reference point notes, with
reference point type along the columns and stocks along the rows
}
}

\section{Newer versions (v4.40 onwards) also contains tables of individual most-used time series}{

\itemize{
\item tb.data: Total Biomass
\item ssb.data: Spawning Stock Biomass
\item tn.data: Total Abundance
\item r.data: Recruits
\item tc.data: Total Catch
\item tl.data: Total Landings
\item recc.data: Recreational Catch
\item f.data: Fishing Mortality
\item er.data: Exploitation Rate
\item divtb.data: TB/TBmsy
\item divssb.data: SSB/SSBmsy
\item ivf.data: F/Fmsy
\item diver.data: ER/ERmsy
\item divbpref.data: B/Bmsypref
\item divupref.data: U/Umsypref
\item tbbest.data: TBbest
\item tcbest.data: TCbest
\item erbest.data: ERbest
\item divtb.mgt.data: TB/TBmgt
\item divssb.mgt.data: SSB/SSBmgt
\item divf.mgt.data: F/Fmgt
\item diver.mgt.data: ER/ERmgt
\item divbpref.mgt.data: B/Bmgtpref
\item divupref.mgt.data: U/Umgtpref
\item cpair.data: Cpair
\item tac.data: TAC
\item cadv.data: Cadvised
\item survb.data: survB
\item cpue.data: CPUE
\item effort.data: EFFORT
}
}

\examples{
\donttest{
\dontshow{Sys.setenv(RAM_HOME = tempfile())}
# first download version 4.44 of the database
download_ramlegacy(version = "4.44")

# get a list containing area and bioparams tables
# from version 4.44 database
load_ramlegacy(version = "4.44", tables = c("area", "bioparams"))
}
}
\seealso{
Other ramlegacy functions: \code{\link{download_ramlegacy}},
  \code{\link{ram_dir}}
}
\concept{ramlegacy functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{ram_dir}
\alias{ram_dir}
\title{Output OS-independent path to the rappdirs directory on user's computer where
the RAM Legacy database is downloaded by default}
\usage{
ram_dir(vers = NULL)
}
\arguments{
\item{vers}{character, version number of the database. As of writing this
package, the available versions are "1.0", "2.0", "2.5", "3.0", "4.3","4.40",
"4.41", and "4.44". If version is not specified the \code{ram_dir()}
returns the path to the rappdirs top-level directory which stores
all the version subdirectories.}
}
\description{
Provides the download location for \code{\link{download_ramlegacy}}
in an OS independent manner. This is also the location from where
\code{\link{load_ramlegacy}} loads the database from.
}
\examples{
# return the path to the rappdirs directory where
# all version subdirectories are stored
ram_dir()

# Returns the path of the subdirectory where v4.3
# of the database is downloaded to and read from.
ram_dir(vers = "4.3")
}
\seealso{
Other ramlegacy functions: \code{\link{download_ramlegacy}},
  \code{\link{load_ramlegacy}}
}
\concept{ramlegacy functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/download_ramlegacy.R
\name{download_ramlegacy}
\alias{download_ramlegacy}
\title{Download RAM Legacy Excel Database}
\usage{
download_ramlegacy(version = NULL, ram_path = NULL,
  ram_url = "https://doi.org/10.5281/zenodo.2542918",
  overwrite = FALSE, quiet = FALSE)
}
\arguments{
\item{version}{A character vector of length 1 specifying the version number
of the database that should be downloaded. As of writing this package, the
available versions are "1.0","2.0", "2.5", "3.0", "4.3", "4.40", "4.41",
and "4.44". If the version argument is not specified then it defaults
to latest version (currently latest version is "4.44").}

\item{ram_path}{A string specifying the path of the local directory where
database will be downloaded. By default this path is set to the location
provided by \pkg{rappdirs} package and can be viewed by calling
\code{\link{ram_dir}}. Although this is not the \strong{recommended}
approach \code{download_ramlegacy} supports downloading to a user-specified path.}

\item{ram_url}{A string. By default it is set to the Zenodo url of the database.
Please \strong{do not pass} in any other url to \code{ram_url}.}

\item{overwrite}{This argument is only relevant if you are trying to
re-download a version that is already present locally in the rappdirs
directory. If overwrite = TRUE then user will not encounter the interactive prompt
that confirms whether to overwrite the version present locally.}

\item{quiet}{If TRUE, suppress status messages}
}
\description{
Downloads a specified version of RAM Legacy Stock Assessment
Excel Database and as an RDS object to a local directory specified by
\code{\link{ram_dir}}. The function will check if the version requested
already exists on the user's computer and if it does then it prompts the user
to download it again unless `overwrite = TRUE` in which case the function will
download the version without displaying the prompt. The function also supports
downloading all the older versions (1.0, 2.0, 2.5, 3.0, 4.3) from
[a github repo](www.github.com/kshtzgupta1/ramlegacy-assets)
}
\examples{
\donttest{
\dontshow{
Sys.setenv(RAM_HOME = tempfile())
}
# download version 4.44
download_ramlegacy(version = "4.44")

# download version 1.0
download_ramlegacy(version = "1.0")
\dontshow{
Sys.setenv(RAM_HOME = tempfile())
}
# If version not specified then default
# to latest version (currently 4.44)
download_ramlegacy()
}
}
\seealso{
Other ramlegacy functions: \code{\link{load_ramlegacy}},
  \code{\link{ram_dir}}
}
\concept{ramlegacy functions}
