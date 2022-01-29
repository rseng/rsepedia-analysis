
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
