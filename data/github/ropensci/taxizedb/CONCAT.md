taxizedb
========



[![cran checks](https://cranchecks.info/badges/worst/taxizedb)](https://cranchecks.info/pkgs/taxizedb)
[![R-check](https://github.com/ropensci/taxizedb/workflows/R-check/badge.svg)](https://github.com/ropensci/taxizedb/actions?query=workflow%3AR-check)
[![CircleCI](https://circleci.com/gh/ropensci/taxizedb.svg?style=svg)](https://circleci.com/gh/ropensci/taxizedb)
[![codecov](https://codecov.io/gh/ropensci/taxizedb/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/taxizedb)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/taxizedb)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/taxizedb)](https://cran.r-project.org/package=taxizedb)
[![DOI](https://zenodo.org/badge/53961466.svg)](https://zenodo.org/badge/latestdoi/53961466)

`taxizedb` - Tools for Working with Taxonomic Databases on your machine

Docs: <https://docs.ropensci.org/taxizedb/>

[taxize](https://github.com/ropensci/taxize) is a heavily used taxonomic toolbelt
package in R - However, it makes web requests for nearly all methods. That is fine
for most cases, but when the user has many, many names it is much more efficient
to do requests to a local SQL database.

## Data sources

Not all taxonomic databases are publicly available, or possible to mash into a SQLized
version. Taxonomic DB's supported:

- NCBI: text files are provided by NCBI, which we stitch into a sqlite db
- ITIS: they provide a sqlite dump, which we use here
- The PlantList: created from stitching together csv files. this
 source is no longer updated as far as we can tell. they say they've
 moved focus to the World Flora Online
- Catalogue of Life: created from Darwin Core Archive dump.
- GBIF: created from Darwin Core Archive dump. right now we only have
 the taxonomy table (called gbif), but will add the other tables in the
 darwin core archive later
- Wikidata: aggregated taxonomy of Open Tree of Life, GLoBI and Wikidata. 
 On Zenodo, created by Joritt Poelen of GLOBI.
- World Flora Online: http://www.worldfloraonline.org/

Update schedule for databases:

- NCBI: since `db_download_ncbi` creates the database when the function
is called, it's updated whenever you run the function
- ITIS: since ITIS provides the sqlite database as a download, you can
delete the old file and run `db_download_itis` to get a new dump;
they I think update the dumps every month or so
- The PlantList: no longer updated, so you shouldn't need to download
this after the first download. hosted on Amazon S3
- Catalogue of Life: a GitHub Actions job runs once a day at 00:00 UTC,
building the lastest COL data into a SQLite database thats hosted on
Amazon S3
- GBIF: a GitHub Actions job runs once a day at 00:00 UTC,
building the lastest GBIF data into a SQLite database thats hosted on
Amazon S3
- Wikidata: last updated April 6, 2018. Scripts are available to 
update the data if you prefer to do it yourself.
- World Flora Online: since `db_download_wfo` creates the database when
the function is called, it's updated whenever you run the function

 Links:

- NCBI: ftp://ftp.ncbi.nih.gov/pub/taxonomy/
- ITIS: https://www.itis.gov/downloads/index.html
- The PlantList - http://www.theplantlist.org/
- Catalogue of Life:
  - latest monthly edition via http://www.catalogueoflife.org/DCA_Export/archive.php
- GBIF: http://rs.gbif.org/datasets/backbone/
- Wikidata: https://zenodo.org/record/1213477
- World Flora Online: http://www.worldfloraonline.org/

Get in touch [in the issues](https://github.com/ropensci/taxizedb/issues) with
any ideas on new data sources.

All databases are SQLite.

## Package API

This package for each data sources performs the following tasks:

* Downloaded taxonomic databases `db_download_*`
* Create `dplyr` SQL backend via `dbplyr::src_dbi` - `src_*` 
* Query and get data back into a data.frame - `sql_collect`
* Manage cached database files - `tdb_cache`
* Retrieve immediate descendents of a taxon - `children`
* Retrieve the taxonomic hierarchies from local database - `classification`
* Retrieve all taxa descending from a vector of taxa - `downstream`
* Convert species names to taxon IDs - `name2taxid`
* Convert taxon IDs to species names - `taxid2name`
* Convert taxon IDs to ranks - `taxid2rank`

You can use the `src` connections with `dplyr`, etc. to do operations downstream. Or use the database connection to do raw SQL queries.

## install

cran version


```r
install.packages("taxizedb")
```

dev version


```r
remotes::install_github("ropensci/taxizedb")
```

## Contributors

* [Scott Chamberlain](https://github.com/sckott)
* [Zebulun Arendsee](https://github.com/arendsee)
* [Tamora James](https://github.com/tdjames1)


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/taxizedb/issues).
* License: MIT
* Get citation information for `taxizedb` in R doing `citation(package = 'taxizedb')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.

[![ropensci](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
taxizedb 0.3.0
==============

## NEW FEATURES

* `db_download()` gains new parameter `overwrite` (logical): used to state that you want to overwrite an existing database on disk. before this you would have to manually delete an older database file (#34)
* new function added `taxa_at()` for getting taxa at specific scientific ranks. For example, your known taxon is the family Lachnospiraceae with NCBI identifier of 186803. You want information on the phylum which the Lachnospiraceae family is in. This function can do that for you.  (#51)

## BUG FIXES

* fixed problem in internal function `txdb_rr()`: in older verions of R (e.g., 3.6) we were creating a data.frame in this function without settings `stringsAsFactors=FALSE`, resulting in different behavior in R v3 vs. R v4 given the change in `stringsAsFactors` behavior in R v4 onward (#54)


taxizedb 0.2.2
==============

## BUG FIXES

* fix failing tests (#50)


taxizedb 0.2.0
==============

## NEW FEATURES

* gains 3 new data sources: NCBI taxonomy, World Flora Online, Wikidata (#18) (#49) (#37)
* gains ports of `taxize` functions to `taxizedb` (NCBI & ITIS supported): `children`, `classification`, `downstream`. beware when both `taxize` and `taxizedb` loaded in the same R session to namespace calls to these three functions (#19) (#25) (#44) (#48)
* gains mapping functions: `name2taxid` (scientific or common name to taxonomy ID); `taxid2name` (taxonomy ID to scientific name); `taxid2rank` (taxonomy ID to rank) (#41) (#42)
* intro vignette added (#17)
* GBIF and COL data sources are not updated daily in the repos https://github.com/ropenscilabs/gbif-backbone-sql and https://github.com/ropenscilabs/col-sql/ via GithHub Actions. See those repos for details (#26)
* update package level manual file (`?taxizedb-package`) with details on each data source, their update schedules, and examples
* all data sources now use SQLite as the database storage engine. passwords/ports/usernames/etc are no longer needed! note that some `db_download*` functions download already created SQLite databases, whereas for other data sources the database is built locally on your machine from other data formats downloaded (see also #36, #46)
* a copy of the taxonomic ranks information from taxize package was ported over for internal use to be able to make `downstream()` work for most data sources

## MINOR IMPROVEMENTS

* remove check for whether SQLite is installed (#5 #29)
* all `src_*` functions now only have two paramters: `path` and `...`. where path by default figures out the path for you using the function `db_path()`, and `...` allows the user to pass on parameters to `DBI::dbConnect`

## DEFUNCT

* `db_load()` is now defunct. Now just use `db_download*` then `src*` for your data source (see also #43)


taxizedb 0.1.4
==============

## BUG FIXES

* Fixes to SQL database connection functions for changes in `dplyr`, 
which now requires `dbplyr` package - also `DBI` now imported (#16)


taxizedb 0.1.0
==============

## NEW FEATURES

* Released to CRAN
## Test environments

* local OS X install, R 4.0.3 patched
* ubuntu 16.04 (on GitHub Actions), R 4.0.3
* win-builder (devel and release)

## R CMD check results

No notes

## Reverse dependencies

There are no reverse dependencies.

---

This version fixes a bug, adds a new function, and one function gains a new parameter.

Thanks!
Scott Chamberlain
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
# CONTRIBUTING #

### Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/taxizedb/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/taxizedb.git`
* Make sure to track progress upstream (i.e., on our version of `taxizedb` at `ropensci/taxizedb`) by doing `git remote add upstream https://github.com/ropensci/taxizedb.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/taxizedb`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Prefer to Email? Get in touch: [myrmecocystus@gmail.com](mailto:myrmecocystus@gmail.com)
<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/) If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
