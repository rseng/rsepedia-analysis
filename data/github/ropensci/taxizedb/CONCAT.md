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
taxizedb
========

```{r echo=FALSE}
knitr::opts_chunk$set(
  warning = FALSE,
  message = FALSE,
  collapse = TRUE,
  comment = "#>"
)
```

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

```{r eval=FALSE}
install.packages("taxizedb")
```

dev version

```{r eval=FALSE}
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
---
title: "taxizedb"
author: Scott Chamberlain
date: "2021-01-14"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{taxizedb}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---




```r
library("taxizedb")
library("dplyr")
```

## Download DBs

ITIS


```r
db_download_itis()
```

The Plant List (TPL)


```r
db_download_tpl()
```

Catalogue of Life (COL)


```r
db_download_col()
```

## connect to the DBs

By default `src_*` functions use a path to the cached database file.
You can alternatively pass in your own path if you've put it 
somewhere else.

ITIS


```r
src_itis <- src_itis()
```

TPL


```r
src_tpl <- src_tpl()
```

COL


```r
src_col <- src_col()
```

## query with SQL syntax


```r
sql_collect(src_itis, "select * from hierarchy limit 5")
#> # A tibble: 5 x 5
#>                     hierarchy_string    tsn parent_tsn level childrencount
#> *                              <chr>  <int>      <int> <int>         <int>
#> 1                             202422 202422          0     0        154282
#> 2                      202422-846491 846491     202422     1          2666
#> 3               202422-846491-660046 660046     846491     2          2654
#> 4        202422-846491-660046-846497 846497     660046     3             7
#> 5 202422-846491-660046-846497-846508 846508     846497     4             6
```


```r
# or pipe the src to sql_collect
src_itis %>% sql_collect("select * from hierarchy limit 5")
```

## use dplyr verbs

get a `tbl`


```r
hiers <- src_itis %>% tbl("hierarchy")
#> # Source:   table<hierarchy> [?? x 5]
#> # Database: postgres 9.6.0 [sacmac@localhost:5432/ITIS]
#>                                              hierarchy_string    tsn parent_tsn level childrencount
#>                                                         <chr>  <int>      <int> <int>         <int>
#>  1                                                     202422 202422          0     0        154282
#>  2                                              202422-846491 846491     202422     1          2666
#>  3                                       202422-846491-660046 660046     846491     2          2654
#>  4                                202422-846491-660046-846497 846497     660046     3             7
#>  5                         202422-846491-660046-846497-846508 846508     846497     4             6
#>  6                  202422-846491-660046-846497-846508-846553 846553     846508     5             5
#>  7           202422-846491-660046-846497-846508-846553-954935 954935     846553     6             3
#>  8      202422-846491-660046-846497-846508-846553-954935-5549   5549     954935     7             2
#>  9 202422-846491-660046-846497-846508-846553-954935-5549-5550   5550       5549     8             0
#> 10           202422-846491-660046-846497-846508-846553-954936 954936     846553     6             0
#> # ... with more rows
```

select certain fields


```r
hiers %>% select(TSN, level)
#> # Source:   lazy query [?? x 2]
#> # Database: postgres 9.6.0 [sacmac@localhost:5432/ITIS]
#>       tsn level
#>     <int> <int>
#>  1 202422     0
#>  2 846491     1
#>  3 660046     2
#>  4 846497     3
#>  5 846508     4
#>  6 846553     5
#>  7 954935     6
#>  8   5549     7
#>  9   5550     8
#> 10 954936     6
#> # ... with more rows
```

## Local versions of `taxize` functions

A few of the key functions from `taxize` have been ported to `taxizedb`.
Support is currently limited to the NCBI taxonomy database.

`children` accesses the nodes immediately descending from a given taxon


```r
children(3701, db='ncbi')
#> $`3701`
#>    childtaxa_id                                                     childtaxa_name childtaxa_rank
#> 1       1837063                         Arabidopsis thaliana x Arabidopsis halleri        species
#> 2       1547872                                              Arabidopsis umezawana        species
#> 3       1328956 (Arabidopsis thaliana x Arabidopsis arenosa) x Arabidopsis suecica        species
#> 4       1240361                         Arabidopsis thaliana x Arabidopsis arenosa        species
#> 5        869750                          Arabidopsis thaliana x Arabidopsis lyrata        species
#> 6        412662                                            Arabidopsis pedemontana        species
#> 7        378006                         Arabidopsis arenosa x Arabidopsis thaliana        species
#> 8        347883                                              Arabidopsis arenicola        species
#> 9        302551                                              Arabidopsis petrogena        species
#> 10        97980                                               Arabidopsis croatica        species
#> 11        97979                                            Arabidopsis cebennensis        species
#> 12        81970                                                Arabidopsis halleri        species
#> 13        59690                                             Arabidopsis kamchatica        species
#> 14        59689                                                 Arabidopsis lyrata        species
#> 15        45251                                               Arabidopsis neglecta        species
#> 16        45249                                                Arabidopsis suecica        species
#> 17        38785                                                Arabidopsis arenosa        species
#> 18         3702                                               Arabidopsis thaliana        species
#> 
#> attr(,"class")
#> [1] "children"
#> attr(,"db")
#> [1] "ncbi"
```

`classification` finds the lineage of a taxon


```r
classification(3702, db='ncbi')
#> $`3702`
#>                    name         rank      id
#> 1    cellular organisms      no rank  131567
#> 2             Eukaryota superkingdom    2759
#> 3         Viridiplantae      kingdom   33090
#> 4          Streptophyta       phylum   35493
#> 5        Streptophytina    subphylum  131221
#> 6           Embryophyta      no rank    3193
#> 7          Tracheophyta      no rank   58023
#> 8         Euphyllophyta      no rank   78536
#> 9         Spermatophyta      no rank   58024
#> 10        Magnoliophyta      no rank    3398
#> 11      Mesangiospermae      no rank 1437183
#> 12       eudicotyledons      no rank   71240
#> 13           Gunneridae      no rank   91827
#> 14         Pentapetalae      no rank 1437201
#> 15               rosids     subclass   71275
#> 16              malvids      no rank   91836
#> 17          Brassicales        order    3699
#> 18         Brassicaceae       family    3700
#> 19           Camelineae        tribe  980083
#> 20          Arabidopsis        genus    3701
#> 21 Arabidopsis thaliana      species    3702
#> 
#> attr(,"class")
#> [1] "classification"
#> attr(,"db")
#> [1] "ncbi"
```

`downstream` finds all taxa descending from a taxon


```r
downstream(3700, db='ncbi')
#> $`3700`
#>      childtaxa_id                                                     childtaxa_name     rank
#> 1         2071891                                                     Draba taylorii  species
#> 2         2071524                                                  Rorippa tenerrima  species
#> 3         2071523                                                Rorippa crystallina  species
#> 4         2071509                                                   Physaria calderi  species
#> 5         2071468                                                 Erysimum arenicola  species
#> 6         2071452                                                   Draba yukonensis  species
#> 7         2071451                                                   Draba thompsonii  species
#> ...
#> 326       1492251                                                 Erysimum lilacinum  species
#> 327       1492250                                              Erysimum leucanthemum  species
#> 328       1492249                                               Erysimum leptostylum  species
#> 329       1492248                                              Erysimum leptophyllum  species
#> 330       1492247                                               Erysimum leptocarpum  species
#> 331       1492246                                                Erysimum ledebourii  species
#> 332       1492245                                                Erysimum laxiflorum  species
#> 333       1492244                                                  Erysimum kurdicum  species
#>  [ reached getOption("max.print") -- omitted 2880 rows ]
#>
#> attr(,"class")
#> [1] "downstream"
#> attr(,"db")
#> [1] "ncbi"
```

All of these functions run very fast. It only takes a few seconds to find all
bacterial taxa and count them: 


```r
downstream(2, db='ncbi')[[1]] %>%
    dplyr::group_by(rank) %>%
    dplyr::count()
#> #> [1] 138695
#> # A tibble: 18 x 2
#> # Groups:   rank [18]
#> rank                 n
#> <chr>            <int>
#>  1 class               83
#>  2 family             483
#>  3 forma                4
#>  4 genus             3497
#>  5 no rank          37140
#>  6 order              198
#>  7 phylum             134
#>  8 species          97031
#>  9 species group       68
#> 10 species subgroup    10
#> 11 subclass             3
#> 12 subfamily            1
#> 13 subgenus             1
#> 14 suborder             8
#> 15 subphylum            1
#> 16 subspecies          10
#> 17 tribe                2
#> 18 varietas            21
```

## Mapping functions

Several mapping functions are available for the NCBI taxonomy database:


```r
# Map scientific or common names to taxonomy IDs
name2taxid("pig")
#> [1] "9823"

# Map taxonomy IDs to scientific names
taxid2name(9823)
#> [1] "Sus scrofa"

# Map taxonomy IDs to rank
taxid2rank(2)
#> [1] "superkingdom"
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/src.R
\name{src_taxizedb}
\alias{src_taxizedb}
\alias{src_itis}
\alias{src_tpl}
\alias{src_col}
\alias{src_gbif}
\alias{src_ncbi}
\alias{src_wikidata}
\alias{src_wfo}
\title{src - dplyr src objects}
\usage{
src_itis(path = db_path("itis"), ...)

src_tpl(path = db_path("tpl"), ...)

src_col(path = db_path("col"), ...)

src_gbif(path = db_path("gbif"), ...)

src_ncbi(path = db_path("ncbi"), ...)

src_wikidata(path = db_path("wikidata"), ...)

src_wfo(path = db_path("wfo"), ...)
}
\arguments{
\item{path}{(character) path to SQLite database. by default
we use the function \code{\link[=db_path]{db_path()}} to get the path}

\item{...}{Further args passed on to \code{\link[DBI:dbConnect]{DBI::dbConnect()}}}
}
\value{
an src object
}
\description{
src - dplyr src objects
}
\examples{
\dontrun{
# src_itis()
# src_tpl()
# src_col()
# src_gbif()
# src_ncbi()
# src_wikidata()
# src_wfo()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ap_taxid2rank.R
\name{taxid2rank}
\alias{taxid2rank}
\title{Convert taxon IDs to scientific ranks}
\usage{
taxid2rank(x, db = "ncbi", verbose = TRUE, warn = TRUE, ...)
}
\arguments{
\item{x}{(character) Vector of taxon keys (name or id) for the given
database}

\item{db}{(character) The database to search, one of ncbi, itis, gbif,
col, or wfo}

\item{verbose}{(logical) Print verbose messages}

\item{warn}{(logical) If \code{TRUE}, raise a warning if any taxon IDs can not
be found}

\item{...}{Additional arguments passed to database specific classification
functions}
}
\value{
character vector of ranks in the same order as the inputs
}
\description{
Convert taxon IDs to scientific ranks
}
\examples{
\dontrun{
taxid2rank(c(3701, 9606))
taxid2rank(c(154395, 154357, 23041, 154396), db="itis")
taxid2rank(c('wfo-4000032377', 'wfo-0000541830'), db="wfo")
taxid2rank("wfo-7000000057", db="wfo")
taxid2rank(2877951, db="gbif")
taxid2rank(c(2877951, 5386), db="gbif")
taxid2rank(c(3960765, 3953606, 3953010), db="col")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_path.R
\name{db_path}
\alias{db_path}
\title{database path}
\usage{
db_path(db)
}
\arguments{
\item{db}{(character) db name. one of: itis, tpl, col, gbif,
ncbi, wikidata, wfo. required}
}
\description{
database path
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_load.R
\name{db_load-defunct}
\alias{db_load-defunct}
\alias{db_load_itis}
\alias{db_load_tpl}
\alias{db_load_col}
\alias{db_load_gbif}
\alias{db_load_ncbi}
\alias{db_load_wikidata}
\title{Load taxonomic databases - NO LONGER NEEDED}
\usage{
db_load_itis(...)

db_load_tpl(...)

db_load_col(...)

db_load_gbif(...)

db_load_ncbi(...)

db_load_wikidata(...)
}
\arguments{
\item{...}{ignored}
}
\description{
Use \link{db_download} then \link{src_taxizedb}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caching.R
\name{tdb_cache}
\alias{tdb_cache}
\title{Caching}
\description{
Manage cached taxizedb files with \pkg{hoardr}
}
\details{
\code{cache_delete} only accepts 1 file name, while
\code{cache_delete_all} doesn't accept any names, but deletes all files.
For deleting many specific files, use \code{cache_delete} in a \code{\link[=lapply]{lapply()}}
type call
}
\section{Useful user functions}{

\itemize{
\item \code{tdb_cache$cache_path_get()} get cache path
\item \code{tdb_cache$cache_path_set()} set cache path
\item \code{tdb_cache$list()} returns a character vector of full
path file names
\item \code{tdb_cache$files()} returns file objects with metadata
\item \code{tdb_cache$details()} returns files with details
\item \code{tdb_cache$delete()} delete specific files
\item \code{tdb_cache$delete_all()} delete all files, returns nothing
}
}

\examples{
\dontrun{
tdb_cache

# list files in cache
tdb_cache$list()

# delete certain database files
# tdb_cache$delete("file path")
# tdb_cache$list()

# delete all files in cache
# tdb_cache$delete_all()
# tdb_cache$list()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ap_children.R
\name{children}
\alias{children}
\title{Retrieve immediate descendents of a taxon}
\usage{
children(x, db = "ncbi", verbose = TRUE, ...)
}
\arguments{
\item{x}{(character) Vector of taxon keys for the given database}

\item{db}{(character) The database to search, one of ncbi, itis, gbif,
col, or wfo}

\item{verbose}{(logical) Print verbose messages}

\item{...}{Additional arguments passed to database specific function.}
}
\value{
list of tibbles with the columns: id, name, rank. This is exactly
equivalent to the output of \code{taxize::children()}
}
\description{
Retrieve immediate descendents of a taxon
}
\examples{
\dontrun{
children(c(3700, 2))
children(c(154395, 154357), db="itis")
children("wfo-4000032377", db="wfo")
children(2877951, db="gbif")
children(3960765, db="col") # Abies
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_download.R
\name{db_download}
\alias{db_download}
\alias{db_download_ncbi}
\alias{db_download_itis}
\alias{db_download_tpl}
\alias{db_download_wfo}
\alias{db_download_col}
\alias{db_download_gbif}
\alias{db_download_wikidata}
\title{Download taxonomic databases}
\usage{
db_download_ncbi(verbose = TRUE, overwrite = FALSE)

db_download_itis(verbose = TRUE, overwrite = FALSE)

db_download_tpl(verbose = TRUE, overwrite = FALSE)

db_download_wfo(verbose = TRUE, overwrite = FALSE)

db_download_col(verbose = TRUE, overwrite = FALSE)

db_download_gbif(verbose = TRUE, overwrite = FALSE)

db_download_wikidata(verbose = TRUE, overwrite = FALSE)
}
\arguments{
\item{verbose}{(logical) Print messages. Default: \code{TRUE}}

\item{overwrite}{(logical) If \code{TRUE} force an update by overwriting
previously downloaded data. Default: \code{FALSE}}
}
\value{
(character) path to the downloaded SQL database
}
\description{
Download taxonomic databases
}
\details{
Downloads sql database, cleans up unneeded files, returns path
to sql file
}
\examples{
\dontrun{
# ITIS
# db_download_itis()
# src_itis()

# Plantlist
# db_download_tpl()
# db_download_tpl(overwrite=TRUE) # overwrite - download again
# src_tpl()

# COL
# db_download_col()
# src_col()

# GBIF
# db_download_gbif()
# src_gbif()

# NCBI
# db_download_ncbi()
# src_ncbi()

# Wikidata
# db_download_wikidata()
# db_download_wikidata(overwrite=TRUE) # overwrite - download again
# src_wikidata()

# World Flora Online
# db_download_wfo()
# src_wfo()
}
}
\seealso{
\link{tdb_cache}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ap_downstream.R
\name{downstream}
\alias{downstream}
\title{Retrieve all taxa descending from a vector of taxa}
\usage{
downstream(x, db = "ncbi", verbose = TRUE, ...)
}
\arguments{
\item{x}{(character) Vector of taxon keys for the given database}

\item{db}{(character) The database to search, one of ncbi, itis,
gbif, col, or wfo}

\item{verbose}{(logical) Print verbose messages}

\item{...}{Additional arguments passed to database specific downstream
functions}
}
\value{
list of data.frames with the columns: childtaxa_id, childtaxa_name,
and rank. This is exactly equivalent to the output of \code{taxize::downstream()}
}
\description{
This function is nearly equivalent to the \code{taxize::downstream()} function
}
\examples{
\dontrun{
# get descendents from all ranks
# downstream(c(3700, 9605)) # takes a while

# limit results to species
downstream(c(3700, 9605), downto='species')

# allow ambiguous nodes but no ambiguous species
downstream(
  c(3700, 9605),
  downto='species',
  ambiguous_nodes=FALSE,
  ambiguous_species=TRUE
)

# ITIS
id <- name2taxid('Aves', db = "itis")
downstream(id, db = "itis", downto = "family")
downstream(id, db = "itis", downto = "genus")
id <- name2taxid('Bombus', db = "itis")
downstream(id, db = "itis", downto = "species")

# COL
id <- name2taxid('Chordata', db = "col")
downstream(id, db = "col", downto = "family")

# GBIF
id <- name2taxid('Pinaceae', db = "gbif")
downstream(id, db = "gbif", downto = "genus")

# World Flora Online
id <- name2taxid('Pinaceae', db = "wfo")
downstream(id, db = "wfo", downto = "species")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pipe.R
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
% Please edit documentation in R/taxa_at.R
\name{taxa_at}
\alias{taxa_at}
\title{Get taxa at specific scientific ranks}
\usage{
taxa_at(
  x,
  rank,
  db = "ncbi",
  missing = "lower",
  verbose = TRUE,
  warn = TRUE,
  ...
)
}
\arguments{
\item{x}{(character) Vector of taxon keys (ids) for the given
database. required}

\item{rank}{(character) A target rank for which to fetch data. required}

\item{db}{(character) The database to search, one of ncbi, itis, gbif,
col, or wfo}

\item{missing}{(character) if no data found at the given rank and input key,
should we get the next closest lower than that given in \code{rank}, or higher.
one of lower (default), higher.}

\item{verbose}{(logical) Print verbose messages}

\item{warn}{(logical) If \code{TRUE}, raise a warning if any taxon IDs can not
be found}

\item{...}{Additional arguments passed to database specific classification
functions}
}
\value{
list of data.frame's for each input taxon key, where each data.frame
has fields: name, rank, id. When no results found, an empty data.frame
}
\description{
Get taxa at specific scientific ranks
}
\examples{
\dontrun{
taxa_at(186803, rank = "order", db="ncbi", missing = "lower")
taxa_at(c(186803, 541000, 216572, 186804, 31979,  186806),
 rank = "family", missing = "lower")
taxa_at(c(154395, 154357, 23041, 154396), rank = "family", db="itis")
taxa_at(c('wfo-4000032377', 'wfo-0000541830'), rank = "family", db="wfo")
taxa_at("wfo-7000000057", rank = "order", db="wfo")
taxa_at(2877951, rank = "phylum", db="gbif")
taxa_at(c(2877951, 5386), rank = "family", db="gbif")
taxa_at(c(3960765, 3953606, 3953010), rank = "family", db="col")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ap_classification.R
\name{classification}
\alias{classification}
\title{Retrieve the taxonomic hierarchies from local database}
\usage{
classification(x, db = "ncbi", verbose = TRUE, ...)
}
\arguments{
\item{x}{character) Vector of taxon keys for the given database}

\item{db}{character) The database to search, one of ncbi, itis,
gbif, col, or wfo}

\item{verbose}{(logical) Print verbose messages}

\item{...}{Additional arguments passed to database specific classification
functions.}
}
\value{
list of data.frames with the columns: name, rank, and id. This is
exactly equivalent to the output of \code{taxize::classification()}
}
\description{
This function is equivalent to the \code{taxize::classification()} function,
except that it uses a local database (so is much faster). The output is
identical to \code{taxize::classification()}
}
\examples{
\dontrun{
classification(c(3702, 9606))
classification(c(154395, 154357), db="itis")
classification(c("wfo-0000291463", "wfo-7000000057"), db="wfo")
classification(2878586, db="gbif")
classification(c(2878586, 2704179), db="gbif")
classification(3960765, db="col") # Abies
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sql.R
\name{sql_collect}
\alias{sql_collect}
\title{Query and get data back into a data.frame}
\usage{
sql_collect(src, query, ...)
}
\arguments{
\item{src}{(src) An \code{src} object, result of calling \code{\link[=src_itis]{src_itis()}},
\code{\link[=src_col]{src_col()}}, or \code{\link[=src_tpl]{src_tpl()}}}

\item{query}{(character) A SQL query}

\item{...}{further args passed on to \code{\link[dplyr:tbl]{dplyr::tbl()}}}
}
\description{
Query and get data back into a data.frame
}
\details{
we run \code{\link[dplyr:tbl]{dplyr::tbl()}}, then \code{\link[dplyr:compute]{dplyr::collect()}}
}
\examples{
\dontrun{
src <- src_itis()
sql_collect(src, "select * from hierarchy limit 5")
## or pipe the src to sql_collect
src \%>\% sql_collect("select * from hierarchy limit 5")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ap_name2taxid.R
\name{name2taxid}
\alias{name2taxid}
\title{Convert species names to taxon IDs}
\usage{
name2taxid(x, db = "ncbi", verbose = TRUE, out_type = c("uid", "summary"), ...)
}
\arguments{
\item{x}{(character) Vector of taxon keys for the given database}

\item{db}{(character) The database to search, one of ncbi, itis, gbif,
wfo, or tpl}

\item{verbose}{(logical) Print verbose messages}

\item{out_type}{(logical) character "uid" for an ID vector, "summary" for a
table with columns 'tax_id' and 'tax_name'.}

\item{...}{Additional arguments passed to database specific classification
functions.}
}
\description{
\code{name2taxid()} returns a vector and dies if there are any ambiguous
names. \code{name2taxid_map()} returns a data.frame mapping names to ids
}
\section{NCBI database}{


The NCBI taxonomy database includes common names, synonyms and misspellings.
However, the database is a little inconsistent. For some species, such as
Arabidopsis thaliana, the misspelling Arabidopsis_thaliana is included, but
the same is NOT done for humans. However, underscores are supported when
querying through entrez, as is done in taxize, which implies entrez is
replacing underscores with spaces. So I do the same. A corner case appears
when an organism uses underscores as part of the name, not just a standin
for space ("haloarchaeon 3A1_DGR"). To deal with this case, we replace
underscores with spaces ONLY if there are not spaces in the original name.
}

\examples{
\dontrun{
name2taxid(c('Arabidopsis thaliana', 'pig'))
name2taxid(c('Arabidopsis thaliana', 'pig'), out_type="summary")
name2taxid(x=c('Arabidopsis thaliana', 'Apis mellifera'), db = "itis")
name2taxid(x=c('Arabidopsis thaliana', 'Apis mellifera'), db = "itis",
 out_type="summary")
name2taxid(x=c('Arabidopsis thaliana', 'Quercus kelloggii'), db = "wfo")
name2taxid(x=c('Arabidopsis thaliana', 'Quercus kelloggii'), db = "wfo",
 out_type="summary")
name2taxid("Austrobaileyaceae", db = "wfo")
name2taxid("Quercus kelloggii", db = "gbif")
name2taxid(c("Quercus", "Fabaceae", "Animalia"), db = "gbif")
name2taxid(c("Abies", "Pinales", "Tracheophyta"), db = "col")
name2taxid(c("Abies mangifica", "Acanthopale aethiogermanica",
  "Acanthopale albosetulosa"), db = "tpl")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxizedb-package.R
\docType{package}
\name{taxizedb-package}
\alias{taxizedb-package}
\alias{taxizedb}
\title{taxizedb}
\description{
Taxonomic databases interface
}
\section{Supported data sources and database structure}{

All are using SQLite as the database
\itemize{
\item NCBI: text files are provided by NCBI, which we stitch into a sqlite db
\item ITIS: they provide a sqlite dump, which we use here
\item The PlantList: created from stitching together csv files. this
source is no longer updated as far as we can tell. they say they've
moved focus to the World Flora Online
\item Catalogue of Life: created from Darwin Core Archive dump. Using the
latest monthly edition via
http://www.catalogueoflife.org/DCA_Export/archive.php
\item GBIF: created from Darwin Core Archive dump. right now we only have
the taxonomy table (called gbif), but will add the other tables in the
darwin core archive later
\item Wikidata: aggregated taxonomy of Open Tree of Life, GLoBI and Wikidata.
On Zenodo, created by Joritt Poelen of GLOBI.
\item World Flora Online: http://www.worldfloraonline.org/
}
}

\section{Update schedule for databases}{

\itemize{
\item NCBI: since \code{db_download_ncbi} creates the database when the function
is called, it's updated whenever you run the function
\item ITIS: since ITIS provides the sqlite database as a download, you can
delete the old file and run \code{db_download_itis} to get a new dump;
they I think update the dumps every month or so
\item The PlantList: no longer updated, so you shouldn't need to download
this after the first download
\item Catalogue of Life: a GitHub Actions job runs once a day at 00:00 UTC,
building the lastest COL data into a SQLite database thats hosted on
Amazon S3
\item GBIF: a GitHub Actions job runs once a day at 00:00 UTC,
building the lastest COL data into a SQLite database thats hosted on
Amazon S3
\item Wikidata: last updated April 6, 2018. Scripts are available to
update the data if you prefer to do it yourself.
\item World Flora Online: since \code{db_download_wfo} creates the database when
the function is called, it's updated whenever you run the function
}
}

\section{Links}{

\itemize{
\item NCBI: ftp://ftp.ncbi.nih.gov/pub/taxonomy/
\item ITIS: https://www.itis.gov/downloads/index.html
\item The PlantList - http://www.theplantlist.org/
\item Catalogue of Life:
via http://www.catalogueoflife.org/content/annual-checklist-archive
\item GBIF: http://rs.gbif.org/datasets/backbone/
\item Wikidata: https://zenodo.org/record/1213477
\item World Flora Online: http://www.worldfloraonline.org/
}
}

\examples{
\dontrun{
library(dplyr)

# data source: NCBI
db_download_ncbi()
src <- src_ncbi()
df <- tbl(src, "names")
filter(df, name_class == "scientific name")

# data source: ITIS
## download ITIS database
db_download_itis()
## connect to the ITIS database
src <- src_itis()
## use SQL syntax
sql_collect(src, "select * from hierarchy limit 5")
### or pipe the src to sql_collect
src \%>\% sql_collect("select * from hierarchy limit 5")
## use dplyr verbs
src \%>\%
  tbl("hierarchy") \%>\%
  filter(ChildrenCount > 1000)
## or create tbl object for repeated use
hiers <- src \%>\% tbl("hierarchy")
hiers \%>\% select(TSN, level)

# data source: The PlantList
## download tpl datababase
db_download_tpl()
## connecto the tpl database
src <- src_tpl()
## do queries
tpl <- tbl(src, "tpl")
filter(tpl, Family == "Pinaceae")

# data source: Catalogue of Life
## download col datababase
db_download_col()
## connec to the col database
src <- src_col()
## do queries
names <- tbl(src, "taxa")
select(names, taxonID, scientificName)

# data source: GBIF
## download gbif datababase
db_download_gbif()
## connecto the gbif database
src <- src_gbif()
## do queries
df <- tbl(src, "gbif")
select(df, taxonID, scientificName)

# data source: Wikidata
db_download_wikidata()
src <- src_wikidata()
df <- tbl(src, "wikidata")
filter(df, rank_id == "Q7432")

# data source: World Flora Online
db_download_wfo()
src <- src_wfo()
df <- tbl(src, "wfo")
filter(df, taxonID == "wfo-0000000010")
}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ap_taxid2name.R
\name{taxid2name}
\alias{taxid2name}
\title{Convert taxon IDs to scientific names}
\usage{
taxid2name(x, db = "ncbi", verbose = TRUE, warn = TRUE, ...)
}
\arguments{
\item{x}{(character) Vector of taxon keys for the given database}

\item{db}{(character) The database to search, one of ncbi, itis, gbif,
col, wfo, or tpl}

\item{verbose}{(logical) Print verbose messages}

\item{warn}{(logical) If \code{TRUE}, raise a warning if any taxon IDs can not
be found}

\item{...}{Additional arguments passed to database specific classification
functions}
}
\value{
character vector of scientific names
}
\description{
Convert taxon IDs to scientific names
}
\examples{
\dontrun{
taxid2name(c(3702, 9606))
taxid2name(c(154395, 154357, 23041, 154396), db = "itis")
taxid2name(c('wfo-0000541830', 'wfo-0000291463'), db = "wfo")
taxid2name("wfo-7000000057", db="wfo")
taxid2name(2877951, db="gbif")
taxid2name(c(2877951, 5386), db="gbif")
taxid2name(c(3960765, 3953606, 3953010), db="col")
taxid2name(c("kew-2614538", "kew-2895433", "kew-2615007"), db="tpl")
}
}
