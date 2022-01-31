
# arkdb <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![R build
status](https://github.com/ropensci/arkdb/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/arkdb/actions)
[![Travis build
status](https://travis-ci.org/ropensci/arkdb.svg?branch=master)](https://travis-ci.org/ropensci/arkdb)
[![Coverage
status](https://codecov.io/gh/ropensci/arkdb/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/arkdb?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/arkdb)](https://cran.r-project.org/package=arkdb)
[![](https://badges.ropensci.org/224_status.svg)](https://github.com/ropensci/software-review/issues/224)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
[![CRAN RStudio mirror
downloads](http://cranlogs.r-pkg.org/badges/grand-total/arkdb)](https://CRAN.R-project.org/package=arkdb)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1343943.svg)](https://doi.org/10.5281/zenodo.1343943)
<!-- badges: end -->

<!-- README.md is generated from README.Rmd. Please edit that file -->

The goal of `arkdb` is to provide a convenient way to move data from
large compressed text files (tsv, csv, etc) into any DBI-compliant
database connection (e.g. MYSQL, Postgres, SQLite; see
[DBI](https://db.rstudio.com/dbi/)), and move tables out of such
databases into text files. The key feature of `arkdb` is that files are
moved between databases and text files in chunks of a fixed size,
allowing the package functions to work with tables that would be much
too large to read into memory all at once. There is also functionality
for filtering and applying transformation to data as it is extracted
from the database.

The `arkdb` package is easily extended to use custom read and write
methods allowing you to dictate your own output formats. See
`R/streamable_table.R` for examples that include using:

-   Base c/tsv
-   Apache arrow’s parquet
-   The `readr` package for c/tsv

## Links

-   A more detailed introduction to package design and use can be found
    in the package
    [Vignette](https://docs.ropensci.org/arkdb/articles/arkdb.html)
-   [Online versions of package
    documentation](https://docs.ropensci.org/arkdb/)

## Installation

You can install arkdb from GitHub with:

``` r
# install.packages("devtools")
devtools::install_github("cboettig/arkdb")
```

# Basic use

``` r
library(arkdb)

# additional libraries just for this demo
library(dbplyr)
library(dplyr)
library(fs)
```

## Creating an archive of a database

Consider the `nycflights` database in SQLite:

``` r
tmp <- tempdir() # Or can be your working directory, "."
db <- dbplyr::nycflights13_sqlite(tmp)
#> Caching nycflights db at /tmp/RtmpKGu2Ay/nycflights13.sqlite
#> Creating table: airlines
#> Creating table: airports
#> Creating table: flights
#> Creating table: planes
#> Creating table: weather
```

Create an archive of the database:

``` r
dir <- fs::dir_create(fs::path(tmp, "nycflights"))
ark(db, dir, lines = 50000)
#> Exporting airlines in 50000 line chunks:
#>  ...Done! (in 0.006583929 secs)
#> Exporting airports in 50000 line chunks:
#>  ...Done! (in 0.02108455 secs)
#> Exporting flights in 50000 line chunks:
#>  ...Done! (in 8.810824 secs)
#> Exporting planes in 50000 line chunks:
#>  ...Done! (in 0.02794719 secs)
#> Exporting weather in 50000 line chunks:
#>  ...Done! (in 0.6644697 secs)
```

## Unarchive

Import a list of compressed tabular files (i.e. `*.csv.bz2`) into a
local SQLite database:

``` r
files <- fs::dir_ls(dir)
new_db <- DBI::dbConnect(RSQLite::SQLite(), fs::path(tmp, "local.sqlite"))

unark(files, new_db, lines = 50000)
#> Importing /tmp/RtmpKGu2Ay/nycflights/airlines.tsv.bz2 in 50000 line chunks:
#>  ...Done! (in 0.0117662 secs)
#> Importing /tmp/RtmpKGu2Ay/nycflights/airports.tsv.bz2 in 50000 line chunks:
#>  ...Done! (in 0.02637362 secs)
#> Importing /tmp/RtmpKGu2Ay/nycflights/flights.tsv.bz2 in 50000 line chunks:
#>  ...Done! (in 6.802646 secs)
#> Importing /tmp/RtmpKGu2Ay/nycflights/planes.tsv.bz2 in 50000 line chunks:
#>  ...Done! (in 0.03848696 secs)
#> Importing /tmp/RtmpKGu2Ay/nycflights/weather.tsv.bz2 in 50000 line chunks:
#>  ...Done! (in 0.3772023 secs)
```

## Using filters

This package can also be used to generate slices of data that are
required for analytical or operational purposes. In the example below we
archive to disk only the flight data that occured in the month of
December. It is recommended to use filters on a single table at a time.

``` r
ark(db, dir, lines = 50000, tables = "flights", filter_statement = "WHERE month = 12")
```

## Using callbacks

It is possible to use a callback to perform just-in-time data
transformations before ark writes your data object to disk in your
preferred format. In the example below, we write a simple transformation
to convert the flights data `arr_delay` field, from minutes, to hours.
It is recommended to use callbacks on a single table at a time. A
callback function can be anything you can imagine so long as it returns
a data.frame that can be written to disk.

``` r
mins_to_hours <- function(data) {
  data$arr_delay <- data$arr_delay/60
  data
}

ark(db, dir, lines = 50000, tables = "flights", callback = mins_to_hours)
```

## ark() in parallel

There are two strategies for using `ark` in parallel. One is to loop
over the tables, re-using the ark function per table in parallel. The
other, introduced in 0.0.15, is to use the “window-parallel” method
which loops over chunks of your table. This is particularly useful if
your tables are very large and can speed up the process significantly.

Note: `window-parallel` currently only works in conjunction with
`streamable_parquet`

``` r
# Strategy 1: Parallel over tables
library(arkdb)
library(future.apply)

plan(multisession)

# Any streamable_table method is acceptable
future_lapply(vector_of_tables, function(x) ark(db, dir, lines, tables = x))

# Strategy 2: Parallel over chunks of a table
library(arkdb)
library(future.apply)

plan(multisession)

ark(
  db, 
  dir, 
  streamable_table = streamable_parquet(), # required for window-parallel
  lines = 50000, 
  tables = "flights", 
  method = "window-parallel"
)

# Strategy 3: Parallel over tables and chunks of tables
library(arkdb)
library(future.apply)
# 16 core machine for example
plan(list(tweak(multisession, n = 4), tweak(multisession, n = 4)))

# 4 tables at a time, 4 threads per table
future_lapply(vector_of_tables, function(x) { 
  ark(
    db, 
    dir, 
    streamable_table = streamable_parquet(), # required for window-parallel
    lines = 50000, 
    tables = x, 
    method = "window-parallel")
  }
)
```

## ETLs with arkdb

The `arkdb` package can also be used to create a number of ETL pipelines
involving text archives or databases given its ability to filter, and
use callbacks. In the example below, we leverage `duckdb` to read a
fictional folder of files by US state, filter by `var_filtered`, apply a
callback transformation `transform_fun` to `var_transformed` save as
parquet, and then load a folder of parquet files for analysis with
Apache Arrow.

``` r
library(arrow)
library(duckdb)

db <- dbConnect(duckdb::duckdb())

transform_fun <- function(data) {
  data$var_transformed <- sqrt(data$var_transformed)
  data
}

for(state in c("DC", state.abb)) {
  path <- paste0("path/to/archives/", state, ".gz")
  
  ark(
    db,
    dir = paste0("output/", state),
    streamable_table = streamable_parquet(), # parquet files of nline rows
    lines = 100000,
    # See: https://duckdb.org/docs/data/csv
    tables = sprintf("read_csv_auto('%s')", path), 
    compress = "none", # Compression meaningless for parquet as it's already compressed
    overwrite = T, 
    filenames = state, # Overload tablename
    filter_statement = "WHERE var_filtered = 1",
    callback = transform_fun
  )
}

# The result is trivial to read in with arrow 
ds <- open_dataset("output", partitioning = "state")
```

------------------------------------------------------------------------

Please note that this project is released with a [Contributor Code of
Conduct](https://ropensci.org/code-of-conduct/). By participating in
this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)
# arkdb 0.0.15

- Added window-parallel option for ark'ing large tables in parallel

# arkdb 0.0.14

- Patch for test suite for Solaris. `arrow` package installs on Solaris, but
  functions do not actually run correctly since the C++ libraries have not
  been set up properly on Solaris. 


# arkdb 0.0.13

- Added ability to name output files directly.
- Add warning when users specify compression for parquet files.
- Added callback functionality to the `ark` function. Allowing users to perform 
  transformations or recodes before chunked data.frames are saved to disk.
- Added ability to filter databases by allowing users to specify a "WHERE" clause. 
- Added parquet as an streamable_table format, allowing users to `ark` to parquet 
  instead of a text format. 

# arkdb 0.0.12

- Bugfix for arkdb

# arkdb 0.0.11

- Make cached connection opt-out instead of applying only to read_only.  This
  allows cache to work on read-write connections by default.  This also avoids
  the condition of a connection being garbage-collected when functions call
  local_db internally.

# arkdb 0.0.10

- Better handling of read_only vs read_write connections.  Only caches
  read_only connections.  
- Includes optional support for MonetDBLite

# arkdb 0.0.8

- Bugfix for dplyr 2.0.0 release


# arkdb 0.0.7

- Bugfix for upcoming dplyr 2.0.0 release

# arkdb 0.0.6

- Support vroom as an opt-in streamable table
- Export `process_chunks`
- Add mechanism to attempt a bulk importer, when available (#27)
- Bugfix for case when text contains `#` characters in base parser (#28)
- Lighten core dependencies.  Fully recursive dependencies include only 4
  non-base packages now, as `progress` is now optional.
- Use "magic numbers" instead of extensions to guess compression type.
  (NOTE: requires that file is local and not a URL)
- Now that `duckdb` is on CRAN and `MonetDBLite` isn't, drop built-in
  support for `MonetDBLite` in favor of `duckdb` alone.

# arkdb 0.0.5 2018-10-31

- `ark()`'s default `keep-open` method would cut off header names for
   Postgres connections (due to variation in the behavior of SQL queries
   with `LIMIT 0`.)  The issue is now resolved by accessing the header in
   a more robust, general way.

# arkdb 0.0.4 2018-09-27

- `unark()` will strip out non-compliant characters in table names by default.
- `unark()` gains the optional argument `tablenames`, allowing the user to
   specify the corresponding table names manually, rather than enforcing
   they correspond with the incoming file names. 
   [#18](https://github.com/ropensci/arkdb/issues/18)
-  `unark()` gains the argument `encoding`, allowing users to directly set
   the encoding of incoming files.  Previously this could only be set by
   setting `options(encoding)`, which will still work as well. See
  `FAO.R` example in `examples` for an illustration.  
- `unark()` will now attempt to guess which streaming parser to use 
   (e.g `csv` or `tsv`) based on the file extension pattern, rather than
   defaulting to a `tsv` parser.  (`ark()` still defaults to exporting in
   the more portable `tsv` format).

# arkdb 0.0.3 2018-09-11

* Remove dependency on utils::askYesNo for backward compatibility, [#17](https://github.com/ropensci/arkdb/issues/17)

# arkdb 0.0.2 2018-08-20 (First release to CRAN)

* Ensure the suggested dependency MonetDBLite is available before running unit test using it.

# arkdb 0.0.1 2018-08-20

* Overwrite existing tables of same name (with warning and
  interactive proceed) in both DB and text-files to avoid
  appending.

# arkdb 0.0.0.9000 2018-08-11

* Added a `NEWS.md` file to track changes to the package.
* Log messages improved as suggested by @richfitz
* Improved mechanism for windowing in most DBs, from @krlmlr [#8](https://github.com/ropensci/arkdb/pull/8)
* Support pluggable I/O, based on @richfitz suggestions [#3](https://github.com/ropensci/arkdb/issues/3), [#10](https://github.com/ropensci/arkdb/pull/10)

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
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/
Dear CRAN maintainers,

Changes in this release are described in NEWS.md.

This release (0.0.14) follows right behind the earlier release (0.0.13), 
because our previous release introduced new optional features that leverage
the arrow library, and these tests fail on Solaris. We apologize for not
catching this, but our tests did assert that they should only be run if 
arrow was available.  However, the arrow package is not properly installed 
and configured on solaris, as described in the error message that results:
as described in https://arrow.apache.org/docs/r/articles/install.html 

Until the CRAN Solaris machine is properly set up to use the arrow package,
we now gracefully skip the tests on that machine. 
 





Thanks!

Note that winbuilder will always throw a NOTE on this package due to the continued
use (as SUGGESTED only) of a CRAN package that is not in the standard repositories.  

winbuilder may also show a NOTE regarding a possible 304 ("Not Modified") code on
https://www.iana.org/assignments/media-types/text/tab-separated-values, the 
canonical IANA link defining this popular media-type format.  I do not believe there
is a preferable link here, but if the 304 code is problematic to CRAN maintainers,
I would be willing to merely remove the hyperlink and provide the URL only in 
plain text.  Just advise on your preference.


Carl




---
output: github_document
---

# arkdb <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![R build status](https://github.com/ropensci/arkdb/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/arkdb/actions)
[![Travis build status](https://travis-ci.org/ropensci/arkdb.svg?branch=master)](https://travis-ci.org/ropensci/arkdb)
[![Coverage status](https://codecov.io/gh/ropensci/arkdb/branch/master/graph/badge.svg)](https://codecov.io/github/ropensci/arkdb?branch=master)
[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/arkdb)](https://cran.r-project.org/package=arkdb)
[![](https://badges.ropensci.org/224_status.svg)](https://github.com/ropensci/software-review/issues/224)
[![lifecycle](https://img.shields.io/badge/lifecycle-stable-brightgreen.svg)](https://lifecycle.r-lib.org/articles/stages.html)
 [![CRAN RStudio mirror downloads](http://cranlogs.r-pkg.org/badges/grand-total/arkdb)](https://CRAN.R-project.org/package=arkdb) 
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.1343943.svg)](https://doi.org/10.5281/zenodo.1343943)
  <!-- badges: end -->
  
<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "README-"
)
```

The goal of `arkdb` is to provide a convenient way to move data from large compressed text files (tsv, csv, etc) into any DBI-compliant database connection (e.g. MYSQL, Postgres, SQLite; see [DBI](https://db.rstudio.com/dbi/)), and move tables out of such databases into text files. The key feature of `arkdb` is that files are moved between databases and text files in chunks of a fixed size, allowing the package functions to work with tables that would be much too large to read into memory all at once. There is also functionality for filtering and applying transformation to data as it is extracted from the database.  

The `arkdb` package is easily extended to use custom read and write methods allowing you to dictate your own output formats. See `R/streamable_table.R` for examples that include using: 

- Base c/tsv
- Apache arrow's parquet
- The `readr` package for c/tsv

## Links

- A more detailed introduction to package design and use can be found in the package [Vignette](https://docs.ropensci.org/arkdb/articles/arkdb.html)
- [Online versions of package documentation](https://docs.ropensci.org/arkdb/)

## Installation

You can install arkdb from GitHub with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("cboettig/arkdb")
```


# Basic use

```{r message = FALSE}
library(arkdb)

# additional libraries just for this demo
library(dbplyr)
library(dplyr)
library(fs)
```

## Creating an archive of a database

Consider the `nycflights` database in SQLite:

```{r example}
tmp <- tempdir() # Or can be your working directory, "."
db <- dbplyr::nycflights13_sqlite(tmp)
```

Create an archive of the database: 

```{r}
dir <- fs::dir_create(fs::path(tmp, "nycflights"))
ark(db, dir, lines = 50000)
```

## Unarchive

Import a list of compressed tabular files (i.e. `*.csv.bz2`) into a local SQLite database:


```{r}
files <- fs::dir_ls(dir)
new_db <- DBI::dbConnect(RSQLite::SQLite(), fs::path(tmp, "local.sqlite"))

unark(files, new_db, lines = 50000)
```



```{r include=FALSE}
disconnect <- function(db){
  ## Cleanup 
  if(inherits(db, "DBIConnection")){
    DBI::dbDisconnect(db)
  } else {
    DBI::dbDisconnect(db$con)
  }
}

DBI::dbDisconnect(db)
DBI::dbDisconnect(new_db)

codemeta::write_codemeta()
```

## Using filters

This package can also be used to generate slices of data that are required for analytical or operational purposes. In the example below we archive to disk only the flight data that occured in the month of December. It is recommended to use filters on a single table at a time. 

```{r, eval=FALSE}
ark(db, dir, lines = 50000, tables = "flights", filter_statement = "WHERE month = 12")
```

## Using callbacks

It is possible to use a callback to perform just-in-time data transformations before ark writes your data object to disk in your preferred format. In the example below, we write a simple transformation to convert the flights data `arr_delay` field, from minutes, to hours. It is recommended to use callbacks on a single table at a time. A callback function can be anything you can imagine so long as it returns a data.frame that can be written to disk.

```{r, eval=FALSE}
mins_to_hours <- function(data) {
  data$arr_delay <- data$arr_delay/60
  data
}

ark(db, dir, lines = 50000, tables = "flights", callback = mins_to_hours)
```

## ark() in parallel

There are two strategies for using `ark` in parallel. One is to loop over the tables, re-using the ark function per table in parallel. The other, introduced in 0.0.15, is to use the "window-parallel" method which loops over chunks of your table. This is particularly useful if your tables are very large and can speed up the process significantly. 

Note: `window-parallel` currently only works in conjunction with `streamable_parquet`

```{r, eval = FALSE}
# Strategy 1: Parallel over tables
library(arkdb)
library(future.apply)

plan(multisession)

# Any streamable_table method is acceptable
future_lapply(vector_of_tables, function(x) ark(db, dir, lines, tables = x))

# Strategy 2: Parallel over chunks of a table
library(arkdb)
library(future.apply)

plan(multisession)

ark(
  db, 
  dir, 
  streamable_table = streamable_parquet(), # required for window-parallel
  lines = 50000, 
  tables = "flights", 
  method = "window-parallel"
)

# Strategy 3: Parallel over tables and chunks of tables
library(arkdb)
library(future.apply)
# 16 core machine for example
plan(list(tweak(multisession, n = 4), tweak(multisession, n = 4)))

# 4 tables at a time, 4 threads per table
future_lapply(vector_of_tables, function(x) { 
  ark(
    db, 
    dir, 
    streamable_table = streamable_parquet(), # required for window-parallel
    lines = 50000, 
    tables = x, 
    method = "window-parallel")
  }
)

```

## ETLs with arkdb

The `arkdb` package can also be used to create a number of ETL pipelines involving text archives or databases given its ability to filter, and use callbacks. In the example below, we leverage `duckdb` to read a fictional folder of files by US state, filter by `var_filtered`, apply a callback transformation `transform_fun` to `var_transformed` save as parquet, and then load a folder of parquet files for analysis with Apache Arrow. 

```{r, eval = FALSE}
library(arrow)
library(duckdb)

db <- dbConnect(duckdb::duckdb())

transform_fun <- function(data) {
  data$var_transformed <- sqrt(data$var_transformed)
  data
}

for(state in c("DC", state.abb)) {
  path <- paste0("path/to/archives/", state, ".gz")
  
  ark(
    db,
    dir = paste0("output/", state),
    streamable_table = streamable_parquet(), # parquet files of nline rows
    lines = 100000,
    # See: https://duckdb.org/docs/data/csv
    tables = sprintf("read_csv_auto('%s')", path), 
    compress = "none", # Compression meaningless for parquet as it's already compressed
    overwrite = T, 
    filenames = state, # Overload tablename
    filter_statement = "WHERE var_filtered = 1",
    callback = transform_fun
  )
}

# The result is trivial to read in with arrow 
ds <- open_dataset("output", partitioning = "state")
```

-----

Please note that this project is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/).
By participating in this project you agree to abide by its terms.


[![ropensci_footer](https://ropensci.org/public_images/ropensci_footer.png)](https://ropensci.org)

---
title: "Introduction to arkdb"
author: "Carl Boettiger"
date: "`r Sys.Date()`"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{arkdb}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

```{r setup, include = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)
```

# arkdb


## Package rationale

Increasing data sizes create challenges for the fundamental tasks of publishing, distributing, and preserving data.  Despite (or perhaps because of) the diverse and ever-expanding number of database and file formats, the humble plain text file such as  comma or tab-separated-values (e.g. `.csv` or `.tsv` files) remains the gold standard for data archiving and distribution.  These files can read on almost any platform or tool and can be efficiently compressed using long-standing and widely available standard open source libraries like `gzip` or `bzip2`.  In contrast, database storage formats and dumps are usually particular to the database platform used to generate them, and will likely not be compatible between different database engines (e.g. PostgreSQL -> SQLite) or even between different versions of the same engine. Researchers unfamiliar with these databases will have difficulty accessing such data, and these dumps may also be in formats that are less efficient to compress.    

Working with tables that are too large for working memory on most machines by using external relational database stores is now a common R practice, thanks to ever-rising availability of data and increasing support and popularity of packages such as `DBI`, `dplyr`, and `dbplyr`.  Working with plain text files becomes increasingly difficult in this context.  Many R users will not have sufficient RAM to simply read in a 10 GB `.tsv` file into R.  Similarly, moving a 10 GB database out of a relational data file and into a plain text file for archiving and distribution is similarly challenging from R. While most relational database back-ends implement some form of `COPY` or `IMPORT` that allows them to read in and export out plain text files directly, these methods are not consistent across database types and not part of the standard SQL interface.  Most importantly for our case, they also cannot be called directly from R, but require a separate stand-alone installation of the database client.  `arkdb` provides a simple solution to these two tasks. 
 

The goal of `arkdb` is to provide a convenient way to move data from large compressed text files (e.g. `*.tsv.bz2`) into any DBI-compliant database connection (see [DBI](https://db.rstudio.com/dbi/)), and move tables out of such databases into text files. The key feature of `arkdb` is that files are moved between databases and text files in chunks of a fixed size, allowing the package functions to work with tables that would be much to large to read into memory all at once.  This will be slower than reading the file into memory at one go, but can be scaled to larger data and larger data with no additional memory requirement. 


## Installation

You can install `arkdb` from GitHub with:

```{r gh-installation, eval = FALSE}
# install.packages("devtools")
devtools::install_github("cboettig/arkdb")
```


# Tutorial

```{r message = FALSE}
library(arkdb)

# additional libraries just for this demo
library(dbplyr)
library(dplyr)
library(nycflights13)
library(fs)
```

## Creating an archive of an existing database

First, we'll need an example database to work with.  Conveniently, there is a nice example using the NYC flights data built into the `dbplyr` package.

```{r example}
tmp <- tempdir() # Or can be your working directory, "."
db <- dbplyr::nycflights13_sqlite(tmp)
```

To create an archive, we just give `ark` the connection to the database and tell it where we want the `*.tsv.bz2` files to be archived.   We can also set the chunk size as the number of `lines` read in a single chunk.  More lines per chunk usually means faster run time at the cost of higher memory requirements. 

```{r}
dir <- fs::dir_create(fs::path(tmp, "nycflights"))
ark(db, dir, lines = 50000)


```

We can take a look and confirm the files have been written. Note that we can use `fs::dir_info` to get a nice snapshot of the file sizes.  Compare the compressed sizes to the original database:

```{r}
fs::dir_info(dir) %>% 
  select(path, size) %>%
  mutate(path = fs::path_file(path))

fs::file_info(fs::path(tmp,"nycflights13.sqlite")) %>% 
  pull(size)


```


## Unarchive

Now that we've gotten all the database into (compressed) plain text files, let's get them back out.  We simply need to pass `unark` a list of these compressed files and a connection to the database.  Here we create a new local SQLite database.  Note that this design means that it is also easy to use `arkdb` to move data between databases.  


```{r}
files <- fs::dir_ls(dir, glob = "*.tsv.bz2")
new_db <- DBI::dbConnect(RSQLite::SQLite(), fs::path(tmp, "local.sqlite"))

```

As with `ark`, we can set the chunk size to control the memory footprint required:

```{r}
unark(files, new_db, lines = 50000)  
```

`unark` returns a `dplyr` database connection that we can use in the usual way:

```{r}
tbl(new_db, "flights")
```



```{r}
# Remove example files we created.
DBI::dbDisconnect(new_db)
unlink(dir, TRUE)
unlink(fs::path(tmp, "local.sqlite"))
```




## Pluggable text formats

 

By default, `arkdb` uses `tsv` format, implemented in base tools, as the text-based serialization.  The `tsv` standard is particularly attractive because it side-steps some of the ambiguities present in the CSV format due to string quoting.  The [IANA Standard for TSV](https://www.iana.org/assignments/media-types/text/tab-separated-values) neatly avoids this for tab-separated values by insisting that a tab can only ever be a separator.

`arkdb` provides a pluggable mechanism for changing the back end utility used to write text files. For instance, if we need to read in or export in `.csv` format, we can simply swap in a `csv` based reader in both `ark()` and `unark()` methods, as illustrated here:

```{r}
dir <- fs::dir_create(fs::path(tmp, "nycflights"))

ark(db, dir, 
    streamable_table = streamable_base_csv())
```




```{r}
files <- fs::dir_ls(dir, glob = "*.csv.bz2")
new_db <- DBI::dbConnect(RSQLite::SQLite(), fs::path(tmp, "local.sqlite"))

unark(files, new_db,
      streamable_table = streamable_base_csv())
```



`arkdb` also provides the function `streamable_table()` to facilitate users creating their own streaming table interfaces.  For instance, if you would prefer to use `readr` methods to read and write `tsv` files, we could construct the table as follows (`streamable_readr_tsv()` and `streamable_readr_csv()` are also shipped inside `arkdb` for convenience):

```{r}
stream <- 
   streamable_table(
     function(file, ...) readr::read_tsv(file, ...),
     function(x, path, omit_header)
       readr::write_tsv(x = x, path = path, append = omit_header),
     "tsv")

```

and we can then pass such a streamable table directly to `ark()` and `unark()`, like so:

```{r}
ark(db, dir, 
    streamable_table = stream)
```

Note several constraints on this design. The write method must be able to take a generic R `connection` object (which will allow it to handle the compression methods used, if any), and the read method must be able to take a `textConnection` object.  `readr` functions handle these cases out of the box, so the above method is easy to write.  Also note that the write method must be able to `append`, i.e. it should use a header if `append=TRUE`, but omit when it is `FALSE`.  See the built-in methods for more examples.


## A note on compression

`unark` can read from a variety of compression formats recognized by base R: `bzip2`, `gzip`, `zip`, and `xz`, and `ark` can choose any of these as the compression algorithm.  Note that there is some trade-off between speed of compression and efficiency (i.e. the final file size).  `ark` uses the `bz2` compression algorithm by default, supported in base R, to compress `tsv` files.  The  `bz2` offers excellent compression levels, but is considerably slower to compress than `gzip` or `zip`.  It is comparably fast to uncompress.  For faster archiving when maximum file size reduction is not critical, `gzip` will give nearly as effective compression in significantly less time.  Compression can also be turned off, e.g. by using `ark()` with `compress="none"` and `unark()` with files that have no compression suffix (e.g. `*.tsv` instead of `*.tsv.gz`). 


## Distributing data

Once you have archived your database files with `ark`, consider sharing them privately or publicly as part of your project GitHub repo using the [`piggyback` R package](https://github.com/ropensci/piggyback). For more permanent, versioned, and citable data archiving, upload your `*.tsv.bz2` files to a data repository like [Zenodo.org](https://zenodo.org).  




```{r include=FALSE}

disconnect <- function(db){
  ## Cleanup 
  if(inherits(db, "DBIConnection")){
    DBI::dbDisconnect(db)
  } else {
    DBI::dbDisconnect(db$con)
  }
}
disconnect(db)
DBI::dbDisconnect(new_db)
unlink(dir, TRUE)
```
---
title: "Working with medium-sized data"
author: "Carl Boettiger"
date: "9/22/2018"
output: pdf_document
vignette: >
  %\VignetteIndexEntry{working_with_data}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---

## *DRAFT POST*

Over the past summer, I have written two small-ish R packages to address challenges I frequently run up against during the course of my research.  Both are challenges with what I will refer to as medium-sized data -- not the kind of petabyte scale "big data" which precludes analysis on standard hardware or existing methodology, but large enough that the size alone starts creating problems for certain bits of a typical workflow.  More precisely, I will take *medium-sized* to refer to data that is too large to comfortably fit in memory on most laptops (e.g. on the order of several GB), or data that is merely too large to commit to GitHub.  By *typical workflow*, I mean easily being able to share all parts of analysis publicly or privately with collaborators (or merely different machines, such as my laptop and cloud server) who should be able to reproduce the results with minimal fuss and configuration.

For data too large to fit into memory, there's already a well-established solution of using an external database, to store the data.  Thanks to `dplyr`'s database backends, many R users can adapt their workflow relatively seamlessly to move from `dplyr` commands that call in-memory data frames to identical or nearly identical commands that call a database. This all works pretty well when your data *is already in a database*, but getting it onto a database, and then moving the data around so that other people/machines can access it is not nearly so straight forward. So far, this part of the problem has received relatively little attention.

The reason for that is because the usual response to this problem is "you're doing it wrong."  The standard practice in this context is simply not to move the data at all.  A central database server, usually with access controlled by password or other credential, can allow multiple users to all query the same database directly.  Thanks to the magical abstractions of SQL queries such as the `DBI` package, the user (aka client), doesn't need to care about the details of where the database is located, or even what particular backend is used. Moving all that data around can be slow and expensive. Arbitrarily large data can be housed in a central/cloud location and provisioned with enough resources to store everything and process complex queries. Consequently, just about every database backend not only to provides a mechanism for doing your `SQL` / `dplyr` querying, filtering, joining etc on data that cannot fit into memory all at once, but also nearly every such backend provides *server* abilities to do so over a network connection, handling secure logins and so forth.  Why would you want to do anything else?

The problem with the usual response is that it is often at odds with our original objectives and typical scientific workflows.  Setting up a database server can be non-trivial; by which I mean: difficult to automate in a portable/cross-platform manner when working entirely from R.  More importantly, it reflects a use-case more typical of industry context than scientific practice.  Individual researchers need to make data available to a global community of scientists who can reproduce results years or decades later; not just to a handful of employees who can be granted authenticated access to a central database.  Archiving data as static text files is far more *scalable*, more *cost-effective* (storing static files is much cheaper than keeping a database server running), more *future-proof* (rapid evolution in database technology is not always backwards compatible) and simplifies or *avoids most security issues* involved in maintaining a public server. In the scientific context, it almost always makes more sense to move the data after all.  

Scientific data repositories are already built on precisely this model: providing long term storage of files that can be downloaded and analyzed locally. For smaller `.csv` files, this works pretty well.



We just wanted to access data that was a bit larger than our active memory.  There is in fact a widely-used solution to this case: the `Lite` flavors of databases like `SQLite`, or my new favorite, `MonetDBLite`, which provide the disk-based storage but not the support for network connection server model.  Using the corresponding R packages, these databases can be easily deployed to store & query data on our local disk.  



% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/process_chunks.R
\name{process_chunks}
\alias{process_chunks}
\title{process a table in chunks}
\usage{
process_chunks(
  file,
  process_fn,
  streamable_table = NULL,
  lines = 50000L,
  encoding = Sys.getenv("encoding", "UTF-8"),
  ...
)
}
\arguments{
\item{file}{path to a file}

\item{process_fn}{a function of a \code{chunk}}

\item{streamable_table}{interface for serializing/deserializing in chunks}

\item{lines}{number of lines to read in a chunk.}

\item{encoding}{encoding to be assumed for input files.}

\item{...}{additional arguments to \code{streamable_table$read} method.}
}
\description{
process a table in chunks
}
\examples{
con <- system.file("extdata/mtcars.tsv.gz", package = "arkdb")
dummy <- function(x) message(paste(dim(x), collapse = " x "))
process_chunks(con, dummy, lines = 8)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_db.R
\name{arkdb_delete_db}
\alias{arkdb_delete_db}
\title{delete the local arkdb database}
\usage{
arkdb_delete_db(db_dir = arkdb_dir(), ask = interactive())
}
\arguments{
\item{db_dir}{neon database location}

\item{ask}{Ask for confirmation first?}
}
\description{
delete the local arkdb database
}
\details{
Just a helper function that deletes the database
files.  Usually unnecessary but can be
helpful in resetting a corrupt database.
}
\examples{

# Create a db
dir <- tempfile()
db <- local_db(dir)

# Delete it
arkdb_delete_db(dir, ask = FALSE)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/arkdb.R
\docType{package}
\name{arkdb-package}
\alias{arkdb}
\alias{arkdb-package}
\title{arkdb: Archive and Unarchive Databases Using Flat Files}
\description{
Flat text files provide a more robust, compressible,
and portable way to store tables.  This package provides convenient
functions for exporting tables from relational database connections
into compressed text files and streaming those text files back into
a database without requiring the whole table to fit in working memory.
}
\details{
It has two functions:
\itemize{
\item \code{\link[=ark]{ark()}}: archive a database into flat files, chunk by chunk.
\item \code{\link[=unark]{unark()}}: Unarchive flat files back int a database connection.
}

arkdb will work with any \code{DBI} supported connection.  This makes it
a convenient and robust way to migrate between different databases
as well.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://github.com/ropensci/arkdb}
  \item Report bugs at \url{https://github.com/ropensci/arkdb/issues}
}

}
\author{
\strong{Maintainer}: Carl Boettiger \email{cboettig@gmail.com} (\href{https://orcid.org/0000-0002-1642-628X}{ORCID}) [copyright holder]

Other contributors:
\itemize{
  \item Richard FitzJohn [contributor]
  \item Brandon Bertelsen \email{brandon@bertelsen.ca} [contributor]
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streamable_table.R
\name{streamable_table}
\alias{streamable_table}
\title{streamable table}
\usage{
streamable_table(read, write, extension)
}
\arguments{
\item{read}{read function. Arguments should be "\code{file}"
(must be able to take a \code{\link[=connection]{connection()}} object) and "\code{...}" (for)
additional arguments.}

\item{write}{write function. Arguments should be "\code{data}" (a data.frame),
\code{file} (must be able to take a \code{\link[=connection]{connection()}} object), and "\code{omit_header}"
logical, include header (initial write) or not (for appending subsequent
chunks)}

\item{extension}{file extension to use (e.g. "tsv", "csv")}
}
\value{
a \code{streamable_table} object (S3)
}
\description{
streamable table
}
\details{
Note several constraints on this design. The write method must be able
to take a generic R \code{connection} object (which will allow it to handle
the compression methods used, if any), and the read method must be able
to take a \code{textConnection} object.  \code{readr} functions handle these cases
out of the box, so the above method is easy to write.  Also note that
the write method must be able to \code{omit_header}. See the built-in methods
for more examples.
}
\examples{

streamable_readr_tsv <- function() {
  streamable_table(
    function(file, ...) readr::read_tsv(file, ...),
    function(x, path, omit_header) {
      readr::write_tsv(x = x, path = path, omit_header = omit_header)
    },
    "tsv"
  )
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_db.R
\name{local_db}
\alias{local_db}
\title{Connect to a local stand-alone database}
\usage{
local_db(
  dbdir = arkdb_dir(),
  driver = Sys.getenv("ARKDB_DRIVER", "duckdb"),
  readonly = FALSE,
  cache_connection = TRUE,
  memory_limit = getOption("duckdb_memory_limit", NA),
  ...
)
}
\arguments{
\item{dbdir}{Path to the database.}

\item{driver}{Default driver, one of "duckdb", "MonetDBLite", "RSQLite".
It will select the first one of those it finds available if a
driver is not set. This fallback can be overwritten either by explicit
argument or by setting the environmental variable \code{ARKDB_DRIVER}.}

\item{readonly}{Should the database be opened read-only? (duckdb only).
This allows multiple concurrent connections (e.g. from different R sessions)}

\item{cache_connection}{should we preserve a cache of the connection? allows
faster load times and prevents connection from being garbage-collected.  However,
keeping open a read-write connection to duckdb or MonetDBLite will block access of
other R sessions to the database.}

\item{memory_limit}{Set a memory limit for duckdb, in GB.  This can
also be set for the session by using options, e.g.
\code{options(duckdb_memory_limit=10)} for a limit of 10GB.  On most systems
duckdb will automatically set a limit to 80\% of machine capacity if not
set explicitly.}

\item{...}{additional arguments (not used at this time)}
}
\value{
Returns a \verb{[DBIconnection]} connection to the default database
}
\description{
This function will provide a connection to the best available database.
This function is a drop-in replacement for \verb{[DBI::dbConnect]} with behaviour
that makes it more subtle for R packages that need a database backend with
minimal complexity, as described in details.
}
\details{
This function provides several abstractions to \verb{[DBI::dbConnect]} to
provide a seamless backend for use inside other R packages.

First, this  provides a generic method that allows the use of a \verb{[RSQLite::SQLite]`` connection if nothing else is available, while being able to automatically select a much faster, more powerful backend from }\link[duckdb:duckdb]{duckdb::duckdb}`
if available.  An argument or environmental variable can be used to override this
to manually set a database endpoint for testing purposes.

Second, this function will cache the database connection in an R environment and
load that cache.  That means you can call \code{local_db()} as fast/frequently as you
like without causing errors that would occur by rapid calls to \verb{[DBI::dbConnect]}

Third, this function defaults to persistent storage location set by \verb{[tools::R_user_dir]}
and configurable by setting the environmental variable \code{ARKDB_HOME}.  This allows
a package to provide persistent storage out-of-the-box, and easily switch that storage
to a temporary directory (e.g. for testing purposes, or custom user configuration) without
having to edit database calls directly.
}
\examples{
\donttest{
## OPTIONAL: you can first set an alternative home location,
## such as a temporary directory:
Sys.setenv(ARKDB_HOME = tempdir())

## Connect to the database:
db <- local_db()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streamable_table.R
\name{streamable_vroom}
\alias{streamable_vroom}
\title{streamable tables using \code{vroom}}
\usage{
streamable_vroom()
}
\value{
a \code{streamable_table} object (S3)
}
\description{
streamable tables using \code{vroom}
}
\seealso{
\code{\link[readr:read_delim]{readr::read_tsv()}}, \code{\link[readr:write_delim]{readr::write_tsv()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streamable_table.R
\name{streamable_base_tsv}
\alias{streamable_base_tsv}
\title{streamable tsv using base R functions}
\usage{
streamable_base_tsv()
}
\value{
a \code{streamable_table} object (S3)
}
\description{
streamable tsv using base R functions
}
\details{
Follows the tab-separate-values standard using \code{\link[utils:read.table]{utils::read.table()}},
see IANA specification at:
\url{https://www.iana.org/assignments/media-types/text/tab-separated-values}
}
\seealso{
\code{\link[utils:read.table]{utils::read.table()}}, \code{\link[utils:write.table]{utils::write.table()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streamable_table.R
\name{streamable_parquet}
\alias{streamable_parquet}
\title{streamable chunked parquet using \code{arrow}}
\usage{
streamable_parquet()
}
\value{
a \code{streamable_table} object (S3)
}
\description{
streamable chunked parquet using \code{arrow}
}
\details{
Parquet files are streamed to disk by breaking them into chunks that are
equal to the \code{nlines} parameter in the initial call to \code{ark}. For each \code{tablename}, a
folder is created and the chunks are placed in the folder in the form \verb{part-000000.parquet}.
The software looks at the folder, and increments the name appropriately for the next
chunk. This is done intentionally so that users can take advantage of \code{arrow::open_dataset}
in the future, when coming back to review or perform analysis of these data.
}
\seealso{
\code{\link[arrow:read_parquet]{arrow::read_parquet()}}, \code{\link[arrow:write_parquet]{arrow::write_parquet()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streamable_table.R
\name{streamable_base_csv}
\alias{streamable_base_csv}
\title{streamable csv using base R functions}
\usage{
streamable_base_csv()
}
\value{
a \code{streamable_table} object (S3)
}
\description{
streamable csv using base R functions
}
\details{
Follows the comma-separate-values standard using \code{\link[utils:read.table]{utils::read.table()}}
}
\seealso{
\code{\link[utils:read.table]{utils::read.table()}}, \code{\link[utils:write.table]{utils::write.table()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streamable_table.R
\name{streamable_readr_csv}
\alias{streamable_readr_csv}
\title{streamable csv using \code{readr}}
\usage{
streamable_readr_csv()
}
\value{
a \code{streamable_table} object (S3)
}
\description{
streamable csv using \code{readr}
}
\seealso{
\code{\link[readr:read_delim]{readr::read_csv()}}, \code{\link[readr:write_delim]{readr::write_csv()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/streamable_table.R
\name{streamable_readr_tsv}
\alias{streamable_readr_tsv}
\title{streamable tsv using \code{readr}}
\usage{
streamable_readr_tsv()
}
\value{
a \code{streamable_table} object (S3)
}
\description{
streamable tsv using \code{readr}
}
\seealso{
\code{\link[readr:read_delim]{readr::read_tsv()}}, \code{\link[readr:write_delim]{readr::write_tsv()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ark.R
\name{ark}
\alias{ark}
\title{Archive tables from a database as flat files}
\usage{
ark(
  db_con,
  dir,
  streamable_table = streamable_base_tsv(),
  lines = 50000L,
  compress = c("bzip2", "gzip", "xz", "none"),
  tables = list_tables(db_con),
  method = c("keep-open", "window", "window-parallel", "sql-window"),
  overwrite = "ask",
  filter_statement = NULL,
  filenames = NULL,
  callback = NULL
)
}
\arguments{
\item{db_con}{a database connection}

\item{dir}{a directory where we will write the compressed text files output}

\item{streamable_table}{interface for serializing/deserializing in chunks}

\item{lines}{the number of lines to use in each single chunk}

\item{compress}{file compression algorithm. Should be one of "bzip2" (default),
"gzip" (faster write times, a bit less compression), "xz", or "none", for
no compression.}

\item{tables}{a list of tables from the database that should be
archived.  By default, will archive all tables. Table list should specify
schema if appropriate, see examples.}

\item{method}{method to use to query the database, see details.}

\item{overwrite}{should any existing text files of the same name be overwritten?
default is "ask", which will ask for confirmation in an interactive session, and
overwrite in a non-interactive script.  TRUE will always overwrite, FALSE will
always skip such tables.}

\item{filter_statement}{Typically an SQL "WHERE" clause, specific to your
dataset. (e.g., \verb{WHERE year = 2013})}

\item{filenames}{An optional vector of names that will be used to name the
files instead of using the tablename from the \code{tables} parameter.}

\item{callback}{An optional function that acts on the data.frame before it is
written to disk by \code{streamable_table}. It is recommended to use this on a single
table at a time. Callback functions must return a data.frame.}
}
\value{
the path to \code{dir} where output files are created (invisibly), for piping.
}
\description{
Archive tables from a database as flat files
}
\details{
\code{ark} will archive tables from a database as (compressed) tsv files.
Or other formats that have a \verb{streamtable_table method}, like parquet.
\code{ark} does this by reading only chunks at a time into memory, allowing it to
process tables that would be too large to read into memory all at once (which
is probably why you are using a database in the first place!)  Compressed
text files will likely take up much less space, making them easier to store and
transfer over networks.  Compressed plain-text files are also more archival
friendly, as they rely on widely available and long-established open source compression
algorithms and plain text, making them less vulnerable to loss by changes in
database technology and formats.

In almost all cases, the default method should be the best choice.
If the \code{\link[DBI:dbSendQuery]{DBI::dbSendQuery()}} implementation for your database platform returns the
full results to the client immediately rather than supporting chunking with \code{n}
parameter, you may want to use "window" method, which is the most generic.  The
"sql-window" method provides a faster alternative for databases like PostgreSQL that
support windowing natively (i.e. \code{BETWEEN} queries). Note that "window-parallel"
only works with \code{streamable_parquet}.
}
\examples{
\donttest{
# setup
library(dplyr)
dir <- tempdir()
db <- dbplyr::nycflights13_sqlite(tempdir())

## And here we go:
ark(db, dir)
}
\dontrun{

## For a Postgres DB with schema, we can append schema names first
## to each of the table names, like so:
schema_tables <- dbGetQuery(db, sqlInterpolate(db,
  "SELECT table_name FROM information_schema.tables
WHERE table_schema = ?schema",
  schema = "schema_name"
))

ark(db, dir, tables = paste0("schema_name", ".", schema_tables$table_name))
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/unark.R
\name{unark}
\alias{unark}
\title{Unarchive a list of compressed tsv files into a database}
\usage{
unark(
  files,
  db_con,
  streamable_table = NULL,
  lines = 50000L,
  overwrite = "ask",
  encoding = Sys.getenv("encoding", "UTF-8"),
  tablenames = NULL,
  try_native = TRUE,
  ...
)
}
\arguments{
\item{files}{vector of filenames to be read in. Must be \code{tsv}
format, optionally compressed using \code{bzip2}, \code{gzip}, \code{zip},
or \code{xz} format at present.}

\item{db_con}{a database src (\code{src_dbi} object from \code{dplyr})}

\item{streamable_table}{interface for serializing/deserializing in chunks}

\item{lines}{number of lines to read in a chunk.}

\item{overwrite}{should any existing text files of the same name be overwritten?
default is "ask", which will ask for confirmation in an interactive session, and
overwrite in a non-interactive script.  TRUE will always overwrite, FALSE will
always skip such tables.}

\item{encoding}{encoding to be assumed for input files.}

\item{tablenames}{vector of tablenames to be used for corresponding files.
By default, tables will be named using lowercase names from file basename with
special characters replaced with underscores (for SQL compatibility).}

\item{try_native}{logical, default TRUE. Should we try to use a native bulk
import method for the database connection?  This can substantially speed up
read times and will fall back on the DBI method for any table that fails
to import.  Currently only MonetDBLite connections support this.}

\item{...}{additional arguments to \code{streamable_table$read} method.}
}
\value{
the database connection (invisibly)
}
\description{
Unarchive a list of compressed tsv files into a database
}
\details{
\code{unark} will read in a files in chunks and
write them into a database.  This is essential for processing
large compressed tables which may be too large to read into
memory before writing into a database.  In general, increasing
the \code{lines} parameter will result in a faster total transfer
but require more free memory for working with these larger chunks.

If using \code{readr}-based streamable-table, you can suppress the progress bar
by using \code{options(readr.show_progress = FALSE)} when reading in large
files.
}
\examples{
\donttest{
## Setup: create an archive.
library(dplyr)
dir <- tempdir()
db <- dbplyr::nycflights13_sqlite(tempdir())

## database -> .tsv.bz2
ark(db, dir)

## list all files in archive (full paths)
files <- list.files(dir, "bz2$", full.names = TRUE)

## Read archived files into a new database (another sqlite in this case)
new_db <- DBI::dbConnect(RSQLite::SQLite())
unark(files, new_db)

## Prove table is returned successfully.
tbl(new_db, "flights")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/local_db.R
\name{local_db_disconnect}
\alias{local_db_disconnect}
\title{Disconnect from the arkdb database.}
\usage{
local_db_disconnect(db = local_db(), env = arkdb_cache)
}
\arguments{
\item{db}{a DBI connection. By default, will call \link{local_db} for the default connection.}

\item{env}{The environment where the function looks for a connection.}
}
\description{
Disconnect from the arkdb database.
}
\details{
This function manually closes a connection to the \code{arkdb} database.
}
\examples{
\donttest{

## Disconnect from the database:
local_db_disconnect()
}
}
