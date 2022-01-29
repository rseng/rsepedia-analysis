
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




