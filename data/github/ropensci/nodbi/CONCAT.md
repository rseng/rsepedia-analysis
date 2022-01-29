
# nodbi

[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran
checks](https://cranchecks.info/badges/worst/nodbi)](https://cranchecks.info/pkgs/nodbi)
[![R-CMD-check](https://github.com/ropensci/nodbi/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/nodbi/actions?query=workflow%3AR-CMD-check)
[![codecov](https://codecov.io/gh/ropensci/nodbi/branch/master/graph/badge.svg)](https://codecov.io/gh/ropensci/nodbi)
[![rstudio mirror
downloads](http://cranlogs.r-pkg.org/badges/nodbi)](https://github.com/r-hub/cranlogs.app)
[![cran
version](https://www.r-pkg.org/badges/version/nodbi)](https://cran.r-project.org/package=nodbi)

`nodbi` is an R package that provides a single interface for several
NoSQL databases, with the same function parameters and return values
across backends.

Currently, `nodbi` supports the following database backends:

-   MongoDB
-   SQLite
-   Elasticsearch
-   CouchDB
-   PostgreSQL (new in `nodbi` version 0.6.0)

for an `R` object of any of these data types:

-   data.frame
-   list
-   JSON string
-   a file name or URL with NDJSON records\*

and for executing the following operations:

-   Create
-   Exists
-   Get
-   Query\*\*
-   Update\*\*
-   Delete
-   List

across all database backends. \* `docdb_create`. \*\*Only simple
(e.g. equality for a single field) queries (and updates) are supported
for Elasticsearch at the moment. Only root fields can be specified for
CouchDB, whereas fields with subitems (in dot notation) can be specified
for Elasticsearch, MongoDB, SQLite and PostgreSQL.

For capabilities in terms of operations and parameter combinations
across any of the database backends, see section
[walk-through](#walk-through) below and see the main file for package
testing, here: [core-nodbi.R](./tests/testthat/core-nodbi.R).

## Install

CRAN version

``` r
install.packages("nodbi")
```

Development version

``` r
remotes::install_github("ropensci/nodbi")
```

Load package from library

``` r
library("nodbi")
```

## API overview

Parameters for `docdb_*()` functions are the same across database
backends.

| Purpose                                                                                                                         | Function call                                                                                                     |
|---------------------------------------------------------------------------------------------------------------------------------|-------------------------------------------------------------------------------------------------------------------|
| Create database connection (see below)                                                                                          | `src <- nodbi::src_{mongo, sqlite, couchdb, elastic}(<see below for parameters>)`                                 |
| Load `myData` (a data frame, a list, a JSON string, or a file name or url with NDJSON records) into database, container `myTbl` | `nodbi::docdb_create(src = src, key = "myTbl", value = myData)`                                                   |
| Get all documents back into a data frame                                                                                        | `nodbi::docdb_get(src = src, key = "myTbl")`                                                                      |
| Get documents selected with query (as MongoDB-compatible JSON) into a data frame                                                | `nodbi::docdb_query(src = src, key = "myTbl", query = '{"age": 20}')`                                             |
| Get selected fields (in MongoDB compatible JSON) from documents selected query                                                  | `nodbi::docdb_query(src = src, key = "myTbl", query = '{"age": 20}', fields = '{"name": 1, "_id": 0, "age": 1}')` |
| Update (patch) selected documents with new data in data frame, list or JSON string from `myData`                                | `nodbi::docdb_update(src = src, key = "myTbl", value = myData, query = '{"age": 20}')`                            |
| Check if container exists                                                                                                       | `nodbi::docdb_exists(src = src, key = "myTbl")`                                                                   |
| List all containers in database                                                                                                 | `nodbi::docdb_list(src = src)`                                                                                    |
| Delete document(s) in container                                                                                                 | `nodbi::docdb_delete(src = src, key = "myTbl", query = '{"age": 20}')`                                            |
| Delete container                                                                                                                | `nodbi::docdb_delete(src = src, key = "myTbl")`                                                                   |
| Close and remove database connection                                                                                            | `rm(src)`                                                                                                         |

## Database connections

Overview on parameters that are specific to the database backend. These
are only needed once, for `src_*()` to create a connection object for
use with `nodbi`.

### MongoDB

Note that only MongoDB requires to specify the container already in the
`src_*()` function. “Container” refers to a MongoDB collection.

``` r
src_mongo(
  collection = "test", db = "test",
  url = "mongodb://localhost", ...)
```

### SQLite

The functionality to process JSON is based on the SQLite extension
[JSON1](https://www.sqlite.org/json1.html), available in RSQLite.
“Container” refers to an SQLite table.

``` r
src_sqlite(dbname = ":memory:", ...)
```

### CouchDB

``` r
src_couchdb(
  host = "127.0.0.1", port = 5984L, path = NULL,
  transport = "http", user = NULL, pwd = NULL, headers = NULL)
```

### Elasticsearch

``` r
src_elastic(
  host = "127.0.0.1", port = 9200L, path = NULL,
  transport_schema = "http", user = NULL, pwd = NULL, force = FALSE, ...)
```

### PostgreSQL

With this database, the order of variables in data frames returned by
`docdb_get()` and `docdb_query()` can differ from the order in which
they were in `docdb_create()`.

``` r
src_postgres(
  dbname = "test", host = "127.0.0.1", port = 5432L, ...)
```

## Walk-through

This example is meant to show how functional `nodbi` is at this time.

``` r
# connect database backend
src <- src_sqlite()

# load data (here data frame, alternatively list or JSON)
docdb_create(src, key = "myTbl", value = mtcars)
#> [1] 32

# load additionally 98 NDJSON records
docdb_create(src, key = "myTbl", "http://httpbin.org/stream/98")
#> Note: container 'myTbl' already exists
#> [1] 98

# load additionally contacts JSON data, from package nodbi
docdb_create(src, key = "myTbl", contacts)
#> Note: container 'myTbl' already exists
#> [1] 5

# get all documents, irrespective of schema
dplyr::tibble(docdb_get(src, "myTbl"))
#> # A tibble: 37 × 22
#>    `_id`  isActive balance    age eyeColor name  email about  registered
#>    <chr>  <lgl>    <chr>    <int> <chr>    <chr> <chr> <chr>  <chr>     
#>  1 5cd67… TRUE     $2,412.…    20 blue     Kris… kris… Sint … 2017-07-1…
#>  2 5cd67… FALSE    $3,400.…    20 brown    Rae … raec… Nisi … 2018-12-1…
# ...
#>  9 Chrys… NA       NA          NA NA       NA    NA    NA     NA        
#> 10 Datsu… NA       NA          NA NA       NA    NA    NA     NA        
#> # … with 27 more rows, and 13 more variables: tags <list>,
#> #   friends <list>, mpg <dbl>, cyl <dbl>, disp <dbl>, hp <dbl>,
#> #   drat <dbl>, wt <dbl>, qsec <dbl>, vs <dbl>, am <dbl>, gear <dbl>,
#> #   carb <dbl>

# query some documents
docdb_query(src, "myTbl", query = '{"mpg": {"$gte": 30}}')
#>              _id mpg cyl disp  hp drat  wt qsec vs am gear carb
#> 1       Fiat 128  32   4   79  66  4.1 2.2   19  1  1    4    1
#> 2    Honda Civic  30   4   76  52  4.9 1.6   19  1  1    4    2
#> 3   Lotus Europa  30   4   95 113  3.8 1.5   17  9  1    5    2
#> 4 Toyota Corolla  34   4   71  65  4.2 1.8   20  1  1    4    1

# query some fields from some documents; 'query' is a mandatory 
# parameter and is used here in its position in the signature
docdb_query(src, "myTbl", '{"mpg": {"$gte": 30}}', fields = '{"wt": 1, "mpg": 1}')
#>    wt mpg
#> 1 2.2  32
#> 2 1.6  30
#> 3 1.5  30
#> 4 1.8  34

# query some subitem fields from some documents
# (only simple queries so far implemented for Elasticsearch)
# (only root, not subitems so far implemented for CouchDB)
str(docdb_query(src, "myTbl", '
                {"$or": [{"age": {"$lt": 21}}, 
                         {"friends.name": {"$regex": "^B[a-z]{3,6}.*"}}]}', 
                fields = '{"age": 1, "friends.name": 1}'))
#> 'data.frame':    4 obs. of  2 variables:
#>  $ age    : int  20 20 22 23
#>  $ friends:'data.frame': 4 obs. of  1 variable:
#>   ..$ name:List of 4
#>   .. ..$ : chr "Pace Bell"
#>   .. ..$ : chr  "Yang Yates" "Lacy Chen"
#>   .. ..$ : chr  "Baird Keller" "Francesca Reese" "Dona Bartlett"
#>   .. ..$ : chr  "Wooten Goodwin" "Brandie Woodward" "Angelique Britt"

# such queries can also be used for updating (patching) selected documents 
# with a new 'value'(s) from a JSON string, a data frame or a list
docdb_update(src, "myTbl", value = '{"vs": 9, "xy": [1, 2]}', query = '{"carb": 3}')
#> [1] 3
docdb_query(src, "myTbl", '{"carb": {"$in": [1,3]}}', fields = '{"vs": 1}')[[1]]
#> [1] 1 1 1 1 9 9 9 1 1 1
docdb_get(src, "myTbl")[126:130, c(1, 27, 28)]
#>                  _id carb   xy
#> 126        Merc 280C    4 NULL
#> 127       Merc 450SE    3 1, 2
#> 128       Merc 450SL    3 1, 2
#> 129      Merc 450SLC    3 1, 2
#> 130 Pontiac Firebird    2 NULL

# use with dplyr
library("dplyr")
docdb_get(src, "myTbl") %>%
  group_by(gear) %>%
  summarise(mean_mpg = mean(mpg))
# # A tibble: 4 × 2
#    gear mean_mpg
#   <dbl>    <dbl>
# 1     3     16.1
# 2     4     24.5
# 3     5     21.4
# 4    NA     NA  

# delete documents; query is optional parameter and has to be 
# specified for deleting documents instead of deleting the container
docdb_delete(src, "myTbl", query = '{"$or": {"gear": 5, "age": {"$gte": 22}}}')
#> TRUE
nrow(docdb_get(src, "myTbl"))
#> [1] 127

# delete container from database
docdb_delete(src, "myTbl")
#> [1] TRUE
```

## Benchmark

``` r
library("nodbi")

COUCHDB_TEST_USER <- Sys.getenv("COUCHDB_TEST_USER")
COUCHDB_TEST_PWD <- Sys.getenv("COUCHDB_TEST_PWD")

srcMongo <- src_mongo()
srcSqlite <- src_sqlite() # here, in-memory database
srcElastic <- src_elastic()
srcCouchdb <- src_couchdb(user = COUCHDB_TEST_USER, pwd = COUCHDB_TEST_PWD)
srcPostgres <- src_postgres()
key <- "test"
query <- '{"clarity": "SI1"}'
fields <- '{"cut": 1, "_id": 1, "clarity": "1"}'
value <- '{"clarity": "XYZ", "new": ["ABC", "DEF"]}'
data <- as.data.frame(diamonds)[1:12000, ]

testFunction <- function(src, key, value, query, fields) {
  docdb_create(src, key, data)
  # Elasticsearch needs a delay to process the data
  docdb_create(src, key, "http://httpbin.org/stream/89")
  # Elasticsearch needs a delay to process the data
  if (inherits(src, "src_elastic")) Sys.sleep(1)
  head(docdb_get(src, key))
  docdb_query(src, key, query = query, fields = fields)
  docdb_update(src, key, value = value, query = query)
  docdb_delete(src, key)
}

rbenchmark::benchmark(
  MongoDB = testFunction(src = srcMongo, key, value, query, fields),
  RSQLite = testFunction(src = srcSqlite, key, value, query, fields),
  Elastic = testFunction(src = srcElastic, key, value, query, fields),
  CouchDB = testFunction(src = srcCouchdb, key, value, query, fields),
  PostgreSQL = testFunction(src = srcPostgres, key, value, query, fields),
  replications = 10L,
  columns = c('test', 'replications', 'elapsed')
)
#> on 2015 mobile computer
#>         test replications elapsed
#> 4    CouchDB           10    1795
#> 3    Elastic           10     196 # 10s subtracted
#> 5 PostgreSQL           10      85
#> 2    RSQLite           10      76
#> 1    MongoDB           10      75
```

## Testing

``` r
testthat::test_local()
```

## Notes

-   Please [report any issues or
    bugs](https://github.com/ropensci/nodbi/issues).
-   License: MIT
-   Get citation information for `nodbi` in R doing
    `citation(package = 'nodbi')`
-   Please note that this package is released with a [Contributor Code
    of Conduct](https://ropensci.org/code-of-conduct/). By contributing
    to this project, you agree to abide by its terms.
-   Support for redis has been removed since version 0.5, because no way
    was found to query and update specific documents in a container.
nodbi 0.7.0.9000
================

### Changes
* new development version

nodbi 0.7.0
===========

### IMPROVEMENTS
* `docdb_create` now supports file names and http urls as argument `value` for importing data
* `docdb_create` (and thus `docdb_update`) now supports quantifiers (e.g., '[a-z]{2,3}') in regular expressions

### BUG FIXES
* for SQLite, return `FALSE` like other backends when using `docdb_delete` for a non-existing container (table, in the case of SQLite)
* better handle special characters and encodings under Windows

nodbi 0.6.0
===========

### IMPROVEMENTS
* full support for PostgreSQL (using jsonb)

### BUG FIXES
* for SQLite add closing file references also on exit

nodbi 0.5.1
===========

### BUG FIXES
* for SQLite under Windows ensure handling of special characters (avoiding encoding conversions with file operations that stream out / in NDJSON)

nodbi 0.5.0
===========

### IMPROVEMENTS
* identical API for `docdb_*()` functions so that `query` and `fields` parameters can be used across database backends
* identical return values across database backends

### UNDER THE HOOD
* re-factored recently added functions for RSQLite
* re-factored most functions to provide identical API
* performance (timing and memory use) profiled and optimised as far as possible

### OTHER CHANGES
* testing now uses the same test file across databases
* currently, no more support for redis (no way was found to query and update specific documents in a container)
* `docdb_list()` added as function to list container in database

### NOTES
* Support for complex queries not yet implemented for Elasticsearch
* Only root fields (no subitems) returned by Elasticsearch and CouchDB

nodbi 0.4.4
===========

### MINOR IMPROVEMENTS
* made remaining `docdb_*()` functions return a logical indicating the success of the function (`docdb_create`, `docdb_delete`), or a data frame (`docdb_get`, `docdb_query`), or the number of documents affected by the function (`docdb_update`)
  
### BUG FIXES
* `docdb_get()` to not return '_id' field for `src_{sqlite,mongo}` since already used for row names

### OTHER
* change testing approach

nodbi 0.4.3
===========

### MINOR IMPROVEMENTS
* `docdb_query.src_sqlite()` now handles JSON objects, returning nested lists (#40)
* `src_sqlite()` now uses transactions for relevant functions (#39)
* `docdb_update.src_mongo()` now returns the number of upserted or matched documents, irrespective of whether
  they were updated or not
  
### BUG FIXES
* `docdb_get()` to not return '_id' field for `src_{sqlite,mongo}` since already used for row names

### OTHER
* change of maintainer agreed

nodbi 0.4.2
===========

### BUG FIXES

* fix for `src_couchdb()`: we were not setting user and password correctly internally, was causing issues in CouchDB v3 (#35) thanks to @drtagkim for the pull request

nodbi 0.4.0
===========

### MINOR IMPROVEMENTS

* in `docdb_query()` and `docdb_get()`, for sqlite source, use a connection instead of a regular file path to avoid certain errors on Windows (#33) work by @rfhb
* in `docdb_query()` and `docdb_create()` for sqlite source, fix to handle mixed values of different types (#34) work by @rfhb
* some Sys.sleep's added to Elasticserch eg's to make sure data is available after creation, and before a data request

nodbi 0.3.0
===========

### NEW FEATURES

* new author Ralf Herold, with contribution of new functions for working with SQLite/json1. new functions: `src_sqlite`, `print.src_sqlite`, `docdb_create.src_sqlite`, `docdb_delete.src_sqlite`, `docdb_exists.src_sqlite`, `docdb_get.src_sqlite`, `docdb_query.src_sqlite`, and `docdb_update.src_sqlite`. includes new dataset `contacts` (#25) (#27) (#28) (#29) (#30) (#31)
* `docdb_update` gains method for working with MongoDB, via (#27)

### MINOR IMPROVEMENTS

* added `.github` files in the source repository to facilitate contributions
* `src_mongo` changes, improved behavior, via (#27)

### DEFUNCT

* `etcd` (via the `etseed` package) integration has been removed from this package as etcd doesn't really fit the main goal of the pkg. functions now defunct are: `src_etcd`, `docdb_create.src_etcd`, `docdb_delete.src_etcd`, `docdb_exists.src_etcd`, `docdb_get.src_etcd`, and `print.src_etcd` (#26)


nodbi 0.2.0
===========

### NEW FEATURES

* `docdb_get()` gains `limit` parameter to do pagination, for CouchDB, 
Elasticsearch and MongoDB only (#17) (#23)
* gains function `docdb_query()` to send queries to each backend (#18) (#22)
* gains function `docdb_exists()` to check if a database or equivalent exists (#21) (#22)

### MINOR IMPROVEMENTS

* Updated package for new version of `elastic`, which has slightly different
setup for connecting to the Elasticsearch instance (#20)


nodbi 0.1.0
===========

### NEW FEATURES

* released to CRAN
## Test environments

* Local: macOS, R 4.1.2 and R 3.6.3; with CouchDB, Elasticsearch, MongoDB, SQLite, PostgreSQL databases
* Github Actions: Ubuntu 20.04; R release and R devel
* win-builder: R Under development (unstable) (2022-01-03 r81439 ucrt)
* R-hub builder: Debian Linux, R-devel, clang, ISO-8859-15 locale; Windows Server 2022, R-devel, 64 bit
* macOS builder: r-release-macosx-arm64|4.1.1|macosx|macOS 11.5.2 (20G95)

## R CMD check results

0 errors | 0 warnings | 0 notes

## IMPROVEMENTS
* `docdb_create` now supports file names and http urls as argument `value` for importing data
* `docdb_create` (and thus `docdb_update`) now supports quantifiers (e.g., '[a-z]{2,3}') in regular expressions

## BUG FIXES
* for SQLite, return `FALSE` like other backends when using `docdb_delete` for a non-existing container (table, in the case of SQLite)
* better handle special characters and encodings under Windows

## revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

--------

Thank you -
Ralf Herold
<!-- Please use a feature branch (i.e., put your work in a new branch that has a name that reflects the feature you are working on; https://docs.gitlab.com/ee/workflow/workflow.html) -->

<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this pull request - most likely the maintainer will have their own equivalent key -->

<!-- If you've updated a file in the man-roxygen directory, make sure to update the man/ files by running devtools::document() or similar as .Rd files should be affected by your change -->

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

### Fixing typos

Small typos or grammatical errors in documentation may be edited directly using
the GitHub web interface, so long as the changes are made in the _source_ file.

*  YES: you edit a roxygen comment in a `.R` file below `R/`.
*  NO: you edit an `.Rd` file below `man/`.

### Prerequisites

Before you make a substantial pull request, you should always file an issue and
make sure someone from the team agrees that it’s a problem. If you’ve found a
bug, create an associated issue and illustrate the bug with a minimal 
[reprex](https://www.tidyverse.org/help/#reprex).

### Pull request process

*  We recommend that you create a Git branch for each pull request (PR).  
*  Look at the Travis and AppVeyor build status before and after making changes.
The `README` should contain badges for any continuous integration services used
by the package.  
*  We have some style tips in our [developer guide](https://devguide.ropensci.org/building.html#code-style)
*  We use [roxygen2](https://cran.r-project.org/package=roxygen2).  
*  We use [testthat](https://cran.r-project.org/package=testthat). Contributions with test cases included are easier to accept.  
*  For user-facing changes, add a bullet to the top of `NEWS.md` below the
current development version header describing the changes made followed by your
GitHub username, and links to relevant issue(s)/PR(s).

### Code of Conduct

Please note that the nodbi project is released with a
[Contributor Code of Conduct](CODE_OF_CONDUCT.md). By contributing to this
project you agree to abide by its terms.

### See rOpenSci [contributing guide](https://devguide.ropensci.org/contributingguide.html) for further details.

### Discussion forum

Check out our [discussion forum](https://discuss.ropensci.org) if you think your issue requires a longer form discussion.
<!-- If authentication is involved: do not share your username/password, or api keys/tokens in this issue - most likely the maintainer will have their own equivalent key -->

<!-- Do not share screen shots of code. Share actual code in text format. -->

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below, AND include a reproducible example (consider using a "reprex" https://cran.rstudio.com/web/packages/reprex/). If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
nodbi
=====

```{r echo=FALSE}
library("knitr")
hook_output <- knitr::knit_hooks$get("output")
knitr::knit_hooks$set(output = function(x, options) {
   lines <- options$output.lines
   if (is.null(lines)) {
     return(hook_output(x, options))  # pass to default hook
   }
   x <- unlist(strsplit(x, "\n"))
   more <- "..."
   if (length(lines)==1) {        # first n lines
     if (length(x) > lines) {
       # truncate the output, but add ....
       x <- c(head(x, lines), more)
     }
   } else {
     x <- c(if (abs(lines[1])>1) more else NULL,
            x[lines],
            if (length(x)>lines[abs(length(lines))]) more else NULL
           )
   }
   # paste these lines together
   x <- paste(c(x, ""), collapse = "\n")
   hook_output(x, options)
})

knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```


`nodbi` provides a single UI for interacting with many NoSQL databases. 

So far we support the following DBs:

* MongoDB
* Redis (server and serverless)
* CouchDB
* Elasticsearch
* etcd
* Riak

Currently we have support for data.frame's for the following operations

* Create - all DBs
* Get - all DBs
* Delete - all DBs
* Update - just CouchDB (others coming)

`sofa`, `mongolite`, `elastic`, and `etseed` are on CRAN

`RedisAPI`, `rrlite`, `reeack` are not on CRAN

## Install

```{r eval=FALSE}
install.packages("devtools")
devtools::install_github("ropensci/nodbi")
```

```{r}
library("nodbi")
```

## Initialize connections

Start CouchDB in your shell of choice by, e.g.: `couchdb`

```{r}
src_couchdb()
```

Start Elaticsearch in your shell of choice by, e.g.:

```sh
cd /usr/local/elasticsearch && bin/elasticsearch
```

```{r}
src_elasticsearch()
```

Start etcd in your shell of choice after installing etcd (https://github.com/coreos/etcd/releases/tag/v2.2.0) by, e.g.: `etcd`

```{r}
src_etcd()
```

Start MongoDB in your shell of choice by: `mongod`

```{r}
src_mongo()
```

If you want to use classic Redis server, we do that through the [RedisAPi][redisapi] 
package, and you'll need to start up Redis by e.g,. `redis-server` in your shell. 

```{r output.lines=1:10}
src_redis()
```

But if you want to use serverless Redis via [rlite][rlite], we do that through using 
with the [rrlite][rrlite] R package - and no need to start a server, of course.

```{r output.lines=1:10}
src_rlite()
```

Start your Riak server, then:

```{r}
src_riak()
```

## CouchDB

```{r}
src <- src_couchdb()
docout <- docdb_create(src, key = "mtcars", value = mtcars)
head( docdb_get(src, "mtcars") )
```

## etcd

```{r echo=FALSE}
src <- src_etcd()
invisible(docdb_delete(src, "/mtcars"))
```

```{r}
src <- src_etcd()
ff <- docdb_create(src, "/mtcars", mtcars)
head( docdb_get(src, "/mtcars") )
```

## Elasticsearch

Put the `iris` dataset into ES

```{r echo=FALSE}
src <- src_elasticsearch()
if (elastic::index_exists("iris")) invisible(docdb_delete(src, "iris"))
```

```{r eval=FALSE}
src <- src_elasticsearch()
ff <- docdb_create(src, "iris", iris)
head( docdb_get(src, "iris") )
#>   Sepal.Length Sepal.Width Petal.Length Petal.Width Species
#>          5.0         3.6          1.4         0.2  setosa
#>          4.9         3.1          1.5         0.1  setosa
#>          4.8         3.4          1.6         0.2  setosa
#>          5.4         3.9          1.3         0.4  setosa
#>          5.1         3.3          1.7         0.5  setosa
#>          5.2         3.4          1.4         0.2  setosa
```

## MongoDB

```{r output.lines=1:10}
library("ggplot2")
src <- src_mongo(verbose = FALSE)
ff <- docdb_create(src, "diamonds", diamonds)
docdb_get(src, "diamonds")
```

## Redis

```{r}
src <- src_rlite()
docdb_create(src, "mtcars", mtcars)
```

```{r output.lines=1:10}
docdb_get(src, "mtcars")
```

## Riak

```{r echo=FALSE}
src <- src_riak()
invisible(docdb_delete(src, "iris"))
```

```{r}
src <- src_riak()
docdb_create(src, "iris", iris)
```


## Use with dplyr

```{r}
library("dplyr")
src <- src_mongo(verbose = FALSE)
```

```{r}
docdb_get(src, "diamonds") %>%
  group_by(cut) %>%
  summarise(mean_depth = mean(depth), mean_price = mean(price))
```

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/nodbi/issues).
* License: MIT
* Get citation information for `nodbi` in R doing `citation(package = 'nodbi')`
* Please note that this project is released with a [Contributor Code of Conduct](CONDUCT.md). By participating in this project you agree to abide by its terms.

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[rlite]: https://github.com/seppo0010/rlite
[redisapi]: https://github.com/ropensci/RedisAPI
[rrlite]: https://github.com/ropensci/rrlite

```r
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

# Benchmarking nodbi


```r
library(nodbi)
library(microbenchmark)
```

make connections for use below


```r
src_m <- src_mongo()
src_red <- src_redis()
src_rl <- src_rlite()
src_c <- src_couchdb()
src_e <- src_elasticsearch()
src_ri <- src_riak()
#src_et <- src_etcd()
```

delete any datasets to be used below

to do ...


## initialize connection


```r
microbenchmark(
  mongo = src_mongo(),
  redis = src_redis(),
  rlite = src_rlite(),
  couch = src_couchdb(),
  elasticsearch = src_elasticsearch(),
  riak = src_riak(),
  etcd = src_etcd(),
  times = 100
)
#> Unit: milliseconds
#>           expr       min        lq      mean    median        uq       max
#>          mongo  1.082408  1.294793  2.040935  1.416415  1.705144  17.53605
#>          redis  2.497329  2.859878  4.767820  3.287540  4.492701  46.12102
#>          rlite  1.745421  2.008409  2.895608  2.190136  2.778588  18.95398
#>          couch 12.900579 14.179201 30.828498 16.013078 20.005475 321.11350
#>  elasticsearch 12.998076 14.447466 23.544338 17.785322 22.305072 203.46777
#>           riak 23.888188 26.996847 40.102676 29.710338 36.934067 185.78502
#>           etcd  5.672403  6.093102 11.164707  6.591032  9.857925 100.61295
#>  neval
#>    100
#>    100
#>    100
#>    100
#>    100
#>    100
#>    100
```

## create


```r
microbenchmark(
  mongo = docdb_create(src_m, paste0("iris", sample(1:1000, 1)), iris),
  redis = docdb_create(src_red, paste0("iris", sample(1:1000, 1)), iris),
  rlite = docdb_create(src_rl, paste0("iris", sample(1:1000, 1)), iris),
  couch = docdb_create(src_c, paste0("iris", sample(1:1000, 1)), iris),
  elasticsearch = docdb_create(src_e, paste0("iris", sample(1:1000, 1)), iris),
  riak = docdb_create(src_ri, paste0("iris", sample(1:1000, 1)), iris),
  #etcd = docdb_create(src_et, paste0("/iris", sample(1:1000, 1)), iris),
  times = 10
)
#> 
Complete! Processed total of 150 rows.
#> 
Complete! Processed total of 150 rows.
#> Error: length(attr(src, "dbs")) == 1 is not TRUE
```

## get

create some data that won't be affected by above


```r
try_del_create <- function(src, key) {
  invisible(tryCatch(docdb_delete(src, key), error = function(e) e))
  invisible(docdb_create(src, key, iris))
}
```


```r
try_del_create(src_m, "iris_get")
#> 
Complete! Processed total of 150 rows.
try_del_create(src_red, "iris_get")
try_del_create(src_rl, "iris_get")
try_del_create(src_c, "iris_get")
try_del_create(src_e, "iris_get")
#> 
  |                                                                       
  |                                                                 |   0%
  |                                                                       
  |=================================================================| 100%
try_del_create(src_ri, "iris_get")
#> Error: length(attr(src, "dbs")) == 1 is not TRUE
#try_del_create(src_et, "/iris_get")
```

benchmark


```r
microbenchmark(
  mongo = docdb_get(src_m, "iris_get"),
  redis = docdb_get(src_red, "iris_get"),
  rlite = docdb_get(src_rl, "iris_get"),
  couch = docdb_get(src_c, "iris_get"),
  elasticsearch = docdb_get(src_e, "iris_get"),
  riak = docdb_get(src_ri, "iris_get"),
  #etcd = docdb_get(src_et, "/iris_get"),
  times = 100
)
#> Error: Not Found (HTTP 404).
```

## delete

to do ...

## cleanup


```r
docdb_delete(src_m, "iris_get")
#> [1] TRUE
docdb_delete(src_red, "iris_get")
#> [1] 1
docdb_delete(src_rl, "iris_get")
#> [1] 1
docdb_delete(src_c, "iris_get")
#> $ok
#> [1] TRUE
docdb_delete(src_e, "iris_get")
#> $acknowledged
#> [1] TRUE
docdb_delete(src_ri, "iris_get")
#> Error: Not Found (HTTP 404).
#docdb_delete(src_et, "/iris_get")
```
## revdepcheck results

We checked 2 reverse dependencies, comparing R CMD check results across CRAN and dev versions of this package.

 * We saw 0 new problems
 * We failed to check 0 packages

# Platform

|field    |value                                                                       |
|:--------|:---------------------------------------------------------------------------|
|version  |R version 4.1.2 (2021-11-01)                                                |
|os       |macOS Big Sur 11.6.2                                                        |
|system   |x86_64, darwin17.0                                                          |
|ui       |RStudio                                                                     |
|language |(EN)                                                                        |
|collate  |de_DE.UTF-8                                                                 |
|ctype    |de_DE.UTF-8                                                                 |
|tz       |Europe/Berlin                                                               |
|date     |2021-12-27                                                                  |
|rstudio  |2021.09.0+351 Ghost Orchid (desktop)                                        |
|pandoc   |2.14.0.3 @ /Applications/RStudio.app/Contents/MacOS/pandoc/ (via rmarkdown) |

# Dependencies

|package |old   |new        |Δ  |
|:-------|:-----|:----------|:--|
|nodbi   |0.6.0 |0.6.0.9001 |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*