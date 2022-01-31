sofa
====



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/sofa)](https://cranchecks.info/pkgs/sofa)
[![R-check](https://github.com/ropensci/sofa/workflows/R-check/badge.svg)](https://github.com/ropensci/sofa/actions)
[![codecov.io](https://codecov.io/github/ropensci/sofa/coverage.svg?branch=master)](https://codecov.io/github/ropensci/sofa?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/sofa?color=ff69b4)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/sofa)](https://cran.r-project.org/package=sofa)

__An easy interface to CouchDB from R__

Note: Check out [*R4couchdb*](https://github.com/wactbprot/R4CouchDB), another R
package to interact with CouchDB.

sofa docs: https://docs.ropensci.org/sofa/

## CouchDB versions

`sofa` works with CouchDB v2 and v3. See [the builds](https://github.com/ropensci/sofa/actions?query=workflow%3AR-check) for checks on various CouchDB versions.

## CouchDB Info

* Docs: <http://docs.couchdb.org/en/latest/index.html>
* Installation: <http://docs.couchdb.org/en/latest/install/index.html>

## Connect to CouchDB

This may be starting it on your terminal/shell

```sh
couchdb
```

Or opening the CouchDB app on your machine, or running it in Docker. Whatever it
is, start it up.

## Install sofa

From CRAN


```r
install.packages("sofa")
```

Development version from GitHub


```r
remotes::install_github("ropensci/sofa")
```


```r
library('sofa')
```

## Cushions

Cushions? What? Since it's couch we gotta use `cushions` somehow. `cushions` are a
connection class containing all connection info to a CouchDB instance.
See `?Cushion` for help.

As an example, connecting to a Cloudant couch:


```r
z <- Cushion$new(
  host = "stuff.cloudant.com",
  transport = 'https',
  port = NULL,
  user = 'foobar',
  pwd = 'things'
)
```

Break down of parameters:

* `host`: the base url, without the transport (`http`/`https`)
* `path`: context path that is appended to the end of the url
* `transport`: `http` or `https`
* `port`: The port to connect to. Default: 5984. For Cloudant, have to set to `NULL`
* `user`: User name for the service.
* `pwd`: Password for the service, if any.
* `headers`: headers to pass in all requests

If you call `Cushion$new()` with no arguments you get a cushion set up for local
use on your machine, with all defaults used.


```r
x <- Cushion$new()
```

Ping the server


```r
x$ping()
```

Nice, it's working.

## More

See the docs https://docs.ropensci.org/sofa/ for more.


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/sofa/issues).
* License: MIT
* Get citation information for `sofa` in R doing `citation(package = 'sofa')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/sofa/blob/master/CODE_OF_CONDUCT.md
## query examples

docs: <https://docs.cloudant.com/cloudant_query.html>

field `smell` with values greater than 5

```
curl -v -XPOST -H 'Content-Type: application/json' 'http://localhost:5984/farts/_find' -d '{
  "selector": {
    "smell": {
      "$gt": 5
    }
  }
}'
```

field `smell` with values greater than 5, and return only `_id` and `smell` fields

```
curl -v -XPOST -H 'Content-Type: application/json' 'http://localhost:5984/farts/_find' -d '{
  "selector": {
    "smell": {
      "$gt": 5
    }
  },
  "fields": [
    "_id",
    "smell"
  ]
}'
```


curl -v -XGET -H 'Content-Type: application/json' 'http://localhost:5984/_stats'
sofa 0.4.0
==========

### NEW FEATURES

* new function `doc_upsert()`: updates an existing document or creates it if it doesn't yet exist (#69) work by @critichu

CouchDB v3 related changes

* made sure sofa works with v3; all examples/tests updated to use username/password  (#73)
* new function `db_bulk_get()` for the `/{db}/_bulk_get` route  (#73)
* fixed `design_search_many()`: in couch v2.2 and greater there's a new route `/{db}/_design/{ddoc}/_view/{view}/queries`, which is used in this fxn now instead of using the `/{db}/_design/{ddoc}/_view/{view}` route (#75)
* Cushion class gains new method `$version()` to get the CouchDB version you're using as a numeric (to enable progammatic couch version checking)
* `db_query()` changes: some new parameters added: `r`, `bookmark`, `update`, `stable`, `stale`, and `execution_stats` (#74)

### DEFUNCT

* `attach_get()` is now defunct, use `doc_attach_get()` (#76)

### MINOR IMPROVEMENTS

* added more tests (#61)
* `design_search()` now allows more possible values for start and end keys: `startkey_docid`, `start_key_doc_id`, `startkey`, `start_key`, `endkey_docid`, `end_key_doc_id`, `endkey`, `end_key` (#62)
* add title to vignettes (#71)
* for `docs_create()` internally support using user's setting for the R option `digits` to pass on to `jsonlite::toJSON` to control number of digits after decimal place (#66)

### BUG FIXES

* fixed authorization problems in `$ping()` method in Cushion; now separate `ping()` function calls `$ping()` method in Cushion (#72)

sofa 0.3.0
==========

### NEW FEATURES

* Gains new functions `db_index`, `db_index_create`, and `db_index_delete` 
for getting an index, creating one, and deleting one
* Gains function `design_search_many` to do many queries at once 
in a `POST` request (#56)
* `design_search` reworked to allow user to do a `GET` request or 
`POST` request depending on if they use `params` parameter or
`body` parameter - many parameters removed in the function 
definition, and are now to be passed to `params` or `body` (#56)
* `db_alldocs` gains new parameter `disk` to optionally
write data to disk instead of into the R session - should help
when data is very large (if disk is used fxn returns a file path) (#64)

### MINOR IMPROVEMENTS

* fix minor issues in vignette, and updated for working with CouchDB v2 and greater (#53) (#54) (#47)
* replace `httr` with `crul` for HTTP requests (#52)
* `design_copy` removed temporarily (#20) (#60)
* new issue and pull request template

### BUG FIXES

* Fix to docs for `design_search` (#57) thanks @michellymenezes
* Fix to `db_query` to make a single field passed to `fields` parameter
work (#63) thanks @gtumuluri
* Fix error in `doc_attach_get` (#58) thanks @gtumuluri


sofa 0.2.0
==========

### NEW FEATURES

* released to CRAN
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
(https://contributor-covenant.org), version 1.0.0, available at 
https://contributor-covenant.org/version/1/0/0/
## Test environments

* local macOS install, R 4.0.2
* ubuntu 16.04 (on travis-ci), R 4.0.2
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

There are no reverse dependencies.

---

This version adds some new functions and makes some bug fixes. This is a re-submission of the same version fixing a non-file package-anchored link in the documentation.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/sofa/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/sofa.git`
* Make sure to track progress upstream (i.e., on our version of `sofa` at `ropensci/sofa`) by doing `git remote add upstream https://github.com/ropensci/sofa.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/sofa`

### Also, check out our [discussion forum](https://discuss.ropensci.org)

### Prefer to Email? Get in touch: [myrmecocystus@gmail.com](mailto:myrmecocystus@gmail.com)

<!-- If this issue relates to usage of the package, whether a question, bug or similar, along with your query, please paste your devtools::session_info() or sessionInfo() into the code block below. If not, delete all this and proceed :) -->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
sofa
====

```{r, echo=FALSE}
knitr::opts_chunk$set(
  collapse=TRUE,
  comment="#>",
  warning=FALSE,
  message=FALSE
)
```

[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![cran checks](https://cranchecks.info/badges/worst/sofa)](https://cranchecks.info/pkgs/sofa)
[![R-check](https://github.com/ropensci/sofa/workflows/R-check/badge.svg)](https://github.com/ropensci/sofa/actions)
[![codecov.io](https://codecov.io/github/ropensci/sofa/coverage.svg?branch=master)](https://codecov.io/github/ropensci/sofa?branch=master)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/sofa?color=ff69b4)](https://github.com/metacran/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/sofa)](https://cran.r-project.org/package=sofa)

__An easy interface to CouchDB from R__

Note: Check out [*R4couchdb*](https://github.com/wactbprot/R4CouchDB), another R
package to interact with CouchDB.

sofa docs: https://docs.ropensci.org/sofa/

## CouchDB versions

`sofa` works with CouchDB v2 and v3. See [the builds](https://github.com/ropensci/sofa/actions?query=workflow%3AR-check) for checks on various CouchDB versions.

## CouchDB Info

* Docs: <http://docs.couchdb.org/en/latest/index.html>
* Installation: <http://docs.couchdb.org/en/latest/install/index.html>

## Connect to CouchDB

This may be starting it on your terminal/shell

```sh
couchdb
```

Or opening the CouchDB app on your machine, or running it in Docker. Whatever it
is, start it up.

## Install sofa

From CRAN

```{r eval=FALSE}
install.packages("sofa")
```

Development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/sofa")
```

```{r}
library('sofa')
```

## Cushions

Cushions? What? Since it's couch we gotta use `cushions` somehow. `cushions` are a
connection class containing all connection info to a CouchDB instance.
See `?Cushion` for help.

As an example, connecting to a Cloudant couch:

```{r eval=FALSE}
z <- Cushion$new(
  host = "stuff.cloudant.com",
  transport = 'https',
  port = NULL,
  user = 'foobar',
  pwd = 'things'
)
```

Break down of parameters:

* `host`: the base url, without the transport (`http`/`https`)
* `path`: context path that is appended to the end of the url
* `transport`: `http` or `https`
* `port`: The port to connect to. Default: 5984. For Cloudant, have to set to `NULL`
* `user`: User name for the service.
* `pwd`: Password for the service, if any.
* `headers`: headers to pass in all requests

If you call `Cushion$new()` with no arguments you get a cushion set up for local
use on your machine, with all defaults used.

```{r eval=FALSE}
x <- Cushion$new()
```

Ping the server

```{r eval=FALSE}
x$ping()
```

Nice, it's working.

## More

See the docs https://docs.ropensci.org/sofa/ for more.


## Meta

* Please [report any issues or bugs](https://github.com/ropensci/sofa/issues).
* License: MIT
* Get citation information for `sofa` in R doing `citation(package = 'sofa')`
* Please note that this project is released with a [Contributor Code of Conduct][coc]. By participating in this project you agree to abide by its terms.

[![ropensci_footer](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)

[coc]: https://github.com/ropensci/sofa/blob/master/CODE_OF_CONDUCT.md
---
title: "CouchDB Queries"
author: "Scott Chamberlain"
date: "`r Sys.Date()`"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{CouchDB Queries}
  \usepackage[utf8]{inputenc}
---

```{r echo=FALSE}
knitr::opts_chunk$set(
	comment = "#>",
	collapse = TRUE,
	warning = FALSE,
	message = FALSE
)
```

## CouchDB query table

| Operator Type | Operator | Argument             | Purpose                                                                                                                                                                                                                                                                                                                                                                                 |
| ------------- | -------- | -----------          | ------------------------------                                                                                                                                                                                                                                                                                                                                                          |
| (In)equality  | $lt      | Any JSON             | The field is less than the argument                                                                                                                                                                                                                                                                                                                                                     |
|               | $lte     | Any JSON             | The field is less than or equal to the argument                                                                                                                                                                                                                                                                                                                                         |
|               | $eq      | Any JSON             | The field is equal to the argument                                                                                                                                                                                                                                                                                                                                                      |
|               | $ne      | Any JSON             | The field is not equal to the argument                                                                                                                                                                                                                                                                                                                                                  |
|               | $gte     | Any JSON             | The field is greater than or equal to the argument                                                                                                                                                                                                                                                                                                                                      |
|               | $gt      | Any JSON             | The field is greater than the to the argument                                                                                                                                                                                                                                                                                                                                           |
| Object        | $exists  | Boolean              | Check whether the field exists or not, regardless of its value                                                                                                                                                                                                                                                                                                                          |
|               | $type    | String               | Check the document field’s type. Valid values are "null", "boolean", "number", "string", "array", and "object"                                                                                                                                                                                                                                                                          |
| Array         | $in      | Array of JSON values | The document field must exist in the list provided                                                                                                                                                                                                                                                                                                                                      |
|               | $nin     | Array of JSON values | The document field not must exist in the list provided                                                                                                                                                                                                                                                                                                                                  |
|               | $size    | Integer              | Special condition to match the length of an array field in a document. Non-array fields cannot match this condition                                                                                                                                                                                                                                                                     |
| Miscellaneous | $mod     | [Divisor, Remainder] | Divisor and Remainder are both positive or negative integers. Non-integer values result in a 404. Matches documents where field % Divisor == Remainder is true, and only when the document field is an integer                                                                                                                                                                          |
|               | $regex   | String               | A regular expression pattern to match against the document field. Only matches when the field is a string value and matches the supplied regular expression. The matching algorithms are based on the Perl Compatible Regular Expression (PCRE) library. For more information about what is implemented, see the see the [Erlang Regular Expression](http://erlang.org/doc/man/re.html) |
---
title: "sofa introduction"
author: "Scott Chamberlain"
date: "2020-06-25"
output:
  html_document:
    toc: true
    toc_float: true
    theme: readable
vignette: >
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{sofa introduction}
  %\VignetteEncoding{UTF-8}
---



Make sure your CouchDB installation is running.

## Install sofa

Stable version


```r
install.packages("sofa")
```

Development version


```r
remotes::install_github("ropensci/sofa")
```

Load library


```r
library(sofa)
```

## sofa package API

The following is a breakdown of the major groups of functions - note that not all are included.

__create a CouchDB client connection__

* `Cushion`

__work with databases__

* `db_alldocs`
* `db_changes`
* `db_compact`
* `db_create`
* `db_delete`
* `db_explain`
* `db_info`
* `db_list`
* `db_query`
* `db_replicate`
* `db_revisions`
* `db_index`
* `db_index_create`
* `db_index_delete`

__work with views/design documents__

* `design_create`
* `design_create_`
* `design_delete`
* `design_get`
* `design_head`
* `design_info`
* `design_search`
* `design_search_many`

__work with documents__

* `doc_create`
* `doc_delete`
* `doc_get`
* `doc_head`
* `doc_update`
* `db_bulk_create`
* `db_bulk_update`
* `doc_attach_create`
* `doc_attach_delete`
* `doc_attach_get`
* `doc_attach_info`

## Create a connection client

If your CouchDB instance requires username and password make sure to pass those to `Cushion$new`


```r
(x <- Cushion$new(user="admin", pwd="password"))
#> <sofa - cushion> 
#>   transport: http
#>   host: 127.0.0.1
#>   port: 5984
#>   path: 
#>   type: 
#>   user: admin
#>   pwd: <secret>
```

## Ping your server


```r
x$ping()
#> $couchdb
#> [1] "Welcome"
#> 
#> $version
#> [1] "3.1.0"
#> 
#> $git_sha
#> [1] "ff0feea20"
#> 
#> $uuid
#> [1] "30ed570659e8b72d688cfab563811c53"
#> 
#> $features
#> $features[[1]]
#> [1] "access-ready"
#> 
#> $features[[2]]
#> [1] "partitioned"
#> 
#> $features[[3]]
#> [1] "pluggable-storage-engines"
#> 
#> $features[[4]]
#> [1] "reshard"
#> 
#> $features[[5]]
#> [1] "scheduler"
#> 
#> 
#> $vendor
#> $vendor$name
#> [1] "The Apache Software Foundation"
```

## Create a new database




```r
db_create(x, 'cats')
#> $ok
#> [1] TRUE
```

## List databases


```r
db_list(x)
#> [1] "cats"
```

## Create a document


```r
doc1 <- '{"name": "leo", "color": "blue", "furry": true, "size": 1}'
doc_create(x, dbname = "cats", doc1, docid = "bluecat")
#> $ok
#> [1] TRUE
#> 
#> $id
#> [1] "bluecat"
#> 
#> $rev
#> [1] "1-41784f190c466d990684003a958c9f39"
```

and another!


```r
doc2 <- '{"name": "samson", "color": "red", "furry": false, "size": 3}'
doc_create(x, dbname = "cats", doc2)
#> $ok
#> [1] TRUE
#> 
#> $id
#> [1] "3e80ffb4c86ccbf35d2c3b3314000bfa"
#> 
#> $rev
#> [1] "1-08aef850a23f5ff95869c9cf5d9604dc"
```

and one more, cause 3's company


```r
doc3 <- '{"name": "matilda", "color": "green", "furry": false, "size": 5, "age": 2}'
doc_create(x, dbname = "cats", doc3)
#> $ok
#> [1] TRUE
#> 
#> $id
#> [1] "3e80ffb4c86ccbf35d2c3b3314001a3a"
#> 
#> $rev
#> [1] "1-953d3cfbbebb977fb75940c2bb0c93a1"
```

Note how we used a document id in the first document creation, but
not in the second and third. Using a document id is optional.

Also note that the third document has an additional field "age".

## Changes feed


```r
db_changes(x, "cats")
#> $results
#> $results[[1]]
#> $results[[1]]$seq
#> [1] "1-g1AAAABteJzLYWBgYMpgTmHgzcvPy09JdcjLz8gvLskBCScyJNX___8_K4M5kTEXKMBubJZqnpacjK4Yh_Y8FiDJ0ACk_oNMSWTIAgD59SI-"
#> 
#> $results[[1]]$id
#> [1] "3e80ffb4c86ccbf35d2c3b3314000bfa"
#> 
#> $results[[1]]$changes
#> $results[[1]]$changes[[1]]
#> $results[[1]]$changes[[1]]$rev
#> [1] "1-08aef850a23f5ff95869c9cf5d9604dc"
#> 
#> 
#> 
#> 
#> $results[[2]]
#> $results[[2]]$seq
#> [1] "2-g1AAAABteJzLYWBgYMpgTmHgzcvPy09JdcjLz8gvLskBCScyJNX___8_K4M5kSkXKMBubJZqnpacjK4Yh_Y8FiDJ0ACk_oNMSWTIAgD6OyI_"
#> 
#> $results[[2]]$id
#> [1] "3e80ffb4c86ccbf35d2c3b3314001a3a"
#> 
#> $results[[2]]$changes
#> $results[[2]]$changes[[1]]
#> $results[[2]]$changes[[1]]$rev
#> [1] "1-953d3cfbbebb977fb75940c2bb0c93a1"
#> 
#> 
#> 
#> 
#> $results[[3]]
#> $results[[3]]$seq
#> [1] "3-g1AAAACLeJzLYWBgYMpgTmHgzcvPy09JdcjLz8gvLskBCScyJNX___8_K4M5kSkXKMBubJZqnpacjK4Yh_Y8FiDJ0ACk_kNNYQSbkpZiZmGaaIauJwsAaKQq8g"
#> 
#> $results[[3]]$id
#> [1] "bluecat"
#> 
#> $results[[3]]$changes
#> $results[[3]]$changes[[1]]
#> $results[[3]]$changes[[1]]$rev
#> [1] "1-41784f190c466d990684003a958c9f39"
#> 
#> 
#> 
#> 
#> 
#> $last_seq
#> [1] "3-g1AAAACLeJzLYWBgYMpgTmHgzcvPy09JdcjLz8gvLskBCScyJNX___8_K4M5kSkXKMBubJZqnpacjK4Yh_Y8FiDJ0ACk_kNNYQSbkpZiZmGaaIauJwsAaKQq8g"
#> 
#> $pending
#> [1] 0
```

## Search

The simplest search just returns the documents.


```r
db_query(x, dbname = "cats", selector = list(`_id` = list(`$gt` = NULL)))$docs
#> [[1]]
#> [[1]]$`_id`
#> [1] "3e80ffb4c86ccbf35d2c3b3314000bfa"
#> 
#> [[1]]$`_rev`
#> [1] "1-08aef850a23f5ff95869c9cf5d9604dc"
#> 
#> [[1]]$name
#> [1] "samson"
#> 
#> [[1]]$color
#> [1] "red"
#> 
#> [[1]]$furry
#> [1] FALSE
#> 
#> [[1]]$size
#> [1] 3
#> 
#> 
#> [[2]]
#> [[2]]$`_id`
#> [1] "3e80ffb4c86ccbf35d2c3b3314001a3a"
#> 
#> [[2]]$`_rev`
#> [1] "1-953d3cfbbebb977fb75940c2bb0c93a1"
#> 
#> [[2]]$name
#> [1] "matilda"
#> 
#> [[2]]$color
#> [1] "green"
#> 
#> [[2]]$furry
#> [1] FALSE
#> 
#> [[2]]$size
#> [1] 5
#> 
#> [[2]]$age
#> [1] 2
#> 
#> 
#> [[3]]
#> [[3]]$`_id`
#> [1] "bluecat"
#> 
#> [[3]]$`_rev`
#> [1] "1-41784f190c466d990684003a958c9f39"
#> 
#> [[3]]$name
#> [1] "leo"
#> 
#> [[3]]$color
#> [1] "blue"
#> 
#> [[3]]$furry
#> [1] TRUE
#> 
#> [[3]]$size
#> [1] 1
```

Search for cats that are red


```r
db_query(x, dbname = "cats", selector = list(color = "red"))$docs
#> [[1]]
#> [[1]]$`_id`
#> [1] "3e80ffb4c86ccbf35d2c3b3314000bfa"
#> 
#> [[1]]$`_rev`
#> [1] "1-08aef850a23f5ff95869c9cf5d9604dc"
#> 
#> [[1]]$name
#> [1] "samson"
#> 
#> [[1]]$color
#> [1] "red"
#> 
#> [[1]]$furry
#> [1] FALSE
#> 
#> [[1]]$size
#> [1] 3
```

Search for cats that are furry


```r
db_query(x, dbname = "cats", selector = list(furry = TRUE))$docs
#> [[1]]
#> [[1]]$`_id`
#> [1] "bluecat"
#> 
#> [[1]]$`_rev`
#> [1] "1-41784f190c466d990684003a958c9f39"
#> 
#> [[1]]$name
#> [1] "leo"
#> 
#> [[1]]$color
#> [1] "blue"
#> 
#> [[1]]$furry
#> [1] TRUE
#> 
#> [[1]]$size
#> [1] 1
```

Return only certain fields


```r
db_query(x, dbname = "cats",
         selector = list(size = list(`$gt` = 2)),
         fields = c("name", "color"))$docs
#> [[1]]
#> [[1]]$name
#> [1] "samson"
#> 
#> [[1]]$color
#> [1] "red"
#> 
#> 
#> [[2]]
#> [[2]]$name
#> [1] "matilda"
#> 
#> [[2]]$color
#> [1] "green"
```

Convert the result of a query into a data.frame using `jsonlite`


```r
library('jsonlite')
res <- db_query(x, dbname = "cats",
                 selector = list(`_id` = list(`$gt` = NULL)),
                 fields = c("name", "color", "furry", "size", "age"),
                 as = "json")

fromJSON(res)$docs
#>      name color furry size age
#> 1  samson   red FALSE    3  NA
#> 2 matilda green FALSE    5   2
#> 3     leo  blue  TRUE    1  NA
```
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_alldocs.R
\name{db_alldocs}
\alias{db_alldocs}
\title{List all docs in a given database.}
\usage{
db_alldocs(
  cushion,
  dbname,
  descending = NULL,
  startkey = NULL,
  endkey = NULL,
  limit = NULL,
  include_docs = FALSE,
  as = "list",
  disk = NULL,
  ...
)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name. (character)}

\item{descending}{Return in descending order? (logical)}

\item{startkey}{Document ID to start at. (character)}

\item{endkey}{Document ID to end at. (character)}

\item{limit}{Number document IDs to return. (numeric)}

\item{include_docs}{(logical) If \code{TRUE}, returns docs themselves,
in addition to IDs. Default: \code{FALSE}}

\item{as}{(character) One of list (default) or json}

\item{disk}{write to disk or not. By default, data is in the R session;
if you give a file path, we'll write data to disk and you'll get back
the file path. by default, we save in the R session}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
List all docs in a given database.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("leothelion" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="leothelion"))
}
db_create(x, dbname='leothelion')
db_bulk_create(x, mtcars, dbname="leothelion")

db_alldocs(x, dbname="leothelion")
db_alldocs(x, dbname="leothelion", as='json')
db_alldocs(x, dbname="leothelion", limit=2)
db_alldocs(x, dbname="leothelion", limit=2, include_docs=TRUE)

# curl options
res <- db_alldocs(x, dbname="leothelion", verbose = TRUE)

# write data to disk - useful when data is very large
## create omdb dataset first
file <- system.file("examples/omdb.json", package = "sofa")
strs <- readLines(file)
if ("omdb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="omdb"))
}
db_create(x, dbname='omdb')
invisible(db_bulk_create(x, "omdb", strs))

## get all docs, writing them to disk
res <- db_alldocs(x, dbname="omdb", disk = (f <- tempfile(fileext=".json")))
res
readLines(res, n = 10)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_info.r
\name{db_info}
\alias{db_info}
\title{List database info.}
\usage{
db_info(cushion, dbname, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
List database info.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("sofadb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="sofadb"))
}
db_create(x, dbname='sofadb')

db_info(x, dbname="sofadb")
db_info(x, dbname="sofadb", as='json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_bulk_create.R
\name{db_bulk_create}
\alias{db_bulk_create}
\title{Create documents via the bulk API}
\usage{
db_bulk_create(
  cushion,
  dbname,
  doc,
  docid = NULL,
  how = "rows",
  as = "list",
  ...
)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. Required.}

\item{doc}{A data.frame, list, or JSON as a character string. Required.}

\item{docid}{Document IDs, ignored for now, eventually, you can pass in a
list, or vector to be the ids for each document created. Has to be the same
length as the number of documents.}

\item{how}{(character) One of rows (default) or columns. If rows, each row
becomes a separate document; if columns, each column becomes a separate
document.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
Either a list or json (depending on \code{as} parameter), with
each element an array of key:value pairs:
\itemize{
\item ok - whether creation was successful
\item id - the document id
\item rev - the revision id
}
}
\description{
Create documents via the bulk API
}
\details{
Note that row.names are dropped from data.frame inputs.
}
\examples{
\dontrun{
# initialize a CouchDB connection
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

# From a data.frame
if ("bulktest" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="bulktest"))
}
db_create(x, dbname="bulktest")
db_bulk_create(x, "bulktest", mtcars)

if ("bulktest2" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="bulktest2"))
}
db_create(x, dbname="bulktest2")
db_bulk_create(x, "bulktest2", iris)

# data.frame with 1 or more columns as neseted lists
mtcars$stuff <- list("hello_world")
mtcars$stuff2 <- list("hello_world","things")
if ("bulktest3" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="bulktest3"))
}
db_create(x, dbname="bulktest3")
db_bulk_create(x, "bulktest3", mtcars)

# From a json character string, or more likely, many json character strings
library("jsonlite")
strs <- as.character(parse_df(mtcars, "columns"))
if ("bulkfromchr" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="bulkfromchr"))
}
db_create(x, dbname="bulkfromchr")
db_bulk_create(x, "bulkfromchr", strs)

# From a list of lists
library("jsonlite")
lst <- parse_df(mtcars, tojson=FALSE)
if ("bulkfromchr" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="bulkfromchr"))
}
db_create(x, dbname="bulkfromchr")
db_bulk_create(x, "bulkfromchr", lst)

# iris dataset - by rows
if ("irisrows" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="irisrows"))
}
db_create(x, dbname="irisrows")
db_bulk_create(x, "irisrows", apply(iris, 1, as.list))

# iris dataset - by columns - doesn't quite work yet
# if ("iriscolumns" \%in\% db_list(x)) {
#   invisible(db_delete(x, dbname="iriscolumns"))
# }
# db_create(x, dbname="iriscolumns")
# db_bulk_create(x, "iriscolumns", parse_df(iris, "columns", tojson=FALSE), how="columns")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/session.R
\name{session}
\alias{session}
\title{session}
\usage{
session(cushion, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
session
}
\examples{
\dontrun{
# Create a CouchDB connection client
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

session(x)
session(x, as = 'json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/membership.R
\name{membership}
\alias{membership}
\title{membership}
\usage{
membership(cushion, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
membership
}
\examples{
\dontrun{
# Create a CouchDB connection client
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

membership(x)
membership(x, as = 'json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_bulk_get.R
\name{db_bulk_get}
\alias{db_bulk_get}
\title{Query many documents at once}
\usage{
db_bulk_get(cushion, dbname, docid_rev, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. Required.}

\item{docid_rev}{A list of named lists, each of which must have the
slot \code{id}, and optionally \code{rev} for the revision of the document id}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
Either a list or json (depending on \code{as} parameter)
}
\description{
Query many documents at once
}
\examples{
\dontrun{
# initialize a CouchDB connection
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("bulkgettest" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="bulkgettest"))
}
db_create(x, dbname="bulkgettest")
db_bulk_create(x, "bulkgettest", mtcars)
res <- db_query(x, dbname = "bulkgettest", selector = list(cyl = 8))

# with ids only
ids <- vapply(res$docs, "[[", "", "_id")
ids_only <- lapply(ids[1:5], function(w) list(id = w))
db_bulk_get(x, "bulkgettest", docid_rev = ids_only)

# with ids and revs
ids_rev <- lapply(res$docs[1:3],
  function(w) list(id = w$`_id`, rev = w$`_rev`))
db_bulk_get(x, "bulkgettest", docid_rev = ids_rev)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/getattach.r
\name{attach_get}
\alias{attach_get}
\title{Get an attachment}
\usage{
attach_get(...)
}
\arguments{
\item{...}{ignored}
}
\description{
This function is defunct. See \code{\link[=doc_attach_get]{doc_attach_get()}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc_delete.r
\name{doc_delete}
\alias{doc_delete}
\title{Delete a document in a database.}
\usage{
doc_delete(cushion, dbname, docid, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name. (character)}

\item{docid}{Document ID (character)}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Delete a document in a database.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

# create a database
if ("sofadb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="sofadb"))
}
db_create(x, dbname='sofadb')

doc3 <- "<top><a/><b/><c><d/><e>bob</e></c></top>"
doc_create(x, dbname="sofadb", doc=doc3, docid="newnewxml")
doc_delete(x, dbname="sofadb", docid="newnewxml")

# wrong docid name
doc_create(x, dbname="sofadb", doc=doc3, docid="newxml")
# doc_delete(x, dbname="sofadb", docid="wrongname")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/attach.r
\name{attachments}
\alias{attachments}
\alias{doc_attach_create}
\alias{doc_attach_info}
\alias{doc_attach_get}
\alias{doc_attach_delete}
\title{Work with attachments}
\usage{
doc_attach_create(
  cushion,
  dbname,
  docid,
  attachment,
  attname,
  as = "list",
  ...
)

doc_attach_info(cushion, dbname, docid, attname, ...)

doc_attach_get(cushion, dbname, docid, attname = NULL, type = "raw", ...)

doc_attach_delete(cushion, dbname, docid, attname, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. Required.}

\item{docid}{(character) Document ID. Required.}

\item{attachment}{(character) A file name. Required.}

\item{attname}{(character) Attachment name. Required.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}

\item{type}{(character) one of raw (default) or text. required.}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Work with attachments
}
\details{
Methods:
\itemize{
\item \code{doc_attach_create} - create an attachment
\item \code{doc_attach_info} - get info (headers) for an attachment
\item \code{doc_attach_get} - get an attachment. this method does not attempt
to read the object into R, but only gets the raw bytes or plain
text. See examples for how to read some attachment types
\item \code{doc_attach_delete} - delete and attachment
}
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("foodb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="foodb"))
}
db_create(x, dbname='foodb')

# create an attachment on an existing document
## create a document first
doc <- '{"name":"stuff", "drink":"soda"}'
doc_create(x, dbname="foodb", doc=doc, docid="asoda")

## create a csv attachment
row.names(mtcars) <- NULL
file <- tempfile(fileext = ".csv")
write.csv(mtcars, file = file, row.names = FALSE)
doc_attach_create(x, dbname="foodb", docid="asoda",
  attachment=file, attname="mtcarstable.csv")

## create a binary (png) attachment
file <- tempfile(fileext = ".png")
png(file)
plot(1:10)
dev.off()
doc_attach_create(x, dbname="foodb", docid="asoda",
  attachment=file, attname="img.png")

## create a binary (pdf) attachment
file <- tempfile(fileext = ".pdf")
pdf(file)
plot(1:10)
dev.off()
doc_attach_create(x, dbname="foodb", docid="asoda",
  attachment=file, attname="plot.pdf")

# get info for an attachment (HEAD request)
doc_attach_info(x, "foodb", docid="asoda", attname="mtcarstable.csv")
doc_attach_info(x, "foodb", docid="asoda", attname="img.png")
doc_attach_info(x, "foodb", docid="asoda", attname="plot.pdf")

# get an attachment (GET request)
res <- doc_attach_get(x, "foodb", docid="asoda",
  attname="mtcarstable.csv", type = "text")
read.csv(text = res)
doc_attach_get(x, "foodb", docid="asoda", attname="img.png")
doc_attach_get(x, "foodb", docid="asoda", attname="plot.pdf")
## OR, don't specify an attachment and list the attachments
(attchms <- doc_attach_get(x, "foodb", docid="asoda", type="text"))
jsonlite::fromJSON(attchms)

# delete an attachment
doc_attach_delete(x, "foodb", docid="asoda", attname="mtcarstable.csv")
doc_attach_delete(x, "foodb", docid="asoda", attname="img.png")
doc_attach_delete(x, "foodb", docid="asoda", attname="plot.pdf")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_compact.R
\name{db_compact}
\alias{db_compact}
\title{Request compaction of the specified database}
\usage{
db_compact(cushion, dbname, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name. Required.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Request compaction of the specified database
}
\details{
Compaction compresses the disk database file by performing the following
operations:
\itemize{
\item Writes a new, optimised, version of the database file, removing any unused
sections from the new version during write. Because a new file is temporarily
created for this purpose, you may require up to twice the current storage space
of the specified database in order for the compaction routine to complete.
\item Removes old revisions of documents from the database, up to the per-database
limit specified by the _revs_limit database parameter.
}

Compaction can only be requested on an individual database; you cannot compact all
the databases for a CouchDB instance. The compaction process runs as a background
process. You can determine if the compaction process is operating on a database by
obtaining the database meta information, the compact_running value of the returned
database structure will be set to true. See GET /{db}. You can also obtain a list of
running processes to determine whether compaction is currently running.
See "/_active_tasks"
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))
# db_compact(x, dbname = "iris")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_changes.R
\name{db_changes}
\alias{db_changes}
\title{List changes to a database.}
\usage{
db_changes(
  cushion,
  dbname,
  descending = NULL,
  startkey = NULL,
  endkey = NULL,
  since = NULL,
  limit = NULL,
  include_docs = NULL,
  feed = "normal",
  heartbeat = NULL,
  filter = NULL,
  as = "list",
  ...
)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name. (character)}

\item{descending}{Return in descending order? (logical)}

\item{startkey}{Document ID to start at. (character)}

\item{endkey}{Document ID to end at. (character)}

\item{since}{Start the results from the change immediately after the given
sequence number.}

\item{limit}{Number document IDs to return. (numeric)}

\item{include_docs}{(character) If "true", returns docs themselves, in
addition to IDs}

\item{feed}{Select the type of feed. One of normal, longpoll, or continuous.
See description. (character)}

\item{heartbeat}{Period in milliseconds after which an empty line is sent in
the results. Only applicable for longpoll or continuous feeds. Overrides any
timeout to keep the feed alive indefinitely. (numeric (milliseconds))}

\item{filter}{Reference a filter function from a design document to
selectively get updates.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
Either a list of json (depending on \code{as} parameter), with
keys:
\itemize{
\item results - Changes made to a database, length 0 if no changes.
Each of these has:
\itemize{
\item changes - List of document`s leafs with single field rev
\item id - Document ID
\item seq - Update sequence
}
\item last_seq - Last change update sequence
\item pending - Count of remaining items in the feed
}
}
\description{
Of course it doesn't make much sense to use certain options in
_changes. For example, using feed=longpoll or continuous doesn't make much
sense within R itself.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("leoalion" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="leoalion"))
}
db_create(x, dbname='leoalion')

# no changes
res <- db_changes(x, dbname="leoalion")
res$results

# create a document
doc1 <- '{"name": "drink", "type": "water", "score": 5}'
doc_create(x, dbname="leoalion", doc1, docid="awater")

# now there's changes
res <- db_changes(x, dbname="leoalion")
res$results

# as JSON
db_changes(x, dbname="leoalion", as='json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/restart.R
\name{restart}
\alias{restart}
\title{Restart your Couchdb instance}
\usage{
restart(cushion = "localhost", as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Restart your Couchdb instance
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

# restart(x)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ping.r
\name{ping}
\alias{ping}
\title{Ping a CouchDB server}
\usage{
ping(cushion, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Ping a CouchDB server
}
\examples{
\dontrun{
# initialize a CouchDB connection
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

# call ping on the cushion object, or pass the cushion to ping()
x$ping()
ping(x)
ping(x, as = "json")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/uuids.r
\name{uuids}
\alias{uuids}
\title{Get uuids.}
\usage{
uuids(cushion, count = 1, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{count}{(numeric) Number of uuids to return. Default: 1}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Get uuids.
}
\examples{
\dontrun{
# Create a CouchDB connection client
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

uuids(x)
uuids(x, as = 'json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cushion.R
\name{Cushion}
\alias{Cushion}
\title{sofa connection client}
\value{
An object of class \code{Cushion}, with variables accessible for
host, port, path, transport, user, pwd, and headers. Functions are callable
to get headers, and to make the base url sent with all requests.
}
\description{
sofa connection client

sofa connection client
}
\section{CouchDB versions}{

\pkg{sofa} was built assuming CouchDB version 2 or greater. Some
functionality of this package will work with versions < 2, while
some may not (mango queries, see \code{\link[=db_query]{db_query()}}). I don't
plan to support older CouchDB versions per se.
}

\examples{
\dontrun{
# Create a CouchDB connection client
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

## metadata
x$host
x$path
x$port
x$type

## ping the CouchDB server
x$ping()

## get CouchDB version
x$version()

# create database
if (!"stuff" \%in\% db_list(x)) {
  db_create(x, "stuff")
}

# add documents to a database
if (!"sofadb" \%in\% db_list(x)) {
  db_create(x, "sofadb")
}
doc1 <- '{"name": "drink", "beer": "IPA", "score": 5}'
doc_create(x, dbname="sofadb", docid="abeer", doc1)

# bulk create
if (!"mymtcars" \%in\% db_list(x)) {
  db_create(x, "mymtcars")
}
db_bulk_create(x, dbname="mymtcars", doc = mtcars)
db_list(x)

## database info
db_info(x, "mymtcars")

## list dbs
db_list(x)

## all docs
db_alldocs(x, "mymtcars", limit = 3)

## changes
db_changes(x, "mymtcars")

# With auth
# x <- Cushion$new(user = 'sckott', pwd = 'sckott')

# Using Cloudant
# z <- Cushion$new(host = "ropensci.cloudant.com", transport = 'https', port = NULL,
#   user = 'ropensci', pwd = Sys.getenv('CLOUDANT_PWD'))
# z
# db_list(z)
# db_create(z, "stuff2")
# db_info(z, "stuff2")
# db_alldocs(z, "foobar")
}
}
\section{Public fields}{
\if{html}{\out{<div class="r6-fields">}}
\describe{
\item{\code{host}}{(character) host}

\item{\code{port}}{(integer) port}

\item{\code{path}}{(character) url path, if any}

\item{\code{transport}}{(character) transport schema, (http|https)}

\item{\code{user}}{(character) username}

\item{\code{pwd}}{(character) password}

\item{\code{headers}}{(list) named list of headers}
}
\if{html}{\out{</div>}}
}
\section{Methods}{
\subsection{Public methods}{
\itemize{
\item \href{#method-new}{\code{Cushion$new()}}
\item \href{#method-print}{\code{Cushion$print()}}
\item \href{#method-ping}{\code{Cushion$ping()}}
\item \href{#method-make_url}{\code{Cushion$make_url()}}
\item \href{#method-get_headers}{\code{Cushion$get_headers()}}
\item \href{#method-get_auth}{\code{Cushion$get_auth()}}
\item \href{#method-version}{\code{Cushion$version()}}
\item \href{#method-clone}{\code{Cushion$clone()}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-new"></a>}}
\if{latex}{\out{\hypertarget{method-new}{}}}
\subsection{Method \code{new()}}{
Create a new \code{Cushion} object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cushion$new(host, port, path, transport, user, pwd, headers)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{host}}{(character) A base URL (without the transport), e.g.,
\code{localhost}, \verb{127.0.0.1}, or \code{foobar.cloudant.com}}

\item{\code{port}}{(numeric) Port. Remember that if you don't want a port set,
set this parameter to \code{NULL}. Default: \code{5984}}

\item{\code{path}}{(character) context path that is appended to the end of
the url. e.g., \code{bar} in \verb{http://foo.com/bar}. Default: \code{NULL},
ignored}

\item{\code{transport}}{(character) http or https. Default: http}

\item{\code{user, pwd}}{(character) user name, and password. these are used in all
requests. if absent, they are not passed to requests}

\item{\code{headers}}{A named list of headers. These headers are used in all
requests. To use headers in individual requests and not others, pass
in headers via \code{...} in a function call.}
}
\if{html}{\out{</div>}}
}
\subsection{Returns}{
A new \code{Cushion} object
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-print"></a>}}
\if{latex}{\out{\hypertarget{method-print}{}}}
\subsection{Method \code{print()}}{
print method for \code{Cushion}
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cushion$print()}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{x}}{self}

\item{\code{...}}{ignored}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-ping"></a>}}
\if{latex}{\out{\hypertarget{method-ping}{}}}
\subsection{Method \code{ping()}}{
Ping the CouchDB server
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cushion$ping(as = "list", ...)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{as}}{(character) One of list (default) or json}

\item{\code{...}}{curl options passed to \link[crul:verb-GET]{crul::verb-GET}}
}
\if{html}{\out{</div>}}
}
}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-make_url"></a>}}
\if{latex}{\out{\hypertarget{method-make_url}{}}}
\subsection{Method \code{make_url()}}{
Construct full base URL from the pieces in the
connection object
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cushion$make_url()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_headers"></a>}}
\if{latex}{\out{\hypertarget{method-get_headers}{}}}
\subsection{Method \code{get_headers()}}{
Get list of headers that will be sent with
each request
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cushion$get_headers()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-get_auth"></a>}}
\if{latex}{\out{\hypertarget{method-get_auth}{}}}
\subsection{Method \code{get_auth()}}{
Get list of auth values, user and pwd
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cushion$get_auth()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-version"></a>}}
\if{latex}{\out{\hypertarget{method-version}{}}}
\subsection{Method \code{version()}}{
Get the CouchDB version as a numeric
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cushion$version()}\if{html}{\out{</div>}}
}

}
\if{html}{\out{<hr>}}
\if{html}{\out{<a id="method-clone"></a>}}
\if{latex}{\out{\hypertarget{method-clone}{}}}
\subsection{Method \code{clone()}}{
The objects of this class are cloneable with this method.
\subsection{Usage}{
\if{html}{\out{<div class="r">}}\preformatted{Cushion$clone(deep = FALSE)}\if{html}{\out{</div>}}
}

\subsection{Arguments}{
\if{html}{\out{<div class="arguments">}}
\describe{
\item{\code{deep}}{Whether to make a deep clone.}
}
\if{html}{\out{</div>}}
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc_head.r
\name{doc_head}
\alias{doc_head}
\title{Get header info for a document}
\usage{
doc_head(cushion, dbname, docid, ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. Required.}

\item{docid}{(character) Document ID. Required.}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Get header info for a document
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

# create a database
if ("sofadb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="sofadb"))
}
db_create(x, dbname='sofadb')

# create a document
doc1 <- '{"name": "drink", "beer": "IPA", "score": 5}'
doc_create(x, dbname="sofadb", doc1, docid="abeer")

# run doc_head
doc_head(x, dbname="sofadb", docid="abeer")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_index.R
\name{db_index}
\alias{db_index}
\alias{db_index_create}
\alias{db_index_delete}
\title{Create and get database indexes}
\usage{
db_index(cushion, dbname, as = "list", ...)

db_index_create(cushion, dbname, body, as = "list", ...)

db_index_delete(cushion, dbname, design, index_name, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name, required}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}

\item{body}{(named list) index fields, required}

\item{design}{(character) Design document name}

\item{index_name}{(character) index name}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Create and get database indexes
}
\section{Body parameters}{

\itemize{
\item index (json) - JSON object describing the index to create.
\item ddoc (string) - Name of the design document in which the index will be
created. By default, each index will be created in its own design
document. Indexes can be grouped into design documents for efficiency.
However, a change to one index in a design document will invalidate all
other indexes in the same document (similar to views). Optional
\item name (string) - Name of the index. If no name is provided, a name will
be generated automatically. Optional
\item type (string) - Can be "json" or "text". Defaults to json. Geospatial
indexes will be supported in the future. Optional Text indexes are
supported via a third party library Optional
\item partial_filter_selector (json) - A selector to apply to documents at
indexing time, creating a partial index. Optional
}
}

\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

# create a database first
if ("testing" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="testing"))
}
db_create(x, "testing")

# get indexes
db_index(x, "testing")

# create indexes
body <- list(index = list(fields = I("foo")), name = "foo-index", type = "json")
db_index_create(x, "testing", body = body)

# get indexes, after creating another index
db_index(x, "testing")

# delete an index
res <- db_index(x, "testing")
db_index_delete(x, "testing", res$indexes[[2]]$ddoc, res$indexes[[2]]$name)
## and it's gone
db_index(x, "testing")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc_update.r
\name{doc_update}
\alias{doc_update}
\title{Update a document.}
\usage{
doc_update(cushion, dbname, doc, docid, rev, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. Required.}

\item{doc}{(character) Document content. Required.}

\item{docid}{(character) Document ID. Required.}

\item{rev}{(character) Revision id. Required.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Update a document.
}
\details{
Internally, this function adds in the docid and revision id,
required to do a document update
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("sofadb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="sofadb"))
}
db_create(x, dbname='sofadb')

doc1 <- '{"name":"drink","beer":"IPA"}'
doc_create(x, dbname="sofadb", doc=doc1, docid="b_beer")
doc_get(x, dbname = "sofadb", docid = "b_beer")
revs <- db_revisions(x, dbname = "sofadb", docid = "b_beer")
doc2 <- '{"name":"drink","beer":"IPA","note":"yummy","note2":"yay"}'
doc_update(x, dbname="sofadb", doc=doc2, docid="b_beer", rev=revs[1])
db_revisions(x, dbname = "sofadb", docid = "b_beer")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc_upsert.R
\name{doc_upsert}
\alias{doc_upsert}
\title{Create a new document or update an existing one}
\usage{
doc_upsert(cushion, dbname, doc, docid)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. Required.}

\item{doc}{(character) Document content. Required.}

\item{docid}{(character) Document ID. Required.}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Create a new document or update an existing one
}
\details{
Internally, this function attempts to update a document with the given name. \cr
If the document does not exist, it is created
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("sofadb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="sofadb"))
}
db_create(x, 'sofadb')

# create a document
doc1 <- '{"name": "drink", "beer": "IPA", "score": 5}'
doc_upsert(x, dbname="sofadb", doc1, docid="abeer")

#update the document
doc2 <- '{"name": "drink", "beer": "lager", "score": 6}'
doc_upsert(x, dbname="sofadb", doc2, docid="abeer")


doc_get(x, dbname = "sofadb", docid = "abeer")
}
}
\author{
George Kritikos
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/active_tasks.R
\name{active_tasks}
\alias{active_tasks}
\title{active tasks}
\usage{
active_tasks(cushion, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
active tasks
}
\examples{
\dontrun{
# Create a CouchDB connection client
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

active_tasks(x)
active_tasks(x, as = 'json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_delete.r
\name{db_delete}
\alias{db_delete}
\title{Delete a database.}
\usage{
db_delete(cushion, dbname, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Delete a database.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

# local databasees
## create database first, then delete
db_create(x, dbname='newdb')
db_delete(x, dbname='newdb')

## with curl info while doing request
library('crul')
db_create(x, 'newdb')
db_delete(x, 'newdb', verbose = TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/revisions.r
\name{db_revisions}
\alias{db_revisions}
\title{Get document revisions.}
\usage{
db_revisions(cushion, dbname, docid, simplify = TRUE, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name}

\item{docid}{Document ID}

\item{simplify}{(logical) Simplify to character vector of revision ids.
If \code{FALSE}, gives back availability info too. Default: \code{TRUE}}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Get document revisions.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("sofadb" \%in\% db_list(x)) {
 db_delete(x, dbname = "sofadb")
}
db_create(x, dbname = "sofadb")

doc1 <- '{"name": "drink", "beer": "IPA", "score": 5}'
doc_create(x, dbname="sofadb", doc1, docid="abeer")
doc_create(x, dbname="sofadb", doc1, docid="morebeer", as='json')

db_revisions(x, dbname="sofadb", docid="abeer")
db_revisions(x, dbname="sofadb", docid="abeer", simplify=FALSE)
db_revisions(x, dbname="sofadb", docid="abeer", as='json')
db_revisions(x, dbname="sofadb", docid="abeer", simplify=FALSE, as='json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_replicate.r
\name{db_replicate}
\alias{db_replicate}
\title{Upload (replicate) a local database to a remote database server,
e.g., Cloudant, Iriscouch}
\usage{
db_replicate(from, to, dbname, createdb = FALSE, as = "list", ...)
}
\arguments{
\item{from}{Couch to replicate from. An object of class \link{Cushion}.
Required.}

\item{to}{Remote couch to replicate to. An object of class \link{Cushion}.
Required.}

\item{dbname}{(character) Database name. Required.}

\item{createdb}{If \code{TRUE}, the function creates the db on the remote
server before uploading. The db has to exist before uploading, so either
you do it separately or this fxn can do it for you. Default: \code{FALSE}}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Upload (replicate) a local database to a remote database server,
e.g., Cloudant, Iriscouch
}
\examples{
\dontrun{
if (interactive()) {
## create a connection
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

# Create a database locally
db_list(x)
if ("hello_earth" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="hello_earth"))
}
db_create(x, 'hello_earth')

## replicate to a remote server
z <- Cushion$new(host = "ropensci.cloudant.com", transport = 'https',
  port = NULL, user = 'ropensci', pwd = Sys.getenv('CLOUDANT_PWD'))

## do the replication
db_replicate(x, z, dbname = "hello_earth")

## check changes on the remote
db_list(z)
db_changes(z, dbname = "hello_earth")

## make some changes on the remote
doc_create(z, dbname = "hello_earth",
  '{"language":"python","library":"requests"}', 'stuff')
changes(z, dbname = "hello_earth")

## create another document, and try to get it
doc_create(z, dbname = "hello_earth", doc = '{"language":"R"}',
  docid="R_rules")
doc_get(z, dbname = "hello_earth", docid='R_rules')

## cleanup - delete the database
db_delete(z, 'hello_earth')
}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc_get.r
\name{doc_get}
\alias{doc_get}
\title{Get a document from a database.}
\usage{
doc_get(
  cushion,
  dbname,
  docid,
  rev = NULL,
  attachments = FALSE,
  deleted = FALSE,
  revs = FALSE,
  revs_info = FALSE,
  conflicts = FALSE,
  deleted_conflicts = FALSE,
  local_seq = FALSE,
  as = "list",
  ...
)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name}

\item{docid}{Document ID}

\item{rev}{Revision id of the document to get. If NULL, gets current revision}

\item{attachments}{(logical) Whether to include _attachments field.}

\item{deleted}{(logical) Whether to include _deleted field.}

\item{revs}{(logical) Whether to include _revisions field.}

\item{revs_info}{(logical) Whether to include _revs_info field.}

\item{conflicts}{(logical) Whether to include _conflicts field.}

\item{deleted_conflicts}{(logical) Whether to include _deleted_conflicts field.}

\item{local_seq}{(logical) Whether to include _local_seq field.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Get a document from a database.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("sofadb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="sofadb"))
}
db_create(x, dbname="sofadb")

# create a document
doc1 <- '{"name": "drink", "type": "drink", "score": 5}'
doc_create(x, dbname="sofadb", doc1, docid="asoda")

doc_get(x, dbname="sofadb", docid="asoda")
revs <- db_revisions(x, dbname="sofadb", docid="asoda")
doc_get(x, dbname="sofadb", docid="asoda", rev=revs[1])
doc_get(x, dbname="sofadb", docid="asoda", as='json')
doc_get(x, dbname="sofadb", docid="asoda", revs=TRUE)
doc_get(x, dbname="sofadb", docid="asoda", revs=TRUE, local_seq=TRUE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_query.R
\name{db_query}
\alias{db_query}
\title{Query a database.}
\usage{
db_query(
  cushion,
  dbname,
  query = NULL,
  selector = NULL,
  limit = NULL,
  skip = NULL,
  sort = NULL,
  fields = NULL,
  use_index = NULL,
  r = NULL,
  bookmark = NULL,
  update = NULL,
  stable = NULL,
  stale = NULL,
  execution_stats = FALSE,
  as = "list",
  ...
)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name}

\item{query}{(character) instead of using the other parameters, you can
compose one R list or json blob here}

\item{selector}{(list/json) - JSON object describing criteria used to select
documents. More information provided in the section on selector syntax.
See the \code{query_tutorial} in this package, and the selectors docs
\url{http://docs.couchdb.org/en/2.0.0/api/database/find.html#find-selectors}}

\item{limit}{(numeric) - Maximum number of results returned. Default is 25.
Optional}

\item{skip}{(numeric) - Skip the first 'n' results, where 'n' is the value
specified. Optional}

\item{sort}{(json) - JSON array following sort syntax. Optional.
See \url{http://docs.couchdb.org/en/2.0.0/api/database/find.html#find-sort}
For some reason, sort doesn't work often, not sure why.}

\item{fields}{(json) - JSON array specifying which fields of each object
should be returned. If it is omitted, the entire object is returned. More
information provided in the section on filtering fields. Optional
See \url{http://docs.couchdb.org/en/2.0.0/api/database/find.html#find-filter}}

\item{use_index}{(json) - Instruct a query to use a specific index.
Specified either as \verb{<design_document>} or \verb{["<design_document>", "<index_name>"]}. Optional}

\item{r}{(numeric) Read quorum needed for the result. This defaults to 1,
in which case the document found in the index is returned. If set to a
higher value, each document is read from at least that many replicas before
it is returned in the results. This is likely to take more time than using
only the document stored locally with the index. Optional, default: 1}

\item{bookmark}{(character) A string that enables you to specify which
page of results you require. Used for paging through result sets. Every
query returns an opaque string under the bookmark key that can then be
passed back in a query to get the next page of results. If any part of
the selector query changes between requests, the results are undefined.
Optional, default: \code{NULL}}

\item{update}{(logical) Whether to update the index prior to returning the
result. Default is true. Optional}

\item{stable}{(logical) Whether or not the view results should be returned
from a “stable” set of shards. Optional}

\item{stale}{(character) Combination of \code{update=FALSE} and \code{stable=TRUE}
options. Possible options: "ok", "FALSE" (default). Optional}

\item{execution_stats}{(logical) Include execution statistics in the query
response. Optional, default: \code{FALSE}}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Query a database.
}
\examples{
\dontrun{
## create a connection
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

file <- system.file("examples/omdb.json", package = "sofa")
strs <- readLines(file)

## create a database
if ("omdb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="omdb"))
}
db_create(x, dbname='omdb')

## add some documents
invisible(db_bulk_create(x, "omdb", strs))

## query all in one json blob
db_query(x, dbname = "omdb", query = '{
  "selector": {
    "_id": {
      "$gt": null
    }
  }
}')

## query with each parameter
db_query(x, dbname = "omdb",
  selector = list(`_id` = list(`$gt` = NULL)))

db_query(x, dbname = "omdb",
  selector = list(`_id` = list(`$gt` = NULL)), limit = 3)

# fields
## single field works
db_query(x, dbname = "omdb",
  selector = list(`_id` = list(`$gt` = NULL)), limit = 3,
  fields = c('_id', 'Actors', 'imdbRating'))

## as well as many fields
db_query(x, dbname = "omdb",
  selector = list(`_id` = list(`$gt` = NULL)), limit = 3,
  fields = '_id')

## other queries
db_query(x, dbname = "omdb",
  selector = list(Year = list(`$gt` = "2013")))

db_query(x, dbname = "omdb", selector = list(Rated = "R"))

db_query(x, dbname = "omdb",
  selector = list(Rated = "PG", Language = "English"))

db_query(x, dbname = "omdb", selector = list(
  `$or` = list(
    list(Director = "Jon Favreau"),
    list(Director = "Spike Lee")
  )
), fields = c("_id", "Director"))

## when selector vars are of same name, use a JSON string
## b/c R doesn't let you have a list with same name slots
db_query(x, dbname = "omdb", query = '{
  "selector": {
    "Year": {"$gte": "1990"},
    "Year": {"$lte": "2000"},
    "$not": {"Year": "1998"}
  },
  "fields": ["_id", "Year"]
}')

## regex
db_query(x, dbname = "omdb", selector = list(
  Director = list(`$regex` = "^R")
), fields = c("_id", "Director"))

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_explain.R
\name{db_explain}
\alias{db_explain}
\title{Explain API}
\usage{
db_explain(
  cushion,
  dbname,
  query = NULL,
  selector = NULL,
  limit = NULL,
  skip = NULL,
  sort = NULL,
  fields = NULL,
  use_index = NULL,
  as = "list",
  ...
)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name}

\item{query}{(character) instead of using the other parameters, you can
compose one R list or json blob here}

\item{selector}{(json) - JSON object describing criteria used to select
documents. More information provided in the section on selector syntax.
See the \code{query_tutorial} in this package, and the selectors docs
\url{http://docs.couchdb.org/en/2.0.0/api/database/find.html#find-selectors}}

\item{limit}{(number) - Maximum number of results returned. Default: 25
Optional}

\item{skip}{(number) - Skip the first 'n' results, where 'n' is the value
specified. Optional}

\item{sort}{(json) - JSON array following sort syntax. Optional.
See \url{http://docs.couchdb.org/en/2.0.0/api/database/find.html#find-sort}
For some reason, sort doesn't work often, not sure why.}

\item{fields}{(json) - JSON array specifying which fields of each object
should be returned. If it is omitted, the entire object is returned. More
information provided in the section on filtering fields. Optional
See \url{http://docs.couchdb.org/en/2.0.0/api/database/find.html#find-filter}}

\item{use_index}{(json) - Instruct a query to use a specific index.
Specified either as \verb{<design_document>} or \verb{["<design_document>", "<index_name>"]}. Optional}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Explain API
}
\examples{
\dontrun{
## create a connection
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

file <- system.file("examples/omdb.json", package = "sofa")
strs <- readLines(file)

## create a database
if ("omdb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="omdb"))
}
db_create(x, dbname='omdb')

## add some documents
invisible(db_bulk_create(x, "omdb", strs))

## query all in one json blob
db_explain(x, dbname = "omdb", query = '{
  "selector": {
    "_id": {
      "$gt": null
    }
  }
}')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_bulk_update.R
\name{db_bulk_update}
\alias{db_bulk_update}
\title{Create documents via the bulk API}
\usage{
db_bulk_update(
  cushion,
  dbname,
  doc,
  docid = NULL,
  how = "rows",
  as = "list",
  ...
)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. Required.}

\item{doc}{For now, a data.frame only. Required.}

\item{docid}{Document IDs, ignored for now, eventually, you can pass in a
list, or vector to be the ids for each document created. Has to be the same
length as the number of documents.}

\item{how}{(character) One of rows (default) or columns. If rows, each row
becomes a separate document; if columns, each column becomes a separate
document.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
Either a list or json (depending on \code{as} parameter), with
each element an array of key:value pairs:
\itemize{
\item ok - whether creation was successful
\item id - the document id
\item rev - the revision id
}
}
\description{
Create documents via the bulk API
}
\examples{
\dontrun{
# initialize a CouchDB connection
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

row.names(mtcars) <- NULL

if ("bulktest" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="bulktest"))
}
db_create(x, dbname="bulktest")
db_bulk_create(x, mtcars, dbname="bulktest")

# modify mtcars
mtcars$letter <- sample(letters, NROW(mtcars), replace = TRUE)
db_bulk_update(x, "bulktest", mtcars)

# change again
mtcars$num <- 89
db_bulk_update(x, "bulktest", mtcars)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/databases.r
\name{databases}
\alias{databases}
\title{Work with databases in your CouchDB's.}
\description{
Work with databases in your CouchDB's.
}
\details{
There are the following functions for working with databases:
\itemize{
\item \code{\link[=db_create]{db_create()}} - Create a database
\item \code{\link[=db_delete]{db_delete()}} - Delete a database
\item \code{\link[=db_info]{db_info()}} - Get info for a database
\item \code{\link[=db_list]{db_list()}} - List databases
\item \code{\link[=db_replicate]{db_replicate()}} - Replicate a database from one couch to another
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sofa-package.r
\docType{package}
\name{sofa-package}
\alias{sofa-package}
\alias{sofa}
\title{R client for CouchDB.}
\description{
Relax.
}
\section{About sofa}{

\pkg{sofa} provides an interface to the NoSQL database CouchDB
(\url{http://couchdb.apache.org}). Methods are provided for managing
databases within CouchDB, including creating/deleting/updating/transferring,
and managing documents within databases. One can connect with a local
CouchDB instance, or a remote CouchDB databases such as Cloudant
(\url{https://cloudant.com}). Documents can be inserted directly from
vectors, lists, data.frames, and JSON.
}

\section{Client connections}{

All functions take as their first parameter a client connection object,
or a \strong{cushion}. Create the object with \link{Cushion}. You
can have multiple connection objects in an R session.
}

\section{CouchDB versions}{

\pkg{sofa} was built assuming CouchDB version 2 or greater. Some
functionality of this package will work with versions < 2, while
some may not (mango queries, see \code{\link[=db_query]{db_query()}}). I don't
plan to support older CouchDB versions per se.
}

\section{Digits after the decimal}{

If you have any concern about number of digits after the decimal
in your documents, make sure to look at \code{digits} in your R options.
The default value is 7 (see \link{options} for more informnation). You
can set the value you like with e.g., \code{options(digits = 10)}, and
get what \code{digits} is set to with \code{getOption("digits")}.

Note that in \code{\link[=doc_create]{doc_create()}} we convert your document to JSON with
\code{jsonlite::toJSON()} if given as a list, which has a \code{digits} parameter.
We pass \code{getOption("digits")} to the \code{digits} parameter in
\code{jsonlite::toJSON()}.
}

\section{Defunct functions}{

\itemize{
\item \link{attach_get}
}
}

\author{
Scott Chamberlain \email{myrmecocystus@gmail.com}
}
\keyword{package}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/doc_create.r
\name{doc_create}
\alias{doc_create}
\title{Create documents to a database.}
\usage{
doc_create(cushion, dbname, doc, docid = NULL, how = "rows", as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name3}

\item{doc}{Document content, can be character string or a list.
The character type can be XML as well, if embedded in JSON. When
the document is retrieved via \code{\link[=doc_get]{doc_get()}}, the XML is given back and
you can parse it as normal.}

\item{docid}{Document ID}

\item{how}{(character) One of rows (default) or columns. If rows, each row
becomes a separate document; if columns, each column becomes a separate
document.}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Create documents to a database.
}
\details{
Documents can have attachments just like email. There are two ways
to use attachments: the first one is via a separate REST call
(see \code{\link[=doc_attach_create]{doc_attach_create()}}); the second is inline within your
document, you can do so with this fxn. See
https://docs.couchdb.org/en/latest/api/document/attachments.html for help
on formatting json appropriately.

Note that you can create documents from a data.frame with this function,
where each row or column is a separate document. However, this function
does not use the bulk API
\url{https://couchdb.readthedocs.org/en/latest/api/database/bulk-api.html#db-bulk-docs}
\itemize{
\item see \code{\link[=db_bulk_create]{db_bulk_create()}} and \code{\link[=db_bulk_update]{db_bulk_update()}} to
create or update documents with the bulk API - which should be much faster
for a large number of documents.
}
}
\section{Digits after the decimal}{

If you have any concern about number of digits after the decimal
in your documents, make sure to look at \code{digits} in your R options.
The default value is 7 (see \link{options} for more informnation). You
can set the value you like with e.g., \code{options(digits = 10)}, and
get what \code{digits} is set to with \code{getOption("digits")}.

Note that in \code{\link[=doc_create]{doc_create()}} we convert your document to JSON with
\code{jsonlite::toJSON()} if given as a list, which has a \code{digits} parameter.
We pass \code{getOption("digits")} to the \code{digits} parameter in
\code{jsonlite::toJSON()}
}

\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("sofadb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="sofadb"))
}
db_create(x, 'sofadb')

# write a document WITH a name (uses PUT)
doc1 <- '{"name": "drink", "beer": "IPA", "score": 5}'
doc_create(x, dbname="sofadb", doc1, docid="abeer")
doc_create(x, dbname="sofadb", doc1, docid="morebeer", as='json')
doc_get(x, dbname = "sofadb", docid = "abeer")
## with factor class values
doc2 <- list(name = as.factor("drink"), beer = "stout", score = 4)
doc_create(x, doc2, dbname="sofadb", docid="nextbeer", as='json')
doc_get(x, dbname = "sofadb", docid = "nextbeer")

# write a json document WITHOUT a name (uses POST)
doc2 <- '{"name": "food", "icecream": "rocky road"}'
doc_create(x, doc2, dbname="sofadb")
doc3 <- '{"planet": "mars", "size": "smallish"}'
doc_create(x, doc3, dbname="sofadb")
## assigns a UUID instead of a user given name
db_alldocs(x, dbname = "sofadb")

# write an xml document WITH a name (uses PUT). xml is written as xml in
# couchdb, just wrapped in json, when you get it out it will be as xml
doc4 <- "<top><a/><b/><c><d/><e>bob</e></c></top>"
doc_create(x, doc4, dbname="sofadb", docid="somexml")
doc_get(x, dbname = "sofadb", docid = "somexml")

# You can pass in lists that autoconvert to json internally
doc1 <- list(name = "drink", type = "soda", score = 9)
doc_create(x, dbname="sofadb", doc1, docid="gooddrink")

# Write directly from a data.frame
## Each row or column becomes a separate document
### by rows
if ("test" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="test"))
}
db_create(x, dbname = "test")
doc_create(x, mtcars, dbname="test", how="rows")
doc_create(x, mtcars, dbname="test", how="columns")

if ("testiris" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="testiris"))
}
db_create(x, dbname = "testiris")
head(iris)
doc_create(x, iris, dbname = "testiris")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parse_df.R
\name{parse_df}
\alias{parse_df}
\title{Parse data.frame to json or list by row or column}
\usage{
parse_df(dat, how = "rows", tojson = TRUE, ...)
}
\arguments{
\item{dat}{(data.frame) A data.frame, matrix, or tbl_df}

\item{how}{(character) One of rows (default) or columns. If rows, each
row becomes a separate document; if columns, each column becomes a
separate document.}

\item{tojson}{(logical) If \code{TRUE} (default) convert to json - if \code{FALSE},
to lists}

\item{...}{Further args passed on to \code{jsonlite::toJSON()}}
}
\description{
Parse data.frame to json or list by row or column
}
\details{
Parse data.frame to get either rows or columns, each as a list
or json string
}
\examples{
\dontrun{
parse_df(mtcars, how="rows")
parse_df(mtcars, how="columns")
parse_df(mtcars, how="rows", tojson=FALSE)
parse_df(mtcars, how="columns", tojson=FALSE)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_list.r
\name{db_list}
\alias{db_list}
\title{List all databases.}
\usage{
db_list(cushion, simplify = TRUE, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{simplify}{(logical) Simplify to character vector, ignored
if \code{as="json"}}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
List all databases.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

db_list(x)
db_list(x, as = 'json')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design.R
\name{design}
\alias{design}
\alias{design_create}
\alias{design_create_}
\alias{design_delete}
\alias{design_get}
\alias{design_head}
\alias{design_info}
\title{Work with design documents}
\usage{
design_create(
  cushion,
  dbname,
  design,
  fxnname,
  key = "null",
  value = "doc",
  as = "list",
  ...
)

design_create_(cushion, dbname, design, fxnname, fxn, as = "list", ...)

design_delete(cushion, dbname, design, as = "list", ...)

design_get(cushion, dbname, design, as = "list", ...)

design_head(cushion, dbname, design, ...)

design_info(cushion, dbname, design, ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. required.}

\item{design}{(character) Design document name. this is the design name
without \verb{_design/}, which is prepended internally. required.}

\item{fxnname}{(character) A function name. required for \code{view_put}
and \code{view_put_}}

\item{key, value}{(character) a key and value, see Examples and Details}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}

\item{fxn}{(character) a javascript function. required for \code{view_put_}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Work with design documents
}
\details{
\code{design_create} is a slightly easier interface to creating
design documents; it just asks for a function name, the key and a
value, then we create the function for you internally. TO have more
flexibility use \code{view_put_} (with underscore on the end) to write the
function yourself.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

file <- system.file("examples/omdb.json", package = "sofa")
strs <- readLines(file)

## create a database
if ("omdb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="omdb"))
}
db_create(x, dbname='omdb')

## add the documents
invisible(db_bulk_create(x, "omdb", strs))

# Create a view, the easy way, but less flexible
design_create(x, dbname='omdb', design='view1', fxnname="foobar1")
design_create(x, dbname='omdb', design='view2', fxnname="foobar2",
  value="doc.Country")
design_create(x, dbname='omdb', design='view5', fxnname="foobar3",
  value="[doc.Country,doc.imdbRating]")

# the harder way, write your own function, but more flexible
design_create_(x, dbname='omdb', design='view22',
  fxnname = "stuffthings", fxn = "function(doc){emit(null,doc.Country)}")

# Delete a view
design_delete(x, dbname='omdb', design='view1')

# Get info on a design document
## HEAD request, returns just response headers
design_head(x, dbname='omdb', design='view2')
design_head(x, dbname='omdb', design='view5')
## GET request, returns information about the design document
design_info(x, dbname='omdb', design='view2')
design_info(x, dbname='omdb', design='view5')

# Get a design document (GET request)
design_get(x, dbname='omdb', design='view2')
design_get(x, dbname='omdb', design='view5')

# Search using a view
res <- design_search(x, dbname='omdb', design='view2', view='foobar2')
head(
  do.call(
    "rbind.data.frame",
    lapply(res$rows, function(x) Filter(length, x))
  )
)

res <- design_search(x, dbname='omdb', design='view5', view='foobar3')
head(
  structure(do.call(
    "rbind.data.frame",
    lapply(res$rows, function(x) x$value)
  ), .Names = c('Country', 'imdbRating'))
)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/design_search.R
\name{design_search}
\alias{design_search}
\alias{design_search_many}
\title{Search design documents}
\usage{
design_search(
  cushion,
  dbname,
  design,
  view,
  params = list(),
  body = list(),
  as = "list",
  ...
)

design_search_many(cushion, dbname, design, view, queries, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{(character) Database name. required.}

\item{design}{(character) Design document name. this is the design name
without \verb{_design/}, which is prepended internally. required.}

\item{view}{(character) a view, same as \code{fxn} param in
\code{\link[=design_create_]{design_create_()}}. required.}

\item{params}{query parameters. a named list}

\item{body}{same as \code{params}, but if any given, a POST request is
sent (if body non-NULL, \code{params} also sent). a named list}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}

\item{queries}{a list of named lists of queries}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Search design documents
}
\section{Options to pass to \code{params}, \code{body}, or \code{queries} params}{

\itemize{
\item conflicts (logical) Includes conflicts information in response.
Ignored if include_docs isn't \code{TRUE}. Default: \code{FALSE}
\item descending (logical) Return the documents in descending by key
order. Default: \code{FALSE}
\item endkey,end_key (list) Stop returning records when the specified key is
reached. Optional. \code{end_key} is an alias for \code{endkey}
\item endkey_docid,end_key_doc_id (character) Stop returning records when
the specified document ID is reached. Requires endkey to be specified for
this to have any effect. Optional. \code{end_key_doc_id} is an alias for
\code{endkey_docid}
\item group (logical) Group the results using the reduce function to a
group or single row. Default: \code{FALSE}
\item group_level (integer) Specify the group level to be used. Optional
\item include_docs (logical) Include the associated document with each
row. Default: \code{FALSE}.
\item attachments (logical) Include the Base64-encoded content of
attachments in the documents that are included if include_docs is
\code{TRUE}. Ignored if include_docs isn't \code{TRUE}.
Default: \code{FALSE}
\item att_encoding_info (logical) Include encoding information in
attachment stubs if include_docs is \code{TRUE} and the particular
attachment is compressed. Ignored if include_docs isn't \code{TRUE}.
Default: \code{FALSE}.
\item inclusive_end (logical) Specifies whether the specified end key
should be included in the result. Default: \code{TRUE}
\item key (list) Return only documents that match the specified
key. Optional
\item keys (list) Return only documents where the key matches one of the
keys specified in the array. Optional
\item limit (integer) Limit the number of the returned documents to the
specified number. Optional
\item reduce (logical)  Use the reduction function. Default: \code{TRUE}
\item skip (integer)  Skip this number of records before starting to
return the results. Default: 0
\item sorted (logical)  Sort returned rows (see Sorting Returned Rows).
Setting this to \code{FALSE} offers a performance boost. The total_rows
and offset fields are not available when this is set to \code{FALSE}.
Default: \code{TRUE}
\item stale (character) Allow the results from a stale view to be used.
Supported values: ok and update_after. Optional
\item startkey,start_key (list) Return records starting with the specified
key. Optional. \code{start_key} is an alias for startkey
\item startkey_docid,start_key_doc_id (character) Return records starting
with the specified document ID. Requires startkey to be specified for this
to have any effect. Optional. \code{start_key_doc_id} is an alias for
\code{startkey_docid}
\item update_seq (logical) Response includes an update_seq value
indicating which sequence id of the database the view reflects.
Default: \code{FALSE}
}
}

\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

file <- system.file("examples/omdb.json", package = "sofa")
strs <- readLines(file)

## create a database
if ("omdb" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="omdb"))
}
db_create(x, dbname='omdb')

## add the documents
invisible(db_bulk_create(x, "omdb", strs))

# Create a view, the easy way, but less flexible
design_create(x, dbname='omdb', design='view1', fxnname="foobar1")
design_create(x, dbname='omdb', design='view2', fxnname="foobar2",
  value="doc.Country")
design_create(x, dbname='omdb', design='view5', fxnname="foobar3",
  value="[doc.Country,doc.imdbRating]")
design_create_(x, dbname='omdb', design='view6', fxnname="foobar4",
  fxn = "function(doc){emit(doc._id,doc.Country)}")

# Search using a view
compact <- function(l) Filter(Negate(is.null), l)
res <- design_search(x, dbname='omdb', design='view2', view ='foobar2')
head(
  do.call(
    "rbind.data.frame",
    Filter(
      function(z) length(z) == 2,
      lapply(res$rows, function(x) compact(x[names(x) \%in\% c('id', 'value')]))
    )
  )
)

res <- design_search(x, dbname='omdb', design='view5', view = 'foobar3')
head(
  structure(do.call(
    "rbind.data.frame",
    lapply(res$rows, function(x) x$value)
  ), .Names = c('Country', 'imdbRating'))
)

# query parameters
## limit
design_search(x, dbname='omdb', design='view5', view = 'foobar3',
  params = list(limit = 5))
## limit and skip
design_search(x, dbname='omdb', design='view5', view = 'foobar3',
  params = list(limit = 5, skip = 3))
## with start and end keys
### important: the key strings have to be in JSON, so here e.g., 
###  need to add escaped double quotes
res <- design_search(
  cushion = x,
  dbname = 'omdb',
  design = 'view6',
  view = 'foobar4',
  params = list(
    startkey = "\"c25bbf4fef99408b3e1115374a03f642\"",
    endkey = "\"c25bbf4fef99408b3e1115374a040f11\""
  )
)

# POST request
ids <- vapply(db_alldocs(x, dbname='omdb')$rows[1:3], "[[", "", "id")
res <- design_search(x, dbname='omdb', design='view6', view = 'foobar4',
  body = list(keys = ids), verbose = TRUE)
res

# Many queries at once in a POST request
queries <- list(
  list(keys = ids),
  list(limit = 3, skip = 2)
)
design_search_many(x, 'omdb', 'view6', 'foobar4', queries)
}
}
\references{
https://docs.couchdb.org/en/latest/api/ddoc/views.html
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documents.r
\name{documents}
\alias{documents}
\title{Work with documents in your CouchDB's.}
\description{
Work with documents in your CouchDB's.
}
\details{
If you are writing a complicated javascript function, better to do
that in the Futon CouchDB interface or otherwise.

There are the following functions for working with documents:
\itemize{
\item \code{doc_create} - Create a document, with or without an ID
\item \code{doc_update} - Update a document
\item \code{doc_get} - Get a document
\item \code{doc_delete} - Delete a document
\item \code{doc_head} - Get headers for a document
\item \code{doc_attach_create} - Attach something to a document
\item \code{doc_attach_info} - Get info on an attachment
\item \code{doc_attach_get} - Fetch an attachment
\item \code{doc_attach_delete} - Delete an attachment
\item \code{db_alldocs} - Get all documents
\item \code{db_revisions} - Get revisions for a document
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/db_create.r
\name{db_create}
\alias{db_create}
\title{Create a database.}
\usage{
db_create(cushion, dbname, delifexists = FALSE, as = "list", ...)
}
\arguments{
\item{cushion}{A \code{\link{Cushion}} object. Required.}

\item{dbname}{Database name}

\item{delifexists}{If \code{TRUE}, delete any database of the same name before
creating it. This is useful for testing. Default: \code{FALSE}}

\item{as}{(character) One of list (default) or json}

\item{...}{Curl args passed on to \code{\link[crul]{HttpClient}}}
}
\value{
JSON as a character string or a list (determined by the
\code{as} parameter)
}
\description{
Create a database.
}
\examples{
\dontrun{
user <- Sys.getenv("COUCHDB_TEST_USER")
pwd <- Sys.getenv("COUCHDB_TEST_PWD")
(x <- Cushion$new(user=user, pwd=pwd))

if ("leothetiger" \%in\% db_list(x)) {
  invisible(db_delete(x, dbname="leothetiger"))
}
db_create(x, dbname='leothetiger')

## see if its there now
db_list(x)
}
}
