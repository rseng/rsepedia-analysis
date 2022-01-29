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
