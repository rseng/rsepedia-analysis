elastic
=======



[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![R-check](https://github.com/ropensci/elastic/workflows/R-check/badge.svg)](https://github.com/ropensci/elastic/actions?query=workflow%3AR-check)
[![cran checks](https://cranchecks.info/badges/worst/elastic)](https://cranchecks.info/pkgs/elastic)
[![rstudio mirror downloads](https://cranlogs.r-pkg.org/badges/elastic?color=E664A4)](https://github.com/r-hub/cranlogs.app)
[![cran version](https://www.r-pkg.org/badges/version/elastic)](https://cran.r-project.org/package=elastic)
<!-- [![codecov.io](https://codecov.io/github/ropensci/elastic/coverage.svg?branch=master)](https://codecov.io/github/ropensci/elastic?branch=master) -->

**A general purpose R interface to [Elasticsearch](https://www.elastic.co/elasticsearch/)**


## Elasticsearch info

* [Elasticsearch home page](https://www.elastic.co/elasticsearch/)
* [API docs](https://www.elastic.co/guide/en/elasticsearch/reference/current/index.html)


## Compatibility

This client is developed following the latest stable releases, currently `v7.10.0`. It is generally compatible with older versions of Elasticsearch. Unlike the [Python client](https://github.com/elastic/elasticsearch-py#compatibility), we try to keep as much compatibility as possible within a single version of this client, as that's an easier setup in R world.

## Security

You're fine running ES locally on your machine, but be careful just throwing up ES on a server with a public IP address - make sure to think about security.

* Elastic has paid products - but probably only applicable to enterprise users
* DIY security - there are a variety of techniques for securing your Elasticsearch installation. A number of resources are collected in a [blog post](https://recology.info/2015/02/secure-elasticsearch/) - tools include putting your ES behind something like Nginx, putting basic auth on top of it, using https, etc.

## Installation

Stable version from CRAN


```r
install.packages("elastic")
```

Development version from GitHub


```r
remotes::install_github("ropensci/elastic")
```


```r
library('elastic')
```

## Install Elasticsearch

* [Elasticsearch installation help](https://www.elastic.co/guide/en/elasticsearch/reference/current/install-elasticsearch.html)

__w/ Docker__

Pull the official elasticsearch image

```
# elasticsearch needs to have a version tag. We're pulling 7.10.1 here
docker pull elasticsearch:7.10.1
```

Then start up a container

```
docker run -d -p 9200:9200 elasticsearch:7.10.1
```

Then elasticsearch should be available on port 9200, try `curl localhost:9200` and you should get the familiar message indicating ES is on.

If you're using boot2docker, you'll need to use the IP address in place of localhost. Get it by doing `boot2docker ip`.

__on OSX__

+ Download zip or tar file from Elasticsearch [see here for download](https://www.elastic.co/downloads), e.g., `curl -L -O https://artifacts.elastic.co/downloads/elasticsearch/elasticsearch-7.10.0-darwin-x86_64.tar.gz`
+ Extract: `tar -zxvf elasticsearch-7.10.0-darwin-x86_64.tar.gz`
+ Move it: `sudo mv elasticsearch-7.10.0 /usr/local`
+ Navigate to /usr/local: `cd /usr/local`
+ Delete symlinked `elasticsearch` directory: `rm -rf elasticsearch`
+ Add shortcut: `sudo ln -s elasticsearch-7.10.0 elasticsearch` (replace version with your version)

You can also install via Homebrew: `brew install elasticsearch`

> Note: for the 1.6 and greater upgrades of Elasticsearch, they want you to have java 8 or greater. I downloaded Java 8 from here http://www.oracle.com/technetwork/java/javase/downloads/jdk8-downloads-2133151.html and it seemed to work great.

## Upgrading Elasticsearch

I am not totally clear on best practice here, but from what I understand, when you upgrade to a new version of Elasticsearch, place old `elasticsearch/data` and `elasticsearch/config` directories into the new installation (`elasticsearch/` dir). The new elasticsearch instance with replaced data and config directories should automatically update data to the new version and start working. Maybe if you use homebrew on a Mac to upgrade it takes care of this for you - not sure.

Obviously, upgrading Elasticsearch while keeping it running is a different thing ([some help here from Elastic](https://www.elastic.co/guide/en/elasticsearch/reference/current/setup-upgrade.html)).

## Start Elasticsearch

* Navigate to elasticsearch: `cd /usr/local/elasticsearch`
* Start elasticsearch: `bin/elasticsearch`

I create a little bash shortcut called `es` that does both of the above commands in one step (`cd /usr/local/elasticsearch && bin/elasticsearch`).

## Initialization

The function `connect()` is used before doing anything else to set the connection details to your remote or local elasticsearch store. The details created by `connect()` are written to your options for the current session, and are used by `elastic` functions.


```r
x <- connect(port = 9200)
```

> If you're following along here with a local instance of Elasticsearch, you'll use `x` below to 
do more stuff.

For AWS hosted elasticsearch, make sure to specify path = "" and the correct port - transport schema pair.


```r
connect(host = <aws_es_endpoint>, path = "", port = 80, transport_schema  = "http")
  # or
connect(host = <aws_es_endpoint>, path = "", port = 443, transport_schema  = "https")
```

If you are using Elastic Cloud or an installation with authentication (X-pack), make sure to specify path = "", user = "", pwd = "" and the correct port - transport schema pair.


```r
connect(host = <ec_endpoint>, path = "", user="test", pwd = "1234", port = 9243, transport_schema  = "https")
```

<br>

## Get some data

Elasticsearch has a bulk load API to load data in fast. The format is pretty weird though. It's sort of JSON, but would pass no JSON linter. I include a few data sets in `elastic` so it's easy to get up and running, and so when you run examples in this package they'll actually run the same way (hopefully).

I have prepare a non-exported function useful for preparing the weird format that Elasticsearch wants for bulk data loads, that is somewhat specific to PLOS data (See below), but you could modify for your purposes. See `make_bulk_plos()` and `make_bulk_gbif()` [here](https://github.com/ropensci/elastic/blob/master/R/docs_bulk.r).

### Shakespeare data

Elasticsearch provides some data on Shakespeare plays. I've provided a subset of this data in this package. Get the path for the file specific to your machine:




```r
shakespeare <- system.file("examples", "shakespeare_data.json", package = "elastic")
# If you're on Elastic v6 or greater, use this one
shakespeare <- system.file("examples", "shakespeare_data_.json", package = "elastic")
shakespeare <- type_remover(shakespeare)
```

Then load the data into Elasticsearch:

> make sure to create your connection object with `connect()`


```r
# x <- connect()  # do this now if you didn't do this above
invisible(docs_bulk(x, shakespeare))
```

If you need some big data to play with, the shakespeare dataset is a good one to start with. You can get the whole thing and pop it into Elasticsearch (beware, may take up to 10 minutes or so.):

```sh
curl -XGET https://download.elastic.co/demos/kibana/gettingstarted/shakespeare_6.0.json > shakespeare.json
curl -XPUT localhost:9200/_bulk --data-binary @shakespeare.json
```

### Public Library of Science (PLOS) data

A dataset inluded in the `elastic` package is metadata for PLOS scholarly articles. Get the file path, then load:


```r
if (index_exists(x, "plos")) index_delete(x, "plos")
plosdat <- system.file("examples", "plos_data.json", package = "elastic")
plosdat <- type_remover(plosdat)
invisible(docs_bulk(x, plosdat))
```

### Global Biodiversity Information Facility (GBIF) data

A dataset inluded in the `elastic` package is data for GBIF species occurrence records. Get the file path, then load:


```r
if (index_exists(x, "gbif")) index_delete(x, "gbif")
gbifdat <- system.file("examples", "gbif_data.json", package = "elastic")
gbifdat <- type_remover(gbifdat)
invisible(docs_bulk(x, gbifdat))
```

GBIF geo data with a coordinates element to allow `geo_shape` queries


```r
if (index_exists(x, "gbifgeo")) index_delete(x, "gbifgeo")
gbifgeo <- system.file("examples", "gbif_geo.json", package = "elastic")
gbifgeo <- type_remover(gbifgeo)
invisible(docs_bulk(x, gbifgeo))
```

### More data sets

There are more datasets formatted for bulk loading in the `sckott/elastic_data` GitHub repository. Find it at <https://github.com/sckott/elastic_data>


<br>

## Search

Search the `plos` index and only return 1 result


```r
Search(x, index = "plos", size = 1)$hits$hits
#> [[1]]
#> [[1]]$`_index`
#> [1] "plos"
#> 
#> [[1]]$`_type`
#> [1] "_doc"
#> 
#> [[1]]$`_id`
#> [1] "0"
#> 
#> [[1]]$`_score`
#> [1] 1
#> 
#> [[1]]$`_source`
#> [[1]]$`_source`$id
#> [1] "10.1371/journal.pone.0007737"
#> 
#> [[1]]$`_source`$title
#> [1] "Phospholipase C-\u03b24 Is Essential for the Progression of the Normal Sleep Sequence and Ultradian Body Temperature Rhythms in Mice"
```

Search the `plos` index, and query for _antibody_, limit to 1 result


```r
Search(x, index = "plos", q = "antibody", size = 1)$hits$hits
#> [[1]]
#> [[1]]$`_index`
#> [1] "plos"
#> 
#> [[1]]$`_type`
#> [1] "_doc"
#> 
#> [[1]]$`_id`
#> [1] "813"
#> 
#> [[1]]$`_score`
#> [1] 5.18676
#> 
#> [[1]]$`_source`
#> [[1]]$`_source`$id
#> [1] "10.1371/journal.pone.0107638"
#> 
#> [[1]]$`_source`$title
#> [1] "Sortase A Induces Th17-Mediated and Antibody-Independent Immunity to Heterologous Serotypes of Group A Streptococci"
```

## Get documents

Get document with id=4


```r
docs_get(x, index = 'plos', id = 4)
#> $`_index`
#> [1] "plos"
#> 
#> $`_type`
#> [1] "_doc"
#> 
#> $`_id`
#> [1] "4"
#> 
#> $`_version`
#> [1] 1
#> 
#> $`_seq_no`
#> [1] 4
#> 
#> $`_primary_term`
#> [1] 1
#> 
#> $found
#> [1] TRUE
#> 
#> $`_source`
#> $`_source`$id
#> [1] "10.1371/journal.pone.0107758"
#> 
#> $`_source`$title
#> [1] "Lactobacilli Inactivate Chlamydia trachomatis through Lactic Acid but Not H2O2"
```

Get certain fields


```r
docs_get(x, index = 'plos', id = 4, fields = 'id')
#> $`_index`
#> [1] "plos"
#> 
#> $`_type`
#> [1] "_doc"
#> 
#> $`_id`
#> [1] "4"
#> 
#> $`_version`
#> [1] 1
#> 
#> $`_seq_no`
#> [1] 4
#> 
#> $`_primary_term`
#> [1] 1
#> 
#> $found
#> [1] TRUE
```


## Get multiple documents via the multiget API

Same index and different document ids


```r
docs_mget(x, index = "plos", id = 1:2)
#> $docs
#> $docs[[1]]
#> $docs[[1]]$`_index`
#> [1] "plos"
#> 
#> $docs[[1]]$`_type`
#> [1] "_doc"
#> 
#> $docs[[1]]$`_id`
#> [1] "1"
#> 
#> $docs[[1]]$`_version`
#> [1] 1
#> 
#> $docs[[1]]$`_seq_no`
#> [1] 1
#> 
#> $docs[[1]]$`_primary_term`
#> [1] 1
#> 
#> $docs[[1]]$found
#> [1] TRUE
#> 
#> $docs[[1]]$`_source`
#> $docs[[1]]$`_source`$id
#> [1] "10.1371/journal.pone.0098602"
#> 
#> $docs[[1]]$`_source`$title
#> [1] "Population Genetic Structure of a Sandstone Specialist and a Generalist Heath Species at Two Levels of Sandstone Patchiness across the Strait of Gibraltar"
#> 
#> 
#> 
#> $docs[[2]]
#> $docs[[2]]$`_index`
#> [1] "plos"
#> 
#> $docs[[2]]$`_type`
#> [1] "_doc"
#> 
#> $docs[[2]]$`_id`
#> [1] "2"
#> 
#> $docs[[2]]$`_version`
#> [1] 1
#> 
#> $docs[[2]]$`_seq_no`
#> [1] 2
#> 
#> $docs[[2]]$`_primary_term`
#> [1] 1
#> 
#> $docs[[2]]$found
#> [1] TRUE
#> 
#> $docs[[2]]$`_source`
#> $docs[[2]]$`_source`$id
#> [1] "10.1371/journal.pone.0107757"
#> 
#> $docs[[2]]$`_source`$title
#> [1] "Cigarette Smoke Extract Induces a Phenotypic Shift in Epithelial Cells; Involvement of HIF1\u03b1 in Mesenchymal Transition"
```

## Parsing

You can optionally get back raw `json` from `Search()`, `docs_get()`, and `docs_mget()` setting parameter `raw=TRUE`.

For example:


```r
(out <- docs_mget(x, index = "plos", id = 1:2, raw = TRUE))
#> [1] "{\"docs\":[{\"_index\":\"plos\",\"_type\":\"_doc\",\"_id\":\"1\",\"_version\":1,\"_seq_no\":1,\"_primary_term\":1,\"found\":true,\"_source\":{\"id\":\"10.1371/journal.pone.0098602\",\"title\":\"Population Genetic Structure of a Sandstone Specialist and a Generalist Heath Species at Two Levels of Sandstone Patchiness across the Strait of Gibraltar\"}},{\"_index\":\"plos\",\"_type\":\"_doc\",\"_id\":\"2\",\"_version\":1,\"_seq_no\":2,\"_primary_term\":1,\"found\":true,\"_source\":{\"id\":\"10.1371/journal.pone.0107757\",\"title\":\"Cigarette Smoke Extract Induces a Phenotypic Shift in Epithelial Cells; Involvement of HIF1\u03b1 in Mesenchymal Transition\"}}]}"
#> attr(,"class")
#> [1] "elastic_mget"
```

Then parse


```r
jsonlite::fromJSON(out)
#> $docs
#>   _index _type _id _version _seq_no _primary_term found
#> 1   plos  _doc   1        1       1             1  TRUE
#> 2   plos  _doc   2        1       2             1  TRUE
#>                     _source.id
#> 1 10.1371/journal.pone.0098602
#> 2 10.1371/journal.pone.0107757
#>                                                                                                                                                _source.title
#> 1 Population Genetic Structure of a Sandstone Specialist and a Generalist Heath Species at Two Levels of Sandstone Patchiness across the Strait of Gibraltar
#> 2                                Cigarette Smoke Extract Induces a Phenotypic Shift in Epithelial Cells; Involvement of HIF1\u03b1 in Mesenchymal Transition
```

## Known pain points

* On secure Elasticsearch servers:
  * `HEAD` requests don't seem to work, not sure why
  * If you allow only `GET` requests, a number of functions that require
  `POST` requests obviously then won't work. A big one is `Search()`, but
  you can use `Search_uri()` to get around this, which uses `GET` instead
  of `POST`, but you can't pass a more complicated query via the body

## Screencast

A screencast introducing the package: <a href="https://vimeo.com/124659179">vimeo.com/124659179</a>

## Meta

* Please [report any issues or bugs](https://github.com/ropensci/elastic/issues)
* License: MIT
* Get citation information for `elastic` in R doing `citation(package = 'elastic')`
* Please note that this package is released with a [Contributor Code of Conduct](https://ropensci.org/code-of-conduct/). By contributing to this project, you agree to abide by its terms.
elastic 1.2.0
=============

### NEW FEATURES

* `Search()` and `Search_uri()` gain new parameter `ignore_unavailable` to determine what happens if an index name does not exist (#273)
* `connect()` gains new parameter `ignore_version`. Internally, `elastic` sometimes checks the Elasticsearch version that the user is connected to to determine what to do. may be useful when it's not possible to check the Elasticsearch version, e.g., when its not possible to ping the root route of the API  (#275)
* all docs bulk functions gain parameter `digits` that is passed down to `jsonlite::toJSON() used internally`. thus, `digits` will control the number of decimal digits used in the JSON the package creates to be bulk loaded into Elasticsearch (#279)

### MINOR IMPROVEMENTS

* fix README instructions on installing Elasticsearch from docker; there's no latest tag, so use a specific version (#277) thanks @ColinFay


elastic 1.1.0
=============

### NEW FEATURES

* types were deprecated in Elasticsearch v7 and greater, and will be removed in Elasticsearch v8 and greater. this version makes type optional in all/most functions so that users with older Elasticsearch versions can still use them, but users with v7 or v8 installations don't have to use them  (#251) (#270)
* gains new method `index_shrink()` for index shrinking (#192)
* through fixing functionality in `docs_bulk()` to allow pipline attachments to work, all `docs_bulk` methods that do http requests (i.e, not prep fxns) gain the parameter `query` to pass through query parameters to the http request, including for example `pipeline`, `_source` etc. (#253)
* `Search()` and `Search_uri()` gain the parameter `track_total_hits` (default: `TRUE`) (#262) thanks @orenov

### MINOR IMPROVEMENTS

* the `warn` parameter in `connect()` was not being used across the entire package; now all methods should capture any warnings returned in the Elasticsearch HTTP API headers  (#261)
* clarify in docs that `connect()` does not create a DBI like connection object (#265)
* fix warning in `index_analyze()` function where as is method `I()` should only be applied if the input parameter is not `NULL` - to avoid a warning (#269)

### BUG FIXES

* fix to `docs_bulk_update()`: subsetting data.frame's was not working correctly when data.frame's had only 1 column; fixed (#260)
* fix to internal method `es_ver()` in the `Elasticsearch` class to be more flexible in capturing Elasticsearch version (#268)
* require newest `crul` version, helps fix a problem with passing along authentication details (#267)


elastic 1.0.0
=============

### BREAKING CHANGE

(#87) The `connect()` function is essentially the same, with some changes, but now you pass the connection object to each function all. This indeed will break code. That's why this is a major version bump. 

There is one very big downside to this: breaks existing code. That's the big one. I do apologize for this, but I believe that is outweighed by the upsides: passing the connection object matches behavior in similar R packages (e.g., all the SQL database clients); you can now manage as many different connection objects as you like in the same R session; having the connection object as an R6 class allows us to have some simple methods on that object to ping the server, etc. In addition, all functions will error with an informative message if you don't pass the connection object as the first thing.

### NEW FEATURES

* gains new ingest functions `pipeline_create`, `pipeline_delete`, `pipeline_get`, `pipeline_simulate`, and `pipeline_attachment()` (#191) (#226) 
* gains new function `docs_delete_by_query()` and `docs_update_by_query()` to delete or update multiple documents at once, respectively; and new function `reindex()` to reindex all documents from one index to another (#237) (#195)
* now using `crul` for HTTP requests. this only should matter with respect to passing in curl options  (#168)
* recent versions of Elasticsearch are starting to include warnings in response headers for deprecations and other things. These can now be turned on or off with `connect()` (#241)
* gains new functions for the bulk API: `docs_bulk_create()`, `docs_bulk_delete()`, `docs_bulk_index()`. each of which are tailored to doing the operation in the function name: creating docs, deleting docs, or indexing docs (#183)
* gains new function `type_remover()` as a utility function to help users remove types from their files to use for bulk loading; could be used on example files in this package or user supplied files (#180)
* gains function `alias_rename()` to rename aliases

### MINOR IMPROVEMENTS

* fixed `scroll()` example that wasn't working (#228)
* rework `alias_create()` (#230)
* move initialize Elasticsearch connection section of README higher up to emphasize it in the right place (#231) thanks @mbannert
* whether you want "simple" or "complete" errors no longer sets env vars internally, but is passed through the internal error checker so that choices about type of errors for different connection objects do not affect one another (#242)
* `docs_get` gains new parameters `source_includes` and `source_excludes` to include or exclude certain fields in the returned document (#246) thanks @Jensxy
* added more examples to `index_create()` (#211)
* add examples to `Search()` and `Search_uri()` docs of how to use profiles (https://www.elastic.co/guide/en/elasticsearch/reference/current/search-profile.html) (#194)
* additional example added to `docs_bulk_prep()` for doing a mix of actions (i.e., delete, create, etc.)
* improved examples throughout package docs so that examples are more self-contained
* add `include_type_name` param in mappings fxns (#250)

### BUG FIXES

* `docs_bulk_update()` was not handling boolean values correctly. now fixed (#239) (#240) thanks to @dpmccabe

### DEPRECATED AND DEFUNCT

* the `info()` method has been moved inside of the connection object. after calling `x = connect()` you can call `x$info()`
* the `ping()` method has been marked as deprecated; instead, call `ping()` on the connection object created by a call to `connect()`


elastic 0.8.4
=============

### NEW FEATURES

* Gains new function `docs_bulk_update()` to do bulk updates to documents (#169)

### MINOR IMPROVEMENTS

* Vignettes weren't showing up on CRAN, fixed (#205)
* Added an example of using WKT in a query (#215)
* using markdown docs (#209)
* `id` is now optional in `docs_create()` - if you don't pass a document identifier Elasticsearch generates one for you (#216) thanks @jbrant
* `docs_bulk()` gains new parameter `quiet` to optionally turn off the progress bar (#202)

### BUG FIXES

* Fix to `docs_bulk()` for encoding in different locales (#223) (#224) thanks @Lchiffon
* Fix for `index_get()`: you can now only pass in one value to the `features` parameter (one of settings, mappings, or aliases) (#218) thanks @happyshows
* Fix to `index_create()` to handle a list body, in addition to a JSON body (#214) thanks @emillykkejensen
* Fix to `docs_bulk()` for document IDs as factors (#212) thanks @AMR-KELEG
* Temporary files created when using `docs_bulk()` (and taking up disk space) are cleaned up now (deleted), though if you pass in your own file paths you have to clean them up (#208) thanks @emillykkejensen


elastic 0.8.0
=============

### Scroll changes

* changed to S3 setup, with methods for `character` and 
`list`.
* first parameter of `scroll()` and `scroll_clear()` is now `x`, should 
only matter if you specified the parameter name for the first parameter
* `scroll` parameter in `scroll()` function is now `time_scroll`
* Added `asdf` (for "as data.frame") to `scroll()` to give back a
data.frame (#163)
* streaming option added to `scroll()`, see parameter `stream_opts` in the
docs and examples (#160)
* general docs improvements (#182)


### NEW FEATURES

* New functions `tasks` and `tasks_cancel` for the tasks API (#145)
* streaming option added to `Search()`, see parameter `stream_opts` in the
docs and examples. `scroll` parameter in `Search()` is now `time_scroll` 
(#160)
* New function `field_caps` (for field capabilities) - in ES v5.4 and 
greater
* New function `reindex` for the reindex ES API (#134)
* New functions `index_template_get`, `index_template_put`, 
`index_template_exists`, and `index_template_delete` for the indices 
templates ES API (#133)
* New function `index_forcemerge` for the ES index `_forcemerge`
route (#176)

### MINOR IMPROVEMENTS

* Added examples to docs for `Search` and `Search_uri` for how 
to show progress bar (#162)
* Small docs fix to `docs_bulk` to clarify what's allowed as first 
parameter input (#173)
* `docs_bulk` change to internal JSON preparation to use 
`na = "null"` and `auto_unbox = TRUE` in the `jsonlite::toJSON` 
call. This means that `NA`'s in R become `null` in the JSON 
and atomic vectors are unboxed (#174) thanks @pieterprovoost
* `mapping_create` gains `update_all_types` parameter; and new man 
file to explain how to enable fielddata if sorting needed (#164)
* `suggest` is used through query DSL instead of a route, added
example to `Search` (#102)
* Now caching internal `ping()` calls - so that after the first one
we used the cached version if called again within the same R session. 
Should help speed up some code with respect to http calls (#184) 
thanks @henfiber
* Fixes to percolate functions and docs for differences in percolate 
functionality pre v5 and post v5 (#176)
* All http requests now contain `content-type` headers, for the most part
`application/json` (#197), though functions that work with the bulk API
use `application/x-ndjson` (#186)
* docs fix to `mapping_create` egs (#199)
* README now includes example of how to connect when your ES is using X-pack 
(#185) thanks @ugosan

### BUG FIXES

* fixes for normalizing url paths (#181)
* fix to `type_exists` to work on ES versions less to and greater than 
v5 (#189)
* fix to `field_stats` to indicate that its no longer avail. in 
ES v5.4 and above - and that the `fields` parameter in ES >= v5 is 
gone (#190)



elastic 0.7.8
=============

### NEW FEATURES

* New function `docs_update()` to do partial document updates (#152)
* New function `docs_bulk_prep()` to prepare bulk format files
that you can use to load into Elasticsearch with this package, on the 
command line, or in any other context (Python, Ruby, etc.) (#154)

### MINOR IMPROVEMENTS

* We're no longer running a check that your ES server is up before
every request to the server. This makes request faster, but may lead to 
less informative errors when your server is down or in some other state
than fully operational (#149)
* Tweaks here and there to make sure `elastic` works with Elasticsearch
v5. Note that not all v5 features are included here yet. (#153)

### BUG FIXES

* `docs_bulk()` was not working on single column data.frame's. now is
working. (#151) thanks @gustavobio
* `docs_*` functions now support ids with whitespace in them. (#155)
* fixes to `docs_mget()` to fix requesting certain fields back.


elastic 0.7.6
=============

### BUG FIXES

* Allow usage of `es_base` parameter in `connect()` - Now, instead of 
`stop()` on `es_base` usage, we use its value for `es_host`. Only 
pass in one or the other of `es_base` and `es_host`, not both. 
(#146) thanks @MarcinKosinski


elastic 0.7.4
=============

### NEW FEATURES

* package gains new set of functions for working with search templates:
`Search_template()`, `Search_template_register()`, `Search_template_get()`, 
`Search_template_delete()`, and `Search_template_render()`  (#101)

### MINOR IMPROVEMENTS

* Improved documentation for `docs_delete`, `docs_get` and `docs_create` 
to list correctly that numeric and character values are accepted for 
the id parameter - before stated that numeric values allowed only (#144)
thanks @dominoFire
* Added tests for illegal characters in index names.

### BUG FIXES

* Fixed bug introduced into `Search` and related functions where 
wildcards in indeces didn't work. Turned out we url escaped twice
unintentionally. Fixed now, and more tests added for wildcards. 
(#143) thanks @martijnvanbeers

elastic 0.7.2
=============

### MINOR IMPROVEMENTS

* Changed `docs_bulk()` to always return a list, whether it's given a file,
data.frame, or list. For a file, a named list is returned, while for a 
data.frame or list an unnamed list is returned as many chunks can be processed
and we don't attempt to wrangle the list output. Inputs of data.frame and list
used to return `NULL` as we didn't return anything from the internal for loop. 
You can wrap `docs_bulk` in `invisible()` if you don't want the list printed 
(#142)

### BUG FIXES

* Fixed bug in `docs_bulk()` and `msearch()` in which base URL construction
was not done correctly (#141) thanks @steeled !

elastic 0.7.0
=============

### NEW FEATURES

* New function `scroll_clear()` to clear search contexts created when
using `scroll()` (#140)
* New function `ping()` to ping an Elasticsearch server to see if
it is up (#138)
* `connect()` gains new parameter `es_path` to specify a context path, 
e.g., the `bar` in `http://foo.com/bar` (#137)

### MINOR IMPROVEMENTS

* Change all `httr::content()` calls to parse to plain text
and UTF-8 encoding (#118)
* Added note to docs that when using `scroll()` all scores are
zero b/c scores are not calculated/tracked (#127)
* `connect()` no longer pings the ES server when run, but can
now be done separately with `ping()` (#139)
* Let http request headers be sent with all requests - set with 
`connect()` (#129)
* Added `transport_schema` param to `connect()` to specify 
http or https (#130)
* By default use UUIDs with bulk API with `docs_bulk()` (#125)

### BUG FIXES

* Fix to fail well on empty body sent by user (#119)
* Fix to `docs_bulk()` function so that user supplied `doc_ids` 
are not changed at all now (#123)

elastic 0.6.0
=============

Compatibility for many Elasticsearch versions has improved. We've tested on ES versions
from the current (`v2.1.1`) back to `v1.0.0`, and `elastic` works with all versions.
There are some functions that stop with a message with some ES versions simply 
because older versions may not have had particular ES features. Please do let us 
know if you have problems with older versions of ES, so we can improve compatibility.

### NEW FEATURES

* Added `index_settings_update()` function to allow updating index settings (#66)
* All errors from the Elasticsearch server are now given back as `JSON`. 
Error parsing has thus changed in `elastic`. We now have two levels of error
behavior: 'simple' and 'complete'. These can be set in `connect()` with the 
`errors` parameter. Simple errors give back often just that there was an error,
sometimes a message with explanation is supplied. Complete errors give 
more explanation and even the ES stack trace if supplied in the ES error 
response (#92) (#93)
* New function `msearch()` to do multi-searches. This works by defining queries 
in a file, much like is done for a file to be used in bulk loading. (#103)
* New function `validate()` to validate a search. (#105)
* New suite of functions to work with the percolator service: `percolate_count()`, 
`percolate_delete()`, `percolate_list()`, `percolate_match()`, `percolate_register()`. 
The percolator works by first storing queries into an index and then you define 
documents in order to retrieve these queries. (#106)
* New function `field_stats()` to find statistical properties of a field without 
executing a search (#107)
* Added a Code of Conduct
* New function `cat_nodeattrs()`
* New function `index_recreate()` as a convenience function that detects if an 
index exists, and if so, deletes it first, then creates it again.

### MINOR IMPROVEMENTS

* `docs_bulk()` now supports passing in document ids (to the `_id` field) 
via the parameter `doc_ids` for each input data.frame or list & supports using ids
already in data.frame's or lists (#83)
* `cat_*()` functions cleaned up. previously, some functions had parameters
that were essentially silently ignored. Those parameters dropped now
from the functions. (#96)
* Elasticsearch had for a while 'search exists' functionality (via `/_search/exists`), 
but have removed that in favor of using regular `_search` with `size=0` and 
`terminate_after=1` instead. (#104)
* New parameter `lenient` in `Search()` and `Search_uri` to allow format based 
failures to be ignored, or not ignored.
* Better error handling for `docs_get()` when gthe document isn't found

### BUG FIXES

* Fixed problems in `docs_bulk()` in the use case where users use 
the function in a for loop, for example, and indexing started over, 
replacing documents with the same id (#83)
* Fixed bug in `cat_()` functions in which they sometimes failed 
when `parse=TRUE` (#88)
* Fixed bug in `docs_bulk()` in which user supplied document IDs weren't being 
passed correctly internally (#90)
* Fixed bug in `Search()` and `Search_uri()` where multiple indices weren't 
supported, whereas they should have been - supported now (#115)

### DEFUNCT

* The following functions are now defunct: `mlt()`, `nodes_shutdown()`, `index_status()`, 
and `mapping_delete()` (#94) (#98) (#99) (#110)

elastic 0.5.0
===============

### NEW FEATURES

* Added `index_settings_update()` function to allow updating index settings (#66)

### MINOR IMPROVEMENTS

* Replace `RCurl::curlEscape()` with `curl::curl_escape()` (#81)
* Explicitly import non-base R functions (#80)

### BUG FIXES

* Fixed problems introduced with `v1` of `httr`


elastic 0.4.0
===============

### NEW FEATURES

* New function `Search_uri()` where the search is defined entirely in the URL itself. 
Especially useful for cases in which `POST` requests are forbidden, e.g, on a server
that prevents `POST` requests (which the function `Search()` uses). (#58)
* New function `nodes_shutdown()` (#23)
* `docs_bulk()` gains ability to push data into Elasticsearch via the bulk http API 
from data.frame or list objects. Previously, this function only would accept a file
formatted correctly. In addition, gains new parameters: `index` - The index name to use. 
`type` - The type name to use. `chunk_size` - Size of each chunk. (#60) (#67) (#68)

### MINOR IMPROVEMENTS

* `cat_*()` functions gain new parameters: `h` to specify what fields to return; `help` to 
output available columns, and their meanings; `bytes` to give numbers back machine 
friendly; `parse` Parse to a data.frame or not
* `cat_*()` functions can now optionally capture data returned in to a data.frame (#64)
* `Search()` gains new parameter `search_path` to set the path that is used for searching. 
The default is `_search`, but sometimes in your configuration you've setup so that 
you don't need that path, or it's a different path. (023d28762e7e1028fcb0ad17867f08b5e2c92f93)

### BUG FIXES

* In `docs_mget()` added internal checker to make sure user passes in the right combination of 
`index`, `type`, and `id` parameters, or `index` and `type_id`, or just `index_type_id` (#42)
* Made `index`, `type`, and `id` parameters required in the function `docs_get()` (#43)
* Fixed bug in `scroll()` to allow long `scroll_id`'s by passing scroll ids in the body instead 
of as query parameter (#44)
* In `Search()` function, in the `error_parser()` error parser function, check to see if 
`error` element returned in response body from Elasticsearch, and if so, parse error, if not, 
pass on body (likely empty) (#45)
* In `Search()` function, added helper function to check size and from parameter
values passed in to make sure they are numbers. (#46)
* Across all functions where `index` and `type` parameters used, now using `RCurl::curlEscape()`
to URL escape. Other parameters passed in are go through `httr` CRUD methods, and do URL escaping
for us. (#49)
* Fixed links to development repo in DESCRIPTION file

elastic 0.3.0
===============

First version to go to CRAN.

### NEW FEATURES

* Added a function `scroll()` and a `scroll` parameter to the `Search()` function (#36)
* Added the function `explain()` to easily get at explanation of search results.
* Added a help file added to help explain timem and distance units. See `?units-time` and 
`?units=distance`
* New help file added to list and explain the various search functions. See `?searchapis`
* New function `tokenizer_set()` to set tokenizers
* `connect()` run on package load to set default base url of `localhost` and port of `9200` - 
you can override this by running that fxn yourself, or storing `es_base`, `es_port`, etc. 
in your `.Rprofile` file.

### IMPROVEMENTS

* Made CouchDB river plugin functions not exported for now, may bring back later. 
* Added vignettes for an intro and for search details and examples (#2)
* `es_search()` changed to `Search()`.
* More datasets included in the package for bulk data load (#16)
* All examples wrapped in `\dontrun` instead of `\donttest` so they don't fail on CRAN checks.
* `es_search_body()` removed - body based queries using the query DSL moved to the `Search()` 
function, passed into the `body` parameter.

elastic 0.2.0
===============

### IMPROVEMENTS

* Remoworked package API. Almost all functions have new names. Sorry for this major change
but it needed to be done. This brings `elastic` more in line with the official Elasticsearch
Python client (http://elasticsearch-py.readthedocs.org/en/master/).
* Similar functions are grouped together in the same manual file now to make finder related
functions easier. For example, all functions that work with indices are on the `index` manual
page, and all functions prefixed with `index_()`. Thematic manual files are: `index`, `cat`,
`cluster`, `alias`, `cdbriver`, `connect`, `documents`, `mapping`, `nodes`, and `search`.
* Note that the function `es_cat()` was changed to `cat_()` - we avoided `cat()` because as 
you know there is already a widely used function in base R, see `base::cat()`.
* We changed `cat` functions to separate functions for each command, instead of passing 
the command in as an argument. For example, `cat('aliases')` becomes `cat_aliases()`.
* The `es_` prefix remains only for `es_search()`, as we have to avoid conflict with 
`base::search()`. 
* Removed `assertthat` package import, using `stopifnot()` instead (#14)

elastic 0.1.0
===============

### NEW FEATURES

* First version.
## Test environments

* local macOS install, R 4.0.4 patched
* ubuntu 16.04 (on github actions), R 4.0.4
* win-builder (devel and release)

## R CMD check results

0 errors | 0 warnings | 0 notes

## Reverse dependencies

* I have run R CMD check on the 2 downstream dependencies, with no problems 
(<https://github.com/ropensci/elastic/blob/master/revdep/README.md>).

-------

This version adds new parameters and improves documentation.

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

* Submit an issue on the [Issues page](https://github.com/ropensci/elastic/issues)

### Code contributions

* Fork this repo to your Github account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/elastic.git`
* Make sure to track progress upstream (i.e., on our version of `elastic` at `ropensci/elastic`) by doing `git remote add upstream https://github.com/ropensci/elastic.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Push up to your account
* Submit a pull request to home base at `ropensci/elastic`

### Also, check out our [discussion forum](https://discuss.ropensci.org)
<!-- If this issue relates to usage of the package, whether a
question, bug or similar, along with your query, please paste
your devtools::session_info() or sessionInfo() into the code
block below. If not, delete all this and proceed :) 

Please also include your Elasticsearch version as behavior of this package varies with the Elasticsearch version you're connecting to.
-->

<!-- make sure to include code examples if applicable and highlight
your code in code blocks like below (remove it if you aren't
showing any code examples) 

Make absolutely sure not to include any secrets in your issue!
-->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
# Platform

|field    |value                                       |
|:--------|:-------------------------------------------|
|version  |R version 4.0.4 Patched (2021-02-17 r80031) |
|os       |macOS Big Sur 10.16                         |
|system   |x86_64, darwin17.0                          |
|ui       |X11                                         |
|language |(EN)                                        |
|collate  |en_US.UTF-8                                 |
|ctype    |en_US.UTF-8                                 |
|tz       |US/Pacific                                  |
|date     |2021-03-16                                  |

# Dependencies

|package |old   |new      |Δ  |
|:-------|:-----|:--------|:--|
|elastic |1.1.0 |1.1.9.91 |*  |
|crul    |NA    |1.1.0    |*  |

# Revdeps

*Wow, no problems at all. :)**Wow, no problems at all. :)*