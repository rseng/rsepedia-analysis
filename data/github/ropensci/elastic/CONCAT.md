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

*Wow, no problems at all. :)**Wow, no problems at all. :)*elastic
=======

```{r echo=FALSE}
knitr::opts_chunk$set(
  comment = "#>",
  collapse = TRUE,
  warning = FALSE,
  message = FALSE
)
```

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

```{r eval=FALSE}
install.packages("elastic")
```

Development version from GitHub

```{r eval=FALSE}
remotes::install_github("ropensci/elastic")
```

```{r}
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

```{r}
x <- connect(port = 9200)
```

> If you're following along here with a local instance of Elasticsearch, you'll use `x` below to 
do more stuff.

For AWS hosted elasticsearch, make sure to specify path = "" and the correct port - transport schema pair.

```{r eval=FALSE}
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

```{r echo=FALSE}
library(elastic)
x <- connect()
if (x$es_ver() < 600) {
  shakespeare <- system.file("examples", "shakespeare_data.json", package = "elastic")
} else {
  shakespeare <- system.file("examples", "shakespeare_data_.json", package = "elastic")
  shakespeare <- type_remover(shakespeare)
}
```

```{r eval=FALSE}
shakespeare <- system.file("examples", "shakespeare_data.json", package = "elastic")
# If you're on Elastic v6 or greater, use this one
shakespeare <- system.file("examples", "shakespeare_data_.json", package = "elastic")
shakespeare <- type_remover(shakespeare)
```

Then load the data into Elasticsearch:

> make sure to create your connection object with `connect()`

```{r eval=FALSE}
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

```{r}
if (index_exists(x, "plos")) index_delete(x, "plos")
plosdat <- system.file("examples", "plos_data.json", package = "elastic")
plosdat <- type_remover(plosdat)
invisible(docs_bulk(x, plosdat))
```

### Global Biodiversity Information Facility (GBIF) data

A dataset inluded in the `elastic` package is data for GBIF species occurrence records. Get the file path, then load:

```{r}
if (index_exists(x, "gbif")) index_delete(x, "gbif")
gbifdat <- system.file("examples", "gbif_data.json", package = "elastic")
gbifdat <- type_remover(gbifdat)
invisible(docs_bulk(x, gbifdat))
```

GBIF geo data with a coordinates element to allow `geo_shape` queries

```{r}
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

```{r}
Search(x, index = "plos", size = 1)$hits$hits
```

Search the `plos` index, and query for _antibody_, limit to 1 result

```{r}
Search(x, index = "plos", q = "antibody", size = 1)$hits$hits
```

## Get documents

Get document with id=4

```{r}
docs_get(x, index = 'plos', id = 4)
```

Get certain fields

```{r}
docs_get(x, index = 'plos', id = 4, fields = 'id')
```


## Get multiple documents via the multiget API

Same index and different document ids

```{r}
docs_mget(x, index = "plos", id = 1:2)
```

## Parsing

You can optionally get back raw `json` from `Search()`, `docs_get()`, and `docs_mget()` setting parameter `raw=TRUE`.

For example:

```{r}
(out <- docs_mget(x, index = "plos", id = 1:2, raw = TRUE))
```

Then parse

```{r}
jsonlite::fromJSON(out)
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
---
title: elastic introduction
author: Scott Chamberlain
date: "2020-07-27"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{elastic introduction}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



`elastic` is an R client for [Elasticsearch](https://www.elastic.co/elasticsearch). This vignette is an introduction to the package, while other vignettes dive into the details of various topics.

## Installation

You can install from CRAN (once the package is up there)


```r
install.packages("elastic")
```

Or the development version from GitHub


```r
install.packages("remotes")
remotes::install_github("ropensci/elastic")
```

Then load the package


```r
library("elastic")
```

## Elasticsearch info

+ [Elasticsearch home page](https://www.elastic.co/)
+ [API docs](https://www.elastic.co/guide/en/elasticsearch/reference/current/index.html)

## Install Elasticsearch

* [Elasticsearch installation help](https://www.elastic.co/guide/en/elasticsearch/reference/current/install-elasticsearch.html)

__Unix (linux/osx)__

Replace `6.5.3` with the version you are working with.

+ Download zip or tar file from Elasticsearch [see here for download](https://www.elastic.co/downloads), e.g., `curl -L -O https://artifacts.elastic.co/downloads/elasticsearch/elasticsearch-6.5.3.tar.gz`
+ Extract: `tar -zxvf elasticsearch-6.5.3.tar.gz`
+ Move it: `sudo mv elasticsearch-6.5.3 /usr/local`
+ Navigate to /usr/local: `cd /usr/local`
+ Delete symlinked `elasticsearch` directory: `rm -rf elasticsearch`
+ Add shortcut: `sudo ln -s elasticsearch-6.5.3 elasticsearch` (replace version with your version)

On OSX, you can install via Homebrew: `brew install elasticsearch`

__Windows__

Windows users can follow the above, but unzip the zip file instead of uncompressing the tar file.

## Start Elasticsearch

* Navigate to elasticsearch: `cd /usr/local/elasticsearch`
* Start elasticsearch: `bin/elasticsearch`

I create a little bash shortcut called `es` that does both of the above commands in one step (`cd /usr/local/elasticsearch && bin/elasticsearch`).

__Note:__ Windows users should run the `elasticsearch.bat` file

## Initialize connection

The function `connect()` is used before doing anything else to set the connection details to your remote or local elasticsearch store. The details created by `connect()` are written to your options for the current session, and are used by `elastic` functions.


```r
x <- connect()
x
```

```
#> <Elasticsearch Connection> 
#>   transport:  http 
#>   host:       127.0.0.1 
#>   port:       9200 
#>   path:       NULL 
#>   username:   NULL 
#>   password:   NULL 
#>   errors:     simple 
#>   headers (names):   
#>   cainfo:  NULL
```

On package load, your base url and port are set to `http://127.0.0.1` and `9200`, respectively. You can of course override these settings per session or for all sessions.

## Get some data

Elasticsearch has a bulk load API to load data in fast. The format is pretty weird though. It's sort of JSON, but would pass no JSON linter. I include a few data sets in `elastic` so it's easy to get up and running, and so when you run examples in this package they'll actually run the same way (hopefully).

I have prepared a non-exported function useful for preparing the weird format that Elasticsearch wants for bulk data loads (see below). See `elastic:::make_bulk_plos` and `elastic:::make_bulk_gbif`.

### Shakespeare data

Elasticsearch provides some data on Shakespeare plays. I've provided a subset of this data in this package. Get the path for the file specific to your machine:


```r
shakespeare <- system.file("examples", "shakespeare_data.json", package = "elastic")
```

Then load the data into Elasticsearch:


```r
docs_bulk(x, shakespeare)
```

If you need some big data to play with, the shakespeare dataset is a good one to start with. You can get the whole thing and pop it into Elasticsearch (beware, may take up to 10 minutes or so.):

```sh
curl -XGET https://www.elastic.co/guide/en/kibana/3.0/snippets/shakespeare.json > shakespeare.json
curl -XPUT localhost:9200/_bulk --data-binary @shakespeare.json
```

### Public Library of Science (PLOS) data

A dataset inluded in the `elastic` package is metadata for PLOS scholarly articles. Get the file path, then load:


```r
plosdat <- system.file("examples", "plos_data.json", package = "elastic")
docs_bulk(x, plosdat)
```

### Global Biodiversity Information Facility (GBIF) data

A dataset inluded in the `elastic` package is data for GBIF species occurrence records. Get the file path, then load:


```r
gbifdat <- system.file("examples", "gbif_data.json", package = "elastic")
docs_bulk(x, gbifdat)
```

GBIF geo data with a coordinates element to allow `geo_shape` queries


```r
gbifgeo <- system.file("examples", "gbif_geo.json", package = "elastic")
docs_bulk(x, gbifgeo)
```

### More data sets

There are more datasets formatted for bulk loading in the `sckott/elastic_data` GitHub repository. Find it at [https://github.com/sckott/elastic_data](https://github.com/sckott/elastic_data)

## Search

Search the `plos` index and only return 1 result


```r
Search(x, index="plos", size=1)$hits$hits
```

```
#> [[1]]
#> [[1]]$`_index`
#> [1] "plos"
#> 
#> [[1]]$`_type`
#> [1] "article"
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
#> [1] "Phospholipase C-β4 Is Essential for the Progression of the Normal Sleep Sequence and Ultradian Body Temperature Rhythms in Mice"
```

Search the `plos` index, and the `article` document type, and query for _antibody_, limit to 1 result


```r
Search(x, index="plos", type="article", q="antibody", size=1)$hits$hits
```

```
#> [[1]]
#> [[1]]$`_index`
#> [1] "plos"
#> 
#> [[1]]$`_type`
#> [1] "article"
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

Get document with `id=1`


```r
docs_get(x, index='plos', type='article', id=1)
```

```
#> $`_index`
#> [1] "plos"
#> 
#> $`_type`
#> [1] "article"
#> 
#> $`_id`
#> [1] "1"
#> 
#> $`_version`
#> [1] 1
#> 
#> $`_seq_no`
#> [1] 1
#> 
#> $`_primary_term`
#> [1] 1
#> 
#> $found
#> [1] TRUE
#> 
#> $`_source`
#> $`_source`$id
#> [1] "10.1371/journal.pone.0098602"
#> 
#> $`_source`$title
#> [1] "Population Genetic Structure of a Sandstone Specialist and a Generalist Heath Species at Two Levels of Sandstone Patchiness across the Strait of Gibraltar"
```

Get certain fields


```r
docs_get(x, index='plos', type='article', id=1, fields='id')
```

```
#> $`_index`
#> [1] "plos"
#> 
#> $`_type`
#> [1] "article"
#> 
#> $`_id`
#> [1] "1"
#> 
#> $`_version`
#> [1] 1
#> 
#> $`_seq_no`
#> [1] 1
#> 
#> $`_primary_term`
#> [1] 1
#> 
#> $found
#> [1] TRUE
```

## Get multiple documents at once

Same index and type, different document ids


```r
docs_mget(x, index="plos", type="article", id=3:4)
```

```
#> $docs
#> $docs[[1]]
#> $docs[[1]]$`_index`
#> [1] "plos"
#> 
#> $docs[[1]]$`_type`
#> [1] "article"
#> 
#> $docs[[1]]$`_id`
#> [1] "3"
#> 
#> $docs[[1]]$`_version`
#> [1] 1
#> 
#> $docs[[1]]$`_seq_no`
#> [1] 3
#> 
#> $docs[[1]]$`_primary_term`
#> [1] 1
#> 
#> $docs[[1]]$found
#> [1] TRUE
#> 
#> $docs[[1]]$`_source`
#> $docs[[1]]$`_source`$id
#> [1] "10.1371/journal.pone.0107756"
#> 
#> $docs[[1]]$`_source`$title
#> [1] "The Effect of S-Adenosylmethionine on Cognitive Performance in Mice: An Animal Model Meta-Analysis"
#> 
#> 
#> 
#> $docs[[2]]
#> $docs[[2]]$`_index`
#> [1] "plos"
#> 
#> $docs[[2]]$`_type`
#> [1] "article"
#> 
#> $docs[[2]]$`_id`
#> [1] "4"
#> 
#> $docs[[2]]$`_version`
#> [1] 1
#> 
#> $docs[[2]]$`_seq_no`
#> [1] 4
#> 
#> $docs[[2]]$`_primary_term`
#> [1] 1
#> 
#> $docs[[2]]$found
#> [1] TRUE
#> 
#> $docs[[2]]$`_source`
#> $docs[[2]]$`_source`$id
#> [1] "10.1371/journal.pone.0107758"
#> 
#> $docs[[2]]$`_source`$title
#> [1] "Lactobacilli Inactivate Chlamydia trachomatis through Lactic Acid but Not H2O2"
```

Different indeces, types, and ids


```r
docs_mget(x, index_type_id=list(c("plos","article",1), c("gbif","record",1)))$docs[[1]]
```

```
#> $`_index`
#> [1] "plos"
#> 
#> $`_type`
#> [1] "article"
#> 
#> $`_id`
#> [1] "1"
#> 
#> $`_version`
#> [1] 1
#> 
#> $`_seq_no`
#> [1] 1
#> 
#> $`_primary_term`
#> [1] 1
#> 
#> $found
#> [1] TRUE
#> 
#> $`_source`
#> $`_source`$id
#> [1] "10.1371/journal.pone.0098602"
#> 
#> $`_source`$title
#> [1] "Population Genetic Structure of a Sandstone Specialist and a Generalist Heath Species at Two Levels of Sandstone Patchiness across the Strait of Gibraltar"
```

## Raw JSON data

You can optionally get back raw `json` from `Search()`, `docs_get()`, and `docs_mget()` setting parameter `raw=TRUE`.

For example:


```r
(out <- docs_mget(x, index="plos", type="article", id=5:6, raw=TRUE))
```

```
#> [1] "{\"docs\":[{\"_index\":\"plos\",\"_type\":\"article\",\"_id\":\"5\",\"_version\":1,\"_seq_no\":5,\"_primary_term\":1,\"found\":true,\"_source\":{\"id\":\"10.1371/journal.pone.0085123\",\"title\":\"MiR-21 Is under Control of STAT5 but Is Dispensable for Mammary Development and Lactation\"}},{\"_index\":\"plos\",\"_type\":\"article\",\"_id\":\"6\",\"_version\":1,\"_seq_no\":6,\"_primary_term\":1,\"found\":true,\"_source\":{\"id\":\"10.1371/journal.pone.0098600\",\"title\":\"Correction: Designing Mixed Species Tree Plantations for the Tropics: Balancing Ecological Attributes of Species with Landholder Preferences in the Philippines\"}}]}"
#> attr(,"class")
#> [1] "elastic_mget"
```

Then parse


```r
jsonlite::fromJSON(out)
```

```
#> $docs
#>   _index   _type _id _version _seq_no _primary_term found
#> 1   plos article   5        1       5             1  TRUE
#> 2   plos article   6        1       6             1  TRUE
#>                     _source.id
#> 1 10.1371/journal.pone.0085123
#> 2 10.1371/journal.pone.0098600
#>                                                                                                                                                     _source.title
#> 1                                                                       MiR-21 Is under Control of STAT5 but Is Dispensable for Mammary Development and Lactation
#> 2 Correction: Designing Mixed Species Tree Plantations for the Tropics: Balancing Ecological Attributes of Species with Landholder Preferences in the Philippines
```
---
title: elastic searching
author: Scott Chamberlain
date: "2020-07-27"
output: rmarkdown::html_vignette
vignette: >
    %\VignetteIndexEntry{elastic searching}
    %\VignetteEngine{knitr::rmarkdown}
    %\VignetteEncoding{UTF-8}
---



## Load elastic


```r
library("elastic")
```

## The Search function

The main interface to searching documents in your Elasticsearch store is the function `Search()`. I nearly always develop R software using all lowercase, but R has a function called `search()`, and I wanted to avoid collision with that function.

`Search()` is an interface to both the HTTP search API (in which queries are passed in the URI of the request, meaning queries have to be relatively simple), as well as the POST API, or the Query DSL, in which queries are passed in the body of the request (so can be much more complex).

There are a huge amount of ways you can search Elasticsearch documents - this tutorial covers some of them, and highlights the ways in which you interact with the R outputs.


```r
x <- connect()
```

### Search an index


```r
out <- Search(x, index="shakespeare")
out$hits$total
```

```
#> $value
#> [1] 5000
#> 
#> $relation
#> [1] "eq"
```


```r
out$hits$hits[[1]]
```

```
#> $`_index`
#> [1] "shakespeare"
#> 
#> $`_type`
#> [1] "_doc"
#> 
#> $`_id`
#> [1] "0"
#> 
#> $`_score`
#> [1] 1
#> 
#> $`_source`
#> $`_source`$line_id
#> [1] 1
#> 
#> $`_source`$play_name
#> [1] "Henry IV"
#> 
#> $`_source`$line_number
#> [1] ""
#> 
#> $`_source`$speaker
#> [1] ""
#> 
#> $`_source`$text_entry
#> [1] "ACT I"
```

### Search an index by type


```r
Search(x, index = "shakespeare")$hits$hits[[1]]
```

```
#> $`_index`
#> [1] "shakespeare"
#> 
#> $`_type`
#> [1] "_doc"
#> 
#> $`_id`
#> [1] "0"
#> 
#> $`_score`
#> [1] 1
#> 
#> $`_source`
#> $`_source`$line_id
#> [1] 1
#> 
#> $`_source`$play_name
#> [1] "Henry IV"
#> 
#> $`_source`$line_number
#> [1] ""
#> 
#> $`_source`$speaker
#> [1] ""
#> 
#> $`_source`$text_entry
#> [1] "ACT I"
```

### Return certain fields


```r
Search(x, index = "shakespeare", body = '{
  "_source": ["play_name", "speaker"]
}')$hits$hits[[1]]
```

```
#> $`_index`
#> [1] "shakespeare"
#> 
#> $`_type`
#> [1] "_doc"
#> 
#> $`_id`
#> [1] "0"
#> 
#> $`_score`
#> [1] 1
#> 
#> $`_source`
#> $`_source`$play_name
#> [1] "Henry IV"
#> 
#> $`_source`$speaker
#> [1] ""
```


### Paging


```r
Search(x, index="shakespeare", size=1, from=1)$hits
```

```
#> $total
#> $total$value
#> [1] 5000
#> 
#> $total$relation
#> [1] "eq"
#> 
#> 
#> $max_score
#> [1] 1
#> 
#> $hits
#> $hits[[1]]
#> $hits[[1]]$`_index`
#> [1] "shakespeare"
#> 
#> $hits[[1]]$`_type`
#> [1] "_doc"
#> 
#> $hits[[1]]$`_id`
#> [1] "1"
#> 
#> $hits[[1]]$`_score`
#> [1] 1
#> 
#> $hits[[1]]$`_source`
#> $hits[[1]]$`_source`$line_id
#> [1] 2
#> 
#> $hits[[1]]$`_source`$play_name
#> [1] "Henry IV"
#> 
#> $hits[[1]]$`_source`$line_number
#> [1] ""
#> 
#> $hits[[1]]$`_source`$speaker
#> [1] ""
#> 
#> $hits[[1]]$`_source`$text_entry
#> [1] "SCENE I. London. The palace."
```

### Queries

Using the `q` parameter you can pass in a query, which gets passed in the URI of the query. This type of query is less powerful than the below query passed in the body of the request, using the `body` parameter.


```r
Search(x, index="shakespeare", q="speaker:KING HENRY IV")$hits$total
```

```
#> $value
#> [1] 5000
#> 
#> $relation
#> [1] "eq"
```

#### More complex queries

Here, query for values from 10 to 20 in the field `line_id`


```r
Search(x, index="shakespeare", q="line_id:[10 TO 20]")$hits$total
```

```
#> $value
#> [1] 11
#> 
#> $relation
#> [1] "eq"
```

### Get version number for each document

Version number usually is not returned.


```r
sapply(Search(x, index="shakespeare", version=TRUE, size=2)$hits$hits, "[[", "_version")
```

```
#> [1] 1 1
```

### Get raw data


```r
Search(x, index="shakespeare", raw=TRUE)
```

```
#> [1] "{\"took\":0,\"timed_out\":false,\"_shards\":{\"total\":1,\"successful\":1,\"skipped\":0,\"failed\":0},\"hits\":{\"total\":{\"value\":5000,\"relation\":\"eq\"},\"max_score\":1.0,\"hits\":[{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"0\",\"_score\":1.0,\"_source\":{\"line_id\":1,\"play_name\":\"Henry IV\",\"line_number\":\"\",\"speaker\":\"\",\"text_entry\":\"ACT I\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"1\",\"_score\":1.0,\"_source\":{\"line_id\":2,\"play_name\":\"Henry IV\",\"line_number\":\"\",\"speaker\":\"\",\"text_entry\":\"SCENE I. London. The palace.\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"2\",\"_score\":1.0,\"_source\":{\"line_id\":3,\"play_name\":\"Henry IV\",\"line_number\":\"\",\"speaker\":\"\",\"text_entry\":\"Enter KING HENRY, LORD JOHN OF LANCASTER, the EARL of WESTMORELAND, SIR WALTER BLUNT, and others\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"3\",\"_score\":1.0,\"_source\":{\"line_id\":4,\"play_name\":\"Henry IV\",\"speech_number\":1,\"line_number\":\"1.1.1\",\"speaker\":\"KING HENRY IV\",\"text_entry\":\"So shaken as we are, so wan with care,\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"4\",\"_score\":1.0,\"_source\":{\"line_id\":5,\"play_name\":\"Henry IV\",\"speech_number\":1,\"line_number\":\"1.1.2\",\"speaker\":\"KING HENRY IV\",\"text_entry\":\"Find we a time for frighted peace to pant,\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"5\",\"_score\":1.0,\"_source\":{\"line_id\":6,\"play_name\":\"Henry IV\",\"speech_number\":1,\"line_number\":\"1.1.3\",\"speaker\":\"KING HENRY IV\",\"text_entry\":\"And breathe short-winded accents of new broils\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"6\",\"_score\":1.0,\"_source\":{\"line_id\":7,\"play_name\":\"Henry IV\",\"speech_number\":1,\"line_number\":\"1.1.4\",\"speaker\":\"KING HENRY IV\",\"text_entry\":\"To be commenced in strands afar remote.\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"7\",\"_score\":1.0,\"_source\":{\"line_id\":8,\"play_name\":\"Henry IV\",\"speech_number\":1,\"line_number\":\"1.1.5\",\"speaker\":\"KING HENRY IV\",\"text_entry\":\"No more the thirsty entrance of this soil\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"8\",\"_score\":1.0,\"_source\":{\"line_id\":9,\"play_name\":\"Henry IV\",\"speech_number\":1,\"line_number\":\"1.1.6\",\"speaker\":\"KING HENRY IV\",\"text_entry\":\"Shall daub her lips with her own childrens blood;\"}},{\"_index\":\"shakespeare\",\"_type\":\"_doc\",\"_id\":\"9\",\"_score\":1.0,\"_source\":{\"line_id\":10,\"play_name\":\"Henry IV\",\"speech_number\":1,\"line_number\":\"1.1.7\",\"speaker\":\"KING HENRY IV\",\"text_entry\":\"Nor more shall trenching war channel her fields,\"}}]}}"
```

### Curl debugging

Common options are `verbose=TRUE`, `timeout_ms=1`, `followlocation=TRUE`.


```r
out <- Search(x, index="shakespeare", verbose = TRUE)
```

### Query DSL searches - queries sent in the body of the request

Pass in as an R list


```r
mapping_create(x, "shakespeare", update_all_types = TRUE, body = '{
   "properties": {
     "text_entry": {
       "type":     "text",
       "fielddata": true
    }
  }
}')
```

```
#> $acknowledged
#> [1] TRUE
```

```r
aggs <- list(aggs = list(stats = list(terms = list(field = "text_entry"))))
Search(x, index="shakespeare", body=aggs)$hits$hits[[1]]
```

```
#> $`_index`
#> [1] "shakespeare"
#> 
#> $`_type`
#> [1] "_doc"
#> 
#> $`_id`
#> [1] "0"
#> 
#> $`_score`
#> [1] 1
#> 
#> $`_source`
#> $`_source`$line_id
#> [1] 1
#> 
#> $`_source`$play_name
#> [1] "Henry IV"
#> 
#> $`_source`$line_number
#> [1] ""
#> 
#> $`_source`$speaker
#> [1] ""
#> 
#> $`_source`$text_entry
#> [1] "ACT I"
```

Or pass in as json query with newlines, easy to read


```r
aggs <- '{
    "aggs": {
        "stats" : {
            "terms" : {
                "field" : "text_entry"
            }
        }
    }
}'
Search(x, index="shakespeare", body=aggs)$hits$hits[[1]]
```

```
#> $`_index`
#> [1] "shakespeare"
#> 
#> $`_type`
#> [1] "_doc"
#> 
#> $`_id`
#> [1] "0"
#> 
#> $`_score`
#> [1] 1
#> 
#> $`_source`
#> $`_source`$line_id
#> [1] 1
#> 
#> $`_source`$play_name
#> [1] "Henry IV"
#> 
#> $`_source`$line_number
#> [1] ""
#> 
#> $`_source`$speaker
#> [1] ""
#> 
#> $`_source`$text_entry
#> [1] "ACT I"
```

Or pass in collapsed json string


```r
aggs <- '{"aggs":{"stats":{"terms":{"field":"text_entry"}}}}'
Search(x, index="shakespeare", body=aggs)$hits$hits[[1]]
```

```
#> $`_index`
#> [1] "shakespeare"
#> 
#> $`_type`
#> [1] "_doc"
#> 
#> $`_id`
#> [1] "0"
#> 
#> $`_score`
#> [1] 1
#> 
#> $`_source`
#> $`_source`$line_id
#> [1] 1
#> 
#> $`_source`$play_name
#> [1] "Henry IV"
#> 
#> $`_source`$line_number
#> [1] ""
#> 
#> $`_source`$speaker
#> [1] ""
#> 
#> $`_source`$text_entry
#> [1] "ACT I"
```

### Aggregations

Histograms


```r
aggs <- '{
    "aggs": {
        "latbuckets" : {
           "histogram" : {
               "field" : "decimalLatitude",
               "interval" : 5
           }
        }
    }
}'
Search(x, index="gbif", body=aggs, size=0)$aggregations$latbuckets$buckets[1:3]
```

```
#> [[1]]
#> [[1]]$key
#> [1] -35
#> 
#> [[1]]$doc_count
#> [1] 1
#> 
#> 
#> [[2]]
#> [[2]]$key
#> [1] -30
#> 
#> [[2]]$doc_count
#> [1] 0
#> 
#> 
#> [[3]]
#> [[3]]$key
#> [1] -25
#> 
#> [[3]]$doc_count
#> [1] 0
```

### A bool query


```r
mmatch <- '{
 "query": {
   "bool" : {
     "must_not" : {
       "range" : {
         "speech_number" : {
           "from" : 1, "to": 5
}}}}}}'
sapply(Search(x, index="shakespeare", body=mmatch)$hits$hits, function(x) x$`_source`$speech_number)
```

```
#> [[1]]
#> NULL
#> 
#> [[2]]
#> NULL
#> 
#> [[3]]
#> NULL
#> 
#> [[4]]
#> [1] 6
#> 
#> [[5]]
#> [1] 6
#> 
#> [[6]]
#> [1] 7
#> 
#> [[7]]
#> [1] 7
#> 
#> [[8]]
#> [1] 7
#> 
#> [[9]]
#> [1] 7
#> 
#> [[10]]
#> [1] 7
```

### Fuzzy query

Fuzzy query on numerics


```r
fuzzy <- list(query = list(fuzzy = list(text_entry = "arms")))
Search(x, index="shakespeare", body = fuzzy)$hits$total
```

```
#> $value
#> [1] 49
#> 
#> $relation
#> [1] "eq"
```


```r
fuzzy <- list(query = list(fuzzy = list(text_entry = list(value = "arms", fuzziness = 4))))
Search(x, index="shakespeare", body=fuzzy)$hits$total
```

```
#> $value
#> [1] 618
#> 
#> $relation
#> [1] "eq"
```

### Range query

With numeric


```r
body <- list(query=list(range=list(decimalLongitude=list(gte=1, lte=3))))
Search(x, 'gbif', body=body)$hits$total
```

```
#> $value
#> [1] 24
#> 
#> $relation
#> [1] "eq"
```


```r
body <- list(query=list(range=list(decimalLongitude=list(gte=2.9, lte=10))))
Search(x, 'gbif', body=body)$hits$total
```

```
#> $value
#> [1] 126
#> 
#> $relation
#> [1] "eq"
```

With dates


```r
body <- list(query=list(range=list(eventDate=list(gte="2012-01-01", lte="now"))))
Search(x, 'gbif', body=body)$hits$total
```

```
#> $value
#> [1] 301
#> 
#> $relation
#> [1] "eq"
```


```r
body <- list(query=list(range=list(eventDate=list(gte="2014-01-01", lte="now"))))
Search(x, 'gbif', body=body)$hits$total
```

```
#> $value
#> [1] 292
#> 
#> $relation
#> [1] "eq"
```

### More-like-this query (more_like_this can be shortened to mlt)


```r
body <- '{
 "query": {
   "more_like_this": {
     "fields": ["abstract","title"],
     "like": "and then",
     "min_term_freq": 1,
     "max_query_terms": 12
   }
 }
}'
Search(x, 'plos', body=body)$hits$total
```

```
#> $value
#> [1] 488
#> 
#> $relation
#> [1] "eq"
```


```r
body <- '{
 "query": {
   "more_like_this": {
     "fields": ["abstract","title"],
     "like": "cell",
     "min_term_freq": 1,
     "max_query_terms": 12
   }
 }
}'
Search(x, 'plos', body=body)$hits$total
```

```
#> $value
#> [1] 58
#> 
#> $relation
#> [1] "eq"
```


### Highlighting


```r
body <- '{
 "query": {
   "query_string": {
     "query" : "cell"
   }
 },
 "highlight": {
   "fields": {
     "title": {"number_of_fragments": 2}
   }
 }
}'
out <- Search(x, 'plos', body=body)
out$hits$total
```

```
#> $value
#> [1] 58
#> 
#> $relation
#> [1] "eq"
```


```r
sapply(out$hits$hits, function(x) x$highlight$title[[1]])[8:10]
```

```
#> [1] "Functional Analysis of the Drosophila Embryonic Germ <em>Cell</em> Transcriptome by RNA Interference"
#> [2] "Diversin Is Overexpressed in Breast Cancer and Accelerates <em>Cell</em> Proliferation and Invasion" 
#> [3] "c-FLIP Protects Eosinophils from TNF-α-Mediated <em>Cell</em> Death In Vivo"
```

### Scrolling search - instead of paging


```r
Search(x, 'shakespeare', q="a*")$hits$total
```

```
#> $value
#> [1] 2747
#> 
#> $relation
#> [1] "eq"
```

```r
res <- Search(x, index = 'shakespeare', q="a*", time_scroll = "1m")
length(scroll(x, res$`_scroll_id`, time_scroll = "1m")$hits$hits)
```

```
#> [1] 10
```


```r
res <- Search(x, index = 'shakespeare', q = "a*", time_scroll = "5m")
out <- res$hits$hits
hits <- 1
while (hits != 0) {
  res <- scroll(x, res$`_scroll_id`)
  hits <- length(res$hits$hits)
  if (hits > 0)
    out <- c(out, res$hits$hits)
}
length(out)
```

```
#> [1] 2747
```

```r
res$hits$total
```

```
#> $value
#> [1] 2747
#> 
#> $relation
#> [1] "eq"
```

Woohoo! Collected all 2747 documents in very little time.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ingest.R
\name{ingest}
\alias{ingest}
\alias{pipeline_create}
\alias{pipeline_attachment}
\alias{pipeline_get}
\alias{pipeline_delete}
\alias{pipeline_simulate}
\title{Ingest API operations}
\usage{
pipeline_create(conn, id, body, ...)

pipeline_attachment(conn, index, id, pipeline, body, type = NULL, ...)

pipeline_get(conn, id, filter_path = NULL, ...)

pipeline_delete(conn, id, body, ...)

pipeline_simulate(conn, body, id = NULL, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{id}{(character) one or more pipeline id's. with delete, you can use
a wildcard match}

\item{body}{body describing pipeline, see examples and Elasticsearch docs}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}, \link[crul:verb-GET]{crul::verb-GET},
\link[crul:verb-PUT]{crul::verb-PUT}, or \link[crul:verb-DELETE]{crul::verb-DELETE}}

\item{index}{(character) an index. only used in \code{pipeline_attachment}}

\item{pipeline}{(character) a pipeline name. only used in \code{pipeline_attachment}}

\item{type}{(character) a type. only used in \code{pipeline_attachment}. by default
ths is set to \code{NULL} - optional in ES <= v6.3; not allowed in ES >= v6.4}

\item{filter_path}{(character) fields to return. deafults to all if not given}
}
\value{
a named list
}
\description{
Ingest API operations
}
\details{
ingest/pipeline functions available in Elasticsearch v5 and
greater
}
\section{Attachments}{

See https://www.elastic.co/guide/en/elasticsearch/plugins/current/ingest-attachment.html
You need to install the attachment processor plugin to be able to use
attachments in pipelines
}

\examples{
\dontrun{
# connection setup
(x <- connect())

# create
body1 <- '{
  "description" : "do a thing",
  "version" : 123,
  "processors" : [
    {
      "set" : {
        "field": "foo",
        "value": "bar"
      }
    }
  ]
}'
body2 <- '{
  "description" : "do another thing",
  "processors" : [
    {
      "set" : {
        "field": "stuff",
        "value": "things"
      }
    }
  ]
}'
pipeline_create(x, id = 'foo', body = body1)
pipeline_create(x, id = 'bar', body = body2)

# get
pipeline_get(x, id = 'foo')
pipeline_get(x, id = 'bar')
pipeline_get(x, id = 'foo', filter_path = "*.version")
pipeline_get(x, id = c('foo', 'bar')) # get >1

# delete
pipeline_delete(x, id = 'foo')

# simulate
## with pipeline included
body <- '{
  "pipeline" : {
    "description" : "do another thing",
    "processors" : [
      {
        "set" : {
          "field": "stuff",
          "value": "things"
        }
      }
    ]
  },
  "docs" : [
    { "_source": {"foo": "bar"} },
    { "_source": {"foo": "world"} }
  ]
}'
pipeline_simulate(x, body)

## referencing existing pipeline
body <- '{
  "docs" : [
    { "_source": {"foo": "bar"} },
    { "_source": {"foo": "world"} }
  ]
}'
pipeline_simulate(x, body, id = "foo")

# attchments - Note: you need the attachment plugin for this, see above
body1 <- '{
  "description" : "do a thing",
  "version" : 123,
  "processors" : [
    {
      "attachment" : {
        "field" : "data"
      }
    }
  ]
}'
pipeline_create(x, "baz", body1)
body_attach <- '{
  "data": "e1xydGYxXGFuc2kNCkxvcmVtIGlwc3VtIGRvbG9yIHNpdCBhbWV0DQpccGFyIH0="
}'
if (!index_exists(x, "boomarang")) index_create(x, "boomarang")
docs_create(x, 'boomarang', id = 1, body = list(title = "New title"))
pipeline_attachment(x, "boomarang", "1", "baz", body_attach)
pipeline_get(x, id = 'baz')
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/ingest-apis.html},
\url{https://www.elastic.co/guide/en/elasticsearch/plugins/current/using-ingest-attachment.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/elastic-package.r
\docType{package}
\name{elastic}
\alias{elastic}
\alias{elastic-package}
\title{elastic}
\description{
An Elasticsearch R client.
}
\section{About}{


This package gives you access to local or remote Elasticsearch databases.
}

\section{Quick start}{


If you're connecting to a Elasticsearch server already running, skip ahead to \strong{Search}

Install Elasticsearch (on OSX)
\itemize{
\item Download zip or tar file from Elasticsearch see here for download:
\url{https://www.elastic.co/downloads/elasticsearch}
\item Unzip it: \verb{untar elasticsearch-2.3.5.tar.gz}
\item Move it: \verb{sudo mv elasticsearch-2.3.5 /usr/local}
(replace version with your version)
\item Navigate to /usr/local: \code{cd /usr/local}
\item Add shortcut: \verb{sudo ln -s elasticsearch-2.3.5 elasticsearch}
(replace version with your version)
}

For help on other platforms, see
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/install-elasticsearch.html}

\strong{Start Elasticsearch}
\itemize{
\item Navigate to elasticsearch: \code{cd /usr/local/elasticsearch}
\item Start elasticsearch: \code{bin/elasticsearch}
}

\strong{Initialization}

The function \code{\link[=connect]{connect()}} is used before doing anything else to set
the connection details to your remote or local elasticsearch store. The
details created by \code{\link[=connect]{connect()}} are written to your options for the
current session, and are used by \code{elastic} functions.

\strong{Search}

The main way to search Elasticsearch is via the \code{\link[=Search]{Search()}} function. E.g.:

\code{Search()}
}

\section{Security}{


Elasticsearch is insecure out of the box! If you are running Elasticsearch
locally on your own machine without exposing a port to the outside world, no
worries, but if you install on a server with a public IP address, take the
necessary precautions. There are a few options:
\itemize{
\item Shield - A paid product - so probably only applicable to enterprise users
\item DIY security - there are a variety of techniques for securing your
Elasticsearch. I collected a number of resources in a blog post at
\url{https://recology.info/2015/02/secure-elasticsearch/}
}
}

\section{Elasticsearch changes}{

As of Elasticsearch v2:
\itemize{
\item You can no longer create fields with dots in the name.
\item Type names may not start with a dot (other than the special \code{.percolator} type)
\item Type names may not be longer than 255 characters
\item Types may no longer be deleted
\item Queries and filters have been merged - all filter clauses are now query clauses.
Instead, query clauses can now be used in query context or in filter context. See
examples in \code{\link[=Search]{Search()}} or \code{\link[=Search_uri]{Search_uri()}}
}
}

\section{index names}{

The following are illegal characters, and can not be used in index names or types:
\verb{\\\\}, \code{/}, \code{*}, \verb{?}, \code{<}, \code{>}, \code{|}, \verb{,} (comma). double quote and whitespace are
also illegal.
}

\author{
Scott Chamberlain
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cluster.R
\name{cluster}
\alias{cluster}
\alias{cluster_settings}
\alias{cluster_health}
\alias{cluster_state}
\alias{cluster_stats}
\alias{cluster_reroute}
\alias{cluster_pending_tasks}
\title{Elasticsearch cluster endpoints}
\usage{
cluster_settings(
  conn,
  index = NULL,
  raw = FALSE,
  callopts = list(),
  verbose = TRUE,
  ...
)

cluster_health(
  conn,
  index = NULL,
  level = NULL,
  wait_for_status = NULL,
  wait_for_relocating_shards = NULL,
  wait_for_active_shards = NULL,
  wait_for_nodes = NULL,
  timeout = NULL,
  raw = FALSE,
  callopts = list(),
  verbose = TRUE,
  ...
)

cluster_state(
  conn,
  index = NULL,
  metrics = NULL,
  raw = FALSE,
  callopts = list(),
  verbose = TRUE,
  ...
)

cluster_stats(
  conn,
  index = NULL,
  raw = FALSE,
  callopts = list(),
  verbose = TRUE,
  ...
)

cluster_reroute(conn, body, raw = FALSE, callopts = list(), ...)

cluster_pending_tasks(
  conn,
  index = NULL,
  raw = FALSE,
  callopts = list(),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{Index}

\item{raw}{If \code{TRUE} (default), data is parsed to list. If \code{FALSE}, then raw JSON.}

\item{callopts}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}

\item{verbose}{If \code{TRUE} (default) the url call used printed to console.}

\item{...}{Further args passed on to elastic search HTTP API as parameters.}

\item{level}{Can be one of cluster, indices or shards. Controls the details level of the
health information returned. Defaults to cluster.}

\item{wait_for_status}{One of green, yellow or red. Will wait (until the timeout
provided) until the status of the cluster changes to the one provided or better, i.e.
green > yellow > red. By default, will not wait for any status.}

\item{wait_for_relocating_shards}{A number controlling to how many relocating shards
to wait for. Usually will be 0 to indicate to wait till all relocations have happened.
Defaults to not wait.}

\item{wait_for_active_shards}{A number controlling to how many active shards to wait for.
Defaults to not wait.}

\item{wait_for_nodes}{The request waits until the specified number N of nodes is
available. It also accepts >=N, <=N, >N and <N. Alternatively, it is possible to use
ge(N), le(N), gt(N) and lt(N) notation.}

\item{timeout}{A time based parameter controlling how long to wait if one of the
wait_for_XXX are provided. Defaults to 30s.}

\item{metrics}{One or more of version, master_node, nodes, routing_table,
metadata, and blocks. See Details.}

\item{body}{Query, either a list or json.}
}
\description{
Elasticsearch cluster endpoints
}
\details{
metrics param options:
\itemize{
\item version Shows the cluster state version.
\item master_node Shows the elected master_node part of the response
\item nodes Shows the nodes part of the response
\item routing_table Shows the routing_table part of the response. If you supply
a comma separated list of indices, the returned output will only contain the
indices listed.
\item metadata Shows the metadata part of the response. If you supply a comma
separated list of indices, the returned output will only contain the indices
listed.
\item blocks Shows the blocks part of the response
}

Additional parameters that can be passed in:
\itemize{
\item metric A comma-separated list of metrics to display. Possible values: '_all',
'completion', 'docs', 'fielddata', 'filter_cache', 'flush', 'get', 'id_cache', 'indexing',
'merge', 'percolate', 'refresh', 'search', 'segments', 'store', 'warmer'
\item completion_fields A comma-separated list of fields for completion metric (supports
wildcards)
\item fielddata_fields A comma-separated list of fields for fielddata metric (supports
wildcards)
\item fields A comma-separated list of fields for fielddata and completion metric (supports
wildcards)
\item groups A comma-separated list of search groups for search statistics
\item allow_no_indices Whether to ignore if a wildcard indices expression resolves into no
concrete indices. (This includes _all string or when no indices have been specified)
\item expand_wildcards Whether to expand wildcard expression to concrete indices that are
open, closed or both.
\item ignore_indices When performed on multiple indices, allows to ignore missing ones
(default: none)
\item ignore_unavailable Whether specified concrete indices should be ignored when unavailable
(missing or closed)
\item human Whether to return time and byte values in human-readable format.
\item level Return stats aggregated at cluster, index or shard level. ('cluster', 'indices'
or 'shards', default: 'indices')
\item types A comma-separated list of document types for the indexing index metric
}
}
\examples{
\dontrun{
# connection setup
(x <- connect())

cluster_settings(x)
cluster_health(x)

cluster_state(x)
cluster_state(x, metrics = "version")
cluster_state(x, metrics = "nodes")
cluster_state(x, metrics = c("version", "nodes"))
cluster_state(x, metrics = c("version", "nodes", 'blocks'))
cluster_state(x, "shakespeare", metrics = "metadata")
cluster_state(x, c("shakespeare", "flights"), metrics = "metadata")

cluster_stats(x)
cluster_pending_tasks(x)

body <- '{
  "commands": [ 
    {
      "move": {
        "index" : "test", "shard" : 0,
        "from_node" : "node1", "to_node" : "node2"
      }
    },
    {
      "allocate_replica" : {
        "index" : "test", "shard" : 1, "node" : "node3"
      }
    }
  ]
}'
# cluster_reroute(x, body =  body)

cluster_health(x)
# cluster_health(x, wait_for_status = "yellow", timeout = "3s")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tasks.R
\name{tasks}
\alias{tasks}
\alias{tasks_cancel}
\title{Elasticsearch tasks endpoints}
\usage{
tasks(
  conn,
  task_id = NULL,
  nodes = NULL,
  actions = NULL,
  parent_task_id = NULL,
  detailed = FALSE,
  group_by = NULL,
  wait_for_completion = FALSE,
  timeout = NULL,
  raw = FALSE,
  ...
)

tasks_cancel(
  conn,
  node_id = NULL,
  task_id = NULL,
  nodes = NULL,
  actions = NULL,
  parent_task_id = NULL,
  detailed = FALSE,
  group_by = NULL,
  wait_for_completion = FALSE,
  timeout = NULL,
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{task_id}{a task id}

\item{nodes}{(character) The nodes}

\item{actions}{(character) Actions}

\item{parent_task_id}{(character) A parent task ID}

\item{detailed}{(character) get detailed results. Default: \code{FALSE}}

\item{group_by}{(character) "nodes" (default, i.e., NULL) or "parents"}

\item{wait_for_completion}{(logical) wait for completion. Default: \code{FALSE}}

\item{timeout}{(integer) timeout time}

\item{raw}{If \code{TRUE} (default), data is parsed to list. If \code{FALSE}, then
raw JSON.}

\item{...}{Curl args passed on to \link[crul:verb-GET]{crul::verb-GET} or
\link[crul:verb-POST]{crul::verb-POST}}

\item{node_id}{a node id}
}
\description{
Elasticsearch tasks endpoints
}
\examples{
\dontrun{
x <- connect()

tasks(x)
# tasks(x, parent_task_id = "1234")

# delete a task
# tasks_cancel(x)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/tasks.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_delete.R
\name{docs_delete}
\alias{docs_delete}
\title{Delete a document}
\usage{
docs_delete(
  conn,
  index,
  id,
  type = NULL,
  refresh = NULL,
  routing = NULL,
  timeout = NULL,
  version = NULL,
  version_type = NULL,
  callopts = list(),
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) The name of the index. Required}

\item{id}{(numeric/character) The document ID. Can be numeric or character.
Required}

\item{type}{(character) The type of the document. optional}

\item{refresh}{(logical) Refresh the index after performing the operation}

\item{routing}{(character) Specific routing value}

\item{timeout}{(character) Explicit operation timeout, e.g,. 5m (for 5
minutes)}

\item{version}{(character) Explicit version number for concurrency control}

\item{version_type}{(character) Specific version type. One of internal
or external}

\item{callopts}{Curl args passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{...}{Further args to query DSL}
}
\description{
Delete a document
}
\examples{
\dontrun{
(x <- connect())
x$ping()

if (!index_exists(x, "plos")) {
 plosdat <- system.file("examples", "plos_data.json",
    package = "elastic")
 plosdat <- type_remover(plosdat)
 docs_bulk(x, plosdat)
}

# delete a document
if (!docs_get(x, index='plos', id=36, exists=TRUE)) {
  docs_create(x, index='plos', id=36, 
    body = list(id="12345", title="New title")
  )
}
docs_get(x, index='plos', id=36)
docs_delete(x, index='plos', id=36)
# docs_get(x, index='plos', id=36) # and the document is gone
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-delete.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tokenizer_set.R
\name{tokenizer_set}
\alias{tokenizer_set}
\title{Tokenizer operations}
\usage{
tokenizer_set(conn, index, body, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) A character vector of index names}

\item{body}{Query, either a list or json.}

\item{...}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Tokenizer operations
}
\examples{
\dontrun{
# connection setup
(x <- connect())

# set tokenizer

## NGram tokenizer
body <- '{
        "settings" : {
             "analysis" : {
                 "analyzer" : {
                     "my_ngram_analyzer" : {
                         "tokenizer" : "my_ngram_tokenizer"
                     }
                 },
                 "tokenizer" : {
                     "my_ngram_tokenizer" : {
                         "type" : "nGram",
                         "min_gram" : "2",
                         "max_gram" : "3",
                         "token_chars": [ "letter", "digit" ]
                     }
                 }
             }
      }
}'
if (index_exists('test1')) index_delete('test1')
tokenizer_set(index = "test1", body=body)
index_analyze(text = "hello world", index = "test1", 
  analyzer='my_ngram_analyzer')
}
}
\references{
https://www.elastic.co/guide/en/elasticsearch/reference/current/analysis-tokenizers.html
}
\author{
Scott Chamberlain \href{mailto:myrmecocystus@gmail.com}{myrmecocystus@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nodes.R
\name{nodes}
\alias{nodes}
\alias{nodes_stats}
\alias{nodes_info}
\alias{nodes_hot_threads}
\title{Elasticsearch nodes endpoints.}
\usage{
nodes_stats(conn, node = NULL, metric = NULL, raw = FALSE, fields = NULL, ...)

nodes_info(conn, node = NULL, metric = NULL, raw = FALSE, ...)

nodes_hot_threads(
  conn,
  node = NULL,
  metric = NULL,
  threads = 3,
  interval = "500ms",
  type = NULL,
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{node}{The node}

\item{metric}{A metric to get. See Details.}

\item{raw}{If \code{TRUE} (default), data is parsed to list. If \code{FALSE}, then
raw JSON.}

\item{fields}{You can get information about field data memory usage on
node level or on index level}

\item{...}{Curl args passed on to \link[crul:verb-GET]{crul::verb-GET}}

\item{threads}{(character) Number of hot threads to provide. Default: 3}

\item{interval}{(character) The interval to do the second sampling of
threads. Default: 500ms}

\item{type}{(character) The type to sample, defaults to cpu, but supports
wait and block to see hot threads that are in wait or block state.}
}
\description{
Elasticsearch nodes endpoints.
}
\details{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/cluster-nodes-stats.html}

By default, all stats are returned. You can limit this by combining any of
indices, os, process, jvm, network, transport, http, fs, breaker and
thread_pool. With the metric parameter you can select zero or more of:
\itemize{
\item indices Indices stats about size, document count, indexing and
deletion times, search times, field cache size, merges and flushes
\item os retrieve information that concern the operating system
\item fs File system information, data path, free disk space,
read/write stats
\item http HTTP connection information
\item jvm JVM stats, memory pool information, garbage collection,
buffer pools
\item network TCP information
\item os Operating system stats, load average, cpu, mem, swap
\item process Process statistics, memory consumption, cpu usage, open
file descriptors
\item thread_pool Statistics about each thread pool, including current
size, queue and rejected tasks
\item transport Transport statistics about sent and received bytes in
cluster communication
\item breaker Statistics about the field data circuit breaker
}

\code{\link[=nodes_hot_threads]{nodes_hot_threads()}} returns plain text, so \code{\link[base:cat]{base::cat()}}
is used to print to the console.
}
\examples{
\dontrun{
# connection setup
(x <- connect())

(out <- nodes_stats(x))
nodes_stats(x, node = names(out$nodes))
nodes_stats(x, metric='get')
nodes_stats(x, metric='jvm')
nodes_stats(x, metric=c('os','process'))
nodes_info(x)
nodes_info(x, metric='process')
nodes_info(x, metric='jvm')
nodes_info(x, metric='http')
nodes_info(x, metric='network')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field_caps.R
\name{field_caps}
\alias{field_caps}
\title{Field capabilities}
\usage{
field_caps(conn, fields, index = NULL, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{fields}{A list of fields to compute stats for. required}

\item{index}{Index name, one or more}

\item{...}{Curl args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
The field capabilities API allows to retrieve the capabilities of fields
among multiple indices.
}
\examples{
\dontrun{
x <- connect()
x$ping()

if (x$es_ver() >= 540) {
  field_caps(x, fields = "speaker", index = "shakespeare")
}

}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/search-field-caps.html}
}
\seealso{
\code{\link[=field_stats]{field_stats()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{index_status}
\alias{index_status}
\title{This function is defunct}
\usage{
index_status(...)
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/preference.R
\name{preference}
\alias{preference}
\title{Preferences.}
\description{
Preferences.
}
\details{
\itemize{
\item _primary The operation will go and be executed only on the primary shards.
\item _primary_first The operation will go and be executed on the primary shard, and if
not available (failover), will execute on other shards.
\item _local The operation will prefer to be executed on a local allocated shard if possible.
\item _only_node:xyz Restricts the search to execute only on a node with the provided
node id (xyz in this case).
\item _prefer_node:xyz Prefers execution on the node with the provided node
id (xyz in this case) if applicable.
\item _shards:2,3 Restricts the operation to the specified shards. (2 and 3 in this case).
This preference can be combined with other preferences but it has to appear
first: _shards:2,3;_primary
\item Custom (string) value A custom value will be used to guarantee that the same shards
will be used for the same custom value. This can help with "jumping values" when hitting
different shards in different refresh states. A sample value can be something like the web
session id, or the user name.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_bulk_prep.R
\name{docs_bulk_prep}
\alias{docs_bulk_prep}
\title{Use the bulk API to prepare bulk format data}
\usage{
docs_bulk_prep(
  x,
  index,
  path,
  type = NULL,
  chunk_size = 1000,
  doc_ids = NULL,
  quiet = FALSE,
  digits = NA
)
}
\arguments{
\item{x}{A data.frame or a list. required.}

\item{index}{(character) The index name. required.}

\item{path}{(character) Path to the file. If data is broken into chunks,
we'll use this path as the prefix, and suffix each file path with a number.
required.}

\item{type}{(character) The type. default: \code{NULL}. Note that \code{type} is
deprecated in Elasticsearch v7 and greater, and removed in Elasticsearch v8}

\item{chunk_size}{(integer) Size of each chunk. If your data.frame is smaller
thank \code{chunk_size}, this parameter is essentially ignored. We write in
chunks because at some point, depending on size of each document, and
Elasticsearch setup, writing a very large number of documents in one go
becomes slow, so chunking can help. This parameter is ignored if you
pass a file name. Default: 1000}

\item{doc_ids}{An optional vector (character or numeric/integer) of document
ids to use. This vector has to equal the size of the documents you are
passing in, and will error if not. If you pass a factor we convert to
character. Default: not passed}

\item{quiet}{(logical) Suppress progress bar. Default: \code{FALSE}}

\item{digits}{digits used by the parameter of the same name by
\code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}} to convert data to JSON before being submitted to
your ES instance. default: \code{NA}}
}
\value{
File path(s). By default we use temporary files; these are cleaned
up at the end of a session
}
\description{
Use the bulk API to prepare bulk format data
}
\section{Tempfiles}{

In \code{docs_bulk} we create temporary files in some cases, and delete
those before the function exits. However, we don't clean up those files
in this function because the point of the function is to create the
newline delimited JSON files that you need. Tempfiles are cleaned up
when you R session ends though - be aware of that. If you want to
keep the files make sure to move them outside of the temp directory.
}

\examples{
\dontrun{
# From a data.frame
ff <- tempfile(fileext = ".json")
docs_bulk_prep(mtcars, index = "hello", path = ff)
readLines(ff)

## field names cannot contain dots
names(iris) <- gsub("\\\\.", "_", names(iris))
docs_bulk_prep(iris, "iris", path = tempfile(fileext = ".json"))

## type can be missing, but index can not
docs_bulk_prep(iris, "flowers", path = tempfile(fileext = ".json"))

# From a list
docs_bulk_prep(apply(iris, 1, as.list), index="iris",
   path = tempfile(fileext = ".json"))
docs_bulk_prep(apply(USArrests, 1, as.list), index="arrests",
   path = tempfile(fileext = ".json"))

# when chunking
## multiple files created, one for each chunk
bigiris <- do.call("rbind", replicate(30, iris, FALSE))
docs_bulk_prep(bigiris, index = "big", path = tempfile(fileext = ".json"))

# When using in a loop
## We internally get last _id counter to know where to start on next bulk
## insert but you need to sleep in between docs_bulk_prep calls, longer the
## bigger the data is
files <- c(system.file("examples", "test1.csv", package = "elastic"),
           system.file("examples", "test2.csv", package = "elastic"),
           system.file("examples", "test3.csv", package = "elastic"))
paths <- vector("list", length = length(files))
for (i in seq_along(files)) {
  d <- read.csv(files[[i]])
  paths[i] <- docs_bulk_prep(d, index = "stuff",
     path = tempfile(fileext = ".json"))
}
unlist(paths)

# You can include your own document id numbers
## Either pass in as an argument
files <- c(system.file("examples", "test1.csv", package = "elastic"),
           system.file("examples", "test2.csv", package = "elastic"),
           system.file("examples", "test3.csv", package = "elastic"))
tt <- vapply(files, function(z) NROW(read.csv(z)), numeric(1))
ids <- list(1:tt[1],
           (tt[1] + 1):(tt[1] + tt[2]),
           (tt[1] + tt[2] + 1):sum(tt))
paths <- vector("list", length = length(files))
for (i in seq_along(files)) {
  d <- read.csv(files[[i]])
  paths[i] <- docs_bulk_prep(d, index = "testes",
    doc_ids = ids[[i]], path = tempfile(fileext = ".json"))
}
unlist(paths)

## or include in the input data
### from data.frame's
files <- c(system.file("examples", "test1_id.csv", package = "elastic"),
           system.file("examples", "test2_id.csv", package = "elastic"),
           system.file("examples", "test3_id.csv", package = "elastic"))
paths <- vector("list", length = length(files))
for (i in seq_along(files)) {
  d <- read.csv(files[[i]])
  paths[i] <- docs_bulk_prep(d, index = "testes",
     path = tempfile(fileext = ".json"))
}
unlist(paths)

### from lists via file inputs
paths <- vector("list", length = length(files))
for (i in seq_along(files)) {
  d <- read.csv(files[[i]])
  d <- apply(d, 1, as.list)
  paths[i] <- docs_bulk_prep(d, index = "testes",
      path = tempfile(fileext = ".json"))
}
unlist(paths)


# A mix of actions
## make sure you use a column named 'es_action' or this won't work
## if you need to delete or update you need document IDs
if (index_exists(x, "baz")) index_delete(x, "baz")
df <- data.frame(a = 1:5, b = 6:10, c = letters[1:5], stringsAsFactors = FALSE) 
f <- tempfile(fileext = ".json")
invisible(docs_bulk_prep(df, "baz", f))
cat(readLines(f), sep = "\n")
docs_bulk(x, f)
Sys.sleep(2)
(res <- Search(x, 'baz', asdf=TRUE)$hits$hits)

df[1, "a"] <- 99
df[1, "c"] <- "aa"
df[3, "c"] <- 33
df[3, "c"] <- "cc"
df$es_action <- c('update', 'delete', 'update', 'delete', 'delete')
df$id <- res$`_id`
df
f <- tempfile(fileext = ".json")
invisible(docs_bulk_prep(df, "baz", path = f, doc_ids = df$id))
cat(readLines(f), sep = "\n")
docs_bulk(x, f)


# suppress progress bar
docs_bulk_prep(mtcars, index = "hello",
  path = tempfile(fileext = ".json"), quiet = TRUE)
## vs. 
docs_bulk_prep(mtcars, index = "hello",
  path = tempfile(fileext = ".json"), quiet = FALSE)
}
}
\seealso{
Other bulk-functions: 
\code{\link{docs_bulk_create}()},
\code{\link{docs_bulk_delete}()},
\code{\link{docs_bulk_index}()},
\code{\link{docs_bulk_update}()},
\code{\link{docs_bulk}()}
}
\concept{bulk-functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/fielddata.R
\name{fielddata}
\alias{fielddata}
\title{fielddata}
\description{
Deep dive on fielddata details
}
\details{
Most fields are indexed by default, which makes them searchable. Sorting,
aggregations, and accessing field values in scripts, however, requires a
different access pattern from search.

Text fields use a query-time in-memory data structure called fielddata.
This data structure is built on demand the first time that a field is
used for aggregations, sorting, or in a script. It is built by reading
the entire inverted index for each segment from disk, inverting the
term-document relationship, and storing the result in memory, in the
JVM heap.

fielddata is disabled on text fields by default. Fielddata can consume a
lot of heap space, especially when loading high cardinality text fields.
Once fielddata has been loaded into the heap, it remains there for the
lifetime of the segment. Also, loading fielddata is an expensive process
which can cause users to experience latency hits. This is why fielddata
is disabled by default. If you try to sort, aggregate, or access values
from a script on a text field, you will see this exception:

"Fielddata is disabled on text fields by default. Set fielddata=true on
\code{your_field_name} in order to load fielddata in memory by uninverting
the inverted index. Note that this can however use significant memory."

To enable fielddata on a text field use the PUT mapping API, for example
\code{mapping_create("shakespeare", body = '{
  "properties": {
    "speaker": { 
      "type":     "text",
      "fielddata": true
    }
  }
}')}

You may get an error about \code{update_all_types}, in which case set
\code{update_all_types=TRUE} in \code{mapping_create}, e.g.,

\code{mapping_create("shakespeare", update_all_types=TRUE, body = '{
  "properties": {
    "speaker": { 
      "type":     "text",
      "fielddata": true
    }
  }
}')}

See \url{https://www.elastic.co/guide/en/elasticsearch/reference/current/fielddata.html#_enabling_fielddata_on_literal_text_literal_fields}
for more information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/msearch.R
\name{msearch}
\alias{msearch}
\title{Multi-search}
\usage{
msearch(conn, x, raw = FALSE, asdf = FALSE, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{x}{(character) A file path}

\item{raw}{(logical) Get raw JSON back or not.}

\item{asdf}{(logical) If \code{TRUE}, use \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}}
to parse JSON directly to a data.frame. If \code{FALSE} (Default), list
output is given.}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}
}
\description{
Performs multiple searches, defined in a file
}
\details{
This function behaves similarly to \code{\link[=docs_bulk]{docs_bulk()}} -
performs searches based on queries defined in a file.
}
\examples{
\dontrun{
x <- connect()

msearch1 <- system.file("examples", "msearch_eg1.json", package = "elastic")
readLines(msearch1)
msearch(x, msearch1)

tf <- tempfile(fileext = ".json")
cat('{"index" : "shakespeare"}', file = tf, sep = "\n")
cat('{"query" : {"match_all" : {}}, "from" : 0, "size" : 5}',  sep = "\n",
   file = tf, append = TRUE)
readLines(tf)
msearch(x, tf)
}
}
\seealso{
\code{\link[=Search_uri]{Search_uri()}} \code{\link[=Search]{Search()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index.R
\name{indices}
\alias{indices}
\alias{index_get}
\alias{index_exists}
\alias{index_delete}
\alias{index_create}
\alias{index_recreate}
\alias{index_close}
\alias{index_open}
\alias{index_stats}
\alias{index_settings}
\alias{index_settings_update}
\alias{index_segments}
\alias{index_recovery}
\alias{index_optimize}
\alias{index_forcemerge}
\alias{index_upgrade}
\alias{index_analyze}
\alias{index_flush}
\alias{index_clear_cache}
\alias{index_shrink}
\title{Index API operations}
\usage{
index_get(
  conn,
  index = NULL,
  features = NULL,
  raw = FALSE,
  verbose = TRUE,
  ...
)

index_exists(conn, index, ...)

index_delete(conn, index, raw = FALSE, verbose = TRUE, ...)

index_create(conn, index = NULL, body = NULL, raw = FALSE, verbose = TRUE, ...)

index_recreate(
  conn,
  index = NULL,
  body = NULL,
  raw = FALSE,
  verbose = TRUE,
  ...
)

index_close(conn, index, ...)

index_open(conn, index, ...)

index_stats(
  conn,
  index = NULL,
  metric = NULL,
  completion_fields = NULL,
  fielddata_fields = NULL,
  fields = NULL,
  groups = NULL,
  level = "indices",
  ...
)

index_settings(conn, index = "_all", ...)

index_settings_update(conn, index = NULL, body, ...)

index_segments(conn, index = NULL, ...)

index_recovery(conn, index = NULL, detailed = FALSE, active_only = FALSE, ...)

index_optimize(
  conn,
  index = NULL,
  max_num_segments = NULL,
  only_expunge_deletes = FALSE,
  flush = TRUE,
  wait_for_merge = TRUE,
  ...
)

index_forcemerge(
  conn,
  index = NULL,
  max_num_segments = NULL,
  only_expunge_deletes = FALSE,
  flush = TRUE,
  ...
)

index_upgrade(conn, index = NULL, wait_for_completion = FALSE, ...)

index_analyze(
  conn,
  text = NULL,
  field = NULL,
  index = NULL,
  analyzer = NULL,
  tokenizer = NULL,
  filters = NULL,
  char_filters = NULL,
  body = list(),
  ...
)

index_flush(
  conn,
  index = NULL,
  force = FALSE,
  full = FALSE,
  wait_if_ongoing = FALSE,
  ...
)

index_clear_cache(
  conn,
  index = NULL,
  filter = FALSE,
  filter_keys = NULL,
  fielddata = FALSE,
  query_cache = FALSE,
  id_cache = FALSE,
  ...
)

index_shrink(conn, index, index_new, body = NULL, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) A character vector of index names}

\item{features}{(character) A single feature. One of settings, mappings, or
aliases}

\item{raw}{If \code{TRUE} (default), data is parsed to list. If FALSE, then raw JSON.}

\item{verbose}{If \code{TRUE} (default) the url call used printed to console.}

\item{...}{Curl args passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{body}{Query, either a list or json.}

\item{metric}{(character) A character vector of metrics to display. Possible
values: "_all", "completion", "docs", "fielddata", "filter_cache", "flush",
"get", "id_cache", "indexing", "merge", "percolate", "refresh", "search",
"segments", "store", "warmer".}

\item{completion_fields}{(character) A character vector of fields for completion metric
(supports wildcards)}

\item{fielddata_fields}{(character) A character vector of fields for fielddata metric
(supports wildcards)}

\item{fields}{(character) Fields to add.}

\item{groups}{(character) A character vector of search groups for search statistics.}

\item{level}{(character) Return stats aggregated on "cluster", "indices" (default) or "shards"}

\item{detailed}{(logical) Whether to display detailed information about shard recovery.
Default: \code{FALSE}}

\item{active_only}{(logical) Display only those recoveries that are currently on-going.
Default: \code{FALSE}}

\item{max_num_segments}{(character) The number of segments the index should be merged into.
Default: "dynamic"}

\item{only_expunge_deletes}{(logical) Specify whether the operation should only expunge
deleted documents}

\item{flush}{(logical) Specify whether the index should be flushed after performing the
operation. Default: \code{TRUE}}

\item{wait_for_merge}{(logical) Specify whether the request should block until the merge
process is finished. Default: \code{TRUE}}

\item{wait_for_completion}{(logical) Should the request wait for the upgrade to complete.
Default: \code{FALSE}}

\item{text}{The text on which the analysis should be performed (when request body is not used)}

\item{field}{Use the analyzer configured for this field (instead of passing the analyzer name)}

\item{analyzer}{The name of the analyzer to use}

\item{tokenizer}{The name of the tokenizer to use for the analysis}

\item{filters}{A character vector of filters to use for the analysis}

\item{char_filters}{A character vector of character filters to use for the analysis}

\item{force}{(logical) Whether a flush should be forced even if it is not necessarily needed
ie. if no changes will be committed to the index.}

\item{full}{(logical) If set to TRUE a new index writer is created and settings that have been
changed related to the index writer will be refreshed.}

\item{wait_if_ongoing}{If TRUE, the flush operation will block until the flush can be executed
if another flush operation is already executing. The default is false and will cause an
exception to be thrown on the shard level if another flush operation is already running.}

\item{filter}{(logical) Clear filter caches}

\item{filter_keys}{(character) A vector of keys to clear when using the \code{filter_cache}
parameter (default: all)}

\item{fielddata}{(logical) Clear field data}

\item{query_cache}{(logical) Clear query caches}

\item{id_cache}{(logical) Clear ID caches for parent/child}

\item{index_new}{(character) an index name, required. only applies to
index_shrink method}
}
\description{
Index API operations
}
\details{
\strong{index_analyze}:
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-analyze.html}
This method can accept a string of text in the body, but this function passes it as a
parameter in a GET request to simplify.

\strong{index_flush}:
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-flush.html}
From the ES website: The flush process of an index basically frees memory from the index by
flushing data to the index storage and clearing the internal transaction log. By default,
Elasticsearch uses memory heuristics in order to automatically trigger flush operations as
required in order to clear memory.

\strong{index_status}: The API endpoint for this function was deprecated in
Elasticsearch \code{v1.2.0}, and will likely be removed soon. Use \code{\link[=index_recovery]{index_recovery()}}
instead.

\strong{index_settings_update}: There are a lot of options you can change with this
function. See
https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-update-settings.html
for all the options.

\strong{index settings}: See
https://www.elastic.co/guide/en/elasticsearch/reference/current/index-modules.html
for the \emph{static} and \emph{dynamic} settings you can set on indices.
}
\section{Mappings}{

The "keyword" type is not supported in Elasticsearch < v5. If you do use a mapping
with "keyword" type in Elasticsearch < v5 \code{\link[=index_create]{index_create()}} should fail.
}

\examples{
\dontrun{
# connection setup
(x <- connect())

# get information on an index
index_get(x, index='shakespeare')
## this one is the same as running index_settings('shakespeare')
index_get(x, index='shakespeare', features='settings')
index_get(x, index='shakespeare', features='mappings')
index_get(x, index='shakespeare', features='alias')

# check for index existence
index_exists(x, index='shakespeare')
index_exists(x, index='plos')

# create an index
if (index_exists(x, 'twitter')) index_delete(x, 'twitter')
index_create(x, index='twitter')
if (index_exists(x, 'things')) index_delete(x, 'things')
index_create(x, index='things')
if (index_exists(x, 'plos')) index_delete(x, 'plos')
index_create(x, index='plos')

# re-create an index
index_recreate(x, "deer")
index_recreate(x, "deer", verbose = FALSE)

# delete an index
if (index_exists(x, 'plos')) index_delete(x, index='plos')

## with a body
body <- '{
 "settings" : {
  "index" : {
    "number_of_shards" : 3,
    "number_of_replicas" : 2
   }
 }
}'
if (index_exists(x, 'alsothat')) index_delete(x, 'alsothat')
index_create(x, index='alsothat', body = body)
## with read only
body <- '{
 "settings" : {
  "index" : {
    "blocks" : {
      "read_only" : true
    }
   }
 }
}'
# index_create(x, index='myindex', body = body)
# then this delete call should fail with something like:
## > Error: 403 - blocked by: [FORBIDDEN/5/index read-only (api)]
# index_delete(x, index='myindex')

## with mappings
body <- '{
 "mappings": {
   "properties": {
     "location" : {"type" : "geo_point"}
   }
 }
}'
if (!index_exists(x, 'gbifnewgeo')) index_create(x, index='gbifnewgeo', body=body)
gbifgeo <- system.file("examples", "gbif_geosmall.json", package = "elastic")
docs_bulk(x, gbifgeo)

# close an index
index_create(x, 'plos')
index_close(x, 'plos')

# open an index
index_open(x, 'plos')

# Get stats on an index
index_stats(x, 'plos')
index_stats(x, c('plos','gbif'))
index_stats(x, c('plos','gbif'), metric='refresh')
index_stats(x, metric = "indexing")
index_stats(x, 'shakespeare', metric='completion')
index_stats(x, 'shakespeare', metric='completion', completion_fields = "completion")
index_stats(x, 'shakespeare', metric='fielddata')
index_stats(x, 'shakespeare', metric='fielddata', fielddata_fields = "evictions")
index_stats(x, 'plos', level="indices")
index_stats(x, 'plos', level="cluster")
index_stats(x, 'plos', level="shards")

# Get segments information that a Lucene index (shard level) is built with
index_segments(x)
index_segments(x, 'plos')
index_segments(x, c('plos','gbif'))

# Get recovery information that provides insight into on-going index shard recoveries
index_recovery(x)
index_recovery(x, 'plos')
index_recovery(x, c('plos','gbif'))
index_recovery(x, "plos", detailed = TRUE)
index_recovery(x, "plos", active_only = TRUE)

# Optimize an index, or many indices
if (x$es_ver() < 500) {
  ### ES < v5 - use optimize
  index_optimize(x, 'plos')
  index_optimize(x, c('plos','gbif'))
  index_optimize(x, 'plos')
} else {
  ### ES > v5 - use forcemerge
  index_forcemerge(x, 'plos')
}

# Upgrade one or more indices to the latest format. The upgrade process converts any
# segments written with previous formats.
if (x$es_ver() < 500) {
  index_upgrade(x, 'plos')
  index_upgrade(x, c('plos','gbif'))
}

# Performs the analysis process on a text and return the tokens breakdown
# of the text
index_analyze(x, text = 'this is a test', analyzer='standard')
index_analyze(x, text = 'this is a test', analyzer='whitespace')
index_analyze(x, text = 'this is a test', analyzer='stop')
index_analyze(x, text = 'this is a test', tokenizer='keyword',
  filters='lowercase')
index_analyze(x, text = 'this is a test', tokenizer='keyword',
  filters='lowercase', char_filters='html_strip')
index_analyze(x, text = 'this is a test', index = 'plos',
  analyzer="standard")
index_analyze(x, text = 'this is a test', index = 'shakespeare',
  analyzer="standard")

## NGram tokenizer
body <- '{
        "settings" : {
             "analysis" : {
                 "analyzer" : {
                     "my_ngram_analyzer" : {
                         "tokenizer" : "my_ngram_tokenizer"
                     }
                 },
                 "tokenizer" : {
                     "my_ngram_tokenizer" : {
                         "type" : "nGram",
                         "min_gram" : "2",
                         "max_gram" : "3",
                         "token_chars": [ "letter", "digit" ]
                     }
                 }
             }
      }
}'
if (index_exists(x, "shakespeare2")) index_delete(x, "shakespeare2")
tokenizer_set(x, index = "shakespeare2", body=body)
index_analyze(x, text = "art thouh", index = "shakespeare2",
  analyzer='my_ngram_analyzer')

# Explicitly flush one or more indices.
index_flush(x, index = "plos")
index_flush(x, index = "shakespeare")
index_flush(x, index = c("plos","shakespeare"))
index_flush(x, index = "plos", wait_if_ongoing = TRUE)
index_flush(x, index = "plos", verbose = TRUE)

# Clear either all caches or specific cached associated with one ore more indices.
index_clear_cache(x)
index_clear_cache(x, index = "plos")
index_clear_cache(x, index = "shakespeare")
index_clear_cache(x, index = c("plos","shakespeare"))
index_clear_cache(x, filter = TRUE)

# Index settings
## get settings
index_settings(x)
index_settings(x, "_all")
index_settings(x, 'gbif')
index_settings(x, c('gbif','plos'))
index_settings(x, '*s')
## update settings
if (index_exists(x, 'foobar')) index_delete(x, 'foobar')
index_create(x, "foobar")
settings <- list(index = list(number_of_replicas = 4))
index_settings_update(x, "foobar", body = settings)
index_get(x, "foobar")$foobar$settings

# Shrink index - Can only shrink an index if it has >1 shard
## index must be read only, a copy of every shard in the index must
## reside on the same node, and the cluster health status must be green
### index_settings_update call to change these
settings <- list(
  index.routing.allocation.require._name = "shrink_node_name",
  index.blocks.write = "true"
)
if (index_exists(x, 'barbarbar')) index_delete(x, 'barbarbar')
index_create(x, "barbarbar")
index_settings_update(x, "barbarbar", body = settings)
cat_recovery(x, index='barbarbar')
# index_shrink(x, "barbarbar", "barfoobbar")
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices.html}
}
\author{
Scott Chamberlain \href{mailto:myrmecocystus@gmail.com}{myrmecocystus@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/index_templates.R
\name{index_template}
\alias{index_template}
\alias{index_template_put}
\alias{index_template_get}
\alias{index_template_exists}
\alias{index_template_delete}
\title{Index templates}
\usage{
index_template_put(
  conn,
  name,
  body = NULL,
  create = NULL,
  flat_settings = NULL,
  master_timeout = NULL,
  order = NULL,
  timeout = NULL,
  ...
)

index_template_get(conn, name = NULL, filter_path = NULL, ...)

index_template_exists(conn, name, ...)

index_template_delete(conn, name, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{name}{(character) The name of the template}

\item{body}{(character/list) The template definition}

\item{create}{(logical) Whether the index template should only be added
if new or can also replace an existing one. Default: \code{FALSE}}

\item{flat_settings}{(logical) Return settings in flat format.
Default: \code{FALSE}}

\item{master_timeout}{(integer) Specify timeout for connection to master}

\item{order}{(integer) The order for this template when merging
multiple matching ones (higher numbers are merged later, overriding the
lower numbers)}

\item{timeout}{(integer) Explicit operation timeout}

\item{...}{Curl options. Or in \code{percolate_list} function, further
args passed on to \code{\link[=Search]{Search()}}}

\item{filter_path}{(character) a regex for filtering output path,
see example}
}
\description{
Index templates allow you to define templates that
will automatically be applied when new indices are created
}
\examples{
\dontrun{
(x <- connect())

body <- '{
  "template": "te*",
  "settings": {
    "number_of_shards": 1
  },
  "mappings": {
    "type1": {
      "_source": {
        "enabled": false
      },
      "properties": {
        "host_name": {
          "type": "keyword"
        },
        "created_at": {
          "type": "date",
          "format": "EEE MMM dd HH:mm:ss Z YYYY"
        }
      }
    }
  }
}'
index_template_put(x, "template_1", body = body)

# get templates
index_template_get(x)
index_template_get(x, "template_1")
index_template_get(x, c("template_1", "template_2"))
index_template_get(x, "template_*")
## filter path
index_template_get(x, "template_1", filter_path = "*.template")

# template exists
index_template_exists(x, "template_1")
index_template_exists(x, "foobar")

# delete a template
index_template_delete(x, "template_1")
index_template_exists(x, "template_1")
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-templates.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{mlt}
\alias{mlt}
\title{This function is defunct}
\usage{
mlt(...)
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_update.R
\name{docs_update}
\alias{docs_update}
\title{Update a document}
\usage{
docs_update(
  conn,
  index,
  id,
  body,
  type = NULL,
  fields = NULL,
  source = NULL,
  version = NULL,
  version_type = NULL,
  routing = NULL,
  parent = NULL,
  timestamp = NULL,
  ttl = NULL,
  refresh = NULL,
  timeout = NULL,
  retry_on_conflict = NULL,
  wait_for_active_shards = NULL,
  detect_noop = NULL,
  callopts = list(),
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) The name of the index. Required}

\item{id}{(numeric/character) The document ID. Can be numeric or character.
Required}

\item{body}{The document, either a list or json}

\item{type}{(character) The type of the document. optional}

\item{fields}{A comma-separated list of fields to return in the response}

\item{source}{Allows to control if and how the updated source should be
returned in the response. By default the updated source is not returned.}

\item{version}{(character) Explicit version number for concurrency control}

\item{version_type}{(character) Specific version type. One of internal,
external, external_gte, or force}

\item{routing}{(character) Specific routing value}

\item{parent}{ID of the parent document. Is is only used for routing and
when for the upsert request}

\item{timestamp}{(date) Explicit timestamp for the document}

\item{ttl}{(aka \dQuote{time to live}) Expiration time for the document.
Expired documents will be expunged automatically. The expiration date that
will be set for a document with a provided ttl is relative to the timestamp
of the document,  meaning it can be based on the time of indexing or on
any time provided. The provided ttl must be strictly positive and can be
a number (in milliseconds) or any valid time value (e.g, 86400000, 1d).}

\item{refresh}{Refresh the index after performing the operation.}

\item{timeout}{(character) Explicit operation timeout, e.g,. 5m (for
5 minutes)}

\item{retry_on_conflict}{Specify how many times should the operation be
retried when a conflict occurs (default: 0)}

\item{wait_for_active_shards}{The number of shard copies required to be
active before proceeding with the update operation.}

\item{detect_noop}{(logical) Specifying \code{TRUE} will cause Elasticsearch
to check if there are changes and, if there aren't, turn the update request
into a noop.}

\item{callopts}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{...}{Further args to query DSL}
}
\description{
Update a document
}
\examples{
\dontrun{
(x <- connect())
if (!index_exists(x, 'plos')) {
  plosdat <- system.file("examples", "plos_data.json",
    package = "elastic")
  plosdat <- type_remover(plosdat)
  invisible(docs_bulk(x, plosdat))
}

docs_create(x, index='plos', id=1002,
  body=list(id="12345", title="New title"))
# and the document is there now
docs_get(x, index='plos', id=1002)
# update the document
docs_update(x, index='plos', id=1002,
  body = list(doc = list(title = "Even newer title again")))
# get it again, notice changes
docs_get(x, index='plos', id=1002)

if (!index_exists(x, 'stuffthings')) {
  index_create(x, "stuffthings")
}
docs_create(x, index='stuffthings', id=1,
  body=list(name = "foo", what = "bar"))
docs_update(x, index='stuffthings', id=1,
  body = list(doc = list(name = "hello", what = "bar")),
  source = 'name')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/type_remover.R
\name{type_remover}
\alias{type_remover}
\title{Utility function to remove 'type' from bulk load files}
\usage{
type_remover(file)
}
\arguments{
\item{file}{(character) a file path, required}
}
\value{
a file path for a temporary file with the types removed
}
\description{
Types are being removed from Elasticsearch. This little function
aims to help remove "_type" fields from bulk newline-delimited JSON
files. See Details.
}
\details{
Looks for any lines that have an "index" key, then drops
any "_type" keys in the hash given by the "index" key.

You can of course manually modify these files as an alternative,
in a text editor or with command line tools like sed, etc.
}
\examples{
\dontrun{
z <- system.file("examples/omdb.json", package = "elastic")
readLines(z, 6)
ff <- type_remover(z)
readLines(ff, 6)
unlink(ff)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ping.r
\name{ping}
\alias{ping}
\title{Ping an Elasticsearch server.}
\usage{
ping(conn, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{...}{Curl args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Ping an Elasticsearch server.
}
\examples{
\dontrun{
x <- connect()
ping(x)
# ideally call ping on the connetion object itself
x$ping()
}
}
\seealso{
\code{\link[=connect]{connect()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mapping.R
\name{mapping}
\alias{mapping}
\alias{mapping_create}
\alias{mapping_get}
\alias{field_mapping_get}
\alias{type_exists}
\title{Mapping management}
\usage{
mapping_create(
  conn,
  index,
  body,
  type = NULL,
  update_all_types = FALSE,
  include_type_name = NULL,
  ...
)

mapping_get(conn, index = NULL, type = NULL, include_type_name = NULL, ...)

field_mapping_get(
  conn,
  index = NULL,
  type = NULL,
  field,
  include_defaults = FALSE,
  include_type_name = NULL,
  ...
)

type_exists(conn, index, type, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) An index}

\item{body}{(list) Either a list or json, representing the query.}

\item{type}{(character) A document type}

\item{update_all_types}{(logical) update all types. default: \code{FALSE}.
This parameter is deprecated in ES v6.3.0 and higher, see
https://github.com/elastic/elasticsearch/pull/28284}

\item{include_type_name}{(logical) If set to \code{TRUE}, you can include a type
name, if not an error will occur. default: not set. See Details.}

\item{...}{Curl options passed on to \link[crul:verb-PUT]{crul::verb-PUT}, \link[crul:verb-GET]{crul::verb-GET},
or \link[crul:verb-HEAD]{crul::verb-HEAD}}

\item{field}{(character) One or more field names}

\item{include_defaults}{(logical) Whether to return default values}
}
\description{
Mapping management
}
\details{
Find documentation for each function at:
\itemize{
\item \code{mapping_create} -
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-put-mapping.html}
\item \code{type_exists} -
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-types-exists.html}
\item \code{mapping_delete} - FUNCTION DEFUNCT - instead of deleting mapping, delete
index and recreate index with new mapping
\item \code{mapping_get} -
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-get-mapping.html}
\item \code{field_mapping_get} -
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-get-field-mapping.html}
}

See \url{https://www.elastic.co/guide/en/elasticsearch/reference/current/removal-of-types.html}
for information on type removal
}
\examples{
\dontrun{
# connection setup
(x <- connect())

# Used to check if a type/types exists in an index/indices
type_exists(x, index = "plos", type = "article")
type_exists(x, index = "plos", type = "articles")
type_exists(x, index = "shakespeare", type = "line")

# The put mapping API allows to register specific mapping definition for a specific type.
## a good mapping body
body <- list(properties = list(
 journal = list(type="text"),
 year = list(type="long")
))
if (!index_exists(x, "plos")) index_create(x, "plos")
mapping_create(x, index = "plos", type = "citation", body=body)
## OR if above fails, try
mapping_create(x, index = "plos", type = "citation", body=body,
  include_type_name=TRUE)
## ES >= 7, no type
mapping_create(x, index = "plos", body=body)

### or as json
body <- '{
  "properties": {
    "journal": { "type": "text" },
      "year": { "type": "long" }
}}'
mapping_create(x, index = "plos", type = "citation", body=body)
mapping_get(x, "plos", "citation")

## A bad mapping body
body <- list(things = list(properties = list(
  journal = list("text")
)))
# mapping_create(x, index = "plos", type = "things", body=body)

# Get mappings
mapping_get(x, '_all')
mapping_get(x, index = "plos")
mapping_get(x, index = c("shakespeare","plos"))
# mapping_get(x, index = "shakespeare", type = "act")
# mapping_get(x, index = "shakespeare", type = c("act","line"))

# Get field mappings
plosdat <- system.file("examples", "plos_data.json",
  package = "elastic")
plosdat <- type_remover(plosdat)
invisible(docs_bulk(x, plosdat))
field_mapping_get(x, index = "_all", field = "text")
field_mapping_get(x, index = "plos", field = "title")
field_mapping_get(x, index = "plos", field = "*")
field_mapping_get(x, index = "plos", field = "title", include_defaults = TRUE)
field_mapping_get(x, type = c("article","record"), field = c("title","class"))
field_mapping_get(x, type = "a*", field = "t*")

# Create geospatial mapping
if (index_exists(x, "gbifgeopoint")) index_delete(x, "gbifgeopoint")
file <- system.file("examples", "gbif_geopoint.json",
  package = "elastic")
file <- type_remover(file)
index_create(x, "gbifgeopoint")
body <- '{
 "properties" : {
   "location" : { "type" : "geo_point" }
 }
}'
mapping_create(x, "gbifgeopoint", body = body)
invisible(docs_bulk(x, file))

# update_all_fields, see also ?fielddata
if (x$es_ver() < 603) {
 mapping_create(x, "shakespeare", "record", update_all_types=TRUE, body = '{
   "properties": {
     "speaker": { 
       "type":     "text",
       "fielddata": true
     }
   }
 }')
} else {
 index_create(x, 'brownchair')
 mapping_create(x, 'brownchair', body = '{
   "properties": {
     "foo": { 
       "type":     "text",
       "fielddata": true
     }
   }
 }')
}

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/validate.R
\name{validate}
\alias{validate}
\title{Validate a search}
\usage{
validate(conn, index, type = NULL, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{Index name. Required.}

\item{type}{Document type. Optional.}

\item{...}{Additional args passed on to \code{\link[=Search]{Search()}}}
}
\description{
Validate a search
}
\examples{
\dontrun{
x <- connect()

if (!index_exists(x, "twitter")) index_create(x, "twitter")
docs_create(x, 'twitter', id=1, body = list(
   "user" = "foobar", 
   "post_date" = "2014-01-03",
   "message" = "trying out Elasticsearch"
 )
)
validate(x, "twitter", q='user:foobar')
validate(x, "twitter", q='user:foobar')

body <- '{
"query" : {
  "bool" : {
    "must" : {
      "query_string" : {
        "query" : "*:*"
      }
    },
    "filter" : {
      "term" : { "user" : "kimchy" }
    }
  }
}
}'
validate(x, "twitter", body = body)
}
}
\seealso{
\code{\link[=Search]{Search()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search.r
\name{Search}
\alias{Search}
\title{Full text search of Elasticsearch}
\usage{
Search(
  conn,
  index = NULL,
  type = NULL,
  q = NULL,
  df = NULL,
  analyzer = NULL,
  default_operator = NULL,
  explain = NULL,
  source = NULL,
  fields = NULL,
  sort = NULL,
  track_scores = NULL,
  timeout = NULL,
  terminate_after = NULL,
  from = NULL,
  size = NULL,
  search_type = NULL,
  lowercase_expanded_terms = NULL,
  analyze_wildcard = NULL,
  version = NULL,
  lenient = NULL,
  body = list(),
  raw = FALSE,
  asdf = FALSE,
  track_total_hits = TRUE,
  time_scroll = NULL,
  search_path = "_search",
  stream_opts = list(),
  ignore_unavailable = FALSE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link{connect}}}

\item{index}{Index name, one or more}

\item{type}{Document type. Note that \code{type} is deprecated in
Elasticsearch v7 and greater, and removed in Elasticsearch v8. We will
strive to support types for folks using older ES versions}

\item{q}{The query string (maps to the query_string query, see Query String
Query for more details). See
https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-query-string-query.html
for documentation and examples.}

\item{df}{(character) The default field to use when no field prefix is
defined within the query.}

\item{analyzer}{(character) The analyzer name to be used when analyzing the
query string.}

\item{default_operator}{(character) The default operator to be used, can be
\code{AND} or \code{OR}. Default: \code{OR}}

\item{explain}{(logical) For each hit, contain an explanation of how
scoring of the hits was computed. Default: \code{FALSE}}

\item{source}{(logical) Set to \code{FALSE} to disable retrieval of the
\code{_source} field. You can also retrieve part of the document by
using \code{_source_include} & \code{_source_exclude} (see the \code{body}
documentation for more details). You can also include a comma-delimited
string of fields from the source document that you want back. See also
the \strong{fields} parameter}

\item{fields}{(character) The selective stored fields of the document to
return for each hit. Not specifying any value will cause no fields to return.
Note that in Elasticsearch v5 and greater, \strong{fields} parameter has
changed to \strong{stored_fields}, which is not on by default. You can
however, pass fields to \strong{source} parameter}

\item{sort}{(character) Sorting to perform. Can either be in the form of
fieldName, or \code{fieldName:asc}/\code{fieldName:desc}. The fieldName
can either be an actual field within the document, or the special
\code{_score} name to indicate sorting based on scores. There can be several
sort parameters (order is important).}

\item{track_scores}{(logical) When sorting, set to \code{TRUE} in order to
still track scores and return them as part of each hit.}

\item{timeout}{(numeric) A search timeout, bounding the search request to
be executed within the specified time value and bail with the hits
accumulated up to that point when expired. Default: no timeout.}

\item{terminate_after}{(numeric) The maximum number of documents to collect
for each shard, upon reaching which the query execution will terminate
early. If set, the response will have a boolean field terminated_early to
indicate whether the query execution has actually terminated_early.
Default: no terminate_after}

\item{from}{(character) The starting from index of the hits to return.
Pass in as a character string to avoid problems with large number
conversion to scientific notation. Default: 0}

\item{size}{(character) The number of hits to return. Pass in as a
character string to avoid problems with large number conversion to
scientific notation. Default: 10. The default maximum is 10,000 - however,
you can change this default maximum by changing the
\code{index.max_result_window} index level parameter.}

\item{search_type}{(character) The type of the search operation to perform.
Can be \code{query_then_fetch} (default) or \code{dfs_query_then_fetch}.
Types \code{scan} and \code{count} are deprecated.
See Elasticsearch docs for more details on the different types of
search that can be performed.}

\item{lowercase_expanded_terms}{(logical) Should terms be automatically
lowercased or not. Default: \code{TRUE}.}

\item{analyze_wildcard}{(logical) Should wildcard and prefix queries be
analyzed or not. Default: \code{FALSE}.}

\item{version}{(logical) Print the document version with each document.}

\item{lenient}{(logical) If \code{TRUE} will cause format based failures (like
providing text to a numeric field) to be ignored. Default: \code{NULL}}

\item{body}{Query, either a list or json.}

\item{raw}{(logical) If \code{FALSE} (default), data is parsed to list.
If \code{TRUE}, then raw JSON returned}

\item{asdf}{(logical) If \code{TRUE}, use \code{\link[jsonlite]{fromJSON}}
to parse JSON directly to a data.frame. If \code{FALSE} (Default), list
output is given.}

\item{track_total_hits}{(logical, numeric) If \code{TRUE} will always track
the number of hits that match the query accurately. If \code{FALSE} will
count documents accurately up to 10000 documents. If \code{is.integer} will
count documents accurately up to the number. Default: \code{TRUE}}

\item{time_scroll}{(character) Specify how long a consistent view of the
index should be maintained for scrolled search, e.g., "30s", "1m". See
\link{units-time}}

\item{search_path}{(character) The path to use for searching. Default
to \verb{_search}, but in some cases you may already have that in the base
url set using \code{\link[=connect]{connect()}}, in which case you can set this
to \code{NULL}}

\item{stream_opts}{(list) A list of options passed to
\code{\link[jsonlite]{stream_out}} - Except that you can't pass \code{x} as
that's the data that's streamed out, and pass a file path instead of a
connection to \code{con}. \code{pagesize} param doesn't do much as
that's more or less controlled by paging with ES.}

\item{ignore_unavailable}{(logical) What to do if an specified index name
doesn't exist. If set to \code{TRUE} then those indices are ignored.}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}}}
}
\description{
Full text search of Elasticsearch
}
\details{
This function name has the "S" capitalized to avoid conflict with the function
\code{base::search}. I hate mixing cases, as I think it confuses users, but in this case
it seems neccessary.
}
\section{profile}{

The Profile API provides detailed timing information about the execution of
individual components in a search request. See
https://www.elastic.co/guide/en/elasticsearch/reference/current/search-profile.html
for more information

In a body query, you can set to \code{profile: true} to enable profiling
results. e.g.

\preformatted{
{
  "profile": true,
  "query" : {
    "match" : { "message" : "some number" }
  }
}
}
}

\examples{
\dontrun{
# make connection object
(x <- connect())

# load some data
if (!index_exists(x, "shakespeare")) {
  shakespeare <- system.file("examples", "shakespeare_data.json",
    package = "elastic")
  shakespeare <- type_remover(shakespeare)
  invisible(docs_bulk(x, shakespeare))
}
if (!index_exists(x, "gbif")) {
  gbif <- system.file("examples", "gbif_data.json",
    package = "elastic")
  gbif <- type_remover(gbif)
  invisible(docs_bulk(x, gbif))
}
if (!index_exists(x, "plos")) {
  plos <- system.file("examples", "plos_data.json",
    package = "elastic")
  plos <- type_remover(plos)
  invisible(docs_bulk(x, plos))
}


# URI string queries
Search(x, index="shakespeare")
## if you're using an older ES version, you may have types
if (gsub("\\\\.", "", x$ping()$version$number) < 700) {
  Search(x, index="shakespeare", type="act")
  Search(x, index="shakespeare", type="scene")
  Search(x, index="shakespeare", type="line")
}

## Return certain fields
if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
  ### ES < v5
  Search(x, index="shakespeare", fields=c('play_name','speaker'))
} else {
  ### ES > v5
  Search(x, index="shakespeare", body = '{
   "_source": ["play_name", "speaker"]
  }')
}

## Search multiple indices
Search(x, index = "gbif")$hits$total$value
Search(x, index = "shakespeare")$hits$total$value
Search(x, index = c("gbif", "shakespeare"))$hits$total$value

## search_type
Search(x, index="shakespeare", search_type = "query_then_fetch")
Search(x, index="shakespeare", search_type = "dfs_query_then_fetch")
### search type "scan" is gone - use time_scroll instead
Search(x, index="shakespeare", time_scroll = "2m")
### search type "count" is gone - use size=0 instead
Search(x, index="shakespeare", size = 0)$hits$total$value

## search exists check
### use size set to 0 and terminate_after set to 1
### if there are > 0 hits, then there are matching documents
Search(x, index="shakespeare", size = 0, terminate_after = 1)

## sorting
### if ES >5, we need to make sure fielddata is turned on for a field 
### before using it for sort 
if (gsub("\\\\.", "", x$ping()$version$number) >= 500) {
 if (index_exists(x, "shakespeare")) index_delete(x, "shakespeare")
 index_create(x, "shakespeare")
 mapping_create(x, "shakespeare", body = '{
    "properties": {
      "speaker": { 
        "type":     "text",
        "fielddata": true
      }
    }
  }'
 )
 shakespeare <- system.file("examples", "shakespeare_data.json",
   package = "elastic")
 shakespeare <- type_remover(shakespeare)
 invisible(docs_bulk(x, shakespeare))
 z <- Search(x, index="shakespeare", sort="speaker", size = 30)
 vapply(z$hits$hits, function(w) w$`_source`$speaker, "")
}

if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
  Search(x, index="shakespeare", type="line", sort="speaker:desc", 
    fields='speaker')
  Search(x, index="shakespeare", type="line",
    sort=c("speaker:desc","play_name:asc"), fields=c('speaker','play_name'))
}


## pagination
Search(x, index="shakespeare", size=1)$hits$hits
Search(x, index="shakespeare", size=1, from=1)$hits$hits

## queries
### Search in all fields
Search(x, index="shakespeare", q="york")

### Searchin specific fields
Search(x, index="shakespeare", q="speaker:KING HENRY IV")$hits$total$value

### Exact phrase search by wrapping in quotes
Search(x, index="shakespeare", q='speaker:"KING HENRY IV"')$hits$total$value

### can specify operators between multiple words parenthetically
Search(x, index="shakespeare", q="speaker:(HENRY OR ARCHBISHOP)")$hits$total$value

### where the field line_number has no value (or is missing)
Search(x, index="shakespeare", q="_missing_:line_number")$hits$total$value

### where the field line_number has any non-null value
Search(x, index="shakespeare", q="_exists_:line_number")$hits$total$value

### wildcards, either * or ?
Search(x, index="shakespeare", q="*ay")$hits$total$value
Search(x, index="shakespeare", q="m?y")$hits$total$value

### regular expressions, wrapped in forward slashes
Search(x, index="shakespeare", q="text_entry:/[a-z]/")$hits$total$value

### fuzziness
Search(x, index="shakespeare", q="text_entry:ma~")$hits$total$value
Search(x, index="shakespeare", q="text_entry:the~2")$hits$total$value
Search(x, index="shakespeare", q="text_entry:the~1")$hits$total$value

### Proximity searches
Search(x, index="shakespeare", q='text_entry:"as hath"~5')$hits$total$value
Search(x, index="shakespeare", q='text_entry:"as hath"~10')$hits$total$value

### Ranges, here where line_id value is between 10 and 20
Search(x, index="shakespeare", q="line_id:[10 TO 20]")$hits$total$value

### Grouping
Search(x, index="shakespeare", q="(hath OR as) AND the")$hits$total$value

# Limit number of hits returned with the size parameter
Search(x, index="shakespeare", size=1)

# Give explanation of search in result
Search(x, index="shakespeare", size=1, explain=TRUE)

## terminate query after x documents found
## setting to 1 gives back one document for each shard
Search(x, index="shakespeare", terminate_after=1)
## or set to other number
Search(x, index="shakespeare", terminate_after=2)

## Get version number for each document
Search(x, index="shakespeare", version=TRUE, size=2)

## Get raw data
Search(x, index="shakespeare", raw = TRUE)

## Curl options 
### verbose 
out <- Search(x, index="shakespeare", verbose = TRUE)


# Query DSL searches - queries sent in the body of the request
## Pass in as an R list

### if ES >5, we need to make sure fielddata is turned on for a field 
### before using it for aggregations 
if (gsub("\\\\.", "", x$ping()$version$number) >= 500) {
  mapping_create(x, "shakespeare", update_all_types = TRUE, body = '{
    "properties": {
      "text_entry": { 
        "type":     "text",
        "fielddata": true
     }
   }
 }')
 aggs <- list(aggs = list(stats = list(terms = list(field = "text_entry"))))
 Search(x, index="shakespeare", body=aggs)
}

### if ES >5, you don't need to worry about fielddata
if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
   aggs <- list(aggs = list(stats = list(terms = list(field = "text_entry"))))
   Search(x, index="shakespeare", body=aggs)
}

## or pass in as json query with newlines, easy to read
aggs <- '{
    "aggs": {
        "stats" : {
            "terms" : {
                "field" : "speaker"
            }
        }
    }
}'
Search(x, index="shakespeare", body=aggs, asdf=TRUE, size = 0)

## or pass in collapsed json string
aggs <- '{"aggs":{"stats":{"terms":{"field":"text_entry"}}}}'
Search(x, index="shakespeare", body=aggs)


## Aggregations
### Histograms
aggs <- '{
    "aggs": {
        "latbuckets" : {
           "histogram" : {
               "field" : "decimalLatitude",
               "interval" : 5
           }
        }
    }
}'
Search(x, index="gbif", body=aggs, size=0)

### Histograms w/ more options
aggs <- '{
    "aggs": {
        "latbuckets" : {
           "histogram" : {
               "field" : "decimalLatitude",
               "interval" : 5,
               "min_doc_count" : 0,
               "extended_bounds" : {
                   "min" : -90,
                   "max" : 90
               }
           }
        }
    }
}'
Search(x, index="gbif", body=aggs, size=0)

### Ordering the buckets by their doc_count - ascending:
aggs <- '{
    "aggs": {
        "latbuckets" : {
           "histogram" : {
               "field" : "decimalLatitude",
               "interval" : 5,
               "min_doc_count" : 0,
               "extended_bounds" : {
                   "min" : -90,
                   "max" : 90
               },
               "order" : {
                   "_count" : "desc"
               }
           }
        }
    }
}'
out <- Search(x, index="gbif", body=aggs, size=0)
lapply(out$aggregations$latbuckets$buckets, data.frame)

### By default, the buckets are returned as an ordered array. It is also possible to
### request the response as a hash instead keyed by the buckets keys:
aggs <- '{
    "aggs": {
        "latbuckets" : {
           "histogram" : {
               "field" : "decimalLatitude",
               "interval" : 10,
               "keyed" : true
           }
        }
    }
}'
Search(x, index="gbif", body=aggs, size=0)

# match query
match <- '{"query": {"match" : {"text_entry" : "Two Gentlemen"}}}'
Search(x, index="shakespeare", body=match)

# multi-match (multiple fields that is) query
mmatch <- '{"query": {"multi_match" : {"query" : "henry", "fields": ["text_entry","play_name"]}}}'
Search(x, index="shakespeare", body=mmatch)

# bool query
mmatch <- '{
 "query": {
   "bool" : {
     "must_not" : {
       "range" : {
         "speech_number" : {
           "from" : 1, "to": 5
}}}}}}'
Search(x, index="shakespeare", body=mmatch)

# Boosting query
boost <- '{
 "query" : {
  "boosting" : {
      "positive" : {
          "term" : {
              "play_name" : "henry"
          }
      },
      "negative" : {
          "term" : {
              "text_entry" : "thou"
          }
      },
      "negative_boost" : 0.8
    }
 }
}'
Search(x, index="shakespeare", body=boost)

# Fuzzy query
## fuzzy query on numerics
fuzzy <- list(query = list(fuzzy = list(text_entry = "arms")))
Search(x, index="shakespeare", body=fuzzy)$hits$total$value
fuzzy <- list(query = list(fuzzy = list(text_entry = list(value = "arms", fuzziness = 4))))
Search(x, index="shakespeare", body=fuzzy)$hits$total$value

# geoshape query
## not working yets
geo <- list(query = list(geo_shape = list(location = list(shape = list(type = "envelope",
   coordinates = "[[2,10],[10,20]]")))))
geo <- '{
 "query": {
   "geo_shape": {
     "location": {
       "point": {
         "type": "envelope",
         "coordinates": [[2,0],[2.93,100]]
       }
     }
   }
 }
}'
# Search(x, index="gbifnewgeo", body=geo)

# range query
## with numeric
body <- list(query=list(range=list(decimalLongitude=list(gte=1, lte=3))))
Search(x, 'gbif', body=body)$hits$total$value

body <- list(query=list(range=list(decimalLongitude=list(gte=2.9, lte=10))))
Search(x, 'gbif', body=body)$hits$total$value

## with dates
body <- list(query=list(range=list(eventDate=list(gte="2012-01-01", lte="now"))))
Search(x, 'gbif', body=body)$hits$total$value

body <- list(query=list(range=list(eventDate=list(gte="2014-01-01", lte="now"))))
Search(x, 'gbif', body=body)$hits$total$value

# more like this query (more_like_this can be shortened to mlt)
body <- '{
 "query": {
   "more_like_this": {
     "fields": ["title"],
     "like": "and then",
     "min_term_freq": 1,
     "max_query_terms": 12
   }
 }
}'
Search(x, 'plos', body=body)$hits$total$value

body <- '{
 "query": {
   "more_like_this": {
     "fields": ["abstract","title"],
     "like": "cell",
     "min_term_freq": 1,
     "max_query_terms": 12
   }
 }
}'
Search(x, 'plos', body=body)$hits$total$value

# Highlighting
body <- '{
 "query": {
   "query_string": {
     "query" : "cell"
   }
 },
 "highlight": {
   "fields": {
     "title": {"number_of_fragments": 2}
   }
 }
}'
out <- Search(x, 'plos', body=body)
out$hits$total$value
sapply(out$hits$hits, function(x) x$`_source`$title[[1]])

### Common terms query
body <- '{
 "query" : {
   "match": {
      "text_entry": {
         "query": "this is"
      }
   }
 }
}'
Search(x, 'shakespeare', body=body)

## Scrolling search - instead of paging
res <- Search(x, index = 'shakespeare', q="a*", time_scroll="1m")
scroll(x, res$`_scroll_id`)

res <- Search(x, index = 'shakespeare', q="a*", time_scroll="5m")
out <- list()
hits <- 1
while(hits != 0){
  res <- scroll(x, res$`_scroll_id`)
  hits <- length(res$hits$hits)
  if(hits > 0)
    out <- c(out, res$hits$hits)
}

### Sliced scrolling
#### For scroll queries that return a lot of documents it is possible to 
#### split the scroll in multiple slices which can be consumed independently
body1 <- '{
  "slice": {
    "id": 0, 
    "max": 2 
  },
  "query": {
    "match" : {
      "text_entry" : "a*"
    }
  }
}'

body2 <- '{
  "slice": {
    "id": 1, 
    "max": 2 
  },
  "query": {
    "match" : {
      "text_entry" : "a*"
    }
  }
}'

res1 <- Search(x, index = 'shakespeare', time_scroll="1m", body = body1)
res2 <- Search(x, index = 'shakespeare', time_scroll="1m", body = body2)
scroll(x, res1$`_scroll_id`)
scroll(x, res2$`_scroll_id`)

out1 <- list()
hits <- 1
while(hits != 0){
  tmp1 <- scroll(x, res1$`_scroll_id`)
  hits <- length(tmp1$hits$hits)
  if(hits > 0)
    out1 <- c(out1, tmp1$hits$hits)
}

out2 <- list()
hits <- 1
while(hits != 0) {
  tmp2 <- scroll(x, res2$`_scroll_id`)
  hits <- length(tmp2$hits$hits)
  if(hits > 0)
    out2 <- c(out2, tmp2$hits$hits)
}

c(
 lapply(out1, "[[", "_source"),
 lapply(out2, "[[", "_source")
) 



# Using filters
## A bool filter
body <- '{
 "query":{
   "bool": {
     "must_not" : {
       "range" : {
         "year" : { "from" : 2011, "to" : 2012 }
       }
     }
   }
 }
}'
Search(x, 'gbif', body = body)$hits$total$value

## Geo filters - fun!
### Note that filers have many geospatial filter options, but queries 
### have fewer, andrequire a geo_shape mapping

body <- '{
 "mappings": {
     "properties": {
         "location" : {"type" : "geo_point"}
      }
   }
}'
index_recreate(x, index='gbifgeopoint', body=body)
path <- system.file("examples", "gbif_geopoint.json",
  package = "elastic")
path <- type_remover(path)
invisible(docs_bulk(x, path))

### Points within a bounding box
body <- '{
 "query":{
   "bool" : {
     "must" : {
       "match_all" : {}
     },
     "filter":{
        "geo_bounding_box" : {
          "location" : {
            "top_left" : {
              "lat" : 60,
              "lon" : 1
            },
            "bottom_right" : {
              "lat" : 40,
              "lon" : 14
            }
          }
       }
     }
   }
 }
}'
out <- Search(x, 'gbifgeopoint', body = body, size = 300)
out$hits$total$value
do.call(rbind, lapply(out$hits$hits, function(x) x$`_source`$location))

### Points within distance of a point
body <- '{
"query": {
  "bool" : {
    "must" : {
      "match_all" : {}
    },
   "filter" : {
     "geo_distance" : {
       "distance" : "200km",
       "location" : {
         "lon" : 4,
         "lat" : 50
       }
     }
  }
}}}'
out <- Search(x, 'gbifgeopoint', body = body)
out$hits$total$value
do.call(rbind, lapply(out$hits$hits, function(x) x$`_source`$location))

### Points within distance range of a point
body <- '{
 "aggs":{
   "points_within_dist" : {
     "geo_distance" : {
        "field": "location",
        "origin" : "4, 50",
        "ranges": [ 
          {"from" : 200},
          {"to" : 400}
         ]
     }
   }
 }
}'
out <- Search(x, 'gbifgeopoint', body = body)
out$hits$total$value
do.call(rbind, lapply(out$hits$hits, function(x) x$`_source`$location))

### Points within a polygon
body <- '{
 "query":{
   "bool" : {
     "must" : {
       "match_all" : {}
     },
     "filter":{
        "geo_polygon" : {
          "location" : {
             "points" : [
               [80.0, -20.0], [-80.0, -20.0], [-80.0, 60.0], [40.0, 60.0], [80.0, -20.0]
             ]
           }
         }
       }
     }
   }
}'
out <- Search(x, 'gbifgeopoint', body = body)
out$hits$total$value
do.call(rbind, lapply(out$hits$hits, function(x) x$`_source`$location))

### Geoshape filters using queries instead of filters
#### Get data with geojson type location data loaded first
body <- '{
 "mappings": {
     "properties": {
         "location" : {"type" : "geo_shape"}
      }
   }
}'
index_recreate(x, index='geoshape', body=body)
path <- system.file("examples", "gbif_geoshape.json",
  package = "elastic")
path <- type_remover(path)
invisible(docs_bulk(x, path))

#### Get data with a square envelope, w/ point defining upper left and the other
#### defining the lower right
body <- '{
 "query":{
   "geo_shape" : {
     "location" : {
         "shape" : {
           "type": "envelope",
            "coordinates": [[-30, 50],[30, 0]]
         }
       }
     }
   }
}'
out <- Search(x, 'geoshape', body = body)
out$hits$total$value

#### Get data with a circle, w/ point defining center, and radius
body <- '{
 "query":{
   "geo_shape" : {
     "location" : {
         "shape" : {
           "type": "circle",
           "coordinates": [-10, 45],
           "radius": "2000km"
         }
       }
     }
   }
}'
out <- Search(x, 'geoshape', body = body)
out$hits$total$value

#### Use a polygon, w/ point defining center, and radius
body <- '{
 "query":{
   "geo_shape" : {
     "location" : {
         "shape" : {
           "type": "polygon",
           "coordinates":  [
              [ [80.0, -20.0], [-80.0, -20.0], [-80.0, 60.0], [40.0, 60.0], [80.0, -20.0] ]
           ]
         }
       }
     }
   }
}'
out <- Search(x, 'geoshape', body = body)
out$hits$total$value


# Geofilter with WKT
# format follows "BBOX (minlon, maxlon, maxlat, minlat)"
body <- '{
    "query": {
        "bool" : {
            "must" : {
                "match_all" : {}
            },
            "filter" : {
                "geo_bounding_box" : {
                    "location" : {
                        "wkt" : "BBOX (1, 14, 60, 40)"
                    }
                }
            }
        }
    }
}'
out <- Search(x, 'gbifgeopoint', body = body)
out$hits$total$value



# Missing filter
if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
  ### ES < v5
  body <- '{
   "query":{
     "constant_score" : {
       "filter" : {
         "missing" : { "field" : "play_name" }
       }
     }
    }
  }'
  Search(x, "shakespeare", body = body)
} else {
  ### ES => v5
  body <- '{
   "query":{
     "bool" : {
       "must_not" : {
         "exists" : { 
           "field" : "play_name" 
         }
       }
    }
   }
  }'
  Search(x, "shakespeare", body = body)
}

# prefix filter
body <- '{
 "query": {
   "bool": {
     "must": {
       "prefix" : {
         "speaker" : "we"
       }
     }
   }
 }
}'
z <- Search(x, "shakespeare", body = body)
z$hits$total$value
vapply(z$hits$hits, "[[", "", c("_source", "speaker"))


# ids filter
if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
  ### ES < v5
  body <- '{
   "query":{
     "bool": {
       "must": {
         "ids" : {
           "values": ["1","2","10","2000"]
        }
      }
    }
   }
  }'
  z <- Search(x, "shakespeare", body = body)
  z$hits$total$value
  identical(
   c("1","2","10","2000"),
   vapply(z$hits$hits, "[[", "", "_id")
  )
} else {
  body <- '{
   "query":{
     "ids" : {
       "values": ["1","2","10","2000"]
     }
   }
  }'
  z <- Search(x, "shakespeare", body = body)
  z$hits$total$value
  identical(
   c("1","2","10","2000"),
   vapply(z$hits$hits, "[[", "", "_id")
  )
}

# combined prefix and ids filters
if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
  ### ES < v5
  body <- '{
   "query":{
     "bool" : {
       "should" : {
         "or": [{
           "ids" : {
             "values": ["1","2","3","10","2000"]
           }
         }, {
         "prefix" : {
           "speaker" : "we"
         }
        }
      ]
     }
    }
   }
  }'
  z <- Search(x, "shakespeare", body = body)
  z$hits$total$value
} else {
  ### ES => v5
  body <- '{
   "query":{
     "bool" : {
       "should" : [
         {
           "ids" : {
             "values": ["1","2","3","10","2000"]
           }
         }, 
         {
           "prefix" : {
             "speaker" : "we"
           }
         }
      ]
     }
    }
  }'
  z <- Search(x, "shakespeare", body = body)
  z$hits$total$value
}

# Suggestions
sugg <- '{
 "query" : {
    "match" : {
      "text_entry" : "late"
     }
 },  
 "suggest" : {
   "sugg" : {
     "text" : "late",
     "term" : {
         "field" : "text_entry"
      }
    }
  }
}'
Search(x, index = "shakespeare", body = sugg, 
  asdf = TRUE, size = 0)$suggest$sugg$options



# stream data out using jsonlite::stream_out
file <- tempfile()
res <- Search(x, "shakespeare", size = 1000, stream_opts = list(file = file))
head(df <- jsonlite::stream_in(file(file)))
NROW(df)
unlink(file)


# get profile data
body <- '{
  "profile": true,
  "query" : {
    "match" : { "text_entry" : "war" }
  }
}'
res <- Search(x, "shakespeare", body = body)
res$profile
# time in nanoseconds across each of the shards
vapply(res$profile$shards, function(w) {
  w$searches[[1]]$query[[1]]$time_in_nanos
}, 1)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/search-search.html}
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl.html}
}
\seealso{
\code{\link[=Search_uri]{Search_uri()}} \code{\link[=Search_template]{Search_template()}} \code{\link[=scroll]{scroll()}} \code{\link[=count]{count()}}
\code{\link[=validate]{validate()}} \code{\link[=fielddata]{fielddata()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/scroll.R
\name{scroll}
\alias{scroll}
\alias{scroll_clear}
\title{Scroll search function}
\usage{
scroll(
  conn,
  x,
  time_scroll = "1m",
  raw = FALSE,
  asdf = FALSE,
  stream_opts = list(),
  ...
)

scroll_clear(conn, x = NULL, all = FALSE, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{x}{(character) For \code{scroll}, a single scroll id; for
\code{scroll_clear}, one or more scroll id's}

\item{time_scroll}{(character) Specify how long a consistent view of the
index should be maintained for scrolled search, e.g., "30s", "1m".
See \link{units-time}.}

\item{raw}{(logical) If \code{FALSE} (default), data is parsed to list.
If \code{TRUE}, then raw JSON.}

\item{asdf}{(logical) If \code{TRUE}, use \code{\link[jsonlite:fromJSON]{jsonlite::fromJSON()}}
to parse JSON directly to a data.frame. If \code{FALSE} (Default), list
output is given.}

\item{stream_opts}{(list) A list of options passed to
\code{\link[jsonlite:stream_in]{jsonlite::stream_out()}} - Except that you can't pass \code{x} as
that's the data that's streamed out, and pass a file path sinstead of a
connection to \code{con}. \code{pagesize} param doesn't do much as
that's more or less controlled by paging with ES.}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}

\item{all}{(logical) If \code{TRUE} (default) then all search contexts
cleared.  If \code{FALSE}, scroll id's must be passed to \code{x}}
}
\value{
\code{scroll()} returns a list, identical to what
\code{\link[=Search]{Search()}} returns. With attribute \code{scroll} that is the
scroll value set via the \code{time_scroll} parameter

\code{scroll_clear()} returns a boolean (\code{TRUE} on success)
}
\description{
Scroll search function
}
\section{Scores}{

Scores will be the same for all documents that are returned from a
scroll request. Dems da rules.
}

\section{Inputs}{

Inputs to \code{scroll()} can be one of:
\itemize{
\item list - This usually will be the output of \code{\link[=Search]{Search()}}, but
you could in theory make a list yourself with the appropriate elements
\item character - A scroll ID - this is typically the scroll id output
from a call to \code{\link[=Search]{Search()}}, accessed like \code{res$`_scroll_id`}
}

All other classes passed to \code{scroll()} will fail with message

Lists passed to \code{scroll()} without a \verb{_scroll_id} element will
trigger an error.

From lists output form \code{\link[=Search]{Search()}} there should be an attribute
("scroll") that is the \code{scroll} value set in the \code{\link[=Search]{Search()}}
request - if that attribute is missing from the list, we'll attempt to
use the \code{time_scroll} parameter value set in the
\code{scroll()} function call

The output of \code{scroll()} has the scroll time value as an attribute so
the output can be passed back into \code{scroll()} to continue.
}

\section{Clear scroll}{

Search context are automatically removed when the scroll timeout has
been exceeded.  Keeping scrolls open has a cost, so scrolls should be
explicitly cleared as soon  as the scroll is not being used anymore
using \code{scroll_clear}
}

\section{Sliced scrolling}{

For scroll queries that return a lot of documents it is possible to split
the scroll in multiple slices which can be consumed independently.

See the example in this man file.
}

\section{Aggregations}{

If the request specifies aggregations, only the initial search response
will contain the aggregations results.
}

\examples{
\dontrun{
# connection setup
(con <- connect())

# Basic usage - can use across all indices
res <- Search(con, time_scroll="1m")
scroll(con, res)$`_scroll_id`

# use on a specific index - and specify a query
res <- Search(con, index = 'shakespeare', q="a*", time_scroll="1m")
res$`_scroll_id`

# Setting "sort=_doc" to turn off sorting of results - faster
res <- Search(con, index = 'shakespeare', q="a*", time_scroll="1m",
  body = '{"sort": ["_doc"]}')
res$`_scroll_id`

# Pass scroll_id to scroll function
scroll(con, res$`_scroll_id`)

# Get all results - one approach is to use a while loop
res <- Search(con, index = 'shakespeare', q="a*", time_scroll="5m",
  body = '{"sort": ["_doc"]}')
out <- res$hits$hits
hits <- 1
while(hits != 0){
  res <- scroll(con, res$`_scroll_id`, time_scroll="5m")
  hits <- length(res$hits$hits)
  if(hits > 0)
    out <- c(out, res$hits$hits)
}
length(out)
res$hits$total
out[[1]]

# clear scroll
## individual scroll id
res <- Search(con, index = 'shakespeare', q="a*", time_scroll="5m",
  body = '{"sort": ["_doc"]}')
scroll_clear(con, res$`_scroll_id`)

## many scroll ids
res1 <- Search(con, index = 'shakespeare', q="c*", time_scroll="5m",
  body = '{"sort": ["_doc"]}')
res2 <- Search(con, index = 'shakespeare', q="d*", time_scroll="5m",
  body = '{"sort": ["_doc"]}')
nodes_stats(con, metric = "indices")$nodes[[1]]$indices$search$open_contexts
scroll_clear(con, c(res1$`_scroll_id`, res2$`_scroll_id`))
nodes_stats(con, metric = "indices")$nodes[[1]]$indices$search$open_contexts

## all scroll ids
res1 <- Search(con, index = 'shakespeare', q="f*", time_scroll="1m",
  body = '{"sort": ["_doc"]}')
res2 <- Search(con, index = 'shakespeare', q="g*", time_scroll="1m",
  body = '{"sort": ["_doc"]}')
res3 <- Search(con, index = 'shakespeare', q="k*", time_scroll="1m",
  body = '{"sort": ["_doc"]}')
scroll_clear(con, all = TRUE)

## sliced scrolling
body1 <- '{
  "slice": {
    "id": 0,
    "max": 2
  },
  "query": {
    "match" : {
      "text_entry" : "a*"
    }
  }
}'

body2 <- '{
  "slice": {
    "id": 1,
    "max": 2
  },
  "query": {
    "match" : {
      "text_entry" : "a*"
    }
  }
}'

res1 <- Search(con, index = 'shakespeare', time_scroll="1m", body = body1)
res2 <- Search(con, index = 'shakespeare', time_scroll="1m", body = body2)
scroll(con, res1$`_scroll_id`)
scroll(con, res2$`_scroll_id`)

out1 <- list()
hits <- 1
while(hits != 0){
  tmp1 <- scroll(con, res1$`_scroll_id`)
  hits <- length(tmp1$hits$hits)
  if(hits > 0)
    out1 <- c(out1, tmp1$hits$hits)
}

out2 <- list()
hits <- 1
while(hits != 0){
  tmp2 <- scroll(con, res2$`_scroll_id`)
  hits <- length(tmp2$hits$hits)
  if(hits > 0)
    out2 <- c(out2, tmp2$hits$hits)
}

c(
 lapply(out1, "[[", "_source"),
 lapply(out2, "[[", "_source")
)


# using jsonlite::stream_out
res <- Search(con, time_scroll = "1m")
file <- tempfile()
scroll(con, 
  x = res$`_scroll_id`,
  stream_opts = list(file = file)
)
jsonlite::stream_in(file(file))
unlink(file)

## stream_out and while loop
(file <- tempfile())
res <- Search(con, index = "shakespeare", time_scroll = "5m",
  size = 1000, stream_opts = list(file = file))
while(!inherits(res, "warning")) {
  res <- tryCatch(scroll(
    conn = con,
    x = res$`_scroll_id`,
    time_scroll = "5m",
    stream_opts = list(file = file)
  ), warning = function(w) w)
}
NROW(df <- jsonlite::stream_in(file(file)))
head(df)
}
}
\seealso{
\code{\link[=Search]{Search()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/count.r
\name{count}
\alias{count}
\title{Get counts of the number of records per index.}
\usage{
count(conn, index = NULL, type = NULL, callopts = list(), verbose = TRUE, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{Index, defaults to all indices}

\item{type}{Document type, optional}

\item{callopts}{Curl args passed on to \link[crul:verb-GET]{crul::verb-GET}}

\item{verbose}{If \code{TRUE} (default) the url call used printed to console.}

\item{...}{Further args passed on to elastic search HTTP API as parameters.}
}
\description{
Get counts of the number of records per index.
}
\details{
See docs for the count API here
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/search-count.html}

You can also get a count of documents using \code{\link[=Search]{Search()}} or
\code{\link[=Search_uri]{Search_uri()}} and setting \code{size = 0}
}
\examples{
\dontrun{
# connection setup
(x <- connect())

if (!index_exists(x, "plos")) {
  plosdat <- system.file("examples", "plos_data.json",
    package = "elastic")
  plosdat <- type_remover(plosdat)
  invisible(docs_bulk(x, plosdat))
}
if (!index_exists(x, "shakespeare")) {
  shake <- system.file("examples", "shakespeare_data_.json", 
    package = "elastic")
  invisible(docs_bulk(x, shake))
}

count(x)
count(x, index='plos')
count(x, index='shakespeare')
count(x, index=c('plos','shakespeare'), q="a*")
count(x, index=c('plos','shakespeare'), q="z*")

# Curl options
count(x, callopts = list(verbose = TRUE))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_update_by_query.R
\name{docs_update_by_query}
\alias{docs_update_by_query}
\title{Update documents by query}
\usage{
docs_update_by_query(
  conn,
  index,
  body = NULL,
  type = NULL,
  conflicts = NULL,
  routing = NULL,
  scroll_size = NULL,
  refresh = NULL,
  wait_for_completion = NULL,
  wait_for_active_shards = NULL,
  timeout = NULL,
  scroll = NULL,
  requests_per_second = NULL,
  pipeline = NULL,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) The name of the index. Required}

\item{body}{(character/json) query to be passed on to POST request body}

\item{type}{(character) The type of the document. optional}

\item{conflicts}{(character) If you’d like to count version conflicts
rather than cause them to abort then set \code{conflicts=proceed}}

\item{routing}{(character) Specific routing value}

\item{scroll_size}{(integer) By default uses scroll batches of 1000.
Change batch size with this parameter.}

\item{refresh}{(logical) Refresh the index after performing the operation}

\item{wait_for_completion}{(logical) If \code{wait_for_completion=FALSE} then
Elasticsearch will perform some preflight checks, launch the request, and
then return a task which can be used with Tasks APIs to cancel or get the
status of the task. Elasticsearch will also create a record of this task
as a document at .tasks/task/${taskId}. This is yours to keep or remove
as you see fit. When you are done with it, delete it so Elasticsearch
can reclaim the space it uses. Default: \code{TRUE}}

\item{wait_for_active_shards}{(logical) controls how many copies of a
shard must be active before proceeding with the request.}

\item{timeout}{(character) Explicit operation timeout, e.g,. 5m (for 5
minutes)}

\item{scroll}{(integer) control how long the "search context" is kept
alive, eg \code{scroll='10m'}, by default it’s 5 minutes (\verb{5m})}

\item{requests_per_second}{(integer) any positive decimal number
(1.4, 6, 1000, etc); throttles rate at which \verb{_delete_by_query} issues
batches of delete operations by padding each batch with a wait time.
The throttling can be disabled by setting \code{requests_per_second=-1}}

\item{pipeline}{(character) a pipeline name}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}
}
\description{
update documents by query via a POST request
}
\examples{
\dontrun{
(x <- connect())
x$ping()

omdb <- system.file("examples", "omdb.json", package = "elastic")
omdb <- type_remover(omdb)
if (!index_exists(x, "omdb")) invisible(docs_bulk(x, omdb))

# can be sent without a body
docs_update_by_query(x, index='omdb')

# update
## note this works with imdbRating, a float, but didn't seem to work
## with Metascore, a long
## See link above for Painless API reference
body <- '{
  "script": {
    "source": "ctx._source.imdbRating++",
    "lang": "painless"
  },
  "query": {
    "match": {
      "Rated": "R"
    }
  }
}'
Search(x, "omdb", q = "Rated:\"R\"", asdf=TRUE,
  source = c("Title", "Rated", "imdbRating"))$hits$hits
docs_update_by_query(x, index='omdb', body = body)
Search(x, "omdb", q = "Rated:\"R\"", asdf=TRUE,
  source = c("Title", "Rated", "imdbRating"))$hits$hits
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-update-by-query.html}
\url{https://www.elastic.co/guide/en/elasticsearch/painless/current/painless-api-reference.html}
}
\seealso{
\code{\link[=docs_delete_by_query]{docs_delete_by_query()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/parsers.r
\name{es_parse}
\alias{es_parse}
\alias{es_parse.es_GET}
\alias{es_parse.index_delete}
\alias{es_parse.bulk_make}
\alias{es_parse.elastic_mget}
\alias{es_parse.elastic_search}
\alias{es_parse.elastic_status}
\alias{es_parse.elastic_stats}
\alias{es_parse.elastic_cluster_health}
\alias{es_parse.elastic_cluster_state}
\alias{es_parse.elastic_cluster_settings}
\alias{es_parse.elastic_cluster_stats}
\alias{es_parse.elastic_cluster_pending_tasks}
\alias{es_parse.elastic_nodes_stats}
\alias{es_parse.elastic_nodes_info}
\title{Parse raw data from es_get, es_mget, or es_search.}
\usage{
es_parse(input, parsetype, verbose)

\method{es_parse}{es_GET}(input, parsetype = "list", verbose = FALSE)

\method{es_parse}{index_delete}(input, parsetype = "list", verbose = FALSE)

\method{es_parse}{bulk_make}(input, parsetype = "list", verbose = FALSE)

\method{es_parse}{elastic_mget}(input, parsetype = "list", verbose = FALSE)

\method{es_parse}{elastic_search}(input, parsetype = "list", verbose = FALSE)

\method{es_parse}{elastic_status}(input, parsetype = "list", verbose = FALSE)

\method{es_parse}{elastic_stats}(input, parsetype = "list", verbose = FALSE)

\method{es_parse}{elastic_cluster_health}(input, parsetype = "list", verbose = TRUE)

\method{es_parse}{elastic_cluster_health}(input, parsetype = "list", verbose = TRUE)

\method{es_parse}{elastic_cluster_state}(input, parsetype = "list", verbose = TRUE)

\method{es_parse}{elastic_cluster_settings}(input, parsetype = "list", verbose = TRUE)

\method{es_parse}{elastic_cluster_stats}(input, parsetype = "list", verbose = TRUE)

\method{es_parse}{elastic_cluster_pending_tasks}(input, parsetype = "list", verbose = TRUE)

\method{es_parse}{elastic_nodes_stats}(input, parsetype = "list", verbose = TRUE)

\method{es_parse}{elastic_nodes_info}(input, parsetype = "list", verbose = TRUE)
}
\arguments{
\item{input}{Output from solr_facet}

\item{parsetype}{One of 'list' or 'df' (data.frame). Only list possible for now.}

\item{verbose}{Print messages or not (default: FALSE).}
}
\description{
Parse raw data from es_get, es_mget, or es_search.
}
\details{
This is the parser used internally in es_get, es_mget, and es_search,
but if you output raw data from es_* functions using raw=TRUE, then you can use this
function to parse that data (a es_* S3 object) after the fact to a list of
data.frame's for easier consumption.
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_bulk_update.R
\name{docs_bulk_update}
\alias{docs_bulk_update}
\title{Use the bulk API to update documents}
\usage{
docs_bulk_update(
  conn,
  x,
  index = NULL,
  type = NULL,
  chunk_size = 1000,
  doc_ids = NULL,
  raw = FALSE,
  quiet = FALSE,
  query = list(),
  digits = NA,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{x}{A list, data.frame, or character path to a file. required.}

\item{index}{(character) The index name to use. Required for data.frame
input, but optional for file inputs.}

\item{type}{(character) The type. default: \code{NULL}. Note that \code{type} is
deprecated in Elasticsearch v7 and greater, and removed in Elasticsearch v8}

\item{chunk_size}{(integer) Size of each chunk. If your data.frame is smaller
thank \code{chunk_size}, this parameter is essentially ignored. We write in
chunks because at some point, depending on size of each document, and
Elasticsearch setup, writing a very large number of documents in one go
becomes slow, so chunking can help. This parameter is ignored if you
pass a file name. Default: 1000}

\item{doc_ids}{An optional vector (character or numeric/integer) of document
ids to use. This vector has to equal the size of the documents you are
passing in, and will error if not. If you pass a factor we convert to
character. Default: not passed}

\item{raw}{(logical) Get raw JSON back or not. If \code{TRUE}
you get JSON; if \code{FALSE} you get a list. Default: \code{FALSE}}

\item{quiet}{(logical) Suppress progress bar. Default: \code{FALSE}}

\item{query}{(list) a named list of query parameters. optional.
options include: pipeline, refresh, routing, _source, _source_excludes,
_source_includes, timeout, wait_for_active_shards. See the docs bulk
ES page for details}

\item{digits}{digits used by the parameter of the same name by
\code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}} to convert data to JSON before being submitted to
your ES instance. default: \code{NA}}

\item{...}{Pass on curl options to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Use the bulk API to update documents
}
\details{
\itemize{
\item \code{doc_as_upsert} - is set to \code{TRUE} for all records
}

For doing updates with a file already prepared for the bulk API,
see \code{\link[=docs_bulk]{docs_bulk()}}

Only data.frame's are supported for now.
}
\examples{
\dontrun{
x <- connect()
if (index_exists(x, "foobar")) index_delete(x, "foobar")

df <- data.frame(name = letters[1:3], size = 1:3, id = 100:102)
invisible(docs_bulk(x, df, 'foobar', es_ids = FALSE))

# add new rows in existing fields
(df2 <- data.frame(size = c(45, 56), id = 100:101))
(df2 <- data.frame(size = c(45, 56)))
df2$`_id` <- 100:101
df2
Search(x, "foobar", asdf = TRUE)$hits$hits
invisible(docs_bulk_update(x, df2, index = 'foobar'))
Search(x, "foobar", asdf = TRUE)$hits$hits

# add new fields (and new rows by extension)
(df3 <- data.frame(color = c("blue", "red", "green"), id = 100:102))
Search(x, "foobar", asdf = TRUE)$hits$hits
invisible(docs_bulk_update(x, df3, index = 'foobar'))
Sys.sleep(2) # wait for a few sec to make sure you see changes reflected
Search(x, "foobar", asdf = TRUE)$hits$hits
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-bulk.html#bulk-update}
}
\seealso{
Other bulk-functions: 
\code{\link{docs_bulk_create}()},
\code{\link{docs_bulk_delete}()},
\code{\link{docs_bulk_index}()},
\code{\link{docs_bulk_prep}()},
\code{\link{docs_bulk}()}
}
\concept{bulk-functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/mtermvectors.R
\name{mtermvectors}
\alias{mtermvectors}
\title{Multi Termvectors}
\usage{
mtermvectors(
  conn,
  index = NULL,
  type = NULL,
  ids = NULL,
  body = list(),
  pretty = TRUE,
  field_statistics = TRUE,
  fields = NULL,
  offsets = TRUE,
  parent = NULL,
  payloads = TRUE,
  positions = TRUE,
  preference = "random",
  realtime = TRUE,
  routing = NULL,
  term_statistics = FALSE,
  version = NULL,
  version_type = NULL,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) The index in which the document resides.}

\item{type}{(character) The type of the document.}

\item{ids}{(character) One or more document ids}

\item{body}{(character) Define parameters and or supply a document to get
termvectors for}

\item{pretty}{(logical) pretty print. Default: \code{TRUE}}

\item{field_statistics}{(character) Specifies if document count, sum of
document frequencies and sum of total term frequencies should be returned.
Default: \code{TRUE}}

\item{fields}{(character) A comma-separated list of fields to return.}

\item{offsets}{(character) Specifies if term offsets should be returned.
Default: \code{TRUE}}

\item{parent}{(character) Parent id of documents.}

\item{payloads}{(character) Specifies if term payloads should be returned.
Default: \code{TRUE}}

\item{positions}{(character) Specifies if term positions should be returned.
Default: \code{TRUE}}

\item{preference}{(character) Specify the node or shard the operation
should be performed on (Default: \code{random}).}

\item{realtime}{(character) Specifies if request is real-time as opposed to
near-real-time (Default: \code{TRUE}).}

\item{routing}{(character) Specific routing value.}

\item{term_statistics}{(character) Specifies if total term frequency and
document frequency should be returned. Default: \code{FALSE}}

\item{version}{(character) Explicit version number for concurrency control}

\item{version_type}{(character) Specific version type, valid choices are:
'internal', 'external', 'external_gte', 'force'}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}
}
\description{
Multi Termvectors
}
\details{
Multi termvectors API allows to get multiple termvectors based on an
index, type and id.
}
\examples{
\dontrun{
x <- connect()

if (index_exists(x, 'omdb')) index_delete(x, "omdb")
omdb <- system.file("examples", "omdb.json", package = "elastic")
omdb <- type_remover(omdb)
invisible(docs_bulk(x, omdb))
out <- Search(x, "omdb", size = 2)$hits$hits
ids <- vapply(out, "[[", "", "_id")

# no index
body <- '{
   "docs": [
      {
         "_index": "omdb",
         "_id": "\%s",
         "term_statistics": true
      },
      {
         "_index": "omdb",
         "_id": "\%s",
         "fields": [
            "Plot"
         ]
      }
   ]
}'
mtermvectors(x, body = sprintf(body, ids[1], ids[2]))

# index given
body <- '{
   "docs": [
      {
         "_id": "\%s",
         "fields": [
            "Plot"
         ],
         "term_statistics": true
      },
      {
         "_id": "\%s",
         "fields": [
            "Title"
         ]
      }
   ]
}'
mtermvectors(x, 'omdb', body = sprintf(body, ids[1], ids[2]))

# parameters same for both documents, so can simplify
body <- '{
    "ids" : ["\%s", "\%s"],
    "parameters": {
        "fields": [
            "Plot"
        ],
        "term_statistics": true
    }
}'
mtermvectors(x, 'omdb', body = sprintf(body, ids[1], ids[2]))

# you can give user provided documents via the 'docs' parameter
## though you have to give index and type that exist in your Elasticsearch 
## instance
body <- '{
   "docs": [
      {
         "_index": "omdb",
         "doc" : {
            "Director" : "John Doe",
            "Plot" : "twitter test test test"
         }
      },
      {
         "_index": "omdb",
         "doc" : {
           "Director" : "Jane Doe",
           "Plot" : "Another twitter test ..."
         }
      }
   ]
}'
mtermvectors(x, body = body)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-multi-termvectors.html}
}
\seealso{
\code{\link[=termvectors]{termvectors()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/alias.R
\name{alias}
\alias{alias}
\alias{alias_get}
\alias{aliases_get}
\alias{alias_exists}
\alias{alias_create}
\alias{alias_rename}
\alias{alias_delete}
\title{Elasticsearch alias APIs}
\usage{
alias_get(conn, index = NULL, alias = NULL, ignore_unavailable = FALSE, ...)

aliases_get(conn, index = NULL, alias = NULL, ignore_unavailable = FALSE, ...)

alias_exists(conn, index = NULL, alias = NULL, ...)

alias_create(
  conn,
  index,
  alias,
  filter = NULL,
  routing = NULL,
  search_routing = NULL,
  index_routing = NULL,
  ...
)

alias_rename(conn, index, alias, alias_new, ...)

alias_delete(conn, index = NULL, alias, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) An index name}

\item{alias}{(character) An alias name}

\item{ignore_unavailable}{(logical) What to do if an specified index name
doesn't exist. If set to \code{TRUE} then those indices are ignored.}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}, \link[crul:verb-GET]{crul::verb-GET},
\link[crul:verb-HEAD]{crul::verb-HEAD}, or \link[crul:verb-DELETE]{crul::verb-DELETE}}

\item{filter}{(named list) provides an easy way to create different "views" of
the same index. Defined using Query DSL and is applied to all Search, Count,
Delete By Query and More Like This operations with this alias. See
examples}

\item{routing, search_routing, index_routing}{(character) Associate a routing
value with an alias}

\item{alias_new}{(character) A new alias name, used in rename only}
}
\description{
Elasticsearch alias APIs
}
\details{
Note that you can also create aliases when you create indices
by putting the directive in the request body. See the Elasticsearch
docs link
}
\examples{
\dontrun{
# connection setup
(x <- connect())

if (!index_exists(x, "plos")) {
  plosdat <- system.file("examples", "plos_data.json", package = "elastic")
  invisible(docs_bulk(x, plosdat))
}
if (!index_exists(x, "shakespeare")) {
  shake <- system.file("examples", "shakespeare_data_.json", package = "elastic")
  invisible(docs_bulk(x, shake))
}

# Create/update an alias
alias_create(x, index = "plos", alias = "candles")
## more than one alias
alias_create(x, index = "plos", alias = c("tables", "chairs"))

# associate an alias with two multiple different indices
alias_create(x, index = c("plos", "shakespeare"), alias = "stools")

# Retrieve a specified alias
alias_get(x, index="plos")
alias_get(x, alias="tables")
alias_get(x, alias="stools")
aliases_get(x)

# rename an alias
aliases_get(x, "plos")
alias_rename(x, index = 'plos', alias = "stools", alias_new = "plates")
aliases_get(x, "plos")

# filtered aliases
alias_create(x, index = "plos", alias = "candles", 
  filter = list(wildcard = list(title = "cell")))
## a search with the alias should give titles with cell in them
(titles <- Search(x, "candles", asdf = TRUE)$hits$hits$`_source.title`)
grepl("cell", titles, ignore.case = TRUE)

# routing
alias_create(x, index = "plos", alias = "candles", 
  routing = "1")

# Check for alias existence
alias_exists(x, index = "plos")
alias_exists(x, alias = "tables")
alias_exists(x, alias = "adsfasdf")

# Delete an alias
alias_delete(x, index = "plos", alias = "tables")
alias_exists(x, alias = "tables")

# Curl options
alias_create(x, index = "plos", alias = "tables")
aliases_get(x, alias = "tables", verbose = TRUE)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/indices-aliases.html}
}
\author{
Scott Chamberlain \href{mailto:myrmecocystus@gmail.com}{myrmecocystus@gmail.com}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_bulk_delete.R
\name{docs_bulk_delete}
\alias{docs_bulk_delete}
\title{Use the bulk API to delete documents}
\usage{
docs_bulk_delete(
  conn,
  x,
  index = NULL,
  type = NULL,
  chunk_size = 1000,
  doc_ids = NULL,
  raw = FALSE,
  quiet = FALSE,
  query = list(),
  digits = NA,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{x}{A list, data.frame, or character path to a file. required.}

\item{index}{(character) The index name to use. Required for data.frame
input, but optional for file inputs.}

\item{type}{(character) The type. default: \code{NULL}. Note that \code{type} is
deprecated in Elasticsearch v7 and greater, and removed in Elasticsearch v8}

\item{chunk_size}{(integer) Size of each chunk. If your data.frame is smaller
thank \code{chunk_size}, this parameter is essentially ignored. We write in
chunks because at some point, depending on size of each document, and
Elasticsearch setup, writing a very large number of documents in one go
becomes slow, so chunking can help. This parameter is ignored if you
pass a file name. Default: 1000}

\item{doc_ids}{An optional vector (character or numeric/integer) of document
ids to use. This vector has to equal the size of the documents you are
passing in, and will error if not. If you pass a factor we convert to
character. Default: not passed}

\item{raw}{(logical) Get raw JSON back or not. If \code{TRUE}
you get JSON; if \code{FALSE} you get a list. Default: \code{FALSE}}

\item{quiet}{(logical) Suppress progress bar. Default: \code{FALSE}}

\item{query}{(list) a named list of query parameters. optional.
options include: pipeline, refresh, routing, _source, _source_excludes,
_source_includes, timeout, wait_for_active_shards. See the docs bulk
ES page for details}

\item{digits}{ignored, used in other docs bulk functions but not used here}

\item{...}{Pass on curl options to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Use the bulk API to delete documents
}
\details{
For doing deletes with a file already prepared for the bulk API,
see \code{\link[=docs_bulk]{docs_bulk()}}

Only data.frame's are supported for now.
}
\examples{
\dontrun{
x <- connect()
if (index_exists(x, "foobar")) index_delete(x, "foobar")

df <- data.frame(name = letters[1:3], size = 1:3, id = 100:102)
invisible(docs_bulk(x, df, 'foobar', es_ids = FALSE))
Search(x, "foobar", asdf = TRUE)$hits$hits

# delete using doc ids from the data.frame you used to create
invisible(docs_bulk_delete(x, df, index = 'foobar'))
Search(x, "foobar", asdf = TRUE)$hits$total$value

# delete by passing in doc ids
## recreate data first
if (index_exists(x, "foobar")) index_delete(x, "foobar")
df <- data.frame(name = letters[1:3], size = 1:3, id = 100:102)
invisible(docs_bulk(x, df, 'foobar', es_ids = FALSE))
docs_bulk_delete(x, df, index = 'foobar', doc_ids = df$id)
Search(x, "foobar", asdf = TRUE)$hits$total$value
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-bulk.html}
}
\seealso{
Other bulk-functions: 
\code{\link{docs_bulk_create}()},
\code{\link{docs_bulk_index}()},
\code{\link{docs_bulk_prep}()},
\code{\link{docs_bulk_update}()},
\code{\link{docs_bulk}()}
}
\concept{bulk-functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpdocs.R
\name{units-distance}
\alias{units-distance}
\title{Distance units}
\description{
Wherever distances need to be specified, such as the distance parameter in the
Geo Distance Filter), the default unit if none is specified is the meter. Distances
can be specified in other units, such as "1km" or "2mi" (2 miles).
}
\details{
\tabular{ll}{
mi or miles \tab Mile \cr
yd or yards \tab Yard \cr
ft or feet \tab Feet \cr
in or inch \tab Inch \cr
km or kilometers \tab Kilometer \cr
m or meters \tab Meter \cr
cm or centimeters \tab Centimeter \cr
mm or millimeters \tab Millimeter \cr
NM, nmi or nauticalmiles \tab Nautical mile \cr
}

The precision parameter in the Geohash Cell Filter accepts distances with the above
units, but if no unit is specified, then the precision is interpreted as the length
of the geohash.
}
\seealso{
\link{units-time}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_delete_by_query.R
\name{docs_delete_by_query}
\alias{docs_delete_by_query}
\title{Delete documents by query}
\usage{
docs_delete_by_query(
  conn,
  index,
  body,
  type = NULL,
  conflicts = NULL,
  routing = NULL,
  scroll_size = NULL,
  refresh = NULL,
  wait_for_completion = NULL,
  wait_for_active_shards = NULL,
  timeout = NULL,
  scroll = NULL,
  requests_per_second = NULL,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) The name of the index. Required}

\item{body}{(character/json) query to be passed on to POST request body}

\item{type}{(character) The type of the document. optional}

\item{conflicts}{(character) If you’d like to count version conflicts
rather than cause them to abort then set \code{conflicts=proceed}}

\item{routing}{(character) Specific routing value}

\item{scroll_size}{(integer) By default uses scroll batches of 1000.
Change batch size with this parameter.}

\item{refresh}{(logical) Refresh the index after performing the operation}

\item{wait_for_completion}{(logical) If \code{wait_for_completion=FALSE} then
Elasticsearch will perform some preflight checks, launch the request, and
then return a task which can be used with Tasks APIs to cancel or get the
status of the task. Elasticsearch will also create a record of this task
as a document at .tasks/task/${taskId}. This is yours to keep or remove
as you see fit. When you are done with it, delete it so Elasticsearch
can reclaim the space it uses. Default: \code{TRUE}}

\item{wait_for_active_shards}{(logical) controls how many copies of a
shard must be active before proceeding with the request.}

\item{timeout}{(character) Explicit operation timeout, e.g,. 5m (for 5
minutes)}

\item{scroll}{(integer) control how long the "search context" is kept
alive, eg \code{scroll='10m'}, by default it’s 5 minutes (\verb{5m})}

\item{requests_per_second}{(integer) any positive decimal number
(1.4, 6, 1000, etc); throttles rate at which \verb{_delete_by_query} issues
batches of delete operations by padding each batch with a wait time.
The throttling can be disabled by setting \code{requests_per_second=-1}}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}
}
\description{
delete documents by query via a POST request
}
\examples{
\dontrun{
(x <- connect())
x$ping()

plosdat <- system.file("examples", "plos_data.json",
  package = "elastic")
plosdat <- type_remover(plosdat)
if (!index_exists(x, "plos")) invisible(docs_bulk(x, plosdat))

# delete with fuzzy matching
body <- '{
  "query": { 
    "match": {
      "title": {
        "query": "cells",
        "fuzziness": 1
      }
    }
  }
}'
docs_delete_by_query(x, index='plos', body = body) 

# delete with no fuzziness
if (index_exists(x, "plos")) index_delete(x, 'plos')
invisible(docs_bulk(x, plosdat))
count(x, "plos")
body <- '{
  "query": { 
    "match": {
      "title": {
        "query": "cells",
        "fuzziness": 0
      }
    }
  }
}'
docs_delete_by_query(x, index='plos', body = body)

# delete all docs with match_all query
if (index_exists(x, "plos")) index_delete(x, 'plos')
invisible(docs_bulk(x, plosdat))
body <- '{
  "query": { 
    "match_all": {}
  }
}'
docs_delete_by_query(x, index='plos', body = body)

# put plos back in 
if (index_exists(x, "plos")) index_delete(x, 'plos')
invisible(docs_bulk(x, plosdat))

# delete docs from more than one index
foo <- system.file("examples/foo.json", package = "elastic")
if (!index_exists(x, "foo")) invisible(docs_bulk(x, foo))
bar <- system.file("examples/bar.json", package = "elastic")
if (!index_exists(x, "bar")) invisible(docs_bulk(x, bar))

body <- '{
  "query": { 
    "match_all": {}
  }
}'
docs_delete_by_query(x, index=c('foo','bar'), 
  body = body, verbose = TRUE)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-delete-by-query.html}
}
\seealso{
\code{\link[=docs_update_by_query]{docs_update_by_query()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_shards.R
\name{search_shards}
\alias{search_shards}
\title{Search shards}
\usage{
search_shards(
  conn,
  index = NULL,
  raw = FALSE,
  routing = NULL,
  preference = NULL,
  local = NULL,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{One or more indeces}

\item{raw}{If \code{TRUE} (default), data is parsed to list. If \code{FALSE}, then
raw JSON}

\item{routing}{A character vector of routing values to take into account
when determining which shards a request would be executed against.}

\item{preference}{Controls a preference of which shard replicas to execute
the search request on. By default, the operation is randomized between the
shard replicas. See \link{preference} for a list of all acceptable
values.}

\item{local}{(logical) Whether to read the cluster state locally in order
to determine where shards are allocated instead of using the Master node's
cluster state.}

\item{...}{Curl args passed on to \link[crul:verb-GET]{crul::verb-GET}}
}
\description{
Search shards
}
\examples{
\dontrun{
# connection setup
(x <- connect())

search_shards(x, index = "plos")
search_shards(x, index = c("plos","gbif"))
search_shards(x, index = "plos", preference='_primary')
search_shards(x, index = "plos", preference='_shards:2')

# curl options
search_shards(x, index = "plos", verbose = TRUE)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/search-shards.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Search_uri.R
\name{Search_uri}
\alias{Search_uri}
\title{Full text search of Elasticsearch with URI search}
\usage{
Search_uri(
  conn,
  index = NULL,
  type = NULL,
  q = NULL,
  df = NULL,
  analyzer = NULL,
  default_operator = NULL,
  explain = NULL,
  source = NULL,
  fields = NULL,
  sort = NULL,
  track_scores = NULL,
  timeout = NULL,
  terminate_after = NULL,
  from = NULL,
  size = NULL,
  search_type = NULL,
  lowercase_expanded_terms = NULL,
  analyze_wildcard = NULL,
  version = NULL,
  lenient = NULL,
  raw = FALSE,
  asdf = FALSE,
  track_total_hits = TRUE,
  search_path = "_search",
  stream_opts = list(),
  ignore_unavailable = FALSE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link{connect}}}

\item{index}{Index name, one or more}

\item{type}{Document type. Note that \code{type} is deprecated in
Elasticsearch v7 and greater, and removed in Elasticsearch v8. We will
strive to support types for folks using older ES versions}

\item{q}{The query string (maps to the query_string query, see Query String
Query for more details). See
https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-query-string-query.html
for documentation and examples.}

\item{df}{(character) The default field to use when no field prefix is
defined within the query.}

\item{analyzer}{(character) The analyzer name to be used when analyzing the
query string.}

\item{default_operator}{(character) The default operator to be used, can be
\code{AND} or \code{OR}. Default: \code{OR}}

\item{explain}{(logical) For each hit, contain an explanation of how
scoring of the hits was computed. Default: \code{FALSE}}

\item{source}{(logical) Set to \code{FALSE} to disable retrieval of the
\code{_source} field. You can also retrieve part of the document by
using \code{_source_include} & \code{_source_exclude} (see the \code{body}
documentation for more details). You can also include a comma-delimited
string of fields from the source document that you want back. See also
the \strong{fields} parameter}

\item{fields}{(character) The selective stored fields of the document to
return for each hit. Not specifying any value will cause no fields to return.
Note that in Elasticsearch v5 and greater, \strong{fields} parameter has
changed to \strong{stored_fields}, which is not on by default. You can
however, pass fields to \strong{source} parameter}

\item{sort}{(character) Sorting to perform. Can either be in the form of
fieldName, or \code{fieldName:asc}/\code{fieldName:desc}. The fieldName
can either be an actual field within the document, or the special
\code{_score} name to indicate sorting based on scores. There can be several
sort parameters (order is important).}

\item{track_scores}{(logical) When sorting, set to \code{TRUE} in order to
still track scores and return them as part of each hit.}

\item{timeout}{(numeric) A search timeout, bounding the search request to
be executed within the specified time value and bail with the hits
accumulated up to that point when expired. Default: no timeout.}

\item{terminate_after}{(numeric) The maximum number of documents to collect
for each shard, upon reaching which the query execution will terminate
early. If set, the response will have a boolean field terminated_early to
indicate whether the query execution has actually terminated_early.
Default: no terminate_after}

\item{from}{(character) The starting from index of the hits to return.
Pass in as a character string to avoid problems with large number
conversion to scientific notation. Default: 0}

\item{size}{(character) The number of hits to return. Pass in as a
character string to avoid problems with large number conversion to
scientific notation. Default: 10. The default maximum is 10,000 - however,
you can change this default maximum by changing the
\code{index.max_result_window} index level parameter.}

\item{search_type}{(character) The type of the search operation to perform.
Can be \code{query_then_fetch} (default) or \code{dfs_query_then_fetch}.
Types \code{scan} and \code{count} are deprecated.
See Elasticsearch docs for more details on the different types of
search that can be performed.}

\item{lowercase_expanded_terms}{(logical) Should terms be automatically
lowercased or not. Default: \code{TRUE}.}

\item{analyze_wildcard}{(logical) Should wildcard and prefix queries be
analyzed or not. Default: \code{FALSE}.}

\item{version}{(logical) Print the document version with each document.}

\item{lenient}{(logical) If \code{TRUE} will cause format based failures (like
providing text to a numeric field) to be ignored. Default: \code{NULL}}

\item{raw}{(logical) If \code{FALSE} (default), data is parsed to list.
If \code{TRUE}, then raw JSON returned}

\item{asdf}{(logical) If \code{TRUE}, use \code{\link[jsonlite]{fromJSON}}
to parse JSON directly to a data.frame. If \code{FALSE} (Default), list
output is given.}

\item{track_total_hits}{(logical, numeric) If \code{TRUE} will always track
the number of hits that match the query accurately. If \code{FALSE} will
count documents accurately up to 10000 documents. If \code{is.integer} will
count documents accurately up to the number. Default: \code{TRUE}}

\item{search_path}{(character) The path to use for searching. Default
to \verb{_search}, but in some cases you may already have that in the base
url set using \code{\link[=connect]{connect()}}, in which case you can set this
to \code{NULL}}

\item{stream_opts}{(list) A list of options passed to
\code{\link[jsonlite]{stream_out}} - Except that you can't pass \code{x} as
that's the data that's streamed out, and pass a file path instead of a
connection to \code{con}. \code{pagesize} param doesn't do much as
that's more or less controlled by paging with ES.}

\item{ignore_unavailable}{(logical) What to do if an specified index name
doesn't exist. If set to \code{TRUE} then those indices are ignored.}

\item{...}{Curl args passed on to \code{\link[crul]{verb-POST}}}
}
\description{
Full text search of Elasticsearch with URI search
}
\examples{
\dontrun{
# connection setup
(x <- connect())

# URI string queries
Search_uri(x, index="shakespeare")
## if you're using an older ES version, you may have types
if (gsub("\\\\.", "", x$ping()$version$number) < 700) {
Search_uri(x, index="shakespeare", type="act")
Search_uri(x, index="shakespeare", type="scene")
Search_uri(x, index="shakespeare", type="line")
}

## Return certain fields
if (gsub("\\\\.", "", ping()$version$number) < 500) {
  ### ES < v5
  Search_uri(x, index="shakespeare", fields=c('play_name','speaker'))
} else {
  ### ES > v5
  Search_uri(x, index="shakespeare", source=c('play_name','speaker'))
}

## Search many indices
Search_uri(x, index = "gbif")$hits$total$value
Search_uri(x, index = "shakespeare")$hits$total$value
Search_uri(x, index = c("gbif", "shakespeare"))$hits$total$value

## search_type
## NOTE: If you're in ES V5 or greater, see \code{?fielddata}
Search_uri(x, index="shakespeare", search_type = "query_then_fetch")
Search_uri(x, index="shakespeare", search_type = "dfs_query_then_fetch")
# Search_uri(x, index="shakespeare", search_type = "scan") # only when scrolling

## sorting
Search_uri(x, index="shakespeare", sort="text_entry")
if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
  Search_uri(x, index="shakespeare", sort="speaker:desc", fields='speaker')
  Search_uri(x, index="shakespeare", sort=c("speaker:desc","play_name:asc"),
    fields=c('speaker','play_name'))
}

## pagination
Search_uri(x, index="shakespeare", size=1)$hits$hits
Search_uri(x, index="shakespeare", size=1, from=1)$hits$hits

## queries
### Search in all fields
Search_uri(x, index="shakespeare", q="york")

### Searchin specific fields
Search_uri(x, index="shakespeare", q="speaker:KING HENRY IV")$hits$total$value

### Exact phrase search by wrapping in quotes
Search_uri(x, index="shakespeare", q='speaker:"KING HENRY IV"')$hits$total$value

### can specify operators between multiple words parenthetically
Search_uri(x, index="shakespeare", q="speaker:(HENRY OR ARCHBISHOP)")$hits$total$value

### where the field line_number has no value (or is missing)
Search_uri(x, index="shakespeare", q="_missing_:line_number")$hits$total$value

### where the field line_number has any non-null value
Search_uri(x, index="shakespeare", q="_exists_:line_number")$hits$total$value

### wildcards, either * or ?
Search_uri(x, index="shakespeare", q="*ay")$hits$total$value
Search_uri(x, index="shakespeare", q="m?y")$hits$total$value

### regular expressions, wrapped in forward slashes
Search_uri(x, index="shakespeare", q="text_entry:/[a-z]/")$hits$total$value

### fuzziness
Search_uri(x, index="shakespeare", q="text_entry:ma~")$hits$total$value
Search_uri(x, index="shakespeare", q="text_entry:the~2")$hits$total$value
Search_uri(x, index="shakespeare", q="text_entry:the~1")$hits$total$value

### Proximity searches
Search_uri(x, index="shakespeare", q='text_entry:"as hath"~5')$hits$total$value
Search_uri(x, index="shakespeare", q='text_entry:"as hath"~10')$hits$total$value

### Ranges, here where line_id value is between 10 and 20
Search_uri(x, index="shakespeare", q="line_id:[10 TO 20]")$hits$total$value

### Grouping
Search_uri(x, index="shakespeare", q="(hath OR as) AND the")$hits$total$value

# Limit number of hits returned with the size parameter
Search_uri(x, index="shakespeare", size=1)

# Give explanation of search in result
Search_uri(x, index="shakespeare", size=1, explain=TRUE)

## terminate query after x documents found
## setting to 1 gives back one document for each shard
Search_uri(x, index="shakespeare", terminate_after=1)
## or set to other number
Search_uri(x, index="shakespeare", terminate_after=2)

## Get version number for each document
Search_uri(x, index="shakespeare", version=TRUE, size=2)

## Get raw data
Search_uri(x, index="shakespeare", raw=TRUE)

## Curl options
### verbose
out <- Search_uri(x, index="shakespeare", verbose = TRUE)
}
}
\seealso{
\code{\link[=fielddata]{fielddata()}}

\code{\link[=Search]{Search()}} \code{\link[=Search_template]{Search_template()}} \code{\link[=count]{count()}} \code{\link[=fielddata]{fielddata()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_create.R
\name{docs_create}
\alias{docs_create}
\title{Create a document}
\usage{
docs_create(
  conn,
  index,
  body,
  type = NULL,
  id = NULL,
  version = NULL,
  version_type = NULL,
  op_type = NULL,
  routing = NULL,
  parent = NULL,
  timestamp = NULL,
  ttl = NULL,
  refresh = NULL,
  timeout = NULL,
  callopts = list(),
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) The name of the index. Required}

\item{body}{The document}

\item{type}{(character) The type of the document. optional}

\item{id}{(numeric/character) The document ID. Can be numeric or character.
Optional. if not provided, Elasticsearch creates the ID for you as a UUID.}

\item{version}{(character) Explicit version number for concurrency control}

\item{version_type}{(character) Specific version type. One of internal,
external, external_gte, or force}

\item{op_type}{(character) Operation type. One of create, or ...}

\item{routing}{(character) Specific routing value}

\item{parent}{(numeric) A parent document ID}

\item{timestamp}{(date) Explicit timestamp for the document}

\item{ttl}{(aka \dQuote{time to live}) Expiration time for the document.
Expired documents will be expunged automatically. The expiration date that
will be set for a document with a provided ttl is relative to the timestamp
of the document,  meaning it can be based on the time of indexing or on
any time provided. The provided ttl must be strictly positive and can be
a number (in milliseconds) or any valid time value (e.g, 86400000, 1d).}

\item{refresh}{(logical) Refresh the index after performing the operation}

\item{timeout}{(character) Explicit operation timeout, e.g,. 5m (for
5 minutes)}

\item{callopts}{Curl options passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{...}{Further args to query DSL}
}
\description{
Create a document
}
\examples{
\dontrun{
(x <- connect())

if (!index_exists(x, 'plos')) {
  plosdat <- system.file("examples", "plos_data.json",
    package = "elastic")
  plosdat <- type_remover(plosdat)
  invisible(docs_bulk(x, plosdat))
}

# give a document id
z <- docs_create(x, index = 'plos', id = 1002,
  body = list(id = "12345", title = "New title"))
z
# and the document is there now
docs_get(x, index = 'plos', id = 1002)

# let Elasticsearch create the document id for you
z <- docs_create(x, index='plos', body=list(id="6789", title="Some title"))
z
# and the document is there now
docs_get(x, index='plos', id=z$`_id`)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-index_.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cat.r
\name{cat}
\alias{cat}
\alias{cat_}
\alias{cat_aliases}
\alias{cat_allocation}
\alias{cat_count}
\alias{cat_segments}
\alias{cat_health}
\alias{cat_indices}
\alias{cat_master}
\alias{cat_nodes}
\alias{cat_nodeattrs}
\alias{cat_pending_tasks}
\alias{cat_plugins}
\alias{cat_recovery}
\alias{cat_thread_pool}
\alias{cat_shards}
\alias{cat_fielddata}
\title{Use the cat Elasticsearch api.}
\usage{
cat_(conn, parse = FALSE, ...)

cat_aliases(
  conn,
  verbose = FALSE,
  index = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  expand_wildcards = "all",
  ...
)

cat_allocation(
  conn,
  verbose = FALSE,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_count(
  conn,
  verbose = FALSE,
  index = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_segments(
  conn,
  verbose = FALSE,
  index = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_health(
  conn,
  verbose = FALSE,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_indices(
  conn,
  verbose = FALSE,
  index = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_master(
  conn,
  verbose = FALSE,
  index = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_nodes(
  conn,
  verbose = FALSE,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_nodeattrs(
  conn,
  verbose = FALSE,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_pending_tasks(
  conn,
  verbose = FALSE,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_plugins(
  conn,
  verbose = FALSE,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_recovery(
  conn,
  verbose = FALSE,
  index = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_thread_pool(
  conn,
  verbose = FALSE,
  index = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_shards(
  conn,
  verbose = FALSE,
  index = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)

cat_fielddata(
  conn,
  verbose = FALSE,
  index = NULL,
  fields = NULL,
  h = NULL,
  help = FALSE,
  bytes = FALSE,
  parse = FALSE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{parse}{(logical) Parse to a data.frame or not. Default: \code{FALSE}}

\item{...}{Curl args passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{verbose}{(logical) If \code{TRUE} (default) the url call used printed to
console}

\item{index}{(character) Index name}

\item{h}{(character) Fields to return}

\item{help}{(logical) Output available columns, and their meanings}

\item{bytes}{(logical) Give numbers back machine friendly. Default: \code{FALSE}}

\item{expand_wildcards}{(character) Whether to expand wildcard expression
to concrete indices that are open, closed or both.  Valid choices: 'open',
'closed', 'hidden', 'none', 'all'. default: 'all'. Available in ES >= v7.7}

\item{fields}{(character) Fields to return, only used with \code{fielddata}}
}
\description{
Use the cat Elasticsearch api.
}
\details{
See \url{https://www.elastic.co/guide/en/elasticsearch/reference/current/cat.html}
for the cat API documentation.

Note how \code{\link[=cat_]{cat_()}} has an underscore at the end to avoid conflict with the
function \code{\link[base:cat]{base::cat()}} in base R.
}
\examples{
\dontrun{
# connection setup
(x <- connect())

# list Elasticsearch cat endpoints
cat_(x)

# Do other cat operations
cat_aliases(x)
alias_create(x, index = "plos", alias = c("tables", "chairs"))
cat_aliases(x, expand_wildcards='open')
cat_aliases(x, expand_wildcards='all')
cat_allocation(x)
cat_allocation(x, verbose=TRUE)
cat_count(x)
cat_count(x, index='plos')
cat_count(x, index='gbif')
cat_segments(x)
cat_segments(x, index='gbif')
cat_health(x)
cat_indices(x)
cat_master(x)
cat_nodes(x)
# cat_nodeattrs(x) # not available in older ES versions
cat_pending_tasks(x)
cat_plugins(x)
cat_recovery(x, verbose=TRUE)
cat_recovery(x, index='gbif')
cat_thread_pool(x)
cat_thread_pool(x, verbose=TRUE)
cat_shards(x)
cat_fielddata(x)
cat_fielddata(x, fields='body')

# capture cat data into a data.frame
cat_(x, parse = TRUE)
cat_indices(x, parse = TRUE)
cat_indices(x, parse = TRUE, verbose = TRUE)
cat_count(x, parse = TRUE)
cat_count(x, parse = TRUE, verbose = TRUE)
cat_health(x, parse = TRUE)
cat_health(x, parse = TRUE, verbose = TRUE)

# Get help - what does each column mean
head(cat_indices(x, help = TRUE, parse = TRUE))
cat_health(x, help = TRUE, parse = TRUE)
head(cat_nodes(x, help = TRUE, parse = TRUE))

# Get back only certain fields
cat_nodes(x)
cat_nodes(x, h = c('ip','port','heapPercent','name'))
cat_nodes(x, h = c('id', 'ip', 'port', 'v', 'm'))
cat_indices(x, verbose = TRUE)
cat_indices(x, verbose = TRUE, h = c('index','docs.count','store.size'))

# Get back machine friendly numbers instead of the normal human friendly
cat_indices(x, verbose = TRUE, bytes = TRUE)

# Curl options
# cat_count(x, timeout_ms = 1)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Elasticsearch.R
\name{connect}
\alias{connect}
\title{Set connection details to an Elasticsearch engine.}
\usage{
connect(
  host = "127.0.0.1",
  port = 9200,
  path = NULL,
  transport_schema = "http",
  user = NULL,
  pwd = NULL,
  headers = NULL,
  cainfo = NULL,
  force = FALSE,
  errors = "simple",
  warn = TRUE,
  ignore_version = FALSE,
  ...
)
}
\arguments{
\item{host}{(character) The base host, defaults to \verb{127.0.0.1}}

\item{port}{(character) port to connect to, defaults to \code{9200}
(optional)}

\item{path}{(character) context path that is appended to the end of the
url. Default: \code{NULL}, ignored}

\item{transport_schema}{(character) http or https. Default: \code{http}}

\item{user}{(character) User name, if required for the connection. You
can specify,  but ignored for now.}

\item{pwd}{(character) Password, if required for the connection. You
can specify, but ignored for now.}

\item{headers}{named list of headers. These headers are used in all requests}

\item{cainfo}{(character) path to a crt bundle, passed to curl option
\code{cainfo}}

\item{force}{(logical) Force re-load of connection details.
Default: \code{FALSE}}

\item{errors}{(character) One of simple (Default) or complete. Simple gives
http code and  error message on an error, while complete gives both http
code and error message,  and stack trace, if available.}

\item{warn}{(logical) whether to throw warnings from the Elasticsearch
server when provided. Pulls warnings from response headers when given.
default: \code{TRUE}. To turn these off, you can set \code{warn=FALSE} or
wrap function calls in \code{\link[=suppressWarnings]{suppressWarnings()}}. You can also see warnings in
headers by using curl verbose.}

\item{ignore_version}{(logical) ignore Elasticsearch version checks?
default: \code{FALSE}. Setting this to \code{TRUE} may cause some problems, it
has not been fully tested yet. You may want to set this to \code{TRUE} if
it's not possible to ping the root route of the Elasticsearch instance,
which has the Elasticsearch version. We use the version to do
alter what request is sent as different Elasticsearch versions allow
different parameters.}

\item{...}{additional curl options to be passed in ALL http requests}
}
\description{
Set connection details to an Elasticsearch engine.
}
\details{
The default configuration is set up for localhost access on port
9200, with no username or password.

Running this connection method doesn't ping the ES server, but only prints
your connection details.

All connection details are stored within the returned object. We used to
store them in various env vars, but are now contained within the object
so you can have any number of connection objects and they shouldn't
conflict with one another.
}
\section{What is the connection object?}{

Creating a connection object with \code{connect()} does not create
a DBI-like connection object. DBI-like objects have externalptr, etc.,
while \code{connect()} simply holds details about your Elasticsearch
instance (host, port, authentication, etc.) that is used by other
methods in this package to interact with your instances' ES API.
\code{connect()} is more or less a fancy list.

You can connect to different Elasticsearch intances within the same
R session by creating a separate connection object for each instance;
then pass the appropriate connection object to each \code{elastic} method.
}

\examples{
\dontrun{
# the default is set to 127.0.0.1 (i.e., localhost) and port 9200
(x <- connect())
x$make_url()
x$ping()

# pass connection object to function calls
Search(x, q = "*:*")

# set username/password (hidden in print method)
connect(user = "me", pwd = "stuff")

# set a different host
# connect(host = '162.243.152.53')
# => http://162.243.152.53:9200

# set a different port
# connect(port = 8000)
# => http://localhost:8000

# set a different context path
# connect(path = 'foo_bar')
# => http://localhost:9200/foo_bar

# set to https
# connect(transport_schema = 'https')
# => https://localhost:9200

# set headers
connect(headers = list(a = 'foobar'))

# set cainfo path (hidden in print method)
connect(cainfo = '/some/path/bundle.crt')
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{nodes_shutdown}
\alias{nodes_shutdown}
\title{This function is defunct}
\usage{
nodes_shutdown(...)
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/percolater.R
\name{percolate}
\alias{percolate}
\alias{percolate_register}
\alias{percolate_match}
\alias{percolate_list}
\alias{percolate_count}
\alias{percolate_delete}
\title{Percolater}
\usage{
percolate_register(
  conn,
  index,
  id,
  type = NULL,
  body = list(),
  routing = NULL,
  preference = NULL,
  ignore_unavailable = NULL,
  percolate_format = NULL,
  refresh = NULL,
  ...
)

percolate_match(
  conn,
  index,
  type = NULL,
  body,
  routing = NULL,
  preference = NULL,
  ignore_unavailable = NULL,
  percolate_format = NULL,
  ...
)

percolate_list(conn, index, ...)

percolate_count(conn, index, type, body, ...)

percolate_delete(conn, index, id, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{Index name. Required}

\item{id}{A precolator id. Required}

\item{type}{Document type. Required}

\item{body}{Body json, or R list.}

\item{routing}{(character) In case the percolate queries are partitioned by a custom
routing value, that routing option makes sure that the percolate request only gets
executed on the shard where the routing value is partitioned to. This means that the
percolate request only gets executed on one shard instead of all shards. Multiple values
can be specified as a comma separated string, in that case the request can be be executed
on more than one shard.}

\item{preference}{(character) Controls which shard replicas are preferred to execute
the request on. Works the same as in the search API.}

\item{ignore_unavailable}{(logical) Controls if missing concrete indices should
silently be ignored. Same as is in the search API.}

\item{percolate_format}{(character) If ids is specified then the matches array in the
percolate response will contain a string array of the matching ids instead of an
array of objects. This can be useful to reduce the amount of data being send back to
the client. Obviously if there are two percolator queries with same id from different
indices there is no way to find out which percolator query belongs to what index. Any
other value to percolate_format will be ignored.}

\item{refresh}{If \code{TRUE} then refresh the affected shards to make this
operation visible to search, if "wait_for" then wait for a refresh to
make this operation visible to search, if \code{FALSE} (default) then do
nothing with refreshes. Valid choices: \code{TRUE}, \code{FALSE}, "wait_for"}

\item{...}{Curl options. Or in \code{percolate_list} function, further args
passed on to \code{\link[=Search]{Search()}}}
}
\description{
Store queries into an index then, via the percolate API, define
documents to retrieve these queries.
}
\details{
Additional body options, pass those in the body. These aren't query string
parameters:
\itemize{
\item filter - Reduces the number queries to execute during percolating. Only the
percolator queries that match with the filter will be included in the percolate
execution. The filter option works in near realtime, so a refresh needs to have
occurred for the filter to included the latest percolate queries.
\item query - Same as the filter option, but also the score is computed. The
computed scores can then be used by the track_scores and sort option.
\item size - Defines to maximum number of matches (percolate queries) to be returned.
Defaults to unlimited.
\item track_scores - Whether the _score is included for each match. The _score is
based on the query and represents how the query matched the percolate query's
metadata, not how the document (that is being percolated) matched the query. The query
option is required for this option. Defaults to false.
\item sort - Define a sort specification like in the search API. Currently only
sorting _score reverse (default relevancy) is supported. Other sort fields will
throw an exception. The size and query option are required for this setting. Like
track_score the score is based on the query and represents how the query matched
to the percolate query's metadata and not how the document being percolated matched
to the query.
\item aggs - Allows aggregation definitions to be included. The aggregations are
based on the matching percolator queries, look at the aggregation documentation on
how to define aggregations.
\item highlight - Allows highlight definitions to be included. The document being
percolated is being highlight for each matching query. This allows you to see how
each match is highlighting the document being percolated. See highlight documentation
on how to define highlights. The size option is required for highlighting, the
performance of highlighting in the percolate API depends of how many matches are
being highlighted.
}
}
\section{The Elasticsearch v5 split}{

In Elasticsearch < v5, there's a certain set of percolate APIs available,
while in Elasticsearch >= v5, there's a different set of APIs available.

Internally within these percolate functions we detect your Elasticsearch
version, then use the appropriate APIs
}

\examples{
\dontrun{
x <- connect(errors = "complete")

##### Elasticsearch < v5
if (x$es_ver() < 500) {
# typical usage
## create an index first
if (index_exists(x, "myindex")) index_delete(x, "myindex")
mapping <- '{
  "mappings": {
    "mytype": {
      "properties": {
        "message": {
           "type": "text"
        },
        "query": {
           "type": "percolator"
        }
      }
    }
  }
}'
index_create(x, "myindex", body = mapping)

## register a percolator
perc_body = '{
 "query" : {
    "match" : {
      "message" : "bonsai tree"
    }
 }
}'
percolate_register(x, index = "myindex", type = "mytype", 
  id = 1, body = perc_body)

## register another
perc_body2 <- '{
  "query" : {
    "match" : {
      "message" : "jane doe"
    }
  }
}'
percolate_register(x, index = "myindex", type = "mytype", 
  id = 2, body = perc_body2)

## match a document to a percolator
doc <- '{
 "query": {
   "percolate": {
     "field": "query",
     "document": {
       "message" : "A new bonsai tree in the office"
     }
   }
 }
}'
percolate_match(x, index = "myindex", type = "mytype", body = doc)

## List percolators - for an index, no type, can't do across indices
percolate_list(x, index = "myindex")$hits$hits

## Percolate counter
percolate_count(x, index = "myindex", type = "mytype", body = doc)$total

## delete a percolator
percolate_delete(x, index = "myindex", id = 2)
} # end ES < 5


##### Elasticsearch >= v5
if (x$es_ver() >= 500 && x$es_ver() <= 700) {
if (index_exists(x, "myindex")) index_delete(x, "myindex")

body <- '{
  "mappings": {
    "mytype": {
      "properties": {
        "message": {
           "type": "text"
        },
        "query": {
           "type": "percolator"
        }
      }
    }
  }
}'

# create the index with mapping
index_create(x, "myindex", body = body)

## register a percolator
z <- '{
  "query" : {
     "match" : {
       "message" : "bonsai tree"
     }
  }
}'
percolate_register(x, index = "myindex", type = "mytype", id = 1, body = z)

## register another
x2 <- '{
  "query" : {
    "match" : {
      "message" : "the office"
    }
  }
}'
percolate_register(x, index = "myindex", type = "mytype", id = 2, body = x2)

## match a document to a percolator
query <- '{
  "query" : {
    "percolate" : {
      "field": "query",
      "document": {
        "message": "A new bonsai tree in the office"
      }
    }
  }
}'
percolate_match(x, index = "myindex", body = query)
} # end ES >= 5




##### Elasticsearch >= v7
if (x$es_ver() >= 700) {
if (index_exists(x, "myindex")) index_delete(x, "myindex")

body <- '{
  "mappings": {
    "properties": {
      "message": {
        "type": "text"
      },
      "query": {
        "type": "percolator"
      }
    }
  }
}'

# create the index with mapping
index_create(x, "myindex", body = body)

## register a percolator
z <- '{
  "query" : {
     "match" : {
       "message" : "bonsai tree"
     }
  }
}'
percolate_register(x, index = "myindex", id = 1, body = z)

## register another
x2 <- '{
  "query" : {
    "match" : {
      "message" : "the office"
    }
  }
}'
percolate_register(x, index = "myindex", id = 2, body = x2)

## match a document to a percolator
query <- '{
  "query" : {
    "percolate" : {
      "field": "query",
      "document": {
        "message": "A new bonsai tree in the office"
      }
    }
  }
}'
percolate_match(x, index = "myindex", body = query)
} # end ES >= 7


}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/query-dsl-percolate-query.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/explain.R
\name{explain}
\alias{explain}
\title{Explain a search query.}
\usage{
explain(
  conn,
  index,
  id,
  type = NULL,
  source2 = NULL,
  fields = NULL,
  routing = NULL,
  parent = NULL,
  preference = NULL,
  source = NULL,
  q = NULL,
  df = NULL,
  analyzer = NULL,
  analyze_wildcard = NULL,
  lowercase_expanded_terms = NULL,
  lenient = NULL,
  default_operator = NULL,
  source_exclude = NULL,
  source_include = NULL,
  body = NULL,
  raw = FALSE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{Only one index. Required}

\item{id}{Document id, only one. Required}

\item{type}{Only one document type, optional}

\item{source2}{(logical) Set to TRUE to retrieve the _source of the document
explained. You can also retrieve part of the document by using
source_include & source_exclude (see Get API for more details). This
matches the \verb{_source} term, but we want to avoid the leading underscore.}

\item{fields}{Allows to control which stored fields to return as part of
the document explained.}

\item{routing}{Controls the routing in the case the routing was used during
indexing.}

\item{parent}{Same effect as setting the routing parameter.}

\item{preference}{Controls on which shard the explain is executed.}

\item{source}{Allows the data of the request to be put in the query string
of the url.}

\item{q}{The query string (maps to the query_string query).}

\item{df}{The default field to use when no field prefix is defined within
the query. Defaults to _all field.}

\item{analyzer}{The analyzer name to be used when analyzing the query
string. Defaults to the analyzer of the _all field.}

\item{analyze_wildcard}{(logical) Should wildcard and prefix queries be
analyzed or not. Default: \code{FALSE}}

\item{lowercase_expanded_terms}{Should terms be automatically lowercased
or not. Default: \code{TRUE}}

\item{lenient}{If set to true will cause format based failures (like
providing text to a numeric field) to be ignored. Default: \code{FALSE}}

\item{default_operator}{The default operator to be used, can be AND or OR.
Defaults to OR.}

\item{source_exclude}{A vector of fields to exclude from the returned
source2 field}

\item{source_include}{A vector of fields to extract and return from the
source2 field}

\item{body}{The query definition using the Query DSL. This is passed in the
body of the request.}

\item{raw}{If \code{TRUE} (default), data is parsed to list. If \code{FALSE}, then
raw JSON.}

\item{...}{Curl args passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Explain a search query.
}
\examples{
\dontrun{
(x <- connect())

explain(x, index = "plos", id = 14, q = "title:Germ")

body <- '{
 "query": {
   "match": { "title": "Germ" }
 }
}'
explain(x, index = "plos", id = 14, body=body)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/search-explain.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helpdocs.R
\name{units-time}
\alias{units-time}
\title{Time units}
\description{
Whenever durations need to be specified, eg for a timeout parameter, the duration can
be specified as a whole number representing time in milliseconds, or as a time value
like 2d for 2 days. The supported units are:
}
\details{
\tabular{ll}{
y \tab Year \cr
M \tab Month \cr
w \tab Week \cr
d \tab Day \cr
h \tab Hour \cr
m \tab Minute \cr
s \tab Second \cr
}
}
\seealso{
\link{units-distance}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_get.r
\name{docs_get}
\alias{docs_get}
\title{Get documents}
\usage{
docs_get(
  conn,
  index,
  id,
  type = NULL,
  source = NULL,
  fields = NULL,
  source_includes = NULL,
  source_excludes = NULL,
  exists = FALSE,
  raw = FALSE,
  callopts = list(),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) The name of the index. Required}

\item{id}{(numeric/character) The document ID. Can be numeric or character.
Required}

\item{type}{(character) The type of the document. optional}

\item{source}{(logical) If \code{TRUE} (default), return source. note that
it is actually set to \code{NULL} in the function definition, but within
Elasticsearch, it returns the source by default. alternatively,
you can pass a vector of field names to return.}

\item{fields}{Fields to return from the response object.}

\item{source_includes, source_excludes}{(character) fields to include in the
returned document, or to exclude. a character vector}

\item{exists}{(logical) Only return a logical as to whether the document
exists or not.}

\item{raw}{If \code{TRUE} (default), data is parsed to list. If \code{FALSE}, then raw
JSON.}

\item{callopts}{Curl args passed on to \link[crul:HttpClient]{crul::HttpClient}}

\item{verbose}{If TRUE (default) the url call used printed to console.}

\item{...}{Further args passed on to elastic search HTTP API as parameters.}
}
\description{
Get documents
}
\examples{
\dontrun{
(x <- connect())

if (!index_exists(x, "shakespeare")) {
  shakespeare <- system.file("examples", "shakespeare_data_.json",
    package = "elastic")
  shakespeare <- type_remover(shakespeare)
  invisible(docs_bulk(x, shakespeare))
}

docs_get(x, index='shakespeare', id=10)
docs_get(x, index='shakespeare', id=12)
docs_get(x, index='shakespeare', id=12, source=TRUE)

# Get certain fields
if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
  ### ES < v5
  docs_get(x, index='shakespeare', id=10, fields='play_name')
  docs_get(x, index='shakespeare', id=10, fields=c('play_name','speaker'))
} else {
  ### ES > v5
  docs_get(x, index='shakespeare', id=10, source='play_name')
  docs_get(x, index='shakespeare', id=10, source=c('play_name','speaker'))
}

# Just test for existence of the document
docs_get(x, index='plos', id=1, exists=TRUE)
docs_get(x, index='plos', id=123456, exists=TRUE)

# source includes / excludes
docs_get(x, index='shakespeare', id=10, source_includes = "play_name")
docs_get(x, index='shakespeare', id=10, source_excludes = "play_name")
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-get.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{elastic-defunct}
\alias{elastic-defunct}
\title{Defunct functions in elastic}
\description{
\itemize{
\item \code{\link[=mlt]{mlt()}}: The MLT API has been removed, use More Like This Query
via \code{\link[=Search]{Search()}}
\item \code{\link[=nodes_shutdown]{nodes_shutdown()}}: The _shutdown API has been removed. Instead,
setup Elasticsearch to run as a service (see Running as a Service on Linux
(\url{https://www.elastic.co/guide/en/elasticsearch/reference/2.0/setup-service.html}) or
Running as a Service on Windows
(\url{https://www.elastic.co/guide/en/elasticsearch/reference/2.0/setup-service-win.html}))
or use the -p command line option to write the PID to a file.
\item \code{\link[=index_status]{index_status()}}: _status route for the index API has been removed.
Replaced with the Indices Stats and Indices Recovery APIs.
\item \code{\link[=mapping_delete]{mapping_delete()}}: Elasticsearch dropped this route in their API. Instead
of deleting a mapping, delete the index and recreate with a new mapping.
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{info}
\alias{info}
\title{This function is defunct}
\usage{
info(...)
}
\description{
This function is defunct
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/reindex.R
\name{reindex}
\alias{reindex}
\title{Reindex}
\usage{
reindex(
  conn,
  body,
  refresh = NULL,
  requests_per_second = NULL,
  slices = NULL,
  timeout = NULL,
  wait_for_active_shards = NULL,
  wait_for_completion = NULL,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{body}{(list/character/json) The search definition using the Query DSL
and the prototype for the index request.}

\item{refresh}{(logical) Should the effected indexes be refreshed?}

\item{requests_per_second}{(integer) The throttle to set on this request in
sub-requests per second. - 1 means no throttle. Default: 0}

\item{slices}{(integer) The number of slices this task should be divided
into. Defaults to 1 meaning the task isn't sliced into subtasks. Default: 1}

\item{timeout}{(character) Time each individual bulk request should wait
for shards that are unavailable. Default: '1m'}

\item{wait_for_active_shards}{(integer) Sets the number of shard copies that
must be active before proceeding with the reindex operation. Defaults to 1,
meaning the primary shard only. Set to all for all shard copies, otherwise
set to any non-negative value less than or equal to the total number of
copies for the shard (number of replicas + 1)}

\item{wait_for_completion}{(logical) Should the request block until the
reindex is complete? Default: \code{TRUE}}

\item{...}{Curl options, passed on to \link[crul:verb-POST]{crul::verb-POST}}
}
\description{
Reindex all documents from one index to another.
}
\examples{
\dontrun{
x <- connect()

if (!index_exists(x, "twitter")) index_create(x, "twitter")
if (!index_exists(x, "new_twitter")) index_create(x, "new_twitter")
body <- '{
  "source": {
    "index": "twitter"
  },
  "dest": {
    "index": "new_twitter"
  }
}'
reindex(x, body = body)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-reindex.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_mget.r
\name{docs_mget}
\alias{docs_mget}
\title{Get multiple documents via the multiple get API}
\usage{
docs_mget(
  conn,
  index = NULL,
  type = NULL,
  ids = NULL,
  type_id = NULL,
  index_type_id = NULL,
  source = NULL,
  fields = NULL,
  raw = FALSE,
  callopts = list(),
  verbose = TRUE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{Index. Required.}

\item{type}{Document type. Required.}

\item{ids}{More than one document id, see examples.}

\item{type_id}{List of vectors of length 2, each with an element for
type and id.}

\item{index_type_id}{List of vectors of length 3, each with an element for
index, type, and id.}

\item{source}{(logical) If \code{TRUE}, return source.}

\item{fields}{Fields to return from the response object.}

\item{raw}{If TRUE (default), data is parsed to list. If FALSE, then raw JSON.}

\item{callopts}{Curl args passed on to \code{\link[crul]{HttpClient}}}

\item{verbose}{If TRUE (default) the url call used printed to console.}

\item{...}{Further args passed on to elastic search HTTP API as parameters.}
}
\description{
Get multiple documents via the multiple get API
}
\details{
You can pass in one of three combinations of parameters:
\itemize{
\item Pass in something for \code{index}, \code{type}, and \code{id}.
This is the simplest, allowing retrieval from the same index, same type,
and many ids.
\item Pass in only \code{index} and \code{type_id} - this allows you to
get multiple documents from the same index, but from different types.
\item Pass in only \code{index_type_id} - this is so that you can get
multiple documents from different indexes and different types.
}
}
\examples{
\dontrun{
(x <- connect())

if (!index_exists(x, 'plos')) {
  plosdat <- system.file("examples", "plos_data.json",
    package = "elastic")
  plosdat <- type_remover(plosdat)
  invisible(docs_bulk(x, plosdat))
}

# same index, many ids
docs_mget(x, index="plos", ids=c(9,10))

# Same index and type
docs_mget(x, index="plos", type="_doc", ids=c(9,10))

tmp <- docs_mget(x, index="plos", ids=c(9, 10), raw=TRUE)
es_parse(tmp)
docs_mget(x, index="plos", ids=c(9, 10), source='title')
docs_mget(x, index="plos", ids=c(14, 19), source=TRUE)

# curl options
docs_mget(x, index="plos", ids=1:2, callopts=list(verbose=TRUE))

# Same index, but different types
if (index_exists(x, 'shakespeare')) index_delete(x, 'shakespeare')
shakedat <- system.file("examples", "shakespeare_data.json",
  package = "elastic")
invisible(docs_bulk(x, shakedat))

docs_mget(x, index="shakespeare", type_id=list(c("scene",1), c("line",20)))
docs_mget(x, index="shakespeare", type_id=list(c("scene",1), c("line",20)),
  source='play_name')

# Different indices and different types pass in separately
docs_mget(x, index_type_id = list(
  c("shakespeare", "line", 20),
  c("plos", "article", 1)
 )
)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-multi-get.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/defunct.R
\name{mapping_delete}
\alias{mapping_delete}
\title{Mapping delete}
\usage{
mapping_delete(...)
}
\description{
Mapping delete
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_body.R
\name{search_body}
\alias{search_body}
\title{Full text search of Elasticsearch - body requests.}
\usage{
search_body(
  conn,
  index = NULL,
  type = NULL,
  raw = FALSE,
  callopts = list(),
  query = list(),
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{Index name}

\item{type}{Document type}

\item{raw}{If \code{TRUE} (default), data is parsed to list. If \code{FALSE}, then
raw JSON.}

\item{callopts}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}

\item{query}{Query, either a list or json.}

\item{...}{Further args passed on to elastic search HTTP API as parameters.
Not used right now.}
}
\description{
Full text search of Elasticsearch - body requests.
}
\examples{
\dontrun{
# connection setup
# x <- connect()
# x$ping()

# pass in as an R list
aggs <- list(aggs = list(stats = list(terms = list(field = "text_entry"))))
# search_body(x, index="shakespeare", query=aggs)

# or pass in as json query with newlines, easy to read
aggs <- '{
    "aggs": {
        "stats" : {
            "terms" : {
                "field" : "text_entry"
            }
        }
    }
}'
# search_body(x, index="shakespeare", query=aggs)


# or pass in collapsed json string
aggs <- '{"aggs":{"stats":{"terms":{"field":"text_entry"}}}}'
# search_body(x, index="shakespeare", query=aggs)

# match query
match <- '{"query": {"match" : {"text_entry" : "Two Gentlemen"}}}'
# search_body(x, index="shakespeare", query=match)

# multi-match (multiple fields that is) query
mmatch <- '{"query": {"multi_match" : {"query" : "henry", "fields": 
["text_entry","play_name"]}}}'
# search_body(x, index="shakespeare", query=mmatch)

# bool query
mmatch <- '{
 "query": {
   "bool" : {
     "must_not" : {
       "range" : {
         "speech_number" : {
           "from" : 1, "to": 5
}}}}}}'
# search_body(x, index="shakespeare", query=mmatch)

# Boosting query
boost <- '{
 "query" : {
  "boosting" : {
      "positive" : {
          "term" : {
              "play_name" : "henry"
          }
      },
      "negative" : {
          "term" : {
              "text_entry" : "thou"
          }
      },
      "negative_boost" : 0.2
    }
 }
}'
# search_body(x, index="shakespeare", query=mmatch)
}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_bulk_create.R
\name{docs_bulk_create}
\alias{docs_bulk_create}
\title{Use the bulk API to create documents}
\usage{
docs_bulk_create(
  conn,
  x,
  index = NULL,
  type = NULL,
  chunk_size = 1000,
  doc_ids = NULL,
  es_ids = TRUE,
  raw = FALSE,
  quiet = FALSE,
  query = list(),
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{x}{A list, data.frame, or character path to a file. required.}

\item{index}{(character) The index name to use. Required for data.frame
input, but optional for file inputs.}

\item{type}{(character) The type. default: \code{NULL}. Note that \code{type} is
deprecated in Elasticsearch v7 and greater, and removed in Elasticsearch v8}

\item{chunk_size}{(integer) Size of each chunk. If your data.frame is smaller
thank \code{chunk_size}, this parameter is essentially ignored. We write in
chunks because at some point, depending on size of each document, and
Elasticsearch setup, writing a very large number of documents in one go
becomes slow, so chunking can help. This parameter is ignored if you
pass a file name. Default: 1000}

\item{doc_ids}{An optional vector (character or numeric/integer) of document
ids to use. This vector has to equal the size of the documents you are
passing in, and will error if not. If you pass a factor we convert to
character. Default: not passed}

\item{es_ids}{(boolean) Let Elasticsearch assign document IDs as UUIDs.
These are sequential, so there is order to the IDs they assign.
If \code{TRUE}, \code{doc_ids} is ignored. Default: \code{TRUE}}

\item{raw}{(logical) Get raw JSON back or not. If \code{TRUE}
you get JSON; if \code{FALSE} you get a list. Default: \code{FALSE}}

\item{quiet}{(logical) Suppress progress bar. Default: \code{FALSE}}

\item{query}{(list) a named list of query parameters. optional.
options include: pipeline, refresh, routing, _source, _source_excludes,
_source_includes, timeout, wait_for_active_shards. See the docs bulk
ES page for details}

\item{...}{Pass on curl options to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Use the bulk API to create documents
}
\details{
For doing create with a file already prepared for the bulk API,
see \code{\link[=docs_bulk]{docs_bulk()}}

Only data.frame's are supported for now.
}
\examples{
\dontrun{
x <- connect()
if (index_exists(x, "foobar")) index_delete(x, "foobar")

df <- data.frame(name = letters[1:3], size = 1:3, id = 100:102)
docs_bulk_create(x, df, 'foobar', es_ids = FALSE)
Search(x, "foobar", asdf = TRUE)$hits$hits

# more examples
docs_bulk_create(x, mtcars, index = "hello")
## field names cannot contain dots
names(iris) <- gsub("\\\\.", "_", names(iris))
docs_bulk_create(x, iris, "iris")
## type can be missing, but index can not
docs_bulk_create(x, iris, "flowers")
## big data.frame, 53K rows, load ggplot2 package first
# res <- docs_bulk_create(x, diamonds, "diam")
# Search(x, "diam")$hits$total$value
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-bulk.html}
}
\seealso{
Other bulk-functions: 
\code{\link{docs_bulk_delete}()},
\code{\link{docs_bulk_index}()},
\code{\link{docs_bulk_prep}()},
\code{\link{docs_bulk_update}()},
\code{\link{docs_bulk}()}
}
\concept{bulk-functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_bulk.r
\name{docs_bulk}
\alias{docs_bulk}
\title{Use the bulk API to create, index, update, or delete documents.}
\usage{
docs_bulk(
  conn,
  x,
  index = NULL,
  type = NULL,
  chunk_size = 1000,
  doc_ids = NULL,
  es_ids = TRUE,
  raw = FALSE,
  quiet = FALSE,
  query = list(),
  digits = NA,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{x}{A list, data.frame, or character path to a file. required.}

\item{index}{(character) The index name to use. Required for data.frame
input, but optional for file inputs.}

\item{type}{(character) The type. default: \code{NULL}. Note that \code{type} is
deprecated in Elasticsearch v7 and greater, and removed in Elasticsearch v8}

\item{chunk_size}{(integer) Size of each chunk. If your data.frame is smaller
thank \code{chunk_size}, this parameter is essentially ignored. We write in
chunks because at some point, depending on size of each document, and
Elasticsearch setup, writing a very large number of documents in one go
becomes slow, so chunking can help. This parameter is ignored if you
pass a file name. Default: 1000}

\item{doc_ids}{An optional vector (character or numeric/integer) of document
ids to use. This vector has to equal the size of the documents you are
passing in, and will error if not. If you pass a factor we convert to
character. Default: not passed}

\item{es_ids}{(boolean) Let Elasticsearch assign document IDs as UUIDs.
These are sequential, so there is order to the IDs they assign.
If \code{TRUE}, \code{doc_ids} is ignored. Default: \code{TRUE}}

\item{raw}{(logical) Get raw JSON back or not. If \code{TRUE}
you get JSON; if \code{FALSE} you get a list. Default: \code{FALSE}}

\item{quiet}{(logical) Suppress progress bar. Default: \code{FALSE}}

\item{query}{(list) a named list of query parameters. optional.
options include: pipeline, refresh, routing, _source, _source_excludes,
_source_includes, timeout, wait_for_active_shards. See the docs bulk
ES page for details}

\item{digits}{digits used by the parameter of the same name by
\code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}} to convert data to JSON before being submitted to
your ES instance. default: \code{NA}}

\item{...}{Pass on curl options to \link[crul:HttpClient]{crul::HttpClient}}
}
\value{
A list
}
\description{
Use the bulk API to create, index, update, or delete documents.
}
\details{
More on the Bulk API:
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-bulk.html}

This function dispatches on data.frame or character input. Character input
has to be a file name or the function stops with an error message.

If you pass a data.frame to this function, we by default do an index
operation, that is, create the record in the index given by those
parameters to the function. Down the road perhaps we will try to support
other operations on the bulk API. if you pass a file, of course in that
file, you can specify any operations you want.

Row names are dropped from data.frame, and top level names for a list
are dropped as well.

A progress bar gives the progress for data.frames and lists - the progress
bar is based around a for loop, where progress indicates progress along
the iterations of the for loop, where each iteration is a chunk of data
that's converted to bulk format, then pushed into Elasticsearch. The
\code{character} method has no for loop, so no progress bar.
}
\section{Document IDs}{

Document IDs can be passed in via the \code{doc_ids} paramater when passing
in data.frame or list, but not with files. If ids are not passed to
\code{doc_ids}, we assign document IDs from 1 to length of the object
(rows of a data.frame, or length of a list). In the future we may allow the
user to select whether they want to assign sequential numeric IDs or
to allow Elasticsearch to assign IDs, which are UUIDs that are actually
sequential, so you still can determine an order of your documents.
}

\section{Document IDs and Factors}{

If you pass in ids that are of class factor, we coerce them to character
with \code{as.character}. This applies to both data.frame and list inputs, but
not to file inputs.
}

\section{Large numbers for document IDs}{

Until recently, if you had very large integers for document IDs,
\code{docs_bulk} failed. It should be fixed now. Let us know if not.
}

\section{Missing data}{

As of \pkg{elastic} version \verb{0.7.8.9515} we convert \code{NA} to
\code{null} before loading into Elasticsearch. Previously, fields that
had an \code{NA} were dropped - but when you read data back from
Elasticsearch into R, you retain those missing values as \pkg{jsonlite}
fills those in for you. Now, fields with \code{NA}'s are made into
\code{null}, and are not dropped in Elasticsearch.

Note also that null values can not be indexed or searched
\url{https://www.elastic.co/guide/en/elasticsearch/reference/5.3/null-value.html}
}

\section{Tips}{

This function returns the response from Elasticsearch, but you'll likely
not be that interested in the response. If not, wrap your call to
\code{docs_bulk} in \code{\link[=invisible]{invisible()}}, like so: \code{invisible(docs_bulk(...))}
}

\section{Connections/Files}{

We create temporary files, and connections to those files, when data.frame's
and lists are passed in to \code{docs_bulk()} (not when a file is passed in
since we don't need to create a file). After inserting data into your
Elasticsearch instance, we close the connections and delete the temporary files.

There are some exceptions though. When you pass in your own file, whether a
tempfile or not, we don't delete those files after using them - in case
you need those files again. Your own tempfile's will be cleaned up/delete
when the R session ends. Non-tempfile's won't be cleaned up/deleted after
the R session ends.
}

\section{Elasticsearch versions that don't support type}{

See the \code{\link[=type_remover]{type_remover()}} function.
}

\examples{
\dontrun{
# connection setup
(x <- connect())

# From a file already in newline delimited JSON format
plosdat <- system.file("examples", "plos_data.json", package = "elastic")
docs_bulk(x, plosdat)
aliases_get(x)
index_delete(x, index='plos')
aliases_get(x)

# From a data.frame
docs_bulk(x, mtcars, index = "hello")
## field names cannot contain dots
names(iris) <- gsub("\\\\.", "_", names(iris))
docs_bulk(x, iris, "iris")
## type can be missing, but index can not
docs_bulk(x, iris, "flowers")
## big data.frame, 53K rows, load ggplot2 package first
# res <- docs_bulk(x, diamonds, "diam")
# Search(x, "diam")$hits$total

# From a list
docs_bulk(x, apply(iris, 1, as.list), index="iris")
docs_bulk(x, apply(USArrests, 1, as.list), index="arrests")
# dim_list <- apply(diamonds, 1, as.list)
# out <- docs_bulk(x, dim_list, index="diamfromlist")

# When using in a loop
## We internally get last _id counter to know where to start on next bulk
## insert but you need to sleep in between docs_bulk calls, longer the
## bigger the data is
files <- c(system.file("examples", "test1.csv", package = "elastic"),
           system.file("examples", "test2.csv", package = "elastic"),
           system.file("examples", "test3.csv", package = "elastic"))
for (i in seq_along(files)) {
  d <- read.csv(files[[i]])
  docs_bulk(x, d, index = "testes")
  Sys.sleep(1)
}
count(x, "testes")
index_delete(x, "testes")

# You can include your own document id numbers
## Either pass in as an argument
index_create(x, "testes")
files <- c(system.file("examples", "test1.csv", package = "elastic"),
           system.file("examples", "test2.csv", package = "elastic"),
           system.file("examples", "test3.csv", package = "elastic"))
tt <- vapply(files, function(z) NROW(read.csv(z)), numeric(1))
ids <- list(1:tt[1],
           (tt[1] + 1):(tt[1] + tt[2]),
           (tt[1] + tt[2] + 1):sum(tt))
for (i in seq_along(files)) {
  d <- read.csv(files[[i]])
  docs_bulk(x, d, index = "testes", doc_ids = ids[[i]],
    es_ids = FALSE)
}
count(x, "testes")
index_delete(x, "testes")

## or include in the input data
### from data.frame's
index_create(x, "testes")
files <- c(system.file("examples", "test1_id.csv", package = "elastic"),
           system.file("examples", "test2_id.csv", package = "elastic"),
           system.file("examples", "test3_id.csv", package = "elastic"))
readLines(files[[1]])
for (i in seq_along(files)) {
  d <- read.csv(files[[i]])
  docs_bulk(x, d, index = "testes")
}
count(x, "testes")
index_delete(x, "testes")

### from lists via file inputs
index_create(x, "testes")
for (i in seq_along(files)) {
  d <- read.csv(files[[i]])
  d <- apply(d, 1, as.list)
  docs_bulk(x, d, index = "testes")
}
count(x, "testes")
index_delete(x, "testes")

# data.frame's with a single column
## this didn't use to work, but now should work
db <- paste0(sample(letters, 10), collapse = "")
index_create(x, db)
res <- data.frame(foo = 1:10)
out <- docs_bulk(x, res, index = db)
count(x, db)
index_delete(x, db)


# data.frame with a mix of actions
## make sure you use a column named 'es_action' or this won't work
## if you need to delete or update you need document IDs
if (index_exists(x, "baz")) index_delete(x, "baz")
df <- data.frame(a = 1:5, b = 6:10, c = letters[1:5], stringsAsFactors = FALSE) 
invisible(docs_bulk(x, df, "baz"))
Sys.sleep(3)
(res <- Search(x, 'baz', asdf=TRUE)$hits$hits)
df[1, "a"] <- 99
df[1, "c"] <- "aa"
df[3, "c"] <- 33
df[3, "c"] <- "cc"
df$es_action <- c('update', 'delete', 'update', 'delete', 'delete')
df$id <- res$`_id`
df
invisible(docs_bulk(x, df, "baz", es_ids = FALSE))
### or es_ids = FALSE and pass in document ids to doc_ids
# invisible(docs_bulk(df, "baz", es_ids = FALSE, doc_ids = df$id))
Search(x, 'baz', asdf=TRUE)$hits$hits


# Curl options
plosdat <- system.file("examples", "plos_data.json",
  package = "elastic")
plosdat <- type_remover(plosdat)
invisible(docs_bulk(x, plosdat, verbose = TRUE))


# suppress progress bar
invisible(docs_bulk(x, mtcars, index = "hello", quiet = TRUE))
## vs. 
invisible(docs_bulk(x, mtcars, index = "hello", quiet = FALSE))
}
}
\seealso{
Other bulk-functions: 
\code{\link{docs_bulk_create}()},
\code{\link{docs_bulk_delete}()},
\code{\link{docs_bulk_index}()},
\code{\link{docs_bulk_prep}()},
\code{\link{docs_bulk_update}()}
}
\concept{bulk-functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_template.R
\name{Search_template}
\alias{Search_template}
\alias{Search_template_register}
\alias{Search_template_get}
\alias{Search_template_delete}
\alias{Search_template_render}
\title{Search or validate templates}
\usage{
Search_template(conn, body = list(), raw = FALSE, ...)

Search_template_register(conn, template, body = list(), raw = FALSE, ...)

Search_template_get(conn, template, ...)

Search_template_delete(conn, template, ...)

Search_template_render(conn, body = list(), raw = FALSE, ...)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{body}{Query, either a list or json.}

\item{raw}{(logical) If \code{FALSE} (default), data is parsed to list.
If \code{TRUE}, then raw JSON returned}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}

\item{template}{(character) a template name}
}
\description{
Search or validate templates
}
\section{Template search}{

With \code{Search_template} you can search with a template, using
mustache templating. Added in Elasticsearch v1.1
}

\section{Template render}{

With \code{Search_template_render} you validate a template without
conducting the search. Added in Elasticsearch v2.0
}

\section{Pre-registered templates}{

Register a template with \code{Search_template_register}. You can get
the template with \code{Search_template_get} and delete the template
with \code{Search_template_delete}

You can also pre-register search templates by storing them in the
\code{config/scripts} directory, in a file using the .mustache
extension. In order to execute the stored template, reference it
by it's name under the template key, like
\verb{"file": "templateName", ...}
}

\examples{
\dontrun{
# connection setup
(x <- connect())

if (!index_exists(x, "iris")) {
  invisible(docs_bulk(x, iris, "iris"))
}

body1 <- '{
  "inline" : {
    "query": { "match" : { "{{my_field}}" : "{{my_value}}" } },
    "size" : "{{my_size}}"
  },
  "params" : {
    "my_field" : "Species",
    "my_value" : "setosa",
    "my_size" : 3
  }
}'
Search_template(x, body = body1)

body2 <- '{
 "inline": {
   "query": {
      "match": {
          "Species": "{{query_string}}"
      }
   }
 },
 "params": {
   "query_string": "versicolor"
 }
}'
Search_template(x, body = body2)

# pass in a list
mylist <- list(
  inline = list(query = list(match = list(`{{my_field}}` = "{{my_value}}"))),
  params = list(my_field = "Species", my_value = "setosa", my_size = 3L)
)
Search_template(x, body = mylist)

## Validating templates w/ Search_template_render()
Search_template_render(x, body = body1)
Search_template_render(x, body = body2)

## pre-registered templates
### register a template
if (x$es_ver() <= 520) {
  body3 <- '{
    "template": {
       "query": {
           "match": {
               "Species": "{{query_string}}"
           }
       }
     }
  }'
  Search_template_register(x, 'foobar', body = body3)
} else {
  body3 <- '{
   "script": {
     "lang": "mustache",
       "source": {
         "query": {
           "match": {
             "Species": "{{query_string}}"
           }
         }
       }
     }
  }'
  Search_template_register(x, 'foobar', body = body3)
}

### get template
Search_template_get(x, 'foobar')

### use the template
body4 <- '{
 "id": "foobar",
  	"params": {
      "query_string": "setosa"
  }
}'
Search_template(x, body = body4)

### delete the template
Search_template_delete(x, 'foobar')
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/search-template.html}
}
\seealso{
\code{\link[=Search]{Search()}}, \code{\link[=Search_uri]{Search_uri()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/searchapis.R
\name{searchapis}
\alias{searchapis}
\title{Overview of search functions}
\description{
Overview of search functions
}
\details{
Elasticsearch search APIs include the following functions:
\itemize{
\item \code{\link[=Search]{Search()}} - Search using the Query DSL via the body of the request.
\item \code{\link[=Search_uri]{Search_uri()}} - Search using the URI search API only. This may be
needed for servers that block POST requests for security, or maybe you don't need
complicated requests, in which case URI only requests are suffice.
\item \code{\link[=msearch]{msearch()}} - Multi Search - execute several search requests defined
in a file passed to \code{msearch}
\item \code{\link[=search_shards]{search_shards()}} - Search shards.
\item \code{\link[=count]{count()}} - Get counts for various searches.
\item \code{\link[=explain]{explain()}} - Computes a score explanation for a query and a specific
document. This can give useful feedback whether a document matches or didn't match
a specific query.
\item \code{\link[=validate]{validate()}} - Validate a search
\item \code{\link[=field_stats]{field_stats()}} - Search field statistics
\item \code{\link[=percolate]{percolate()}} - Store queries into an index then, via the percolate API,
define documents to retrieve these queries.
}

More will be added soon.
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/search.html}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/termvectors.R
\name{termvectors}
\alias{termvectors}
\title{Termvectors}
\usage{
termvectors(
  conn,
  index,
  type = NULL,
  id = NULL,
  body = list(),
  pretty = TRUE,
  field_statistics = TRUE,
  fields = NULL,
  offsets = TRUE,
  parent = NULL,
  payloads = TRUE,
  positions = TRUE,
  realtime = TRUE,
  preference = "random",
  routing = NULL,
  term_statistics = FALSE,
  version = NULL,
  version_type = NULL,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{index}{(character) The index in which the document resides.}

\item{type}{(character) The type of the document. optional}

\item{id}{(character) The id of the document, when not specified a doc
param should be supplied.}

\item{body}{(character) Define parameters and or supply a document to get
termvectors for}

\item{pretty}{(logical) pretty print. Default: \code{TRUE}}

\item{field_statistics}{(character) Specifies if document count, sum
of document frequencies and sum of total term frequencies should be
returned. Default: \code{TRUE}}

\item{fields}{(character) A comma-separated list of fields to return.}

\item{offsets}{(character) Specifies if term offsets should be returned.
Default: \code{TRUE}}

\item{parent}{(character) Parent id of documents.}

\item{payloads}{(character) Specifies if term payloads should be returned.
Default: \code{TRUE}}

\item{positions}{(character) Specifies if term positions should be returned.
Default: \code{TRUE}}

\item{realtime}{(character) Specifies if request is real-time as opposed to
near-real-time (Default: \code{TRUE}).}

\item{preference}{(character) Specify the node or shard the operation
should be performed on (Default: \code{random}).}

\item{routing}{(character) Specific routing value.}

\item{term_statistics}{(character) Specifies if total term frequency and
document frequency should be returned. Default: \code{FALSE}}

\item{version}{(character) Explicit version number for concurrency control}

\item{version_type}{(character) Specific version type, valid choices are:
'internal', 'external', 'external_gte', 'force'}

\item{...}{Curl args passed on to \link[crul:verb-POST]{crul::verb-POST}}
}
\description{
Termvectors
}
\details{
Returns information and statistics on terms in the fields of a
particular document. The document could be stored in the index or
artificially provided by the user (Added in 1.4). Note that for
documents stored in the index, this is a near realtime API as the term
vectors are not available until the next refresh.
}
\examples{
\dontrun{
x <- connect()

if (!index_exists(x, 'plos')) {
  plosdat <- system.file("examples", "plos_data.json",
    package = "elastic")
  plosdat <- type_remover(plosdat)
  invisible(docs_bulk(x, plosdat))
}
if (!index_exists(x, 'omdb')) {
  omdb <- system.file("examples", "omdb.json", package = "elastic")
  omdb <- type_remover(omdb)
  invisible(docs_bulk(x, omdb))
}

body <- '{
  "fields" : ["title"],
  "offsets" : true,
  "positions" : true,
  "term_statistics" : true,
  "field_statistics" : true
}'
termvectors(x, 'plos', id = 29, body = body)

body <- '{
  "fields" : ["Plot"],
  "offsets" : true,
  "positions" : true,
  "term_statistics" : true,
  "field_statistics" : true
}'
termvectors(x, 'omdb', id = Search(x, "omdb", size=1)$hits$hits[[1]]$`_id`,
body = body)
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-termvectors.html}
}
\seealso{
\code{\link[=mtermvectors]{mtermvectors()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docs_bulk_index.R
\name{docs_bulk_index}
\alias{docs_bulk_index}
\title{Use the bulk API to index documents}
\usage{
docs_bulk_index(
  conn,
  x,
  index = NULL,
  type = NULL,
  chunk_size = 1000,
  doc_ids = NULL,
  es_ids = TRUE,
  raw = FALSE,
  quiet = FALSE,
  query = list(),
  digits = NA,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{x}{A list, data.frame, or character path to a file. required.}

\item{index}{(character) The index name to use. Required for data.frame
input, but optional for file inputs.}

\item{type}{(character) The type. default: \code{NULL}. Note that \code{type} is
deprecated in Elasticsearch v7 and greater, and removed in Elasticsearch v8}

\item{chunk_size}{(integer) Size of each chunk. If your data.frame is smaller
thank \code{chunk_size}, this parameter is essentially ignored. We write in
chunks because at some point, depending on size of each document, and
Elasticsearch setup, writing a very large number of documents in one go
becomes slow, so chunking can help. This parameter is ignored if you
pass a file name. Default: 1000}

\item{doc_ids}{An optional vector (character or numeric/integer) of document
ids to use. This vector has to equal the size of the documents you are
passing in, and will error if not. If you pass a factor we convert to
character. Default: not passed}

\item{es_ids}{(boolean) Let Elasticsearch assign document IDs as UUIDs.
These are sequential, so there is order to the IDs they assign.
If \code{TRUE}, \code{doc_ids} is ignored. Default: \code{TRUE}}

\item{raw}{(logical) Get raw JSON back or not. If \code{TRUE}
you get JSON; if \code{FALSE} you get a list. Default: \code{FALSE}}

\item{quiet}{(logical) Suppress progress bar. Default: \code{FALSE}}

\item{query}{(list) a named list of query parameters. optional.
options include: pipeline, refresh, routing, _source, _source_excludes,
_source_includes, timeout, wait_for_active_shards. See the docs bulk
ES page for details}

\item{digits}{digits used by the parameter of the same name by
\code{\link[jsonlite:fromJSON]{jsonlite::toJSON()}} to convert data to JSON before being submitted to
your ES instance. default: \code{NA}}

\item{...}{Pass on curl options to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Use the bulk API to index documents
}
\details{
For doing index with a file already prepared for the bulk API,
see \code{\link[=docs_bulk]{docs_bulk()}}

Only data.frame's are supported for now.
}
\examples{
\dontrun{
x <- connect()
if (index_exists(x, "foobar")) index_delete(x, "foobar")

df <- data.frame(name = letters[1:3], size = 1:3, id = 100:102)
docs_bulk_index(x, df, 'foobar')
docs_bulk_index(x, df, 'foobar', es_ids = FALSE)
Search(x, "foobar", asdf = TRUE)$hits$hits

# more examples
docs_bulk_index(x, mtcars, index = "hello")
## field names cannot contain dots
names(iris) <- gsub("\\\\.", "_", names(iris))
docs_bulk_index(x, iris, "iris")
## type can be missing, but index can not
docs_bulk_index(x, iris, "flowers")
## big data.frame, 53K rows, load ggplot2 package first
# res <- docs_bulk_index(x, diamonds, "diam")
# Search(x, "diam")$hits$total$value
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/current/docs-bulk.html}
}
\seealso{
Other bulk-functions: 
\code{\link{docs_bulk_create}()},
\code{\link{docs_bulk_delete}()},
\code{\link{docs_bulk_prep}()},
\code{\link{docs_bulk_update}()},
\code{\link{docs_bulk}()}
}
\concept{bulk-functions}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/documents.R
\name{documents}
\alias{documents}
\title{Elasticsearch documents functions.}
\description{
Elasticsearch documents functions.
}
\details{
There are five functions to work directly with documents.
\itemize{
\item \code{\link[=docs_get]{docs_get()}}
\item \code{\link[=docs_mget]{docs_mget()}}
\item \code{\link[=docs_create]{docs_create()}}
\item \code{\link[=docs_delete]{docs_delete()}}
\item \code{\link[=docs_bulk]{docs_bulk()}}
}
}
\examples{
\dontrun{
# Get a document
# docs_get(index='plos', type='article', id=1)

# Get multiple documents
# docs_mget(index="shakespeare", type="line", id=c(9,10))

# Create a document
# docs_create(index='plos', type='article', id=35, body=list(id="12345", title="New title"))

# Delete a document
# docs_delete(index='plos', type='article', id=35)

# Bulk load documents
# plosdat <- system.file("examples", "plos_data.json", package = "elastic")
# docs_bulk(plosdat)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/field_stats.R
\name{field_stats}
\alias{field_stats}
\title{Search field statistics}
\usage{
field_stats(
  conn,
  fields = NULL,
  index = NULL,
  level = "cluster",
  body = list(),
  raw = FALSE,
  asdf = FALSE,
  ...
)
}
\arguments{
\item{conn}{an Elasticsearch connection object, see \code{\link[=connect]{connect()}}}

\item{fields}{A list of fields to compute stats for. optional}

\item{index}{Index name, one or more}

\item{level}{Defines if field stats should be returned on a per index level
or on a cluster wide level. Valid values are 'indices' and 'cluster'
(default)}

\item{body}{Query, either a list or json}

\item{raw}{(logical) Get raw JSON back or not}

\item{asdf}{(logical) If \code{TRUE}, use \code{\link[jsonlite]{fromJSON}}
to parse JSON directly to a data.frame. If \code{FALSE} (Default), list
output is given.}

\item{...}{Curl args passed on to \link[crul:HttpClient]{crul::HttpClient}}
}
\description{
Search field statistics
}
\details{
The field stats api allows you to get statistical properties of a
field without executing a search, but looking up measurements that are
natively available in the Lucene index. This can be useful to explore a
dataset which you don't know much about. For example, this allows creating
a histogram aggregation with meaningful intervals based on the min/max range
of values.

The field stats api by defaults executes on all indices, but can execute on
specific indices too.
}
\note{
Deprecated in Elasticsearch versions equal to/greater than 5.4.0
}
\examples{
\dontrun{
x <- connect()

if (gsub("\\\\.", "", x$ping()$version$number) < 500) {
  field_stats(x, body = '{ "fields": ["speaker"] }', index = "shakespeare")
  ff <- c("scientificName", "continent", "decimalLatitude", "play_name", 
    "speech_number")
  field_stats(x, "play_name")
  field_stats(x, "play_name", level = "cluster")
  field_stats(x, ff, level = "indices")
  field_stats(x, ff)
  field_stats(x, ff, index = c("gbif", "shakespeare"))

  # can also pass a body, just as with Search()
  # field_stats(x, body = list(fields = "rating")) # doesn't work
  field_stats(x, body = '{ "fields": ["scientificName"] }', index = "gbif")

  body <- '{
    "fields" : ["scientificName", "decimalLatitude"]
  }'
  field_stats(x, body = body, level = "indices", index = "gbif")
}
}
}
\references{
\url{https://www.elastic.co/guide/en/elasticsearch/reference/5.6/search-field-stats.html}
}
\seealso{
\code{\link[=field_caps]{field_caps()}}
}
