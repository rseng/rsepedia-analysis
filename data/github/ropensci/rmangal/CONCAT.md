# rmangal :package: - an R Client for the Mangal database <img src="man/figures/logo.png" width="130" align="right"/>

[![](https://badges.ropensci.org/332_status.svg)](https://github.com/ropensci/software-review/issues/332)
[![R CMD Check](https://github.com/ropensci/rmangal/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/ropensci/rmangal/actions/workflows/R-CMD-check.yaml)
[![codecov](https://app.codecov.io/gh/ropensci/rmangal/branch/master/graph/badge.svg?token=lGqUVLM2o3)](https://app.codecov.io/gh/ropensci/ropensci/rmangal)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)
[![CRAN status](https://www.r-pkg.org:443/badges/version/rmangal)](https://CRAN.R-project.org/package=rmangal)


## Context

[Mangal](https://mangal.io/#/) -- a global ecological interactions database --
serializes ecological interaction matrices into nodes (e.g. taxon, individuals
or population) and interactions (i.e. edges). For each network, Mangal offers
the opportunity to store study context such as the location, sampling
environment, inventory date and informations pertaining to the original
publication. For every nodes involved in the ecological networks, Mangal
references unique taxonomic identifiers such as Encyclopedia of Life (EOL),
Catalogue of Life (COL), Global Biodiversity Information Facility (GBIF) etc.
and can extend nodes informations to individual traits.

**rmangal** is an R client to the Mangal database and provides various functions
to explore his content through search functions. It offers methods to retrieve
networks structured as `mgNetwork` or `mgNetworksCollection` S3 objects and
methods to convert `mgNetwork` to other class objects in order to analyze and
visualize networks properties: [`igraph`](https://igraph.org/r/),
[`tidygraph`](https://github.com/thomasp85/tidygraph), and
[`ggraph`](https://github.com/thomasp85/ggraph).


## Installation

So far, only the development version is available and can be installed via the [remotes](https://CRAN.R-project.org/package=remotes) :package:

```r
R> remotes::install_github("ropensci/rmangal")
R> library("rmangal")
```


## How to use `rmangal`

There are [seven `search_*()` functions](https://docs.ropensci.org/rmangal/reference/index.html#section-explore-database) to explore the content of Mangal, for
instance `search_datasets()`:

```r
R> mgs <- search_datasets("lagoon")
Found 2 datasets
```

Once this first step achieved, networks found can be retrieved with the `get_collection()` function.

```r
R> mgn <- get_collection(mgs)
```

`get_collection()` returns an object `mgNetwork` if there is one network
returned, otherwise an object `mgNetworkCollection`, which is a list of
`mgNetwork` objects.


```r
R> class(mgn)
[1] "mgNetworksCollection"
R> mgn
A collection of 3 networks

* Network # from data set #
* Description: Dietary matrix of the Huizache–Caimanero lagoon
* Includes 189 edges and 26 nodes
* Current taxonomic IDs coverage for nodes of this network:
  --> ITIS: 81%, BOLD: 81%, EOL: 85%, COL: 81%, GBIF: 0%, NCBI: 85%
* Published in ref # DOI:10.1016/s0272-7714(02)00410-9

* Network # from data set #
* Description: Food web of the Brackish lagoon
* Includes 27 edges and 11 nodes
* Current taxonomic IDs coverage for nodes of this network:
  --> ITIS: 45%, BOLD: 45%, EOL: 45%, COL: 45%, GBIF: 18%, NCBI: 45%
* Published in ref # DOI:NA

* Network # from data set #
* Description: Food web of the Costal lagoon
* Includes 34 edges and 13 nodes
* Current taxonomic IDs coverage for nodes of this network:
  --> ITIS: 54%, BOLD: 54%, EOL: 54%, COL: 54%, GBIF: 15%, NCBI: 54%
* Published in ref # DOI:NA
```

[`igraph`](https://igraph.org/r/) and
[`tidygraph`](https://github.com/thomasp85/tidygraph) offer powerful features to
analyze networks and **rmangal** provides functions to convert `mgNetwork` to
`igraph` and `tbl_graph` so that the user can easily benefit from those
packages.

```r
R> ig <- as.igraph(mgn[[1]])
R> class(ig)
[1] "igraph"
R> library(tidygraph)
R> tg <- as_tbl_graph(mgn[[1]])
R> class(tg)
[1] "tbl_graph" "igraph"
```

:book: Note that the vignette ["Get started with
rmangal"](https://docs.ropensci.org/rmangal/articles/rmangal.html) will guide
the reader through several examples and provide further details about **rmangal** features.

## How to publish ecological networks

We are working on that part. The networks publication process will be
facilitated with structured objects and tests suite to maintain data integrity
and quality.Comments and suggestions are welcome, feel free to open issues.

## `rmangal` vs `rglobi`

Those interested only in pairwise interactions among taxons may consider using
`rglobi`, an R package that provides an interface to the [GloBi
infrastructure](https://www.globalbioticinteractions.org/about.html). GloBi
provides open access to aggregated interactions from heterogeneous sources. In
contrast, Mangal gives access to the original networks and open the gate to
study ecological networks properties (i.e. connectance, degree etc.) along large
environmental gradients, which wasn't possible using the GloBi infrastructure.


## Acknowledgment

We are grateful to [Noam Ross](https://github.com/noamross) for acting as an editor during the review process. We also thank [Anna Willoughby](https://github.com/arw36) and [Thomas Lin Petersen](https://github.com/thomasp85) for reviewing the package. Their comments strongly contributed to improve the quality of rmangal.


## Code of conduct

Please note that the `rmangal` project is released with a [Contributor Code of Conduct](https://mangal.io/doc/r/CODE_OF_CONDUCT.html). By contributing to this project, you agree to abide by its terms.

## Meta

* Get citation information for `rmangal` in R doing `citation(package = 'rmangal')`
* Please [report any issues or bugs](https://github.com/ropensci/rmangal/issues).

[![rofooter](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# rmangal 2.1.0

* All examples are within the `\donttest` tag (see #100).
* `get_collection()` methods always return an object of class `mgNetworksCollection` (see #100).
* `get_network_by_id()` gains an argument `force_collection` to force the class collection (see #100).
* Vignette now precomputed (see #100).
* Tests now use `vcr` (see #100).
* Travis and Appveyor removed, use GitHub Actions (see #100).
* `avail_type()` is no longer exported.

# rmangal 2.0.2

* Fix a minor bug `search_datasets()` related to absent networks attached on a dataset (see #97 and #98);
* Update Travis CI environment test (`travis.yml`).

# rmangal 2.0.1

* Fix a minor bug in the print method for `mgNetwork` objects see #94;
* Fix broken URIs in README;
* Remove mapview from vignette (CRAN issue with missing PhantomJS).

# rmangal 2.0.0

* Revisions see https://github.com/ropensci/software-review/issues/332;
* add summary method [#87];
* `mg_to_igraph` is now `as.igraph()`;
* `search_references()` has been rewritten [#85];
* vignette now includes examples to use `tigygraph` and `ggraph`;
* pkgdown website is now deployed by Travis CI [#86];
* `geom` column has been removed from `mgSearchInteractions` objects;
* `sf` features are only used in `search_networks_sf()` and when argument `as_sf` is set to `TRUE` [#89];
* query with spatial (`sf`) objects are handle in `query_networks_sf()` that is now exported.

# rmangal 1.9.0.9000

* Version submitted to ROpenSci for review;
* Added a `NEWS.md` file to track changes to the package.


# Previous version rmangal

* See https://github.com/mangal-wg/rmangal-v1 for the first version. Note that due to changes in the RESTful API, there is no backward compatibility.## Submission of rmangal v2.1.0

Dear CRAN, 

rmangal was archived on 2020-11-02 for policy violation on internet access. 
This minor release addresses this concern: tests now use `vcr` and the vignette has been precomputed. Minor bugs have been squashed along the way. 

## CRAN comment 


### Comment

Missing Rd-tags:
     avail_type.Rd: \value
     clear_cache_rmangal.Rd: \value


### Answer

1. `avail_type` is no longer exported, instead type of interactions are detailed in the documentation of `search_interactions()`. 

2. \value tag has been added for `clear_cache_rmangal`



## Test environments

  * GitHub Actions, Ubuntu 20.04: R-release,
  * GitHub Actions, Ubuntu 20.04: R-devel,
  * GitHub Actions, macOS 11.6: R-release,
  * GitHub Actions, Microsoft Windows Server 2019 10.0.17763: R-release,
  * win-builder (R-oldrelease, R-release and R-devel),
  * local Debian 11 (Kernel: 5.14.0-2-amd64 x86_64), R-4.1.1.


## R CMD check results

0 ERRORs | 0 WARNINGs | 0 NOTES.


## Downstream dependencies

There are currently no downstream dependencies for this package.
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
# CONTRIBUTING

## We love collaboration!

## Bugs?

* Submit an issue on the [Issues page](https://github.com/ropensci/rmangal/issues)
- be sure to include R session information and a reproducible example (repex).


## Code contributions

### Broad overview of contributing workflow

* Fork this repo to your GitHub account
* Clone your version on your account down to your machine from your account, e.g,. `git clone https://github.com/<yourgithubusername>/rmangal.git`
* Make sure to track progress upstream (i.e., on our version of `rmangal` at `ropensci/rmangal`) by doing `git remote add upstream https://github.com/ropensci/rmangal.git`. Before making changes make sure to pull changes in from upstream by doing either `git fetch upstream` then merge later or `git pull upstream` to fetch and merge in one step
* Make your changes (bonus points for making changes on a new feature branch)
* Please do write a test(s) for your changes if they affect code and not just docs (see Tests below)
* Push up to your account
* Submit a pull request to home base at `ropensci/rmangal`

### Tests

To add tests, go to the folder `tests/testthat/`. Tests are generally organized as individual files for each exported function from the package (that is, listed as an export in the `NAMESPACE` file). If you are adding a new exported function, add a new test file. If you are changing an existing function, work in the tests file for that function, unless it doesn't have tests, in which case make a new test file.

The book R packages book provides [a chapter on testing in general](http://r-pkgs.had.co.nz/tests.html). Do consult that first if you aren't familiar with testing in R.

The easiest set up to run tests is from within an R session:

```r
library(devtools)
library(testthat)
# loads the package
load_all()
```

To test an individual test file

```r
test_file("tests/testthat/test-foobar.R")
```

To run all tests

```r
devtools::test()
```

If you are running tests that have `skip_on_cran()` in them, set `Sys.setenv(NOT_CRAN = "true")` prior to running tests.


### Making changes

In addition to changing the code, do make sure to update the documentation if applicable. The R packages book book has a [chapter on documentation](http://r-pkgs.had.co.nz/man.html) you should read if you aren't familiar.

After code and documentation has been changed, update documentation by running either `devtools::document()` or `roxygen2::roxygenise()`.

Make sure if you change what packages or even functions within packages are
imported, most likely add the package to Imports in the DESCRIPTION file.

Be conservative about adding new dependencies.


### Style

* Make sure code, documentation, and comments are no more than 80 characters in width;
* Use `<-` instead of `=` for assignment;
* Always use `snake_case` (and all lowercase) instead of `camelCase`.



## Thanks for contributing!
<!--
If this issue relates to usage of the package, whether a question, bug or similar,
along with your query, please paste your devtools::session_info() or sessionInfo()
into the code block below. If not, delete all this and proceed :)
-->

<details> <summary><strong>Session Info</strong></summary>

```r

```
</details>
---
title: "Get started with rmangal"
author:
  - name: Steve Vissault & Kevin Cazelles
bibliography:
  - ../inst/bib/main.bib
date: "2021-11-09"
output: rmarkdown::html_vignette
vignette: >
  %\VignetteIndexEntry{Get started with rmangal}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
---






# Context


## The Mangal project


[The Mangal project](https://mangal.io) aims at archiving published ecological networks and at easing their retrieval. To do so, Mangal:


1. uses a data specification for ecological networks [described in @poisot_mangal_2016];


2. archives ecological networks in a [PostgreSQL](https://www.postgresql.org/) database;


3. provides:
  - [a data explorer](https://mangal.io/#/) to visualize and download data available;
  - a [RESTful Application Programming Interface (API)](https://mangal-interactions.github.io/mangal-api/);
  - a client library for Julia: [Mangal.jl](https://github.com/EcoJulia/Mangal.jl);
  - a client of this API for R: the **rmangal** package described below.


Currently, 172 datasets are including in the database representing over [1300 ecological
networks](https://mangal.io/#/network). In 2016, the first paper describing the project was
published and introduced the first release of **rmangal** [@poisot_mangal_2016]. Since then, the
structure of the database has been improved (new fields have been added), several ecological
networks added and the API entirely rewritten. Consequently, [the first release of the
**rmangal**](https://github.com/mangal-interactions/rmangal-v1) is obsolete (and archived) and we introduce
**rmangal v2.0** in this vignette.


## Data structure


<br>


<div class = "row">
<div class = "col-md-6">
<img src="img/data_structure.png" title="plot of chunk unnamed-chunk-1" alt="plot of chunk unnamed-chunk-1" width="100%" />
</div>
<div class = "col-md-6">

The diagram on the left side represents the structure of the Mangal database. All *references*
included in Mangal correspond to a specific publication that includes one or several *dataset(s)*. This dataset is
basically a collection of ecological *networks* whose *nodes* and *interactions* (edges) are stored
in separate tables. Below, we briefly describe the content of each table.


**References** -- Information pertaining to a reference (scientific article, book, online website,
etc.) characterizing an original collection of ecological networks. URLs of data and publication
sources are included as well as persistent identifiers (when available) such as digital object
identifiers (DOIs). This allows the user to retrieve more details about original publications using
appropriate R packages such as [crossref](https://docs.ropensci.org/rcrossref/).


**Datasets** -- Metadata of the datasets attached to a reference. It includes a general description
of the networks.


**Networks** -- Metadata of the networks attached to a dataset. It provides the sampling location, date and
specific description of the network.


**Nodes** -- Information on the population, taxa or individu in the network. Each node has the
original taxon name documented and taxonomic backbone provided by all services embedded in taxize
[@chamberlain_2019].

**Interactions** -- Information on the interaction type (e.g. mutualism,
predation, etc.), the strength, and the direction of the interaction between
two nodes.


</div>
</div>



## Authentication


So far, the `rmangal` package provides methods to get access to the data store. Data requests
(performed via `httr::GET()`) do not require any authentication.


A bearer authentication strategy using [ORCID](https://orcid.org/) credentials
(as a third-party services) has been implemented on all `POST`, `DELETE`, `PUT`
API operations to allow the user to add and delete new ecological to the data
base. These features are not currently included in the **rmangal** package, but
are under consideration for future major releases.




# How to use **rmangal**


## Overall approach


In order to efficiently retrieve networks from the database, **rmangal**
includes 7 search functions querying the 5 tables described above as well as a table dedicated to the taxonomy backbone.


1. `search_references()`: search in the reference table, for instance the user can look for a specific `doi`;
2. `search_datasets()`: search among datasets using a keyword;
3. `search_networks()` and `search_networks_sf()`: search networks based on a keyword or a geographical area;
4. `search_interactions()`: list all networks containing a specific interaction type;
5. `search_nodes()`: identify nodes based on nodes information;
6. `search_taxonomy()`: identify nodes based on taxonomic names and unique identifiers.



All of these functions return specific class objects with the information needed
to retrieve the corresponding set of ecological networks with
`get_collection()`. Hence, the user can easily retrieve data in two steps:


```r
networks <- search_*() %>% get_collection()
```


Note that if there is only one network to be retrieved, `get_collection()`
returns a `mgNetwork` object, otherwise it returns an object of class
`mgNetworksCollection` which is a collection (a list) of `mgNetwork` objects.
Below, we exemplify how to use the search functions, how to get a collection of
networks and how to use other packages to carry out specific analyses.

## Search functions


In **rmangal**, every functions queries a specific table and allow only one
query at a time (see section [Batch analysis](#batch-analysis) to learn
how to perform more than one query). All the functions offer two ways to query
the corresponding table:

1. a keyword: in this case, the entries returned are the partial or full keyword match of any strings contained in the table;
2. a custom query: in this case, entries returned are exact matches.


Let's load **rmangal** as well as two helper packages:


```r
library(rmangal)
library(magrittr) # for the pipe %>%
library(tibble) # to use tibble (enhanced data frames)
```

### Search and list available datasets


Let's assume we are looking for ecological networks including species living in
lagoons. If we have no idea about any existing data set, the best starting point
is then to query the `dataset` table with `lagoon` as a keyword:



```r
lagoon <- search_datasets(query = "lagoon")
class(lagoon)
#> [1] "tbl_df"           "tbl"              "data.frame"       "mgSearchDatasets"
lagoon
#> # A tibble: 2 × 10
#>      id name        description                   public created_at        updated_at        ref_id user_id references  networks  
#>   <int> <chr>       <chr>                         <lgl>  <chr>             <chr>              <int>   <int> <list>      <list>    
#> 1    22 zetina_2003 Dietary matrix of the Huizac… TRUE   2019-02-23T17:04… 2019-02-23T17:04…     22       3 <df [1 × 1… <df [1 × …
#> 2    52 yanez_1978  Food web of the Guerrero lag… TRUE   2019-02-24T23:42… 2019-02-24T23:42…     53       3 <df [1 × 1… <df [2 × …
```


If the Mangal reference id containing the lagoon networks was known, we could build a custom query as follow:



```r
lagoon_zetina <- search_datasets(list(ref_id = 22))
lagoon_zetina
#> # A tibble: 1 × 10
#>      id name        description                   public created_at        updated_at        ref_id user_id references  networks  
#>   <int> <chr>       <chr>                         <lgl>  <chr>             <chr>              <int>   <int> <list>      <list>    
#> 1    22 zetina_2003 Dietary matrix of the Huizac… TRUE   2019-02-23T17:04… 2019-02-23T17:04…     22       3 <df [1 × 1… <df [1 × …
```


Note that if an empty character is passed, i.e. `""`, all entries are returned. We can use this behavior to list all datasets available:



```r
all_datasets <- search_datasets("", verbose = FALSE)
glimpse(all_datasets)
#> Rows: 175
#> Columns: 10
#> $ id          <int> 2, 7, 9, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38,…
#> $ name        <chr> "howking_1968", "lundgren_olesen_2005", "elberling_olesen_1999", "johnston_1956", "havens_1992", "kemp_1977"…
#> $ description <chr> "Insect activity recorded on flower at Lake Hazen, Ellesmere Island, N.W.T., Canada", "Pollnator activity re…
#> $ public      <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, …
#> $ created_at  <chr> "2019-02-22T15:39:00.427Z", "2019-02-22T20:04:25.322Z", "2019-02-22T20:09:17.994Z", "2019-02-22T21:10:45.269…
#> $ updated_at  <chr> "2019-02-22T15:39:00.427Z", "2019-02-22T20:04:25.322Z", "2019-02-22T20:09:17.994Z", "2019-02-22T21:10:45.269…
#> $ ref_id      <int> 2, 7, 9, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 38, 39, 40,…
#> $ user_id     <int> 2, 2, 2, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, 3, …
#> $ references  <list> [<data.frame[1 x 11]>], [<data.frame[1 x 11]>], [<data.frame[1 x 11]>], [<data.frame[1 x 11]>], [<data.fram…
#> $ networks    <list> [<data.frame[1 x 13]>], [<data.frame[1 x 13]>], [<data.frame[1 x 13]>], [<data.frame[1 x 13]>], [<data.fram…
```


As shown in the diagram above, a dataset comes from a specific reference and `search_references()`
queries the reference table directly. A handy argument of this function is `doi` as it
allows to pass a Digital Object Identifier and so to retrieve all datasets attached to a specific
publication.




```r
zetina_2003 <- search_references(doi = "10.1016/s0272-7714(02)00410-9")
```


### Finding a specific network


We can also search by keyword across all networks.


```r
insect_coll <- search_networks(query="insect%")
glimpse(insect_coll)
#> Rows: 14
#> Error: Input must be a vector, not a <data.frame/mgSearchNetworks> object.
```

It is also possible to retrieve all networks based on interaction types involved:



```r
# List all interaction types available
avail_type()
#>  [1] "competition"  "amensalism"   "neutralism"   "commensalism" "mutualism"    "parasitism"   "predation"    "herbivory"   
#>  [9] "symbiosis"    "scavenger"    "detritivore"  "unspecified"
comp_interac <- search_interactions(type="competition")
# Number of competition interactions in mangal
nrow(comp_interac)
#> [1] 12
```

`search_networks_sf()` handles spatial queries: argument `query_sf` takes a
[`sf`](https://cran.r-project.org/package=sf) object as input and returns all
networks included in the spatial extent of this object. For instance, one can
retrieve all Californian networks included in Mangal like so:



```r
library(sf)
library(mapview)
library(USAboundaries)

area <- us_states(state = "california")
in_CA <- search_networks_sf(area, verbose = FALSE)
```


```r
mapView(st_geometry(area), color = "red", legend = FALSE, col.regions = "#FF000033") + mapView(in_CA, legend = FALSE) 
#> Error in path.expand(path): invalid 'path' argument
```

<img src="img/map1.png" title="plot of chunk unnamed-chunk-9" alt="plot of chunk unnamed-chunk-9" width="60%" />

### Search for a specific taxon


The user can easily identify networks including a specific taxonomic entity
with `search_taxonomy()`:


```r
sr_ficus <- search_taxonomy("Ficus")
```


This function allows to search for a specific taxonomic entity using it's validated
name or unique identifiers, i.e. EOL, TSN, GBIF, COL, BOLD and NCBI IDs.
Taxon names of the `taxonomy` table were validated with
TNRS (see <https://tnrs.biendata.org/> and/or GNR (see <https://resolver.globalnames.org/>). The taxon names in this table
might not be the taxon name documented in the original publication.
In order to identify relevant networks with the original name, use
[search_nodes()].

The validation of taxon names was performed by an automated
procedure using taxize [@chamberlain_2019] and if there is any doubt, the original names recorded
by authors should be regarded as the most reliable information. Please
report any issue related to taxonomy at <https://github.com/mangal-interactions/contribute/issues/new/choose>.



```r
glimpse(search_taxonomy(tsn = 28749))
#> Rows: 1
#> Error: Input must be a vector, not a <data.frame/mgSearchTaxonomy> object.
glimpse(search_taxonomy(eol = 583069))
#> Rows: 1
#> Error: Input must be a vector, not a <data.frame/mgSearchTaxonomy> object.
```


Note that in some case, one may need to find a dataset based on the original name included in the
publication, in such case, `search_nodes()` must be used:



```r
sr_ficus2 <- search_nodes("Ficus")
```


## Get networks associated with a `search_*` object


Once the search performed, ecological networks are accessible from the object
returned with `get_collection()`:



```r
nets_lagoons <- lagoon %>% get_collection
nets_in_CA <- in_CA %>% get_collection
nets_competition <- comp_interac %>% get_collection
```



```r
nets_lagoons
#> A collection of 3 networks
#> 
#> * Network #86 included in dataset #22
#> * Description: Dietary matrix of the Huizache–Caimanero lagoon
#> * Includes 189 edges and 26 nodes 
#> * Current taxonomic IDs coverage for nodes of this network: 
#>   --> ITIS: 81%, BOLD: 81%, EOL: 85%, COL: 81%, GBIF: 0%, NCBI: 85%
#> * Published in ref # DOI:10.1016/s0272-7714(02)00410-9
#> 
#> * Network #927 included in dataset #52
#> * Description: Food web of the Brackish lagoon
#> * Includes 27 edges and 11 nodes 
#> * Current taxonomic IDs coverage for nodes of this network: 
#>   --> ITIS: 45%, BOLD: 45%, EOL: 45%, COL: 45%, GBIF: 18%, NCBI: 45%
#> * Published in ref # DOI:NA
#> 
#> * Network #926 included in dataset #52
#> * Description: Food web of the Costal lagoon
#> * Includes 34 edges and 13 nodes 
#> * Current taxonomic IDs coverage for nodes of this network: 
#>   --> ITIS: 54%, BOLD: 54%, EOL: 54%, COL: 54%, GBIF: 15%, NCBI: 54%
#> * Published in ref # DOI:NA
class(nets_lagoons)
#> [1] "mgNetworksCollection"
```


Note that `mgNetworksCollection` objects are lists of `mgNetwork` object which are a list of five datasets reflecting the 5 tables presented in the diagram in the first section:



```r
names(nets_lagoons[[1]])
#> [1] "network"      "nodes"        "interactions" "dataset"      "reference"
glimpse(nets_lagoons[[1]]$network)
#> Rows: 1
#> Columns: 13
#> $ network_id       <int> 86
#> $ name             <chr> "zetina_2003_20030101_86"
#> $ date             <chr> "2003-01-01T00:00:00.000Z"
#> $ description      <chr> "Dietary matrix of the Huizache–Caimanero lagoon"
#> $ public           <lgl> TRUE
#> $ all_interactions <lgl> FALSE
#> $ created_at       <chr> "2019-02-23T17:04:34.046Z"
#> $ updated_at       <chr> "2019-02-23T17:04:34.046Z"
#> $ dataset_id       <int> 22
#> $ user_id          <int> 3
#> $ geom_type        <chr> "Point"
#> $ geom_lon         <list> -106.1099
#> $ geom_lat         <list> 22.98531
glimpse(nets_lagoons[[1]]$nodes)
#> Rows: 26
#> Columns: 19
#> $ node_id             <int> 4904, 4905, 4906, 4907, 4908, 4909, 4910, 4911, 4912, 4913, 4914, 4915, 4916, 4917, 4918, 4919, 4920…
#> $ original_name       <chr> "Scianids", "Elopids", "Lutjanids", "Carangids", "Centropomids", "Ariids", "Haemulids", "Pleuronecto…
#> $ node_level          <chr> "taxon", "taxon", "taxon", "taxon", "taxon", "taxon", "taxon", "taxon", "taxon", "taxon", "taxon", "…
#> $ network_id          <int> 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, …
#> $ taxonomy_id         <int> 4363, 4364, 4365, 4366, 4367, 4368, 4369, 4355, 3823, 4370, 4371, 4372, 4373, 4374, 4279, 4375, 4376…
#> $ created_at          <chr> "2019-02-23T17:04:42.505Z", "2019-02-23T17:04:42.571Z", "2019-02-23T17:04:42.622Z", "2019-02-23T17:0…
#> $ updated_at          <chr> "2019-02-23T17:04:42.505Z", "2019-02-23T17:04:42.571Z", "2019-02-23T17:04:42.622Z", "2019-02-23T17:0…
#> $ taxonomy.id         <int> 4363, 4364, 4365, 4366, 4367, 4368, 4369, 4355, 3823, 4370, 4371, 4372, 4373, 4374, 4279, 4375, 4376…
#> $ taxonomy.name       <chr> "Sciaenidae", "Elops", "Lutjanidae", "Carangidae", "Centropomidae", "Ariidae", "Haemulidae", "Pleuro…
#> $ taxonomy.ncbi       <int> 30870, 7927, 30850, 8157, 8184, 31017, 30840, 8256, 6762, 94935, 55118, 274463, 8079, 8219, 8189, 66…
#> $ taxonomy.tsn        <int> 169237, 28630, 168845, 168584, 167642, 43998, 169055, 172859, 13951, 165546, NA, 169013, 165876, 171…
#> $ taxonomy.eol        <int> 5211, 46561210, 5294, 5361, 5355, 5115, 5317, 5173, 46508442, 8246, 6893, 5321, 5517, 46575119, 5287…
#> $ taxonomy.bold       <int> 1856, 4061, 1858, 1851, 586, 1313, 1855, 1126, 4985, 1326, 1252, 422, 1259, NA, 1863, 1504, 28100, 7…
#> $ taxonomy.gbif       <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
#> $ taxonomy.col        <chr> "81a86c329909d507edb5c296906ef3f4", "94532a14786adeb25bcec244a53aadc1", "7150078b7dd31a5f7575240f1b7…
#> $ taxonomy.rank       <chr> "family", "genus", "family", "family", "family", "family", "family", "family", "genus", "family", "f…
#> $ taxonomy.created_at <chr> "2019-02-23T17:04:35.620Z", "2019-02-23T17:04:35.744Z", "2019-02-23T17:04:35.870Z", "2019-02-23T17:0…
#> $ taxonomy.updated_at <chr> "2019-06-14T15:25:46.438Z", "2019-06-14T15:25:46.492Z", "2019-06-14T15:25:46.546Z", "2019-06-14T15:2…
#> $ taxonomy            <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, …
glimpse(nets_lagoons[[1]]$interactions)
#> Rows: 189
#> Columns: 20
#> $ interaction_id        <int> 48376, 48377, 48378, 48379, 48380, 48381, 48382, 48383, 48384, 48385, 48388, 48389, 48390, 48391, …
#> $ node_from             <int> 4912, 4912, 4912, 4912, 4912, 4912, 4912, 4912, 4912, 4912, 4913, 4913, 4913, 4913, 4913, 4913, 49…
#> $ node_to               <int> 4912, 4914, 4915, 4918, 4919, 4920, 4921, 4922, 4925, 4926, 4914, 4916, 4917, 4919, 4920, 4922, 49…
#> $ date                  <chr> "2003-01-01T00:00:00.000Z", "2003-01-01T00:00:00.000Z", "2003-01-01T00:00:00.000Z", "2003-01-01T00…
#> $ direction             <chr> "directed", "directed", "directed", "directed", "directed", "directed", "directed", "directed", "d…
#> $ type                  <chr> "predation", "predation", "predation", "predation", "predation", "predation", "predation", "predat…
#> $ method                <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
#> $ attr_id               <int> 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12…
#> $ value                 <dbl> 0.026, 0.025, 0.003, 0.009, 0.009, 0.016, 0.284, 0.231, 0.079, 0.090, 0.100, 0.002, 0.004, 0.006, …
#> $ geom                  <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
#> $ public                <lgl> TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TRUE, TR…
#> $ network_id            <int> 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86, 86…
#> $ created_at            <chr> "2019-02-23T17:05:45.061Z", "2019-02-23T17:05:45.131Z", "2019-02-23T17:05:45.193Z", "2019-02-23T17…
#> $ updated_at            <chr> "2019-02-23T17:05:45.061Z", "2019-02-23T17:05:45.131Z", "2019-02-23T17:05:45.193Z", "2019-02-23T17…
#> $ attribute.id          <int> 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12, 12…
#> $ attribute.name        <chr> "dietary matrix", "dietary matrix", "dietary matrix", "dietary matrix", "dietary matrix", "dietary…
#> $ attribute.description <chr> "Proportions of the consumer diets made up by the prey.", "Proportions of the consumer diets made …
#> $ attribute.unit        <lgl> NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA…
#> $ attribute.created_at  <chr> "2019-02-23T17:04:25.350Z", "2019-02-23T17:04:25.350Z", "2019-02-23T17:04:25.350Z", "2019-02-23T17…
#> $ attribute.updated_at  <chr> "2019-02-23T17:04:25.350Z", "2019-02-23T17:04:25.350Z", "2019-02-23T17:04:25.350Z", "2019-02-23T17…
glimpse(nets_lagoons[[1]]$dataset)
#> Rows: 1
#> Columns: 8
#> $ dataset_id  <int> 22
#> $ name        <chr> "zetina_2003"
#> $ description <chr> "Dietary matrix of the Huizache–Caimanero lagoon"
#> $ public      <lgl> TRUE
#> $ created_at  <chr> "2019-02-23T17:04:32.017Z"
#> $ updated_at  <chr> "2019-02-23T17:04:32.017Z"
#> $ ref_id      <int> 22
#> $ user_id     <int> 3
glimpse(nets_lagoons[[1]]$reference)
#> Rows: 1
#> Columns: 11
#> $ ref_id       <int> 22
#> $ doi          <chr> "10.1016/s0272-7714(02)00410-9"
#> $ first_author <chr> "manuel j. zetina-rejon"
#> $ year         <chr> "2003"
#> $ jstor        <lgl> NA
#> $ pmid         <lgl> NA
#> $ bibtex       <chr> "@article{Zetina_Rej_n_2003, doi = {10.1016/s0272-7714(02)00410-9}, url = {https://doi.org/10.1016%2Fs0272-…
#> $ paper_url    <chr> "https://doi.org/10.1016%2Fs0272-7714%2802%2900410-9"
#> $ data_url     <chr> "https://globalwebdb.com/"
#> $ created_at   <chr> "2019-02-23T17:04:28.307Z"
#> $ updated_at   <chr> "2019-02-23T17:04:28.307Z"
```




# Integrated workflow with **rmangal**

## Batch analysis

So far, the search functions of **rmangal** allow the user to perform only a
single search at a time. The simplest way to do more than one search is to loop
over a vector or a list of queries. Below we exemplify how to do so using
`lapply()`:


```r
tsn <- c(837855, 169237)
mgn <- lapply(tsn, function(x) search_taxonomy(tsn = x)) %>%
  lapply(get_collection) %>%
  combine_mgNetworks
mgn
#> A collection of 3 networks
#> 
#> * Network #948 included in dataset #66
#> * Description: Flower and anthophilous insect interactions in the primary cool-temperate subalpine forests and meadows at Mt. Kushigata, Yamanashi Prefecture, Japan
#> * Includes 871 edges and 456 nodes 
#> * Current taxonomic IDs coverage for nodes of this network: 
#>   --> ITIS: 20%, BOLD: 33%, EOL: 46%, COL: 43%, GBIF: 35%, NCBI: 38%
#> * Published in ref # DOI:NA
#> 
#> * Network #86 included in dataset #22
#> * Description: Dietary matrix of the Huizache–Caimanero lagoon
#> * Includes 189 edges and 26 nodes 
#> * Current taxonomic IDs coverage for nodes of this network: 
#>   --> ITIS: 81%, BOLD: 81%, EOL: 85%, COL: 81%, GBIF: 0%, NCBI: 85%
#> * Published in ref # DOI:10.1016/s0272-7714(02)00410-9
#> 
#> * Network #1101 included in dataset #77
#> * Description: Food web of the Angolan fishery landings
#> * Includes 127 edges and 28 nodes 
#> * Current taxonomic IDs coverage for nodes of this network: 
#>   --> ITIS: 61%, BOLD: 50%, EOL: 61%, COL: 54%, GBIF: 4%, NCBI: 57%
#> * Published in ref # DOI:10.3989/scimar.2011.75n2309
```

## Geolocalize Mangal networks with `sf`

The function `get_collection()` has an argument `as_sf` than converts network metadata of mgNetwork objects to `sf` objects, which requires [`sf`](https://cran.r-project.org/package=sf) to be installed. This allows the user to easily
geolocalize the networks retrieved from Mangal.


```r
# assuming sf and mapview are is loaded (as we did above)
mg_lag_sf <- search_datasets(query = 'lagoon') %>% get_collection(as_sf = TRUE)
class(mg_lag_sf[[1]]$network)
#> [1] "sf"         "data.frame"
```


```r
# let's combine all these sf object into a single one
mapView(mg_lag_sf[[1]]$network) + mapView(mg_lag_sf[[2]]$network)
```
<img src="img/map2.png" title="plot of chunk unnamed-chunk-18" alt="plot of chunk unnamed-chunk-18" width="60%" />
## Taxonomic analysis with `taxize`


As Mangal includes taxonomic identifiers, **rmangal** can readily be combined
with `taxize` (see [taxize](https://github.com/ropensci/taxize) for more details about this package):


```r
library(taxize)
tsn_acer <- search_taxonomy("Acer")$taxonomy.tsn
classification(tsn_acer, db = "itis")
#> $`28749`
#>               name          rank     id
#> 1          Plantae       kingdom 202422
#> 2    Viridiplantae    subkingdom 954898
#> 3     Streptophyta  infrakingdom 846494
#> 4      Embryophyta superdivision 954900
#> 5     Tracheophyta      division 846496
#> 6  Spermatophytina   subdivision 846504
#> 7    Magnoliopsida         class  18063
#> 8          Rosanae    superorder 846548
#> 9       Sapindales         order  28643
#> 10     Sapindaceae        family  28657
#> 11            Acer         genus  28727
#> 12    Acer negundo       species  28749
#> 
#> $`28757`
#>                name          rank     id
#> 1           Plantae       kingdom 202422
#> 2     Viridiplantae    subkingdom 954898
#> 3      Streptophyta  infrakingdom 846494
#> 4       Embryophyta superdivision 954900
#> 5      Tracheophyta      division 846496
#> 6   Spermatophytina   subdivision 846504
#> 7     Magnoliopsida         class  18063
#> 8           Rosanae    superorder 846548
#> 9        Sapindales         order  28643
#> 10      Sapindaceae        family  28657
#> 11             Acer         genus  28727
#> 12 Acer saccharinum       species  28757
#> 
#> $<NA>
#> [1] NA
#> 
#> $<NA>
#> [1] NA
#> 
#> $`837855`
#>               name          rank     id
#> 1          Plantae       kingdom 202422
#> 2    Viridiplantae    subkingdom 954898
#> 3     Streptophyta  infrakingdom 846494
#> 4      Embryophyta superdivision 954900
#> 5     Tracheophyta      division 846496
#> 6  Spermatophytina   subdivision 846504
#> 7    Magnoliopsida         class  18063
#> 8          Rosanae    superorder 846548
#> 9       Sapindales         order  28643
#> 10     Sapindaceae        family  28657
#> 11            Acer         genus  28727
#> 12  Acer japonicum       species 837855
#> 
#> attr(,"class")
#> [1] "classification"
#> attr(,"db")
#> [1] "itis"
```

## Network analysis with `igraph`


Once the data are retrieved and a `mgNetwork` or a `mgNetworkCollection` objects
obtained, it is straightforward to convert it as a `igraph` (see the [dedicated
website](https://igraph.org/r/)) object and then to carry out network analysis:



```r
library(igraph)
mg_lagoons <- search_datasets(query = 'lagoon') %>% get_collection
# NB the line below returns a list of igraph objects
ig_lagoons <- as.igraph(mg_lagoons)
## Modularity analysis for the first network
modularity(ig_lagoons[[1]], membership(cluster_walktrap(ig_lagoons[[1]])))
#> [1] 0.04824893
## Degree values for all networks
lapply(ig_lagoons, degree)
#> [[1]]
#> 4904 4905 4906 4907 4908 4909 4910 4911 4912 4913 4914 4915 4916 4917 4918 4919 4920 4921 4922 4924 4925 4926 4927 4923 4929 4928 
#>   17   11   14   13   18   20   14   10   18   14   12   15    7   15   14   12   14   11   26    7   22   15   21   16    5   17 
#> 
#> [[2]]
#> 6459 6460 6461 6463 6464 6465 6458 6462 6466 6456 6457 
#>    4    7    9    3    3    7    4    6    3    4    4 
#> 
#> [[3]]
#> 6445 6447 6448 6449 6450 6452 6453 6454 6446 6451 6455 6443 6444 
#>    6    4    5    5   11    2    5    8    3    5    4    5    5
```


## Network manipulation and visualization with `tidygraph` and `ggraph`

The package [`tidygraph`](https://github.com/thomasp85/tidygraph) treats
networks as two tidy tables (one for the edges and one for the nodes) that can
be modified using the grammar of data manipulation developed in the
[tidyverse](https://www.tidyverse.org/). Moreover, `tidygraph` wraps over most
of the `igraph` functions so that the user can call a vast variety of algorithms
to properly analysis networks. Fortunately, objects of class `mgNetwork` can
readily be converted into `tbl_graph` objects which allows the user to benefit
from all the tools included in `tidygraph`:



```r
library(tidygraph)
# NB the line below would not work with a mgNetworksCollection (use lapply)
tg_lagoons <-  as_tbl_graph(mg_lagoons[[1]]) %>%
  mutate(centrality_dg = centrality_degree(mode = 'in'))
tg_lagoons %E>% as_tibble
#> # A tibble: 189 × 19
#>     from    to interaction_id date    direction type   method attr_id value public network_id created_at  updated_at  attribute.id
#>    <int> <int>          <int> <chr>   <chr>     <chr>  <lgl>    <int> <dbl> <lgl>       <int> <chr>       <chr>              <int>
#>  1     9     9          48376 2003-0… directed  preda… NA          12 0.026 TRUE           86 2019-02-23… 2019-02-23…           12
#>  2     9    11          48377 2003-0… directed  preda… NA          12 0.025 TRUE           86 2019-02-23… 2019-02-23…           12
#>  3     9    12          48378 2003-0… directed  preda… NA          12 0.003 TRUE           86 2019-02-23… 2019-02-23…           12
#>  4     9    15          48379 2003-0… directed  preda… NA          12 0.009 TRUE           86 2019-02-23… 2019-02-23…           12
#>  5     9    16          48380 2003-0… directed  preda… NA          12 0.009 TRUE           86 2019-02-23… 2019-02-23…           12
#>  6     9    17          48381 2003-0… directed  preda… NA          12 0.016 TRUE           86 2019-02-23… 2019-02-23…           12
#>  7     9    18          48382 2003-0… directed  preda… NA          12 0.284 TRUE           86 2019-02-23… 2019-02-23…           12
#>  8     9    19          48383 2003-0… directed  preda… NA          12 0.231 TRUE           86 2019-02-23… 2019-02-23…           12
#>  9     9    21          48384 2003-0… directed  preda… NA          12 0.079 TRUE           86 2019-02-23… 2019-02-23…           12
#> 10     9    22          48385 2003-0… directed  preda… NA          12 0.09  TRUE           86 2019-02-23… 2019-02-23…           12
#> # … with 179 more rows, and 5 more variables: attribute.name <chr>, attribute.description <chr>, attribute.unit <lgl>,
#> #   attribute.created_at <chr>, attribute.updated_at <chr>
tg_lagoons %N>% as_tibble %>%
  select(original_name, taxonomy.tsn, centrality_dg)
#> # A tibble: 26 × 3
#>    original_name  taxonomy.tsn centrality_dg
#>    <chr>                 <int>         <dbl>
#>  1 Scianids             169237             1
#>  2 Elopids               28630             0
#>  3 Lutjanids            168845             1
#>  4 Carangids            168584             2
#>  5 Centropomids         167642             2
#>  6 Ariids                43998             1
#>  7 Haemulids            169055             4
#>  8 Pleuronectoids       172859             3
#>  9 Callinectes           13951             6
#> 10 Belonoids            165546             4
#> # … with 16 more rows
```

Another strong advantage of `tbl_graph` objects is that there are the objects
used by the package [`ggraph`](https://github.com/thomasp85/ggraph) that that
offers various functions (theme, geoms, etc.) to efficiently visualize networks:


```r
library(ggraph)
ggraph(tg_lagoons, layout = "stress") +
  geom_edge_parallel(end_cap = circle(.5), start_cap = circle(.5),
        arrow = arrow(length = unit(1, 'mm'), type = 'closed')) +
  geom_node_point(aes(colour = taxonomy.rank), size = 8) +
  theme_graph(background = "grey40", foreground = NA, text_colour = 'white')
```

![plot of chunk ggraph](img/ggraph-1.png)




## Creating a list references for a set of networks


We can easily print the BibTeX of all publications involved in the networks collection.


```r
   search_datasets(query = 'lagoon') %>%
   get_collection %>% get_citation %>% cat(sep = "\n\n")
#> @article{Zetina_Rej_n_2003, doi = {10.1016/s0272-7714(02)00410-9}, url = {https://doi.org/10.1016%2Fs0272-7714%2802%2900410-9}, year = 2003, month = {aug}, publisher = {Elsevier {BV}}, volume = {57}, number = {5-6}, pages = {803--815}, author = {Manuel J. Zetina-Rejón and Francisco Arreguí-Sánchez and Ernesto A. Chávez}, title = {Trophic structure and flows of energy in the Huizache{	extendash}Caimanero lagoon complex on the Pacific coast of Mexico},journal = {Estuarine, Coastal and Shelf Science}}
#> 
#> @book{yanez_1978, Author = {Yáñez-Arancibia, Alejandro}, Editor = {Universidad Nacional Autónoma de México, Centro de Ciencias del Mar y Limnología. Ciudad Universitaria, México, D.F. -- 1a ed.},Title = {Taxonomía, ecología y estructura de las comunidades de peces en lagunas costeras con bocas efímeras del Pacífico de México}, Year = {1978}}
```




## References
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_nodes.R
\name{search_nodes}
\alias{search_nodes}
\title{Query nodes}
\usage{
search_nodes(query, verbose = TRUE, ...)
}
\arguments{
\item{query}{either a character string including a single keyword or a named list containing a custom query (see details section below).
Note that if an empty character string is passed, then all datasets available are returned.}

\item{verbose}{a \code{logical}. Should extra information be reported on progress?}

\item{...}{further arguments to be passed to \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
An object of class \code{mgSearchNodes}, which basically is a \code{data.frame} object
including taxa that are matching the query and corresponding information.
All networks in which taxa are involved are also attached to the \code{data.frame}.
}
\description{
Search for networks by querying the nodes table.
If the \code{query} is a character string, then all character columns in the table
are searched and the entries for which at least one
partial match was found are returned.
Alternatively, a named list can be used to look for an exact match in a specific column (see Details section)
}
\details{
Names of the list should match one of the column names within the table.
For the \code{networks} table, those are:
\itemize{
\item id: unique identifier of the nodes;
\item original_name: taxonomic name as in the original publication;
\item node_level: either population, taxon or individual;
\item network_id: Mangal network identifier.
}

Note that for lists with more than one element, only the first element is used, the others are ignored. An example is provided below.
}
\examples{
\donttest{
 res_acer <- search_nodes("Acer")
 res_926 <- search_nodes(list(network_id = 926))
}
}
\references{
\itemize{
\item \url{https://mangal.io/#/}
\item \url{https://mangal-interactions.github.io/mangal-api/#nodes}
}
}
\seealso{
\code{\link[=search_taxonomy]{search_taxonomy()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_network_by_id.R
\name{get_network_by_id}
\alias{get_network_by_id}
\alias{get_network_by_id_indiv}
\alias{print.mgNetwork}
\alias{print.mgNetworksCollection}
\alias{summary.mgNetwork}
\alias{summary.mgNetworksCollection}
\title{Retrieve network information, nodes, edges and references for a given set of Mangal network IDs}
\usage{
get_network_by_id(ids, as_sf = FALSE, force_collection = FALSE, verbose = TRUE)

get_network_by_id_indiv(id, as_sf = FALSE, verbose = TRUE)

\method{print}{mgNetwork}(x, ...)

\method{print}{mgNetworksCollection}(x, ...)

\method{summary}{mgNetwork}(object, ...)

\method{summary}{mgNetworksCollection}(object, ...)
}
\arguments{
\item{ids}{a vector of Mangal ID for networks (\code{numeric}).}

\item{as_sf}{a logical. Should networks metadata be converted into an sf object? Note that to use this feature \code{sf} must be installed.}

\item{force_collection}{a logical. Should the output to be of class  \code{mgNetworksCollection} even if it includes only one network.}

\item{verbose}{a logical. Should extra information be reported on progress?}

\item{id}{a single ID network (\code{numeric}).}

\item{x}{an object of class \code{mgNetwork} or \code{mgNetworksCollection}.}

\item{...}{ignored.}

\item{object}{object of of class \code{mgNetwork} or \code{mgNetworksCollection}.}
}
\value{
A \code{mgNetwork} object includes five data frames:
\itemize{
\item network: includes all generic information on the network (if \code{as_sf=TRUE} then it is an object of class \code{sf});
\item nodes: information pertaining to nodes (includes taxonomic information);
\item interactions: includes ecological interactions and their attributes;
\item dataset: information pertaining to the original dataset;
\item reference: details about the original publication.
}

A summary method is available for objects of class \code{mgNetwork} object and returns the following network properties:
\itemize{
\item the number of nodes;
\item the number of edges;
\item the connectance;
\item the linkage density;
\item the degree (in, out an total) and the eigenvector centrality of every nodes.
}
}
\description{
Summarize mgNetwork properties.

Summarize mgNetworksCollection properties.
}
\section{Functions}{
\itemize{
\item \code{get_network_by_id_indiv}: Retrieve a network by its  collection of networks (default).
}}

\examples{
\donttest{
 net18 <- get_network_by_id(id = 18)
 net18_c <- get_network_by_id(id = 18, force_collection = TRUE)  
 nets <- get_network_by_id(id = c(18, 23))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_networks.R
\name{search_networks}
\alias{search_networks}
\alias{search_networks_sf}
\title{Query networks}
\usage{
search_networks(query, verbose = TRUE, ...)

search_networks_sf(query_sf, verbose = TRUE, ...)
}
\arguments{
\item{query}{either a character string including a single keyword or a named list containing a custom query (see details section below), or a spatial object (see the description of \code{query_sf}).
Note that if an empty character string is passed, then all datasets available are returned.}

\item{verbose}{a \code{logical}. Should extra information be reported on progress?}

\item{...}{further arguments to be passed to \code{\link[httr:GET]{httr::GET()}}.}

\item{query_sf}{a spatial object of class \code{sf} used to search in a specific geographical area.}
}
\value{
An object of class \code{mgSearchNetworks}, which is a \code{data.frame} object with all networks informations
}
\description{
Search over all networks using a keyword, a custom query or a spatial object
If the \code{query} is a character string, then all character columns in the table
are searched and the entries for which at least one
partial match was found are returned.
Alternatively, a named list can be used to look for an exact match in a specific column (see Details section)
}
\details{
Names of the list should match one of the column names within the table.
For the \code{networks} table, those are
\itemize{
\item id: unique identifier of the network;
\item all_interactions: false interaction can be considered as real false interaction
\item dataset_id: the identifier of the dataset;
\item public: network publicly available;
}

Note that for lists with more than one element, only the first element is used, the others are ignored. An example is provided below.
}
\section{Functions}{
\itemize{
\item \code{search_networks_sf}: Search networks within a spatial object passed as an argument. Note that \code{sf} must be installed to use this function.
}}

\examples{
\donttest{
 mg_insect <- search_networks(query="insect\%")
 # Retrieve the search results
 nets_insect <- get_collection(mg_insect)
 # Spatial query
 library(sf)
 library(USAboundaries)
 area <- us_states(state="california")
 networks_in_area <- search_networks_sf(area, verbose = FALSE)
 plot(networks_in_area)
 # Retrieve network ID 5013
 net_5013 <- search_networks(query = list(id = 5013))
 # Network(s) of dataset ID 19
 mg_19 <- search_networks(list(dataset_id = 19))
}

}
\references{
\itemize{
\item \url{https://mangal.io/#/}
\item \url{https://mangal-interactions.github.io/mangal-api/#networks}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{get_gen}
\alias{get_gen}
\title{Generic API function to retrieve several entries}
\usage{
get_gen(endpoint, query = NULL, limit = 100, verbose = TRUE, ...)
}
\arguments{
\item{endpoint}{\code{character} API entry point}

\item{query}{\code{list} list of parameters passed to the API}

\item{limit}{\code{integer} number of entries return by the API (max: 1000)}

\item{verbose}{\code{logical} print API code status on error; default: \code{TRUE}}

\item{...}{Further named parameters passed to \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
Object of class \code{mgGetResponses}
}
\description{
Generic API function to retrieve several entries
}
\details{
See endpoints available with \code{endpoints()}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{get_from_fkey}
\alias{get_from_fkey}
\title{Get entries based on foreign key}
\usage{
get_from_fkey(endpoint, verbose = TRUE, ...)
}
\arguments{
\item{endpoint}{\code{character} API entry point}

\item{...}{foreign key column name with the id}
}
\value{
Object returned by \code{\link[=get_gen]{get_gen()}}
}
\description{
Get entries based on foreign key
}
\details{
See endpoints available with \code{endpoints()}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.igraph.R
\name{as.igraph.mgNetwork}
\alias{as.igraph.mgNetwork}
\alias{as.igraph.mgNetworksCollection}
\title{Coerce \code{mgNetworksCollection} or \code{mgNetwork} objects to \code{igraph} objects.}
\usage{
\method{as.igraph}{mgNetwork}(x, ...)

\method{as.igraph}{mgNetworksCollection}(x, ...)
}
\arguments{
\item{x}{either a \code{mgNetworksCollection} or a \code{mgNetwork} object.}

\item{...}{currently ignored.}
}
\value{
An object of class \code{igraph} for a \code{mgNetwork} object and a list of
\code{igraph} objects for \code{mgNetworksCollection}.
}
\description{
Coerce \code{mgNetworksCollection} or \code{mgNetwork} objects to \code{igraph} objects.
}
\section{Methods (by class)}{
\itemize{
\item \code{mgNetwork}: Convert \code{mgNetwork} objects to \code{igraph} objects.

\item \code{mgNetworksCollection}: Convert \code{mgNetworksCollection} objects to a list of \code{igraph} objects.
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/cache.R
\name{clear_cache_rmangal}
\alias{clear_cache_rmangal}
\title{Clear memoise cache}
\usage{
clear_cache_rmangal()
}
\value{
\code{TRUE} when the cache has been reset.
}
\description{
Resets the cache of the memoised function used for http GET queries (see \code{\link[memoise:forget]{memoise::forget()}}).
}
\examples{
clear_cache_rmangal()
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\name{get_singletons}
\alias{get_singletons}
\title{Generic API function to retrieve singletons}
\usage{
get_singletons(endpoint = NULL, ids = NULL, verbose = TRUE, ...)
}
\arguments{
\item{endpoint}{\code{character} API entry point.}

\item{ids}{\code{numeric} vector of ids.}

\item{verbose}{\code{logical} print API code status on error; default: \code{TRUE}}

\item{...}{httr options, see \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
Object of class \code{mgGetResponses}
}
\description{
Generic API function to retrieve singletons
}
\details{
See endpoints available with \code{endpoints()}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_datasets.R
\name{search_datasets}
\alias{search_datasets}
\title{Query datasets}
\usage{
search_datasets(query, verbose = TRUE, ...)
}
\arguments{
\item{query}{either a character string including a single keyword or a named list containing a custom query (see details section below).
Note that if an empty character string is passed, then all datasets available are returned.}

\item{verbose}{a logical. Should extra information be reported on progress?}

\item{...}{further arguments to be passed to \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
An object of class \code{mgSearchDatasets}, which basically is a \code{data.frame}
object including all datasets corresponding to the query. For each dataset
entry,  the networks and the original reference are attached.
}
\description{
Identify relevant datasets using a keyword or a custom query.
If the \code{query} is a character string, then all character columns in the table
are searched and the entries for which at least one
partial match was found are returned.
Alternatively, a named list can be used to look for an exact match in a specific column (see Details section)
}
\details{
If \code{query} is a named list, the name  used should be one of the following:
\itemize{
\item id: unique identifier of the dataset
\item name: name of the dataset
\item date: date (\code{YYYY-mm-dd}) of the corresponding publication
\item description: a brief description of the data set
\item ref_id: the Mangal identifier of the dataset
}

Note that for lists with more than one element, only the first element is used, the others are ignored.
Examples covering custom queries are provided below.
}
\examples{
\donttest{
 # Return all datasets (takes time)
 all_datasets <- search_datasets("")
 all_datasets
 class(all_datasets)
 # Search with keyword
 mg_lagoon <- search_datasets(query = 'lagoon')
 # Search with a custom query (specific column)
 mg_kemp <- search_datasets(query = list(name = 'kemp_1977'))
 mg_16 <- search_datasets(query = list(ref_id = 16))
}
}
\references{
\itemize{
\item \url{https://mangal.io/#/}
\item \url{https://mangal-interactions.github.io/mangal-api/#datasets}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_citations.R
\name{get_citation}
\alias{get_citation}
\alias{get_citation.mgNetwork}
\alias{get_citation.mgNetworksCollection}
\title{Retrieve all references pertaining to the networks collection or individual network}
\usage{
get_citation(x)

\method{get_citation}{mgNetwork}(x)

\method{get_citation}{mgNetworksCollection}(x)
}
\arguments{
\item{x}{an object of class \code{mgNetworksCollection} or \code{mgNetworks}.}
}
\value{
Bibtex entries as a character vector.
}
\description{
Retrieve all references pertaining to the networks collection or individual network
}
\section{Methods (by class)}{
\itemize{
\item \code{mgNetwork}: Get BibTeX entries for the publication associated to the network.

\item \code{mgNetworksCollection}: Get BibTeX entries for the publication associated to the networks.
}}

\examples{
\donttest{
 # network collection
 lagoon_net_collection <- get_collection(search_datasets("lagoon"))
 get_citation(lagoon_net_collection)
 # individual network
 mg_18 <- get_network_by_id(18)
 get_citation(mg_18)
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_interactions.R
\name{search_interactions}
\alias{search_interactions}
\title{Query interactions}
\usage{
search_interactions(
  query,
  type = NULL,
  expand_node = FALSE,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{either a character string including a single keyword or a named list containing a custom query (see details section below).
Note that if an empty character string is passed, then all datasets available are returned.}

\item{type}{a \code{character} one of the interactions type available (see details). Note that \code{query} is ignored if \code{type} is used.}

\item{expand_node}{a logical. Should the function returned extra information pertaining to nodes? Default is set to \code{FALSE}, which means that only the Mangal IDs of nodes are returned.}

\item{verbose}{a \code{logical}. Should extra information be reported on progress?}

\item{...}{further arguments to be passed to \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
An object of class \code{mgSearchInteractions}, i.e. a \code{data.frame} object including interactions.
All networks in which interactions are involved are also attached to the \code{data.frame}.
}
\description{
Search for specific interactions using a keyword or a specific type of
interactions (e.g. mutualism). If the \code{query} is a character string, then all character columns in the table
are searched and the entries for which at least one
partial match was found are returned.
Alternatively, a named list can be used to look for an exact match in a specific column (see Details section)
}
\details{
Names of the list should match one of the column names within the table.
For the \code{interaction} table, those are:
\itemize{
\item id: unique identifier of the interaction;
\item attr_id: identifier of a specific attribute;
\item direction: edge direction ("directed", "undirected" or "unknown");
\item network_id: Mangal network identifier;
\item node_from: node id which the interaction end to;
\item node_to: node to which the interaction end to;
\item type: use argument \code{type} instead.
}

Note that for lists with more than one element, only the first element is
used, the others are ignored. The type of interactions (argument \code{type})
currently available are the following
\itemize{
\item "competition";
\item "amensalism";
\item "neutralism";
\item "commensalism";
\item "mutualism";
\item "parasitism";
\item "predation";
\item "herbivory";
\item "symbiosis";
\item "scavenger";
\item "detritivore".
}
}
\examples{
\donttest{
 df_inter <- search_interactions(type = "competition", verbose = FALSE)
 # Get all networks containing competition
 competition_networks <- get_collection(df_inter, verbose = FALSE)
 df_net_926 <- search_interactions(list(network_id = 926), verbose = FALSE)
}
}
\references{
\itemize{
\item \url{https://mangal.io/#/}
\item \url{https://mangal-interactions.github.io/mangal-api/#taxonomy}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_references.R
\name{search_references}
\alias{search_references}
\title{Query references}
\usage{
search_references(query, doi = NULL, verbose = TRUE, ...)
}
\arguments{
\item{query}{either a character string including a single keyword or a named list containing a custom query (see details section below).
Note that if an empty character string is passed, then all datasets available are returned.}

\item{doi}{\code{character} a Digital Object Identifier  (DOI) of the article. Note that \code{query} is ignored if \code{doi} is specified.}

\item{verbose}{a \code{logical}. Should extra information be reported on progress?}

\item{...}{further arguments to be passed to \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
An object of class \code{mgSearchReferences}, which is a list that includes a
wide range of details associated to the reference, including all datasets
and networks related to the publication that are included in Mangal database.
}
\description{
Search for a specific reference using a key word or a Digital Object Identifier (DOI).
If the \code{query} is a character string, then all character columns in the table
are searched and the entries for which at least one
partial match was found are returned.
Alternatively, a named list can be used to look for an exact match in a specific column (see Details section).
}
\details{
Names of the list should match one of the column names within the table.
For the \code{reference} table, those are:
\itemize{
\item id: unique identifier of the reference
\item first_author: first author
\item doi: use \code{doi} instead
\item jstor: JSTOR identifier
\item year: year of publication.
}

Note that for lists with more than one element, only the first element is used, the others are ignored. An example is provided below.
}
\examples{
\donttest{
 search_references(doi = "10.2307/3225248")
 search_references(list(jstor = 3683041))
 search_references(list(year = 2010))
}
}
\references{
\itemize{
\item \url{https://mangal.io/#/}
\item \url{https://mangal-interactions.github.io/mangal-api/#references}
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/as.igraph.R
\docType{import}
\name{reexports}
\alias{reexports}
\alias{as.igraph}
\title{Objects exported from other packages}
\keyword{internal}
\description{
These objects are imported from other packages. Follow the links
below to see their documentation.

\describe{
  \item{igraph}{\code{\link[igraph]{as.igraph}}}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_collection.R
\name{get_collection}
\alias{get_collection}
\alias{get_collection.default}
\alias{get_collection.mgSearchDatasets}
\alias{get_collection.mgSearchNetworks}
\alias{get_collection.mgSearchReferences}
\alias{get_collection.mgSearchNodes}
\alias{get_collection.mgSearchTaxonomy}
\alias{get_collection.mgSearchInteractions}
\title{Get a collection of networks}
\usage{
get_collection(x, ...)

\method{get_collection}{default}(x, ...)

\method{get_collection}{mgSearchDatasets}(x, ...)

\method{get_collection}{mgSearchNetworks}(x, ...)

\method{get_collection}{mgSearchReferences}(x, ...)

\method{get_collection}{mgSearchNodes}(x, ...)

\method{get_collection}{mgSearchTaxonomy}(x, ...)

\method{get_collection}{mgSearchInteractions}(x, ...)
}
\arguments{
\item{x}{\code{numeric} vector of Mangal network IDs or an object returned by
by one of the \verb{search_*()} functions.}

\item{...}{arguments to be passed on to \code{\link[=get_network_by_id]{get_network_by_id()}}.}
}
\value{
Returns a object of class \code{mgNetworksCollection} which is a collection
(actually, a list) of \code{mgNetwork} objects \code{\link[=get_network_by_id]{get_network_by_id()}}).
}
\description{
Retrieve a set of networks based on the results of one of the \verb{search_*()}
function. The function also accepts a numeric vector of Mangal network IDs.
}
\section{Methods (by class)}{
\itemize{
\item \code{default}: Get a collection of networks (default).

\item \code{mgSearchDatasets}: Get a collection of networks from a \code{mgSearchDatasets} object.

\item \code{mgSearchNetworks}: Get a collection of networks from a \code{mgSearchNetworks} object.

\item \code{mgSearchReferences}: Get a collection of networks from a \code{mgSearchReferences} object.

\item \code{mgSearchNodes}: Get a collection of networks from a \code{mgSearchNodes} object.

\item \code{mgSearchTaxonomy}: Get a collection of networks from a \code{mgSearchTaxa} object.

\item \code{mgSearchInteractions}: Get a collection of networks from a \code{mgSearchTaxa} object.
}}

\examples{
\donttest{
 mg_2 <- get_collection(c(1076:1077), verbose = FALSE)
 mg_anemone <- get_collection(search_networks(query='anemone\%'), verbose = FALSE)
}
}
\seealso{
\code{\link[=search_datasets]{search_datasets()}}, \code{\link[=search_interactions]{search_interactions()}}, \code{\link[=search_networks]{search_networks()}},
\code{\link[=search_nodes]{search_nodes()}}, \code{\link[=search_references]{search_references()}}, \code{\link[=search_taxonomy]{search_taxonomy()}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/zzz.R
\docType{package}
\name{rmangal}
\alias{rmangal}
\alias{rmangal-package}
\title{rmangal}
\description{
A programmatic interface to the Mangal API \url{https://mangal-interactions.github.io/mangal-api/}.
}
\seealso{
Useful links:
\itemize{
  \item \url{https://docs.ropensci.org/rmangal/}
  \item \url{https://mangal.io}
  \item \url{https://github.com/ropensci/rmangal}
  \item Report bugs at \url{https://github.com/ropensci/rmangal/issues}
}

}
\author{
\strong{Maintainer}: Kevin Cazelles \email{kevin.cazelles@gmail.com} (\href{https://orcid.org/0000-0001-6619-9874}{ORCID})

Authors:
\itemize{
  \item Steve Vissault \email{s.vissault@yahoo.fr} (\href{https://orcid.org/0000-0002-0866-4376}{ORCID}) [contributor]
  \item Gabriel Bergeron \email{gabriel.bergeron3@usherbrooke.ca} [contributor]
  \item Benjamin Mercier \email{Benjamin.B.Mercier@usherbrooke.ca} [contributor]
  \item Clément Violet \email{Clement.Violet@etudiant.univ-brest.fr} [contributor]
  \item Dominique Gravel \email{dominique.gravel@usherbrooke.ca}
  \item Timothée Poisot \email{timothee.poisot@umontreal.ca}
}

Other contributors:
\itemize{
  \item Thomas Lin Pedersen (\href{https://orcid.org/0000-0002-5147-4711}{ORCID}) [reviewer]
  \item Anna Willoughby (\href{https://orcid.org/0000-0002-0504-0605}{ORCID}) [reviewer]
}

}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_taxonomy.R
\name{search_taxonomy}
\alias{search_taxonomy}
\title{Query taxonomy}
\usage{
search_taxonomy(
  query,
  tsn = NULL,
  gbif = NULL,
  eol = NULL,
  col = NULL,
  bold = NULL,
  ncbi = NULL,
  verbose = TRUE,
  ...
)
}
\arguments{
\item{query}{a character string including a single keyword.
Note that if an empty character string is passed, then all datasets available are returned.}

\item{tsn}{a \code{numeric}. Unique taxonomic identifier from Integrated Taxonomic Information System (\url{https://www.itis.gov}).}

\item{gbif}{a \code{numeric}. Unique taxonomic identifier from Global Biodiversity Information Facility (\url{https://www.gbif.org}).}

\item{eol}{a \code{numeric}. Unique taxonomic identifier from Encyclopedia of Life (\url{https://eol.org}).}

\item{col}{a \code{numeric}. Unique taxonomic identifier from Catalogue of Life (\url{https://www.catalogueoflife.org}).}

\item{bold}{a \code{numeric}. Unique taxonomic identifier from Barcode of Life (\url{http://www.boldsystems.org}).}

\item{ncbi}{a \code{numeric}. Unique taxonomic identifier from National Center for Biotechnology Information (\url{https://www.ncbi.nlm.nih.gov}).}

\item{verbose}{a \code{logical}. Should extra information be reported on progress?}

\item{...}{further arguments to be passed to \code{\link[httr:GET]{httr::GET()}}.}
}
\value{
An object of class \code{mgSearchTaxonomy}, which is a \code{data.frame} including
all taxa matching the query.
}
\description{
Search network by taxon names and unique taxonomic identifiers.
This function offers the opportunity to retrieve taxon based on (i) known identifier
such as the taxonomic serial number (TSN), GBIF ID etc. or (ii) text search using partial match.
Have a look at the list of arguments to see the complete list of identifiers accessible.
If any unique identifier argument is used (i.e. tsn etc.), then \code{query} is ignored. Moreover,
if several taxonomic identifiers are specified, then only the first one is considered.
}
\details{
Taxon names of the \code{taxonomy} table were validated with
TNRS (see \url{https://tnrs.biendata.org} and/or GNR
might not be the taxon name documented in the original publication.
In order to identify relevant networks with the original name, use
\code{\link[=search_nodes]{search_nodes()}}.

The validation of taxon names was performed by an automated
procedure and if there is any doubt, the original names recorded
by authors should be regarded as the most reliable information. Please
report any issue related to taxonomy at \url{https://github.com/mangal-interactions/contribute/issues/new/choose/}.
}
\examples{
\donttest{
 search_taxonomy("Acer")
 # Retrieve higher classification
 tsn_acer <- search_taxonomy("Acer")$taxonomy.tsn
}
}
\references{
\itemize{
\item \url{https://mangal.io/#/}
\item \url{https://mangal-interactions.github.io/mangal-api/#taxonomy}
}
}
\seealso{
\code{\link[=search_nodes]{search_nodes()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/search_interactions.R
\name{avail_type}
\alias{avail_type}
\title{List interactions type contains in mangal-db}
\usage{
avail_type()
}
\description{
List interactions type contains in mangal-db
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/combine_mgNetworks.R
\name{combine_mgNetworks}
\alias{combine_mgNetworks}
\title{Combine Mangal networks}
\usage{
combine_mgNetworks(...)
}
\arguments{
\item{...}{objects of class \code{mgNetworksCollection} or \code{mgNetwork} or a list #' of objects of these classes.}
}
\value{
An object of class \code{mgNetworksCollection}.
}
\description{
Combine \code{mgNetworksCollection} and \code{mgNetwork} objects into a
\code{mgNetworksCollection} object.
}
\examples{
\donttest{
 mg_random_1071 <- get_collection(c(1071))
 mg_random_1074 <- get_collection(c(1074))
 combine_mgNetworks(mg_random_1071, mg_random_1074)
} 

}
