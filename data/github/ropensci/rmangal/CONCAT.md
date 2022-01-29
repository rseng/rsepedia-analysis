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
