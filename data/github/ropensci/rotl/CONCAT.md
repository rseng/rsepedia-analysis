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

<!-- badges: start -->

[![R build
status](https://github.com/ropensci/rotl/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rotl/actions)
[![codecov.io](https://codecov.io/github/ropensci/rotl/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rotl?branch=master)
[![](https://www.r-pkg.org/badges/version/rotl)](https://www.r-pkg.org/pkg/rotl)
[![CRAN RStudio mirror
downloads](https://cranlogs.r-pkg.org/badges/rotl)](https://www.r-pkg.org/pkg/rotl)
[![](https://badges.ropensci.org/17_status.svg)](https://github.com/ropensci/software-review/issues/17)
[![Project Status: Active – The project has reached a stable, usable
state and is being actively
developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

<!-- badges: end -->

# An R interface to Open Tree API

`rotl` is an R package to interact with the Open Tree of Life data APIs.
It was initially developed as part of the [NESCENT/OpenTree/Arbor
hackathon](https://blog.opentreeoflife.org/2014/06/11/apply-for-tree-for-all-a-hackathon-to-access-opentree-resources/).

Client libraries to interact with the Open Tree of Life API also exists
for [Python](https://github.com/OpenTreeOfLife/pyopentree) and
[Ruby](https://github.com/SpeciesFileGroup/bark).

## Installation

The current stable version is available from CRAN, and can be installed
by typing the following at the prompt in R:

``` r
install.packages("rotl")
```

If you want to test the development version, you first need to install
the `remotes` package.

``` r
install.packages("remotes")
```

Then you can install `rotl` using:

``` r
remotes::install_github("ropensci/rotl")
```

## Vignettes

There are three vignettes:

  - Start by checking out the “How to use `rotl`?” by typing:
    `vignette("rotl", package="rotl")` after installing the package.

  - Then explore how you can use `rotl` with other packages to combine
    your data with trees from the Open Tree of Life project by typing:
    `vignette("data_mashups", package="rotl")`.

  - The vignette “Using the Open Tree Synthesis in a comparative
    analsysis” demonstrates how you can reproduce an analysis of a
    published paper by downloading the tree they used, and data from the
    supplementary material: `vignette("meta-analysis", package="rotl")`.

The vignettes are also available from CRAN: [How to use
`rotl`?](https://cran.r-project.org/package=rotl/vignettes/rotl.html),
[Data
mashups](https://cran.r-project.org/package=rotl/vignettes/data_mashups.html),
and [Using the Open Tree synthesis in a comparative
analysis](https://cran.r-project.org/package=rotl/vignettes/meta-analysis.html).

## Quick start

### Get a little bit of the big Open Tree tree

Taxonomic names are represented in the Open Tree by numeric identifiers,
the `ott_ids` (Open Tree Taxonomy identifiers). To extract a portion of
a tree from the Open Tree, you first need to find `ott_ids` for a set of
names using the `tnrs_match_names` function:

``` r
library(rotl)
apes <- c("Pongo", "Pan", "Gorilla", "Hoolock", "Homo")
(resolved_names <- tnrs_match_names(apes))
```

    ##   search_string unique_name approximate_match ott_id is_synonym          flags
    ## 1         pongo       Pongo             FALSE 417949      FALSE               
    ## 2           pan         Pan             FALSE 417957      FALSE sibling_higher
    ## 3       gorilla     Gorilla             FALSE 417969      FALSE sibling_higher
    ## 4       hoolock     Hoolock             FALSE 712902      FALSE               
    ## 5          homo        Homo             FALSE 770309      FALSE sibling_higher
    ##   number_matches
    ## 1              2
    ## 2              1
    ## 3              1
    ## 4              1
    ## 5              1

Now we can get the tree with just those tips:

``` r
tr <- tol_induced_subtree(ott_ids = ott_id(resolved_names))
```

    ## Warning in collapse_singles(tr, show_progress): Dropping singleton nodes with
    ## labels: mrcaott83926ott6145147, mrcaott83926ott3607728, mrcaott83926ott3607876,
    ## mrcaott83926ott3607873, mrcaott83926ott3607687, mrcaott83926ott3607716,
    ## mrcaott83926ott3607689, mrcaott83926ott3607732, mrcaott770295ott3607719,
    ## mrcaott770295ott3607692, Ponginae ott1082538, Hylobatidae ott166544

``` r
plot(tr)
```

![](https://i.imgur.com/5Fdb927.png)<!-- -->

The code above can be summarized in a single pipe:

``` r
library(magrittr)
```

    ## 
    ## Attaching package: 'magrittr'

    ## The following objects are masked from 'package:testthat':
    ## 
    ##     equals, is_less_than, not

``` r
## or expressed as a pipe:
c("Pongo", "Pan", "Gorilla", "Hoolock", "Homo") %>%
  tnrs_match_names() %>%
  ott_id() %>%
  tol_induced_subtree() %>%
  plot()
```

    ## Warning in collapse_singles(tr, show_progress): Dropping singleton nodes with
    ## labels: mrcaott83926ott6145147, mrcaott83926ott3607728, mrcaott83926ott3607876,
    ## mrcaott83926ott3607873, mrcaott83926ott3607687, mrcaott83926ott3607716,
    ## mrcaott83926ott3607689, mrcaott83926ott3607732, mrcaott770295ott3607719,
    ## mrcaott770295ott3607692, Ponginae ott1082538, Hylobatidae ott166544

![](https://i.imgur.com/43LgNKf.png)<!-- -->

## Citation and Manuscript

To cite `rotl` in publications pleases use:

> Michonneau, F., Brown, J. W. and Winter, D. J. (2016), rotl: an R
> package to interact with the Open Tree of Life data. Methods in
> Ecology and Evolution. 7(12):1476-1481. doi:
> [10.1111/2041-210X.12593](https://doi.org/10.1111/2041-210X.12593)

You may also want to cite the paper for the Open Tree of Life

> Hinchliff, C. E., et al. (2015). Synthesis of phylogeny and taxonomy
> into a comprehensive tree of life. Proceedings of the National Academy
> of Sciences 112.41 (2015): 12764-12769 doi:
> [10.1073/pnas.1423041112](https://doi.org/10.1073/pnas.1423041112)

The manuscript in *Methods in Ecology and Evolution* includes additional
examples on how to use the package. The manuscript and the code it
contains are also hosted on GitHub at:
<https://github.com/fmichonneau/rotl-ms>

## Versioning

Starting with v3.0.0 of the package, the major and minor version numbers
(the first 2 digits of the version number) will be matched to those of
the API. The patch number (the 3rd digit of the version number) will be
used to reflect bug fixes and other changes that are independent from
changes to the API.

`rotl` can be used to access other versions of the API (if they are
available) but most likely the high level functions will not work.
Instead, you will need to parse the output yourself using the “raw”
returns from the unexported low-level functions (all prefixed with a
`.`). For instance to use the `tnrs/match_names` endpoint for `v2` of
the API:

``` r
rotl:::.tnrs_match_names(c("pan", "pango", "gorilla", "hoolock", "homo"), otl_v = "v2")
```

### Code of Conduct

Please note that this project is released with a [Contributor Code of
Conduct](https://github.com/ropensci/rotl/blob/master/CONDUCT.md). By
participating in this project you agree to abide by its terms.

[![](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org)
# rotl 3.0.12

* The default argument `context_name` for the function `tnrs_match_names` was
  changed from `NULL` to `All life`. This changes is made to address what could
  lead to unexpected results. Previously, the context was inferred based on the
  first match when several names were provided (see #134 reported by @LunaSare,
  and https://github.com/OpenTreeOfLife/feedback/issues/528)/
* The "Suggest" dependency `fulltext` was removed following its archival from
  CRAN.

# rotl 3.0.11

## Other changes

* When none of the names provided to `tnrs_match_names` had a match in the Open
  Tree of Life, an error was thrown and nothing was returned. To make the
  behavior of the function more consistent with other situations, when none of
  the names provided have a match, a tibble is returned and a warning is issued.

## Bug Fix

* When attempting to match a name that did not exist in the tree, an error was
  thrown (bug #128, PR #129, @daijiang)

# rotl 3.0.10

* Small fixes following updates to the Open Tree of Life API (no visible change
  for users).
* Updated documentation to reflect new value in output of `tol_node_info()`.

# rotl 3.0.9

* Small fixes following updates to the Open Tree of Life API.

# rotl 3.0.7

* Minor update to vignette to address change to TNRS endpoint (underscores can't
  be included in the taxon names anymore).

# rotl 3.0.6

* Minor update to address warnings seen on CRAN.

# rotl 3.0.5

## New features

* The data types in the data frame returned by `tnrs_match_names` are
  consistent, and remain the same even after using `update()`.
  
## Other changes

* Small internal changes that reflect changes in the data structures returned by
  the API.

# rotl 3.0.4

## New features

* To improve stability of results across releases of the Open Tree Taxonomy, the
  TNRS match with the lowest OTT id is returned instead of the first one in case
  a name is shared across multiple domains (related to #88)
* A warning is issued when users attempt to use TNRS on duplicated names.

## Other changes

* Fix typos and workaround broken API to retrieve supplementary materials in
  vignette.

# rotl 3.0.3

## New features

* The function `get_study_subtree` gains the argument `tip_label` to control the
  formatting of the tip labels, #90, reported by @bomeara
* The new function `is_in_tree` takes a list of OTT ids (i.e., the output of
  `ott_id()`), and returns a vector of logical indiicating whether they are
  included in the synthetic tree (workaround #31).

## Bug fixes

* The function `get_study_subtree` ignored the argument `subtree_id`, #89
  reported by @bomeara

## Other chaanges

* `citation("rotl")` now includes the reference to the Open Tree of Life
  publication.
* The "How to use rotl?" vignette was updated to document the behavior of v3 of
  the OTL API which returns an HTTP error code 400 when the request for induced
  subtree includes taxa that are not in the synthetic tree (fix #84)

# rotl 3.0.1

* Fix tests and vignette to reflect changes accompanying release 6.1 of the
  synthetic tree

* Add section in vignette "How to use rotl?" about how to get the higher
  taxonomy from a given taxon.

* Add `CITATION` file with MEE manuscript information (#82)

# rotl 3.0.0

* `rotl` now interacts with v3.0 of the Open Tree of Life APIs. The
  documentation has been updated to reflect the associated changes. More
  information about the v3.0 of the Open Tree of Life APIs can be found
  [on their wiki](https://github.com/OpenTreeOfLife/germinator/wiki/Open-Tree-of-Life-Web-APIs).


## New features

* New methods: `tax_sources`, `is_suppressed`, `tax_rank`, `unique_name`,
  `name`, `ott_id`, for objects returned by `tnrs_match_names()`,
  `taxonomy_taxon_info()`, `taxonomy_taxon_mrca()`, `tol_node_info()`,
  `tol_about()`, and `tol_mrca()`. Each of these methods have their own class.

* New method `tax_lineage()` to extract the higher taxonomy from an object
  returned by `taxonomy_taxon_info()` (initally suggested by Matt Pennell, #57).

* New method `tol_lineage()` to extract the nodes towards the root of the tree.

* New print methods for `tol_node_info()` and `tol_mrca()`.

* New functions `study_external_IDs()` and `taxon_external_IDs()` that return
  the external identifiers for a study and associated trees (e.g., DOI, TreeBase
  ID); and the identifiers of taxon names in taxonomic databases. The vignette
  "Data mashup" includes an example on how to use it.

* The function `strip_ott_id()` gains the argument `remove_underscores` to remove
  underscores from tips in trees returned by OTL.

## Changes

* Rename method `ott_taxon_name()` to `tax_name()` for consistency.

* Rename method `synth_sources()` and `study_list()` to `source_list()`.

* Refactor how result of query is checked and parsed (invisible to the user).

## Bug fixes

* Fix bug in `studies_find_studies()`, the arguments `verbose` and `exact` were
  ignored.

* The argument `only_current` has been dropped for the methods associated with
  objects returned by `tnrs_match_names()`

* The print method for `tnrs_context()` duplicated some names.

* `inspect()`, `update()` and `synonyms()` methods for `tnrs_match_names()` did
  not work if the query included unmatched taxa.


# rotl 0.5.0

* New vignette: `meta-analysis`

* Added arguments `include_lineage` and `list_terminal_descendants` to
  `taxonomy_taxon()`

* Improve warning and format of the result if one of the taxa requested doesn't
  match anything `tnrs_match_names`.

* In the data frame returned by `tnrs_match_names`, the columns
  `approximate_match`, `is_synonym` and `is_deprecated` are now `logical`
  (instead of `character`) [issue #54]

* New utility function `strip_ott_ids` removes OTT id information from
  a character vector, making it easier to match tip labels in trees returned by
  `tol_induced_subtree` to taxonomic names in other data sources. This function
  can also remove underscores from the taxon names.

* New method `list_trees` returns a list of tree ids associated with
  studies. The function takes the output of `studies_find_studies` or
  `studies_find_trees`.

* `studies_find_studies` and `studies_find_trees` gain argument `detailed`
  (default set to `TRUE`), that produces a data frame summarizing information
  (title of the study, year of publication, DOI, ids of associated trees, ...)
  about the studies matching the search criteria.

* `get_study_tree` gains argument `deduplicate`. When `TRUE`, if the tree
  returned for a given study contains duplicated tip labels, they will be made
  unique before being parsed by NCL by appending a suffix (`_1`, `_2`, `_3`,
  etc.). (#46, reported by @bomeara)

* New method `get_study_year` for objects of class `study_meta` that returns the
  year of publication of the study.

* A more robust approach is used by `get_tree_ids` to identify the tree ids in
  the metadata returned by the API

# rotl 0.4.1

* Initial CRAN release on July, 24th 2015
## Test environments

- local Ubuntu 20.10, R 4.1.1
- Ubuntu 18.04 (GitHub Actions) R 3.4.4, R 3.5.3, R 3.6.3, R 4.0.5, R 4.1.2, and R-devel (2021-11-16 r81199)
- macOS Big Sur 10.16, R 4.1.2
- Windows (GitHub Actions), R 4.1.2, R 3.6.3
- Fedora (clang) on R-hub with R-devel (2019-10-17 r79346)
- Ubuntu (gcc) on R-hub with R 4.0.3


## R CMD check results

- There were no ERRORs, WARNINGs, or NOTEs

## Downstream dependencies

* taxize: no WARNING or NOTE generated.
---
output: github_document
---

```{r setup, eval=TRUE, echo=FALSE}
library(knitr)
# opts_knit$set(upload.fun = image_uri)
opts_knit$set(upload.fun = imgur_upload)
```

<!-- badges: start -->
[![R build status](https://github.com/ropensci/rotl/workflows/R-CMD-check/badge.svg)](https://github.com/ropensci/rotl/actions)
[![codecov.io](https://codecov.io/github/ropensci/rotl/coverage.svg?branch=master)](https://codecov.io/github/ropensci/rotl?branch=master)
[![](https://www.r-pkg.org/badges/version/rotl)](https://www.r-pkg.org/pkg/rotl)
[![CRAN RStudio mirror downloads](https://cranlogs.r-pkg.org/badges/rotl)](https://www.r-pkg.org/pkg/rotl)
[![](https://badges.ropensci.org/17_status.svg)](https://github.com/ropensci/software-review/issues/17)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

<!-- badges: end -->


# An R interface to Open Tree API

`rotl` is an R package to interact with the Open Tree of Life data APIs. It was
initially developed as part of the
[NESCENT/OpenTree/Arbor hackathon](https://blog.opentreeoflife.org/2014/06/11/apply-for-tree-for-all-a-hackathon-to-access-opentree-resources/).

Client libraries to interact with the Open Tree of Life API also exists for
[Python](https://github.com/OpenTreeOfLife/pyopentree)
and [Ruby](https://github.com/SpeciesFileGroup/bark).


## Installation

The current stable version is available from CRAN, and can be installed by
typing the following at the prompt in R:

```{r, eval=FALSE}
install.packages("rotl")
```

If you want to test the development version, you first need to install
the `remotes` package.

```{r, eval=FALSE}
install.packages("remotes")
```

Then you can install `rotl` using:

```{r, eval=FALSE}
remotes::install_github("ropensci/rotl")
```

## Vignettes

There are three vignettes:

- Start by checking out the "How to use `rotl`?" by typing:
  `vignette("rotl", package="rotl")` after installing the
  package.

- Then explore how you can use `rotl` with other packages to combine your data
  with trees from the Open Tree of Life project by typing:
  `vignette("data_mashups", package="rotl")`.

- The vignette "Using the Open Tree Synthesis in a comparative analsysis"
  demonstrates how you can reproduce an analysis of a published paper by
  downloading the tree they used, and data from the supplementary material:
  `vignette("meta-analysis", package="rotl")`.

The vignettes are also available from CRAN:
[How to use `rotl`?](https://cran.r-project.org/package=rotl/vignettes/rotl.html),
[Data mashups](https://cran.r-project.org/package=rotl/vignettes/data_mashups.html),
and
[Using the Open Tree synthesis in a comparative analysis](https://cran.r-project.org/package=rotl/vignettes/meta-analysis.html).

## Quick start

### Get a little bit of the big Open Tree tree

Taxonomic names are represented in the Open Tree by numeric identifiers, the
`ott_ids` (Open Tree Taxonomy identifiers). To extract a portion of a tree from
the Open Tree, you first need to find `ott_ids` for a set of names using the
`tnrs_match_names` function:

```{r resolve}
library(rotl)
apes <- c("Pongo", "Pan", "Gorilla", "Hoolock", "Homo")
(resolved_names <- tnrs_match_names(apes))
```

Now we can get the tree with just those tips:

```{r get_tr}
tr <- tol_induced_subtree(ott_ids = ott_id(resolved_names))
plot(tr)
```

The code above can be summarized in a single pipe:

```{r pipe, plot=FALSE}
library(magrittr)
## or expressed as a pipe:
c("Pongo", "Pan", "Gorilla", "Hoolock", "Homo") %>%
  tnrs_match_names() %>%
  ott_id() %>%
  tol_induced_subtree() %>%
  plot()
```

## Citation and Manuscript

To cite `rotl` in publications pleases use:

>  Michonneau, F., Brown, J. W. and Winter, D. J. (2016), rotl: an R package to
>  interact with the Open Tree of Life data.  Methods in Ecology and
>  Evolution. 7(12):1476-1481. doi:
>  [10.1111/2041-210X.12593](https://doi.org/10.1111/2041-210X.12593)

You may also want to cite the paper for the Open Tree of Life

>  Hinchliff, C. E., et al. (2015). Synthesis of phylogeny and taxonomy into a
>  comprehensive tree of life. Proceedings of the National Academy of Sciences
>  112.41 (2015): 12764-12769
>  doi: [10.1073/pnas.1423041112](https://doi.org/10.1073/pnas.1423041112)

The manuscript in *Methods in Ecology and Evolution* includes additional
examples on how to use the package. The manuscript and the code it contains are
also hosted on GitHub at: https://github.com/fmichonneau/rotl-ms


## Versioning

Starting with v3.0.0 of the package, the major and minor version numbers (the
first 2 digits of the version number) will be matched to those of the API. The
patch number (the 3rd digit of the version number) will be used to reflect
bug fixes and other changes that are independent from changes to the API.

`rotl` can be used to access other versions of the API (if they are available)
but most likely the high level functions will not work. Instead, you will need
to parse the output yourself using the "raw" returns from the unexported
low-level functions (all prefixed with a `.`). For instance to use the
`tnrs/match_names` endpoint for `v2` of the API:

```{r versioning, eval=FALSE}
rotl:::.tnrs_match_names(c("pan", "pango", "gorilla", "hoolock", "homo"), otl_v = "v2")
```


### Code of Conduct

Please note that this project is released with a
[Contributor Code of Conduct](https://github.com/ropensci/rotl/blob/master/CONDUCT.md). By participating in this project you
agree to abide by its terms.

[![](https://ropensci.org/public_images/github_footer.png)](https://ropensci.org/)
---
title: "How to use rotl?"
author: "François Michonneau"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    css: vignette.css
vignette: >
  %\VignetteIndexEntry{How to use rotl?}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  \usepackage[utf8]{inputenc}
---

`rotl` provides an interface to the Open Tree of Life (OTL) API and allows users
to query the API, retrieve parts of the Tree of Life and integrate these parts
with other R packages.

The OTL API provides services to access:

* the **Tree of Life** a.k.a. TOL (the synthetic tree): a single draft tree that is
  a combination of **the OTL taxonomy** and the **source trees** (studies)
* the **Taxonomic name resolution services** a.k.a. TNRS: the methods for
  resolving taxonomic names to the internal identifiers used by the TOL and the
  GOL (the `ott ids`).
* the **Taxonomy** a.k.a. OTT (for Open Tree Taxonomy): which represents the
  synthesis of the different taxonomies used as a backbone of the TOL when no
  studies are available.
* the **Studies** containing the source trees used to build the TOL, and
  extracted from the scientific literature.

In `rotl`, each of these services correspond to functions with different
prefixes:

| Service       | `rotl` prefix |
|---------------|---------------|
| Tree of Life  | `tol_`        |
| TNRS          | `tnrs_`       |
| Taxonomy      | `taxonomy_`   |
| Studies       | `studies_`    |

`rotl` also provides a few other functions and methods that can be used to
extract relevant information from the objects returned by these functions.


## Demonstration of a basic workflow

The most common use for `rotl` is probably to start from a list of species and
get the relevant parts of the tree for these species. This is a two step
process:

1. the species names need to be matched to their `ott_id` (the Open Tree
	Taxonomy identifiers) using the Taxonomic name resolution services (TNRS)
1. these `ott_id` will then be used to retrieve the relevant parts of the Tree
   of Life.

### Step 1: Matching taxonomy to the `ott_id`

Let's start by doing a search on a diverse group of taxa: a tree frog (genus
_Hyla_), a fish (genus _Salmo_), a sea urchin (genus _Diadema_), and a nautilus
(genus _Nautilus_).

```{r}
library(rotl)
taxa <- c("Hyla", "Salmo", "Diadema", "Nautilus")
resolved_names <- tnrs_match_names(taxa)
```

It's always a good idea to check that the resolved names match what you
intended:

`r knitr::kable(resolved_names)`

The column `unique_name` sometimes indicates the higher taxonomic level
associated with the name. The column `number_matches` indicates the number of
`ott_id` that corresponds to a given name. In this example, our search on
_Diadema_ returns 2 matches, and the one returned by default is indeed the sea
urchin that we want for our query. The argument `context_name` allows you to
limit the taxonomic scope of your search. _Diadema_ is also the genus name of a
fungus. To ensure that our search is limited to animal names, we could do:

```{r}
resolved_names <- tnrs_match_names(taxa, context_name = "Animals")
```

If you are trying to build a tree with deeply divergent taxa that the argument
`context_name` cannot fix, see "How to change the ott ids assigned to my taxa?"
in the FAQ below.


### Step 2: Getting the tree corresponding to our taxa

Now that we have the correct `ott_id` for our taxa, we can ask for the tree
using the `tol_induced_subtree()` function. By default, the object returned by
`tol_induced_subtree` is a phylo object (from the
[ape](https://cran.r-project.org/package=ape) package), so we can plot it
directly.

```{r, fig.width=7, fig.height=4}
my_tree <- tol_induced_subtree(ott_ids = resolved_names$ott_id)
plot(my_tree, no.margin = TRUE)
```


## FAQ

### How to change the ott ids assigned to my taxa?

If you realize that `tnrs_match_names` assigns the incorrect taxonomic group to
your name (e.g., because of synonymy) and changing the `context_name` does not
help, you can use the function `inspect`. This function takes the object
resulting from `tnrs_match_names()`, and either the row number, the taxon name
(you used in your search in lowercase), or the `ott_id` returned by the initial
query.

To illustrate this, let's re-use the previous query but this time pretending that
we are interested in the fungus _Diadema_ and not the sea urchin:

```{r}
taxa <- c("Hyla", "Salmo", "Diadema", "Nautilus")
resolved_names <- tnrs_match_names(taxa)
resolved_names
inspect(resolved_names, taxon_name = "diadema")
```

In our case, we want the second row in this data frame to replace the
information that initially matched for _Diadema_. We can now use the `update()`
function, to change to the correct taxa (the fungus not the sea urchin):

```{r}
resolved_names <- update(resolved_names,
  taxon_name = "diadema",
  new_row_number = 2
)

## we could also have used the ott_id to replace this taxon:
## resolved_names <- update(resolved_names, taxon_name = "diadema",
##                          new_ott_id = 4930522)
```

And now our `resolved_names` data frame includes the taxon we want:

`r knitr::kable(resolved_names)`

### How do I know that the taxa I'm asking for is the correct one?

The function `taxonomy_taxon_info()` takes `ott_ids` as arguments and returns
taxonomic information about the taxa. This output can be passed to some helpers
functions to extract the relevant information. Let's illustrate this with our
_Diadema_ example

```{r}
diadema_info <- taxonomy_taxon_info(631176)
tax_rank(diadema_info)
synonyms(diadema_info)
tax_name(diadema_info)
```

In some cases, it might also be useful to investigate the taxonomic tree
descending from an `ott_id` to check that it's the correct taxon and to
determine the species included in the Open Tree Taxonomy:

```{r}
diadema_tax_tree <- taxonomy_subtree(631176)
diadema_tax_tree
```

By default, this function return all taxa (including self, and internal)
descending from this `ott_id` but it also possible to return `phylo` object.

### How do I get the tree for a particular taxonomic group?

If you are looking to get the tree for a particular taxonomic group, you need to
first identify it by its node id or ott id, and then use the `tol_subtree()`
function:

```{r, fig.width=7, fig.height=4}
mono_id <- tnrs_match_names("Monotremata")
mono_tree <- tol_subtree(ott_id = ott_id(mono_id))
plot(mono_tree)
```


### How do I find trees from studies focused on my favourite taxa?

The function `studies_find_trees()` allows the user to search for studies
matching a specific criteria. The function `studies_properties()` returns the
list of properties that can be used in the search.

```{r}
furry_studies <- studies_find_studies(property = "ot:focalCladeOTTTaxonName", value = "Mammalia")
furry_ids <- furry_studies$study_ids
```

Now that we know the `study_id`, we can ask for the meta data information
associated with this study:

```{r}
furry_meta <- get_study_meta("pg_2550")
get_publication(furry_meta) ## The citation for the source of the study
get_tree_ids(furry_meta) ## This study has 10 trees associated with it
candidate_for_synth(furry_meta) ## None of these trees are yet included in the OTL
```

Using `get_study("pg_2550")` would returns a `multiPhylo` object (default) with
all the trees associated with this particular study, while
`get_study_tree("pg_2550", "tree5513")` would return one of these trees.

### The tree returned by the API has duplicated tip labels, how can I work around it?

You may encounter the following error message:

```
Error in rncl(file = file, ...) : Taxon number 39 (coded by the token Pratia
angulata) has already been encountered in this tree. Duplication of taxa in a
tree is prohibited.
```

This message occurs as duplicate labels are not allowed in the NEXUS format and
it is stricly enforced by the part of the code used by `rotl` to import the
trees in memory.

If you use a version of `rotl` more recent than 0.4.1, this should not happen by
default for the function `get_study_tree`. If it happens with another function,
please [let us know](https://github.com/ropensci/rotl/issues).

The easiest way to work around this is to save the tree in a file, and use APE
to read it in memory:

```{r, eval=FALSE}
get_study_tree(
  study_id = "pg_710", tree_id = "tree1277",
  tip_label = "ott_taxon_name", file = "/tmp/tree.tre",
  file_format = "newick"
)
tr <- ape::read.tree(file = "/tmp/tree.tre")
```

### How do I get the higher taxonomy for a given taxa?

If you encounter a taxon name you are not familiar with, it might be useful to
obtain its higher taxonomy to see where it fits in the tree of life. We can
combine several taxonomy methods to extract this information easily.

```{r}
giant_squid <- tnrs_match_names("Architeuthis")
tax_lineage(taxonomy_taxon_info(ott_id(giant_squid), include_lineage = TRUE))
```

### Why are OTT IDs discovered with `rotl` missing from an induced subtree?

Some taxonomic names that can be retrieved through the taxonomic name
resolution service are not part of the Open Tree's synthesis tree. These are
usually traditional higher-level taxa that have been found to be paraphyletic.

For instance, if you wanted to fetch a tree relating the three birds that go
into a [Turkducken](https://en.wikipedia.org/wiki/Turducken) as well as the pork
used for stuffing, you might search for the turkey, duck, chicken, and pork
genera:

```{r}
turducken <- c("Meleagris", "Anas", "Gallus", "Sus")
taxa <- tnrs_match_names(turducken, context = "Animals")
taxa
```

We have the OTT ids for each genus, however, if we tried to get the induced
subtree from these results, we would get an error:

```{r, error=TRUE}
tr <- tol_induced_subtree(ott_id(taxa))
```

As the error message suggests, some of the taxa are not found in the synthetic
tree. This occurs for 2 main reasons: either the taxa is invalid, or it is part
of a group that is not monophyletic in the synthetic tree. There are two ways to
get around this issue: (1) removing the taxa that are not part of the Open Tree;
(2) using the complete species name.

#### Removing the taxa missing from the synthetic tree

To help with this situation, `rotl` provides a way to identify the OTT ids that
are not part of the synthetic tree. The function `is_in_tree()` takes the output
of the `ott_id()` function and returns a vector of logical indicating whether
the taxa are part of the synthetic tree. We can then use to only keep the taxa that appear in the synthetic tree:

```{r}
in_tree <- is_in_tree(ott_id(taxa))
in_tree
tr <- tol_induced_subtree(ott_id(taxa)[in_tree])
```

#### Using the full taxonomic names

The best way to avoid these problems is to specify complete species names
(species being the lowest level of classification in the Open Tree taxonomy they
are guaranteed to be monophyletic):

```{r, fig.width=7, fig.height=4}
turducken_spp <- c("Meleagris gallopavo", "Anas platyrhynchos", "Gallus gallus", "Sus scrofa")
taxa <- tnrs_match_names(turducken_spp, context = "Animals")
tr <- tol_induced_subtree(ott_id(taxa))
plot(tr)
```
---
title: "Connecting data to Open Tree trees"
author: "David Winter"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    css: vignette.css
vignette: >
  %\VignetteIndexEntry{Connecting data to Open Tree trees}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  \usepackage[utf8]{inputenc}
---

## Combining data from OToL and other sources.

One of the major goals of `rotl` is to help users combine data from other
sources with the phylogenetic trees in the Open Tree database. This examples
document describes some of the ways in which a user might connect data to trees
from Open Tree.

## Get Open Tree IDs to match your data.

Let's say you have a dataset where each row represents a measurement taken from
one species, and your goal is to put these measurements in some phylogenetic
context. Here's a small example: the best estimate of the mutation rate for a
set of unicellular Eukaryotes along with some other property of those species
which might explain the mutation rate:

```{r, data}
csv_path <- system.file("extdata", "protist_mutation_rates.csv", package = "rotl")
mu <- read.csv(csv_path, stringsAsFactors = FALSE)
mu
```

If we want to get a tree for these species we need to start by finding the
unique ID for each of these species in the Open Tree database. We can use the
Taxonomic Name Resolution Service (`tnrs`) functions to do this. Before we do
that we should see if any of the taxonomic contexts, which can be used to narrow
a search and avoid conflicts between different codes, apply to our group of species:

```{r, context}
library(rotl)
tnrs_contexts()
```

Hmm, none of those groups contain all of our species. In this case we can
search using the `All life` context and the function `tnrs_match_names`:

```{r, match}
taxon_search <- tnrs_match_names(names = mu$species, context_name = "All life")
knitr::kable(taxon_search)
```

Good, all of our  species are known to Open Tree. Note, though, that one of the names
is a synonym. _Saccharomyces pombe_ is older name for what is now called
_Schizosaccharomyces pombe_. As the name suggests, the Taxonomic Name
Resolution Service is designed to deal with these problems (and similar ones
like misspellings), but it is always a good idea to check the results of
`tnrs_match_names` closely to ensure the results are what you expect.

In this case we have a good ID for each of our species so we can move on. Before
we do that, let's ensure we can match up our original data to the Open Tree
names and IDs by adding them to our `data.frame`:

```{r, munge}
mu$ott_name <- unique_name(taxon_search)
mu$ott_id <- taxon_search$ott_id
```

## Find a tree with your taxa

Now let's find a tree. There are two possible options here: we can search for
published studies that include our taxa or we can use the 'synthetic tree' from
Open Tree. We can try both approaches.

### Published trees

Before we can search for published studies or trees, we should check out the
list of properties we can use to perform such searches:

```{r, properties}
studies_properties()
```

We have `ottIds` for our taxa, so let's use those IDs to search for trees that
contain them.  Starting with our first species _Tetrahymena thermophila_ we can
use `studies_find_trees` to do this search.

```{r taxon_count}
studies_find_trees(property = "ot:ottId", value = as.character(ott_id(taxon_search)[1]))
```

Well... that's not very promising. We can repeat that process for all of the IDs
to see if the other species are better represented.


```{r, all_taxa_count}
hits <- lapply(mu$ott_id, studies_find_trees, property = "ot:ottId", detailed = FALSE)
sapply(hits, function(x) sum(x[["n_matched_trees"]]))
```

OK, most of our species are not in any of the published trees available. You can
help fix this sort of problem by [making sure you submit your published trees to
Open Tree](https://tree.opentreeoflife.org/curator).

### A part of the synthesis tree

Thankfully, we can still use the complete Tree of Life made from the
combined results of all of the published trees and taxonomies that go into Open
Tree. The function `tol_induced_subtree` will fetch a tree relating a set of IDs.

Using the default arguments you can get a tree object into your R session:


```{r subtree,  fig.width=7, fig.height=4}
ott_in_tree <- ott_id(taxon_search)[is_in_tree(ott_id(taxon_search))]
tr <- tol_induced_subtree(ott_ids = ott_in_tree)
plot(tr)
```

### Connect your data to the tips of your tree

Now we have a tree for of our species, how can we use the tree and the data
together?

The package `phylobase` provide an object class called `phylo4d`, which is
designed to represent a phylogeny and data associated with its tips. In oder to
get our tree and data into one of these objects we have to make sure the labels
in the tree and in our data match exactly. That's not quite the case at the
moment (tree labels have underscores and IDs appended):

```{r, match_names}
mu$ott_name[1]
tr$tip.label[4]
```

`rotl` provides a convienence function `strip_ott_ids` to deal with these.

```{r, sub}
tr$tip.label <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
tr$tip.label %in% mu$ott_name
```

Ok, now the tips are together we can make a new dataset. The `phylo4d()`
functions matches tip labels to the row names of a `data.frame`, so let's make
a new dataset that contains just the relevant data and has row names to match
the tree

```{r phylobase}
library(phylobase)
mu_numeric <- mu[, c("mu", "pop.size", "genome.size")]
rownames(mu_numeric) <- mu$ott_name
tree_data <- phylo4d(tr, mu_numeric)
```
And now we can plot the data and the tree together


```{r,  fig.width=7, fig.height=5}
plot(tree_data)
```

## Find external data associated with studies, trees and taxa from Open Tree

In the above example we looked for a tree that related species in another dataset.
Now we will go the other way, and try to find data associated with Open Tree records
in other databases.

### Get external data from a study

Let's imagine you were interested in extending or reproducing the results of a
published study. If that study is included in Open Tree you can find it via
`studies_find_studies` or `studies_find_trees` and retrieve the published trees
with `get_study`. `rotl` will also help you find external. The function
`study_external_IDs` retrieves the DOI for a given study, and uses that to
gather some more data:

```{r}
extra_data <- try(study_external_IDs("pg_1980"), silent = TRUE)
if (!inherits(extra_data, "try-error")) {
  extra_data
}
```

Here the returned object contains an `external_data_url` (in this case a link to
the study in Treebase), a pubmed ID for the paper and a vector IDs for the
NCBI's nuleotide database. The packages `treebase` and `rentrez` provide
functions to make use of these IDs within R.

As an example, let's use `rentrez` to download the first two DNA seqences and
print them.

```{r}
library(rentrez)
seqs <- try(entrez_fetch(db = "nucleotide", id = extra_data$nucleotide_ids[1:2], rettype = "fasta"), silent = TRUE)

if (inherits(seqs, "try-error")) {
  cat("NCBI temporarily down.")
} else {
  cat(seqs)
}
```

You could further process these sequences in R with the function `read.dna` from
`ape` or save them to disk by specifying a file name with `cat`.

### Find a OTT taxon in another taxonomic database

It is also possible map an Open Tree taxon to a record in another taxonomic
database. For instance, if we wanted to search for data about one of the tips of
the sub-tree we fetched in the example above we could do so using
`taxon_external_IDs`:

```{r}
Tt_ids <- taxon_external_IDs(mu$ott_id[2])
Tt_ids
```

A user could then use `rgbif` to find locality records using the gbif ID or
`rentrez` to get genetic or bibliometric data about from the NCBI's databases.


## What next

The demonstration gets you to the point of visualizing your data in a
phylogenetic context. But there's a lot more you do with this sort of data in R.
For instance, you could use packages like `ape`, `caper`, `phytools` and
`mcmcGLMM` to perform phylogenetic comparative analyses of your data. You could
gather more data on your species using packages that connect to
trait databases like `rfishbase`, `AntWeb` or `rnpn` which provides data from
the US National Phenology Network. You could also use `rentrez` to find genetic
data for each of your species, and use that data to generate branch lengths for
the phylogeny.
---
title: "Using the Open Tree synthesis in a comparative analysis"
author: "David Winter"
date: "`r Sys.Date()`"
output:
  rmarkdown::html_vignette:
    css: vignette.css
vignette: >
  %\VignetteIndexEntry{Using the Open Tree synthesis in a comparative analysis}
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteEncoding{UTF-8}
  %\VignetteEncoding{UTF-8}
---

## Phylogenetic Comparative Methods

The development of phylogenetic comparative methods has made phylogenies and
important source of data in fields as diverse as ecology, genomic and medicine.
Comparative  methods can be used to investigate patterns in the evolution of
traits or the diversification of lineages. In other cases a phylogeny is treated
as a "nuisance parameter", allowing with the autocorrelation created by the shared
evolutionary history of the different species included to be controlled for.

In many cases finding a tree that relates the species for which trait data are
available is a rate-limiting step in such comparative analyses. Here we show
how the synthetic tree provided by Open Tree of Life (and made available in R via
`rotl`) can help to fill this gap.

## A phylogenetic meta-analysis

To demonstrate the use of `rotl` in a comparative analysis, we will partially
reproduce the results of [Rutkowska _et al_ 2014](https://doi.org/10.1111/jeb.12282).
Very briefly, this study is a meta-analysis summarising the results of multiple
studies testing for systematic differences in the size of eggs which contain
male and female offspring. Such a difference might mean that birds invest more
heavily in one sex than the other.

Because this study involves data from 51 different species, Rutkowska _et al_
used a phylogenetic comparative approach to account for the shared evolutionary
history among some of the studied-species.

### Gather the data

If we are going to reproduce this analysis, we will first need to gather the
data. Thankfully, the data is available as supplementary material from the
publisher's website. We provide a copy of this data with the package:

```{r load-package}
library(rotl)
```


```{r egg_data}
## This dataset is available from the publisher's study website:
egg_data <- read.csv(system.file("extdata", "egg.csv", package = "rotl"),
  stringsAsFactors = FALSE
)
## }
head(egg_data)
```

The most important variable in this dataset is `Zr`, which is a [normalized
effect size](https://en.wikipedia.org/wiki/Fisher_transformation) for difference
,in size between eggs that contain males and females. Values close to zero come
from studies that found the sex of an egg's inhabitant had little effect in its size,
while large positive or negative values correspond to studies with substantial
sex biases (towards males and females respectively). Since this is a
meta-analysis we should produce the classic [funnel plot](https://en.wikipedia.org/wiki/Funnel_plot)
with effects-size on the y-axis and precision (the inverse of the sample
standard error) on the x-axis. Here we calculate precision from the sample
variance (`Vzr`):

```{r eggs_in_a_funnel, fig.width=6, fig.height=3}
plot(1 / sqrt(egg_data$VZr), egg_data$Zr,
  pch = 16,
  ylab = "Effect size (Zr)",
  xlab = "Precision (1/SE)",
  main = "Effect sizes for sex bias in egg size among 51 brid species"
)
```

In order to use this data later on we need to first convert it to a standard
`data.frame`. We can also convert the `animal` column (the species names) to
lower case, and remove the underscores in their names, which will make it easier to match names later on:

```{r, clean_eggs}
egg_data <- as.data.frame(egg_data)
## Convert taxon names to lower case
egg_data$animal <- tolower(egg_data$animal)
## Let's remove the underscores (_) from the taxon names
egg_data$animal <- gsub("_", " ", egg_data$animal)
```
### Find the species in OTT

We can use the OTL synthesis tree to relate these species. To do so we first need to
find Open Tree Taxonomy (OTT) IDs for each species. We can do that with the
Taxonomic Name Resolution Service function `tnrs_match_names`:

```{r birds}
taxa <- tnrs_match_names(unique(egg_data$animal), context = "Animals")
head(taxa)
```

All of these species are in OTT, but a few of them go by different names in the
Open Tree than we have in our data set. Because the tree `rotl` fetches
will have Open Tree names, we need to create a named vector that maps the names
we have for each species to the names Open Tree uses for them:


```{r bird_map}
taxon_map <- structure(taxa$search_string, names = taxa$unique_name)
```

Now we can use this map to retrieve "data set names" from "OTT names":


```{r odd_duck}
taxon_map["Anser caerulescens"]
```

### Get a tree

Now we can get the tree. There are really too many tips here to show nicely, so
we will leave them out of this plot

```{r birds_in_a_tree, fig.width=5, fig.height=5, fig.align='center'}
tr <- tol_induced_subtree(ott_id(taxa)[is_in_tree(ott_id(taxa))])
plot(tr, show.tip.label = FALSE)
```

There are a few things to note here. First, the tree has no branch lengths.
At present this is true for the whole of the Open Tree synthetic tree. Some
comparative methods require either branch lengths or an ultrametric tree. Before
you can use one of those methods you will need to get a tree with branch
lengths. You could try looking for published trees made available by the Open
Tree with `studies_find_trees`. Alternatively, you could estimate branch lengths
from the toplogy of a phylogeny returned by `tol_induced_subtree`, perhaps by
downloading DNA sequences from the NCBI with `rentrez` or "hanging" the tree on
nodes of known-age using  penalized likelihood method in `ape::chronos`.
In this case, we will use only the topology of the tree as input to our
comparative analysis, so we can skip these steps.

Second, the tip labels contain OTT IDs, which means they will not perfectly
match the species names in our dataset or the taxon map that we created earlier:


```{r tip_lab}
tr$tip.label[1:4]
```

Finally, the tree contains node labels for those nodes that match a higher taxonomic
group, and empty character vectors (`""`) for all other nodes. Some
comparative methods either do no expect node labels at all, or require all
labeled nodes to have a unique name (meaning multiple "empty" labels will cause
and error).

We can deal with all these details easily. `rotl` provides  the convenience
function `strip_ott_ids` to remove the extra information from the tip labels.
With the IDs removed, we can use our taxon map to replace the tip labels in the tree
with the species names from dataset.



```{r clean_tips}
otl_tips <- strip_ott_ids(tr$tip.label, remove_underscores = TRUE)
tr$tip.label <- taxon_map[ otl_tips ]
```

Finally, we can remove the node labels by setting the `node.label` attribute of
the tree to `NULL`.

```{r remove_nodes}
tr$node.label <- NULL
```

```{r match_species_tree}
egg_data <- egg_data[egg_data$animal %in% tr$tip.label, ]
```


### Perform the meta-analysis


Now we have data and a tree, and we know the names in the tree match the ones in
the data. It's time to do the comparative analysis. Rutkowska _et al_. used `MCMCglmm`, a
Bayesian MCMC approach to fitting multi-level models,to perform their meta-analysis,
and we will do the same. Of course, to properly analyse these data you would
take some care in deciding on the appropriate priors to use and inspect the
results carefully. In this case, we are really interested in using this as a
demonstration, so we will just run a simple model.

Specifically we sill fit a model where the only variable that might explain the
values of `Zr` is the random factor `animal`, which corresponds to the
phylogenetic relationships among species. We also provide `Zvr` as the measurement
error variance, effectively adding extra weight to the results of more powerful
studies. Here's how we specify and fit that model with `MCMCglmm`:


```{r model}
set.seed(123)
if (require(MCMCglmm, quietly = TRUE)) {
  pr <- list(
    R = list(V = 1, nu = 0.002),
    G = list(G1 = list(V = 1, nu = 0.002))
  )

  model <- MCMCglmm(Zr ~ 1,
    random = ~animal,
    pedigree = tr,
    mev = egg_data$VZr,
    prior = pr,
    data = egg_data,
    verbose = FALSE
  )
} else {
  model <- readRDS(file = system.file("extdata", "mcmcglmm_model.rds", package = "rotl"))
}
```


Now that we have a result we can find out how much phylogenetic signal exists
for sex-biased differences in egg-size. In a multi-level model we can use variance
components to look at this, specifically the proportion of the total variance
that can be explained by phylogeny is called the phylogenetic reliability, _H_. Let's
calculate the _H_ for this model:


```{r PhyH}
var_comps <- colMeans(model$VCV)
var_comps["animal"] / sum(var_comps)
```

It appears there is almost no phylogenetic signal to the data.
The relationships among species explain much less that one percent of the total
variance in the data. If you were wondering,  Rutkowska _et al_. report a similar result,
even after adding more predictors to their model most of the variance in `Zr`
was left unexplained.

## What other comparative methods can I use in R?

Here we have demonstrated just one comparative analysis that you might do in R.
There are an ever-growing number of packages that allow an ever-growing number
of analysis to performed in R. Some "classics" like ancestral state
reconstruction,  phylogenetic independent contrasts and lineage through time plots
are implemented in `ape`. Packages like `phytools`, `caper` and `diversitree`
provide extensions to these methods.  The [CRAN Phylogenetics Taskview](https://CRAN.R-project.org/view=Phylogenetics)
gives a good idea of the diversity of packages and analyses that can be
completed in R.
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R, R/tol.R
\name{tol_lineage}
\alias{tol_lineage}
\alias{tol_node_info}
\alias{tax_rank.tol_node}
\alias{tax_sources.tol_node}
\alias{unique_name.tol_node}
\alias{tax_name.tol_node}
\alias{ott_id.tol_node}
\alias{source_list.tol_node}
\alias{tax_lineage.tol_node}
\alias{tol_lineage.tol_node}
\title{Node info}
\usage{
tol_lineage(tax, ...)

tol_node_info(ott_id = NULL, node_id = NULL, include_lineage = FALSE, ...)

\method{tax_rank}{tol_node}(tax, ...)

\method{tax_sources}{tol_node}(tax, ...)

\method{unique_name}{tol_node}(tax, ...)

\method{tax_name}{tol_node}(tax, ...)

\method{ott_id}{tol_node}(tax, ...)

\method{source_list}{tol_node}(tax, ...)

\method{tax_lineage}{tol_node}(tax, ...)

\method{tol_lineage}{tol_node}(tax, ...)
}
\arguments{
\item{tax}{an object returned by \code{tol_node_info}.}

\item{...}{additional arguments to customize the API call (see
?rotl for more information)}

\item{ott_id}{Numeric. The OpenTree taxonomic identifier.}

\item{node_id}{Character. The OpenTree node identifier.}

\item{include_lineage}{Logical (default = FALSE). Whether to return the
lineage of the node from the synthetic tree.}
}
\value{
\code{tol_node_info} returns an invisible list of summary
    information about the queried node:

\itemize{

    \item {node_id} {String. The canonical identifier of the node.}

    \item {num_tips} {Numeric. The number of descendent tips.}

    \item {partial_path_of} {List. The edge below this synthetic tree node
    is compatible with the edge below each of these input tree nodes (one
    per tree). Each returned element is reported as sourceid:nodeid.}

   \item {query} { The node id that resolved to this node. This can differ
   from the node_id field if the query id is not canonical. }

    \item {taxon} {A list of taxonomic properties. Only returned if
    the queried node is a taxon. Each source has:}

        \itemize{
            \item {ott_id} {Numeric. The OpenTree Taxonomy ID (ottID).}

            \item {name} {String. The taxonomic name of the queried node.}

            \item {unique_name} {String. The string that uniquely
            identifies the taxon in OTT.}

            \item {rank} {String. The taxonomic rank of the taxon in OTT.}

            \item {tax_sources} {List. A list of identifiers for taxonomic
            sources, such as other taxonomies, that define taxa judged
            equivalent to this taxon.}
        }

    The following properties list support/conflict for the node across
    synthesis source trees. All properties involve sourceid keys and
    nodeid values (see \code{source_id_map} below).

    \item {supported_by} {List. Input tree nodes (one per tree) that support
    this synthetic tree node. Each returned element is reported as
    sourceid:nodeid.}

    \item {terminal} {List. Input tree nodes (one per tree) that are
    equivalent to this synthetic tree node (via an exact mapping, or the
    input tree terminal may be the only terminal descended from this
    synthetic tree node. Each returned element is reported as
    sourceid:nodeid.}

    \item {conflicts_with} {Named list of lists. Names correspond to
    sourceid keys. Each list contains input tree node ids (one or more per
    tree) that conflict with this synthetic node.}
  }

    \code{tol_lineage} and \code{tax_lineage} return data
        frames. \code{tol_lineage} indicate for each ancestor its
        node identifier, the number of tips descending from that
        node, and whether it corresponds to a taxonomic level.
}
\description{
Get summary information about a node in the synthetic tree
}
\details{
Returns summary information about a node in the graph. The
    node of interest may be specified using either a node id or an
    taxon id, but not both. If the specified node or OTT id is not
    in the graph, an error will be returned.

    If the argument \code{include_lineage=TRUE} is used, you can
    use \code{tax_lineage()} or \code{tol_lineage} to return the
    taxonomic information or the node information for all the
    ancestors to this node, down to the root of the tree.
}
\examples{
\dontrun{
birds <- tol_node_info(ott_id=81461, include_lineage=TRUE)
source_list(birds)
tax_rank(birds)
ott_id(birds)
tax_lineage(birds)
tol_lineage(birds)}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studies-methods.R, R/studies.R
\name{get_tree_ids}
\alias{get_tree_ids}
\alias{get_publication}
\alias{candidate_for_synth}
\alias{get_study_year}
\alias{get_tree_ids.study_meta}
\alias{get_publication.study_meta}
\alias{candidate_for_synth.study_meta}
\alias{get_study_year.study_meta}
\alias{get_study_meta}
\title{Study Metadata}
\usage{
get_tree_ids(sm)

get_publication(sm)

candidate_for_synth(sm)

get_study_year(sm)

\method{get_tree_ids}{study_meta}(sm)

\method{get_publication}{study_meta}(sm)

\method{candidate_for_synth}{study_meta}(sm)

\method{get_study_year}{study_meta}(sm)

get_study_meta(study_id, ...)
}
\arguments{
\item{sm}{an object created by \code{get_study_meta}}

\item{study_id}{the study identifier (character)}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
named-list containing the metadata associated with the
    study requested
}
\description{
Retrieve metadata about a study in the Open Tree of Life datastore.
}
\details{
\code{get_study_meta} returns a long list of attributes for the
studies that are contributing to the synthetic tree. To help with
the extraction of relevant information from this list, several
helper functions exists: \itemize{

  \item {get_tree_ids} { The identifiers of the trees
  associated with the study }

  \item {get_publication} { The citation information of the
  publication for the study. The DOI (or URL) for the study is
  available as an attribute to the returned object (i.e.,
  \code{attr(object, "DOI")} ) }.

  \item {candidate_for_synth} { The identifier of the tree(s) from
  the study used in the synthetic tree. This is a subset of the
  result of \code{get_tree_ids}.

  \item {get_study_year} { The year of publication of the study. }

  }
}
}
\examples{
\dontrun{
req <- get_study_meta("pg_719")
get_tree_ids(req)
candidate_for_synth(req)
get_publication(req)
get_study_year(req)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_names.R
\name{inspect.match_names}
\alias{inspect.match_names}
\alias{inspect}
\alias{update.match_names}
\title{Inspect and Update alternative matches for a name returned
    by tnrs_match_names}
\usage{
\method{inspect}{match_names}(response, row_number, taxon_name, ott_id, ...)

inspect(response, ...)

\method{update}{match_names}(object, row_number, taxon_name, ott_id, new_row_number, new_ott_id, ...)
}
\arguments{
\item{response}{an object generated by the
\code{\link{tnrs_match_names}} function}

\item{row_number}{the row number corresponding to the name to
inspect}

\item{taxon_name}{the taxon name corresponding to the name to
inspect}

\item{ott_id}{the ott id corresponding to the name to inspect}

\item{...}{currently ignored}

\item{object}{an object created by \code{\link{tnrs_match_names}}}

\item{new_row_number}{the row number in the output of
\code{\link{inspect}} to replace the taxa specified by
\code{row_number}, \code{taxon_name}, or \code{ott_id}.}

\item{new_ott_id}{the ott id of the taxon to replace the taxa
specified by \code{row_number}, \code{taxon_name}, or
\code{ott_id}.}
}
\value{
a data frame
}
\description{
Taxonomic names may have different meanings in different taxonomic
contexts, as the same genus name can be applied to animals and
plants for instance. Additionally, the meaning of a taxonomic name
may have change throughout its history, and may have referred to a
different taxon in the past. In such cases, a given names might
have multiple matches in the Open Tree Taxonomy. These functions
allow users to inspect (and update) alternative meaning of a given
name and its current taxonomic status according to the Open Tree
Taxonomy.
}
\details{
To inspect alternative taxonomic meanings of a given name, you
need to provide the object resulting from a call to the
tnrs_match_names function, as well as one of either the row number
corresponding to the name in this object, the name itself (as used
in the original query), or the ott_id listed for this name.

To update one of the name, you also need to provide the row number
in which the name to be replaced appear or its ott id.
}
\examples{
  \dontrun{
   matched_names <- tnrs_match_names(c("holothuria", "diadema", "boletus"))
   inspect(matched_names, taxon_name="diadema")
   new_matched_names <- update(matched_names, taxon_name="diadema",
                               new_ott_id = 631176)
   new_matched_names
   }
}
\seealso{
\code{\link{tnrs_match_names}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tnrs.R
\name{tnrs_match_names}
\alias{tnrs_match_names}
\title{Match names to the Open Tree Taxonomy}
\usage{
tnrs_match_names(
  names = NULL,
  context_name = "All life",
  do_approximate_matching = TRUE,
  ids = NULL,
  include_suppressed = FALSE,
  ...
)
}
\arguments{
\item{names}{taxon names to be queried. Currently limited to 10,000 names
for exact matches and 2,500 names for approximate matches (character
vector)}

\item{context_name}{name of the taxonomic context to be searched (length-one
character vector or \code{NULL}). Must match (case sensitive) one of the
values returned by \code{\link{tnrs_contexts}}. Default to "All life".}

\item{do_approximate_matching}{A logical indicating whether or not to
perform approximate string (a.k.a. \dQuote{fuzzy}) matching. Using
\code{FALSE} will greatly improve speed. Default, however, is \code{TRUE}.}

\item{ids}{A vector of ids to use for identifying names. These will be
assigned to each name in the names array. If ids is provided, then ids and
names must be identical in length.}

\item{include_suppressed}{Ordinarily, some quasi-taxa, such as incertae
sedis buckets and other non-OTUs, are suppressed from TNRS results. If
this parameter is true, these quasi-taxa are allowed as possible TNRS
results.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
A data frame summarizing the results of the query. The original
  query output is appended as an attribute to the returned object (and can
  be obtained using \code{attr(object, "original_response")}).
}
\description{
Match taxonomic names to the Open Tree Taxonomy.
}
\details{
Accepts one or more taxonomic names and returns information about
potential matches for these names to known taxa in the Open Tree
Taxonomy.

This service uses taxonomic contexts to disambiguate homonyms and
misspelled names; a context may be specified using the
\code{context_name} argument. If no context is specified, then the
context will be inferred (i.e., the shallowest taxonomic context
that contains all unambiguous names in the input). Taxonomic
contexts are uncontested higher taxa that have been selected to
allow limits to be applied to the scope of TNRS searches
(e.g. 'match names only within flowering plants'). Once a context
has been identified (either user-specified or inferred), all taxon
name matches will performed only against taxa within that
context. For a list of available taxonomic contexts, see
\code{\link{tnrs_contexts}}.

A name is considered unambiguous if it is not a synonym and has
only one exact match to any taxon name in the entire taxonomy.

Several functions listed in the \sQuote{See also} section can be
used to inspect and manipulate the object generated by this
function.
}
\examples{
\dontrun{
 deuterostomes <- tnrs_match_names(names=c("echinodermata", "xenacoelomorpha",
                                            "chordata", "hemichordata"))
}
}
\seealso{
\code{\link{inspect.match_names}},
  \code{\link{update.match_names}}, \code{\link{synonyms.match_names}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxonomy_subtree}
\alias{taxonomy_subtree}
\title{Taxonomy subtree}
\usage{
taxonomy_subtree(
  ott_id = NULL,
  output_format = c("taxa", "newick", "phylo", "raw"),
  label_format = NULL,
  file,
  ...
)
}
\arguments{
\item{ott_id}{The ott id of the taxon of interest.}

\item{output_format}{the format of the object to be returned. See
the \sQuote{Return} section.}

\item{label_format}{Character. Defines the label type; one of
\dQuote{\code{name}}, \dQuote{\code{id}}, or
 \dQuote{\code{name_and_id}} (the default).}

\item{file}{the file name where to save the output of the
function. Ignored unless \code{output_format} is set to
\dQuote{\code{phylo}}.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
If the \code{file} argument is missing: \itemize{

    \item{\dQuote{\code{taxa}}} { a list of the taxa names
    (species) in slot \code{tip_label}, and higher-level taxonomy
    (e.g., families, genera) in slot \code{edge_label}, descending
    from the taxa corresponding to the \code{ott_id} provided. }

    \item{\dQuote{\code{newick}}} { a character vector containing
    the newick formatted string corresponding to the taxonomic
    subtree for the \code{ott_id} provided. }

    \item{\dQuote{\code{phylo}}} { an object of the class
    \code{phylo} from the \code{ape} package. }

    \item{\dQuote{\code{raw}}} { the direct output from the API,
    i.e., a list with an element named \sQuote{newick} that
    contains the subtree as a newick formatted string. }

    }

    If a \code{file} argument is provided (and
    \code{output_format} is set to \dQuote{\code{phylo}}), a
    logical indicating whether the file was successfully created.
}
\description{
Given an ott id, return the inclusive taxonomic subtree descended
from the specified taxon.
}
\details{
If the output of this function is exported to a file, the only
possible value for the \code{output_format} argument is
\dQuote{\code{newick}}. If the file provided already exists, it
will be silently overwritten.
}
\examples{
\dontrun{
req <- taxonomy_subtree(ott_id=515698)
plot(taxonomy_subtree(ott_id=515698, output_format="phylo"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studies.R
\name{studies_find_trees}
\alias{studies_find_trees}
\title{Find Trees}
\usage{
studies_find_trees(
  property = NULL,
  value = NULL,
  verbose = FALSE,
  exact = FALSE,
  detailed = TRUE,
  ...
)
}
\arguments{
\item{property}{The property to be searched on (character)}

\item{value}{The property-value to be searched on (character)}

\item{verbose}{Should the output include all metadata? (logical,
default \code{FALSE})}

\item{exact}{Should exact matching be used for the value?
(logical, default \code{FALSE})}

\item{detailed}{Should a detailed report be provided? If
\code{TRUE} (default), the output will include metadata about
the study that include trees matching the property. Otherwise,
only information about the trees will be provided.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
A data frame that summarizes the trees found (and their
    associated studies) for the requested criteria. If a study has
    more than 5 trees, the \code{tree_ids} of the first ones will
    be shown, followed by \code{...} to indicate that more are
    present.

    If \code{detailed=FALSE}, the data frame will include the
    study ids of the study (\code{study_ids}), the number of trees
    in this study that match the search criteria
    (\code{n_matched_trees}), the tree ids that match the search
    criteria (\code{match_tree_ids}).

    If \code{detailed=TRUE}, in addition of the fields listed
    above, the data frame will also contain the total number of
    trees associated with the study (\code{n_trees}), all the tree
    ids associated with the study (\code{tree_ids}), the tree id
    that is a potential candidate for inclusion in the synthetic
    tree (if any) (\code{candidate}), the year the study was
    published (\code{study_year}), the title of the study
    (\code{title}), the DOI for the study (\code{study_doi}).
}
\description{
Return a list of studies for which trees match a given set of
properties
}
\details{
The list of possible values to be used as values for the argument
\code{property} can be found using the function
\code{\link{studies_properties}}.
}
\examples{
\dontrun{
res <- studies_find_trees(property="ot:ottTaxonName", value="Drosophila",
                          detailed=FALSE)
## summary of the trees and associated studies that match this criterion
res
## With metadata about the studies (default)
res <- studies_find_trees(property="ot:ottTaxonName", value="Drosophila",
                          detailed=TRUE)
## The list of trees for each study that match the search criteria
list_trees(res)
## the trees for a given study
list_trees(res, study_id = "pg_2769")
}
}
\seealso{
\code{\link{studies_properties}} which lists properties
  the studies can be searched on. \code{\link{list_trees}} for
  listing the trees that match the query.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/external_data.R
\name{study_external_IDs}
\alias{study_external_IDs}
\title{Get external identifiers for data associated with an Open Tree study}
\usage{
study_external_IDs(study_id)
}
\arguments{
\item{study_id}{An open tree study ID}
}
\value{
A study_external_data object (which inherits from a list) which
contains some of the following.

doi, character, the DOI for the paper describing this study

external_data_url, character, a URL to an external data repository
(e.g. a treebase entry) if one exists.

pubmed_id character, the unique ID for this study in the NCBI's pubmed database

popset_ids character, vector of IDs for the NCBI's popset database

nucleotide_ids character, vector of IDs for the NCBI's nucleotide database
}
\description{
Data associated with studies contributing to the Open Tree synthesis may
be available from other databases. In particular, trees and alignments
may be available from treebase and DNA sequences and bibliographic
information associated with a given study may be available from the NCBI.
This function retrieves that information for a given study.
}
\examples{
\dontrun{
flies <- studies_find_studies(property="ot:focalCladeOTTTaxonName", value="Drosophilidae")
study_external_IDs(flies[2,]$study_ids)
}
}
\seealso{
studies_find_studies (used to discover study IDs)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_names.R
\name{synonyms.match_names}
\alias{synonyms.match_names}
\title{List the synonyms for a given name}
\usage{
\method{synonyms}{match_names}(tax, row_number, taxon_name, ott_id, ...)
}
\arguments{
\item{tax}{a data frame generated by the
\code{\link{tnrs_match_names}} function}

\item{row_number}{the row number corresponding to the name for
which to list the synonyms}

\item{taxon_name}{the taxon name corresponding to the name for
which to list the synonyms}

\item{ott_id}{the ott id corresponding to the name for which to
list the synonyms}

\item{...}{currently ignored}
}
\value{
a list whose elements are all synonym names (as vectors of
    character) for the taxonomic names that match the query (the
    names of the elements of the list).
}
\description{
When querying the Taxonomic Name Resolution Services for a
particular taxonomic name, the API returns as possible matches all
names that include the queried name as a possible synonym. This
function allows you to explore other synonyms for an accepted
name, and allows you to determine why the name you queried is
returning an accepted synonym.
}
\details{
To list synonyms for a given taxonomic name, you need to provide
the object resulting from a call to the
\code{\link{tnrs_match_names}} function, as well as one of either
the row number corresponding to the name in this object, the name
itself (as used in the original query), or the ott_id listed for
this name. Otherwise, the synonyms for all the currently matched
names are returned.
}
\examples{
\dontrun{
   echino <- tnrs_match_names(c("Diadema", "Acanthaster", "Fromia"))
   ## These 3 calls are identical
   synonyms(echino, taxon_name="Acanthaster")
   synonyms(echino, row_number=2)
   synonyms(echino, ott_id=337928)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studies-methods.R
\name{list_trees}
\alias{list_trees}
\alias{list_trees.matched_studies}
\title{List trees ids in objects returned by
\code{\link{studies_find_studies}} and
\code{\link{studies_find_trees}}.}
\usage{
list_trees(matched_studies, ...)

\method{list_trees}{matched_studies}(matched_studies, study_id, ...)
}
\arguments{
\item{matched_studies}{an object created by
\code{studies_find_trees} or \code{studies_find_studies}.}

\item{...}{Currently unused}

\item{study_id}{a \code{study_id} listed in the object returned by
\code{studies_find_trees}}
}
\value{
\code{list_trees} returns a list of the tree_ids for each
    study that match the requested criteria. If a \code{study_id}
    is provided, then only the trees for this study are returned
    as a vector.
}
\description{
\code{list_trees} returns all trees associated with a particular
study when used on an object returned by
\code{\link{studies_find_studies}}, but only the trees that match
the search criteria when used on objects returned by
\code{\link{studies_find_trees}}.
}
\seealso{
\code{\link{studies_find_studies}} and
    \code{\link{studies_find_trees}}. The help for these functions
    have examples demonstrating the use of \code{list_trees}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R, R/tol.R
\name{source_list}
\alias{source_list}
\alias{source_list.tol_summary}
\title{List of studies used in the Tree of Life}
\usage{
source_list(tax, ...)

\method{source_list}{tol_summary}(tax, ...)
}
\arguments{
\item{tax}{a list containing a \code{source_id_map} slot.}

\item{...}{additional arguments (currently unused)}
}
\value{
a data frame
}
\description{
Retrieve the detailed information for the list of studies used in
the Tree of Life.
}
\details{
This function takes the object resulting from
    \code{tol_about(study_list = TRUE)}, \code{tol_mrca()},
    \code{tol_node_info()}, and returns a data frame listing the
    \code{tree_id}, \code{study_id} and \code{git_sha} for the
    studies currently included in the Tree of Life.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tol.R
\name{tol_induced_subtree}
\alias{tol_induced_subtree}
\title{Subtree from the Open Tree of Life}
\usage{
tol_induced_subtree(
  ott_ids = NULL,
  node_ids = NULL,
  label_format = NULL,
  file,
  ...
)
}
\arguments{
\item{ott_ids}{Numeric vector. OTT ids indicating nodes to be used as tips
in the induced tree.}

\item{node_ids}{Character vector. Node ids indicating nodes to be used as
tips in the induced tree.}

\item{label_format}{Character. Defines the label type; one of
\dQuote{\code{name}}, \dQuote{\code{id}}, or \dQuote{\code{name_and_id}}
(the default).}

\item{file}{If specified, the function will write the subtree to a file in
newick format.}

\item{...}{additional arguments to customize the API call (see
\code{\link{rotl}} for more information).}
}
\value{
If no value is specified to the \code{file} argument
    (default), a phylogenetic tree of class \code{phylo}.

    Otherwise, the function returns invisibly a logical indicating
    whether the file was successfully created.

    Note that the tree returned when specifying a file name with the
    \code{file} argument is the \dQuote{raw} Newick string returned by Open
    Tree of Life. This string contains singleton nodes, and therefore will
    be different from the tree returned as a \code{phylo} object which will
    not contain these singleton nodes.
}
\description{
Return the induced subtree on the synthetic tree that relates a list of nodes.
}
\details{
Return a tree with tips corresponding to the nodes identified in
the input set that is consistent with the topology of the current
synthetic tree. This tree is equivalent to the minimal subtree
induced on the draft tree by the set of identified nodes.
}
\examples{
\dontrun{
## Result as a `phylo` object
res <- tol_induced_subtree(ott_ids = c(292466, 267845, 316878, 102710))

## Raw Newick string from Open Tree of Life
tree_file <- tempfile(fileext = ".tre")
tol_induced_subtree(ott_ids = c(292466, 267845, 316878, 102710),
                    file=tree_file)

}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/is_in_tree.R
\name{is_in_tree}
\alias{is_in_tree}
\title{Check that OTT ids occur in the Synthetic Tree}
\usage{
is_in_tree(ott_ids, ...)
}
\arguments{
\item{ott_ids}{a vector of Open Tree Taxonomy identifiers}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
A named logical vector. \code{TRUE} indicates that the OTT
    id is in the synthetic tree, and \code{FALSE} that it is not.
}
\description{
Some valid taxonomic names do not occur in the Synthetic
Tree. This convenience function allows you to check whether a
given Open Tree Taxonomy identifier (OTT id) is in the tree. A taxonomic
name may not occur in the synthetic tree because (1) it is an
extinct or invalid taxon, or (2) it is part of a group that is not
monophyletic in the tree.
}
\examples{
\dontrun{
  plant_families <- c("Asteraceae", "Solanaceae", "Poaceae", "Amaranthaceae",
                      "Zamiaceae", "Araceae", "Juncaceae")
  matched_names <- tnrs_match_names(plant_families)
  ## This fails because some ott ids are not in the tree
  ## plant_tree <- tol_induced_subtree(ott_id(matched_names))
  ## So let's check which ones are actually in the tree first:
  in_tree <- is_in_tree(ott_id(matched_names))
  ## This now works:
  plant_tree <- tol_induced_subtree(ott_id(matched_names)[in_tree])
}

}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tol.R
\name{tol_mrca}
\alias{tol_mrca}
\alias{tax_sources.tol_mrca}
\alias{unique_name.tol_mrca}
\alias{tax_name.tol_mrca}
\alias{tax_rank.tol_mrca}
\alias{ott_id.tol_mrca}
\alias{source_list.tol_mrca}
\title{MRCA of taxa from the synthetic tree}
\usage{
tol_mrca(ott_ids = NULL, node_ids = NULL, ...)

\method{tax_sources}{tol_mrca}(tax, ...)

\method{unique_name}{tol_mrca}(tax, ...)

\method{tax_name}{tol_mrca}(tax, ...)

\method{tax_rank}{tol_mrca}(tax, ...)

\method{ott_id}{tol_mrca}(tax, ...)

\method{source_list}{tol_mrca}(tax, ...)
}
\arguments{
\item{ott_ids}{Numeric vector. The ott ids for which the MRCA is desired.}

\item{node_ids}{Character vector. The node ids for which the MRCA is desired.}

\item{...}{additional arguments to customize the API call (see
\code{\link{rotl}} for more information).}

\item{tax}{an object returned by \code{tol_mrca()}.}
}
\value{
An invisible list of the MRCA node properties:

\itemize{

    \item {mrca} {List of node properties.}

    \itemize{
        \item {node_id} {String. The canonical identifier of the node.}

        \item {num_tips} {Numeric. The number of descendent tips.}

        \item {taxon} {A list of taxonomic properties. Only returned if
        the queried node is a taxon. (If the node is not a taxon, a
        \code{nearest_taxon} list is returned (see below)).}

            \itemize{
                \item {ott_id} {Numeric. The OpenTree Taxonomy ID (ottID).}

                \item {name} {String. The taxonomic name of the queried node.}

                \item {unique_name} {String. The string that uniquely
                identifies the taxon in OTT.}

                \item {rank} {String. The taxonomic rank of the taxon in OTT.}

               \item {tax_sources} {List. A list of identifiers for taxonomic
                sources, such as other taxonomies, that define taxa judged
                equivalent to this taxon.}
            }

        The following properties list support/conflict for the node across
        synthesis source trees. All properties involve sourceid keys and
        nodeid values (see \code{source_id_map} below) Not all properties are
        are present for every node.

        \item {partial_path_of} {List. The edge below this synthetic tree node
        is compatible with the edge below each of these input tree nodes (one
        per tree). Each returned element is reported as sourceid:nodeid.}

        \item {supported_by} {List. Input tree nodes (one per tree) that support
        this synthetic tree node. Each returned element is reported as
        sourceid:nodeid.}

        \item {terminal} {List. Input tree nodes (one per tree) that are equivalent
        to this synthetic tree node (via an exact mapping, or the input tree
        terminal may be the only terminal descended from this synthetic tree node.
        Each returned element is reported as sourceid:nodeid.}

        \item {conflicts_with} {Named list of lists. Names correspond to
        sourceid keys. Each list contains input tree node ids (one or more per tree)
        that conflict with this synthetic node.}
    }

    \item {nearest_taxon} {A list of taxonomic properties of the nearest rootward
    taxon node to the MRCA node. Only returned if the MRCA node is a not taxon
    (otherwise the \code{taxon} list above is returned).}

        \itemize{
            \item {ott_id} {Numeric. The OpenTree Taxonomy ID (ottID).}

            \item {name} {String. The taxonomic name of the queried node.}

            \item {unique_name} {String. The string that uniquely
            identifies the taxon in OTT.}

            \item {rank} {String. The taxonomic rank of the taxon in OTT.}

           \item {tax_sources} {List. A list of identifiers for taxonomic
            sources, such as other taxonomies, that define taxa judged
            equivalent to this taxon.}
        }

    \item {source_id_map} {Named list of lists. Names correspond to the
    sourceid keys used in the support/conflict properties of the \code{mrca}
    list above. Source trees will have the following properties:}

        \itemize{
            \item {git_sha} {The git SHA identifying a particular source
            version.}

            \item {tree_id} {The tree id associated with the study id used.}

            \item {study_id} {The study identifier. Will typically include
            a prefix ("pg_" or "ot_").}
        }
    The only sourceid that does not correspond to a source tree is the taxonomy,
    which will have the name "ott"+`taxonomy_version`, and the value is the
    ott_id of the taxon in that taxonomy version. "Taxonomy" will only ever
    appear in \code{supported_by}.

   }
}
\description{
Most Recent Common Ancestor for a set of nodes
}
\details{
Get the MRCA of a set of nodes on the current synthetic
    tree. Accepts any combination of node ids and ott ids as
    input. Returns information about the most recent common
    ancestor (MRCA) node as well as the most recent taxonomic
    ancestor (MRTA) node (the closest taxonomic node to the MRCA
    node in the synthetic tree; the MRCA and MRTA may be the same
    node). If they are the same, the taxonomic information will be
    in the \code{mrca} slot, otherwise they will be in the
    \code{nearest_taxon} slot of the list. If any of the specified
    nodes is not in the synthetic tree an error will be returned.

    Taxonomic methods (\code{tax_sources()}, \code{ott_id()},
    \code{unique_name()}, ...) are available on the objects
    returned by \code{tol_mrca()}. If the MRCA node is MRTA, the
    name of the object returned by these methods will start with
    \sQuote{ott}, otherwise it will start with \sQuote{mrca}.
}
\examples{
\dontrun{
birds_mrca <- tol_mrca(ott_ids=c(412129, 119214))
ott_id(birds_mrca)
tax_sources(birds_mrca)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tol.R
\name{tol_about}
\alias{tol_about}
\alias{tax_rank.tol_summary}
\alias{tax_sources.tol_summary}
\alias{unique_name.tol_summary}
\alias{tax_name.tol_summary}
\alias{ott_id.tol_summary}
\title{Information about the Tree of Life}
\usage{
tol_about(include_source_list = FALSE, ...)

\method{tax_rank}{tol_summary}(tax, ...)

\method{tax_sources}{tol_summary}(tax, ...)

\method{unique_name}{tol_summary}(tax, ...)

\method{tax_name}{tol_summary}(tax, ...)

\method{ott_id}{tol_summary}(tax, ...)
}
\arguments{
\item{include_source_list}{Logical (default =
\code{FALSE}). Return an ordered list of source trees.}

\item{...}{additional arguments to customize the API call (see
\code{\link{rotl}} for more information).}

\item{tax}{an object created with a call to \code{tol_about}.}
}
\value{
An invisible list of synthetic tree summary statistics:

\itemize{

    \item {date_created} {String. The creation date of the tree.}

    \item {num_source_studies} {Integer. The number of studies
    (publications)used as sources.}

    \item {num_source_trees} {The number of trees used as sources
    (may be >1 tree per study).}

    \item {taxonomy_version} {The Open Tree Taxonomy version used
    as a source.}

    \item {filtered_flags} {List. Taxa with these taxonomy flags were
    not used in construction of the tree.}

    \item {root} {List. Describes the root node:}
        \itemize{
            \item {node_id} {String. The canonical identifier of the node.}

            \item {num_tips} {Numeric. The number of descendent tips.}

            \item {taxon} {A list of taxonomic properties:}
            \itemize{
                \item {ott_id} {Numeric. The OpenTree Taxonomy ID (ott_id).}

                \item {name} {String. The taxonomic name of the queried node.}

                \item {unique_name} {String. The string that uniquely
                identifies the taxon in OTT.}

                \item {rank} {String. The taxonomic rank of the taxon in OTT.}

                \item {tax_sources} {List. A list of identifiers for taxonomic
                sources, such as other taxonomies, that define taxa judged
                equivalent to this taxon.}
            }
        }

    \item {source_list} {List. Present only if
    \code{include_source_list} is \code{TRUE}. The sourceid
    ordering is the precedence order for synthesis, with
    relationships from earlier trees in the list having priority
    over those from later trees in the list. See
    \code{source_id_map} below for study details.}

    \item {source_id_map} {Named list of lists. Present only if
    \code{include_source_list} is \code{TRUE}. Names correspond to
    the \sQuote{sourceids} used in \code{source_list}
    above. Source trees will have the following properties:}

        \itemize{
            \item {git_sha} {String. The git SHA identifying a particular source
            version.}

            \item {tree_id} {String. The tree id associated with the study id used.}

            \item {study_id} {String. The study identifier. Will typically include
            a prefix ("pg_" or "ot_").}
        }

    \item {synth_id} {The unique string for this version of the tree.}
}
}
\description{
Basic information about the Open Tree of Life (the synthetic tree)
}
\details{
Summary information about the current draft tree of life,
    including information about the list of trees and the taxonomy
    used to build it. The object returned by \code{tol_about} can
    be passed to the taxonomy methods (\code{tax_name()},
    \code{tax_rank()}, \code{tax_sources()}, \code{ott_id}), to
    extract relevant taxonomic information for the root of the
    synthetic tree.
}
\examples{
\dontrun{
res <- tol_about()
tax_sources(res)
ott_id(res)
studies <- source_list(tol_about(include_source_list=TRUE))}
}
\seealso{
\code{\link{source_list}} to explore the list of studies
    used in the synthetic tree (see example).
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/external_data.R
\name{taxon_external_IDs}
\alias{taxon_external_IDs}
\title{Get external identifiers for data associated with an Open Tree taxon}
\usage{
taxon_external_IDs(taxon_id)
}
\arguments{
\item{taxon_id}{An open tree study ID}
}
\value{
a data.frame in which each row represents a unique record in an
external database. The column "source" provides and abbreviated name for the
database, and "id" the unique ID for the record.
}
\description{
The Open Tree taxonomy is a synthesis of multiple reference taxonomies. This
function retrieves identifiers to external taxonomic records that have
contributed the rank, position and definition of a given Open Tree taxon.
}
\examples{
\dontrun{
   gibbon_IDs <- taxon_external_IDs(712902)
}
}
\seealso{
tnrs_matchnames, which can be used to search for taxa by name.

taxonomy_taxon, for more information about a given taxon.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tnrs.R
\name{tnrs_contexts}
\alias{tnrs_contexts}
\title{TNRS contexts}
\usage{
tnrs_contexts(...)
}
\arguments{
\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
Returns invisibly a list for each major clades (e.g.,
    animals, microbes, plants, fungi, life) whose elements
    contains the possible contexts.
}
\description{
This function returns a list of pre-defined taxonomic contexts
(i.e. clades) which can be used to limit the scope of tnrs
queries.
}
\details{
Taxonomic contexts are available to limit the scope of TNRS
searches. These contexts correspond to uncontested higher taxa
such as 'Animals' or 'Land plants'. This service returns a list
containing all available taxonomic context names, which may be
used as input (via the \code{context_name} argument in other
functions) to limit the search scope of other services including
\code{\link{tnrs_match_names}}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tol.R
\name{tol_subtree}
\alias{tol_subtree}
\title{Extract a subtree from the synthetic tree}
\usage{
tol_subtree(ott_id = NULL, node_id = NULL, label_format = NULL, file, ...)
}
\arguments{
\item{ott_id}{Numeric. The ott id of the node in the tree that should
serve as the root of the tree returned.}

\item{node_id}{Character. The node id of the node in the tree that should
serve as the root of the tree returned.}

\item{label_format}{Character. Defines the label type; one of
\dQuote{\code{name}}, \dQuote{\code{id}}, or
 \dQuote{\code{name_and_id}} (the default).}

\item{file}{If specified, the function will write the subtree to a
file in newick format.}

\item{...}{additional arguments to customize the API call (see
\code{\link{rotl}} for more information).}
}
\value{
If no value is specified to the \code{file} argument
    (default), a phylogenetic tree of class \code{phylo}.
    Otherwise, the function returns invisibly a logical indicating
    whether the file was successfully created.
}
\description{
Extract a subtree from the synthetic tree from an Open Tree node id.
}
\details{
Given a node, return the subtree of the synthetic tree descended
    from that node. The start node may be specified using either a node id
    or an ott id, but not both. If the specified node is not in the
    synthetic tree an error will be returned. There is a size limit of
    25000 tips for this method.
}
\examples{
\dontrun{
res <- tol_subtree(ott_id=241841)}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/match_names.R, R/methods.R
\name{ott_id.match_names}
\alias{ott_id.match_names}
\alias{flags.match_names}
\alias{flags}
\title{\code{ott_id} and \code{flags} for taxonomic names matched
    by \code{tnrs_match_names}}
\usage{
\method{ott_id}{match_names}(tax, row_number, taxon_name, ott_id, ...)

\method{flags}{match_names}(tax, row_number, taxon_name, ott_id, ...)

flags(tax, ...)
}
\arguments{
\item{tax}{an object returned by \code{\link{tnrs_match_names}}}

\item{row_number}{the row number corresponding to the name for
which to list the synonyms}

\item{taxon_name}{the taxon name corresponding to the name for
which to list the synonyms}

\item{ott_id}{the ott id corresponding to the name for which to
list the synonyms}

\item{...}{currently ignored}
}
\value{
A list of the ott ids or flags for the taxonomic names
    matched with \code{\link{tnrs_match_names}}, for either one or
    all the names.
}
\description{
\code{rotl} provides a collection of functions that allows users
to extract relevant information from an object generated by
\code{\link{tnrs_match_names}} function.
}
\details{
These methods optionally accept one of the arguments
\code{row_number}, \code{taxon_name} or \code{ott_id} to retrieve
the corresponding information for one of the matches in the object
returned by the \code{\link{tnrs_match_names}} function.

If these arguments are not provided, these methods can return
information for the matches currently listed in the object
returned by \code{\link{tnrs_match_names}}.
}
\examples{
\dontrun{
  rsp <- tnrs_match_names(c("Diadema", "Tyrannosaurus"))
  rsp$ott_id    # ott id for match currently in use
  ott_id(rsp)   # similar as above but elements are named

  ## flags() is useful for instance to determine if a taxon is extinct
  flags(rsp, taxon_name="Tyrannosaurus")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studies.R
\name{studies_properties}
\alias{studies_properties}
\title{Properties of the Studies}
\usage{
studies_properties(...)
}
\arguments{
\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
A list of the study properties that can be used to find
    studies and trees that are contributing to the synthetic tree.
}
\description{
Return the list of study properties that can be used to search
studies and trees used in the synthetic tree.
}
\details{
The list returned has 2 elements \code{tree_properties} and
\code{studies_properties}. Each of these elements lists additional
arguments to customize the API request properties that can be used
to search for trees and studies that are contributing to the
synthetic tree. The definitions of these properties are available
from
\url{https://github.com/OpenTreeOfLife/phylesystem-api/wiki/NexSON}
}
\examples{
\dontrun{
 all_the_properties <- studies_properties()
 unlist(all_the_properties$tree_properties)
}
}
\seealso{
\code{\link{studies_find_trees}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studies.R
\name{get_study_tree}
\alias{get_study_tree}
\title{Study Tree}
\usage{
get_study_tree(
  study_id = NULL,
  tree_id = NULL,
  object_format = c("phylo"),
  tip_label = c("original_label", "ott_id", "ott_taxon_name"),
  file_format,
  file,
  deduplicate = TRUE,
  ...
)
}
\arguments{
\item{study_id}{the identifier of a study (character)}

\item{tree_id}{the identifier of a tree within the study}

\item{object_format}{the class of the object to be returned
(default and currently only possible value \code{phylo} from
the \code{ape} package).}

\item{tip_label}{the format of the tip
labels. \dQuote{\code{original_label}} (default) returns the
original labels as provided in the study,
\dQuote{\code{ott_id}} labels are replaced by their ott IDs,
\dQuote{\code{ott_taxon_name}} labels are replaced by their
Open Tree Taxonomy taxon name.}

\item{file_format}{the format of the file to be generated
(\code{newick} default, \code{nexus}, or \code{json}).}

\item{file}{the file name where the output of the function will be
saved.}

\item{deduplicate}{logical (default \code{TRUE}). If the tree
returned by the study contains duplicated taxon names, should they
be made unique? It is normally illegal for NEXUS/Newick tree
strings to contain duplicated tip names. This is a workaround to
circumvent this requirement. If \code{TRUE}, duplicated tip labels
will be appended \code{_1}, \code{_2}, etc.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
if \code{file_format} is missing, an object of class
    \code{phylo}, otherwise a logical indicating whether the file
    was successfully created.
}
\description{
Returns a specific tree from within a study
}
\examples{
\dontrun{
 tree <- get_study_tree(study_id="pg_1144", tree_id="tree2324")

 ## comparison of the first few tip labels depending on the options used
 head(get_study_tree(study_id="pg_1144", tree_id="tree2324", tip_label="original_label")$tip.label)
 head(get_study_tree(study_id="pg_1144", tree_id="tree2324", tip_label="ott_id")$tip.label)
 head(get_study_tree(study_id="pg_1144", tree_id="tree2324", tip_label="ott_taxon_name")$tip.label)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tnrs.R
\name{tnrs_infer_context}
\alias{tnrs_infer_context}
\title{Infer the taxonomic context from a list of names}
\usage{
tnrs_infer_context(names = NULL, ...)
}
\arguments{
\item{names}{Vector of taxon names.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
A list including the context name, the context ott id and
    possibly the names in the query that have an ambiguous
    taxonomic meaning in the query.
}
\description{
Return a taxonomic context given a list of taxonomic names
}
\details{
Find the least inclusive taxonomic context that includes all the
unambiguous names in the input set. Unambiguous names are names
with exact matches to non-homonym taxa. Ambiguous names (those
without exact matches to non-homonym taxa) are indicated in
results.
}
\examples{
\dontrun{
res <- tnrs_infer_context(names=c("Stellula calliope", "Struthio camelus"))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studies.R
\name{studies_find_studies}
\alias{studies_find_studies}
\title{Find a Study}
\usage{
studies_find_studies(
  property = NULL,
  value = NULL,
  verbose = FALSE,
  exact = FALSE,
  detailed = TRUE,
  ...
)
}
\arguments{
\item{property}{The property to be searched on (character)}

\item{value}{The property value to be searched on (character)}

\item{verbose}{Should the output include all metadata (logical
default \code{FALSE})}

\item{exact}{Should exact matching be used? (logical, default
\code{FALSE})}

\item{detailed}{If \code{TRUE} (default), the function will return
a data frame that summarizes information about the study (see
\sQuote{Value}). Otherwise, it only returns the study
identifiers.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
If \code{detailed=TRUE}, the function returns a data frame
    listing the study id (\code{study_ids}), the number of trees
    associated with this study (\code{n_trees}), the tree ids (at
    most 5) associated with the studies (\code{tree_ids}), the
    tree id that is a candidate for the synthetic tree if any
    (\code{candidate}), the year of publication of the study
    (\code{study_year}), the title of the publication for the
    study (\code{title}), and the DOI (Digital Object Identifier)
    for the study (\code{study_doi}).

    If \code{detailed=FALSE}, the function returns a data frame
    with a single column containing the study identifiers.
}
\description{
Return the identifiers of studies that match given properties
}
\examples{
\dontrun{
## To match a study for which the identifier is already known
one_study <- studies_find_studies(property="ot:studyId", value="pg_719")
list_trees(one_study)

## To find studies pertaining to Mammals
mammals <- studies_find_studies(property="ot:focalCladeOTTTaxonName",
                                value="mammalia")
## To extract the tree identifiers for each of the studies
list_trees(mammals)
## ... or for a given study
list_trees(mammals, "ot_308")

## Just the identifiers without other information about the studies
mammals <- studies_find_studies(property="ot:focalCladeOTTTaxonName",
                                value="mammalia", detailed=FALSE)
}
}
\seealso{
\code{\link{studies_properties}} which lists properties
    against which the studies can be
    searched. \code{\link{list_trees}} that returns a list for all
    tree ids associated with a study.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studies.R
\name{get_study_subtree}
\alias{get_study_subtree}
\title{Study Subtree}
\usage{
get_study_subtree(
  study_id,
  tree_id,
  subtree_id,
  object_format = c("phylo"),
  tip_label = c("original_label", "ott_id", "ott_taxon_name"),
  file_format,
  file,
  deduplicate = TRUE,
  ...
)
}
\arguments{
\item{study_id}{the study identifier (character)}

\item{tree_id}{the tree identifier (character)}

\item{subtree_id, }{either a node id that specifies a subtree or
\dQuote{ingroup} which returns the ingroup for this subtree.}

\item{object_format}{the class of the object returned by the
function (default, and currently only possibility \code{phylo}
from the \code{ape} package)}

\item{tip_label}{the format of the tip
labels. \dQuote{\code{original_label}} (default) returns the
original labels as provided in the study,
\dQuote{\code{ott_id}} labels are replaced by their ott IDs,
\dQuote{\code{ott_taxon_name}} labels are replaced by their
Open Tree Taxonomy taxon name.}

\item{file_format}{character, the file format to use to save the
results of the query (possible values, \sQuote{newick} or
\sQuote{nexus}).}

\item{file}{character, the path and file name where the output
should be written.}

\item{deduplicate}{logical (default \code{TRUE}). If the tree
returned by the study contains duplicated taxon names, should
they be made unique? It is normally illegal for NEXUS/Newick
tree strings to contain duplicated tip names. This is a
workaround to circumvent this requirement. If \code{TRUE},
duplicated tip labels will be appended \code{_1}, \code{_2},
etc.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\description{
Retrieve subtree from a specific tree in the Open Tree of Life data store
}
\examples{
\dontrun{
small_tr <- get_study_subtree(study_id="pg_1144", tree_id="tree5800", subtree_id="node991044")
ingroup  <- get_study_subtree(study_id="pg_1144", tree_id="tree5800", subtree_id="ingroup")
nexus_file <- tempfile(fileext=".nex")
get_study_subtree(study_id="pg_1144", tree_id="tree5800", subtree_id="ingroup", file=nexus_file,
                  file_format="nexus")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxonomy_taxon_info}
\alias{taxonomy_taxon_info}
\alias{tax_rank.taxon_info}
\alias{tax_name.taxon_info}
\alias{unique_name.taxon_info}
\alias{synonyms.taxon_info}
\alias{ott_id.taxon_info}
\alias{tax_sources.taxon_info}
\alias{is_suppressed.taxon_info}
\alias{flags.taxon_info}
\title{Taxon information}
\usage{
taxonomy_taxon_info(
  ott_ids,
  include_children = FALSE,
  include_lineage = FALSE,
  include_terminal_descendants = FALSE,
  ...
)

\method{tax_rank}{taxon_info}(tax, ...)

\method{tax_name}{taxon_info}(tax, ...)

\method{unique_name}{taxon_info}(tax, ...)

\method{synonyms}{taxon_info}(tax, ...)

\method{ott_id}{taxon_info}(tax, ...)

\method{tax_sources}{taxon_info}(tax, ...)

\method{is_suppressed}{taxon_info}(tax, ...)

\method{flags}{taxon_info}(tax, ...)
}
\arguments{
\item{ott_ids}{the ott ids of the taxon of interest (numeric or
character containing only numbers)}

\item{include_children}{whether to include information about all
the children of this taxon. Default \code{FALSE}.}

\item{include_lineage}{whether to include information about all
the higher level taxa that include the \code{ott_ids}.
Default \code{FALSE}.}

\item{include_terminal_descendants}{whether to include the list of
terminal \code{ott_ids} contained in the \code{ott_ids}
provided.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}

\item{tax}{an object generated by the \code{taxonomy_taxon_info}
function}
}
\value{
\code{taxonomy_taxon_info} returns a list detailing
    information about the taxa. \code{tax_rank} and
    \code{tax_name} return a vector. \code{synonyms} returns a
    list whose elements are the synonyms for each of the
    \code{ott_id} requested.
}
\description{
Information about taxa.
}
\details{
Given a vector of ott ids, \code{taxonomy_taxon_info} returns
information about the specified taxa.

The functions \code{tax_rank}, \code{tax_name}, and
\code{synonyms} can extract this information from an object
created by the \code{taxonomy_taxon_info()}.
}
\examples{
\dontrun{
req <- taxonomy_taxon_info(ott_id=515698)
tax_rank(req)
tax_name(req)
synonyms(req)
}
}
\seealso{
\code{\link{tnrs_match_names}} to obtain \code{ott_id}
    from a taxonomic name.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxonomy_about}
\alias{taxonomy_about}
\title{Information about the Open Tree Taxonomy}
\usage{
taxonomy_about(...)
}
\arguments{
\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
A list with the following properties:
\itemize{

    \item {weburl} {String. The release page for this version
    of the taxonomy.}

    \item {author} {String. The author string.}

    \item {name} {String. The name of the taxonomy.}

    \item {source} {String. The full identifying information for
    this version of the taxonomy.}

    \item {version} {String. The version number of the taxonomy.}
}
}
\description{
Summary information about the Open Tree Taxonomy (OTT)
}
\details{
Return metadata and information about the taxonomy
itself. Currently, the available metadata is fairly sparse, but
includes (at least) the version, and the location from which the
complete taxonomy source files can be downloaded.
}
\examples{
\dontrun{
taxonomy_about()
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxonomy.R
\name{taxonomy_mrca}
\alias{taxonomy_mrca}
\alias{tax_rank.taxon_mrca}
\alias{tax_name.taxon_mrca}
\alias{ott_id.taxon_mrca}
\alias{unique_name.taxon_mrca}
\alias{tax_sources.taxon_mrca}
\alias{flags.taxon_mrca}
\alias{is_suppressed.taxon_mrca}
\title{Taxonomic MRCA}
\usage{
taxonomy_mrca(ott_ids = NULL, ...)

\method{tax_rank}{taxon_mrca}(tax, ...)

\method{tax_name}{taxon_mrca}(tax, ...)

\method{ott_id}{taxon_mrca}(tax, ...)

\method{unique_name}{taxon_mrca}(tax, ...)

\method{tax_sources}{taxon_mrca}(tax, ...)

\method{flags}{taxon_mrca}(tax, ...)

\method{is_suppressed}{taxon_mrca}(tax, ...)
}
\arguments{
\item{ott_ids}{a vector of ott ids for the taxa whose MRCA is to
be found (numeric).}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}

\item{tax}{an object generated by the \code{taxonomy_mrca}
function}
}
\value{
\itemize{

    \item{\code{taxonomy_mrca}} { returns a list about the
    taxonomic information relating to the MRCA for the ott_ids
    provided. }

    \item{\code{tax_rank}} { returns a character vector of the
    taxonomic rank for the MRCA. }

    \item{\code{tax_name}} { returns a character vector the
    Open Tree Taxonomy name for the MRCA. }

    \item{\code{ott_id}} { returns a numeric vector of the ott id
    for the MRCA. }

}
}
\description{
Taxonomic Least Inclusive Common Ancestor (MRCA)
}
\details{
Given a set of OTT ids, get the taxon that is the most recent common
ancestor (the MRCA) of all the identified taxa.
}
\examples{
\dontrun{
req <- taxonomy_mrca(ott_ids=c(515698,590452,643717))
tax_rank(req)
tax_name(req)
ott_id(req)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/studies.R
\name{get_study}
\alias{get_study}
\title{Get all the trees associated with a particular study}
\usage{
get_study(
  study_id = NULL,
  object_format = c("phylo", "nexml"),
  file_format,
  file,
  ...
)
}
\arguments{
\item{study_id}{the study ID for the study of interest (character)}

\item{object_format}{the class of the object the query should
return (either \code{phylo} or \code{nexml}). Ignored if
\code{file_format} is specified.}

\item{file_format}{the format of the file to be generated
(\code{newick}, \code{nexus}, \code{nexml} or \code{json}).}

\item{file}{the file name where the output of the function will be
saved.}

\item{...}{additional arguments to customize the API request (see
\code{\link{rotl}} package documentation).}
}
\value{
if \code{file_format} is missing, an object of class
    \code{phylo} or \code{nexml}, otherwise a logical indicating
    whether the file was successfully created.
}
\description{
Returns the trees associated with a given study
}
\details{
If \code{file_format} is missing, the function returns an object
of the class \code{phylo} from the \code{ape} package
(default), or an object of the class \code{nexml} from the
\code{RNeXML} package.

Otherwise \code{file_format} can be either \code{newick},
\code{nexus}, \code{nexml} or \code{json}, and the function will
generate a file of the selected format. In this case, a file name
needs to be provided using the argument \code{file}. If a file
with the same name already exists, it will be silently
overwritten.
}
\examples{
\dontrun{
that_one_study <- get_study(study_id="pg_719", object_format="phylo")
if (require(RNeXML)) { ## if RNeXML is installed get the object directly
   nexml_study <- get_study(study_id="pg_719", object_format="nexml")
} else { ## otherwise write it to a file
   get_study(study_id="pg_719", file_format="nexml", file=tempfile(fileext=".nexml"))
}
}
}
\seealso{
\code{\link{get_study_meta}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R
\name{tax_rank}
\alias{tax_rank}
\alias{ott_id}
\alias{synonyms}
\alias{tax_sources}
\alias{is_suppressed}
\alias{unique_name}
\alias{tax_name}
\title{Methods for Taxonomy}
\usage{
tax_rank(tax, ...)

ott_id(tax, ...)

synonyms(tax, ...)

tax_sources(tax, ...)

is_suppressed(tax, ...)

unique_name(tax, ...)

tax_name(tax, ...)
}
\arguments{
\item{tax}{an object returned by \code{\link{taxonomy_taxon_info}},
\code{\link{taxonomy_mrca}}, or \code{\link{tnrs_match_names}}}

\item{...}{additional arguments (see
\code{\link{tnrs_match_names}})}
}
\description{
Methods for dealing with objects containing taxonomic information
(Taxonomy, TNRS endpoints)
}
\details{
This is the page for the generic methods. See the help pages for
\code{\link{taxonomy_taxon_info}}, \code{\link{taxonomy_mrca}}, and
\code{\link{tnrs_match_names}} for more information.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/methods.R, R/taxonomy.R
\name{tax_lineage}
\alias{tax_lineage}
\alias{tax_lineage.taxon_info}
\title{Lineage of a taxon}
\usage{
tax_lineage(tax, ...)

\method{tax_lineage}{taxon_info}(tax, ...)
}
\arguments{
\item{tax}{an object created by \code{\link{taxonomy_taxon_info}}
using the argument \code{include_lineage=TRUE}.}

\item{...}{additional arguments (currently unused).}
}
\value{
A list with one slot per taxon that contains a data frame
    with 3 columns: the taxonomy rank, the name, and unique name
    for all taxa included in the lineage of the taxon up to the
    root of the tree.
}
\description{
Extract the lineage information (higher taxonomy) from an object
returned by \code{\link{taxonomy_taxon_info}}.
}
\details{
The object passed to this function must have been created using
the argument \code{include_lineage=TRUE}.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rotl-package.R
\docType{package}
\name{rotl}
\alias{rotl}
\title{An Interface to the Open Tree of Life API}
\description{
The Open Tree of Life is an NSF funded project that is generating
an online, comprehensive phylogenetic tree for 1.8 million
species. \code{rotl} provides an interface that allows you to
query and retrieve the parts of the tree of life that is of
interest to you.
}
\details{
\code{rotl} provides function to most of the end points the API
provides. The documentation of the API is available at:
\url{https://github.com/OpenTreeOfLife/opentree/wiki/Open-Tree-of-Life-APIs}
}
\section{Customizing API calls}{


    All functions that use API end points can take 2 arguments to
    customize the API call and are passed as \code{...} arguments.

    \itemize{

    \item{ \code{otl_v} } { This argument controls which version
    of the API your call is using. The default value for this
    argument is a call to the non-exported function
    \code{otl_version()} which returns the current version of the
    Open Tree of Life APIs (v2).}

    \item{ \code{dev_url} } { This argument controls whether to use
    the development version of the API. By default, \code{dev_url}
    is set to \code{FALSE}, using \code{dev_url = TRUE} in your
    function calls will use the development version.}

    }

    For example, to use the development version of the API, you
    could use: \code{tnrs_match_names("anas", dev_url=TRUE)}

    Additional arguments can also be passed to the
    \code{\link[httr]{GET}} and \code{\link[httr]{POST}} methods.
}

\section{Acknowledgments}{


    This package was started during the Open Tree of Life
    \href{https://blog.opentreeoflife.org/2014/06/11/apply-for-tree-for-all-a-hackathon-to-access-opentree-resources/}{Hackathon}
    organized by OpenTree, the NESCent Hackathon Interoperability
    Phylogenetic group, and Arbor.
}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tol.R
\name{strip_ott_ids}
\alias{strip_ott_ids}
\title{Strip OTT ids from tip labels}
\usage{
strip_ott_ids(tip_labels, remove_underscores = FALSE)
}
\arguments{
\item{tip_labels}{a character vector containing tip labels (most
likely the \code{tip.label} element from a tree returned by
\code{\link{tol_induced_subtree}}}

\item{remove_underscores}{logical (defaults to FALSE). If set to
TRUE underscores in tip labels are converted to spaces}
}
\value{
A character vector containing the contents of
    \code{tip_labels} with any OTT ids removed.
}
\description{
Strip OTT ids from tip labels
}
\examples{
\dontrun{
genera <- c("Perdix", "Setophaga", "Cinclus", "Struthio")
tr <- tol_induced_subtree(ott_ids=c(102710, 285198, 267845, 292466))
tr$tip.label \%in\% genera
tr$tip.label <- strip_ott_ids(tr$tip.label)
tr$tip.label \%in\% genera
}
}
