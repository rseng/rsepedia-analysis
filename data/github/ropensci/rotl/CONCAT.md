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
