
<!-- README.md is generated from README.Rmd. Please edit that file -->

# RNeXML <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![DOI](https://zenodo.org/badge/11856817.svg)](https://zenodo.org/badge/latestdoi/11856817)
[![Build
Status](https://api.travis-ci.org/ropensci/RNeXML.png)](https://travis-ci.org/ropensci/RNeXML)
[![Windows build
status](https://ci.appveyor.com/api/projects/status/dhiwp5blx2ns2yba/branch/master?svg=true)](https://ci.appveyor.com/project/cboettig/rnexml/branch/master)
[![CRAN
status](https://www.r-pkg.org/badges/version/RNeXML)](https://cran.r-project.org/package=RNeXML)
[![codecov.io](https://codecov.io/github/ropensci/RNeXML/coverage.svg?branch=master)](https://codecov.io/github/ropensci/RNeXML?branch=master)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/RNeXML)

  - Maintainer: Carl Boettiger
  - Authors: Carl Boettiger, Scott Chamberlain, Hilmar Lapp, Kseniia
    Shumelchyk, Rutger Vos
  - License: BSD-3
  - [Issues](https://github.com/ropensci/RNeXML/issues): Bug reports,
    feature requests, and development discussion.

An extensive and rapidly growing collection of richly annotated
phylogenetics data is now available in the NeXML format. NeXML relies on
state-of-the-art data exchange technology to provide a format that can
be both validated and extended, providing a data quality assurance and
adaptability to the future that is lacking in other formats. See [Vos et
al
2012](http://doi.org/10.1093/sysbio/sys025 "NeXML: Rich, Extensible, and Verifiable Representation of Comparative Data and Metadata.")
for further details on the NeXML format.

## How to cite

RNeXML has been published in the following article:

> Boettiger C, Chamberlain S, Vos R and Lapp H (2016). “RNeXML: A
> Package for Reading and Writing Richly Annotated Phylogenetic,
> Character, and Trait Data in R.” *Methods in Ecology and Evolution*,
> **7**, pp. 352-357.
> [doi:10.1111/2041-210X.12469](http://doi.org/10.1111/2041-210X.12469)

Although the published version of the article is pay-walled, the source
of the manuscript, and a much better rendered PDF, are included in this
package (in the `manuscripts` folder). You can also find it [freely
available on arXiv](http://arxiv.org/abs/1506.02722).

## Installation

The latest stable release of RNeXML is on CRAN, and can be installed
with the usual `install.packages("RNeXML")` command. Some of the more
specialized functionality described in the Vignettes (such as RDF
manipulation) requires additional packages which can be installed using:

``` r
install.packages("RNeXML", deps = TRUE)
```

The development version can be installed using:

``` r
remotes::install_github("ropensci/RNeXML")
```

## Getting Started

See the vignettes below for both a general quick start and an overview
of more specialized features.

  - [A Brief Introduction to
    RNeXML](https://docs.ropensci.org/RNeXML/articles/intro)
  - [RNeXML: A Package for Reading and Writing Richly Annotated
    Phylogenetic, Character, and Trait Data in
    R](https://github.com/ropensci/RNeXML/tree/master/manuscripts)
    (published in MEE).
  - [Handling Metadata in
    RNeXML](https://docs.ropensci.org/RNeXML/articles/metadata)
  - [The `nexml` S4
    Object](https://docs.ropensci.org/RNeXML/articles/S4)
  - [Semantic data & SPARQL with
    RNeXML](https://docs.ropensci.org/RNeXML/articles/sparql)
  - [Extending NeXML: an example based on
    simmap](https://docs.ropensci.org/RNeXML/articles/simmap)

-----

[![ropensci
footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)
NEWS
====

For more fine-grained list of changes or to report a bug, consult 

* [The issues log](https://github.com/ropensci/RNeXML/issues)
* [The commit log](https://github.com/ropensci/RNeXML/commits/master)

v2.4.4
------

Compatibility with upcoming `dplyr` release.

v2.4.3
------

- This update fixes a minor bug in a unit test for compatibility 
  with R-devel (4.0.0) as requested by CRAN (#245)

v2.4.2
------

- This update fixes a minor bug in a unit test which was introduced by a recent change to the geiger package.

v2.4.0
------

- Makes various tests more robust, and uses symbolic address for nexml.org (#238)
- Provides a real `summary()` and improved pretty-print (#237)
- Makes `character(0)` metadata value behave as empty string (#236)
- Fixes detection of having to split matrix by class (#235)
- Switch over to Additional_repositories for CRAN (#229)
- Do not add `ter` namespace by default. (#227)
- Replace taxize with taxald (#226)
- Fixes how metadata arguments are passed on to `add_basic_meta()` (#220)
- Fixes CDAO namespace definition [#219]
- Enables handling of files with `rootedge` [218]


v2.3.0
-------


This release addresses several aspects improving the handling of metadata:

- `add_meta()` now works for trees and characters (#213, PR #217)
- Handles nested meta elements properly (#196, PR #197)

Misc fixes:

- enable handling of `rootEdge` (#207, PR #218)
- Replaces taxize backend with faster alternative `taxadb` method. (#224, PR #226).
 This remains only a suggested package and has much lighter dependencies as well.  
- add hex

v2.2.0
------

- Fixes various (previously broken) aspects of handling polymorphic
  and uncertain states for discrete (non-molecular) and continuous
  characters, including obtaining a character matrix (#174), ensuring
  proper column types (#188), and serializing to NeXML (#192).
- Adds the optional ability to, in addition to the character matrix,
  obtain a concordant formatted matrix of state types (standard,
  polymorphic, uncertain).
- Fixes loss of certain literal-valued metadata when serializing to
  NeXML. #193
- Drops package phylobase as dependency. (Also removes circular
  dependency chain, because phylobase depends on RNeXML.)

v2.1.2
------

- Fix failing checks on CRAN that require a network connection

v2.1.1 
------

- avoid rdf-based tests on solaris architecture, where suggested
  package rdflib is not available. (CRAN request.)

v2.1.0 2018-05-05
------

- `taxize` as Suggests only
- drop `rrdf` in favor of `rdflib`
- drop `Sxslt` in favor of `xslt`


v2.0.8 2017-11-17
------

- patch for compatibility with upcoming release of `testthat`

v2.0.7 2016-06-28
------

- Bugfixes following release of new dplyr and new tidyr dependencies

v2.0.6 2016-03-07
------

- Migrate Additional_repositories to new address for OmegaHat project.

v2.0.5  2015-12-31
-------

- `get_metadata()`, `get_taxa()` now return much richer `data.frames` instead of named vectors. 
  This is potentially a non-backwards compatible change if scripts use the output of these
  functions as lists (#129).  See updated metadata vignette.  This introduces new dependencies
  `dplyr` and `lazyeval`. 
- more robust `nexml_read()` method for URLs, (#123)
- Avoid assuming the namespace prefix `nex` for nexml elements (#51, #124, #126). Includes a
  fix server-side on the NeXML validator as well.
- `nexml_validate()` points to the new validator. (#126)


v2.0.4 2015-10-14
-------

- Fix compatibility issue with recent phytools release.

v2.0.3 2015-05-27
------

- Upgrade tests to be compatible with newest testthat (0.10.0), bumps testthat dependency version up (#119) thanks @hadley

v2.0.2 2015-05-01
------

- Add four new vignettes describing the use of various advanced
  features in the package: the use of SPARQL queries, advanced
	use of metadata features, an example of how to extend NeXML
	with simmap data as the use case, and documentation on the 
	central S4 data structure used in the package.
- Implements the use of Title Case in the package title, as
  requested (on several occasions) by the CRAN maintainers.


v2.0.1 2014-12-26
-------

- Update DESCRIPTION to provide a standard `install.packages()` compatible repository for `rrdf`, as per request from the CRAN team.

v2.0.0   2014-12-06
---------

* add URL and BugReports to Description. [#103](https://github.com/ropensci/RNeXML/issues/103)

* for consistency with other `add_` methods, the `nexml` object is now the _last_, not the _first_, 
argument to `add_basic_meta`.  As this changes the function API, it could break code that does not
explicitly name the arguments, so we release this as 2.0.0


v1.1.3 2014-08-06
------

Minor bugfix

* Fixes typo that caused validator to fail when nexml.org couldn't be reached

v1.1.2  2014-07-19
-------

Less aggressive unit-tests

* nexml_validate now returns NULL if the validation cannot be performed. Unit tests now consider either TRUE or NULL as acceptable.   
* Just skips the uuid unit test if uuid package is not available
* Documented versioning practice in NEWS
* Unit tests relying on the Figshare API are not run (without failing) if authentication to figshare server fails
* Documentation updated to include examples for all functions

v1.1-0 2014-07-18
------

Initial Release 

Contributing Guidelines
=======================


Repository structure
--------------------

This repository is structured as a standard R package
following the conventions outlined in the [Writing R
extensions](http://cran.r-project.org/doc/manuals/R-exts.html) manual.
A few additional files are provided that are not part of the built
R package and are listed in `.Rbuildignore`, such as `.travis.yml`,
which is used for continuous testing and integration.


Code
----

All code for this package is found in `R/`, (except compiled source
code, if used, which is in `/src`).  All functions should be thoroughly
documented with `roxygen2` notation; see Documentation.

Testing
-------

Any new feature or bug-fix should include a unit-test demonstrating the
change.  Unit tests follow the `testthat` framework with files in
`tests/testthat`.  Please make sure that the testing suite passes
before issuing a pull request.  This can be done by running `check()`
from the `devtools` package, which will also check for consistent
documentation, etc.


This package uses the [travis](https://github.com/craigcitro/r-travis)
continuous testing mechanism for R to ensure that the test suite is run
on each push to Github.  An icon at the top of the README.md indicates
whether or not the tests are currently passing.


Documentation
-------------

All of the function documentation is generated automatically.
Please do not edit any of the documentation files in `man/`
or the `NAMESPACE`.  Instead, construct the appropriate
[roxygen2](https://github.com/klutometis/roxygen) documentation in the
function files in `R/` themselves.  The documentation is then generated
by running the `document()` function from the `devtools` package.  Please
consult the [Advanced R programming](http://adv-r.had.co.nz/) guide if
this workflow is unfamiliar to you.  Note that functions should include
examples in the documentation. Please use `\dontrun` for examples that
take more than a few seconds to execute or require an internet connection.

Likewise, the README.md file in the base directory should not be edited
directly.  This file is created automatically from code that runs the
examples shown, helping to ensure that they are functioning as advertised
and consistent with the package README vignette.  Instead, edit the
`README.Rmd` source file in `manuscripts` and run `make` to build
the README.



Manuscript
----------

The manuscript files are built using the dynamic documentation tool
`knitr` from the `.Rmd` versions of the file.  Please do not edit
the `.md` versions since such files are built automatically.
The `.md` versions are built for viewing on Github, and take advantage
of Github's rendering and display of text-based diffs.  The `.md`
files should then be generated using the `Makefile` provided, which
will also handle details such as Github-compatible syntax highlighting
and the embedding of images.

The `Makefile` will also generate a pdf version of the manuscript
using pandoc and the appropriate LaTeX templates.

Text should be hard-wrapped at less than 80 characters width when
possible. This allows git to better track real changes to the files
and impoves the display of line-based diffs.  For this reason,
also avoid re-wrapping text frequently, or changing line end encodings,
etc.


**Embedding images**: Image generation is handled by the markdown
file, which will embed online png images published to imgur for the
`.md` output, and vector pdf graphics for the `.pdf` manuscript.

**Citations**: Citations should be added to the `.Rmd` file using
pandoc markdown notation, with the corresponding bibtex entries
added to `citations.bib`. Citations can also be added as a standard
markdown link.

**Caching** To avoid rerunning potentially slow R code embedded in the
mansucript simply to view changes to the text, results from running the
code are cached in `cache` (see the knitr documentation for details on
how caching is used).  Run `make clean` to erase the cache and clear
your workspace.  Recall that the Makefile will only rerun the relevant
command if the source file has changed.  Consequently, changes to to
the package functions themselves will not automatically cause make to
recompile the manuscript.


<!-- Should add a utility that will generate citation metadata from
the mauscript.Rmd links using knitcitations.

Consider yaml-based citation format instead, see:
http://blog.martinfenner.org/2013/07/30/citeproc-yaml-for-bibliographies/#comment-1046228784
-->


Branches
--------

Please ensure that any pull requests are made to the relevant branch.


Questions or comments?
---------------------

Do not hesistate to open an issue in the issues tracker to raise any
questions or comments about the package or these guidelines.



Dear CRAN maintainers,

This update accomodates upcoming changes in dplyr.

Carl



Instructions for compiling manuscripts
======================================


[![Circle CI](https://circleci.com/gh/ropensci/RNeXML.svg?style=svg)](https://circleci.com/gh/ropensci/RNeXML)


Install dependencies
--------------------

Install the dependencies required for the supplementary examples using the `devtools` package:

```r
install.packages("devtools")
devtools::install_github(c("egonw/rrdf/rrdflibs", "egonw/rrdf/rrdf", "cboettig/Sxslt"))
```

Then install the `RNeXML` R package, including the suggested packages, using the following R command:

```r
install.packages("RNeXML", dependencies=TRUE)
```

Note that `rmarkdown` requires `pandoc` (>= 0.12.3) and `pandoc-citeproc` be installed. These ship with the current version of RStudio (`>=0.98`). Additionally, a LaTeX environment is required to generate the output pdf. 



Build the manuscript
--------------------


Make sure you set the `manuscripts/` as your working directory and then do:

```r
rmarkdown::render("manuscript.Rmd")
```
or use the `knit2pdf` button in your RStudio console. 

Alternately: Using Docker
-------------------------

Instead of installing R packages seperately, you can try out RNeXML
by running RStudio in a Docker container.  This (a) avoids having to install
software dependencies, and (b) avoids altering anything on your local
library. If the above doesn't work, or just for fun, give this a try.

The `RNeXML` package and all dependencies are installed on the [rocker/ropensci](http://registry.hub.docker.com/u/rocker/ropensci) Docker container.  You will still need all
the files from this directory (the `manuscripts` directory on the RNeXML Github repository)
to build the manuscript. Users can decide to run either an R console (accessed through a terminal) 
or an RStudio instance (accessed through the browser) on the container. 



### Docker Installation

In a Mac or Windows machine, this will aslo install boot2docker
(easy point & click install, ~24 MB). On Linux, this installs
natively and we can run everything in the terminal instead of in
the boot2docker window.
([Mac](https://docs.docker.com/installation/mac/),
[Windows](https://docs.docker.com/installation/windows/),
[Linux](https://docs.docker.com/installation)).

### R console

With boot2docker running, run `R` on the `rocker/ropensci` image,
linking the location of your copy of this directory to
`/home/rstudio` on the container, setting the container's
working directory to the same, and setting user as `rstudio`:

```bash
docker run -v /path/to/RNeXML/manuscripts:/home/rstudio \ 
  -w /home/rstudio -u rstudio -ti --rm rocker/ropensci R
```

At the R prompt, you can use `rmarkdown` to render the manuscript PDF from the `Rmd` file:

```r
rmarkdown::render('manuscript.Rmd')
```

`manuscript.pdf` should now be created in the manuscripts directory.  

### Using RStudio

1) From the command line (with boot2docker running on Mac/Windows), do:

```bash
sudo docker run -d -p 8787:8787 rocker/ropensci
```

That will take a while to download the image the first time you run it.

2) Once it is done, try:

```bash
boot2docker ip
```
that should return an ip address you can paste into your browser.

3) Add a `:8787` to the end of this address and paste it into your
browser address bar. (e.g. it's probably `http://92.168.59.103:8787`
but that can change).

4) You should get the RStudio welcome screen.  you should be able to
login with user/password `rstudio/rstudio`.

5) Clone the RNeXML repository using New Project from Version Control (https://github.com/ropensci/RNeXML), switch into the `manuscripts` directory and you should be good to go as above.  

## simmap NeXML definitions 

- Author: Carl Boettiger
- Initial version: 2014-03-21



Definitions of the `simmap` namespace, as defined for the use in `RNeXML`. The prefix `nex:` refers to the [NeXML schema](http://www.nexml.org/2009).


  term               | definition
 ------------------- | -------------
 `simmap:reconstructions`   | A container of one or more stochastic character map reconstructions, as a `meta` child of a the `nex:edge` element to which the contained stochastic character map reconstructions are being assigned.
 `simmap:reconstruction`    | A single stochastic character map reconstruction for a given `nex:edge`. Normally nested within a `simmap:reconstructions` element.
 `simmap:char`              | The id of a character trait, as defined by the `nex:char` element with this value as its `id`. This is a property of a `simmap:reconstruction`.
 `simmap:stateChange`       | A character state assignment to the given `nex:edge` during a specified interval, as a property of a `simmap:reconstruction`. Must have children `simmap:order`, `simmap:length`, and `simmap:state`.
`simmap:order`              | The chronological order (from the root) in which the state is assigned to the edge.  An edge that does not change states still has `simmap:order` 1.   This is a property of a `simmap:stateChange`.  
`simmap:length`             | The duration for which the edge occupies the assigned state, in the same units as the `nex:length` attribute defined on the `nex:edge` being annotated. This is a property of a `simmap:stateChange`.  
 `simmap:state`             | The id of a `nex:state` of the `nex:char` identified by the `simmap:char` property of the `simmap:reconstruction`. This is a property of a `simmap:stateChange`.  
``` r
library("RNeXML")
```

    ## Loading required package: ape

``` r
library("dplyr")
```

    ## 
    ## Attaching package: 'dplyr'
    ## 
    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag
    ## 
    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library("geiger")
knitr::opts_chunk$set(message = FALSE, comment = NA)
```

Let's generate a `NeXML` file using the tree and trait data from the `geiger` package's "primates" data:

``` r
data("primates")
add_trees(primates$phy) %>% 
  add_characters(primates$dat, ., append=TRUE) %>% 
  taxize_nexml() -> nex 
```

    Warning in taxize_nexml(.): ID for otu Alouatta_coibensis not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Aotus_hershkovitzi not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Aotus_miconax not found. Consider
    checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Callicebus_cinerascens not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Callicebus_dubius not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Callicebus_modestus not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Callicebus_oenanthe not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Callicebus_olallae not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Euoticus_pallidus not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Lagothrix_flavicauda not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Leontopithecus_caissara not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Leontopithecus_chrysomela not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Pithecia_aequatorialis not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Pithecia_albicans not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Procolobus_pennantii not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Procolobus_preussi not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Procolobus_rufomitratus not found.
    Consider checking the spelling or alternate classification

    Warning in taxize_nexml(.): ID for otu Tarsius_pumilus not found. Consider
    checking the spelling or alternate classification

(Note that we've used `dplyr`'s cute pipe syntax, but unfortunately our `add_` methods take the `nexml` object as the *second* argument instead of the first, so this isn't as elegant since we need the stupid `.` to show where the piped output should go...)

We now read in the three tables of interest. Note that we tell `get_characters` to give us species labels as there own column, rather than as rownames. The latter is the default only because this plays more nicely with the default format for character matrices that is expected by `geiger` and other phylogenetics packages, but is in general a silly choice for data manipulation.

``` r
otu_meta <- get_metadata(nex, "otus/otu")
taxa <- get_taxa(nex)
char <- get_characters(nex, rownames_as_col = TRUE)
```

We can take a peek at what the tables look like, just to orient ourselves:

``` r
otu_meta
```

    Source: local data frame [215 x 9]

          id property datatype content     xsi.type        rel
       (chr)    (lgl)    (lgl)   (lgl)        (chr)      (chr)
    1     m1       NA       NA      NA ResourceMeta tc:toTaxon
    2     m2       NA       NA      NA ResourceMeta tc:toTaxon
    3     m3       NA       NA      NA ResourceMeta tc:toTaxon
    4     m4       NA       NA      NA ResourceMeta tc:toTaxon
    5     m5       NA       NA      NA ResourceMeta tc:toTaxon
    6     m6       NA       NA      NA ResourceMeta tc:toTaxon
    7     m7       NA       NA      NA ResourceMeta tc:toTaxon
    8     m8       NA       NA      NA ResourceMeta tc:toTaxon
    9     m9       NA       NA      NA ResourceMeta tc:toTaxon
    10   m10       NA       NA      NA ResourceMeta tc:toTaxon
    ..   ...      ...      ...     ...          ...        ...
    Variables not shown: href (chr), otu (chr), otus (chr)

``` r
taxa
```

    Source: local data frame [233 x 5]

          id                       label about xsi.type  otus
       (chr)                       (chr) (chr)    (lgl) (chr)
    1    ou1 Allenopithecus_nigroviridis  #ou1       NA   os1
    2    ou2         Allocebus_trichotis  #ou2       NA   os1
    3    ou3           Alouatta_belzebul  #ou3       NA   os1
    4    ou4             Alouatta_caraya  #ou4       NA   os1
    5    ou5          Alouatta_coibensis  #ou5       NA   os1
    6    ou6              Alouatta_fusca  #ou6       NA   os1
    7    ou7           Alouatta_palliata  #ou7       NA   os1
    8    ou8              Alouatta_pigra  #ou8       NA   os1
    9    ou9               Alouatta_sara  #ou9       NA   os1
    10  ou10          Alouatta_seniculus #ou10       NA   os1
    ..   ...                         ...   ...      ...   ...

``` r
head(char)
```

                             taxa        x
    1 Allenopithecus_nigroviridis 8.465900
    2          Alouatta_seniculus 8.767173
    3               Galago_alleni 5.521461
    4             Galago_gallarum 5.365976
    5            Galago_matschiei 5.267858
    6               Galago_moholi 5.375278

Now that we have nice `data.frame` objects for all our data, it's easy to join them into the desired table with a few obvious `dplyr` commands:

``` r
taxa %>% 
  left_join(char, by = c("label" = "taxa")) %>% 
  left_join(otu_meta, by = c("id" = "otu")) %>%
  select(id, label, x, href)
```

    Warning in left_join_impl(x, y, by$x, by$y): joining factor and character
    vector, coercing into character vector

    Source: local data frame [233 x 4]

          id                       label        x
       (chr)                       (chr)    (dbl)
    1    ou1 Allenopithecus_nigroviridis 8.465900
    2    ou2         Allocebus_trichotis 4.368181
    3    ou3           Alouatta_belzebul 8.729074
    4    ou4             Alouatta_caraya 8.628735
    5    ou5          Alouatta_coibensis 8.764053
    6    ou6              Alouatta_fusca 8.554489
    7    ou7           Alouatta_palliata 8.791790
    8    ou8              Alouatta_pigra 8.881836
    9    ou9               Alouatta_sara 8.796339
    10  ou10          Alouatta_seniculus 8.767173
    ..   ...                         ...      ...
    Variables not shown: href (chr)

Because these are all from the same otus block anyway, we haven't selected that column, but were it of interest it is also available in the taxa table.
---
output: github_document
---

<!-- README.md is generated from README.Rmd. Please edit that file -->

```{r, echo = FALSE}
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.path = "tools/README-"
)
```


# RNeXML  <img src="man/figures/logo.svg" align="right" alt="" width="120" />

[![DOI](https://zenodo.org/badge/11856817.svg)](https://zenodo.org/badge/latestdoi/11856817)
[![Build Status](https://api.travis-ci.org/ropensci/RNeXML.png)](https://travis-ci.org/ropensci/RNeXML)
[![Windows build status](https://ci.appveyor.com/api/projects/status/dhiwp5blx2ns2yba/branch/master?svg=true)](https://ci.appveyor.com/project/cboettig/rnexml/branch/master)
[![CRAN status](https://www.r-pkg.org/badges/version/RNeXML)](https://cran.r-project.org/package=RNeXML)
[![codecov.io](https://codecov.io/github/ropensci/RNeXML/coverage.svg?branch=master)](https://codecov.io/github/ropensci/RNeXML?branch=master)
![downloads](http://cranlogs.r-pkg.org/badges/grand-total/RNeXML)



* Maintainer: Carl Boettiger
* Authors: Carl Boettiger, Scott Chamberlain, Hilmar Lapp, Kseniia Shumelchyk, Rutger Vos
* License: BSD-3 
* [Issues](https://github.com/ropensci/RNeXML/issues): Bug reports, feature requests, and development discussion.

An extensive and rapidly growing collection of richly annotated phylogenetics data is now available in the NeXML format. NeXML relies on state-of-the-art data exchange technology to provide a format that can be both validated and extended, providing a data quality assurance and adaptability to the future that is lacking in other formats. See [Vos et al 2012](http://doi.org/10.1093/sysbio/sys025 "NeXML: Rich, Extensible, and Verifiable Representation of Comparative Data and Metadata.") for further details on the NeXML format.

How to cite
-----------

RNeXML has been published in the following article:

> Boettiger C, Chamberlain S, Vos R and Lapp H (2016). “RNeXML: A Package for
> Reading and Writing Richly Annotated Phylogenetic, Character, and Trait Data
> in R.” _Methods in Ecology and Evolution_, **7**, pp. 352-357.
> [doi:10.1111/2041-210X.12469](http://doi.org/10.1111/2041-210X.12469)

Although the published version of the article is pay-walled, the source of the manuscript, and a much better rendered PDF, are included in this package (in the `manuscripts` folder). You can also find it [freely available on arXiv](http://arxiv.org/abs/1506.02722).


```{r compile-settings, include=FALSE}
library(methods)
```

Installation
---------------

The latest stable release of RNeXML is on CRAN, and can be installed with the usual `install.packages("RNeXML")` command.  Some of the more specialized functionality described in the Vignettes (such as RDF manipulation) requires additional packages which can be installed using:

```{r eval=FALSE}
install.packages("RNeXML", deps = TRUE)
```

The development version can be installed using:

```{r eval=FALSE}
remotes::install_github("ropensci/RNeXML")
```


Getting Started
----------------

See the vignettes below for both a general quick start and an overview of more specialized features.

- [A Brief Introduction to RNeXML](https://docs.ropensci.org/RNeXML/articles/intro)
- [RNeXML: A Package for Reading and Writing Richly Annotated Phylogenetic, Character, and Trait Data in R](https://github.com/ropensci/RNeXML/tree/master/manuscripts) (published in MEE).
- [Handling Metadata in RNeXML](https://docs.ropensci.org/RNeXML/articles/metadata)
- [The `nexml` S4 Object](https://docs.ropensci.org/RNeXML/articles/S4)
- [Semantic data & SPARQL with RNeXML](https://docs.ropensci.org/RNeXML/articles/sparql)
- [Extending NeXML: an example based on simmap](https://docs.ropensci.org/RNeXML/articles/simmap)
 



-----

[![ropensci footer](http://ropensci.org/public_images/github_footer.png)](http://ropensci.org)

---
#layout: preprint
layout: review, 11pt
linenumbers: true
title: "RNeXML: a package for reading and writing richly annotated phylogenetic, character, and trait data in R"
author: 
  - name: Carl Boettiger
    affiliation: ucb
    email: cboettig(at)gmail.com
    footnote: Corresponding author
  - name: Scott Chamberlain
    affiliation: ropensci
  - name: Rutger Vos
    affiliation: NBC
  - name: Hilmar Lapp
    affiliation: dukeplus
address: 
  - code: ucb 
    address: |
      University of California, Berkeley,
      130 Mulford Hall \#3114,
      Berkeley, CA 94720-3114, USA 
  - code: ropensci
    address: |
      University of California, Berkeley, CA, USA
  - code: NBC
    address: | 
      Naturalis Biodiversity Center, Leiden, the Netherlands
  - code: dukeplus
    address: | 
      Center for Genomic and Computational Biology, Duke University,
      and 
      National Evolutionary Synthesis Center, Durham, NC, USA
abstract: | 
      1. NeXML is a powerful and extensible exchange standard
      recently proposed to better meet the expanding needs for
      phylogenetic data and metadata sharing. Here we
      present the RNeXML package, which provides users of
      the R programming language with easy-to-use tools for 
      reading and writing NeXML documents, including rich metadata, in
      a way that interfaces seamlessly with the extensive library of
      phylogenetic tools already available in the R ecosystem.

      2. Wherever possible, we designed RNeXML to map NeXML document
      contents, whose arrangement is influenced by the format's
      XML Schema definition, to their most intuitive or
      useful representation in R. To make NeXML's powerful facility
      for recording semantically rich machine-readable metadata
      accessible to R users, we designed a functional programming
      interface to it that hides the semantic web standards leveraged
      by NeXML from R users who are unfamiliar with them.

      3. RNeXML can read any NeXML document that validates, and it
      generates valid NeXML documents from phylogeny and character
      data in various R representations in use. The metadata
      programming interface at a basic level aids fulfilling data
      documentation best practices, and at an advanced level preserves
      NeXML's nearly limitless extensibility, for which we provide a
      fully working demonstration. Furthermore, to lower the barriers
      to sharing well-documented phylogenetic data, RNeXML has started
      to integrate with taxonomic metadata augmentation services on
      the web, and with online repositories for data archiving.

      4. RNeXML allows R's rich ecosystem to read and write data in
      the NeXML format through an interface that is no more involved
      than reading or writing data from other, less powerful data
      formats. It also provides an interface designed to feel familiar
      to R programmers and to be consistent with recommended practices
      for R package development, yet that retains the full power for
      users to add their own custom data and metadata to the
      phylogenies they work with, without introducing potentially
      incompatible changes to the exchange standard.
      
bibliography: components/references.bib
csl: components/methods-in-ecology-and-evolution.csl
documentclass: components/elsarticle

output: 
  pdf_document:
    template: components/elsarticle.latex
    keep_tex: true
    fig_caption: true


---


```{r compile-settings, include=FALSE}
library("methods")
library("knitr")
opts_chunk$set(tidy = FALSE, warning = FALSE, message = FALSE, 
               cache = 1, comment = NA, verbose = TRUE)
basename <- gsub(".Rmd", "", knitr:::knit_concord$get('infile')) 
opts_chunk$set(fig.path = paste("components/figure/", basename, "-", sep=""),
               cache.path = paste("components/cache/", basename, "/", sep=""))


```

# Introduction

Users of the popular statistical and mathematical computing platform
R [@R] enjoy a wealth of readily installable comparative phylogenetic
methods and tools [@taskview]. Exploiting the opportunities arising from
this wealth for complex and integrative comparative research questions
relies on the ability to reuse and integrate previously generated or
published data and metadata. The expanding data exchange needs of the
evolutionary research community are rapidly outpacing the capabilities of
most current and widely used data exchange standards [@Vos_2012], which
were all developed a decade or more ago. This has resulted in a radiation
of different data representations and exchange standard "flavors" that are
no longer interoperable at the very time when the growth of available data
and methods has made that interoperability most valuable. In response to
the unmet needs for standardized data exchange in phylogenetics, a modern
XML-based exchange standard, called NeXML, has recently been developed
[@Vos_2012]. NeXML comprehensively supports current data exchange needs,
is predictably machine-readable, and is forward compatible.

The exchange problem for phylogenetic data is particularly acute in light
of the challenges in finding and sharing phylogenetic data without the
otherwise common loss of most data and metadata semantics [@Drew_2013;
@Stoltzfus_2012; @Cranston_2014]. For example, the still popular NEXUS
file format [@Maddison_1997] cannot consistently represent horizontal
gene transfer or ambiguity in reading a character (such as a DNA sequence
base pair).  This and other limitations have led to modifications of NEXUS
in different ways for different needs, with the unfortunate result that
NEXUS files generated by one program can be incompatible with another
[@Vos_2012]. Without a formal grammar, software based on NEXUS files
may also make inconsistent assumptions about tokens, quoting, or element
lengths.  @Vos_2012 estimates that as many as 15% of the NEXUS files in
the CIPRES portal contain unrecoverable but hard to diagnose errors.

A detailed account of how the NeXML standard addresses these and
other relevant challenges can be found in @Vos_2012. In brief,
NeXML was designed with the following important properties. First,
NeXML is defined by a precise grammar that can be programmatically
**validated**; i.e., it can be verified whether a file precisely
follows this grammar, and therefore whether it can be read (parsed) 
without errors by software that uses the NeXML grammar (e.g. RNeXML)
is predictable. 
Second, NeXML is **extensible**: a user
can define representations of new, previously unanticipated information
(as we will illustrate) without violating its defining grammar. Third
and most importantly, NeXML is rich in **computable semantics**: it is
designed for expressing metadata such that machines can understand
their meaning and make inferences from it. For example, OTUs in a tree
or character matrix for frog species can be linked to concepts in a
formally defined hierarchy of taxonomic concepts such as the
Vertebrate Taxonomy Ontology [@Midford2013], which enables a machine
to infer that a query for amphibia is to include the frog data in what
is returned. (For a more broader discussion of the value of such
capabilities for evolutionary and biodiversity science we refer the
reader to @Parr2011.)


To make the capabilities of NeXML available to R users in an easy-to-use
form, and to lower the hurdles to adoption of the standard, we present
RNeXML, an R package that aims to provide easy programmatic access to
reading and writing NeXML documents, tailored for the kinds of use-cases
that will be common for users and developers of the wealth of evolutionary
analysis methods within the R ecosystem.

```{r echo=FALSE, results='hide', eval=FALSE}
install.packages("RNeXML", dependencies=TRUE)
```

```{r echo=FALSE, results='hide'}
library(RNeXML)
```

# The RNeXML package

The `RNeXML` package is written entirely in R and available under a
Creative Commons Zero public domain waiver. The current development
version can be found on Github at [https://github.com/ropensci/RNeXML](),
and the stable version can be installed from the CRAN repository.
`RNeXML` is part of the rOpenSci project. Users of `RNeXML` are encouraged
to submit bug reports or feature requests in the issues log on Github,
or the phylogenetics R users group list at `r-sig-phylo@r-project.org`
for help. Vignettes with more detailed examples of specific features
of RNeXML are distributed with the R package and serve as a supplement
to this manuscript. Each of the vignettes can be found at
[http://ropensci.github.io/RNeXML/]().

## Representation of NeXML documents in R

Conceptually, a NeXML document has the following components: (1)
phylogeny topology and branch length data, (2) character or trait data in matrix form,
(3) operational taxonomic units (OTUs), and (4) metadata. To represent
the contents of a NeXML document (currently in memory), `RNeXML` defines
the `nexml` object type. This type therefore holds phylogenetic trees as
well as character or trait matrices, and all metadata, which is similar
to the phylogenetic data object types defined in the `phylobase` package
[@phylobase], but contrasts with the more widely used ones defined in the
`ape` package [@Paradis_2004], which represents trees alone.

When reading and writing NeXML documents, `RNeXML` aims to map
their components to and from, respectively, their most widely used
representations in R. As a result, the types of objects accepted
or returned by the package's methods are the `phylo` and `multiPhylo`
objects from the `ape` package [@Paradis_2004] for phylogenies, and R's
native `data.frame` list structure for data matrices.

## Reading phylogenies and character data

The method `nexml_read()` reads NeXML files, either from a local file, or
from a remote location via its URL, and returns an object of type `nexml`:

```{r}
nex <- nexml_read("components/trees.xml")
```

The method `get_trees_list()` can be used to extract the phylogenies
as an `ape::multiPhylo`
object, which can be treated as a list of `ape::phylo` objects:

```{r}
phy <- get_trees_list(nex)
```

The `get_trees_list()` method is designed for use in scripts, providing
a consistent and predictable return type regardless of the number
of phylogenies a NeXML document contains. For greater convenience in
interactive use, the method `get_trees()` returns the R object most
intuitive given the arrangement of phylogeny data in the source NeXML
document. For example, the method returns an `ape::phylo` object if
the NeXML document contains a single phylogeny, an `ape::multiPhylo`
object if it contains multiple phylogenies arranged in a single `trees`
block, and a list of `ape::multiPhylo` objects if it contains multiple
`trees` blocks (the capability for which NeXML inherits from NEXUS).

If the location parameter with which the `nexml_read()` method is
invoked is recognized as a URL, the method will automatically download
the document to the local working directory and read it from there. This
gives convenient and rapid access to phylogenetic data published in
NeXML format on the web, such as the content of the phylogenetic data
repository TreeBASE [@Piel_2002; @Piel_2009]. For example, the following plots a
tree in TreeBASE (using ape's plot function):

```{r eval=FALSE, fig.show="hide"}
tb_nex <- nexml_read(
"https://raw.github.com/TreeBASE/supertreebase/master/data/treebase/S100.xml")
tb_phy <- get_trees_list(tb_nex)
plot(tb_phy[[1]]) 
```

The method `get_characters()` obtains character data matrices from a
`nexml` object, and returns them as a standard `data.frame` R object
with columns as characters and rows as taxa:

```{r }
nex <- nexml_read("components/comp_analysis.xml")
get_characters(nex)
```

A NeXML data matrix can be of molecular (for molecular sequence
alignments), discrete (for most morphological character data), or
continuous type (for many trait data). To enable strict validation of data
types NeXML allows, and if their data types differ requires multiple data
matrices to be separated into different "blocks". Since the `data.frame`
data structure in R has no such constraints, the `get_characters()`
method combines such blocks as separate columns into a single `data.frame`
object, provided they correspond to the same taxa. Otherwise, a list of
`data.frame`s is returned, with list elements corresponding to characters
blocks. Similar to the methods for obtaining trees, there is also a method
`get_characters_list()`, which always returns a list of `data.frame`s,
one for each character block.

## Writing phylogenies and character data

The method `nexml_write()` generates a NeXML file from its input
parameters. In its simplest invocation, the method writes a tree to
a file:

```{r, results='hide'}
data(bird.orders)
nexml_write(bird.orders, file = "birds.xml")
```

The first argument to `nexml_write()` is either an object of type `nexml`,
or any object that can be coerced to it, such as in the above example an
`ape::phylo` phylogeny. Alternatively, passing a `multiPhylo` object
would write a list of phylogenies to the file.

In addition to trees, the `nexml_write()` method also allows to specify
character data as another parameter. The following example uses data
from the comparative phylogenetics R package `geiger` [@Pennell2014].

```{r, results='hide'}
library("geiger")
data(geospiza)
nexml_write(trees = geospiza$phy, 
            characters = geospiza$dat,
            file="geospiza.xml")
```

Note that the NeXML format is well-suited for incomplete data: for instance,
here it does not assume the character matrix has data for every tip in the tree.

## Validating NeXML

File validation is a central feature of the NeXML format which ensures
that any properly implemented NeXML parser will always be able to read
the NeXML file. The function takes the path to any NeXML file and returns 
`TRUE` to indicate a valid file, or `FALSE` otherwise, along with a
display of any error messages generated by the validator.

```{r }
nexml_validate("geospiza.xml")
```


The `nexml_validate()` function performs this validation
using the online NeXML validator (when a network connection is available),
which performs additional checks not expressed in the NeXML schema itself [@Vos_2012].
If a network connection is not available, the function falls back on the 
schema validation method from the `XML` package [@XML]. 


## Creating and populating `nexml` objects

Instead of packaging the various components for a NeXML file at the
time of writing the file, `RNeXML` also allows users to create and
iteratively populate in-memory `nexml` objects. The methods to do
this are `add_characters()`, `add_trees()`, and `add_meta()`, for
adding characters, trees, and metadata, respectively. Each of these
functions will automatically create a new nexml object if not supplied
with an existing one as the last (optional) argument.

For example, here we use `add_trees()` to first create a `nexml`
object with the phylogeny data, and then add the character data to it:

```{r}
nexObj <- add_trees(geospiza$phy)
nexObj <- add_characters(geospiza$dat, nexObj)
```

The data with which a `nexml` object is populated need not share the
same OTUs.  `RNeXML` automatically adds new, separate OTU blocks into
the NeXML file for each data matrix and tree that uses a different set of OTUs.

Other than storage size, there is no limit to the number of
phylogenies and character matrices that can be included in a
single NeXML document. This allows, for example, to capture samples from a
posterior probability distribution of inferred or simulated phylogenies and
character states in a single NeXML file.

## Data documentation and annotation with built-in metadata

NeXML allows attaching ("_annotating_") metadata to any data element,
and even to metadata themselves. Whether at the level of the document as
a whole or an individual data matrix or phylogeny, metadata can provide
bibliographic and provenance information, for example about the study
as part of which the phylogeny was generated or applied, which data
matrix and which methods were used to generate it. Metadata can also be
attached to very specific elements of the data, such as specific traits,
individual OTUs, nodes, or even edges of the phylogeny.

As described in @Vos_2012, to encode metadata annotations NeXML uses the
"Resource Description Framework in Annotations" (RDFa) [@W3C_2014]. This
standard provides for a strict machine-readable format yet enables future
backwards compatibility with compliant NeXML parsers (and thus `RNeXML`),
because the capacity of a tool to _parse_ annotations is not predicated
on _understanding_ the meaning of annotations it has not seen before.

To lower the barriers to sharing well-documented phylogenetic data,
`RNeXML` aims to make recording useful and machine-readable metadata
easier at several levels.

First, when writing a NeXML file the package adds certain basic metadata
automatically if they are absent, using default values consistent with
recommended best practices [@Cranston_2014]. Currently, this includes
naming the software generating the NeXML, a time-stamp of when a tree
was produced, and an open data license. These are merely default arguments
to `add_basic_meta()` and can be configured.

Second, `RNeXML` provides a simple method, called `add_basic_metadata()`,
to set metadata attributes commonly recommended for inclusion with data
to be publicly archived or shared [@Cranston_2014]. The currently accepted
parameters include `title`, `description`, `creator`, `pubdate`, `rights`,
`publisher`, and `citation`. Behind the scenes the method automatically
anchors these attributes in common vocabularies (such as Dublin Core).

Third, `RNeXML` integrates with the R package `taxize` [@Chamberlain_2013]
to mitigate one of the most common obstacles to reuse of phylogenetic
data, namely the misspellings and inconsistent taxonomic naming with
which OTU labels are often fraught. The `taxize_nexml()` method in
`RNeXML` uses `taxize` to match OTU labels against the NCBI database,
and, where a unique match is found, it annotates the respective OTU with
the matching NCBI identifier.


## Data annotation with custom metadata

The `RNeXML` interface described above for built-in metadata allows
users to create precise and semantically rich annotations without 
confronting any of the complexity of namespaces and ontologies.
Nevertheless, advanced users may desire the explicit control over
these semantic tools that takes full advantage of the flexibility
and extensibility of the NeXML specification [@Vos_2012; @Parr2011].
In this section we detail how to accomplish these more complex uses
in RNeXML.


Using a vocabulary or ontology terms rather than simple text strings to 
describe data is crucial for allowing
machines to not only parse but also interpret and potentially reason
over their semantics. To achieve this benefit for custom metadata
extensions, the user necessarily needs to handle certain technical
details from which the `RNeXML` interface shields her otherwise, in particular
the globally unique identifiers (normally HTTP URIs) of metadata terms
and vocabularies. To be consistent with XML terminology, `RNeXML` calls
vocabulary URIs _namespaces_, and their abbreviations _prefixes_. For
example, the namespace for the Dublin Core Metadata Terms vocabulary
is "http://purl.org/dc/elements/1.1/". Using its common abbreviation
"dc", a metadata property "dc:title" expands to the identifier
"http://purl.org/dc/elements/1.1/title". This URI resolves to a human and
machine-readable (depending on access) definition of precisely what the
term `title` in Dublin Core means. In contrast, just using the text string
"title" could also mean the title of a person, a legal title, the verb
title, etc.  URI identifiers of metadata vocabularies and terms are not
mandated to resolve, but if machines are to derive the maximum benefit
from them, they should resolve to a definition of their semantics in RDF.

`RNeXML` includes methods to obtain and manipulate metadata properties,
values, identifiers, and namespaces. The `get_namespaces()` method
accepts a `nexml` object and returns a named list of namespace prefixes
and their corresponding identifiers known to the object:


```{r}
birds <- nexml_read("birds.xml")
prefixes <- get_namespaces(birds)
prefixes["dc"]
```

The `get_metadata()` method returns, as a named list, the metadata
annotations for a given `nexml` object at a given level, with the whole
NeXML document being the default level (`"all"` extracts all metadata
objects):

```{r}
meta <- get_metadata(birds) 
otu_meta <- get_metadata(birds, level="otu")
```

The returned list does not include the data elements to which the
metadata are attached. Therefore, a different approach, documented in
the metadata vignette, is recommended for accessing the metadata attached
to data elements.

The `meta()` method creates a new metadata object from a property name
and content (value). For example, the following creates a modification
date metadata object, using a property in the PRISM vocabulary:


```{r}
modified <- meta(property = "prism:modificationDate", content = "2013-10-04")
```

Metadata annotations in `NeXML` can be nested within another
annotation, which the `meta()` method accommodates by accepting a
parameter `children`, with the list of nested metadata objects (which
can themselves be nested) as value.

The `add_meta()` function adds metadata objects as annotations to a
`nexml` object at a specified level, with the default level being the
NeXML document as a whole:

```{r}
birds <- add_meta(modified, birds) 
```

If the prefix used by the metadata property is not among the built-in
ones (which can be obtained using `get_namespaces()`), it has to be
provided along with its URI as the `namespaces` parameter. For example,
the following uses the "[Simple Knowledge Organization System](http://www.w3.org/TR/skos-reference/)" (SKOS) 
vocabulary to add a note to the trees in the `nexml` object:

```{r}
history <- meta(property = "skos:historyNote",
  content = "Mapped from the bird.orders data in the ape package using RNeXML")
birds <- add_meta(history, 
                birds, 
                level = "trees",
                namespaces = c(skos = "http://www.w3.org/2004/02/skos/core#"))
```

Alternatively, additional namespaces can also be added in batch using
the `add_namespaces()` method.

By virtue of subsetting the S4 `nexml` object, `RNeXML` also offers
fine control of where a `meta` element is added, for which the package
vignette on S4 subsetting of `nexml` contains examples.

Because NeXML expresses all metadata using the RDF standard, and stores
them compliant with RDFa, they can be extracted as an RDF graph, queried,
analyzed, and mashed up with other RDF data, local or on the web, using
a wealth of off-the-shelf tools for working with RDF (see @W3C_2014
or @Hartig_2012). Examples for these possibilities are included in the
`RNeXML` SPARQL vignette (a recursive acronym for SPARQL Protocol and
RDF Query Language, see http://www.w3.org/TR/rdf-sparql-query/), and
the package also comes with a demonstration that can be run from R using
the following command: `demo("sparql", "RNeXML")`).

## Using metadata to extend the NeXML standard 

NeXML was designed to prevent the need for future non-interoperable
"flavors" of the standard in response to new research directions. Its
solution to this inevitable problem is a highly flexible metadata
system without sacrificing strict validation of syntax and structure.

Here we illustrate how `RNeXML`'s interface to NeXML's metadata system
can be used to record and share a type of phylogenetic data not taken
into account when NeXML was designed, in this case stochastic character
maps [@Huelsenbeck_2003]. Such data assign certain parts (corresponding
to time) of each branch in a time-calibrated phylogeny to a particular
"state" (typically of a morphological characteristic). The current
de-facto format for sharing stochastic character maps, created by
`simmap` [@Bollback_2006], a widely used tool for creating such maps,
is a non-interoperable modification of the standard Newick tree format.
This means that computer programs designed to read Newick or NEXUS formats
may fail when trying to read in a phylogeny that includes `simmap` annotations.

In contrast, by allowing new data types to be added as --- sometimes
complex --- metadata annotations NeXML can accommodate data extensions
without compromise to its grammar and thus syntax In NeXML. To illustrate
how RNeXML facilitates extending the NeXML standard in this way, we
have implemented two functions in the package, `nexml_to_simmap` and
`simmap_to_nexml`.  These functions show
how simmap data can be represented as `meta` annotations on the branch
length elements of a NeXML tree, and provide routines to convert between
this NeXML representation and the extended `ape::phylo` representation
of a `simmap` tree in R that was introduced by @Revell_2012. We encourage
readers interested in this capability to consult the example code in
`simmap_to_nexml` to see how this is implemented. 

Extensions to NeXML must also be defined in the file's namespace 
in order to valid.  This provides a way to ensure that a URI providing
documentation of the extension is always included.  Our examples here
use the prefix, `simmap`, to group the newly introduced
metadata properties in a vocabulary, for which the `add_namespace()`
method can be used to give a URI as an identifier:

```{r}
nex <- add_namespaces(c(simmap = 
  "https://github.com/ropensci/RNeXML/tree/master/inst/simmap.md"))
```

Here the URI does not resolve to a fully machine-readable definition of
the terms and their semantics, but it can nonetheless be used to provide
at least a human-readable informal definition of the terms.


## Publishing NeXML files from R

Data archiving is increasingly required by scientific journals,
including in evolutionary biology, ecology, and biodiversity
(e.g. @Rausher_2010). The effort involved with preparing and
submitting properly annotated data to archives remains a notable
barrier to the broad adoption of data archiving and sharing as a
normal part of the scholarly publication workflow
[@Tenopir_2011; @Stodden_2014]. In particular, the majority of
phylogenetic trees published in the scholarly record are inaccessible
or lost to the research community [@Drew_2013].

One of `RNeXML`'s aims is to promote the archival of well-documented
phylogenetic data in scientific data repositories, in the form of
NeXML files. To this end, the method `nexml_publish()` provides an API
directly from within R that allows data archival to become a step programmed
into data management scripts. Initially, the method supports the data repository Figshare (http://figshare.com):

```{r eval = FALSE}
doi <- nexml_publish(birds, repository="figshare")
```

This method reserves a permanent identifier (DOI) on the figshare
repository that can later be made public through the figshare web
interface.  This also acts as a secure backup of the data to a repository
and a way to share with collaborators prior to public release.

# Conclusions and future directions

`RNeXML` allows R's ecosystem to read and write data in the NeXML
format through an interface that is no more involved than reading or
writing data from other phylogenetic data formats. It also carries
immediate benefits for its users compared to other formats. For
example, comparative analysis R packages and users frequently add
their own metadata annotations to the phylogenies they work with, such
as annotations of species, stochastic character maps, trait values,
model estimates and parameter values. `RNeXML` affords R the
capability to harness machine-readable semantics and an extensible
metadata schema to capture, preserve, and share these and other kinds
of information, all through an API instead of having to understand in
detail the schema underlying the NeXML standard. To assist users in
meeting the rising bar for best practices in data sharing in
phylogenetic research [@Cranston_2014], `RNeXML` captures metadata
information from the R environment to the extent possible, and applies
reasonable defaults.

The goals for continued development of `RNeXML` revolve primarily
around better interoperability with other existing phylogenetic data
representations in R, such as those found in the `phylobase` package
[@phylobase]; and better integration of the rich metadata semantics
found in ontologies defined in the Web Ontology Language (OWL),
including programmatic access to machine reasoning with such metadata.

## Acknowledgements

This project was supported in part by the National Evolutionary
Synthesis Center (NESCent) (NSF #EF-0905606), and grants from the
National Science Foundation (DBI-1306697) and the Alfred P Sloan
Foundation (Grant 2013-6-22). `RNeXML` started as a project idea for
the Google Summer of Code(TM), and we thank Kseniia Shumelchyk for taking
the first steps to implement it. We are grateful to F. Michonneau for
helpful comments on an earlier version of this manuscript, and reviews
by Matthew Pennell, Associate Editor Richard FitzJohn, and an anonymous
reviewer. At their behest, the reviews of FitzJohn and Pennell can be found in this
project's GitHub page at [github.com/ropensci/RNeXML/issues/121](https://github.com/ropensci/RNeXML/issues/121) and [github.com/ropensci/RNeXML/issues/120](https://github.com/ropensci/RNeXML/issues/120), together with our replies and a record of our revisions.

## Data Accessibility

All software, scripts and data used in this paper can be found in the permanent 
data archive Zenodo under the digital object identifier
 doi:10.5281/zenodo.13131 [@zenodo]. This DOI corresponds to a snapshot of the GitHub repository
at [github.com/ropensci/RNeXML](https://github.com/ropensci/RNeXML).


```{r cleanup, include=FALSE, cache=FALSE}
unlink("birds.xml")
unlink("geospiza.xml")
```

# References 
---
output:
 md_document:
   variant: markdown_github
---




```{r}
library("RNeXML")
library("dplyr")
library("geiger")
knitr::opts_chunk$set(message = FALSE, warning=FALSE, comment = NA)
```

Let's generate a `NeXML` file using the tree and trait data from the `geiger` package's "primates" data:

```{r}
data("primates")
add_trees(primates$phy) %>% 
  add_characters(primates$dat, ., append=TRUE) %>% 
  taxize_nexml() -> nex 
```

(Note that we've used `dplyr`'s cute pipe syntax, but unfortunately our `add_` methods take the `nexml` object as the _second_
argument instead of the first, so this isn't as elegant since we need the stupid `.` to show where the piped output should go...)

We now read in the three tables of interest.  Note that we tell `get_characters` to give us species labels as there own column, rather than as rownames.  The latter is the default only because this plays more nicely with the default format for character matrices that is expected by `geiger` and other phylogenetics packages, but is in general a silly choice for data manipulation. 

```{r}
otu_meta <- get_metadata(nex, "otus/otu")
taxa <- get_taxa(nex)
char <- get_characters(nex, rownames_as_col = TRUE)
```


We can take a peek at what the tables look like, just to orient ourselves:

```{r}
otu_meta
taxa
head(char)
```

Now that we have nice `data.frame` objects for all our data, it's easy to join them into the desired table with a few obvious `dplyr` commands:

```{r}
taxa %>% 
  left_join(char, by = c("label" = "taxa")) %>% 
  left_join(otu_meta, by = c("id" = "otu")) %>%
  select(id, label, x, href)
```

Because these are all from the same otus block anyway, we haven't selected that column, but were it of interest it is also available in the taxa table.

---
title: "A Brief Introduction to RNeXML"
author:
  - Carl Boettiger
  - Scott Chamberlain
  - Rutger Vos
  - Hilmar Lapp

output: html_vignette
---

<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{A Brief Introduction to RNeXMLL}
-->




Read in a `nexml` file:


```r
f <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nexml <- nexml_read(f)
```

Pretty-print an overview of the components and metadata that make up the nexml object:

```r
nexml  # this is the same as show(nexml)
```

```
## A nexml object representing:
##  	 1 phylogenetic tree block(s), where:
## 	   block 1 contains 1 phylogenetic tree(s) 
##  	 2 character block(s), where:
## 	   block 1 defines 1 continuous character(s) 
## 	   matrix 1 has 10 row(s)
## 	   block 2 defines 1 standard/discrete character(s), with 2 states each 
## 	    and  0 polymorphic or uncertain state(s) defined 
## 	   matrix 2 has 10 row(s) 
##  	 10 taxonomic units in 1 block(s) 
##   Taxa:	 taxon_1, taxon_2, taxon_3, taxon_4, taxon_5, taxon_6 ...
##   Metadata annotations: 
## 	2 at top level 
## 	0 in block 1 at otu level 
## 	0 in block 1 at char level
## 	0 in block 2 at char level 
## 	0 in block 1 at state level
## 	0 in block 2 at state level
## 
## Author(s): rvosa
## 
## NeXML generated by Bio::Phylo::Project v.0.56 using schema version: 0.9 
## Size: 289.3 Kb
```

Create a summary object of various component and metadata counts (the `show()` method uses this):

```r
summary(nexml)
```

```
## $nblocks
##      trees       otus characters 
##          1          1          2 
## 
## $ncharacters
## block.1 block.2 
##       1       1 
## 
## $nstates
##         block.1 block.2
## Min.         NA       2
## 1st Qu.      NA       2
## Median       NA       2
## Mean         NA       2
## 3rd Qu.      NA       2
## Max.         NA       2
## 
## $nnonstdstatedefs
##         polymorphic uncertain
## block.1          NA        NA
## block.2           0         0
## 
## $nmatrixrows
## block.1 block.2 
##      10      10 
## 
## $ntrees
## block.1 
##       1 
## 
## $notus
## block.1 
##      10 
## 
## $nmeta
## $nmeta$nexml
## [1] 2
## 
## $nmeta$otu
## block.1 
##       0 
## 
## $nmeta$char
## block.1 block.2 
##       0       0 
## 
## $nmeta$state
## block.1 block.2 
##       0       0
```

Extract trees from nexml into the `ape::phylo` format:

```r
tr <- get_trees(nexml) # or: as(nexml, "phylo")
plot(tr)
```

![plot of chunk unnamed-chunk-5](intro-unnamed-chunk-5-1.png)

Write an `ape::phylo` tree into the `nexml` format:


```r
data(bird.orders)
nexml_write(bird.orders, "test.xml", creator = "Carl Boettiger")
```

```
## [1] "test.xml"
```

A key feature of NeXML is the ability to formally validate the construction of the data file against the standard (the lack of such a feature in nexus files had lead to inconsistencies across different software platforms, and some files that cannot be read at all).  While it is difficult to make an invalid NeXML file from `RNeXML`, it never hurts to validate just to be sure:


```r
nexml_validate("test.xml")
```

```
## [1] TRUE
```



Extract metadata from the NeXML file: 


```r
birds <- nexml_read("test.xml")
get_taxa(birds)
```

```
##     otu            label xsi.type otus
## 1  ou37 Struthioniformes       NA  os3
## 2  ou38     Tinamiformes       NA  os3
## 3  ou39      Craciformes       NA  os3
## 4  ou40      Galliformes       NA  os3
## 5  ou41     Anseriformes       NA  os3
## 6  ou42    Turniciformes       NA  os3
## 7  ou43       Piciformes       NA  os3
## 8  ou44    Galbuliformes       NA  os3
## 9  ou45   Bucerotiformes       NA  os3
## 10 ou46      Upupiformes       NA  os3
## 11 ou47    Trogoniformes       NA  os3
## 12 ou48    Coraciiformes       NA  os3
## 13 ou49      Coliiformes       NA  os3
## 14 ou50     Cuculiformes       NA  os3
## 15 ou51   Psittaciformes       NA  os3
## 16 ou52      Apodiformes       NA  os3
## 17 ou53   Trochiliformes       NA  os3
## 18 ou54  Musophagiformes       NA  os3
## 19 ou55     Strigiformes       NA  os3
## 20 ou56    Columbiformes       NA  os3
## 21 ou57       Gruiformes       NA  os3
## 22 ou58    Ciconiiformes       NA  os3
## 23 ou59    Passeriformes       NA  os3
```

```r
get_metadata(birds) 
```

```
##           property   datatype                 content     xsi.type                                              href Meta
## 1       dc:creator xsd:string          Carl Boettiger  LiteralMeta                                              <NA> m278
## 2 dcterms:modified xsd:string 2020-01-28 23:24:19 GMT  LiteralMeta                                              <NA> m279
## 3       cc:license       <NA>                    <NA> ResourceMeta http://creativecommons.org/publicdomain/zero/1.0/ m280
```

--------------------------------------------


Add basic additional metadata:  


```r
  nexml_write(bird.orders, file="meta_example.xml",
              title = "My test title",
              description = "A description of my test",
              creator = "Carl Boettiger <cboettig@gmail.com>",
              publisher = "unpublished data",
              pubdate = "2012-04-01")
```

```
## [1] "meta_example.xml"
```
By default, `RNeXML` adds certain metadata, including the NCBI taxon id numbers for all named taxa.  This acts a check on the spelling and definitions of the taxa as well as providing a link to additional metadata about each taxonomic unit described in the dataset.  


### Advanced annotation


We can also add arbitrary metadata to a NeXML tree by define `meta` objects:


```r
modified <- meta(property = "prism:modificationDate",
                 content = "2013-10-04")
```

Advanced use requires specifying the namespace used.  Metadata follows the RDFa conventions.  Here we indicate the modification date using the prism vocabulary. This namespace is included by default, as it is used for some of the basic metadata shown in the previous example.  We can see from this list:


```r
RNeXML:::nexml_namespaces
```

```
##                                              nex                                              xsi 
##                      "http://www.nexml.org/2009"      "http://www.w3.org/2001/XMLSchema-instance" 
##                                              xml                                             cdao 
##           "http://www.w3.org/XML/1998/namespace"                "http://purl.obolibrary.org/obo/" 
##                                              xsd                                               dc 
##              "http://www.w3.org/2001/XMLSchema#"               "http://purl.org/dc/elements/1.1/" 
##                                          dcterms                                            prism 
##                      "http://purl.org/dc/terms/" "http://prismstandard.org/namespaces/1.2/basic/" 
##                                               cc                                             ncbi 
##                 "http://creativecommons.org/ns#"          "http://www.ncbi.nlm.nih.gov/taxonomy#" 
##                                               tc 
##  "http://rs.tdwg.org/ontology/voc/TaxonConcept#"
```

This next block defines a resource (link), described by the `rel` attribute as a homepage, a term in the `foaf` vocabulary.  Because `foaf` is not a default namespace, we will have to provide its URL in the full definition below. 


```r
website <- meta(href = "http://carlboettiger.info", 
                rel = "foaf:homepage")
```

Here we create a history node using the `skos` namespace.  We can also add id values to any metadata element to make the element easier to reference externally: 


```r
  history <- meta(property = "skos:historyNote", 
                  content = "Mapped from the bird.orders data in the ape package using RNeXML",
                  id = "meta123")
```

For this kind of richer annotation, it is best to build up our NeXML object sequentially. First
we will add `bird.orders` phylogeny to a new phylogenetic object, and then we will add the metadata
elements created above to this object. Finally, we will write the object out as an XML file:


```r
  birds <- add_trees(bird.orders)
  birds <- add_meta(meta = list(history, modified, website),
                    namespaces = c(skos = "http://www.w3.org/2004/02/skos/core#",
                                   foaf = "http://xmlns.com/foaf/0.1/"),
                    nexml=birds)
  nexml_write(birds, 
              file = "example.xml")
```

```
## [1] "example.xml"
```

### Taxonomic identifiers

Add taxonomic identifier metadata to the OTU elements:
<!-- This block relies on a robust internet connection that can occassionally fail.  Also it's a bit slow, so don't run it. After all, this command is tested in the unit tests.-->


```r
nex <- add_trees(bird.orders)
nex <- taxize_nexml(nex)
```



## Working with character data

NeXML also provides a standard exchange format for handling character data.  The R platform is particularly popular in the context of phylogenetic comparative methods, which consider both a given phylogeny and a set of traits.  NeXML provides an ideal tool for handling this metadata.  

### Extracting character data

We can load the library, parse the NeXML file and extract both the characters and the phylogeny.  


```r
library(RNeXML)
nexml <- read.nexml(system.file("examples", "comp_analysis.xml", package="RNeXML"))
traits <- get_characters(nexml)
tree <- get_trees(nexml)
```

(Note that `get_characters` would return both discrete and continuous characters together in the same data.frame, but we use `get_characters_list` to get separate data.frames for the continuous `characters` block and the discrete `characters` block).  

We can then fire up `geiger` and fit, say, a Brownian motion model the continuous data and a Markov transition matrix to the discrete states:  

```r
library(geiger)
fitContinuous(tree, traits[1], ncores=1)
fitDiscrete(tree, traits[2], ncores=1)
```





---
title: "Handling Metadata in RNeXML"
author:
  - Carl Boettiger
  - Scott Chamberlain
  - Rutger Vos
  - Hilmar Lapp

output: html_vignette
---

<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{Handling Metadata in RNeXML}
-->




## Writing NeXML metadata

The `add_basic_meta()` function takes as input an existing `nexml` object
(like the other `add_` functions, if none is provided it will create one), and at the time
of this writing any of the following
parameters: `title`, `description`, `creator`, `pubdate`, `rights`, `publisher`,
`citation`.  Other metadata elements and corresponding parameters may
be added in the future.

Load the packages and data:


```r
library('RNeXML')
data(bird.orders)
```

Create an `nexml` object for the phylogeny `bird.orders` and add appropriate metadata:


```r
birds <- add_trees(bird.orders)
birds <- add_basic_meta(
  title = "Phylogeny of the Orders of Birds From Sibley and Ahlquist",

  description = "This data set describes the phylogenetic relationships of the
     orders of birds as reported by Sibley and Ahlquist (1990). Sibley
     and Ahlquist inferred this phylogeny from an extensive number of
     DNA/DNA hybridization experiments. The ``tapestry'' reported by
     these two authors (more than 1000 species out of the ca. 9000
     extant bird species) generated a lot of debates.

     The present tree is based on the relationships among orders. The
     branch lengths were calculated from the values of Delta T50H as
     found in Sibley and Ahlquist (1990, fig. 353).",

  citation = "Sibley, C. G. and Ahlquist, J. E. (1990) Phylogeny and
     classification of birds: a study in molecular evolution. New
     Haven: Yale University Press.",

  creator = "Sibley, C. G. and Ahlquist, J. E.",
	nexml=birds)
```

Instead of a literal string, citations can also be provided in R's
`bibentry` type, which is the one in which R package citations are obtained:


```r
birds <- add_basic_meta(citation = citation("ape"), nexml = birds)
```

## Taxonomic identifiers

RNeXML now uses `taxald` for fast look-up of taxonomic information.
If a unique match is found, a metadata annotation is added to the taxon
providing the unique URI to the taxonomic unit.



```r
birds <- taxize_nexml(birds, "NCBI")
```

If no perfect match is found, the user is warned to check for possible typographic
errors in the taxonomic labels provided.


## Custom metadata extensions

We can get a list of namespaces along with their prefixes from the `nexml` object: 


```r
prefixes <- get_namespaces(birds)
prefixes["dc"]
```

```
                                dc 
"http://purl.org/dc/elements/1.1/" 
```

We create a `meta` element containing this annotation using the `meta` function:


```r
modified <- meta(property = "prism:modificationDate", content = "2013-10-04")
```

We can add this annotation to our existing `birds` NeXML file using the
`add_meta()` function.  Because we do not specify a level, it is added to
the root node, referring to the NeXML file as a whole.


```r
birds <- add_meta(modified, birds) 
```

The built-in vocabularies are just the tip of the iceberg of established
vocabularies. Here we add an annotation from the `skos` namespace which
describes the history of where the data comes from:


```r
history <- meta(property = "skos:historyNote",
  content = "Mapped from the bird.orders data in the ape package using RNeXML")
```

Because `skos` is not in the current namespace list, we add it with a
url when adding this meta element.  We also specify that this annotation
be placed at the level of the `trees` sub-node in the NeXML file.


```r
birds <- add_meta(history, 
                birds, 
                level = "trees",
                namespaces = c(skos = "http://www.w3.org/2004/02/skos/core#"))
```


For finer control of the level at which a `meta` element is added,
we will manipulate the `nexml` R object directly using S4 sub-setting,
as shown in the supplement.


Much richer metadata annotation is possible. Later we illustrate how
metadata annotation can be used to extend the base NeXML format to
represent new forms of data while maintaining compatibility with any
NeXML parser. The `RNeXML` package can be easily extended to support
helper functions such as `taxize_nexml` to add additional metadata
without imposing a large burden on the user.


## Reading NeXML metadata

A call to the `nexml` object prints some metadata summarizing the data structure: 


```r
birds
```

```
A nexml object representing:
 	 1 phylogenetic tree block(s), where:
	   block 1 contains 1 phylogenetic tree(s) 
 	 0 character block(s),  
 	 23 taxonomic units in 1 block(s) 
  Taxa:	 Struthioniformes, Tinamiformes, Craciformes, Galliformes, Anseriformes, Turniciformes ...
  Metadata annotations: 
	18 at top level 
	0 in block 1 at otu level  

This data set describes the phylogenetic relationships of the
     orders of birds as reported by Sibley and Ahlquist (1990). Sibley
     and Ahlquist inferred this phylogeny from an extensive number of
     DNA/DNA hybridization experiments. The ``tapestry'' reported by
     these two authors (more than 1000 species out of the ca. 9000
     extant bird species) generated a lot of debates.

     The present tree is based on the relationships among orders. The
     branch lengths were calculated from the values of Delta T50H as
     found in Sibley and Ahlquist (1990, fig. 353).

Author(s): Sibley, C. G. and Ahlquist, J. E., cboettig

License: http://creativecommons.org/publicdomain/zero/1.0/

NeXML generated by RNeXML using schema version: 0.9 
Size: 358.3 Kb 
```

We can extract all metadata pertaining to the NeXML document as a whole
(annotations of the XML root node, `<nexml>`) with the command


```r
meta <- get_metadata(birds) 
```

This returns a data.frame of available metadata. We can see the kinds
of metadata recorded from the names:


```r
meta
```

```
                        property   datatype
1                       dc:title xsd:string
2                     dc:creator xsd:string
3                 dc:description xsd:string
4                     cc:license       <NA>
5  dcterms:bibliographicCitation xsd:string
6                     dc:creator xsd:string
7             dcterms:references       <NA>
8         prism:modificationDate xsd:string
9                   prism:volume xsd:string
10                  dc:publisher xsd:string
11         prism:publicationName xsd:string
12              prism:endingPage xsd:string
13            prism:startingPage xsd:string
14         prism:publicationDate xsd:string
15                      dc:title xsd:string
16                dc:contributor xsd:string
17                dc:contributor xsd:string
18 dcterms:bibliographicCitation xsd:string
                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                             content
1                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                          Phylogeny of the Orders of Birds From Sibley and Ahlquist
2                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                  Sibley, C. G. and Ahlquist, J. E.
3  This data set describes the phylogenetic relationships of the\n     orders of birds as reported by Sibley and Ahlquist (1990). Sibley\n     and Ahlquist inferred this phylogeny from an extensive number of\n     DNA/DNA hybridization experiments. The ``tapestry'' reported by\n     these two authors (more than 1000 species out of the ca. 9000\n     extant bird species) generated a lot of debates.\n\n     The present tree is based on the relationships among orders. The\n     branch lengths were calculated from the values of Delta T50H as\n     found in Sibley and Ahlquist (1990, fig. 353).
4                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               <NA>
5                                                                                                                                                                                                                                                                                                                                                                                                                                                      Sibley, C. G. and Ahlquist, J. E. (1990) Phylogeny and\n     classification of birds: a study in molecular evolution. New\n     Haven: Yale University Press.
6                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           cboettig
7                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               <NA>
8                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                         2013-10-04
9                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 35
10                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    Bioinformatics
11                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    Bioinformatics
12                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               528
13                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                               526
14                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                              2018
15                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                 ape 5.0: an environment for modern phylogenetics and evolutionary analyses in {R}
16                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        E. Paradis
17                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                        K. Schliep
18                                                                                                                                                                                                                                                                                                                                                                                                                                                                 Paradis E, Schliep K (2018). "ape 5.0: an environment for modern phylogenetics and evolutionary analyses in R."\n_Bioinformatics_, *35*, 526-528.
       xsi.type                                              href Meta meta
1   LiteralMeta                                              <NA> m544 <NA>
2   LiteralMeta                                              <NA> m545 <NA>
3   LiteralMeta                                              <NA> m546 <NA>
4  ResourceMeta http://creativecommons.org/publicdomain/zero/1.0/ m547 <NA>
5   LiteralMeta                                              <NA> m548 <NA>
6   LiteralMeta                                              <NA> m549 <NA>
7  ResourceMeta                                              <NA> m560 <NA>
8   LiteralMeta                                              <NA> m561 <NA>
9   LiteralMeta                                              <NA> <NA> m560
10  LiteralMeta                                              <NA> <NA> m560
11  LiteralMeta                                              <NA> <NA> m560
12  LiteralMeta                                              <NA> <NA> m560
13  LiteralMeta                                              <NA> <NA> m560
14  LiteralMeta                                              <NA> <NA> m560
15  LiteralMeta                                              <NA> <NA> m560
16  LiteralMeta                                              <NA> <NA> m560
17  LiteralMeta                                              <NA> <NA> m560
18  LiteralMeta                                              <NA> <NA> m560
```

We can also access a table of taxonomic metadata:



```r
get_taxa(birds)
```

```
     otu            label xsi.type otus
1  ou362 Struthioniformes       NA  os8
2  ou363     Tinamiformes       NA  os8
3  ou364      Craciformes       NA  os8
4  ou365      Galliformes       NA  os8
5  ou366     Anseriformes       NA  os8
6  ou367    Turniciformes       NA  os8
7  ou368       Piciformes       NA  os8
8  ou369    Galbuliformes       NA  os8
9  ou370   Bucerotiformes       NA  os8
10 ou371      Upupiformes       NA  os8
11 ou372    Trogoniformes       NA  os8
12 ou373    Coraciiformes       NA  os8
13 ou374      Coliiformes       NA  os8
14 ou375     Cuculiformes       NA  os8
15 ou376   Psittaciformes       NA  os8
16 ou377      Apodiformes       NA  os8
17 ou378   Trochiliformes       NA  os8
18 ou379  Musophagiformes       NA  os8
19 ou380     Strigiformes       NA  os8
20 ou381    Columbiformes       NA  os8
21 ou382       Gruiformes       NA  os8
22 ou383    Ciconiiformes       NA  os8
23 ou384    Passeriformes       NA  os8
```

Which returns text from the otu element labels, typically used to define
taxonomic names, rather than text from explicit meta elements.

We can also access metadata at a specific level (or use `level=all`
to extract all meta elements in a list).  Here we show only the first
few results:


```r
otu_meta <- get_metadata(birds, level="otus/otu")
otu_meta
```

```
# A tibble: 0 x 3
# … with 3 variables: otu <chr>, otus <chr>, Meta <chr>
```


## Merging metadata tables

We often want to combine metadata from multiple tables.  For instance, in this exercise we want to include the taxonomic identifier and id value for each species returned in the character table.  This helps us more precisely identify the species whose traits are described by the table.   



```r
library("RNeXML")
library("dplyr")
library("geiger")
```

To begin, let's generate a `NeXML` file using the tree and trait data from the `geiger` package's "primates" data:


```r
data("primates")
add_trees(primates$phy) %>% 
  add_characters(primates$dat, ., append=TRUE) %>% 
  taxize_nexml("itis") -> nex 
```

(Note that we've used `dplyr`'s cute pipe syntax, but unfortunately our `add_` methods take the `nexml` object as the _second_
argument instead of the first, so this isn't as elegant since we need the stupid `.` to show where the piped output should go...)

We now read in the three tables of interest.  Note that we tell `get_characters` to give us species labels as there own column, rather than as rownames.  The latter is the default only because this plays more nicely with the default format for character matrices that is expected by `geiger` and other phylogenetics packages, but is in general a silly choice for data manipulation. 


```r
otu_meta <- get_metadata(nex, "otus/otu")
taxa <- get_taxa(nex)
char <- get_characters(nex, rownames_as_col = TRUE)
```


We can take a peek at what the tables look like, just to orient ourselves:


```r
otu_meta
```

```
      property                                                                                  href     xsi.type   otu otus
1   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572992 ResourceMeta ou385  os9
2   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572857 ResourceMeta ou386  os9
3   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572935 ResourceMeta ou387  os9
4   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572936 ResourceMeta ou388  os9
5   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945168 ResourceMeta ou389  os9
6   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945164 ResourceMeta ou390  os9
7   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572939 ResourceMeta ou391  os9
8   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572940 ResourceMeta ou392  os9
9   tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572941 ResourceMeta ou393  os9
10  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572942 ResourceMeta ou394  os9
11  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944172 ResourceMeta ou395  os9
12  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572944 ResourceMeta ou396  os9
13  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572949 ResourceMeta ou397  os9
14  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945149 ResourceMeta ou398  os9
15  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572949 ResourceMeta ou399  os9
16  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572947 ResourceMeta ou400  os9
17  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572945 ResourceMeta ou401  os9
18  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572946 ResourceMeta ou402  os9
19  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572951 ResourceMeta ou403  os9
20  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572952 ResourceMeta ou404  os9
21  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572887 ResourceMeta ou405  os9
22  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572888 ResourceMeta ou406  os9
23  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572954 ResourceMeta ou407  os9
24  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572955 ResourceMeta ou408  os9
25  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572957 ResourceMeta ou409  os9
26  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572958 ResourceMeta ou410  os9
27  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572956 ResourceMeta ou411  os9
28  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572953 ResourceMeta ou412  os9
29  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572881 ResourceMeta ou413  os9
30  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572959 ResourceMeta ou414  os9
31  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572983 ResourceMeta ou415  os9
32  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572984 ResourceMeta ou416  os9
33  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025131 ResourceMeta ou417  os9
34  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025134 ResourceMeta ou418  os9
35  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025122 ResourceMeta ou419  os9
36  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025127 ResourceMeta ou420  os9
37  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025112 ResourceMeta ou421  os9
38  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025134 ResourceMeta ou422  os9
39  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025125 ResourceMeta ou423  os9
40  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025116 ResourceMeta ou424  os9
41  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025119 ResourceMeta ou425  os9
42  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025114 ResourceMeta ou426  os9
43  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025115 ResourceMeta ou427  os9
44  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572962 ResourceMeta ou428  os9
45  tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025106 ResourceMeta ou429  os9
46  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572909 ResourceMeta ou430  os9
47  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944125 ResourceMeta ou431  os9
48  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572911 ResourceMeta ou432  os9
49  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572912 ResourceMeta ou433  os9
50  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572913 ResourceMeta ou434  os9
51  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944138 ResourceMeta ou435  os9
52  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572915 ResourceMeta ou436  os9
53  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572916 ResourceMeta ou437  os9
54  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572918 ResourceMeta ou438  os9
55  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944123 ResourceMeta ou439  os9
56  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572975 ResourceMeta ou440  os9
57  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944156 ResourceMeta ou441  os9
58  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572977 ResourceMeta ou442  os9
59  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572978 ResourceMeta ou443  os9
60  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572993 ResourceMeta ou444  os9
61  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572994 ResourceMeta ou445  os9
62  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572995 ResourceMeta ou446  os9
63  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572996 ResourceMeta ou447  os9
64  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572997 ResourceMeta ou448  os9
65  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572998 ResourceMeta ou449  os9
66  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572999 ResourceMeta ou450  os9
67  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573000 ResourceMeta ou451  os9
68  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573001 ResourceMeta ou452  os9
69  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573002 ResourceMeta ou453  os9
70  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573003 ResourceMeta ou454  os9
71  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944219 ResourceMeta ou455  os9
72  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573005 ResourceMeta ou456  os9
73  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573006 ResourceMeta ou457  os9
74  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573007 ResourceMeta ou458  os9
75  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573008 ResourceMeta ou459  os9
76  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573009 ResourceMeta ou460  os9
77  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573010 ResourceMeta ou461  os9
78  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945227 ResourceMeta ou462  os9
79  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573014 ResourceMeta ou463  os9
80  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944221 ResourceMeta ou464  os9
81  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573011 ResourceMeta ou465  os9
82  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572858 ResourceMeta ou466  os9
83  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572859 ResourceMeta ou467  os9
84  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572985 ResourceMeta ou468  os9
85  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572986 ResourceMeta ou469  os9
86  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=552515 ResourceMeta ou470  os9
87  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573035 ResourceMeta ou471  os9
88  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573036 ResourceMeta ou472  os9
89  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573037 ResourceMeta ou473  os9
90  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573038 ResourceMeta ou474  os9
91  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572886 ResourceMeta ou475  os9
92  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573015 ResourceMeta ou476  os9
93  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572865 ResourceMeta ou477  os9
94  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572864 ResourceMeta ou478  os9
95  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572866 ResourceMeta ou479  os9
96  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572867 ResourceMeta ou480  os9
97  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572868 ResourceMeta ou481  os9
98  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572893 ResourceMeta ou482  os9
99  tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572894 ResourceMeta ou483  os9
100 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944100 ResourceMeta ou484  os9
101 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572896 ResourceMeta ou485  os9
102 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572897 ResourceMeta ou486  os9
103 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572898 ResourceMeta ou487  os9
104 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572899 ResourceMeta ou488  os9
105 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572900 ResourceMeta ou489  os9
106 tc:toTaxon http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025177 ResourceMeta ou490  os9
107 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573080 ResourceMeta ou491  os9
108 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572869 ResourceMeta ou492  os9
109 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572870 ResourceMeta ou493  os9
110 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944072 ResourceMeta ou494  os9
111 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=180092 ResourceMeta ou495  os9
112 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573069 ResourceMeta ou496  os9
113 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944286 ResourceMeta ou497  os9
114 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944292 ResourceMeta ou498  os9
115 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945338 ResourceMeta ou499  os9
116 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573073 ResourceMeta ou500  os9
117 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573074 ResourceMeta ou501  os9
118 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944290 ResourceMeta ou502  os9
119 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573076 ResourceMeta ou503  os9
120 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573077 ResourceMeta ou504  os9
121 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573078 ResourceMeta ou505  os9
122 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944293 ResourceMeta ou506  os9
123 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572882 ResourceMeta ou507  os9
124 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572961 ResourceMeta ou508  os9
125 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944187 ResourceMeta ou509  os9
126 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572872 ResourceMeta ou510  os9
127 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572921 ResourceMeta ou511  os9
128 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572920 ResourceMeta ou513  os9
129 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572922 ResourceMeta ou514  os9
130 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572874 ResourceMeta ou515  os9
131 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572875 ResourceMeta ou516  os9
132 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572876 ResourceMeta ou517  os9
133 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572877 ResourceMeta ou518  os9
134 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572878 ResourceMeta ou519  os9
135 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572879 ResourceMeta ou520  os9
136 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572880 ResourceMeta ou521  os9
137 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573016 ResourceMeta ou522  os9
138 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572889 ResourceMeta ou523  os9
139 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573017 ResourceMeta ou524  os9
140 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573018 ResourceMeta ou525  os9
141 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573019 ResourceMeta ou526  os9
142 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=180098 ResourceMeta ou527  os9
143 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=555659 ResourceMeta ou528  os9
144 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573020 ResourceMeta ou529  os9
145 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=180099 ResourceMeta ou530  os9
146 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573021 ResourceMeta ou531  os9
147 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573022 ResourceMeta ou532  os9
148 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573023 ResourceMeta ou533  os9
149 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573024 ResourceMeta ou534  os9
150 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573025 ResourceMeta ou535  os9
151 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573026 ResourceMeta ou536  os9
152 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573027 ResourceMeta ou537  os9
153 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573028 ResourceMeta ou538  os9
154 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573029 ResourceMeta ou539  os9
155 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573030 ResourceMeta ou540  os9
156 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573031 ResourceMeta ou541  os9
157 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944041 ResourceMeta ou542  os9
158 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572861 ResourceMeta ou543  os9
159 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572862 ResourceMeta ou544  os9
160 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573032 ResourceMeta ou545  os9
161 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944263 ResourceMeta ou546  os9
162 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573040 ResourceMeta ou547  os9
163 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572890 ResourceMeta ou548  os9
164 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572891 ResourceMeta ou549  os9
165 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572902 ResourceMeta ou550  os9
166 tc:toTaxon  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572903 ResourceMeta ou551  os9
    Meta
1   m563
2   m564
3   m565
4   m566
5   m567
6   m568
7   m569
8   m570
9   m571
10  m572
11  m573
12  m574
13  m575
14  m576
15  m577
16  m578
17  m579
18  m580
19  m581
20  m582
21  m583
22  m584
23  m585
24  m586
25  m587
26  m588
27  m589
28  m590
29  m591
30  m592
31  m593
32  m594
33  m595
34  m596
35  m597
36  m598
37  m599
38  m600
39  m601
40  m602
41  m603
42  m604
43  m605
44  m606
45  m607
46  m608
47  m609
48  m610
49  m611
50  m612
51  m613
52  m614
53  m615
54  m616
55  m617
56  m618
57  m619
58  m620
59  m621
60  m622
61  m623
62  m624
63  m625
64  m626
65  m627
66  m628
67  m629
68  m630
69  m631
70  m632
71  m633
72  m634
73  m635
74  m636
75  m637
76  m638
77  m639
78  m640
79  m641
80  m642
81  m643
82  m644
83  m645
84  m646
85  m647
86  m648
87  m649
88  m650
89  m651
90  m652
91  m653
92  m654
93  m655
94  m656
95  m657
96  m658
97  m659
98  m660
99  m661
100 m662
101 m663
102 m664
103 m665
104 m666
105 m667
106 m668
107 m669
108 m670
109 m671
110 m672
111 m673
112 m674
113 m675
114 m676
115 m677
116 m678
117 m679
118 m680
119 m681
120 m682
121 m683
122 m684
123 m685
124 m686
125 m687
126 m688
127 m689
128 m690
129 m691
130 m692
131 m693
132 m694
133 m695
134 m696
135 m697
136 m698
137 m699
138 m700
139 m701
140 m702
141 m703
142 m704
143 m705
144 m706
145 m707
146 m708
147 m709
148 m710
149 m711
150 m712
151 m713
152 m714
153 m715
154 m716
155 m717
156 m718
157 m719
158 m720
159 m721
160 m722
161 m723
162 m724
163 m725
164 m726
165 m727
166 m728
 [ reached 'max' / getOption("max.print") -- omitted 66 rows ]
```

```r
taxa
```

```
      otu                        label xsi.type otus
1   ou385  Allenopithecus_nigroviridis       NA  os9
2   ou386          Allocebus_trichotis       NA  os9
3   ou387            Alouatta_belzebul       NA  os9
4   ou388              Alouatta_caraya       NA  os9
5   ou389           Alouatta_coibensis       NA  os9
6   ou390               Alouatta_fusca       NA  os9
7   ou391            Alouatta_palliata       NA  os9
8   ou392               Alouatta_pigra       NA  os9
9   ou393                Alouatta_sara       NA  os9
10  ou394           Alouatta_seniculus       NA  os9
11  ou395                 Aotus_azarai       NA  os9
12  ou396              Aotus_brumbacki       NA  os9
13  ou397           Aotus_hershkovitzi       NA  os9
14  ou398              Aotus_infulatus       NA  os9
15  ou399              Aotus_lemurinus       NA  os9
16  ou400                Aotus_miconax       NA  os9
17  ou401              Aotus_nancymaae       NA  os9
18  ou402              Aotus_nigriceps       NA  os9
19  ou403            Aotus_trivirgatus       NA  os9
20  ou404             Aotus_vociferans       NA  os9
21  ou405            Arctocebus_aureus       NA  os9
22  ou406      Arctocebus_calabarensis       NA  os9
23  ou407             Ateles_belzebuth       NA  os9
24  ou408                Ateles_chamek       NA  os9
25  ou409             Ateles_fusciceps       NA  os9
26  ou410             Ateles_geoffroyi       NA  os9
27  ou411            Ateles_marginatus       NA  os9
28  ou412              Ateles_paniscus       NA  os9
29  ou413                Avahi_laniger       NA  os9
30  ou414      Brachyteles_arachnoides       NA  os9
31  ou415               Cacajao_calvus       NA  os9
32  ou416       Cacajao_melanocephalus       NA  os9
33  ou417          Callicebus_brunneus       NA  os9
34  ou418         Callicebus_caligatus       NA  os9
35  ou419       Callicebus_cinerascens       NA  os9
36  ou420           Callicebus_cupreus       NA  os9
37  ou421      Callicebus_donacophilus       NA  os9
38  ou422            Callicebus_dubius       NA  os9
39  ou423        Callicebus_hoffmannsi       NA  os9
40  ou424          Callicebus_modestus       NA  os9
41  ou425            Callicebus_moloch       NA  os9
42  ou426          Callicebus_oenanthe       NA  os9
43  ou427           Callicebus_olallae       NA  os9
44  ou428        Callicebus_personatus       NA  os9
45  ou429         Callicebus_torquatus       NA  os9
46  ou430            Callimico_goeldii       NA  os9
47  ou431         Callithrix_argentata       NA  os9
48  ou432            Callithrix_aurita       NA  os9
49  ou433         Callithrix_flaviceps       NA  os9
50  ou434         Callithrix_geoffroyi       NA  os9
51  ou435      Callithrix_humeralifera       NA  os9
52  ou436           Callithrix_jacchus       NA  os9
53  ou437            Callithrix_kuhlii       NA  os9
54  ou438       Callithrix_penicillata       NA  os9
55  ou439           Callithrix_pygmaea       NA  os9
56  ou440              Cebus_albifrons       NA  os9
57  ou441                 Cebus_apella       NA  os9
58  ou442              Cebus_capucinus       NA  os9
59  ou443              Cebus_olivaceus       NA  os9
60  ou444            Cercocebus_agilis       NA  os9
61  ou445         Cercocebus_galeritus       NA  os9
62  ou446         Cercocebus_torquatus       NA  os9
63  ou447       Cercopithecus_ascanius       NA  os9
64  ou448      Cercopithecus_campbelli       NA  os9
65  ou449         Cercopithecus_cephus       NA  os9
66  ou450          Cercopithecus_diana       NA  os9
67  ou451          Cercopithecus_dryas       NA  os9
68  ou452  Cercopithecus_erythrogaster       NA  os9
69  ou453     Cercopithecus_erythrotis       NA  os9
70  ou454        Cercopithecus_hamlyni       NA  os9
71  ou455        Cercopithecus_lhoesti       NA  os9
72  ou456          Cercopithecus_mitis       NA  os9
73  ou457           Cercopithecus_mona       NA  os9
74  ou458      Cercopithecus_neglectus       NA  os9
75  ou459      Cercopithecus_nictitans       NA  os9
76  ou460     Cercopithecus_petaurista       NA  os9
77  ou461       Cercopithecus_pogonias       NA  os9
78  ou462        Cercopithecus_preussi       NA  os9
79  ou463       Cercopithecus_sclateri       NA  os9
80  ou464        Cercopithecus_solatus       NA  os9
81  ou465          Cercopithecus_wolfi       NA  os9
82  ou466           Cheirogaleus_major       NA  os9
83  ou467          Cheirogaleus_medius       NA  os9
84  ou468         Chiropotes_albinasus       NA  os9
85  ou469           Chiropotes_satanas       NA  os9
86  ou470         Chlorocebus_aethiops       NA  os9
87  ou471           Colobus_angolensis       NA  os9
88  ou472              Colobus_guereza       NA  os9
89  ou473            Colobus_polykomos       NA  os9
90  ou474              Colobus_satanas       NA  os9
91  ou475 Daubentonia_madagascariensis       NA  os9
92  ou476           Erythrocebus_patas       NA  os9
93  ou477            Eulemur_coronatus       NA  os9
94  ou478               Eulemur_fulvus       NA  os9
95  ou479               Eulemur_macaco       NA  os9
96  ou480               Eulemur_mongoz       NA  os9
97  ou481          Eulemur_rubriventer       NA  os9
98  ou482         Euoticus_elegantulus       NA  os9
99  ou483            Euoticus_pallidus       NA  os9
100 ou484                Galago_alleni       NA  os9
101 ou485              Galago_gallarum       NA  os9
102 ou486             Galago_matschiei       NA  os9
103 ou487                Galago_moholi       NA  os9
104 ou488          Galago_senegalensis       NA  os9
105 ou489          Galagoides_demidoff       NA  os9
106 ou490      Galagoides_zanzibaricus       NA  os9
107 ou491              Gorilla_gorilla       NA  os9
108 ou492             Hapalemur_aureus       NA  os9
109 ou493            Hapalemur_griseus       NA  os9
110 ou494              Hapalemur_simus       NA  os9
111 ou495                 Homo_sapiens       NA  os9
112 ou496             Hylobates_agilis       NA  os9
113 ou497           Hylobates_concolor       NA  os9
114 ou498         Hylobates_gabriellae       NA  os9
115 ou499            Hylobates_hoolock       NA  os9
116 ou500            Hylobates_klossii       NA  os9
117 ou501                Hylobates_lar       NA  os9
118 ou502         Hylobates_leucogenys       NA  os9
119 ou503             Hylobates_moloch       NA  os9
120 ou504           Hylobates_muelleri       NA  os9
121 ou505           Hylobates_pileatus       NA  os9
122 ou506        Hylobates_syndactylus       NA  os9
123 ou507                  Indri_indri       NA  os9
124 ou508         Lagothrix_flavicauda       NA  os9
125 ou509         Lagothrix_lagotricha       NA  os9
126 ou510                  Lemur_catta       NA  os9
127 ou511      Leontopithecus_caissara       NA  os9
128 ou512    Leontopithecus_chrysomela       NA  os9
129 ou513   Leontopithecus_chrysopygus       NA  os9
130 ou514       Leontopithecus_rosalia       NA  os9
131 ou515           Lepilemur_dorsalis       NA  os9
132 ou516           Lepilemur_edwardsi       NA  os9
133 ou517           Lepilemur_leucopus       NA  os9
134 ou518           Lepilemur_microdon       NA  os9
135 ou519         Lepilemur_mustelinus       NA  os9
136 ou520       Lepilemur_ruficaudatus       NA  os9
137 ou521    Lepilemur_septentrionalis       NA  os9
138 ou522          Lophocebus_albigena       NA  os9
139 ou523            Loris_tardigradus       NA  os9
140 ou524             Macaca_arctoides       NA  os9
141 ou525            Macaca_assamensis       NA  os9
142 ou526              Macaca_cyclopis       NA  os9
143 ou527          Macaca_fascicularis       NA  os9
144 ou528               Macaca_fuscata       NA  os9
145 ou529                 Macaca_maura       NA  os9
146 ou530               Macaca_mulatta       NA  os9
147 ou531            Macaca_nemestrina       NA  os9
148 ou532                 Macaca_nigra       NA  os9
149 ou533              Macaca_ochreata       NA  os9
150 ou534               Macaca_radiata       NA  os9
151 ou535               Macaca_silenus       NA  os9
152 ou536                Macaca_sinica       NA  os9
153 ou537              Macaca_sylvanus       NA  os9
154 ou538             Macaca_thibetana       NA  os9
155 ou539              Macaca_tonkeana       NA  os9
156 ou540       Mandrillus_leucophaeus       NA  os9
157 ou541            Mandrillus_sphinx       NA  os9
158 ou542         Microcebus_coquereli       NA  os9
159 ou543           Microcebus_murinus       NA  os9
160 ou544             Microcebus_rufus       NA  os9
161 ou545         Miopithecus_talapoin       NA  os9
162 ou546             Nasalis_concolor       NA  os9
163 ou547             Nasalis_larvatus       NA  os9
164 ou548           Nycticebus_coucang       NA  os9
165 ou549          Nycticebus_pygmaeus       NA  os9
166 ou550      Otolemur_crassicaudatus       NA  os9
167 ou551           Otolemur_garnettii       NA  os9
168 ou552                 Pan_paniscus       NA  os9
169 ou553              Pan_troglodytes       NA  os9
170 ou554              Papio_hamadryas       NA  os9
171 ou555           Perodicticus_potto       NA  os9
172 ou556              Phaner_furcifer       NA  os9
173 ou557       Pithecia_aequatorialis       NA  os9
174 ou558            Pithecia_albicans       NA  os9
175 ou559            Pithecia_irrorata       NA  os9
176 ou560            Pithecia_monachus       NA  os9
177 ou561            Pithecia_pithecia       NA  os9
178 ou562               Pongo_pygmaeus       NA  os9
179 ou563             Presbytis_comata       NA  os9
180 ou564          Presbytis_femoralis       NA  os9
181 ou565           Presbytis_frontata       NA  os9
182 ou566              Presbytis_hosei       NA  os9
183 ou567         Presbytis_melalophos       NA  os9
184 ou568         Presbytis_potenziani       NA  os9
185 ou569          Presbytis_rubicunda       NA  os9
186 ou570            Presbytis_thomasi       NA  os9
187 ou571            Procolobus_badius       NA  os9
188 ou572         Procolobus_pennantii       NA  os9
189 ou573           Procolobus_preussi       NA  os9
190 ou574      Procolobus_rufomitratus       NA  os9
191 ou575             Procolobus_verus       NA  os9
192 ou576          Propithecus_diadema       NA  os9
193 ou577      Propithecus_tattersalli       NA  os9
194 ou578        Propithecus_verreauxi       NA  os9
195 ou579          Pygathrix_avunculus       NA  os9
196 ou580              Pygathrix_bieti       NA  os9
197 ou581           Pygathrix_brelichi       NA  os9
198 ou582            Pygathrix_nemaeus       NA  os9
199 ou583          Pygathrix_roxellana       NA  os9
200 ou584             Saguinus_bicolor       NA  os9
201 ou585         Saguinus_fuscicollis       NA  os9
202 ou586           Saguinus_geoffroyi       NA  os9
203 ou587           Saguinus_imperator       NA  os9
204 ou588             Saguinus_inustus       NA  os9
205 ou589            Saguinus_labiatus       NA  os9
206 ou590            Saguinus_leucopus       NA  os9
207 ou591               Saguinus_midas       NA  os9
208 ou592              Saguinus_mystax       NA  os9
209 ou593         Saguinus_nigricollis       NA  os9
210 ou594             Saguinus_oedipus       NA  os9
211 ou595         Saguinus_tripartitus       NA  os9
212 ou596          Saimiri_boliviensis       NA  os9
213 ou597            Saimiri_oerstedii       NA  os9
214 ou598             Saimiri_sciureus       NA  os9
215 ou599                Saimiri_ustus       NA  os9
216 ou600           Saimiri_vanzolinii       NA  os9
217 ou601       Semnopithecus_entellus       NA  os9
218 ou602             Tarsius_bancanus       NA  os9
219 ou603               Tarsius_dianae       NA  os9
220 ou604              Tarsius_pumilus       NA  os9
221 ou605             Tarsius_spectrum       NA  os9
222 ou606             Tarsius_syrichta       NA  os9
223 ou607         Theropithecus_gelada       NA  os9
224 ou608       Trachypithecus_auratus       NA  os9
225 ou609     Trachypithecus_cristatus       NA  os9
226 ou610     Trachypithecus_francoisi       NA  os9
227 ou611          Trachypithecus_geei       NA  os9
228 ou612        Trachypithecus_johnii       NA  os9
229 ou613      Trachypithecus_obscurus       NA  os9
230 ou614       Trachypithecus_phayrei       NA  os9
231 ou615      Trachypithecus_pileatus       NA  os9
232 ou616       Trachypithecus_vetulus       NA  os9
233 ou617            Varecia_variegata       NA  os9
```

```r
head(char)
```

```
                         taxa        x
1 Allenopithecus_nigroviridis 8.465900
2         Allocebus_trichotis 4.368181
3           Alouatta_belzebul 8.729074
4             Alouatta_caraya 8.628735
5          Alouatta_coibensis 8.764053
6              Alouatta_fusca 8.554489
```

Now that we have nice `data.frame` objects for all our data, it's easy to join them into the desired table with a few obvious `dplyr` commands:


```r
taxa %>% 
  left_join(char, by = c("label" = "taxa")) %>% 
  left_join(otu_meta, by = "otu") %>%
  select(otu, label, x, href)
```

```
      otu                        label         x
1   ou385  Allenopithecus_nigroviridis  8.465900
2   ou386          Allocebus_trichotis  4.368181
3   ou387            Alouatta_belzebul  8.729074
4   ou388              Alouatta_caraya  8.628735
5   ou389           Alouatta_coibensis  8.764053
6   ou390               Alouatta_fusca  8.554489
7   ou391            Alouatta_palliata  8.791790
8   ou392               Alouatta_pigra  8.881836
9   ou393                Alouatta_sara  8.796339
10  ou394           Alouatta_seniculus  8.767173
11  ou395                 Aotus_azarai  7.127694
12  ou396              Aotus_brumbacki  6.401917
13  ou397           Aotus_hershkovitzi  6.684612
14  ou398              Aotus_infulatus  7.122867
15  ou399              Aotus_lemurinus  6.781058
16  ou400                Aotus_miconax  6.684612
17  ou401              Aotus_nancymaae  6.678342
18  ou402              Aotus_nigriceps  6.966024
19  ou403            Aotus_trivirgatus  6.815640
20  ou404             Aotus_vociferans  6.771936
21  ou405            Arctocebus_aureus  5.459586
22  ou406      Arctocebus_calabarensis  5.552960
23  ou407             Ateles_belzebuth  8.812843
24  ou408                Ateles_chamek  8.872067
25  ou409             Ateles_fusciceps  9.112728
26  ou410             Ateles_geoffroyi  8.937218
27  ou411            Ateles_marginatus  8.738735
28  ou412              Ateles_paniscus  9.072227
29  ou413                Avahi_laniger  6.927558
30  ou414      Brachyteles_arachnoides  9.268609
31  ou415               Cacajao_calvus  8.137396
32  ou416       Cacajao_melanocephalus  8.051978
33  ou417          Callicebus_brunneus  6.747587
34  ou418         Callicebus_caligatus  6.899723
35  ou419       Callicebus_cinerascens  6.899723
36  ou420           Callicebus_cupreus  7.021084
37  ou421      Callicebus_donacophilus  6.801283
38  ou422            Callicebus_dubius  6.899723
39  ou423        Callicebus_hoffmannsi  6.975414
40  ou424          Callicebus_modestus  6.899723
41  ou425            Callicebus_moloch  6.872128
42  ou426          Callicebus_oenanthe  6.899723
43  ou427           Callicebus_olallae  6.899723
44  ou428        Callicebus_personatus  7.237059
45  ou429         Callicebus_torquatus  7.106606
46  ou430            Callimico_goeldii  6.324359
47  ou431         Callithrix_argentata  5.910797
48  ou432            Callithrix_aurita  5.961005
49  ou433         Callithrix_flaviceps  5.926926
50  ou434         Callithrix_geoffroyi  5.834811
51  ou435      Callithrix_humeralifera  5.926926
52  ou436           Callithrix_jacchus  5.673323
53  ou437            Callithrix_kuhlii  5.926926
54  ou438       Callithrix_penicillata  5.828946
55  ou439           Callithrix_pygmaea  4.828314
56  ou440              Cebus_albifrons  7.832014
57  ou441                 Cebus_apella  7.922986
58  ou442              Cebus_capucinus  8.009695
59  ou443              Cebus_olivaceus  7.937375
60  ou444            Cercocebus_agilis  8.870663
61  ou445         Cercocebus_galeritus  8.865029
62  ou446         Cercocebus_torquatus  8.894259
63  ou447       Cercopithecus_ascanius  8.171882
64  ou448      Cercopithecus_campbelli  8.199739
65  ou449         Cercopithecus_cephus  8.143227
66  ou450          Cercopithecus_diana  8.382518
67  ou451          Cercopithecus_dryas  7.930206
68  ou452  Cercopithecus_erythrogaster  8.143227
69  ou453     Cercopithecus_erythrotis  8.086410
70  ou454        Cercopithecus_hamlyni  8.438150
71  ou455        Cercopithecus_lhoesti  8.579229
72  ou456          Cercopithecus_mitis  8.525161
73  ou457           Cercopithecus_mona  8.289037
74  ou458      Cercopithecus_neglectus  8.579229
75  ou459      Cercopithecus_nictitans  8.567886
76  ou460     Cercopithecus_petaurista  8.083329
77  ou461       Cercopithecus_pogonias  8.183118
78  ou462        Cercopithecus_preussi  8.544808
79  ou463       Cercopithecus_sclateri  8.029433
80  ou464        Cercopithecus_solatus  8.567886
81  ou465          Cercopithecus_wolfi  8.089482
82  ou466           Cheirogaleus_major  6.104793
83  ou467          Cheirogaleus_medius  5.283204
84  ou468         Chiropotes_albinasus  7.937375
85  ou469           Chiropotes_satanas  7.996317
86  ou470         Chlorocebus_aethiops  8.301522
87  ou471           Colobus_angolensis  9.103868
88  ou472              Colobus_guereza  9.210340
89  ou473            Colobus_polykomos  9.060680
90  ou474              Colobus_satanas  9.121509
91  ou475 Daubentonia_madagascariensis  7.919356
92  ou476           Erythrocebus_patas  8.988446
93  ou477            Eulemur_coronatus  7.438384
94  ou478               Eulemur_fulvus  7.779049
95  ou479               Eulemur_macaco  7.816014
96  ou480               Eulemur_mongoz  7.478735
97  ou481          Eulemur_rubriventer  7.620705
98  ou482         Euoticus_elegantulus  5.690359
99  ou483            Euoticus_pallidus  5.627621
100 ou484                Galago_alleni  5.521461
101 ou485              Galago_gallarum  5.365976
102 ou486             Galago_matschiei  5.267858
103 ou487                Galago_moholi  5.375278
104 ou488          Galago_senegalensis  5.594711
105 ou489          Galagoides_demidoff  4.206184
106 ou490      Galagoides_zanzibaricus  4.997212
107 ou491              Gorilla_gorilla 11.643954
108 ou492             Hapalemur_aureus  7.358831
109 ou493            Hapalemur_griseus  6.862758
110 ou494              Hapalemur_simus  7.620705
111 ou495                 Homo_sapiens 10.980195
112 ou496             Hylobates_agilis  8.675905
113 ou497           Hylobates_concolor  8.768730
114 ou498         Hylobates_gabriellae  9.047821
115 ou499            Hylobates_hoolock  8.809863
116 ou500            Hylobates_klossii  8.674197
117 ou501                Hylobates_lar  8.635865
118 ou502         Hylobates_leucogenys  8.898366
119 ou503             Hylobates_moloch  8.677610
120 ou504           Hylobates_muelleri  8.687779
121 ou505           Hylobates_pileatus  8.630522
122 ou506        Hylobates_syndactylus  9.296518
123 ou507                  Indri_indri  9.065315
124 ou508         Lagothrix_flavicauda  9.020390
125 ou509         Lagothrix_lagotricha  8.743532
126 ou510                  Lemur_catta  7.878534
127 ou511      Leontopithecus_caissara  6.405228
128 ou512    Leontopithecus_chrysomela  6.352629
129 ou513   Leontopithecus_chrysopygus  6.486161
130 ou514       Leontopithecus_rosalia  6.385194
131 ou515           Lepilemur_dorsalis  6.232448
132 ou516           Lepilemur_edwardsi  6.711740
133 ou517           Lepilemur_leucopus  6.396930
134 ou518           Lepilemur_microdon  6.863803
135 ou519         Lepilemur_mustelinus  6.508769
136 ou520       Lepilemur_ruficaudatus  6.637258
137 ou521    Lepilemur_septentrionalis  6.638568
138 ou522          Lophocebus_albigena  8.903815
139 ou523            Loris_tardigradus  5.517453
140 ou524             Macaca_arctoides  9.148465
141 ou525            Macaca_assamensis  9.054855
142 ou526              Macaca_cyclopis  8.658693
143 ou527          Macaca_fascicularis  8.431635
144 ou528               Macaca_fuscata  9.220291
145 ou529                 Macaca_maura  8.894259
146 ou530               Macaca_mulatta  8.771835
147 ou531            Macaca_nemestrina  8.970813
148 ou532                 Macaca_nigra  8.906529
149 ou533              Macaca_ochreata  7.919356
150 ou534               Macaca_radiata  8.517193
151 ou535               Macaca_silenus  8.699515
152 ou536                Macaca_sinica  8.446771
153 ou537              Macaca_sylvanus  9.350102
154 ou538             Macaca_thibetana  9.268609
155 ou539              Macaca_tonkeana  9.220291
156 ou540       Mandrillus_leucophaeus  9.560997
157 ou541            Mandrillus_sphinx  9.729134
158 ou542         Microcebus_coquereli  5.793014
159 ou543           Microcebus_murinus  4.198705
160 ou544             Microcebus_rufus  3.881564
161 ou545         Miopithecus_talapoin  7.130899
162 ou546             Nasalis_concolor  8.979354
163 ou547             Nasalis_larvatus  9.417355
164 ou548           Nycticebus_coucang  6.878326
165 ou549          Nycticebus_pygmaeus  5.840642
166 ou550      Otolemur_crassicaudatus  7.130899
167 ou551           Otolemur_garnettii  6.701960
168 ou552                 Pan_paniscus 10.468801
169 ou553              Pan_troglodytes 10.716638
170 ou554              Papio_hamadryas  9.735069
171 ou555           Perodicticus_potto  6.984716
172 ou556              Phaner_furcifer  6.016157
173 ou557       Pithecia_aequatorialis  7.774856
174 ou558            Pithecia_albicans  7.937375
175 ou559            Pithecia_irrorata  7.749322
176 ou560            Pithecia_monachus  7.654443
177 ou561            Pithecia_pithecia  7.420579
178 ou562               Pongo_pygmaeus 10.860920
179 ou563             Presbytis_comata  9.039552
180 ou564          Presbytis_femoralis  8.797851
181 ou565           Presbytis_frontata  8.853665
182 ou566              Presbytis_hosei  8.725832
183 ou567         Presbytis_melalophos  8.748305
184 ou568         Presbytis_potenziani  8.774931
185 ou569          Presbytis_rubicunda  8.773385
186 ou570            Presbytis_thomasi  8.765615
187 ou571            Procolobus_badius  8.808369
188 ou572         Procolobus_pennantii  9.122601
189 ou573           Procolobus_preussi  9.094930
190 ou574      Procolobus_rufomitratus  8.995909
191 ou575             Procolobus_verus  8.289037
192 ou576          Propithecus_diadema  8.794825
193 ou577      Propithecus_tattersalli  8.171882
194 ou578        Propithecus_verreauxi  8.202482
195 ou579          Pygathrix_avunculus  9.164296
196 ou580              Pygathrix_bieti  9.118225
197 ou581           Pygathrix_brelichi  9.305651
198 ou582            Pygathrix_nemaeus  9.417355
199 ou583          Pygathrix_roxellana  9.510445
200 ou584             Saguinus_bicolor  6.142037
201 ou585         Saguinus_fuscicollis  5.976351
202 ou586           Saguinus_geoffroyi  6.200509
203 ou587           Saguinus_imperator  6.016157
204 ou588             Saguinus_inustus  6.688355
205 ou589            Saguinus_labiatus  6.234411
206 ou590            Saguinus_leucopus  6.124683
207 ou591               Saguinus_midas  6.293419
208 ou592              Saguinus_mystax  6.326149
209 ou593         Saguinus_nigricollis  6.109248
210 ou594             Saguinus_oedipus  6.139885
211 ou595         Saguinus_tripartitus  5.953243
212 ou596          Saimiri_boliviensis  6.690842
213 ou597            Saimiri_oerstedii  6.570883
214 ou598             Saimiri_sciureus  6.620073
215 ou599                Saimiri_ustus  6.791221
216 ou600           Saimiri_vanzolinii  6.664409
217 ou601       Semnopithecus_entellus  9.441452
218 ou602             Tarsius_bancanus  4.736198
219 ou603               Tarsius_dianae  4.709530
220 ou604              Tarsius_pumilus  4.808111
221 ou605             Tarsius_spectrum  5.111988
222 ou606             Tarsius_syrichta  4.753590
223 ou607         Theropithecus_gelada  9.680344
224 ou608       Trachypithecus_auratus  9.181941
225 ou609     Trachypithecus_cristatus  8.872067
226 ou610     Trachypithecus_francoisi  9.008224
227 ou611          Trachypithecus_geei  9.031214
228 ou612        Trachypithecus_johnii  9.268609
229 ou613      Trachypithecus_obscurus  8.870663
230 ou614       Trachypithecus_phayrei  8.948976
231 ou615      Trachypithecus_pileatus  9.323669
232 ou616       Trachypithecus_vetulus  8.969542
233 ou617            Varecia_variegata  8.261010
                                                                                     href
1    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572992
2    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572857
3    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572935
4    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572936
5    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945168
6    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945164
7    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572939
8    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572940
9    http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572941
10   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572942
11   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944172
12   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572944
13   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572949
14   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945149
15   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572949
16   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572947
17   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572945
18   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572946
19   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572951
20   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572952
21   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572887
22   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572888
23   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572954
24   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572955
25   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572957
26   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572958
27   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572956
28   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572953
29   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572881
30   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572959
31   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572983
32   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572984
33  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025131
34  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025134
35  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025122
36  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025127
37  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025112
38  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025134
39  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025125
40  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025116
41  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025119
42  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025114
43  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025115
44   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572962
45  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025106
46   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572909
47   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944125
48   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572911
49   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572912
50   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572913
51   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944138
52   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572915
53   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572916
54   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572918
55   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944123
56   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572975
57   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944156
58   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572977
59   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572978
60   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572993
61   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572994
62   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572995
63   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572996
64   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572997
65   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572998
66   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572999
67   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573000
68   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573001
69   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573002
70   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573003
71   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944219
72   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573005
73   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573006
74   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573007
75   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573008
76   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573009
77   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573010
78   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945227
79   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573014
80   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944221
81   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573011
82   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572858
83   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572859
84   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572985
85   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572986
86   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=552515
87   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573035
88   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573036
89   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573037
90   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573038
91   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572886
92   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573015
93   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572865
94   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572864
95   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572866
96   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572867
97   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572868
98   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572893
99   http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572894
100  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944100
101  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572896
102  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572897
103  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572898
104  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572899
105  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572900
106 http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025177
107  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573080
108  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572869
109  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572870
110  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944072
111  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=180092
112  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573069
113  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944286
114  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944292
115  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945338
116  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573073
117  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573074
118  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944290
119  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573076
120  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573077
121  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573078
122  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944293
123  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572882
124  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572961
125  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944187
126  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572872
127  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572921
128                                                                                  <NA>
129  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572920
130  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572922
131  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572874
132  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572875
133  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572876
134  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572877
135  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572878
136  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572879
137  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572880
138  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573016
139  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572889
140  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573017
141  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573018
142  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573019
143  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=180098
144  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=555659
145  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573020
146  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=180099
147  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573021
148  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573022
149  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573023
150  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573024
151  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573025
152  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573026
153  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573027
154  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573028
155  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573029
156  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573030
157  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573031
158  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944041
159  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572861
160  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572862
161  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573032
162  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944263
163  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573040
164  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572890
165  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572891
166  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572902
167  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572903
168  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573081
169  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573082
170  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573033
171  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572892
172  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572863
173  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572987
174  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572988
175  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572989
176  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572990
177  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572991
178  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573083
179  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573041
180  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573042
181  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573043
182  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573044
183  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573045
184  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573046
185  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573047
186  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573048
187  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944230
188  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944235
189  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944238
190  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944239
191  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573053
192  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572883
193  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572884
194  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572885
195  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944258
196  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944259
197  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944260
198  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573057
199  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945308
200  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572923
201  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=998616
202  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572925
203  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572926
204  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572927
205  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572928
206  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572929
207  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572930
208  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572931
209  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=998617
210  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572933
211 http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=1025098
212  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572979
213  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572980
214  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=180095
215  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572981
216  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572982
217  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573059
218  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945103
219  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944117
220  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572906
221  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944115
222  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=945107
223  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573034
224  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573060
225  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573061
226  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573062
227  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573063
228  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944270
229  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573065
230  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573066
231  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=573067
232  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=944269
233  http://www.itis.gov/servlet/SingleRpt/SingleRpt?search_topic=TSN&search_value=572873
```

Because these are all from the same otus block anyway, we haven't selected that column, but were it of interest it is also available in the taxa table.




---
title: "Extending NeXML: an example based on simmap"
author:
  - Carl Boettiger
  - Scott Chamberlain
  - Rutger Vos
  - Hilmar Lapp
output: html_vignette
bibliography: references.bib
---

<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{Extending NeXML to simmap}
-->





## Extending the NeXML standard through metadata annotation.

Here we illustrate this process using the example of stochastic character
mapping [@Huelsenbeck_2003]. A stochastic character map is simply
an annotation of the branches on a phylogeny, assigning each section
of each branch to a particular "state" (typically of a morphological
characteristic).

@Bollback_2006 provides a widely used stand-alone software implementation
of this method in the software `simmap`, which modified the standard
Newick tree format to express this additional information. This can
break compatibility with other software, and creates a format that
cannot be interpreted without additional information describing this
convention.  By contrast, the NeXML extension is not only backwards
compatible but contains a precise and machine-readable description of
what it is encoding.


In this example, we illustrate how the additional information required
to define a stochastic character mapping (a `simmap` mapping) in NeXML. 

@Revell_2012 describes the `phytools` package for R, which includes
utilities for reading, manipulating, and writing `simmap` files in R.
In this example, we also show how to define `RNeXML` functions that 
map the R representation used by Revell (an extension of the `ape` class)
into the NeXML extension we have defined by using `RNeXML` functions.

Since a stochastic character map simply assigns different states to
parts of a branch (or edge) on the phylogenetic tree, we can create
a NeXML representation by annotating the `edge` elements with appropriate
`meta` elements.  These elements need to describe the character state
being assigned and the duration (in terms of branch-length) that the edge 
spends in that state (Stochastic character maps are specific to time-calibrated
or ultrametric trees).  

NeXML already defines the `characters` element to handle discrete character traits (`nex:char`)
and the states they can assume (`nex:state`).  We will thus reuse the `characters` element for
this purpose, referring to both the character trait and the states by the ids assigned to them
in that element.  (NeXML's convention of referring to everything by id permits a single canonical
definition of each term, making it clear where additional annotation belongs).  For each edge, we 
need to indicate:

- That our annotation contains a stochastic character mapping reconstruction
- Since many reconstructions are possible for a single edge, we give each reconstruction an id
- We indicate for which character trait we are defining the reconstruction 
- We then indicate which states the character assumes on that edge. 
  For each state realized on the edge, that involves stating: 
    + the state assignment
    + the duration (length of time) for which the edge spends in the given state
    + the order in which the state changes happen (Though we could just assume 
      state transitions are listed chronologically, NeXML suggests making all 
      data explicit, rather than relying on the structure of the data file to
      convey information).  

Thus the annotation for an edge that switches from state `s2` to state 
`s1` of character `cr1` would be constructed like this:


```r
 m <- meta("simmap:reconstructions", children = c(
        meta("simmap:reconstruction", children = c(

          meta("simmap:char", "cr1"),
          meta("simmap:stateChange", children = c(
            meta("simmap:order", 1),
            meta("simmap:length", "0.2030"),
            meta("simmap:state", "s2"))),
          
          meta("simmap:char", "cr1"),
          meta("simmap:stateChange", children = c(
            meta("simmap:order", 2),
            meta("simmap:length", "0.0022"),
            meta("simmap:state", "s1")))
          ))))
```

Of course writing out such a definition manually becomes tedious quickly. Because
these are just R commands, we can easily define a function that can loop over an
assignment like this for each edge, extracting the appropriate order, length and
state from an existing R object such as that provided in the `phytools` package.  
Likewise, it is straightforward to define a function that reads this data using
the `RNeXML` utilities and converts it back to the `phytools` package. The full
implementation of this mapping can be seen in the `simmap_to_nexml()` and the
`nexml_to_simmap()` functions provided in the `RNeXML` package.  

As the code indicates, the key step is simply to define the data in meta elements. In 
so doing, we have defined a custom namespace, `simmap`, to hold our variables.  This
allows us to provide a URL with more detailed descriptions of what each of these 
elements mean:


```r
nex <- add_namespaces(c(simmap = "https://github.com/ropensci/RNeXML/tree/master/inst/simmap.md"))
```

At that URL we have posted a simple description of each term. 

Using this convention we can generate NeXML files containing `simmap`
data, read those files into R, and convert them back into the `phytools`
package format. These simple functions serve as further illustration of
how `RNeXML` can be used to extend the NeXML standard.  We illustrate
their use briefly here, starting with  loading a `nexml` object containing
a `simmap` reconstruction into R:



```r
f <- system.file("examples", "simmap_ex.xml", package = "RNeXML")
simmap_ex <- read.nexml(f)
```

The `get_trees()` function can be used to return an `ape::phylo` tree as
usual.  `RNeXML` automatically detects the `simmap` reconstruction data
and returns includes this in a `maps` element of the `ape::phylo` object,
for use with other `phytools` functions.


```r
phy <- nexml_to_simmap(simmap_ex)
```

We can then use various functions from `phytools` designed for `simmap`
objects [@Revell_2012], such as the plotting function:


```r
library("phytools")
plotSimmap(phy)
```

```
no colors provided. using the following legend:
       A        B        C 
 "black"    "red" "green3" 
```

![Stochastic character mapping on a phylogeny, as generated by the phytools package after parsing the simmap-extended NeXML.](simmap-Figure1-1.png)

Likewise, we can convert the object back in the NeXML format and write
it out to file to be read by other users. 


```r
nex <- simmap_to_nexml(phy) 
nexml_write(nex, "simmap.xml")
```

```
[1] "simmap.xml"
```

Though other NeXML parsers (for instance, for Perl or Python) have
not been written explicitly to express `simmap` data, those parsers will
nonetheless be able to successfully parse this file and expose the `simmap`
data to the user.



---
title: "The nexml S4 Object"
author:
  - Carl Boettiger
  - Scott Chamberlain
  - Rutger Vos
  - Hilmar Lapp
output: html_vignette
---

<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{The nexml S4 Object}
-->







## Understanding the `nexml` S4 object

The `RNeXML` package provides many convenient functions to add and extract
information from `nexml` objects in the R environment without requiring
the reader to understand the details of the NeXML data structure and
making it less likely that a user will generate invalid NeXML syntax
that could not be read by other parsers. The `nexml` object we have been using
in all of the examples is built on R's S4 mechanism. Advanced users may
sometimes prefer to interact with the data structure more directly using 
R's S4 class mechanism and subsetting methods. Many R users are more familiar
with the S3 class mechanism (such as in the `ape` package phylo objects)
rather than the S4 class mechanism used in phylogenetics packages such as 
`ouch` and `phylobase`. The `phylobase` vignette provides an excellent introduction
to these data structures.  Users already familiar with subsetting lists and other
S3 objects in R are likely familiar with the use of the `$` operator, such as
`phy$edge`. S4 objects simply use an `@` operator instead (but cannot be subset
using numeric arguments such as `phy[[1]]` or named arguments such as phy[["edge"]]).  


The `nexml` object is an S4 object, as are all of its components (slots).  Its 
hierarchical structure corresponds exactly with the XML tree of a NeXML file, with 
the single exception that both XML attributes and children are represented as slots.  
S4 objects have constructor functions to initialize them.  We create a new `nexml` 
object with the command:


```r
nex <- new("nexml")
```

We can see a list of slots contained in this object with


```r
slotNames(nex)
```

```
 [1] "version"            "generator"          "xsi:schemaLocation" "namespaces"         "otus"              
 [6] "trees"              "characters"         "meta"               "about"              "xsi:type"          
```

Some of these slots have already been populated for us, for instance, the schema version and default namespaces:


```r
nex@version
```

```
[1] "0.9"
```

```r
nex@namespaces
```

```
                                             nex                                              xsi 
                     "http://www.nexml.org/2009"      "http://www.w3.org/2001/XMLSchema-instance" 
                                             xml                                             cdao 
          "http://www.w3.org/XML/1998/namespace"                "http://purl.obolibrary.org/obo/" 
                                             xsd                                               dc 
             "http://www.w3.org/2001/XMLSchema#"               "http://purl.org/dc/elements/1.1/" 
                                         dcterms                                            prism 
                     "http://purl.org/dc/terms/" "http://prismstandard.org/namespaces/1.2/basic/" 
                                              cc                                             ncbi 
                "http://creativecommons.org/ns#"          "http://www.ncbi.nlm.nih.gov/taxonomy#" 
                                              tc                                                  
 "http://rs.tdwg.org/ontology/voc/TaxonConcept#"                      "http://www.nexml.org/2009" 
```

Recognize that `nex@namespaces` serves the same role as `get_namespaces`
function, but provides direct access to the slot data.  For instance,
with this syntax we could also overwrite the existing namespaces with
`nex@namespaces <- NULL`.  Changing the namespace in this way is not
advised.

Some slots can contain multiple elements of the same type, such as
`trees`, `characters`, and `otus`.  For instance, we see that


```r
class(nex@characters)
```

```
[1] "ListOfcharacters"
attr(,"package")
[1] "RNeXML"
```

is an object of class `ListOfcharacters`, and is currently empty,


```r
length(nex@characters)
```

```
[1] 0
```

In order to assign an object to a slot, it must match the class definition
of the slot.  We can create a new element of any given class with the
`new` function,


```r
nex@characters <- new("ListOfcharacters", list(new("characters")))
```

and now we have a length-1 list of character matrices,


```r
length(nex@characters)
```

```
[1] 1
```

and we access the first character matrix using the list notation,
`[[1]]`. Here we check the class is a `characters` object.


```r
class(nex@characters[[1]])
```

```
[1] "characters"
attr(,"package")
[1] "RNeXML"
```

Direct subsetting has two primary use cases: (a) useful in looking up
(and possibly editing) a specific value of an element, or (b) when adding
metadata annotations to specific elements. Consider the example file



```r
f <- system.file("examples", "trees.xml", package="RNeXML")
nex <- nexml_read(f)
```

We can look up the species label of the first `otu` in the first `otus` block:


```r
nex@otus[[1]]@otu[[1]]@label
```

```
      label 
"species 1" 
```

We can add metadata to this particular OTU using this subsetting format


```r
nex@otus[[1]]@otu[[1]]@meta <- 
  c(meta("skos:note", 
          "This species was incorrectly identified"),
         nex@otus[[1]]@otu[[1]]@meta)
```

Here we use the `c` operator to append this element to any existing meta annotations to this otu.  



---
title: "SPARQL with RNeXML"
author:
  - Carl Boettiger
  - Scott Chamberlain
  - Rutger Vos
  - Hilmar Lapp
output: html_vignette
---

<!--
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{SPARQL with RNeXML}
-->









## SPARQL Queries

Rich, semantically meaningful metadata lies at the heart of the NeXML
standard.  R provides a rich environment to unlock this information.
While our previous examples have relied on the user knowing exactly
what metadata they intend to extract (title, publication date, citation
information, and so forth), _semantic_ metadata has meaning that a
computer can make use of, allowing us to make much more conceptually
rich queries than those simple examples.  The SPARQL query language
is a powerful way to make use of such semantic information in making 
complex queries. 

While users should consult a formal introduction to SPARQL for further
background, here we illustrate how SPARQL can be used in combination
with R functions in ways that would be much more tedious to assemble
with only traditional/non-semantic queries.  The SPARQL query language
is provided for the R environment through the `rrdf` package [@Willighagen_2014],
so we start by loading that package.  We will also make use of functions
from `phytools` and `RNeXML`.  


```r
library("rdflib")
library("xml2")
library("phytools")
library("RNeXML")
```

We read in an example file that contains semantic metadata annotations
describing the taxonomic units (OTUs) used in the tree.  


```r
nexml <- nexml_read(system.file("examples/primates.xml", package="RNeXML"))
```

In particular, this example declares the taxon rank, NCBI identifier and parent taxon
for each OTU, such as:

```xml
<otu about="#ou541" id="ou541" label="Alouatta guariba">
      <meta href="http://ncbi.nlm.nih.gov/taxonomy/182256" 
            id="ma20" 
            rel="concept:toTaxon" 
            xsi:type="nex:ResourceMeta"/>
      <meta href="http://rs.tdwg.org/ontology/voc/TaxonRank#Species" 
            id="ma21" 
            rel="concept:rank" 
            xsi:type="nex:ResourceMeta"/>
      <meta href="http://ncbi.nlm.nih.gov/taxonomy/9499" 
            id="ma22" 
            rel="rdfs:subClassOf" 
            xsi:type="nex:ResourceMeta"/>
    </otu>

```

In this example, we will construct a cladogram by using this information
to identify the taxonomic rank of each OTU, and its shared parent
taxonomic rank.  (If this example looks complex, try writing down the
steps to do this without the aid of the SPARQL queries).  These examples
show the manipulation of semantic triples, Unique Resource Identifiers
(URIs) and use of the SPARQL "Join" operator.

Note that this example can be run using `demo("sparql", "RNeXML")` to see the 
code displayed in the R terminal and to avoid character errors that
can occur in having to copy and paste from PDF files.  

We begin by extracting the RDF graph from the NeXML,


```r
rdf <- get_rdf(system.file("examples/primates.xml", package="RNeXML"))
tmp <- tempfile()  # so we must write the XML out first
xml2::write_xml(rdf, tmp) 

graph <- rdf_parse(tmp)
```

We then fetch the NCBI URI for the taxon that has rank 'Order', i.e. the
root of the primates phylogeny. The dot operator `.` between clauses
implies a join, in this case


```r
root <- rdf_query(graph, 
"SELECT ?uri WHERE { 
    ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#rank> <http://rs.tdwg.org/ontology/voc/TaxonRank#Order> . 
    ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> ?uri    
}")
```


This makes use of the SPARQL query language provided by the `rrdf`
package.  We will also define some helper functions that use SPARQL
queries.  Here we define a function to get the name


```r
get_name <- function(id) {
  max <- length(nexml@otus[[1]]@otu)
  for(i in 1:max) {
    if ( nexml@otus[[1]]@otu[[i]]@id == id ) {
      label <- nexml@otus[[1]]@otu[[i]]@label
      label <- gsub(" ","_",label)
      return(label)
    }
  }
}
```


Next, we define a recursive function to build a newick tree from the taxonomic rank information.  


```r
recurse <- function(node){
  
    # fetch the taxonomic rank and id string
    rank_query <- paste0(
        "SELECT ?rank ?id WHERE {
            ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> <",node,"> .
            ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#rank> ?rank
          }")
    result <- rdf_query(graph, rank_query)
    
    # get the local ID, strip URI part
    id <- result$id
    id <- gsub("^.+#", "", id, perl = TRUE)
    
    # if rank is terminal, return the name
    if (result$rank == "http://rs.tdwg.org/ontology/voc/TaxonRank#Species") {
        return(get_name(id))
    }
    
    # recurse deeper
    else {
        child_query <- paste0(
            "SELECT ?uri WHERE {
                ?id <http://www.w3.org/2000/01/rdf-schema#subClassOf> <",node,"> .
                ?id <http://rs.tdwg.org/ontology/voc/TaxonConcept#toTaxon> ?uri
            }")
        children <- rdf_query(graph, child_query)
        
        return(paste("(", 
                     paste(sapply(children$uri, recurse), 
                           sep = ",", collapse = "," ), 
                     ")",  
                     get_name(id), # label interior nodes
                     sep = "", collapse = ""))
    }
}
```


With these functions in place, it is straight forward to build the tree
from the semantic RDFa data and then visualize it



```r
newick <- paste(recurse(root), ";", sep = "", collapse = "")
tree <- read.newick(text = newick)
collapsed <- collapse.singles(tree)
plot(collapsed, 
     type='cladogram', 
     show.tip.label=FALSE, 
     show.node.label=TRUE, 
     cex=0.75, 
     edge.color='grey60', 
     label.offset=-9)
```

![plot of chunk unnamed-chunk-8](sparql-unnamed-chunk-8-1.png)
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_level.R
\name{get_level}
\alias{get_level}
\title{get_level}
\usage{
get_level(nex, level)
}
\arguments{
\item{nex}{a nexml object}

\item{level}{a character vector indicating the class of node, see details}
}
\value{
Returns the attributes of specified class of nodes as a data.frame
}
\description{
get a data.frame of attribute values of a given node
}
\details{
level should be a character vector giving the path to the specified node
group.  For instance, \code{otus}, \code{characters}, and \code{trees} are top-level blocks (e.g.
child nodes of the root nexml block), and can be specified directly.  To get metadata
for all "char" elements from all characters blocks, you must specify that \code{char} nodes
are child nodes to \code{character} nodes: e.g. \code{get_level(nex, "characters/char")},
or similarly for states: \code{get_level(nex, characters/states)}.

The return object is a data frame whose columns are the attribute names of the elements
specified. The column names match the attribute names except for "id" attribute, for which the column
is renamed using the node itself. (Thus \verb{<otus id="os2">} would be rendered in a data.frame with column
called "otus" instead of "id"). Additional columns are
added for each parent element in the path; e.g. \code{get_level(nex, "otus/otu")} would include a column
named "otus" with the id of each otus block.  Even though the method always returns the data frame
for all matching nodes in all blocks, these ids let you see which otu values came from which
otus block.  This is identical to the function call \code{get_taxa()}.
Similarly, \code{get_level(nex, "otus/otu/meta")} would return additional columns 'otus' and
also a column, 'otu', with the otu parent ids of each metadata block.  (This is identical to a
function call to \code{get_metadata}).  This makes it easier to join data.frames as well, see examples
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_meta.R
\name{add_meta}
\alias{add_meta}
\title{Add metadata to a nexml file}
\usage{
add_meta(
  meta,
  nexml = new("nexml"),
  level = c("nexml", "otus", "trees", "characters"),
  namespaces = NULL,
  i = 1,
  at_id = NULL
)
}
\arguments{
\item{meta}{a meta S4 object, e.g. ouput of the function \code{\link{meta}}, or a list of these meta objects}

\item{nexml}{(S4) object}

\item{level}{the level at which the metadata annotation should be added.}

\item{namespaces}{named character string for any additional namespaces that should be defined.}

\item{i}{for otus, trees, characters: if there are multiple such blocks, which one should be annotated?  Default is first/only block.}

\item{at_id}{the id of the element to be annotated.  Optional, advanced use only.}
}
\value{
the updated nexml object
}
\description{
Add metadata to a nexml file
}
\examples{
## Create a new nexml object with a single metadata element: 
modified <- meta(property = "prism:modificationDate", content = "2013-10-04")
nex <- add_meta(modified) # Note: 'prism' is defined in nexml_namespaces by default.  

## Write multiple metadata elements, including a new namespace:  
website <- meta(href = "http://carlboettiger.info", 
                rel = "foaf:homepage")              # meta can be link-style metadata
nex <- add_meta(list(modified,  website), 
                namespaces = c(foaf = "http://xmlns.com/foaf/0.1/"))

## Append more metadata, and specify a level: 
history <- meta(property = "skos:historyNote",
                 content = "Mapped from the bird.orders data in the ape package using RNeXML")
data(bird.orders)
nex <- add_trees(bird.orders) # need to have created a trees block first
nex <- add_meta(history, 
                nexml = nex,
                level = "trees",
                namespaces = c(skos = "http://www.w3.org/2004/02/skos/core#"))

}
\seealso{
\code{\link{meta}} \code{\link{add_trees}} \code{\link{add_characters}} \code{\link{add_basic_meta}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/concatenate_nexml.R
\name{c,nexml-method}
\alias{c,nexml-method}
\title{Concatenate nexml files}
\usage{
\S4method{c}{nexml}(x, ..., recursive = FALSE)
}
\arguments{
\item{x, ...}{nexml objects to be concatenated, e.g. from
\code{\link{write.nexml}} or \code{\link{read.nexml}}.
Must have unique ids on all elements}

\item{recursive}{logical.  If 'recursive = TRUE', the function recursively
descends through lists (and pairlists) combining all their
elements into a vector. (Not implemented).}
}
\value{
a concatenated nexml file
}
\description{
Concatenate nexml files
}
\examples{
\dontrun{
f1 <- system.file("examples", "trees.xml", package="RNeXML")
f2 <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex1 <- read.nexml(f1)
nex2 <- read.nexml(f2)
nex <- c(nex1, nex2)
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.R
\name{slot,ResourceMeta-method}
\alias{slot,ResourceMeta-method}
\alias{slot-ResourceMeta}
\alias{slot<-,ResourceMeta-method}
\title{Access or set slot of S4 object}
\usage{
\S4method{slot}{ResourceMeta}(object, name)

\S4method{slot}{ResourceMeta}(object, name) <- value
}
\arguments{
\item{object}{the object}

\item{name}{name of the slot}

\item{value}{the new value}
}
\description{
See \code{\link[methods:slot]{methods::slot()}}. This version allows using "property" consistently
for both LiteralMeta and ResourceMeta (which internally uses "rel" because
RDFa does), which is easier to program. It also allows using "meta"
as an alias for "children" for ResourceMeta, to be consistent with the
corresponding slot for instances of \code{Annotated}.
}
\seealso{
\code{\link[methods:slot]{methods::slot()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nexml_get.R
\name{nexml_get}
\alias{nexml_get}
\alias{get_item}
\title{Get the desired element from the nexml object}
\usage{
nexml_get(
  nexml,
  element = c("trees", "trees_list", "flat_trees", "metadata", "otu", "taxa",
    "characters", "characters_list", "namespaces"),
  ...
)
}
\arguments{
\item{nexml}{a nexml object (from read_nexml)}

\item{element}{the kind of object desired, see details.}

\item{...}{additional arguments, if applicable to certain elements}
}
\value{
return type depends on the element requested.  See details.
}
\description{
Get the desired element from the nexml object
}
\details{
\itemize{
\item{"tree"}{ an ape::phylo tree, if only one tree is represented.  Otherwise returns a list of lists of multiphylo trees.  To consistently receive the list of lists format (preserving the hierarchical nature of the nexml), use \code{trees} instead.}
\item{"trees"}{ returns a list of lists of multiphylo trees, even if all trees are in the same \code{trees} node (and hence the outer list will be of length 1) or if there is only a single tree (and hence the inner list will also be of length 1.  This ensures a consistent return type regardless of the number of trees present in the nexml file, and also preserves any hierarchy/grouping of trees.  }
\item{"flat_trees"}{ a multiPhylo object (list of ape::phylo objects) Note that this method collapses any hierarchical structure that may have been present as multiple \code{trees} nodes in the original nexml (though such a feature is rarely used).  To preserve that structure, use \code{trees} instead.}
\item{"metadata"}{Get metadata from the specified level (default is top/nexml level) }
\item{"otu"}{ returns a named character vector containing all available metadata.  names indicate \code{property} (or \code{rel} in the case of links/resourceMeta), while values indicate the \code{content} (or \code{href} for links). }
\item{"taxa"}{ alias for otu }
}
For a slightly cleaner interface, each of these elements is also defined as an S4 method
for a nexml object.  So in place of \code{get_item(nexml, "tree")}, one could use \code{get_tree(nexml)},
and so forth for each element type.
}
\examples{
comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex <- nexml_read(comp_analysis)
nexml_get(nex, "trees")
nexml_get(nex, "characters_list")
}
\seealso{
\code{\link{get_trees}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4-utils.R
\name{.methodWithNext}
\alias{.methodWithNext}
\title{Saves the next method in the method meta data}
\usage{
.methodWithNext(method, nextMethod, .cache = FALSE)
}
\arguments{
\item{method}{the \code{MethodDefinition} object to promote}

\item{nextMethod}{the \code{MethodDefinition}
object to record as the next method}

\item{.cache}{whether to cache the promoted method definition object
(using \code{methods::cacheMethod()})}
}
\value{
an instance of \code{MethodWithNext},
which has the next method in the \code{nextMethod} slot
}
\description{
Promotes the given method definition to an instance of
\code{MethodWithNext}, thereby recording the next
method in the \code{nextMethod} slot.
}
\note{
\code{MethodWithNext} objects are normally returned by
\code{methods::addNextMethod()}, but a constructor function for the class
seems missing (or is undocumented?). This provides one.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_metadata.R
\name{get_metadata}
\alias{get_metadata}
\title{get_metadata}
\usage{
get_metadata(nexml, level = "nexml", simplify = TRUE)
}
\arguments{
\item{nexml}{a nexml object}

\item{level}{the name of the level of element desired, see details}

\item{simplify}{logical, see Details}
}
\value{
the requested metadata as a data.frame. Additional columns
indicate the parent element of the return value.
}
\description{
get_metadata
}
\details{
'level' should be either the name of a child element of a NeXML document
(e.g. "otu", "characters"), or a path to the desired element, e.g. 'trees/tree'
will return the metadata for all phylogenies in all trees blocks.

If a metadata element has other metadata elements nested within it, the
nested metadata are returned as well. A column "Meta" will contain the
IDs consolidated from the type-specific LiteralMeta and ResourceMeta
columns, and IDs are generated for meta elements that have nested elements
but do not have an ID ("blank nodes"). A column "meta" contains the
IDs of the parent meta elements for nested ones. This means that the
resulting table can be self-joined on those columns.

If \code{simplify} is \code{FALSE}, the type-specific "LiteralMeta" and "ResourceMeta"
columns will be retained even if a consolidated "Meta" column is present.
Otherwise, only the consolidated column will be included in the result.
Also, if \code{simplify} is \code{TRUE} the values for "property" (LiteralMeta) and
"rel" (ResourceMeta) will be consolidated to "property", and "rel" will be
removed from the result.
}
\examples{
\dontrun{
comp_analysis <- system.file("examples", "primates.xml", package="RNeXML")
nex <- nexml_read(comp_analysis)
get_metadata(nex)
get_metadata(nex, "otus/otu")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{lcapply}
\alias{lcapply}
\title{Compact list then lapply}
\usage{
lcapply(X, ...)
}
\arguments{
\item{X}{the list object}

\item{...}{remaining arguments to \code{lapply()}}
}
\description{
Compacts the list (i.e., removes NULL objects), then calls \code{\link[base:lapply]{lapply()}}
on the result with the remaining parameters.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructors.R
\name{New}
\alias{New}
\title{new with namespaced class name}
\usage{
New(Class, ...)
}
\arguments{
\item{Class}{the name of the S4 class to be instantiated}

\item{...}{additional parameters for \code{methods::new()}}
}
\description{
Convenience function for \code{\link[methods:new]{methods::new()}} that ensures that the provided
class name is namespaced with a package name.
}
\details{
If the provided class name is not already namespaced (see
\code{methods::packageSlot()}), it will be namespaced with this package. This
mechanism is used by \code{new()} to disambiguate if the class name clashes
with a class defined in another package.
}
\note{
This may not completely eliminate messages on standard error about
classes with the same name having been found in different packages. If
they appear, they will most likely have come from the call to the
\code{methods::initialize()} generic that \code{new()} issues at the end.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{Annotated-class}
\alias{Annotated-class}
\title{Class of objects that have metadata as lists of meta elements}
\description{
Class of objects that have metadata as lists of meta elements
}
\section{Slots}{

\describe{
\item{\code{meta}}{list of \code{meta} objects}

\item{\code{about}}{for RDF extraction, the identifier for the resource that this
object is about}
}}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_metadata.R
\name{get_meta}
\alias{get_meta}
\title{Extracts meta objects matching properties}
\usage{
get_meta(nexml, annotated = NULL, props)
}
\arguments{
\item{nexml}{a nexml object}

\item{annotated}{the nexml component object from which to obtain metadata
annotations, or a list of such objects. Defaults to the nexml object itself.}

\item{props}{a character vector of property names for which to extract
metadata annotations}
}
\value{
a named list of the matching meta objects
}
\description{
Extracts the metadata annotations for the given property or properties,
and returns the result as a list of \code{meta} objects.
}
\details{
For matching property identifiers (i.e., URIs), prefixes in the input list
as well as in the \code{annotated} object will be expanded using the namespaces
of the \code{nexml} object. Names in the returned list are mapped to the
(possibly prefixed) form in the input list. The resulting list is flat,
and hence does not retain the nesting hierarchy in the object's annotation.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trees.R
\name{get_flat_trees}
\alias{get_flat_trees}
\title{get_flat_trees}
\usage{
get_flat_trees(nexml)
}
\arguments{
\item{nexml}{a representation of the nexml object from  which the data is to be retrieved}
}
\value{
a multiPhylo object (list of ape::phylo objects).  See details.
}
\description{
extract a single multiPhylo object containing all trees in the nexml
}
\details{
Note that this method collapses any hierarchical structure that may have been present as multiple \code{trees} nodes in the original nexml (though such a feature is rarely used).  To preserve that structure, use \code{\link{get_trees}} instead.
}
\examples{
comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex <- nexml_read(comp_analysis)
get_flat_trees(nex)
}
\seealso{
\code{\link{get_trees}} \code{\link{get_trees}} \code{\link{get_item}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxa.R
\name{get_taxa_list}
\alias{get_taxa_list}
\alias{get_otus_list}
\title{get_taxa_list}
\usage{
get_taxa_list(nexml)
}
\arguments{
\item{nexml}{a nexml object}
}
\value{
the list of taxa
}
\description{
Retrieve names of all species/otus otus (operational taxonomic units) included in the nexml
}
\seealso{
\code{\link{get_item}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxa_meta.R
\name{get_taxa_meta}
\alias{get_taxa_meta}
\title{get_taxa_meta}
\usage{
get_taxa_meta(nexml, what = "href")
}
\arguments{
\item{nexml}{a nexml object}

\item{what}{One of href, rel, id, or xsi:type}
}
\value{
the list of metadata for each taxon
}
\description{
Retrieve metadata of all species/otus otus (operational taxonomic units) included in the nexml
}
\examples{
\dontrun{
data(bird.orders)
birds <- add_trees(bird.orders)
birds <- taxize_nexml(birds, "NCBI")
RNeXML:::get_taxa_meta(birds)
RNeXML:::get_taxa_meta(birds, 'rel')
RNeXML:::get_taxa_meta(birds, 'id')
RNeXML:::get_taxa_meta(birds, 'xsi:type')
 }
}
\seealso{
\code{\link{get_item}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{coalesce_}
\alias{coalesce_}
\title{Front-end to dplyr::coalesce to deal with NULL vectors}
\usage{
coalesce_(...)
}
\arguments{
\item{...}{the vectors to coalesce on NA}
}
\value{
a vector of the same type and length as the last argument
}
\description{
Replaces any NULL argument with a vector of \code{NA}, and casts every vector
to the same type as the last vector. After that, calls \code{\link[dplyr:coalesce]{dplyr::coalesce()}}.
}
\seealso{
\code{\link[dplyr:coalesce]{dplyr::coalesce()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nexml_add.R
\name{nexml_add}
\alias{nexml_add}
\title{add elements to a new or existing nexml object}
\usage{
nexml_add(
  x,
  nexml = new("nexml"),
  type = c("trees", "characters", "meta", "namespaces"),
  ...
)
}
\arguments{
\item{x}{the object to be added}

\item{nexml}{an existing nexml object onto which the object should be appended}

\item{type}{the type of object being provided.}

\item{...}{additional optional arguments to the add functions}
}
\value{
a nexml object with the additional data
}
\description{
add elements to a new or existing nexml object
}
\examples{
library("geiger")
data(geospiza)
geiger_nex <- nexml_add(geospiza$phy, type="trees")
geiger_nex <- nexml_add(geospiza$dat, nexml = geiger_nex, type="characters")
}
\seealso{
\code{\link{add_trees}} \code{\link{add_characters}} \code{\link{add_meta}} \code{\link{add_namespaces}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_namespaces.R
\name{get_namespaces}
\alias{get_namespaces}
\title{get namespaces}
\usage{
get_namespaces(nexml)
}
\arguments{
\item{nexml}{a nexml object}
}
\value{
a named character vector providing the URLs defining each
of the namespaces used in the nexml file.  Names correspond to
the prefix abbreviations of the namespaces.
}
\description{
get namespaces
}
\examples{
comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex <- nexml_read(comp_analysis)
get_namespaces(nex)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/constructors.R, R/meta.R
\name{nexml.tree}
\alias{nexml.tree}
\alias{nexml.trees}
\alias{nexml.node}
\alias{nexml.edge}
\alias{nexml.otu}
\alias{nexml.otus}
\alias{nexml.char}
\alias{nexml.characters}
\alias{nexml.format}
\alias{nexml.state}
\alias{nexml.uncertain_state}
\alias{nexml.states}
\alias{nexml.uncertain_states}
\alias{nexml.polymorphic_states}
\alias{nexml.member}
\alias{nexml.matrix}
\alias{nexml.row}
\alias{nexml.seq}
\alias{nexml.cell}
\alias{nexml.meta_}
\title{Constructor for the respective class}
\usage{
nexml.tree(...)

nexml.trees(...)

nexml.node(...)

nexml.edge(...)

nexml.otu(...)

nexml.otus(...)

nexml.char(...)

nexml.characters(...)

nexml.format(...)

nexml.state(...)

nexml.uncertain_state(...)

nexml.states(...)

nexml.uncertain_states(...)

nexml.polymorphic_states(...)

nexml.member(...)

nexml.matrix(...)

nexml.row(...)

nexml.seq(...)

nexml.cell(...)
}
\arguments{
\item{...}{optionally, parameters passed on to \code{\link[methods:new]{new()}}}
}
\description{
Creates an instance of the class corresponding to the respective NeXML
element, and initializes its slots with the provided parameters, if any.
}
\details{
Usually, users won't need to invoke this directly.
}
\seealso{
\link[=meta]{nexml.meta()} for documentation of \code{nexml.meta()}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trees.R
\name{get_trees_list}
\alias{get_trees_list}
\title{extract all phylogenetic trees in ape format}
\usage{
get_trees_list(nexml)
}
\arguments{
\item{nexml}{a representation of the nexml object from
which the data is to be retrieved}
}
\value{
returns a list of lists of multiphylo trees, even if all trees
are in the same \code{trees} node (and hence the outer list will be of length
\enumerate{
\item or if there is only a single tree (and hence the inner list will also
be of length 1.  This ensures a consistent return type regardless of
the number of trees present in the nexml file, and also preserves any
grouping of trees.
}
}
\description{
extract all phylogenetic trees in ape format
}
\examples{
comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex <- nexml_read(comp_analysis)
get_trees_list(nex)
}
\seealso{
\code{\link{get_trees}} \code{\link{get_flat_trees}} \code{\link{get_item}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nexml_publish.R
\name{nexml_publish}
\alias{nexml_publish}
\title{publish nexml files to the web and receive a DOI}
\usage{
nexml_publish(nexml, ..., repository = "figshare")
}
\arguments{
\item{nexml}{a nexml object (or file path)}

\item{...}{additional arguments, depending on repository. See examples.}

\item{repository}{destination repository}
}
\value{
a digital object identifier to the published data
}
\description{
publish nexml files to the web and receive a DOI
}
\examples{
\dontrun{
data(bird.orders)
birds <- add_trees(bird.orders)
doi <- nexml_publish(birds, visibility = "public", repository="figshare")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trees.R
\name{toPhylo}
\alias{toPhylo}
\title{nexml to phylo}
\usage{
toPhylo(tree, otus)
}
\arguments{
\item{tree}{an nexml tree element}

\item{otus}{a character string of taxonomic labels, named by the otu ids.
e.g. (from get_otu_maps for the otus set matching the relevant trees node.}
}
\value{
phylo object.  If a "reconstructions" annotation is found on the
edges, return simmap maps slot as well.
}
\description{
nexml to phylo coercion
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.R
\name{meta}
\alias{meta}
\alias{nexml.meta}
\title{Constructor function for metadata nodes}
\usage{
meta(
  property = NULL,
  content = NULL,
  rel = NULL,
  href = NULL,
  datatype = NULL,
  id = NULL,
  type = NULL,
  children = list()
)
}
\arguments{
\item{property}{specify the ontological definition together with it's namespace, e.g. dc:title}

\item{content}{content of the metadata field}

\item{rel}{Ontological definition of the reference provided in href}

\item{href}{A link to some reference}

\item{datatype}{optional RDFa field}

\item{id}{optional id element (otherwise id will be automatically generated).}

\item{type}{optional xsi:type.  If not given, will use either "LiteralMeta" or "ResourceMeta" as
determined by the presence of either a property or a href value.}

\item{children}{Optional element containing any valid XML block (XMLInternalElementNode class, see the XML package for details).}
}
\description{
Constructor function for metadata nodes
}
\details{
User must either provide property+content or rel+href.  Mixing these will result in potential garbage.
The datatype attribute will be detected automatically from the class of the content argument.  Maps from R class
to schema datatypes are as follows:
character - xs:string,
Date - xs:date,
integer - xs:integer,
numeric - xs:decimal,
logical - xs:boolean
}
\examples{
meta(content="example", property="dc:title")
}
\seealso{
\code{\link{nexml_write}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/classes.R
\docType{class}
\name{nexml-class}
\alias{nexml-class}
\alias{nexml}
\title{Class representing a NeXML document}
\description{
The \code{nexml} class represents a NeXML document, and is the top of the
class hierarchy defined in this package, corresponding to the root node
of the corresponding XML document.
}
\details{
Normally objects of this type are created by the package as a result of
reading a NeXML file, or of converting from another type, such as
\code{ape::phylo}. Also, interacting directly with the slots of the class is
normally not necessary. Instead, use the \code{get_XXX()} and \code{add_XXX()}
functions in the API.
}
\section{Slots}{

\describe{
\item{\code{trees}}{list, corresponding to the list of \verb{<trees/>} elements in
NeXML. Elements will be of class \code{trees}.}

\item{\code{characters}}{list, corresponding to the list of \verb{<characters/>}
elements in NeXML. Elements will be of class \code{characters}.}

\item{\code{otus}}{list, corresponding to the list of \verb{<otus/>} elements in NeXML.
Elements will be of class \code{otus}.}

\item{\code{about}}{inherited, see \link[=Annotated-class]{Annotated}}

\item{\code{meta}}{inherited, see \link[=Annotated-class]{Annotated}}

\item{\code{xsi:type}}{for internal use}

\item{\code{version}}{NeXML schema version, do not change}

\item{\code{generator}}{name of software generating the XML}

\item{\code{xsi:schemaLocation}}{for internal use, do not change}

\item{\code{namespaces}}{named character vector giving the XML namespaces}
}}

\examples{
nex <- nexml() # a nexml object with no further content
nex <- new("nexml") # accomplishes the same thing
nex@generator
length(nex@trees)

data(bird.orders)
nex <- as(bird.orders, "nexml")
summary(nex)
length(nex@trees)
}
\seealso{
\code{\link[=read.nexml]{read.nexml()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_basic_metadata.R
\name{get_citation}
\alias{get_citation}
\title{Get citation from metadata}
\usage{
get_citation(nexml)
}
\arguments{
\item{nexml}{a nexml object}
}
\value{
the citation if the metadata provides one that is non-empty, and
NA otherwise. If multiple non-empty annotations are found, only the first
one is returned.
}
\description{
Extracts the citation annotation from the metadata annotation of the\code{nexml}
object, and returns its value.
}
\details{
Currently the implementation looks for \code{dcterms:bibliographicCitation}
annotations. (Note that these may be given with any prefix in the metadata
so long as they expand to the same full property URIs.)
}
\seealso{
\code{\link[=get_metadata_values]{get_metadata_values()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/prefixed-uris.R
\name{expand_prefix}
\alias{expand_prefix}
\title{Expand namespace-prefixed string}
\usage{
expand_prefix(x, namespaces = NULL)
}
\arguments{
\item{x}{a character vector of potentially namespace-prefixed strings}

\item{namespaces}{a named vector of namespaces, with namespace prefixes
being the names. A "base" namespace with an empty name can be included.
If not provided, or if empty, the input vector is returned as is.}
}
\value{
a character vector, of the same length as the input vector
}
\description{
Substitutes the namespace prefix in the input vector of strings with
the corresponding namespaces.
}
\details{
Namespace prefixes are expected to be separated by one or more semicolons.
Prefixes that cannot be matched to the vector of namespaces will be left
as is. For strings that do not have a namespace prefix, the vector of
namespaces can contain a base namespace, identified as not having a name,
with which these strings will be expanded.
}
\examples{
uris <- c("cc:license", "dc:title")
ns <- c(dc = "http://purl.org/dc/elements/1.1/",
        dcterms = "http://purl.org/dc/terms/",
        dct = "http://purl.org/dc/terms/",
        cc = "http://creativecommons.org/ns#")
# expansion is vectorized
expand_prefix(uris, ns)

# strings with non-matching or no prefix are left as is
uris <- c(uris, "my:title", "title")
expand_prefix(uris, ns)

# NAs in the input list become NA in the output
uris <- c(uris, NA)
expand_prefix(uris, ns)

# can include a "base" (unnamed) namespace for expanding unprefixed strings
ns <- c(ns, "//local/")
xuris <- expand_prefix(uris, ns)
xuris
xuris[uris == "title"] == paste0("//local/", uris[uris == "title"])

# different prefixes may expand to the same result
expand_prefix("dcterms:modified", ns) == expand_prefix("dct:modified", ns)

# or they may result in different expansions
expand_prefix("dc:title", ns) != expand_prefix("dcterms:title", ns)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nexml_validate.R
\name{nexml_validate}
\alias{nexml_validate}
\title{validate nexml using the online validator tool}
\usage{
nexml_validate(file, schema = CANONICAL_SCHEMA)
}
\arguments{
\item{file}{path to the nexml file to validate}

\item{schema}{URL of schema (for fallback method only, set by default).}
}
\value{
TRUE if the file is valid, FALSE or error message otherwise
}
\description{
validate nexml using the online validator tool
}
\details{
Requires an internet connection.  see http://www.nexml.org/nexml/phylows/validator for more information in debugging invalid files
}
\examples{
\dontrun{
data(bird.orders)
birds <- nexml_write(bird.orders, "birds_orders.xml")
nexml_validate("birds_orders.xml")
unlink("birds_orders.xml") # delete file to clean up
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_rdf.R
\name{get_rdf}
\alias{get_rdf}
\title{Extract rdf-xml from a NeXML file}
\usage{
get_rdf(file)
}
\arguments{
\item{file}{the name of a nexml file, or otherwise a nexml object.}
}
\value{
an RDF-XML object (XMLInternalDocument).  This can be manipulated with
tools from the XML R package, or converted into a triplestore for use with
SPARQL queries from the rdflib R package.
}
\description{
Extract rdf-xml from a NeXML file
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/meta.R
\name{c,meta-method}
\alias{c,meta-method}
\alias{c-meta}
\alias{c,ListOfmeta-method}
\alias{c-ListOfmeta}
\title{Concatenate meta elements into a ListOfmeta}
\usage{
\S4method{c}{meta}(x, ..., recursive = TRUE)

\S4method{c}{ListOfmeta}(x, ..., recursive = TRUE)
}
\arguments{
\item{x, ...}{\code{meta} and \code{ListOfmeta} elements to be concatenated, see \code{\link{meta}}}

\item{recursive}{logical, if 'recursive=TRUE', the function recursively
descends through lists and combines their elements into a flat vector.
This method does not support \code{recursive=FALSE}, use \link[base:list]{list}
instead.}
}
\value{
a ListOfmeta object containing a flat list of meta elements.
}
\description{
Concatenate meta elements into a ListOfmeta

Concatenate ListOfmeta elements into a flat ListOfmeta
}
\examples{
c(meta(content="example", property="dc:title"),
  meta(content="Carl", property="dc:creator"))
metalist <- c(meta(content="example", property="dc:title"),
              meta(content="Carl", property="dc:creator"))
out <- c(metalist, metalist) 
out <- c(metalist, meta(content="a", property="b")) 
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/taxize_nexml.R
\name{taxize_nexml}
\alias{taxize_nexml}
\title{taxize nexml}
\usage{
taxize_nexml(
  nexml,
  type = c("ncbi", "itis", "col", "tpl", "gbif", "wd"),
  warnings = TRUE,
  ...
)
}
\arguments{
\item{nexml}{a nexml object}

\item{type}{the name of the identifier to use}

\item{warnings}{should we show warning messages if no match can be found?}

\item{...}{additional arguments to \verb{[taxadb::get_ids()]}}
}
\description{
Check taxonomic names against the specified service and
add appropriate semantic metadata to the nexml OTU unit
containing the corresponding identifier.
}
\examples{
\dontrun{
data(bird.orders)
birds <- add_trees(bird.orders)
birds <- taxize_nexml(birds, "NCBI")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxa_meta.R
\name{get_taxa_meta_list}
\alias{get_taxa_meta_list}
\title{get_taxa_meta_list}
\usage{
get_taxa_meta_list(nexml, what = "href")
}
\arguments{
\item{nexml}{a nexml object}

\item{what}{One of href, rel, id, or xsi:type}
}
\value{
the list of metadata for each taxon
}
\description{
Retrieve metadata of all species/otus otus (operational taxonomic units) included in the nexml
}
\examples{
\dontrun{
data(bird.orders)
birds <- add_trees(bird.orders)
birds <- taxize_nexml(birds, "NCBI")
RNeXML:::get_taxa_meta_list(birds)
RNeXML:::get_taxa_meta_list(birds, 'rel')
RNeXML:::get_taxa_meta_list(birds, 'id')
RNeXML:::get_taxa_meta_list(birds, 'xsi:type')
}
}
\seealso{
\code{\link{get_item}}
}
\keyword{internal}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_metadata.R
\name{get_metadata_values}
\alias{get_metadata_values}
\title{Get the value(s) for metadata}
\usage{
get_metadata_values(nexml, annotated = NULL, props)
}
\arguments{
\item{nexml}{a nexml object}

\item{annotated}{the nexml component object from which to obtain metadata
annotations, defaults to the nexml object itself}

\item{props}{a character vector of property names for which to extract
metadata annotations}
}
\value{
a named character vector, giving the values and names being the
property names
}
\description{
Extracts the values from the metadata annotations for the given property
or properties, and returns the result.
}
\details{
For matching property identifiers (i.e., URIs), prefixes in the input list
as well as in the \code{annotated} object will be expanded using the namespaces
of the \code{nexml} object. Names in the returned vector are mapped to the
(possibly prefixed) form in the input list.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/utils.R
\name{charzero_as_empty}
\alias{charzero_as_empty}
\title{Treats zero-length character vectors as empty strings}
\usage{
charzero_as_empty(x)
}
\arguments{
\item{x}{the object to be tested for zero-length character vector}
}
\value{
an empty string if \code{x} is a character vector of length zero, and \code{x}
otherwise
}
\description{
If the argument is a zero-length character vector (character(0)), returns
an empty string (which is a character vector of length 1). Otherwise passes
through the argument.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/simmap.R
\name{simmap_to_nexml}
\alias{simmap_to_nexml}
\alias{nexml_to_simmap}
\title{Convert phylo with attached simmap to nexml object}
\usage{
simmap_to_nexml(phy, state_ids = NULL)

nexml_to_simmap(nexml)
}
\arguments{
\item{phy}{a phylo object containing simmap \code{phy$maps} element,
from the phytools package}

\item{state_ids}{a named character vector giving the state
names corresponding to the ids used to refer to each state
in nexml.  If null ids will be generated and states taken from
the phy$states names.}

\item{nexml}{a nexml object}
}
\value{
a nexml representation of the simmap

a simmap object (phylo object with a \verb{$maps} element
for use in phytools functions).
}
\description{
Convert phylo with attached simmap to nexml object

Convert nexml object with simmap to phylo
}
\section{Functions}{
\itemize{
\item \code{nexml_to_simmap}: Convert nexml object with simmap to phylo
}}

\examples{
simmap_ex <- read.nexml(system.file("examples","simmap_ex.xml", package="RNeXML"))
phy <- nexml_to_simmap(simmap_ex)
nex <- simmap_to_nexml(phy) 
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4-utils.R
\name{.callGeneric}
\alias{.callGeneric}
\title{Calls the given generic with the given arguments}
\usage{
.callGeneric(f, ..., .package = NULL)
}
\arguments{
\item{f}{the generic, as a character string or a \code{standardGeneric}
object}

\item{...}{the arguments (named and/or unnamed) with which to call the
matching method}

\item{.package}{the package name for finding the generic (if \code{f} is a character
string); by default the package is determined from the calling environment}
}
\value{
the value returned by the method
}
\description{
Calls the given generic with the given arguments, using the method
whose signature matches the arguments.
}
\details{
Uses \code{methods::selectMethod()} to find the matching method. In theory,
this is at best wholly redundant with what standard S4 generics already
do by themselves. However, the generics dispatch for S4 seems (at least
currently) broken at least if the first argument in the signature is
a class that name-clashes with a class defined in another package. In
that case, whether the standard dispatch works correctly or not can depend
on \code{search()} order, and can change within a session
depending on the order in which packages are loaded.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nexml_read.R
\name{nexml_read}
\alias{nexml_read}
\alias{read.nexml}
\alias{nexml_read.character}
\alias{nexml_read.XMLInternalDocument}
\alias{nexml_read.XMLInternalNode}
\title{Read NeXML files into various R formats}
\usage{
nexml_read(x, ...)

\method{nexml_read}{character}(x, ...)

\method{nexml_read}{XMLInternalDocument}(x, ...)

\method{nexml_read}{XMLInternalNode}(x, ...)
}
\arguments{
\item{x}{Path to the file to be read in. An \code{XML::XMLDocument-class}
or \code{\link[XML]{XMLNode-class}}}

\item{...}{Further arguments passed on to \code{\link[XML]{xmlTreeParse}}}
}
\description{
Read NeXML files into various R formats
}
\examples{
# file
f <- system.file("examples", "trees.xml", package="RNeXML")
nexml_read(f)
\dontrun{ # may take > 5 s
# url
url <- "https://raw.githubusercontent.com/ropensci/RNeXML/master/inst/examples/trees.xml"
nexml_read(url)
# character string of XML
str <- paste0(readLines(f), collapse = "")
nexml_read(str)
# XMLInternalDocument
library("httr")
library("XML")
x <- xmlParse(content(GET(url)))
nexml_read(x)
# XMLInternalNode
nexml_read(xmlRoot(x))
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nexml_publish.R
\name{nexml_figshare}
\alias{nexml_figshare}
\title{publish nexml to figshare}
\usage{
nexml_figshare(
  nexml,
  file = "nexml.xml",
  categories = "Evolutionary Biology",
  tags = list("phylogeny", "NeXML"),
  visibility = c("public", "private", "draft"),
  id = NULL,
  ...
)
}
\arguments{
\item{nexml}{a nexml object (or file path to a nexml file)}

\item{file}{The filename desired for the object, if nexml is not already a file.
if the first argument is already a path, this value is ignored.}

\item{categories}{The figshare categories, must match available set. see \code{fs_add_categories}}

\item{tags}{Any keyword tags you want to add to the data.}

\item{visibility}{whether the results should be published (public), or kept private,
or kept as a draft for further editing before publication.  (New versions can be updated,
but any former versions that was once made public will always be archived and cannot be removed).}

\item{id}{an existing figshare id (e.g. from fs_create), to which this file can be appended.}

\item{...}{additional arguments}
}
\value{
the figshare id of the object
}
\description{
publish nexml to figshare
}
\examples{
\dontrun{
data(bird.orders)
birds <- add_trees(bird.orders)
doi <- nexml_figshare(birds, visibility = "public", repository="figshare")
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_taxa.R
\name{get_taxa}
\alias{get_taxa}
\alias{get_otu}
\title{get_taxa}
\usage{
get_taxa(nexml)
}
\arguments{
\item{nexml}{a nexml object}
}
\value{
the list of taxa
}
\description{
Retrieve names of all species/otus otus (operational taxonomic units) included in the nexml
}
\examples{
comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex <- nexml_read(comp_analysis)
get_taxa(nex)
}
\seealso{
\code{\link{get_item}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/deprecated.R
\name{get_characters_list}
\alias{get_characters_list}
\title{Extract the character matrix}
\usage{
get_characters_list(nexml, rownames_as_col = FALSE)
}
\arguments{
\item{nexml}{nexml object (e.g. from read.nexml)}

\item{rownames_as_col}{option to return character matrix rownames
(with taxon ids) as it's own column in the data.frame. Default is FALSE
for compatibility with geiger and similar packages.}
}
\value{
the list of taxa
}
\description{
Extract the character matrix
}
\examples{
comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex <- nexml_read(comp_analysis)
get_characters_list(nex)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trees.R
\name{get_trees}
\alias{get_trees}
\title{extract a phylogenetic tree from the nexml}
\usage{
get_trees(nexml)
}
\arguments{
\item{nexml}{a representation of the nexml object from
which the data is to be retrieved}
}
\value{
an ape::phylo tree, if only one tree is represented.
Otherwise returns a list of lists of multiphylo trees.
To consistently receive the list of lists format (preserving
the hierarchical nature of the nexml), use \code{\link{get_trees_list}} instead.
}
\description{
extract a phylogenetic tree from the nexml
}
\examples{
comp_analysis <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex <- nexml_read(comp_analysis)
get_trees(nex)
}
\seealso{
\code{\link{get_trees}} \code{\link{get_flat_trees}} \code{\link{get_item}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_metadata.R
\name{get_all_meta}
\alias{get_all_meta}
\title{Get flattened list of meta annotations}
\usage{
get_all_meta(annotated)
}
\arguments{
\item{annotated}{the object from which to extract meta object annotations}
}
\value{
a flat list of \code{meta} objects
}
\description{
Collects recursively (in the case of nested meta annotations) all meta
object annotations for the given object, and returns the result as a flat
list.
}
\details{
Does not check that the input object can actually have meta annotations.
An invalid slot error will be generated if it can't.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4-utils.R
\name{.cacheNextMethod}
\alias{.cacheNextMethod}
\title{Caches next method in the calling environment}
\usage{
.cacheNextMethod()
}
\description{
If the calling environment does not have the next method to be invoked
in the inheritance chain cached yet, this will find the next method
(using \code{findNextMethod()}, and cache it in the calling environment such
that a subsequent call to \code{methods::callNextMethod()} will find and use
it.
}
\details{
As per the description, what this function does would normally already
be done by invoking \code{methods::callNextMethod()}, so in theory this should
be entirely redundant at best. However, \code{methods::addNextMethod()}, which
is invoked by \code{callNextMethod()} if a next method isn't cached yet, is
broken (errors out) if one of the classes in the signature name-clashes
with a class defined in another package. Calling this function prior to
\code{callNextMethod()} is meant to work around that.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_basic_meta.R
\name{add_basic_meta}
\alias{add_basic_meta}
\title{Add basic metadata}
\usage{
add_basic_meta(
  title = NULL,
  description = NULL,
  creator = Sys.getenv("USER"),
  pubdate = NULL,
  rights = "CC0",
  publisher = NULL,
  citation = NULL,
  nexml = new("nexml")
)
}
\arguments{
\item{title}{A title for the dataset}

\item{description}{a description of the dataset}

\item{creator}{name of the data creator.  Can be a string or R person object}

\item{pubdate}{publication date.  Default is current date.}

\item{rights}{the intellectual property rights associated with the data.
The default is Creative Commons Zero (CC0) public domain declaration,
compatible with all other licenses and appropriate for deposition
into the Dryad or figshare repositories.  CC0 is also recommended by the Panton Principles.
Alternatively, any other plain text string can be added and will be provided as the content
attribute to the dc:rights property.}

\item{publisher}{the publisher of the dataset. Usually where a user may go to find the canonical
copy of the dataset: could be a repository, journal, or academic institution.}

\item{citation}{a citation associated with the data.  Usually an academic journal
article that indicates how the data should be cited in an academic context.  Multiple citations
can be included here.
citation can be a plain text object, but is preferably an R \code{citation} or \code{bibentry} object (which
can include multiple citations.  See examples}

\item{nexml}{a nexml object to which metadata should be added.  A new
nexml object will be created if none exists.}
}
\value{
an updated nexml object
}
\description{
adds Dublin Core metadata elements to (top-level) nexml
}
\details{
\code{add_basic_meta()} is just a wrapper for \code{\link{add_meta}} to make it easy to
provide generic metadata without explicitly providing the namespace.  For instance,
\code{add_basic_meta(title="My title", description="a description")} is identical to:
\code{add_meta(list(meta("dc:title", "My title"), meta("dc:description", "a description")))}
Most function arguments are mapped directly to the Dublin Core terms
of the same name, with the exception of \code{rights}, which by default maps
to the Creative Commons namespace when using CC0 license.
}
\examples{
nex <- add_basic_meta(title = "My test title",
             description = "A description of my test",
             creator = "Carl Boettiger <cboettig@gmail.com>",
             publisher = "unpublished data",
             pubdate = "2012-04-01")

 ## Adding citation to an R package:
 nexml <- add_basic_meta(citation=citation("ape"))
\dontrun{
 ## Use knitcitations package to add a citation by DOI:
 library(knitcitations)
 nexml <- add_basic_meta(citation = bib_metadata("10.2307/2408428"))
 }
}
\seealso{
\code{\link{add_trees}} \code{\link{add_characters}} \code{\link{add_meta}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nexml_methods.R
\name{summary,nexml-method}
\alias{summary,nexml-method}
\alias{summary.nexml}
\title{Summary method for nexml objects}
\usage{
\S4method{summary}{nexml}(object)
}
\arguments{
\item{object}{the \link[=nexml-class]{nexml} object}
}
\value{
A list with the following elements:
\itemize{
\item \code{nblocks} the number of trees, otus, and characters blocks
\item \code{ncharacters} the number of characters in each characters block
\item \code{nstates} summary statistics of the number of character states per state set
defined for each characters block
\item \code{nnonstdstatedefs} the number of polymorphic and uncertain states defined
for each character block
\item \code{nmatrixrows} the number of rows in the matrix for each character block
\item \code{ntrees} the number of trees contained in each trees block
\item \code{notus} the number of OTUs defined in each OTUs block
\item \code{nmeta} a list of the number of the number of metadata annotations at
several levels, specifically:
\itemize{
\item \code{nexml} at the top (nexml) level
\item \code{otu} at the OTU level, for each OTUs block
\item \code{char} at the character level, for each characters block
\item \code{state} at the character state level, for each characters block
}
}
}
\description{
Generates a list of various counts of the major elements that comprise a
\link[=nexml-class]{nexml} object, such as number of different kinds of blocks,
characters, states, OTUs (taxa), etc.
}
\details{
The \link[methods:show]{show} method uses this summary for pretty-printing a
summary of the NeXML object, but it can be used on its own as well, in
particular for quick inspection of key properties of a NeXML file.
}
\examples{
nex <- nexml_read(system.file("examples", "comp_analysis.xml", package = "RNeXML"))
s <- summary(nex)
# number of major blocks:
s$nblocks

# each characters block defines 1 character:
s$ncharacters

# summary stats of states per character (for morphological matrices there is
# typically one state set per character)
s$nstates # note that first block is of continuous type, so no stats there

# pretty-printed summary:
nex # this is the same as show(nex)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_trees.R
\name{flatten_multiphylo}
\alias{flatten_multiphylo}
\title{Flatten a multiphylo object}
\usage{
flatten_multiphylo(object)
}
\arguments{
\item{object}{a list of multiphylo objects}
}
\description{
Flatten a multiphylo object
}
\details{
NeXML has the concept of multiple \verb{<trees>} nodes, each with multiple child \verb{<tree>} nodes.
This maps naturally to a list of multiphylo  objects.  Sometimes
this hierarchy conveys important structural information, so it is not discarded by default.
Occasionally it is useful to flatten the structure though, hence this function.  Note that this
discards the original structure, and the nexml file must be parsed again to recover it.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/internal_nexml_id.R
\name{reset_id_counter}
\alias{reset_id_counter}
\title{reset id counter}
\usage{
reset_id_counter()
}
\description{
reset the id counter
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4-utils.R
\name{findNextMethod}
\alias{findNextMethod}
\title{Finds the method that callNextMethod() should chain to}
\usage{
findNextMethod(method, f = NULL, envir = topenv())
}
\arguments{
\item{method}{\code{MethodDefinition}, the method for which to find
the next method}

\item{f}{\code{standardGeneric}, the standard generic for which to find
the next method. By default this will be obtained from \code{method}.}

\item{envir}{the environment in which to find the method}
}
\value{
a \code{MethodDefinition} object that is the next method in the
chain by inheritance
}
\description{
Attempts to find the "next" method in the inheritance chain. This would
(ideally) be the method that \code{methods::callNextMethod()} would chain to,
as a result of the method \code{methods::addNextMethod()} would find (and return
in the \code{nextMethod} slot of the \code{MethodWithNext}
object). Hence, in theory one shouldn't ever need this, but unfortunately
\code{addNextMethod()} is broken (and errors out) if one of the classes in the
signature name-clashes with an S4 class defined in another package that is
loaded.
}
\details{
The next method will be determined by the S4 inheritance chain. However,
this function will walk only the inheritance chain of those arguments in
the signature that are defined in the package of the generic method from
which this function was invoked (directly or indirectly). If there are
no such parameters in the signature, or if there is more than one,
finding the next method is handed off to \code{methods::addNextMethod()}.
}
\note{
In theory a class name clash between packages shouldn't be a problem
because class names can be namespaced, and the \code{MethodDefinition}
object passed to \code{addNextMethod()} has all the necessary namespace
information. Hopefully, at some point this gets fixed in R, and then we
don't need this anymore.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/S4-utils.R
\name{.sigLabel}
\alias{.sigLabel}
\title{Create a label for a method signature}
\usage{
.sigLabel(signature)
}
\arguments{
\item{signature}{the signature for which to create a label, as a vector
or list of strings, or as an instance of \code{signature}.}
}
\value{
a character string
}
\description{
Creates a label for a signature mirroring the result of \code{.sigLabel()}
in the \code{methods} package, which unfortunately does not export the function.
This is needed, for example, for the \code{excluded} slot in the
\code{MethodWithNext} class.
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/nexml_write.R
\name{nexml_write}
\alias{nexml_write}
\alias{write.nexml}
\title{Write nexml files}
\usage{
nexml_write(
  x = nexml(),
  file = NULL,
  trees = NULL,
  characters = NULL,
  meta = NULL,
  ...
)
}
\arguments{
\item{x}{a nexml object, or any phylogeny object (e.g. phylo, phylo4)
that can be coerced into one. Can also be omitted, in which case a new
nexml object will be constructed with the additional parameters specified.}

\item{file}{the name of the file to write out}

\item{trees}{phylogenetic trees to add to the nexml file (if not already given in x)
see \code{\link{add_trees}} for details.}

\item{characters}{additional characters}

\item{meta}{A meta element or list of meta elements, see \code{\link{add_meta}}}

\item{...}{additional arguments to add__basic_meta, such as the title.  See \code{\link{add_basic_meta}}.}
}
\value{
Writes out a nexml file
}
\description{
Write nexml files
}
\examples{
 ## Write an ape tree to nexml, analgous to write.nexus:
 library(ape); data(bird.orders)
 ex <- tempfile(fileext=".xml")
 write.nexml(bird.orders, file=ex)

\dontrun{ # takes > 5s
 ## Assemble a nexml section by section and then write to file:
 library(geiger)
 data(geospiza)
 nexml <- add_trees(geospiza$phy) # creates new nexml
 nexml <- add_characters(geospiza$dat, nexml = nexml) # pass the nexml obj to append character data
 nexml <- add_basic_meta(title="my title", creator = "Carl Boettiger", nexml = nexml)
 nexml <- add_meta(meta("prism:modificationDate", format(Sys.Date())), nexml = nexml)

 ex <- tempfile(fileext=".xml")
 write.nexml(nexml, file=ex)

 ## As above, but in one call (except for add_meta() call).  
 write.nexml(trees = geospiza$phy, 
             characters = geospiza$dat, 
             title = "My title", 
             creator = "Carl Boettiger",
             file = ex)
 
 ## Mix and match: identical to the section by section: 
 nexml <- add_meta(meta("prism:modificationDate", format(Sys.Date())))
 write.nexml(x = nexml,
             trees = geospiza$phy, 
             characters = geospiza$dat, 
             title = "My title", 
             creator = "Carl Boettiger",
             file = ex)

}
}
\seealso{
\code{\link{add_trees}} \code{\link{add_characters}} \code{\link{add_meta}} \code{\link{nexml_read}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/get_characters.R
\name{get_characters}
\alias{get_characters}
\title{Get character data.frame from nexml}
\usage{
get_characters(
  nex,
  rownames_as_col = FALSE,
  otu_id = FALSE,
  otus_id = FALSE,
  include_state_types = FALSE
)
}
\arguments{
\item{nex}{a nexml object}

\item{rownames_as_col}{option to return character matrix rownames (with taxon ids) as it's own column in the
data.frame. Default is FALSE for compatibility with geiger and similar packages.}

\item{otu_id}{logical, default FALSE. return a column with the
otu id (for joining with otu metadata, etc)}

\item{otus_id}{logical, default FALSE. return a column with the
otus block id (for joining with otu metadata, etc)}

\item{include_state_types}{logical, default FALSE. whether to also return a
matrix of state types (with values standard, polymorphic, and uncertain)}
}
\value{
the character matrix as a data.frame, or if \code{include_state_types} is
TRUE a list of two elements, \code{characters} as the character matrix, and
\code{state_types} as a matrix of state types. Both matrices will be in the same
ordering of rows and columns.
}
\description{
Get character data.frame from nexml
}
\details{
RNeXML will attempt to return the matrix using the NeXML taxon (otu) labels to name the rows
and the NeXML char labels to name the traits (columns).  If these are unavailable or not unique, the NeXML
id values for the otus or traits will be used instead.
}
\examples{
\dontrun{
# A simple example with a discrete and a continous trait
f <- system.file("examples", "comp_analysis.xml", package="RNeXML")
nex <- read.nexml(f)
get_characters(nex)

# A more complex example -- currently ignores sequence-type characters
f <- system.file("examples", "characters.xml", package="RNeXML")
nex <- read.nexml(f)
get_characters(nex)

# if polymorphic or uncertain states need special treatment, request state
# types to be returned as well:
f <- system.file("examples", "ontotrace-result.xml", package="RNeXML")
nex <- read.nexml(f)
res <- get_characters(nex, include_state_types = TRUE)
row.has.p <- apply(res$state_types, 1, 
                   function(x) any(x == "polymorphic", na.rm = TRUE))
col.has.p <- apply(res$state_types, 2, 
                   function(x) any(x == "polymorphic", na.rm = TRUE))
res$characters[row.has.p, col.has.p, drop=FALSE] # polymorphic rows and cols
res$characters[!row.has.p, drop=FALSE] # drop taxa with polymorphic states
# replace polymorphic state symbols in matrix with '?'
m1 <- mapply(function(s, s.t) ifelse(s.t == "standard", s, "?"), 
             res$characters, res$state_types)
row.names(m1) <- row.names(res$characters)
m1
}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_namespaces.R
\name{add_namespaces}
\alias{add_namespaces}
\title{Add namespaces}
\usage{
add_namespaces(namespaces, nexml = new("nexml"))
}
\arguments{
\item{namespaces}{a named character vector of namespaces}

\item{nexml}{a nexml object. will create a new one if none is given.}
}
\value{
a nexml object with updated namespaces
}
\description{
Add namespaces and their prefixes as a named vector of URIs, with the
names being the prefixes. Namespaces have most relevance for meta objects'
\code{rel} and \code{property}, and for embedded XML literals.
}
\details{
The implementation attempts to avoid duplication, currently using the
prefix. I.e., namespaces with prefixes already defined will not get added.
Namespaces needed by the NeXML format, and for commonly used metadata
terms, are already included by default, see \code{\link[=get_namespaces]{get_namespaces()}}.
}
\note{
Often a user won't call this directly, but instead provide the
namespace(s) through \code{\link[=add_meta]{add_meta()}}.
}
\examples{
## Write multiple metadata elements, including a new namespace:  
website <- meta(href = "http://carlboettiger.info", 
                rel = "foaf:homepage")     # meta can be link-style metadata
modified <- meta(property = "prism:modificationDate",
                 content = "2013-10-04")
nex <- add_meta(list(modified,  website), 
                namespaces = c(foaf = "http://xmlns.com/foaf/0.1/"))
                # prism prefix already included by default

## Add namespace "by hand" before adding meta:
nex <- add_namespaces(c(skos = "http://www.w3.org/2004/02/skos/core#"),
                      nexml = nex)
history <- meta(property = "skos:historyNote",
                content = "Mapped from the bird.orders data in the ape package using RNeXML")
nex <- add_meta(history, nexml = nex)

}
\seealso{
\code{\link[=meta]{meta()}} \code{\link[=add_meta]{add_meta()}} \code{\link[=get_namespaces]{get_namespaces()}}
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_trees.R
\name{add_trees}
\alias{add_trees}
\title{add_trees}
\usage{
add_trees(phy, nexml = new("nexml"), append_to_existing_otus = FALSE)
}
\arguments{
\item{phy}{a phylo object, multiPhylo object, or list of
multiPhylo to be added to the nexml}

\item{nexml}{a nexml object to which we should append this phylo.
By default, a new nexml object will be created.}

\item{append_to_existing_otus}{logical, indicating if we should
make a new OTU block (default) or append to the existing one.}
}
\value{
a nexml object containing the phy in nexml format.
}
\description{
add_trees
}
\examples{
library("geiger")
data(geospiza)
geiger_nex <- add_trees(geospiza$phy)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/add_characters.R
\name{add_characters}
\alias{add_characters}
\title{Add character data to a nexml object}
\usage{
add_characters(x, nexml = new("nexml"), append_to_existing_otus = FALSE)
}
\arguments{
\item{x}{character data, in which character traits labels are column names
and taxon labels are row names.  x can be in matrix or data.frame
format.}

\item{nexml}{a nexml object, if appending character table to an existing
nexml object.  If omitted will initiate a new nexml object.}

\item{append_to_existing_otus}{logical. If TRUE, will add any new taxa
(taxa not matching any existing otus block) to the existing (first)
otus block.  Otherwise (default), a new otus block is created, even
though it may contain duplicate taxa to those already present.  While
FALSE is the safe option, TRUE may be appropriate when building nexml
files from scratch with both characters and trees.}
}
\description{
Add character data to a nexml object
}
\examples{
library("geiger")
data(geospiza)
geiger_nex <- add_characters(geospiza$dat)
}
