
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
