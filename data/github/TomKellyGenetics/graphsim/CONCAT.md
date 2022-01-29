[![Travis-CI Build Status](https://travis-ci.com/TomKellyGenetics/graphsim.svg?branch=master)](https://travis-ci.com/TomKellyGenetics/graphsim)
[![CircleCI build status](https://circleci.com/gh/TomKellyGenetics/graphsim.svg?style=svg)](https://circleci.com/gh/TomKellyGenetics/graphsim)
[![AppVeyor Build Status](https://ci.appveyor.com/api/projects/status/github/TomKellyGenetics/graphsim?branch=master&svg=true)](https://ci.appveyor.com/project/TomKellyGenetics/graphsim)
[![Project Status: Active – The project has reached a stable, usable state and is being actively developed.](https://www.repostatus.org/badges/latest/active.svg)](https://www.repostatus.org/#active)

[![CRAN_Status_Badge](http://www.r-pkg.org/badges/version/graphsim)](https://cran.r-project.org/package=graphsim)
[![Downloads](https://cranlogs.r-pkg.org/badges/graphsim)](https://CRAN.R-project.org/package=graphsim)
[![Total Downloads](https://cranlogs.r-pkg.org/badges/grand-total/graphsim?color=orange)](https://CRAN.R-project.org/package=graphsim)

[![DOI](https://zenodo.org/badge/67395230.svg)](https://zenodo.org/badge/latestdoi/67395230)
[![bioRxiv](https://img.shields.io/badge/biorXiv-10.1101%2F2020.03.02.972471-blue)](https://doi.org/10.1101/2020.03.02.972471)
[![status](https://joss.theoj.org/papers/96016c6a55d7f74bacebd187c6ededd6/status.svg)](https://joss.theoj.org/papers/10.21105/joss.02161)
[![](https://img.shields.io/badge/Altmetric-72-blue.svg)](https://www.altmetric.com/details/77053356)

--------------------------------------------------

# graphsim

## Version 1.0.2


###  Simulate Expression Data from 'igraph' Networks 

This package provides functions to develop simulated continuous data 
(e.g., gene expression) from a sigma covariance matrix derived from a 
graph structure in 'igraph' objects. Intended to extend 'mvtnorm' to 
take 'igraph' structures rather than sigma matrices as input. This 
allows the use of simulated data that correctly accounts for pathway
relationships and correlations. Here we present a versatile statistical 
framework to simulate correlated gene expression data from biological 
pathways, by sampling from a multivariate normal distribution derived 
from a graph structure. This package allows the simulation of biological
pathways from a graph structure based on a statistical model of 
gene expression, such as simulation of expression profiles that
of log-transformed and normalised data from microarray and RNA-Seq data.

#### Motivation

Network analysis of molecular biological pathways is important
for insights into biology and medical genetics. 
Gene expression profiles capture the regulatory state of a cell
and can be used to analyse complex molecular states with genome-scale data.
Biological pathways are more than simply sets of genes involved in functions,
they are rich in information of relationships defined by pathway structure.

Methods to infer biological pathways and gene regulatory networks from gene
expression data can be tested on  simulated datasets using this framework. This also allows for
pathway structures to be considered as a confounding variable when 
simulating gene expression data to test the performance of genomics analyses.

This package enables the generation of simulated gene expression datasets
containing pathway relationships from a known underlying network.
These simulated datasets can be used to evaluate various bioinformatics
methodologies, including statistical and network inference procedures.

Network analysis techniques have an important role in understanding
of biological pathways and interpretation of genomics studies.
Modelling biological pathways allows the evaluation of gene
regulatory network inference techniques (which so far rely on
experimental validation or resampling). This technique also
enables modelling datasets with correlated pathway-structures
to assess whether other genomics analysis techniques perform
as expected with the background of complex pathways.


## Installation

To install the latest release from CRAN:

```R
install.packages("graphsim")
```

To install the stable release of this package from github:

```R
# install.packages("devtools")
devtools::install_github("TomKellyGenetics/graphsim", ref = "master")
```

To get the development version of this package from github:

```R
# install.packages("devtools")
devtools::install_github("TomKellyGenetics/graphsim", ref = "dev")
```

## Usage

Please see the vignettes for demonstrations of this package on examples of simple simulated networks and the reactome pathway TGF-&beta; receptor signaling activates SMADs (R-HSA-2173789). An [article](https://doi.org/10.21105/joss.02161) with further details has been published in the 
_Journal of Open Source Software_.

A help menu can also be accessed within the R environment:

```
?graphsim
```

```
help("graphsim-package")
```

This will display a help page and link to documentation for each function.

--------------------------------------------------

## Citation

To cite package 'graphsim' in publications use:

>S. Thomas Kelly and Michael A. Black (2020). graphsim: Simulate Expression Data from
>'igraph' Networks. R package version 1.0.2.
>https://github.com/TomKellyGenetics/graphsim doi:10.5281/zenodo.3931288

A BibTeX entry for LaTeX users is:

```
  @Manual{,
    title = {{graphsim}: Simulate Expression Data from 'igraph' Networks },
    author = {S. Thomas Kelly and Michael A. Black},
    year = {2020},
    note = {R package version R package version 1.0.2.},
    url = {https://github.com/TomKellyGenetics/graphsim},
    doi = {10.5281/zenodo.3931288},
  }
```

Please also cite the publication describing use of this package where appropriate.

>Kelly, S.T. and Black, M.A. (2020). graphsim: An R package for simulating gene
>expression data from graph structures of biological pathways.
>_Journal of Open Source Software_, **5**(51), 2161, https://doi.org/10.21105/joss.02161


```
  @article{Kelly2020joss02161,
    doi = {10.21105/joss.02161},
    url ={https://doi.org/10.21105/joss.02161},
    year = {2020},
    publisher = {The Open Journal},
    volume = {5},
    number = {51},
    pages = {2161},
    author = {S. Thomas Kelly and Michael A. Black},
    title = {graphsim: An R package for simulating gene expression data from graph structures of biological pathways},
    journal = {Journal of Open Source Software} }
```

This article is also avaliable as a preprint.

>S. Thomas Kelly, Michael A. Black (2020)
> graphsim: An R package for simulating gene expression data from graph structures of biological pathways
> bioRxiv 2020.03.02.972471; doi:https://doi.org/10.1101/2020.03.02.972471

```
@article {Kelly2020.03.02.972471,
	author = {Kelly, S Thomas and Black, Michael A},
	title = {graphsim: An R package for simulating gene expression data from graph structures of biological pathways},
	elocation-id = {2020.03.02.972471},
	year = {2020},
	doi = {10.1101/2020.03.02.972471},
	publisher = {Cold Spring Harbor Laboratory},
	URL = {https://www.biorxiv.org/content/early/2020/03/04/2020.03.02.972471},
	eprint = {https://www.biorxiv.org/content/early/2020/03/04/2020.03.02.972471.full.pdf},
	journal = {bioRxiv}
}
```

--------------------------------------------------

## Contributions and Bug Reports

Please submit [issues](https://github.com/TomKellyGenetics/graphsim/issues) on GitHub to report
problems or suggest features. [Pull requests](https://github.com/TomKellyGenetics/graphsim/pulls)
to the `dev` branch on GitHub are also welcome to add features or correct problems. Please see
the [contributor guide](https://github.com/TomKellyGenetics/graphsim/blob/master/CONTRIBUTING.md) for more details.


# graphsim 1.0.2

Updates maintainer contact details.

- resolves vignette formatting #11 

- passes updated CRAN checks (links updated)

# graphsim 1.0.1

* Update citation to reflect acceptance at JOSS

* Update documentation (package help page, links and cross-references)

* Critical changes to vignettes to reduce build time (required for regular CRAN checks)

# graphsim 1.0.0

* Major stable release: note changes to results are possible! (legacy code should run without breaking)

* Expanded documentation and examples (consolidate into fewer vignettes for clarity)

* Resolves errors handling inhibiting edges

* Efficiently compute a state matrix from a vector of edge properties from paths

* Enables passing "sd" (standard deviation) to alter covariance of Sigma matrix

* Adds methods for computing using Laplacian matrices

* Adds function to compute simulations directly from an adjacency matrix

* Migrates computing states to sigma (these matrices include inhibitions)

# graphsim 0.1.1

* Initial CRAN release

* Unit testing for all functions

* Checking for compatible inputs

* Passing layout parameters to plotting function

# graphsim 0.1.1

* Full vignette of biological pathway

* Merge plot_directed from function from plot.igraph to avoid github dependencies

* Added example data for reactome pathways

# graphsim 0.1.0

* Initial version of the package for generating simulated data
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
opening an issue or emailing the project lead (<tom.kelly@riken.jp>).

This Code of Conduct is adapted from the Contributor Covenant 
(http://contributor-covenant.org), version 1.0.0, available at 
http://contributor-covenant.org/version/1/0/0/# Contributing to graphsim development

Thank you for helping to make this package better. We value all contributions
and rely on your feedback to identify problems and use cases.

## TL;DR

Send your PR! Thanks!

## More Details

You want to contribute? Awesome! Small changes, like fixing typos in
documentation are completely fine and also most welcome. For bigger
changes, we suggest that you open an issue before you start coding, so that
we can maximize the probability that we can successfully merge in your
code.

The goal of this guide is to help you get up and contributing to graphsim as 
quickly as possible. The guide is divided into two main pieces:
  
  1. Filing a bug report or feature request in an issue.
  
  2. Suggesting a change via a pull request.

Please note that graphsim is released with a [Contributor Code of Conduct](CODE_OF_CONDUCT.md).
By contributing to this project,  you agree to abide by its terms.

## Issues

When filing an issue, the most important thing is to include a minimal 
reproducible example so that we can quickly verify the problem, and then figure 
out how to fix it. There are three things you need to include to make your 
example reproducible: required packages, data, code.

1.  **Packages** should be loaded at the top of the script, so it's easy to
    see which ones the example needs.
  
1.  The easiest way to include **data** is to use `dput()` to generate the R code 
    to recreate it. For example, to recreate the `mtcars` dataset in R,
    I'd perform the following steps:
  
  1. Run `dput(mtcars)` in R
2. Copy the output
3. In my reproducible script, type `mtcars <- ` then paste.

But even better is if you can create a `data.frame()` with just a handful
of rows and columns that still illustrates the problem.

1.  Spend a little bit of time ensuring that your **code** is easy for others to
read:
  
  * make sure you've used spaces and your variable names are concise, but
      informative
  
    * use comments to indicate where your problem lies
  
    * do your best to remove everything that is not related to the problem.  
     The shorter your code is, the easier it is to understand.

You can check you have actually made a reproducible example by starting up a 
fresh R session and pasting your script in.

(Unless you've been specifically asked for it, please don't include the output 
of `sessionInfo()`.)

## Pull requests

To contribute a change to graphsim, you follow these steps:

1. Create a branch in git and make your changes.
1. Push branch to github and issue pull request (PR).
1. Discuss the pull request.
1. Iterate until either we accept the PR or decide that it's not
a good fit for graphsim.

Each of these steps are described in more detail below. This might feel 
overwhelming the first time you get set up, but it gets easier with practice. 
If you get stuck at any point, please reach out for help on the [graphsim-dev](https://groups.google.com/forum/#!forum/graphsim-dev) mailing list.
                                                                                
If you're not familiar with git or github, please start by reading <http://r-pkgs.had.co.nz/git.html>

<!--
* [ ] Motivate the change in one paragraph, and include it in NEWS.
      In parentheses, reference your github user name and this issue:
      `(@hadley, #1234)`
* [ ] Check pull request only includes relevant changes.
* [ ] Use the [official style](http://adv-r.had.co.nz/Style.html).
* [ ] Update documentation and re-run roxygen2
* [ ] Add test, if bug in non-graphical function
* [ ] Add visual test, if bug in graphical function
* [ ] Add minimal example, if new graphical feature

See http://docs.graphsim.org/dev/vignettes/development.html for more details.
--->

Pull requests will be evaluated against a seven point checklist:

1.  __Motivation__. Your pull request should clearly and concisely motivate the
    need for change. Unfortunately I am busy with other projects
    these days, so you need to describe the problem and show
    how your pull request solves it as concisely as possible.

    Also include this motivation in `NEWS` so that when a new release of
    graphsim comes out it's easy for users to see what's changed. Add your
    item at the top of the file and use markdown for formatting. The
    news item should end with `(@yourGithubUsername, #the_issue_number)`.

1.  __Only related changes__. Before you submit your pull request, please
    check to make sure that you haven't accidentally included any unrelated
    changes. These make it harder to see exactly what's changed, and to
    evaluate any unexpected side effects.

    Each PR corresponds to a git branch, so if you expect to submit
    multiple changes make sure to create multiple branches. If you have
    multiple changes that depend on each other, start with the first one
    and don't submit any others until the first one has been processed.
                                                                              
1.  __Use graphsim coding style__. Please follow the
[official tidyverse style](http://style.tidyverse.org). Maintaining
a consistent style across the whole code base makes it much easier to
jump into the code. If you're modifying existing graphsim code that
doesn't follow the style guide, a separate pull request to fix the
style would be greatly appreciated.
                                                                              
1.  If you're adding new parameters or a new function, you'll also need
to document them with [roxygen](https://github.com/klutometis/roxygen).
Make sure to re-run `devtools::document()` on the code before submitting.
                                                                              
Currently, graphsim uses the development version of roxygen2, which you
can get with `install_github("klutometis/roxygen")`. This will be
available on CRAN in the near future.
                                                                              
1.  If fixing a bug or adding a new feature to a non-graphical function,
please add a [testthat](https://github.com/r-lib/testthat) unit test.
                                                                              
1.  If fixing a bug in the visual output, please add a visual test.
(Instructions to follow soon)
                                                                              
1.  If you're adding a new graphical feature, please add a short example
    to the appropriate function.

This seems like a lot of work but don't worry if your pull request isn't perfect.
It's a learning process and members of the graphsim team will be on hand to help you
out. A pull request ("PR") is a process, and unless you've submitted a few in the
past it's unlikely that your pull request will be accepted as is. All PRs require
review and approval from at least one member of the graphsim development team 
before merge.
                                                                              
Please remember that graphsim is package used by other people. 
This means that changing any existing functionality could result in
breaking someone's code (or another package on CRAN). 
Please don't submit pull requests that change existing behaviour. Instead, 
think about how you can add a new feature in a minimally invasive way.

## Making Small Changes

* Please always use the `dev` branch. Choose this branch in your fork. (We
  build the `master` branch from the `dev` branch automatically, to make
  sure that the repo is compatible with the `devtools` R package which uses
  the `master` branch by default.)
* Then look for the file you want to modify.
* Click on the edit symbol (pen) on the upper right corner of the file
  view.
* Make your edits.
* Write a short commit message, less than 65 characters. E.g.  "Fix manual
  page typo" or "Fix degree bug for loops". If needed, elaborate your
  changes below in the "extended description" field.
* Commit your changes.
* Go back to the start page of *your* forked repository. It is at
  `https://github.com/<username>/graphsim`.
* Click on the green button before the branch name to create a pull
  request.
* Click on "Create pull request".
* Provide a more detailed description if you like. Please also indicate
  that you are fine with licensing your contribution under graphsim's license
  (see Legal Stuff below).
* Click on "Create pull request".
* That's it! It is probably a good idea to keep your forked repository
  until the change is accepted into graphsim, in case you need to modify it.
* Now you need to wait for us, unfortunately. Please ping us, if it takes
  long to respond. E.g. a week is considered to be long.
* Once your pull request is accepted, you can delete your forked repository.

## Making More Involved Changes

This is mostly the same as for trivial changes, but you probably want to
edit the sources on your computer, instead of online on Github.

* Open an issue in the issue tracker about the proposed changes.  This is
  not required for smaller things, but I suggest you do it for others. Just
  in case somebody is already working on the same thing, or it is something
  we don't want in graphsim.
* Fork the repository, and clone it to the machine you'll work on.
* We usually build graphsim on OSX, so the `dev` branch is usally fine on
  that platform. It might have problems on other systems. If this happens,
  please open an issue and tell us.
* Make sure you work on the `dev` branch.
* Once ready with your changes, build graphsim, and run the tests. If you use
  the `devtools` package, this (assuming you are in the right directory)
  means running:

  ```R
  library("devtools")
  check()
  install()
  build()
  test()
  ```

* Submit your pull request.
* Now you need to wait for us, unfortunately. Please ping us,
  by email or mentioning the maintainer's username on GitHub
  if it takes longer than a week or so to respond. 

## Writing graphsim Code 

Some tips on writing graphsim code. In general, look at how things are done,
and try to do them similarly. (Unless you think they are not done well, in
which case please tell us.)

### Code Formatting

Look at the style (indentation, braces, etc.) of some recently committed
bigger change, and try to mimic that. The code style within graphsim is not
stricly the same, but we want to keep it reasonably similar. If you are 
unsure on this, we can address this when reviewing the Pull Request so
don't worry about it too much.

### Documentation

Please document your new functions using `roxygen`.

### Test Cases

Unless you change something trivial, please consider adding test cases.
This is important! See the files in the `tests/testthat` directory for
examples.

### Ask Us!

In general, if you are not sure about something, please ask! You can
open an issue on Github or contact the maintainer <tom.kelly@riken.jp>.
We to answer publicly so that others can learn from it, too. There
are silly questions, if you're having trouble others probably are too.

## Legal Stuff

This is a pain to deal with, but we can't avoid it, unfortunately.  So,
graphsim is licensed under the "General Public License (GPL) version 3, or
later". If your contribution is bigger than a typo fix, then please
indicate that you are fine with releasing your code/text under these
licenses.  E.g. adding a sentence that reads as "I'm fine with GPL 3"
is perfectly enough.
## Test environments
* ubuntu 14.04 (on travis-ci), R 3.6.2, R 4.0.0, R 4.0.2
* ubuntu 14.04 (on circle-ci), R 3.6.2
* win-builder (devel and release) Windows Server 2008 R2 SP1, R-devel, 32/64 bit
* Windows Server 2012 R2 x64 on appveyor)  R 3.6.2 Patched 
* rhub (release) Ubuntu Linux 16.04 LTS, R-release, GCC
* Fedora (devel) Linux, R-devel, clang, gfortran
* MacOS 8.6.0 R R 3.6.2
* MacOS Mojave 10.14.6 R 3.6.1, R 3.6.2, R 4.0.2
* MacOS Catalina 10.15.5 R 4.0.2

## Local R CMD check results

── R CMD check results ───────────────────────────────────── graphsim 1.0.2 ────
Duration: 57.9s

0 errors ✓ | 0 warnings ✓ | 0 notes ✓

R CMD check succeeded

## Possible issues

Reactome.org urls have been checked and direct to the correct database. These are only needed for documentation.

graphsim/tests/figs is used exclusively for testing plotting functions with "vdiffr".

Link to local file works on GitHub. Links changed to "\\href"" calls

## Release

This package has also been accepted at the Journal of Open-Source Software for peer-review. This release updates citation information to the publication. 

## Vignettes

Vignettes build with knitr::rmarkdown as used for this CRAN package: https://github.com/yixuan/prettydoc/

Building vignettes that contain real biological datasets are important to
demonstrate functionality of the package (requested by reviewers for JOSS).

This means checks take considerable time to build vignettes. Pre-generated HTMLs
are provided and "eval=FALSE" has been added to the Rmarkdown version of examples
that take a long time to run.

Vignettes are pregenerated to preserve headers, table of contents, and HTML style.

Vignettes in Rmarkdown have been marked as eval=FALSE in some cases (to avoid checks timing out on CRAN).

Pre-generated results so that code and results match will be considered in the future.
This is an urgent release to ensure that checks pass on CRAN (without timing out)
and the (published) package is not archived and remains accessible.

## Minor Release

I will maintain this package at my current address <tom.kelly[at]riken.jp>
---
title: 'graphsim: An R package for simulating gene expression data from graph structures of biological pathways'
output:
 rmarkdown::pdf_document:
   fig_crop: no
   keep_md: TRUE
   #keep_tex: TRUE
   fig_caption: yes
# rmarkdown::html_document:
#  fig_crop: no
#  keep_md: TRUE
#  #keep_tex: TRUE
#  fig_caption: yes
tags:
  - R
  - gene-expression
  - simulation
  - genomics
  - pathway
  - network
authors:
  - name: S. Thomas Kelly
    email: "tom.kelly@postgrad.otago.ac.nz, tom.kelly@riken.jp"
    orcid: 0000-0003-3904-6690
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Michael A. Black
    email: mik.black@otago.ac.nz
    orcid: 0000-0003-1174-6054
    affiliation: "1"
affiliations:
 - name: "Department of Biochemistry, University of Otago, PO Box 56, Dunedin 9054, New Zealand \n"
   index: 1
 - name: "RIKEN Center for Integrative Medical Sciences, Suehiro-cho-1-7-22, Tsurumi Ward, Yokohama, Kanagawa 230-0045, Japan"
   index: 2
date: "24 June 2020"
bibliography: paper.bib
header-includes:
  - \usepackage{caption}
---




### Summary
Transcriptomic analysis is used to capture the molecular state of a cell
or sample in many biological and medical applications. In addition to 
identifying alterations in activity at the level of individual genes, 
understanding changes in the gene networks that regulate fundamental
biological mechanisms is also an important objective of molecular 
analysis. As a result, databases that describe biological pathways 
are increasingly used to assist with the interpretation of results
from large-scale genomics studies. Incorporating information from 
biological pathways and gene regulatory networks into a genomic data
analysis is a popular strategy, and there are many methods that provide
this functionality for gene expression data. When developing or comparing
such methods, it is important to gain an accurate assessment of their 
performance. Simulation-based validation studies are frequently used
for this.
This necessitates the use of simulated data that correctly accounts for
pathway relationships and correlations. Here we present a versatile
statistical framework to simulate correlated gene expression data from
biological pathways, by sampling from a multivariate normal distribution
derived from a graph structure. This procedure has been released as the
\texttt{graphsim} R package on CRAN and GitHub (\url{https://github.com/TomKellyGenetics/graphsim})
 and is compatible with any graph structure that can be described using
the `igraph` package. This package allows the simulation of biological
pathways from a graph structure based on a statistical model of gene expression.


Introduction: inference and modelling of biological networks {#sec:intro}
===============================================================================

Network analysis of molecular biological pathways has the potential to
lead to new insights into biology and medical genetics
[@Barabasi2004; @Hu2016]. Since gene expression profiles capture a
consistent signature of the regulatory state of a cell
[@Perou2000; @Ozsolak2011; @Svensson2018], they can be used to analyse
complex molecular states with genome-scale data. However, biological
pathways are often analysed in a reductionist paradigm as amorphous sets
of genes involved in particular functions, despite the fact that the
relationships defined by pathway structure could further inform gene
expression analyses. In many cases, the pathway relationships are
well-defined, experimentally-validated, and are available in public
databases [@Reactome]. As a result, network analysis techniques can
play an important role in furthering our understanding of biological
pathways and aiding in the interpretation of genomics studies.

Gene networks provide insights into how cells are regulated, by mapping
regulatory interactions between target genes and transcription factors,
enhancers, and sites of epigenetic marks or chromatin structures
[@Barabasi2004; @Yamaguchi2007]. Inference using these regulatory
interactions genomic analysis has the potential to radically
expand the range of candidate biological pathways to be further
explored, or to improve the accuracy of bioinformatics and functional
genomic analysis. A number of methods have been developed to
utilise timecourse gene expression data [@Yamaguchi2007; @Arner2015]
using gene regulatory modules in state-space models and recursive vector
autoregressive models [@Hirose2008; @Shimamura2009]. Various approaches
to gene regulation and networks at the genome-wide scale have led to
novel biological insights [@Arner2015; @Komatsu2013], however, inference
of regulatory networks has thus far primarily relied on experimental
validation or resampling-based approaches to estimate the likelihood
of specific network modules being predicted [@Markowetz2007; @Hawe2019].

Simulating datasets that account for pathway structure are of particular
interest for benchmarking regulatory network inference techniques
and methods being developed for genomics data containing complex biological 
interactions [@Schaffter2011; @Saelens2019]. Dynamical models using
differential equations have been employed, such as by GeneNetWeaver
[@Schaffter2011], to generate simulated datasets
specifically for benchmarking gene regulatory network inference techniques.
There is also renewed interest in modelling biological pathways and simulating
data for benchmarking due to the emergence of single-cell genomics
technologies and the growing number of bioinformatics techniques developed
to use this data [@Zappia2017; @Saelens2019]. Packages such as 'splatter'
[@Zappia2017], which uses the gamma-poisson distribution,
have been developed to model single-cell data.
SERGIO [@Dibaeinia2019] and dyngen [@Cannoodt2020] build on
this by adding gene regulatory networks and multimodality
respectively. These methods have been designed based on known
deterministic relationships or synthetic reaction states,
to which stochasticity is then added.
However, it is computationally-intensive to model
these reactions at scale or run many iterations for benchmarking.
In some cases, it is only necessary to model the statistical
variability and "noise" of RNA-Seq data in order to evaluate
methods in the presence of multivariate correlation structures.

There is a need, therefore, for a systematic framework for statistical
modelling and simulation of gene expression data derived from
hypothetical, inferred or known gene networks. Here we present a
package to achieve this, where samples from a multivariate normal
distribution are used to generate normally-distributed log-expression
data, with correlations between genes derived from the structure of the
underlying pathway or gene regulatory network. This methodology enables
simulation of expression profiles that approximate the log-transformed
and normalised data from microarray studies, as well as bulk or single-cell RNA-Seq
experiments. This procedure has been released as the \texttt{graphsim}
package to enable the generation of simulated gene expression datasets
containing pathway relationships from a known underlying network.
These simulated datasets can be used to evaluate various bioinformatics
methodologies, including
statistical and network inference procedures.

Methodology and software {#sec:methods}
===============================================================================

Here we present a procedure to simulate gene expression data with
correlation structure derived from a known graph structure. This
procedure assumes that transcriptomic data have been generated and
follow a log-normal distribution (i.e.,
$log(X_{ij}) \sim MVN({\bf\mu}, \Sigma)$, where ${\bf\mu}$ and $\Sigma$
are the mean vector and variance-covariance matrix respectively, for
gene expression data derived from a biological pathway) after
appropriate normalisation [@Law2014; @Li2015]. Log-normality of gene
expression matches the assumptions of the popular \texttt{limma} package [@limma], which is
often used for the analysis of intensity-based data from gene expression
microarray studies and count-based data from RNA-Seq experiments. This
approach has also been applied for modelling UMI-based count data from
single-cell RNA-Seq experiments in the \texttt{DESCEND} R package [@Wang2018].

In order to simulate transcriptomic data, a pathway is first constructed
as a graph structure, using the \texttt{igraph} R package [@igraph], with the status of
the edge relationships defined (i.e, whether they activate or inhibit
downstream pathway members). [This procedure uses]{style="color: black"}
a graph structure such as that presented in
Figure [1a](#fig:simple_graph:first){reference-type="ref"
reference="fig:simple_graph:first"}. The graph can be defined by an
adjacency matrix, **$A$** (with elements
$A_{ij}$), where $$A_{ij} = 
\begin{cases}
   1   & \mbox{if genes } i \mbox{ and } j \mbox{ are adjacent} \\
   0   & \mbox{otherwise}
\end{cases}$$

A matrix, **$R$**, with elements
[$R_{ij}$]{style="color: black"}, is calculated based on distance (i.e.,
number of edges contained in the shortest path) between nodes, such that
closer nodes are given more weight than more distant nodes, to define
inter-node relationships. A geometrically-decreasing (relative) distance
weighting is used to achieve this:

[ $$R_{ij} = 
\begin{cases}
   1  & \mbox{if genes } i \mbox{ and } j \mbox{ are adjacent} \\
   (\frac{1}{2})^{d_{ij}}  & \mbox{if a path can be found between genes } i \mbox{ and } j \\
   0  & \mbox{if no path exists between genes } i \mbox{ and } j 
\end{cases}$$]{style="color: black"} where $d_{ij}$ is the length of
the shortest path (i.e., minimum number of edges traversed) between
genes (nodes) $i$ and $j$ in graph $G$. Each more distant node is thus
related by $\frac{1}{2}$ compared to the next nearest, as shown in
Figure [2b](#fig:simulation_activating:second){reference-type="ref"
reference="fig:simulation_activating:second"}. An
arithmetically-decreasing (absolute) distance weighting is also
supported in the \texttt{graphsim} R package which implements this procedure: [ $$R_{ij} = 
\begin{cases}
   1  & \mbox{if genes } i \mbox{ and } j \mbox{ are adjacent} \\
   1-\frac{d_{ij}}{diam(G)}   & \mbox{if a path can be found between genes } i \mbox{ and } j \\
   0  & \mbox{if no path exists between genes } i \mbox{ and } j 
\end{cases}$$ ]{style="color: black"}

Assuming a unit variance for each gene, these values can be used to
derive a $\Sigma$ matrix: $$\Sigma_{ij} = 
\begin{cases}
   1  & \mbox{if } i=j \\
   \rho R_{ij}  & \mbox{otherwise}
\end{cases}$$ where $\rho$ is the correlation between adjacent nodes.
Thus covariances between adjacent nodes are assigned by a correlation
parameter ($\rho$) and the remaining off-diagonal values in the matrix
are based on scaling these correlations by the geometrically weighted
relationship matrix (or the nearest positive definite matrix for
$\Sigma$ with negative correlations).\

Computing the nearest positive definite matrix is necessary to ensure
that the variance-covariance matrix can be inverted when used as a
parameter in multivariate normal simulations, particularly when negative
correlations are included for inhibitions (as shown below). Matrices
that cannot be inverted occur rarely with biologically plausible
graph structures but this approach allows for the computation of a
plausible correlation matrix when the given graph structure is
incomplete or contains loops. When required, the nearest positive
definite matrix is computed using the `nearPD` function of
the \texttt{Matrix} R package [@Matrix] to perform
Higham's algorithm [@Higham2002] on variance-covariance matrices.
The \texttt{graphsim} package gives a warning when this occurs.


Illustrations {#sec:illustrations}
===============================================================================

Generating a Graph Structure {#sec:plot_graph}
-------------------------------------------------------------------------------

The graph structure in
Figure [1a](#fig:simple_graph:first){reference-type="ref"
reference="fig:simple_graph:first"} was used to simulate correlated gene
expression data by sampling from a multivariate normal distribution
using the \texttt{mvtnorm} R package [@Genz2009; @mvtnorm]. The graph structure
visualisation in
Figure [1](#fig:simple_graph){reference-type="ref"
reference="fig:simple_graph"} was specifically developed for (directed)
\texttt{igraph} objects in and is available in the \texttt{graphsim} package. The
\texttt{plot\_directed} function enables customisation of plot parameters for
each node or edge, and mixed (directed) edge types for indicating
activation or inhibition. These inhibition links (which occur frequently
in biological pathways) are demonstrated in
Figure [1b](#fig:simple_graph:second){reference-type="ref"
reference="fig:simple_graph:second"}.

\begin{figure}[!htbp]

{\centering \includegraphics[width=.375\linewidth,height=.375\linewidth]{Plotsimple_graph-1} \includegraphics[width=.375\linewidth,height=.375\linewidth]{Plotsimple_graph-2} 

}

\caption{\textbf{Simulated graph structures}. A constructed graph structure used as an example to demonstrate the simulation procedure in Figures 2 and 3. Activating links are denoted by black arrows and inhibiting links by red edges. }\label{fig:simple_graph}
\end{figure}

A graph structure can be generated and plotted using the following
commands in R:


```r
#install packages required (once per computer)
install.packages("graphsim")
```

```r
#load required packages (once per R instance)
library("graphsim")
#load packages for examples
library("igraph"); library("gplots"); library("scales")
```


```r
#generate graph structure
graph_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
   c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph <- graph.edgelist(graph_edges, directed = TRUE)

#generate parameters for inhibitions for each edge in E(graph)
state <- c(1, 1, -1, 1, 1, 1, 1, -1)
#plot graph structure with inhibitions
plot_directed(graph, state=state, layout = layout.kamada.kawai,
  cex.node = 2, cex.arrow = 4, arrow_clip = 0.2)
```

Generating a Simulated Expression Dataset {#sec:graphsim_demo}
-----------------------------------------

The correlation parameter of $\rho = 0.8$ is used to demonstrate the
inter-correlated datasets using a geometrically-generated relationship
matrix (as used for the example in
Figure [2c](#fig:simulation_activating:third){reference-type="ref"
reference="fig:simulation_activating:third"}). This $\Sigma$ matrix was
then used to sample from a multivariate normal distribution such that
each gene had a mean of $0$, standard deviation $1$, and covariance
within the range $[0,1]$ so that the off-diagonal elements of $\Sigma$
represent correlations. This procedure generated a simulated (continuous
normally-distributed) log-expression profile for each node
(Figure [2e](#fig:simulation_activating:fourth){reference-type="ref"
reference="fig:simulation_activating:fourth"}) with a corresponding
correlation structure (Figure [2d](#fig:simulation_activating:fifth){reference-type="ref"
reference="fig:simulation_activating:fifth"}). The simulated correlation
structure closely resembled the expected correlation structure ($\Sigma$
in Figure [2c](#fig:simulation_activating:third){reference-type="ref"
reference="fig:simulation_activating:third"}) even for the relatively
modest sample size ($N=100$) illustrated in
Figure [2](#fig:simulation_activating){reference-type="ref"
reference="fig:simulation_activating"}. Once a gene expression dataset
comprising multiple pathways has been generated (as in
Figure [2e](#fig:simulation_activating:fourth){reference-type="ref"
reference="fig:simulation_activating:fourth"}), it can then be used to
test procedures designed for analysis of empirical gene expression data
(such as those generated by microarrays or RNA-Seq) that have been
normalised on a log-scale.

\begin{figure}[!htbp]

{\centering \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_activating-1} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_activating-2} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_activating-3} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_activating-4} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_activating-5} 

}

\caption{\textbf{Simulating expression from a graph structure}. An example of a graph structure (a) that has been used to derive a relationship matrix (b), $\Sigma$  matrix (c) and correlation structure (d) from the relative distances between the nodes. Non-negative values are coloured white to red from $0$ to $1$ (e). The $\Sigma$ matrix has been used to generate a simulated expression dataset of 100 samples (coloured blue to red from low to high) via sampling from the multivariate normal distribution. Here genes with closer relationships in the pathway structure show a higher correlation between simulated values.}\label{fig:simulation_activating}
\end{figure}

The simulated dataset can be generated using the following code:


```r
#plot relationship matrix
heatmap.2(make_distance_graph(graph, absolute = FALSE),
  scale = "none", trace = "none", col = colorpanel(50, "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))

#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph, cor = 0.8, absolute = FALSE),
  scale = "none", trace = "none", col = colorpanel(50, "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))

#simulate data
expr <- generate_expression(100, graph, cor = 0.8, mean = 0,
  comm = FALSE, dist =TRUE, absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
```



The simulation procedure
(Figure [2](#fig:simulation_activating){reference-type="ref"
reference="fig:simulation_activating"}) can similarly be used for
pathways containing inhibitory links
(Figure [3](#fig:simulation_inhibiting){reference-type="ref"
reference="fig:simulation_inhibiting"}) with several refinements. With
the inhibitory links
(Figure [3a](#fig:simulation_inhibiting:first){reference-type="ref"
reference="fig:simulation_inhibiting:first"}), distances are calculated
in the same manner as before
(Figure [3b](#fig:simulation_inhibiting:second){reference-type="ref"
reference="fig:simulation_inhibiting:second"}) with inhibitions
accounted for by iteratively multiplying downstream nodes by $-1$ to
form modules with negative correlations between them
(Figures [3c
](#fig:simulation_inhibiting:third){reference-type="ref"
reference="fig:simulation_inhibiting:third"}
and [3d](#fig:simulation_inhibiting:fifth){reference-type="ref"
reference="fig:simulation_inhibiting:fifth"}). A multivariate normal
distribution with these negative correlations can be sampled to generate
simulated data
(Figure [3e](#fig:simulation_inhibiting:fourth){reference-type="ref"
reference="fig:simulation_inhibiting:fourth"}).

The following changes are needed to handle inhibitions:


```r
#generate parameters for inhibitions
state <- c(1, 1, -1, 1, 1, 1, 1, -1)
plot_directed(graph, state=state, layout = layout.kamada.kawai,
  cex.node=2, cex.arrow=4, arrow_clip = 0.2)

#plot sigma matrix for inhibitions
heatmap.2(make_sigma_mat_dist_graph(graph, state, cor = 0.8, absolute = FALSE),
  scale = "none", trace = "none", col = colorpanel(50, "blue", "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))

#simulate data for inhibitions
expr <- generate_expression(100, graph, state, cor = 0.8, mean = 0,
  comm = FALSE, dist =TRUE, absolute = FALSE)
```

The simulation procedure is also demonstrated here
(Figure [4](#fig:simulation_smad){reference-type="ref"
reference="fig:simulation_smad"}) on a pathway structure for a known
biological pathway (Reactome pathway R-HSA-2173789): "TGF-$\beta$ receptor
signaling activates SMADs"
(Figure [4a](#fig:simulation_smad:first){reference-type="ref"
reference="fig:simulation_smad:first"}) derived from the Reactome
database version 52 [@Reactome]. Distances are calculated in the same
manner as before
(Figure [4b](#fig:simulation_smad:second){reference-type="ref"
reference="fig:simulation_smad:second"}) producing blocks of correlated
genes
(Figures [4c](#fig:simulation_smad:third){reference-type="ref"
reference="fig:simulation_smad:third"}
and [4d](#fig:simulation_smad:fifth){reference-type="ref"
reference="fig:simulation_smad:fifth"}). This shows that the
multivariate normal distribution can be sampled to generate simulated
data to represent expression with the complexity of a biological pathway
(Figure [4e](#fig:simulation_smad:fourth){reference-type="ref"
reference="fig:simulation_smad:fourth"}). Here *SMAD7* exhibits
negative correlations with the other SMADs consistent with its
functions as an "inhibitor SMAD" which competitively inhibits *SMAD4*.


\begin{figure}[!htbp]

{\centering \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_inhibiting-1} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_inhibiting-2} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_inhibiting-3} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_inhibiting-4} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_inhibiting-5} 

}

\caption{\textbf{Simulating expression from graph structure with inhibitions}. An example of a graph structure (a), that has been used to derive a relationship matrix (b), $\Sigma$ matrix (c), and correlation structure (d), from the relative distances between the nodes. These values are coloured blue to red from $-1$ to $1$ (e). This has been used to generate a simulated expression dataset of 100 samples (coloured blue to red from low to high) via sampling from the multivariate normal distribution. Here the inhibitory relationships between genes are reflected in negatively correlated simulated  values.}\label{fig:simulation_inhibiting}
\end{figure}


We can import the graph structure into R as follows and run simulations as above:


```r
#import graph from data
graph <- identity(TGFBeta_Smad_graph)
#generate parameters for inhibitions
state <- E(graph)$state

plot_directed(graph, state = state, layout = layout.kamada.kawai,
  border.node=alpha("black", 0.75), fill.node="lightblue",
  col.arrow = c(alpha("navyblue", 0.25), alpha("red", 0.25))[state], 
  cex.node = 1.5, cex.label = 0.8, cex.arrow = 2)
```

These simulated datasets can also be used for simulating gene 
expression data within a graph network to test genomic analysis techniques.
Correlation structure can be included in datasets generated
for testing whether true positive genes or samples can be detected
in a sample with the background of complex pathway structure.


\begin{figure}[!htbp]

{\centering \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_smad-1} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_smad-2} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_smad-3} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_smad-4} \includegraphics[width=.750\linewidth,height=.375\linewidth]{Plotsimulation_smad-5} 

}

\caption{\textbf{Simulating expression from a biological pathway graph structure}. The graph structure (a) of a known biological pathway, "TGF-$\beta$ receptor signaling activates SMADs" (R-HSA-2173789), was used to derive a relationship matrix (b), $\Sigma$ matrix (c) and correlation structure (d) from the relative distances between the nodes. These values are coloured blue to red from $-1$ to $1$ (e). This has been used to generate a simulated expression dataset of 100 samples (coloured blue to red from low to high) via sampling from the multivariate normal distribution. Here modules of genes with correlated expression can be clearly discerned.}\label{fig:simulation_smad}
\end{figure}


Summary and discussion {#sec:summary}
===============================================================================

Biological pathways are of fundamental importance to understanding
molecular biology. In order to translate findings from genomics studies
into real-world applications such as improved healthcare, the roles of
genes must be studied in the context of molecular pathways. Here we
present a statistical framework to simulate gene expression from
biological pathways, and provide the \texttt{graphsim} package in R 
to generate these simulated datasets. This approach is versatile and
can be fine-tuned for modelling existing biological pathways or for
testing whether constructed pathways can be detected by other means.
In particular, methods to infer biological pathways and gene regulatory
networks from gene expression data can be tested on simulated datasets
using this framework. The package also enables simulation of complex gene
expression datasets to test how these pathways impact on statistical
analysis of gene expression data using existing methods or novel
statistical methods being developed for gene expression data analysis.
This approach is intended to be applied to bulk gene expression data
but could in principle be adapted to modelling single-cell or
different modalities such as genome-wide epigenetic data.


Computational details {#computational-details .unnumbered .unnumbered}
===============================================================================

Complete examples of code needed to produce the figures in this paper are
available in the Rmarkdown version in the package GitHub repository
 (\url{https://github.com/TomKellyGenetics/graphsim}).
 Further details are available in the vignettes as well.

The results in this paper were obtained using R 4.0.2 with the \texttt{igraph} 1.2.5
\texttt{Matrix} 1.2-17,
\texttt{matrixcalc} 1.0-3, and \texttt{mvtnorm} 1.1-1 packages.
R itself and all dependent packages
used are available from the Comprehensive Archive Network (CRAN) at
\url{https://CRAN.R-project.org}. The \texttt{graphsim} 1.0.0 package
can be installed from CRAN and the issues can  be reported to
the development version on GitHub.
This package is included in the \texttt{igraph.extensions} library on GitHub (\url{https://github.com/TomKellyGenetics/igraph.extensions})
which installs various
tools for \texttt{igraph} analysis. This software is cross-platform and
compatible with installations on Windows, Mac, and Linux operating
systems. 
Updates to the  package (\texttt{graphsim} 1.0.0) will be released on CRAN.
 
Acknowledgements {#acknowledgements .unnumbered .unnumbered}
===============================================================================

This package was developed as part of a PhD research project funded by
the Postgraduate Tassell Scholarship in Cancer Research Scholarship
awarded to STK. We thank members of the Laboratory of Professor Satoru
Miyano at the University of Tokyo, Institute for Medical Science,
Professor Seiya Imoto, Associate Professor Rui Yamaguchi, and Dr Paul
Sheridan (Assistant Professor at Hirosaki University, CSO at Tupac Bio)
for helpful discussions in this field. We also thank Professor Parry
Guilford at the University of Otago, Professor Cristin Print at the
University of Auckland, and Dr Erik Arner at the RIKEN Center for
Integrative Medical Sciences for their excellent advice during this
project.

Author Contributions  {#contributions .unnumbered .unnumbered}
===============================================================================

S.T.K. and M.A.B. conceived of the presented methodology. S.T.K. developed the theory and performed the computations.
M.A.B. provided guidance throughout the project and gave feedback on the package. All authors discussed the package and contributed to the final manuscript.

# References
title: 'Gala: A Python package for galactic dynamics'
tags:
  - Python
  - astronomy
  - dynamics
  - galactic dynamics
  - milky way
authors:
  - name: Adrian M. Price-Whelan
    orcid: 0000-0003-0872-7098
    affiliation: "1, 2" # (Multiple affiliations must be quoted)
  - name: Author Without ORCID
    affiliation: 2
affiliations:
 - name: Lyman Spitzer, Jr. Fellow, Princeton University
   index: 1
 - name: Institution 2
   index: 2
date: 13 August 2017
bibliography: paper.bib

# Optional fields if submitting to a AAS journal too, see this blog post:
# https://blog.joss.theoj.org/2018/12/a-new-collaboration-with-aas-publishing
aas-doi: 10.3847/xxxxx <- update this with the DOI from AAS once you know it.
aas-journal: Astrophysical Journal <- The name of the AAS journal.
---

# Summary

The forces on stars, galaxies, and dark matter under external gravitational
fields lead to the dynamical evolution of structures in the universe. The orbits
of these bodies are therefore key to understanding the formation, history, and
future state of galaxies. The field of "galactic dynamics," which aims to model
the gravitating components of galaxies to study their structure and evolution,
is now well-established, commonly taught, and frequently used in astronomy.
Aside from toy problems and demonstrations, the majority of problems require
efficient numerical tools, many of which require the same base code (e.g., for
performing numerical orbit integration).

``Gala`` is an Astropy-affiliated Python package for galactic dynamics. Python
enables wrapping low-level languages (e.g., C) for speed without losing
flexibility or ease-of-use in the user-interface. The API for ``Gala`` was
designed to provide a class-based and user-friendly interface to fast (C or
Cython-optimized) implementations of common operations such as gravitational
potential and force evaluation, orbit integration, dynamical transformations,
and chaos indicators for nonlinear dynamics. ``Gala`` also relies heavily on and
interfaces well with the implementations of physical units and astronomical
coordinate systems in the ``Astropy`` package [@astropy] (``astropy.units`` and
``astropy.coordinates``).

``Gala`` was designed to be used by both astronomical researchers and by
students in courses on gravitational dynamics or astronomy. It has already been
used in a number of scientific publications [@Pearson:2017] and has also been
used in graduate courses on Galactic dynamics to, e.g., provide interactive
visualizations of textbook material [@Binney:2008]. The combination of speed,
design, and support for Astropy functionality in ``Gala`` will enable exciting
scientific explorations of forthcoming data releases from the *Gaia* mission
[@gaia] by students and experts alike.

# Mathematics

Single dollars ($) are required for inline mathematics e.g. $f(x) = e^{\pi/x}$

Double dollars make self-standing equations:

$$\Theta(x) = \left\{\begin{array}{l}
0\textrm{ if } x < 0\cr
1\textrm{ else}
\end{array}\right.$$


# Citations

Citations to entries in paper.bib should be in
[rMarkdown](http://rmarkdown.rstudio.com/authoring_bibliographies_and_citations.html)
format.

For a quick reference, the following citation commands can be used:
- `@author:2001`  ->  "Author et al. (2001)"
- `[@author:2001]` -> "(Author et al., 2001)"
- `[@author1:2001; @author2:2001]` -> "(Author1 et al., 2001; Author2 et al., 2002)"

# Figures

Figures can be included like this: ![Example figure.](figure.png)

# Acknowledgements

We acknowledge contributions from Brigitta Sipocz, Syrtis Major, and Semyeong
Oh, and support from Kathryn Johnston during the genesis of this project.

# References
---
name: Report Error
about: Create a report to help us improve the package
title: "[BUG REPORT] with <function> on <data>"
labels: bug
assignees: ''

---

We hope that you are finding the package useful.
Thank you for reporting errors so that we can improve it.
Please describe the problem you are having in as much
 detail as you can so that we can fix it and test our solutions.

**Describe the bug**
Please briefly describe your problem and what output you expect. If you have a question,
please check that it has not already been addressed by a previous issue or
asked before. Please do an internet search for the error first and see if it is 
similar to problems others have had.

Please check if the problem is with the `graphsim` package or any of it's dependencies. 
You can run `rlang::trace_back()` in R to check this.

**To Reproduce**
Please include a minimal reprex. The goal of a reprex is to make it as easy as possible
for me to recreate your problem so that I can fix it. If you've never heard of a reprex
before, start by reading <https://www.tidyverse.org/help/#reprex>, and follow the advice
further down the page. Do NOT include session info unless it's explicitly asked for,
or you've used `reprex::reprex(..., si = TRUE)` to hide it away.  

**Expected behavior**
Please describe what you expect to happen

**Screenshots**
If applicable, give examples of the code you've tried or screenshots
of results  to help explain your problem.

**Additional context**
R is a cross-platform language so there is no need to specify which OS,
version or packages are installed. We will ask for it if it is relevant
such as problems installing the package.

  Delete these instructions once you have read them.

---

Describe your problem

```r
# insert reprex here
```
---
name: Feature request
about: Suggest a new functionality for this project
title: "[NEW FEATURE] Please add <functionality> because <reasons>"
labels: enhancement
assignees: ''

---

We are busy but will consider all requests to add features.
Please describe how it will benefit you and other users 
in as much detail as you can so that we can consider it.

**Is your feature request related to a problem? Please describe.**
A clear and concise description of what the problem is. Ex. I'm always frustrated when [...]

**Describe the solution you'd like**
A clear and concise description of what you want the package to be able to do.

**Describe alternatives you've considered**
A clear and concise description of any alternative solutions or features you've considered.

**Additional context**
Add any other context or screenshots about the feature request here.

  Delete these instructions once you have read them.

---

Describe the feature you wish to be added
---
name: New Issue
about: Open an Issue for discussion
title: "[NEW FEATURE] Please add <functionality> because <reasons>"
labels: question
assignees: ''

---

Please briefly describe your problem and what output you expect. If you have a question,
please check that it has not already been addressed by a previous issue or
asked before.

Please include a minimal reprex. The goal of a reprex is to make it as easy as possible
for me to recreate your problem so that I can fix it. If you've never heard of a reprex
before, start by reading <https://www.tidyverse.org/help/#reprex>, and follow the advice
further down the page. Do NOT include session info unless it's explicitly asked for,
or you've used `reprex::reprex(..., si = TRUE)` to hide it away.  
  
  Delete these instructions once you have read them.

---
  
  Brief description of the problem

```r
# insert reprex here
```