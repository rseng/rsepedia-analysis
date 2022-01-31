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
```---
title: 'graphsim: An R package for simulating gene expression data from graph structures of biological pathways'
output:
  rmarkdown::pdf_document:
    fig_crop: no
    #keep_md: TRUE
    keep_tex: TRUE
    fig_caption: yes
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
 - name: "Department of Biochemistry, University of Otago, PO Box 56, Dunedin 9054, New Zealand"
   index: 1
 - name: "RIKEN Center for Integrative Medical Sciences, Suehiro-cho-1-7-22, Tsurumi Ward, Yokohama, Kanagawa 230-0045, Japan"
   index: 2
date: "`r  format(Sys.time(), '%d %B %Y')`"
bibliography: paper.bib
header-includes:
  - \usepackage{caption}
---


```{r, include = FALSE}
knitr::opts_chunk$set(collapse = TRUE, comment = "#>", width = 68)
knitr::opts_chunk$set(fig.cap = "", fig.path = "Plot")
knitr::opts_chunk$set(fig.pos = "!htbp")
options(width = 68, cli.unicode = FALSE, cli.width = 68)
#par(mai=c(2.82, 2.82, 0.82, 0.82)-0.82)
par(mar=c(7, 10, 4, 2) + 0.1)
par(bty="o")
#captioner::captioner(prefix = "Fig.")
```

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

A graph structure can be generated and plotted using the following
commands in R:

```{r, eval = FALSE}
#install packages required (once per computer)
install.packages("graphsim")
```
```{r, warning=FALSE, results='hide', message=FALSE}
#load required packages (once per R instance)
library("graphsim")
#load packages for examples
library("igraph"); library("gplots"); library("scales")
```

```{r simple_graph_hide, fig.align='center', fig.show='hold', fig.width='1.0\\linewidth', fig.height='1.0\\linewidth', out.height='.375\\linewidth', out.width='.375\\linewidth', fig.retina=10, fig.margin = FALSE, fig.ncol = 2, eval=FALSE}
#generate graph structure
graph_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
   c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph <- graph.edgelist(graph_edges, directed = TRUE)

#plot graph structure (Figure 1a)
plot_directed(graph, state ="activating", layout = layout.kamada.kawai,
  cex.node = 2, cex.arrow = 4, arrow_clip = 0.2)

#generate parameters for inhibitions for each edge in E(graph)
state <- c(1, 1, -1, 1, 1, 1, 1, -1)
#plot graph structure with inhibitions (Figure 1b)
plot_directed(graph, state=state, layout = layout.kamada.kawai,
  cex.node = 2, cex.arrow = 4, arrow_clip = 0.2)
```

```{r simple_graph, fig.align='center', fig.show='hold', fig.width='1.0\\linewidth', fig.height='1.0\\linewidth', out.height='.375\\linewidth', out.width='.375\\linewidth', fig.retina=10, fig.margin = FALSE, fig.ncol = 2, fig.cap = '\\textbf{Simulated graph structures}. A constructed graph structure used as an example to demonstrate the simulation procedure in Figures 2 and 3. Activating links are denoted by black arrows and inhibiting links by red edges.', echo=FALSE}
#generate graph structure
graph_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
   c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph <- graph.edgelist(graph_edges, directed = TRUE)

#plot graph structure (Figure 1a)
plot_directed(graph, state ="activating", layout = layout.kamada.kawai,
  cex.node = 2, cex.arrow = 4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1, line=3.5, at=0.05, adj=0.5, cex=1.75)
box()

#generate parameters for inhibitions for each edge in E(graph)
state <- c(1, 1, -1, 1, 1, 1, 1, -1)
#plot graph structure with inhibitions (Figure 1b)
plot_directed(graph, state=state, layout = layout.kamada.kawai,
  cex.node = 2, cex.arrow = 4, arrow_clip = 0.2)
mtext(text = "(b) Inhibiting pathway structure", side=1, line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
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

The simulated dataset can be generated using the following code:

```{r, include=FALSE, eval=FALSE}
graphics::layout(matrix(c(1:5, 5), nrow = 3, byrow = TRUE))
```

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_activating_hide, fig.align='center', fig.show='hold', fig.width='1.0\\linewidth', fig.height='1.0\\linewidth', out.height='.375\\linewidth', out.width='.375\\linewidth', fig.retina=10, fig.margin = FALSE, fig.ncol = 2, eval=FALSE, warning=FALSE}
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

```{r simulation_activating, fig.align='center', fig.show='hold', fig.width='1.0\\linewidth', fig.height='1.0\\linewidth', out.height='.375\\linewidth', out.width='.375\\linewidth', fig.retina=10, fig.margin = FALSE, fig.ncol = 2, fig.cap = '\\textbf{Simulating expression from a graph structure}. An example of a graph structure (a) that has been used to derive a relationship matrix (b), $\\Sigma$  matrix (c) and correlation structure (d) from the relative distances between the nodes. Non-negative values are coloured white to red from $0$ to $1$ (e). The $\\Sigma$ matrix has been used to generate a simulated expression dataset of 100 samples (coloured blue to red from low to high) via sampling from the multivariate normal distribution. Here genes with closer relationships in the pathway structure show a higher correlation between simulated values.', warning=FALSE, echo=FALSE}
# activating graph
state <- rep(1, length(E(graph)))
plot_directed(graph, state=state, layout = layout.kamada.kawai,
  cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1, line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph, absolute = FALSE),
  scale = "none", trace = "none", key = FALSE,
  col = colorpanel(50, "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
mtext(text = "(b) Relationship matrix", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph, cor = 0.8, absolute = FALSE),
scale = "none", trace = "none", key = FALSE,
col = colorpanel(50, "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#simulate data
expr <- generate_expression(100, graph, cor = 0.8, mean = 0,
comm = FALSE, dist =TRUE, absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none", key = FALSE, 
          col = colorpanel(50, "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
mtext(text = "(d) Simulated correlation", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", key = FALSE, col = bluered(50),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
```

```{r, include=FALSE, warning=FALSE}
pdf("Plotsimulation_activating-5.pdf", width = 12.14, height = 6.072)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", key = FALSE, col = bluered(50),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "", mar = c(6, 6), keysize = 1)
mtext(text = "samples", side=1, line=1.5, at=0.55, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=0.45, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1, line=3.5, at=0.5, adj=0.5, cex=1.75)
dev.off()
#system("sed -i.bak '/Plotsimulation_activating-5/s/\\(.*\\)width=[.]415/\1width=.830/g' paper.md && rm paper.md.bak")
#system("sed -i.bak '/Plotsimulation_activating-5/s/\\(.*\\)width=[.]415/\1width=.830/g' paper.tex && rm paper.tex.bak")
#sed -i '/Plot.*-5/s/\(.*\)width=[.]415/\1width=.830/g' paper.md paper.tex
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

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_inhibiting_hide, fig.align='center', fig.show='hold', fig.width='1.0\\linewidth', fig.height='1.0\\linewidth', out.height='.375\\linewidth', out.width='.375\\linewidth', fig.retina=10, fig.margin = FALSE, fig.ncol = 2, eval=FALSE, warning=FALSE}
#generate parameters for inhibitions
state <- c(1, 1, -1, 1, 1, 1, 1, -1)
plot_directed(graph, state=state, layout = layout.kamada.kawai,
  cex.node=2, cex.arrow=4, arrow_clip = 0.2)

#plot relationship matrix
heatmap.2(make_distance_graph(graph, absolute = FALSE),
  scale = "none", trace = "none", col = colorpanel(50, "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))

#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph, state, cor = 0.8, absolute = FALSE),
  scale = "none", trace = "none", col = colorpanel(50, "blue", "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))

#simulated data
expr <- generate_expression(100, graph, state, cor = 0.8, mean = 0,
  comm = FALSE, dist =TRUE, absolute = FALSE)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
```


```{r simulation_inhibiting, fig.align='center', fig.show='hold', fig.width='1.0\\linewidth', fig.height='1.0\\linewidth', out.height='.375\\linewidth', out.width='.375\\linewidth', fig.retina=10, fig.margin = FALSE, fig.ncol = 2, fig.cap = '\\textbf{Simulating expression from graph structure with inhibitions}. An example of a graph structure (a), that has been used to derive a relationship matrix (b), $\\Sigma$ matrix (c), and correlation structure (d), from the relative distances between the nodes (e). These values are coloured blue to red from $-1$ to $1$. This has been used to generate a simulated expression dataset of 100 samples (coloured blue to red from low to high) via sampling from the multivariate normal distribution. Here the inhibitory relationships between genes are reflected in negatively correlated simulated  values.', warning=FALSE, echo=FALSE}
#generate parameters for inhibitions
state <- c(1, 1, -1, 1, 1, 1, 1, -1)
plot_directed(graph, state=state, layout = layout.kamada.kawai,
  cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Inhibiting pathway structure", side=1, line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph, absolute = FALSE),
  scale = "none", trace = "none", key = FALSE, col = colorpanel(50, "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
mtext(text = "(b) Relationship matrix", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph, state, cor = 0.8, absolute = FALSE),
scale = "none", trace = "none", key = FALSE, col = colorpanel(50, "blue", "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph, state, cor = 0.8, mean = 0,
comm = FALSE, dist =TRUE, absolute = FALSE)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none", key = FALSE, col = colorpanel(50, "blue", "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)))
mtext(text = "(d) Simulated correlation", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", key = FALSE, col = bluered(50),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
```

```{r, include=FALSE, warning=FALSE}
pdf("Plotsimulation_inhibiting-5.pdf", width = 12.14, height = 6.072)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", key = FALSE, col = bluered(50),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "", mar = c(6, 6), keysize = 1)
mtext(text = "samples", side=1, line=1.5, at=0.55, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=0.45, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1, line=3.5, at=0.5, adj=0.5, cex=1.75)
dev.off()
#("sed -i.bak '/Plotsimulation_inhibiting-5/s/\\(.*\\)width=[.]415/\1width=.830/g' paper.md && rm paper.md.bak")
#system("sed -i.bak '/Plotsimulation_inhibiting-5/s/\\(.*\\)width=[.]415/\1width=.830/g' paper.tex && rm paper.tex.bak")
#sed -i '/Plot.*-5/s/\(.*\)width=[.]415/\1width=.830/g' paper.md paper.tex
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

We can import the graph structure into R as follows and run simulations as above:

```{r, include = FALSE}
set.seed(9000)
```

```{r simulation_smad_hide, fig.align='center', fig.show='hold', fig.width='1.0\\linewidth', fig.height='1.0\\linewidth', out.height='.375\\linewidth', out.width='.375\\linewidth', fig.retina=10, fig.margin = FALSE, fig.ncol = 2, eval=FALSE, warning=FALSE}
#import graph from data
graph <- identity(TGFBeta_Smad_graph)
#generate parameters for inhibitions
state <- E(graph)$state

plot_directed(graph, state = state, layout = layout.kamada.kawai,
  border.node=alpha("black", 0.75), fill.node="lightblue",
  col.arrow = c(alpha("navyblue", 0.25), alpha("red", 0.25))[state], 
  cex.node = 1.5, cex.label = 0.8, cex.arrow = 2)

#plot relationship matrix
heatmap.2(make_distance_graph(graph, absolute = FALSE),
  scale = "none", trace = "none", col = colorpanel(50, "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")

#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph, state, cor = 0.8, absolute = FALSE),
  scale = "none", trace = "none", col = colorpanel(50, "blue", "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")

#simulated data
expr <- generate_expression(100, graph, state, cor = 0.8,
  mean = 0,comm = FALSE, dist =TRUE, absolute = FALSE)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none", 
          col = colorpanel(50, "blue", "white", "red"),
  colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
```

```{r simulation_smad, fig.align='center', fig.show='hold', fig.width='1.0\\linewidth', fig.height='1.0\\linewidth', out.height='.375\\linewidth', out.width='.375\\linewidth', fig.retina=10, fig.margin = FALSE, fig.ncol = 2, fig.cap = '\\textbf{Simulating expression from a biological pathway graph structure}. The graph structure (a) of a known biological pathway,  "TGF-$\\beta$ receptor signaling activates SMADs" (R-HSA-2173789), was used to derive a relationship matrix (b), $\\Sigma$ matrix (c) and correlation structure (d) from the relative distances between the nodes. These values are coloured blue to red from $-1$ to $1$ (e). This has been used to generate a simulated expression dataset of 100 samples (coloured blue to red from low to high) via sampling from the multivariate normal distribution. Here modules of genes with correlated expression can be clearly discerned.', warning=FALSE, echo=FALSE}
#import graph from data
graph <- identity(TGFBeta_Smad_graph)
#generate parameters for inhibitions
state <- rep(1, length(E(graph))); pathway <- get.edgelist(graph)
state[pathway[,1] %in% c("SMAD6", "SMAD7", "BAMBI", "SMURF1", "SMURF2", "UCHL5",
  "USP15", "UBB", "UBC", "PMEPA1", "PPP1CA", "PPP1CB", "PPP1CC", "PPP1R15A")] <- 2
state[is.na(state)] <- 1
plot_directed(graph, state = state, layout = layout.kamada.kawai,
  border.node=alpha("black", 0.75), fill.node="lightblue",
  col.arrow = c(alpha("navyblue", 0.25), alpha("red", 0.25))[state], 
  cex.node = 1.5, cex.label = 0.8, cex.arrow = 2, 
  sub = expression(paste("(a) TFG-", beta, " activates SMADs")), cex.sub = 1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph, absolute = FALSE),
  scale = "none", trace = "none", key = FALSE, col = colorpanel(50, "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
mtext(text = "(b) Relationship matrix", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph, state, cor = 0.8, absolute = FALSE),
scale = "none", trace = "none", key = FALSE, col = colorpanel(50, "blue", "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph, state, cor = 0.8,
  mean = 0,comm = FALSE, dist =TRUE, absolute = FALSE)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none", key = FALSE, col = colorpanel(50, "blue", "white", "red"),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
mtext(text = "(d) Simulated correlation", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", key = FALSE, col = bluered(50),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
```


```{r, include=FALSE, warning=FALSE}
pdf("Plotsimulation_smad-5.pdf", width = 12.14, height = 6.072)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
colsep = 1:length(V(graph)), rowsep = 1:length(V(graph)), labCol = "", mar = c(6, 6), keysize = 1)
mtext(text = "samples", side=1, line=1.5, at=0.55, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=0.45, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1, line=3.5, at=0.5, adj=0.5, cex=1.75)
dev.off()
#system("sed -i.bak '/Plotsimulation_smad-5/s/\\(.*\\)width=[.]415/\1width=.830/g' paper.md && rm paper.md.bak")
#system("sed -i.bak '/Plotsimulation_smad-5/s/\\(.*\\)width=[.]415/\1width=.830/g' paper.tex && rm paper.tex.bak")
#sed -i '/Plot.*-5/s/\(.*\)width=[.]415/\1width=.830/g' paper.md paper.tex
```

These simulated datasets can also be used for simulating gene 
expression data within a graph network to test genomic analysis techniques.
Correlation structure can be included in datasets generated
for testing whether true positive genes or samples can be detected
in a sample with the background of complex pathway structure.


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
Sheridan (Assistant Professor at Hirosaki University,CSO at Tupac Bio)
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
---
title: "graphsim: simulating continuous data based on network graph structures"
author: "S. Thomas Kelly^1,2^, Michael A. Black^1^ <br> ^1^ Department of Biochemistry, University of Otago, PO Box 56, Dunedin 9054, New Zealand <br> ^2^ RIKEN Center for Integrative Medical Sciences, Suehiro-cho-1-7-22, Tsurumi Ward, Yokohama"
date: "`r  format(Sys.time(), '%A %d %B %Y')`"
output:
  prettydoc::html_pretty:
       theme: cayman
  #html_document:
       #theme: united
       number_sections: true
       toc: true
       toc_depth: 4
       #toc_float: true
       #code_folding: show
       keep_html: true
toc-title: "Table of Contents"
vignette: >
  %\VignetteIndexEntry{Simulating network graph structure in continuous data}
  %\VignetteEngine{R.rsp::asis}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
#knitr::opts_chunk$set(collapse = TRUE, comment = "#>", width = 68)
knitr::opts_chunk$set(fig.cap = "", fig.path = "Plot")
options(width = 68, cli.unicode = FALSE, cli.width = 68)
#par(mai=c(2.82, 2.82, 0.82, 0.82)-0.82)
par(mar=c(7, 10, 4, 2) + 0.1)
par(bty="o")
#captioner::captioner(prefix = "Fig.")
```

**Summary**

Guide on how to simulate gene expression data from pathway structures defined by graphs from `igraph` objects.

**Package**

graphsim 1.0.0

# Overview of graphsim

This package is designed to balance user-friendliness (by providing sensible defaults and built-in functions) and flexbility (many options are available to use as needed). 

If you have problems or feedback, sumbmitting an issue to the the GitHub repository is encouraged. See the DESCRIPTION and README.md for more details on how to suggest changes to the package.


## Motivations

Pathway and graph structures have a wide array of applications. Here we consider the simulation of (log-normalised) gene expression data from a biological pathway. If you have another use for this software you are welcome to apply it to your problem but please bear in mind that it was designed with this application in mind. In principle, normally-distributed continuous data can be generated based on any defined relationships. This package uses the graph structure to define a ∑ covariance matrix and generate simulated data by sampling from a multivariate normal distribution.

Crucially, this allows the simulation of negative correlations based on inhibitory or repressive relationships, as commonly occur in biology <span class="citation">(Barabási and Oltvai 2004)</span>. A custom plotting function `plot_directed` is provided to visualise these relationships using the "state" parameter. This plotting function has a [dedicated vignette](plots_directed.html).

For more details on the background of this package, see the [paper](https://github.com/TomKellyGenetics/graphsim/blob/master/paper/paper.Rmd) included with the package on GitHub. This vignette provides more detail on the code needed to reproduce the figures presented in the manuscript.
  
# Getting Started
  
## Install dependencies

The package can be installed as follows. Run the following command to install the current release from CRAN (recommended).

```{r, eval = FALSE}
#install packages required (once per computer)
install.packages("graphsim")
```

Run the following command to install the development version from GitHub (advanced users). This will import the latest changes to the package so behaviour may be unstable.

```{r, eval = FALSE}
#install stable release
remotes::install_github("TomKellyGenetics", ref = "master")
#install development version
remotes::install_github("TomKellyGenetics", ref = "dev")
```

Once the required packages are installed, we must load the packages required to use the package functions with `igraph` objects and to generate plots <span class="citation">(Csardi and Nepusz 2006)</span>. Here `igraph` is required to create `igraph` objects and `gplots` is required for plotting heatmaps.

```{r, warning=FALSE, results='hide', message=FALSE}
library("igraph")
library("gplots")
library("graphsim")
library("scales")
```

## Set up simulated graphs

Here we set up a simple graph to demonstrate how connections in the graph structure lead to correlations in the final output. We create a simple pathway of 9 genes with various branches.


```{r, out.width = '50%', out.height  = '50%', fig.align='center', fig.height = 6, fig.width = 6, fig.retina=1.5}
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"),c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
plot_directed(graph_structure, layout = layout.kamada.kawai)
```

# Generating  simulated expression data from graph

## Minimal example 

A simulated dataset can be generated with a single command. This is all that is required to get started.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1,message=FALSE, warning=FALSE}
expr <- generate_expression(100, graph_structure, cor = 0.8, mean = 0,
                            comm = FALSE, dist = TRUE, absolute = FALSE)
heatmap.2(expr, scale = "none", trace = "none",
          col = bluered(50), colsep = 1:4, rowsep = 1:4)
```

Here we've generated a simulated dataset of 100 samples with gene expression data for the genes in the graph shown above. We will show below how these are used to demonstrate what computations are being performed to generate this data from the graph structure given.

Various arguments are supported to alter how the simulated datasets are computed. The other functions used below are called internally by `generate_expression` and are not needed to compute the final dataset in the above heatmap plot. See the documentation for details.

## How it works step-by-step

Here we show the data generated by this graph structure. This demonstrates how several of the available options are used to compute the necessary steps.

### Adjacency matrix

The data can be summarised by an "adjacency matrix" where a one (1) is given between a row `i` and column `j` if there is an edge between genes `i` and `j`. Otherise it is a zero (0) for genes that are not connected. For an undirected graph, edges are shown in a symmetical matrix.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
adj_mat <- make_adjmatrix_graph(graph_structure)
heatmap.2(make_adjmatrix_graph(graph_structure),
          scale = "none", trace = "none",
          col = colorpanel(3, "grey75", "white", "blue"),
          colsep = 1:4, rowsep = 1:4)
```

For a directed graph, the adjacency matrix may be assymetric. A non-zero element `adjmat[i,j]` represents the presence or weight of the edge from gene `i` (matrix rows)  to gene `j`  (matrix columns).

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
heatmap.2(make_adjmatrix_graph(graph_structure, directed = TRUE),
          scale = "none", trace = "none",
          col = colorpanel(3, "grey75", "white", "blue"),
          colsep = 1:4, rowsep = 1:4)
```

We can compute the common links between each pair of genes. This shows how many genes are connected to both genes `i` and `j`.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
comm_mat <- make_commonlink_graph(graph_structure)
heatmap.2(make_commonlink_graph(graph_structure),
          scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

This shows how many edges to a shared neighbour these nodes have between them. The diagonal will therefore reflect vertex degree as all edges are counted.

We define commonlinks between each pair of nodes as how many nodes are mutually connected to both
of the nodes in the adjacency matrix (how many paths of length 2 exist between them). Note that this weights towards genes with a higher vertex degree (as does the Laplacian).


The Laplacian matrix has the same dimensions as the adjancency matrix. For undirected graphs it is a symmetric matrix but for directed graphs it is not. It has the same number of rows and columns as the number of nodes. 

The Laplacian matrix is defined as `laplacian[i,j] = vertex_degree(i)` if `i == j` and  `laplacian[i,j] = -1` if `i != j`. The wieghted Laplacian matrix is defined as `laplacian[i,j] = -wieghts(graph)[i,j]`  for the off-diagonal terms.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
laplacian_mat <- make_laplacian_graph(graph_structure)
heatmap.2(make_laplacian_graph(graph_structure),
          scale = "none", trace = "none", 
          col = bluered(50),colsep = 1:4, rowsep = 1:4)
```

As expected, the off-diagonal terms of the Laplacian are negative integer values and the diagonals reflect the vertex degree.

### Distance matrix

To compute the relationships between each gene by "distance" we first compute the shortest paths between each pair of nodes, using Dijkstra's algorithm <span class="citation">(Prim 1957; Dijkstra 1959)</span>. 

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
shortest.paths(graph_structure)
heatmap.2(shortest.paths(graph_structure),
          scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
           thincolsep = 1:4, rowsep = 1:4)
```

Here we plot the number of edges in the shortest paths between each pair of nodes in the graph (as an integrer value). Relative to the "diameter" (length of the longest shortest path between any 2 nodes), we can show which genes are more similar or different based on the graph structure.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
diam <- diameter(graph_structure)
relative_dist <- (1+diam-shortest.paths(graph_structure))/diam
relative_dist
heatmap.2(relative_dist, scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

These relationships are used to create a distance graph relative to the diameter. A relative geometrically decreasing distance is computed as follows. In this case every connected node is weighted in fractions of the diameter.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
dist_mat <- make_distance_graph(graph_structure, absolute = FALSE)
dist_mat
heatmap.2(dist_mat, scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

An arithmetically decreasing distance is computed as follows. In this case every connected node is by the length of their shortest paths relative to the diameter.


```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
dist_mat <- make_distance_graph(graph_structure, absolute = TRUE)
dist_mat
heatmap.2(dist_mat, scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

### Sigma (Σ) matrix

The Sigma (Σ) covariance matrix defines the relationships between the simulated gene distributions. Where the diagonal is one (1), the covariance terms are correlations between each gene. Where possible these are derived from the distance relationships described above. In cases where this is not compatible, the nearest "positive definite" symmetric matrix is computed.

These can be computed directly from an adjacency matrix.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
#sigma from adj mat
sigma_mat <- make_sigma_mat_graph(graph_structure, 0.8)
sigma_mat
heatmap.2(sigma_mat, scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

A commonlink matrix can also be used to compute a Σ matrix.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
#sigma from comm mat
sigma_mat <- make_sigma_mat_graph(graph_structure, 0.8, comm = TRUE)
sigma_mat
heatmap.2(sigma_mat, scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

It is recommended to compute the distance relationships and use these. This is supported with the built-in functions. For instance Σ from the geometrically computed distances.


```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
# sigma from geometric distance matrix
make_sigma_mat_dist_graph(graph_structure, 0.8, absolute = FALSE)
heatmap.2(make_sigma_mat_dist_graph(graph_structure, 0.8, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

Sigma (Σ) can also be computed for arithmetically computed distances.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
# sigma from absolute distance matrix
sigma_mat <- make_sigma_mat_dist_graph(graph_structure, 0.8, absolute = TRUE)
heatmap.2(sigma_mat, scale = "none", trace = "none", 
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

### Simulated expression and observed correlation

Here we generate the final simulated expression dataset. Note that none of the prior steps are required. These are called internalled as needed.

For example, the adjacency matrix is derived to generate the following dataset. Note that the nearest positive definite matrix is required for the Σ matrix in this case.


```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54, warning=FALSE}
#simulate expression data
#adj mat
expr <- generate_expression(100, graph_structure, cor = 0.8,
                            mean = 0, comm = FALSE)
heatmap.2(expr, scale = "none", trace = "none",
          col = bluered(50), colsep = 1:4, rowsep = 1:4)
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

Here we compute a simulated dataset based on common links shared to other nodes.

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
#comm mat
expr <- generate_expression(100, graph_structure, cor = 0.8, mean = 0, comm =T)
heatmap.2(expr, scale = "none", trace = "none",
          col = bluered(50), colsep = 1:4, rowsep = 1:4)
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

Here we use relative distance (relationships are geometrically weighted to the diameter).

```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54}
# relative dist
expr <- generate_expression(100, graph_structure, cor = 0.8,mean = 0,
                            comm = FALSE, dist = TRUE, absolute = FALSE)
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:4, rowsep = 1:4)
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "grey75", "red"),
          colsep = 1:4, rowsep = 1:4)
```

Here we use absolute distance (relationships are arithmetically weighted to the diameter).

```{r, include = FALSE}
set.seed(9000)
```


```{r, out.width = '50%', out.height  = '50%', fig.align='center', dpi=54, warning=FALSE}
#absolute dist
expr <- generate_expression(100, graph_structure, cor = 0.8, mean = 0,
                            comm = FALSE, dist = TRUE, absolute = TRUE) 
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:4, rowsep = 1:4)
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
           col = bluered(50), colsep = 1:4, rowsep = 1:4)
```


## Summary

In summary, we compute the following expression dataset but on these underlying relationships in the graph structure. Here we use geometrically decreasing correlations between more distant nodes in the network. This code is a complete example to compute the following plots.

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_activating_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, cache=TRUE}
# activating graph
state <- rep(1, length(E(graph_structure)))
plot_directed(graph_structure, state=state,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph_structure, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
sigma_mat <- make_sigma_mat_dist_graph(graph_structure,
               cor = 0.8, absolute = FALSE)
heatmap.2(sigma_mat, scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph_structure, cor = 0.8, mean = 0,
          comm = FALSE, dist =TRUE, absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
```

### Inhibiting relationships

Here we simulate the same graph structure with inhibiting edges but passing the `"state"` parameter. This takes one argument for each each to identify which are inhibitory (as documented).

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_inhibiting_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, cache=TRUE}
# activating graph
state <- c(1, 1, -1, 1, 1, 1, 1, -1)
plot_directed(graph_structure, state=state, layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph_structure, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph_structure, state,
                                    cor = 0.8, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph_structure, state,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

# Examples

## Toy examples

Here we give some toy examples to show how the simulations behave in simple cases. This serves to understand how modules within a larger graph will translate to correlations in the final simulated datasets.

### Diverging branches

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_graph_diverging_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, cache=TRUE}
graph_diverging_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_diverging <- graph.edgelist(graph_diverging_edges, directed = TRUE)

# activating graph
state <- rep(1, length(E(graph_diverging)))
plot_directed(graph_diverging, state=state,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph_diverging, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_diverging)),
          rowsep = 1:length(V(graph_diverging)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph_diverging,
                                    cor = 0.8, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_diverging)),
          rowsep = 1:length(V(graph_diverging)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph_diverging,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_diverging)),
          rowsep = 1:length(V(graph_diverging)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none",col = bluered(50),
          colsep = 1:length(V(graph_diverging)),
          rowsep = 1:length(V(graph_diverging)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

#### Inhibiting relationships

Here we simulate the same graph structure with inhibiting edges but passing the `"state"` parameter. This takes one argument for each each to identify which are inhibitory (as documented).

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_graph_diverging_inhibiting_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, cache=TRUE}
# activating graph
state <- c(1, 1, -1)
plot_directed(graph_diverging, state=state, layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph_diverging, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_diverging)),
          rowsep = 1:length(V(graph_diverging)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph_diverging, state,
                                    cor = 0.8, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(graph_diverging)),
          rowsep = 1:length(V(graph_diverging)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph_diverging, state,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(graph_diverging)),
          rowsep = 1:length(V(graph_diverging)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(graph_diverging)),
          rowsep = 1:length(V(graph_diverging)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

### Converging branches

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_graph_converging_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, cache=TRUE}
graph_converging_edges <- rbind(c("C", "E"), c("D", "E"), c("E", "F"))
graph_converging <- graph.edgelist(graph_converging_edges, directed = TRUE)

# activating graph
state <- rep(1, length(E(graph_converging)))
plot_directed(graph_converging, state=state,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph_converging, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_converging)),
          rowsep = 1:length(V(graph_converging)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph_converging, cor = 0.8,
                                    absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_converging)),
          rowsep = 1:length(V(graph_converging)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph_converging,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_converging)),
          rowsep = 1:length(V(graph_converging)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(graph_converging)),
          rowsep = 1:length(V(graph_converging)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

#### Inhibiting relationships

Here we simulate the same graph structure with inhibiting edges but passing the `"state"` parameter. This takes one argument for each each to identify which are inhibitory (as documented).

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_graph_converging_inhibiting_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, cache=TRUE}
# activating graph
state <- c(-1, 1, -1)
plot_directed(graph_converging, state=state,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph_converging, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_converging)),
          rowsep = 1:length(V(graph_converging)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph_converging, state,
                                    cor = 0.8, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(graph_converging)),
          rowsep = 1:length(V(graph_converging)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph_converging, state,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(graph_converging)),
          rowsep = 1:length(V(graph_converging)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(graph_converging)),
          rowsep = 1:length(V(graph_converging)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

### Reconnecting paths

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_graph_reconnecting_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, cache=TRUE}
graph_reconnecting_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"),
                                  c("C", "E"), c("D", "E"), c("E", "F"))
graph_reconnecting <- graph.edgelist(graph_reconnecting_edges, directed = TRUE)

# activating graph
state <- rep(1, length(E(graph_reconnecting)))
plot_directed(graph_reconnecting, state=state,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph_reconnecting, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_reconnecting)),
          rowsep = 1:length(V(graph_reconnecting)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph_reconnecting, cor = 0.8,
                                    absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_reconnecting)),
          rowsep = 1:length(V(graph_reconnecting)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph_reconnecting,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_reconnecting)),
          rowsep = 1:length(V(graph_reconnecting)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(graph_reconnecting)),
          rowsep = 1:length(V(graph_reconnecting)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

#### Inhibiting relationships

Here we simulate the same graph structure with inhibiting edges but passing the `"state"` parameter. This takes one argument for each each to identify which are inhibitory (as documented).

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_graph_reconnecting_inhibiting_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, cache=TRUE}
# activating graph
state <- c(1, 1, -1, -1, 1, 1, 1, 1)
plot_directed(graph_reconnecting, state=state,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2)
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(graph_reconnecting, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(graph_reconnecting)),
          rowsep = 1:length(V(graph_reconnecting)))
mtext(text = "(b) Relationship matrix", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(graph_reconnecting, state,
                                    cor = 0.8, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(graph_reconnecting)),
          rowsep = 1:length(V(graph_reconnecting)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, graph_reconnecting, state,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(graph_reconnecting)),
          rowsep = 1:length(V(graph_reconnecting)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(graph_reconnecting)),
          rowsep = 1:length(V(graph_reconnecting)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

## Empirical examples

Next we demonstrate the simulation procedure based on real biological pathways from the "Reactome" database <span class="citation">(Croft _et al_. 2014)</span>. We can import these from the `data` directory included with this package. These graphs are given for examples and convenience.

### Kinase pathways

The following pathways are treated as all relationships are activating.

#### RAF/MAP kinase cascade

Here we generate simulated data for the RAF/MAP kinase cascade pathway.

```{r,  fig.align='center', out.width="80%", fig.height = 5, fig.width = 5, eval = FALSE, echo = FALSE, warning=FALSE, message = FALSE}
#RAF_MAP_graph <- identity(RAF_MAP_graph)
#plot_directed(graph, col.arrow = "#00A9FF", fill.node = "lightblue")
```

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_RAF_MAP_graph_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, eval=FALSE}
RAF_MAP_graph <- identity(RAF_MAP_graph)

# activating graph
state <- rep(1, length(E(RAF_MAP_graph)))
plot_directed(RAF_MAP_graph, state=state,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2,
              col.arrow = "#00A9FF", fill.node = "lightblue")
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(RAF_MAP_graph, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(RAF_MAP_graph)),
          rowsep = 1:length(V(RAF_MAP_graph)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(RAF_MAP_graph, cor = 0.8,
                                    absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(RAF_MAP_graph)),
          rowsep = 1:length(V(RAF_MAP_graph)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, RAF_MAP_graph,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(RAF_MAP_graph)),
          rowsep = 1:length(V(RAF_MAP_graph)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(RAF_MAP_graph)),
          rowsep = 1:length(V(RAF_MAP_graph)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

#### Phosphoinositide-3-kinase cascade

Here we generate simulated data for the phosphoinositide-3-kinase (Pi3K) cascade pathway.

```{r, include = FALSE}
set.seed(9000)
```


```{r simulation_Pi3K_graph_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, eval=FALSE}
graph <- identity(Pi3K_graph)

# activating graph
state <- rep(1, length(E(Pi3K_graph)))
plot_directed(Pi3K_graph, state=state,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2,
              col.arrow = "#00A9FF", fill.node = "lightblue")
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(Pi3K_graph, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(Pi3K_graph)),
          rowsep = 1:length(V(Pi3K_graph)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(Pi3K_graph, cor = 0.8,
                                    absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(Pi3K_graph)),
          rowsep = 1:length(V(Pi3K_graph)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, Pi3K_graph, 
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = state)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(Pi3K_graph)),
          rowsep = 1:length(V(Pi3K_graph)))
mtext(text = "(d) Simulated correlation", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(Pi3K_graph)),
          rowsep = 1:length(V(Pi3K_graph)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

### Pathways with repression

#### The Pi3K/AKT pathway

Here we generate simulated data for the phosphoinositide-3-kinase activation of Protein kinase B (PKB) cascade (also known as Pi3k/AKT) pathway. States are imported as edge attributes from the imported graph.


```{r,  fig.align='center', out.width="80%", fig.height = 5, fig.width = 5, warning=FALSE, eval = FALSE, echo = FALSE, message = FALSE}
Pi3K_AKT_graph <- identity(Pi3K_AKT_graph)
edge_properties <- E(Pi3K_AKT_graph)$state
plot_directed(Pi3K_AKT_graph, state = edge_properties,
  col.arrow = c(alpha("navyblue", 0.25), alpha("red", 0.25))[edge_properties],
  fill.node = c("lightblue"))
```

```{r, include = FALSE}
set.seed(9000)
```



```{r simulation_Pi3K_AKT_graph_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, eval = FALSE}
Pi3K_AKT_graph <- identity(Pi3K_AKT_graph)
Pi3K_AKT_graph <- simplify(Pi3K_AKT_graph,
      edge.attr.comb = function(x){
         ifelse(any(x %in% list(-1, 2, "inhibiting", "inhibition")), 2, 1)
        })
Pi3K_AKT_graph <- simplify(Pi3K_AKT_graph, edge.attr.comb = "first")
edge_properties <- E(Pi3K_AKT_graph)$state

# activating graph
#state <- rep(1, length(E(Pi3K_AKT_graph)))
plot_directed(Pi3K_AKT_graph, state = edge_properties,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2,
              col.arrow = c(alpha("navyblue", 0.25),
                            alpha("red", 0.25))[edge_properties],
              fill.node = "lightblue")
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(Pi3K_AKT_graph, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))
mtext(text = "(b) Relationship matrix", side=1, line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(Pi3K_AKT_graph, state = edge_properties,
                                    cor = 0.8, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, Pi3K_AKT_graph, cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE, absolute = FALSE,
                            state = edge_properties)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none",
          col = bluered(50), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

#### The TGFβ-Smad pathway

Here we generate simulated data for the TGFβ-Smad gene regulatory pathway with inhibitions known. States are imported as edge attributes from the imported graph.


```{r,  fig.align='center', out.width="80%", fig.height = 5, fig.width = 5}
TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)
edge_properties <- E(TGFBeta_Smad_graph)$state
plot_directed(TGFBeta_Smad_graph, state = edge_properties,
              col.arrow = c(alpha("navyblue", 0.25),
                            alpha("red", 0.25))[edge_properties],
              fill.node = c("lightblue"))
```

```{r, include = FALSE}
set.seed(9000)
```

```{r simulation_TGFBeta_Smad_graph_hide, fig.show='hold', out.width = '50%', out.height  = '50%', fig.align='center', dpi=36, fig.retina=1, fig.margin = FALSE, fig.ncol = 2, warning=FALSE, message=FALSE, eval=FALSE}
TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)
edge_properties <- E(TGFBeta_Smad_graph)$state

# activating graph
plot_directed(TGFBeta_Smad_graph, state = edge_properties,
              layout = layout.kamada.kawai,
              cex.node=2, cex.arrow=4, arrow_clip = 0.2,
              col.arrow = c(alpha("navyblue", 0.25),
                            alpha("red", 0.25))[edge_properties],
              fill.node = "lightblue")
mtext(text = "(a) Activating pathway structure", side=1,
      line=3.5, at=0.075, adj=0.5, cex=1.75)
box()
#plot relationship matrix
heatmap.2(make_distance_graph(TGFBeta_Smad_graph, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"),
          colsep = 1:length(V(TGFBeta_Smad_graph)),
          rowsep = 1:length(V(TGFBeta_Smad_graph)))
mtext(text = "(b) Relationship matrix", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot sigma matrix
heatmap.2(make_sigma_mat_dist_graph(TGFBeta_Smad_graph,
                                    state = edge_properties,
                                    cor = 0.8, absolute = FALSE),
          scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(TGFBeta_Smad_graph)),
          rowsep = 1:length(V(TGFBeta_Smad_graph)))
mtext(text = expression(paste("(c) ", Sigma, " matrix")), side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#simulated data
expr <- generate_expression(100, TGFBeta_Smad_graph,
                            cor = 0.8, mean = 0,
                            comm = FALSE, dist =TRUE,
                            absolute = FALSE, state = edge_properties)
#plot simulated correlations
heatmap.2(cor(t(expr)), scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"),
          colsep = 1:length(V(TGFBeta_Smad_graph)),
          rowsep = 1:length(V(TGFBeta_Smad_graph)))
mtext(text = "(d) Simulated correlation", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)
#plot simulated expression data
heatmap.2(expr, scale = "none", trace = "none", col = bluered(50),
          colsep = 1:length(V(TGFBeta_Smad_graph)),
          rowsep = 1:length(V(TGFBeta_Smad_graph)), labCol = "")
mtext(text = "samples", side=1, line=1.5, at=0.2, adj=0.5, cex=1.5)
mtext(text = "genes", side=4, line=1, at=-0.4, adj=0.5, cex=1.5)
mtext(text = "(e) Simulated expression data (log scale)", side=1,
      line=3.5, at=0, adj=0.5, cex=1.75)

```

# Summary

Here we have demonstrated that simulated datasets can be generated based on graph structures. These can be either constructed networks for modelling purposes or empirical networks such as those from curated biological databases. 
Various parameters are available and described in the documentation. You can alter these parameters in the examples given here to see the impact they have on the final network. It is encouraged to try different parameters and examine the results carefully, in addition to carefully considering which assumptions are appropriate for your investigations. A model or simulation is never correct, it is a tool to test your assumptions and find weaknesses in your technique, consider which conditions could your method struggle with and model these. Pathway structure in particular should be considered in biological datasets as correlations within a pathway can lead to false positive results and confounding.

The intended application for this package is modelling RNA-Seq gene expression data. However, other applications are encouraged, provided that they require multivariate normal simulations based on relationships in graph structures.

----------------------------------------------------------

## Citation

If you use this package, please cite it where appropriate to recognise the efforts of the developers.

```{r echo=TRUE, message=TRUE, warning=FALSE}
citation("graphsim")
```

## Reporting issues

Please see the GitHub repository for reporting problems and requesting changes. Details on how to contribute can be found in the DESCRIPTION and README.

----------------------------------------------------------

# Session info
Here is the output of `sessionInfo()` on the system on which this document was compiled running pandoc 2.1:

```{r, echo = FALSE}
sessionInfo()
```


----------------------------------------------------------

# References

<p>Barabási, A. L., and Oltvai, Z. N.  2004. “Network Biology: Understanding the Cell’s Functional Organization.” <em>Nat Rev Genet</em> 5 (2): 101–13.</p>

<p>Croft, D., Mundo, A.F., Haw, R., Milacic, M., Weiser, J., Wu, G., Caudy, M., et al. 2014. “The Reactome pathway knowledgebase.” Journal Article. <em>Nucleic Acids Res</em> 42 (database issue): D472–D477. <a href="https://doi.org/10.1093/nar/gkt1102">https://doi.org/10.1093/nar/gkt1102</a>.</p>

<p>Csardi, G., and Nepusz, T. 2006. “The Igraph Software Package for Complex Network Research.” <em>InterJournal</em> Complex Systems: 1695. <a href="https://igraph.org/">https://igraph.org/</a>.</p>

<p>Dijkstra, E.W., 1959. “A note on two problems in connexion with graphs.” <em>Numerische mathematik</em> 1 (1): 269–271.</p>

<p>Prim, R.C., 1957. “Shortest connection networks And some generalizations” <em>Bell System Technical Journal</em> 36 (6): 1389-1401.</p>
---
title: "graphsim: directed plots for igraph objects"
author: "S. Thomas Kelly^1,2^, Michael A. Black^1^ <br> ^1^ Department of Biochemistry, University of Otago, PO Box 56, Dunedin 9054, New Zealand <br> ^2^ RIKEN Center for Integrative Medical Sciences, Suehiro-cho-1-7-22, Tsurumi Ward, Yokohama"
date: "`r  format(Sys.time(), '%A %d %B %Y')`"
output:
  prettydoc::html_pretty:
       theme: cayman
  #html_document:
       #theme: united
       number_sections: true
       toc: true
       toc_depth: 4
       #toc_float: true
       #code_folding: show
       keep_html: true
toc-title: "Table of Contents"
vignette: >
  %\VignetteIndexEntry{Directed plots for igraph objects}
  %\VignetteEngine{R.rsp::asis}
  %\VignetteEncoding{UTF-8}
---

```{r, include = FALSE}
#knitr::opts_chunk$set(collapse = TRUE, comment = "#>", width = 68)
knitr::opts_chunk$set(fig.cap = "", fig.path = "Plot")
knitr::opts_chunk$set(fig.align = "center")
options(width = 68, cli.unicode = FALSE, cli.width = 68)
#par(mai=c(2.82, 2.82, 0.82, 0.82)-0.82)
par(mar=c(7, 10, 4, 2) + 0.1)
par(bty="o")
#captioner::captioner(prefix = "Fig.")
```


<!-- <h1 class="title toc-ignore">Using the graphsim package:  Directed plots for igraph objects</h1> -->

<!-- <h4 class="author">S. Thomas Kelly^1,2^, Michael A. Black^2^</h4> **** -->

<!-- ^1^ Department of Biochemistry, University of Otago, PO Box 56, Dunedin 9054, New Zealand -->

<!-- ^2^ RIKEN Center for Integrative Medical Sciences, Suehiro-cho-1-7-22, Tsurumi Ward, Yokohama -->

<!-- <h4 class="date">`r  format(Sys.time(), '%A %d %B %Y')`</h4> -->

<!-- ---------------------------------------------------------- -->

**Summary**

Guide on how to plot directed graphs from `igraph` objects to display activating and repressive states.

**Package**

graphsim 1.0.0

# Directed graph plots

## Motivations

Here we demonstrate the plotting functions that come built-in with `graphsim`, an alternative to the `plot.igraph` provided by the `igraph` package. Here we provide additional functionality for plotting directed graphs. This draws upon many functions provided by `igraph`, including layout settings <span class="citation">(Csardi and Nepusz 2006)</span>.

In particular, graph and network representations in biology often require displaying edge properties <span class="citation">(Barabási and Oltvai 2004)</span>. Here we have the "state" parameter which can be used to differentiate these, and allow us to represent activating and inhibiting or repressing relationships differently. We use different arrowheads as per convention in molecular biology. 

## Getting started

To generate these plots, the following packages must be imported.

```{r, message=F}
library("igraph")
library("scales")
library("graphsim")
```

# Toy example

Here we demonstrate the plot functions on a small toy graph.

## Set up simulated graph

First we generate a graph object using the `igraph` package.

```{r}
graph_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                     c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph <- graph.edgelist(graph_edges, directed = TRUE)
```

## Plotting

We next demonstrate the plotting function for directed graph objects. `plot_directed` with default settings uses the `layout.fruchterman.reingold` as does built-in plotting function `igraph::plot.igraph`. This function provides additional functionality to displaying directed graphs in particular.

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
plot_directed(graph)
```

Here you can see that the plotting function displays a graph in a similar layout to `plot.igraph` with different aesthetic parameters. We suggest that you choose the function that suits your needs and datasets. We demonstrate the features available for `plot_directed` below.

### Custom aesthetics

We support various aesthetic parameters to control the colour and relative size of nodes and edges.

`plot_directed` supports customised layouts:

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
plot_directed(graph, layout = layout.kamada.kawai)
```

In addition, custom colouts are supported:

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
plot_directed(graph, fill.node = "lightblue", border.node = "royalblue")
```

#### Vectorisation of aesthetics

Colours may also be entered as a vector for each node in `V(graph)`:

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
names(V(graph))
colour_vector <- ifelse(names(V(graph)) %in% c("A", "D", "I"), 1, 2)
plot_directed(graph, fill.node = c("lightblue", "grey")[colour_vector], border.node = c("royalblue", "black")[colour_vector])
```

This functionality allows highlighting of particular groups based on known properties of the graph. For examples `V(graph)$type` for bipartite graphs or partitions from Louvain (`igraph::cluster_louvain`) or Leiden (`leiden::leiden`) clustering algorithms.

### Arrow customisation

The `state` parameter controls whether the links are "activating" or "inhibiting". These can denote activation and repression: foe example, positive and negative regulation of genes or kinase and phosphatase activity of proteins. These may be specified globally as either a character string or numeric:

Activating links are displated with any of the following:

- "activating"
- `1`
- `0`

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
plot_directed(graph, state = "activating")
```

Note that activating states can also be specified as follows:

```{r,  fig.align='center', out.width="80%", fig.height = 6, fig.width = 6, fig.retina=6, eval=FALSE}
plot_directed(graph, state = 1)
plot_directed(graph, state = 0)
```

Inhibiting links are displated with any of the following:

- "inhibiting"
- `-1`
- `2`

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
plot_directed(graph, state = "inhibiting")
```

Note that inhibiting states can also be specified as follows:

```{r,  fig.align='center', out.width="80%", fig.height = 6, fig.width = 6, fig.retina=6, eval=FALSE}
plot_directed(graph, state = -1)
plot_directed(graph, state = 2)
```

#### Vectorisation of edge properties

The state parameter may also be applied as a vector to each edge in `E(graph)` respectively.

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
E(graph)
plot_directed(graph, state = c(1, 1, -1, -1, 1, -1, 1, -1))
```

Note that by default, inhibiting relationships are highlighted with a different `col.arrow` value, which can be controlled by the input parameter.

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
edge_properties <- c(1, 1, -1, -1, 1, -1, 1, -1)/2 + 1.5
plot_directed(graph, state = edge_properties, col.arrow = c("darkgreen", "red")[edge_properties])
```

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
edge_properties <- c(1, 1, -1, -1, 1, -1, 1, -1)/2 + 1.5
ggplot_colours <- c("#F8766D", "#CD9600", "#7CAE00", "#00BE67",
                    "#00BFC4", "#00A9FF", "#C77CFF", "#FF61CC")
plot_directed(graph, state = edge_properties,
              col.arrow = ggplot_colours, fill.node =  ggplot_colours)
```

# Empirical examples

Here we demonstrate using the plotting package to display real biological pathways from the "Reactome" database <span class="citation">(Croft _et al_. 2014)</span>. We can import these from the `data` directory included with this package. These graphs are given for examples and convenience. Any empirical data that consists of a list of directed edges can be imported as an igraph object and handled similarly. Below are some demonstrations.

## RAF/MAP kinase cascade

Here we plot the RAF/MAP kinase cascade pathway.


```{r,  fig.align='center', out.width="80%", fig.height = 5, fig.width = 5, warning=FALSE, message = FALSE}
graph <- identity(RAF_MAP_graph)
plot_directed(graph,col.arrow = alpha("#00A9FF", 0.25),
              fill.node = "lightblue",
              layout = layout.kamada.kawai)
```

## Pi3K cascade


Here we plot the phosphoinositide-3-kinase (Pi3K) cascade pathway.

```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
graph <- identity(Pi3K_graph)
plot_directed(graph, col.arrow = alpha("#00A9FF", 0.25),
              fill.node = "lightblue",
              layout = layout.kamada.kawai)
```


## TGFβ-Smad pathway

Here we plot the TGFβ-Smad pathway with inhibitions known. States are imported as edge attributes from the imported graph.


```{r,  fig.align='center', out.width="80%",fig.height = 6, fig.width = 6, fig.retina=1.5}
graph <- identity(TGFBeta_Smad_graph)
edge_properties <- E(graph)$state
plot_directed(graph, state = edge_properties,
              col.arrow = c(alpha("navyblue", 0.25), alpha("red", 0.25))[edge_properties],
              fill.node = c("lightblue"),
              layout = layout.kamada.kawai)
```


----------------------------------------------------------

# Session info
Here is the output of `sessionInfo()` on the system on which this document was compiled running pandoc 2.1:

```{r, echo = FALSE}
sessionInfo()
```


----------------------------------------------------------

# References

<p>Barabási, A. L., and Oltvai, Z. N.  2004. “Network Biology: Understanding the Cell’s Functional Organization.” <em>Nat Rev Genet</em> 5 (2): 101–13.</p>

<p>Croft, D., Mundo, A.F., Haw, R., Milacic, M., Weiser, J., Wu, G., Caudy, M., et al. 2014. “The Reactome pathway knowledgebase.” Journal Article. <em>Nucleic Acids Res</em> 42 (database issue): D472–D477. <a href="https://doi.org/10.1093/nar/gkt1102">https://doi.org/10.1093/nar/gkt1102</a>.</p>

<p>Csardi, G., and Nepusz, T. 2006. “The Igraph Software Package for Complex Network Research.” <em>InterJournal</em> Complex Systems: 1695. <a href="https://igraph.org/">https://igraph.org/</a>.</p>
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pi3k_graph.R
\docType{data}
\name{Pi3K_graph}
\alias{Pi3K_graph}
\title{PI3K Cascade}
\format{
A graph object of 35 vertices and 251 edges:
\describe{
  \item{V}{gene symbol (human)}
  \item{E}{directed relationship for pathway}
  \item{state}{type of relationship (activating or inhibiting) as edge attribute}
}
}
\source{
PathwayCommons \url{https://reactome.org/content/detail/R-HSA-109704}
}
\usage{
Pi3K_graph
}
\description{
Reactome pathway R-HSA-109704 for the interactions in the phosphoinositide-3-kinase cascade
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plot_directed.R
\name{plot_directed}
\alias{plot_directed}
\alias{plot.directed}
\title{Extensions to igraph for Customising plots}
\usage{
plot_directed(
  graph,
  state = NULL,
  labels = NULL,
  layout = layout.fruchterman.reingold,
  cex.node = 1,
  cex.label = 0.75,
  cex.arrow = 1.25,
  cex.main = 0.8,
  cex.sub = 0.8,
  arrow_clip = 0.075,
  pch = 21,
  border.node = "grey33",
  fill.node = "grey66",
  col.label = NULL,
  col.arrow = NULL,
  main = NULL,
  sub = NULL,
  xlab = "",
  ylab = "",
  frame.plot = F
)
}
\arguments{
\item{graph}{An \code{\link[igraph:aaa-igraph-package]{igraph}} object. Must be directed with known states.}

\item{state}{character or integer. Defaults to "activating" if no "state" edge attribute 
found. May be applied a scalar across all edges or as a vector for each edge respectively. 
Accepts non-integer values for weighted edges provided that the sign indicates whether links
 are activating (positive) or inhibiting (negative). May also be entered as text for 
 "activating" or "inhibiting" or as integers for activating (0,1) or inhibiting (-1,2). 
 Compatible with inputs for make_state_matrix or generate_expression_graph in the graphsim 
 package \url{https://github.com/TomKellyGenetics/graphsim}. Vector input is supported}

\item{labels}{character vector. For labels to plot nodes. Defaults to vertex names in 
graph object. Entering "" would yield unlabelled nodes.}

\item{layout}{function. Layout function as selected from \code{\link[igraph:aaa-igraph-package]{layout_}}. 
Defaults to layout.fruchterman.reingold. Alternatives include layout.kamada.kawai, 
layout.reingold.tilford, layout.sugiyama, and layout.davidson.harel. A 2-column 
layout matrix giving x and y co-ordinates of each node can be given.}

\item{cex.node}{numeric. Defaults to 1.}

\item{cex.label}{numeric. Defaults to 0.75.}

\item{cex.arrow}{numeric Defaults to 1.25. May take a scalar applied to all edges 
or a vector with values for each edge respectively.}

\item{cex.main}{numeric. Defaults to 0.8.}

\item{cex.sub}{numeric. Defaults to 0.8.}

\item{arrow_clip}{numeric Defaults to 0.075 (7.5\%).}

\item{pch}{parameter passed to plot. Defaults to 21. Recommends using selecting 
between 21-25 to preserve colour behaviour. Otherwise entire node will inherit 
border.node as it's colour, in which case a light colour is recommended to see labels.}

\item{border.node}{character. Specifies the colours of node border passed to plot.
Defaults to grey33. Applies to whole node shape if pch has only one colour.}

\item{fill.node}{character. Specfies the colours of node fill passed to plot. 
Defaults to grey66.}

\item{col.label}{character. Specfies the colours of node labels passed to plot. 
Defaults to par("fg").}

\item{col.arrow}{character. Specfies the colours of arrows passed to plot. 
Defaults to par("fg").  May take a scalar applied to all edges or a vector
 with colours for each edge respectively.}

\item{main, sub, xlab, ylab}{Plotting parameters to specify plot titles or axes labels}

\item{frame.plot}{logical. Whether to frame plot with a box. Defaults to FALSE.}
}
\value{
base R graphics
}
\description{
Functions to plot_directed or graph structures including customised colours, layout, states, arrows. Uses graphs functions as an extension of \code{\link[igraph:aaa-igraph-package]{igraph}}. Designed for plotting directed graphs.
}
\examples{

# generate example graphs
library("igraph")
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                           c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)

# plots with igraph defaults
plot(graph_structure, layout = layout.fruchterman.reingold)
plot(graph_structure, layout = layout.kamada.kawai)

# plots with scalar states
plot_directed(graph_structure, state="activating")
plot_directed(graph_structure, state="inhibiting")

# plots with vector states
plot_directed(graph_structure, state = c(1, 1, 1, 1, -1, 1, 1, 1))
plot_directed(graph_structure, state = c(1, 1, -1, 1, -1, 1, -1, 1))
plot_directed(graph_structure, state = c(1, 1, -1, 1, 1, 1, 1, -1))

# plots states with graph attributes
E(graph_structure)$state <- 1
plot_directed(graph_structure)
E(graph_structure)$state <- c(1, 1, -1, 1, -1, 1, -1, 1)
plot_directed(graph_structure)

# plot layout customised
plot_directed(graph_structure, state=c(1, 1, -1, 1, -1, 1, -1, 1), layout = layout.kamada.kawai)

}
\seealso{
See also \code{\link[graphsim]{generate_expression}} for computing the simulated data,
\code{\link[graphsim]{make_sigma}} for computing the Sigma (\eqn{\Sigma}) matrix,
\code{\link[graphsim]{make_distance}} for computing distance from a graph object,
\code{\link[graphsim]{make_state}} for resolving inhibiting states.

See also \code{\link[gplots]{heatmap.2}} for plotting matrices.

See also \code{\link[graphsim]{make_laplacian}}, \code{\link[graphsim]{make_commonlink}}, 
or \code{\link[graphsim]{make_adjmatrix}} for computing input matrices.

See also \code{\link[igraph:aaa-igraph-package]{igraph}} for handling graph objects
and \code{\link[igraph:plot.igraph]{plot.igraph}} for base R \code{\link[base]{plot}} methods.

Other graphsim functions: 
\code{\link{generate_expression}()},
\code{\link{make_adjmatrix}},
\code{\link{make_commonlink}},
\code{\link{make_distance}},
\code{\link{make_laplacian}},
\code{\link{make_sigma}},
\code{\link{make_state}}
}
\author{
Tom Kelly \email{tom.kelly@riken.jp}
}
\concept{graph plotting functions}
\concept{graphsim functions}
\keyword{graph}
\keyword{igraph}
\keyword{plot}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/laplacian.R
\name{make_laplacian}
\alias{make_laplacian}
\alias{make_laplacian_adjmat}
\alias{make_laplacian_graph}
\title{Generate Laplacian Matrix}
\usage{
make_laplacian_adjmat(mat, directed = FALSE)

make_laplacian_graph(graph, directed = FALSE)
}
\arguments{
\item{mat}{precomputed adjacency matrix.}

\item{directed}{logical. Whether directed information is passed to the Laplacian matrix.}

\item{graph}{An \code{\link[igraph:aaa-igraph-package]{igraph}} object. May be directed or weighted.}
}
\value{
An Laplacian matrix compatible with generating an expression matrix
}
\description{
Compute the Laplacian matrix of a (directed) \code{\link[igraph:aaa-igraph-package]{igraph}}
structure, preserving node/column/row names (and direction).
}
\examples{

# construct a synthetic graph module
library("igraph") 
graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
# compute Laplacian matrix for toy example
laplacian_matrix <- make_laplacian_graph(graph_test)
laplacian_matrix

# compute Laplacian matrix from adjacency matrix
adjacency_matrix <- make_adjmatrix_graph(graph_test)
laplacian_matrix <- make_laplacian_adjmat(adjacency_matrix)
laplacian_matrix

# construct a synthetic graph network
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
# compute Laplacian matrix for toy network
graph_structure_laplacian_matrix <- make_laplacian_graph(graph_structure)
graph_structure_laplacian_matrix
 
# import graph from package for reactome pathway
# TGF-\eqn{\Beta} receptor signaling activates SMADs (R-HSA-2173789)
TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)

# compute Laplacian matrix for TGF-\eqn{\Beta} receptor signaling activates SMADs
TGFBeta_Smad_laplacian_matrix <- make_laplacian_graph(TGFBeta_Smad_graph)
dim(TGFBeta_Smad_laplacian_matrix)
TGFBeta_Smad_laplacian_matrix[1:12, 1:12]
# visualise matrix
library("gplots")
heatmap.2(TGFBeta_Smad_laplacian_matrix, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))

}
\seealso{
See also \code{\link[graphsim]{generate_expression}} for computing the simulated data,
\code{\link[graphsim]{make_sigma}} for computing the Sigma (\eqn{\Sigma}) matrix,
\code{\link[graphsim]{make_distance}} for computing distance from a graph object,
\code{\link[graphsim]{make_state}} for resolving inhibiting states.

See also \code{\link[graphsim]{plot_directed}} for plotting graphs or 
\code{\link[gplots]{heatmap.2}} for plotting matrices.

See also \code{\link[graphsim]{make_commonlink}}
or \code{\link[graphsim]{make_adjmatrix}} for computing input matrices.

See also \code{\link[igraph:aaa-igraph-package]{igraph}} for handling graph objects.

Other graphsim functions: 
\code{\link{generate_expression}()},
\code{\link{make_adjmatrix}},
\code{\link{make_commonlink}},
\code{\link{make_distance}},
\code{\link{make_sigma}},
\code{\link{make_state}},
\code{\link{plot_directed}()}

Other graph conversion functions: 
\code{\link{make_adjmatrix}},
\code{\link{make_commonlink}}
}
\author{
Tom Kelly \email{tom.kelly@riken.jp}
}
\concept{graph conversion functions}
\concept{graphsim functions}
\keyword{Laplacian}
\keyword{graph}
\keyword{igraph}
\keyword{network}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/commonlink.R
\name{make_commonlink}
\alias{make_commonlink}
\alias{make_commonlink_adjmat}
\alias{make_commonlink_graph}
\title{Generate Common Link Matrix}
\usage{
make_commonlink_adjmat(adj_mat)

make_commonlink_graph(graph, directed = FALSE)
}
\arguments{
\item{adj_mat}{precomputed adjacency matrix.}

\item{graph}{An \code{\link[igraph:aaa-igraph-package]{igraph}} object. May be directed or weighted.}

\item{directed}{logical. Whether directed information is passed to the adjacency matrix.}
}
\value{
An integer matrix of number of links shared between nodes
}
\description{
Compute the common link matrix of a (directed) \code{\link[igraph:aaa-igraph-package]{igraph}}
structure, preserving node / column / row names (and direction). We can compute the common 
links between each pair of nodes. This shows how many nodes are mutually connected to both
of the nodes in the matrix (how many paths of length 2 exist between them).
}
\examples{

# construct a synthetic graph module
library("igraph")
graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)

# compute adjacency matrix for toy example
adjacency_matrix <- make_adjmatrix_graph(graph_test)
# compute nodes with shared edges to a 3rd node
common_link_matrix <- make_commonlink_adjmat(adjacency_matrix)
common_link_matrix

# construct a synthetic graph network
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
# compute adjacency matrix for toy network
graph_structure_adjacency_matrix <- make_adjmatrix_graph(graph_structure)
# compute nodes with shared edges to a 3rd node
graph_structure_common_link_matrix <- make_commonlink_adjmat(graph_structure_adjacency_matrix)
graph_structure_common_link_matrix

# import graph from package for reactome pathway
# TGF-\eqn{\Beta} receptor signaling activates SMADs (R-HSA-2173789)
TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)
# compute nodes with shared edges to a 3rd node 
TGFBeta_Smad_adjacency_matrix <- make_adjmatrix_graph(TGFBeta_Smad_graph)
TGFBeta_Smad_common_link_matrix <- make_commonlink_adjmat(TGFBeta_Smad_adjacency_matrix)
# we show summary statistics as the graph is large
dim(TGFBeta_Smad_common_link_matrix)
TGFBeta_Smad_common_link_matrix[1:12, 1:12]
# visualise matrix
library("gplots")
heatmap.2(TGFBeta_Smad_common_link_matrix, scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))

}
\seealso{
See also \code{\link[graphsim]{generate_expression}} for computing the simulated data,
\code{\link[graphsim]{make_sigma}} for computing the Sigma (\eqn{\Sigma}) matrix,
\code{\link[graphsim]{make_distance}} for computing distance from a graph object,
\code{\link[graphsim]{make_state}} for resolving inhibiting states.

See also \code{\link[graphsim]{plot_directed}} for plotting graphs or 
\code{\link[gplots]{heatmap.2}} for plotting matrices.

See also \code{\link[graphsim]{make_laplacian}}
or \code{\link[graphsim]{make_adjmatrix}} for computing input matrices.

See also \code{\link[igraph:aaa-igraph-package]{igraph}} for handling graph objects.

Other graphsim functions: 
\code{\link{generate_expression}()},
\code{\link{make_adjmatrix}},
\code{\link{make_distance}},
\code{\link{make_laplacian}},
\code{\link{make_sigma}},
\code{\link{make_state}},
\code{\link{plot_directed}()}

Other graph conversion functions: 
\code{\link{make_adjmatrix}},
\code{\link{make_laplacian}}
}
\author{
Tom Kelly \email{tom.kelly@riken.jp}
}
\concept{graph conversion functions}
\concept{graphsim functions}
\keyword{graph}
\keyword{igraph}
\keyword{neighbourhood}
\keyword{network}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/sigma.R
\name{make_sigma}
\alias{make_sigma}
\alias{make_sigma_mat_adjmat}
\alias{make_sigma_mat_comm}
\alias{make_sigma_mat_laplacian}
\alias{make_sigma_mat_graph}
\alias{make_sigma_mat_dist_adjmat}
\alias{make_sigma_mat_dist_graph}
\title{Generate Sigma (\eqn{\Sigma}) Matrix}
\usage{
make_sigma_mat_adjmat(mat, state = NULL, cor = 0.8, sd = 1)

make_sigma_mat_comm(mat, state = NULL, cor = 0.8, sd = 1)

make_sigma_mat_laplacian(mat, state = NULL, cor = 0.8, sd = 1)

make_sigma_mat_graph(
  graph,
  state = NULL,
  cor = 0.8,
  sd = 1,
  comm = FALSE,
  laplacian = FALSE,
  directed = FALSE
)

make_sigma_mat_dist_adjmat(
  mat,
  state = NULL,
  cor = 0.8,
  sd = 1,
  absolute = FALSE
)

make_sigma_mat_dist_graph(
  graph,
  state = NULL,
  cor = 0.8,
  sd = 1,
  absolute = FALSE
)
}
\arguments{
\item{mat}{precomputed adjacency, laplacian, commonlink, or scaled distance matrix (generated by \code{\link[graphsim]{make_distance}}).}

\item{state}{numeric vector. Vector of length E(graph). Sign used to calculate 
state matrix, may be an integer state or inferred directly from expected correlations
for each edge. May be applied a scalar across all edges or as a vector for each edge
respectively. May also be entered as text for "activating" or "inhibiting" or as
integers for activating (0,1) or inhibiting (-1,2). Compatible with inputs for 
\code{\link[graphsim]{plot_directed}}. Also takes a pre-computed state matrix from
\code{\link[graphsim]{make_state}} if applied to the same graph multiple times.}

\item{cor}{numeric. Simulated maximum correlation/covariance of two adjacent nodes.
Default to 0.8.}

\item{sd}{standard deviations of each gene. Defaults to 1. May be entered as a scalar
applying to all genes or a vector with a separate value for each.}

\item{graph}{An \code{\link[igraph:aaa-igraph-package]{igraph}} object. May be directed or weighted.}

\item{comm}{logical whether a common link matrix is used to compute sigma.
Defaults to FALSE (adjacency matrix).}

\item{laplacian}{logical whether a Laplacian matrix is used to compute sigma.
Defaults to FALSE (adjacency matrix).}

\item{directed}{logical. Whether directed information is passed to the distance matrix.}

\item{absolute}{logical. Whether distances are scaled as the absolute difference from
the diameter (maximum possible). Defaults to TRUE. The alternative is to calculate a
relative difference from the diameter for a geometric decay in distance.}
}
\value{
a numeric covariance matrix of values in the range [-1, 1]
}
\description{
Compute the Sigma (\eqn{\Sigma}) matrix from an \code{\link[igraph:aaa-igraph-package]{igraph}} structure 
or pre-computed matrix. These are compatible with \code{\link[mvtnorm:Mvnorm]{rmvnorm}} and
 \code{\link[graphsim]{generate_expression}}.
By default data is generated with a mean of 0 and standard deviation of 1 for 
each gene (with correlations between derived from the graph structure).
Thus where the Sigma (\eqn{\Sigma}) matrix has diagonals of 1 (for the variance of each gene)
then the symmetric non-diagonal terms (for covariance) determine the correlations
between each gene in the output from \code{\link[graphsim]{generate_expression}}.
}
\examples{

# construct a synthetic graph module
library("igraph")
graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
# compute sigma (\eqn{\Sigma}) matrix for toy example
sigma_matrix <- make_sigma_mat_graph(graph_test, cor = 0.8)
sigma_matrix

# compute sigma (\eqn{\Sigma}) matrix  from adjacency matrix for toy example
adjacency_matrix <- make_adjmatrix_graph(graph_test)
sigma_matrix <- make_sigma_mat_adjmat(adjacency_matrix, cor = 0.8)
sigma_matrix

# compute sigma (\eqn{\Sigma}) matrix from shared edges for toy example
common_link_matrix <- make_commonlink_graph(graph_test)
sigma_matrix <- make_sigma_mat_comm(common_link_matrix, cor = 0.8)
sigma_matrix

# compute sigma (\eqn{\Sigma}) matrix from Laplacian for toy example
laplacian_matrix <- make_laplacian_graph(graph_test)
sigma_matrix <- make_sigma_mat_laplacian(laplacian_matrix, cor = 0.8)
sigma_matrix

# compute sigma (\eqn{\Sigma}) matrix from distance matrix for toy example
distance_matrix <- make_distance_graph(graph_test, absolute = FALSE)
sigma_matrix <- make_sigma_mat_dist_adjmat(distance_matrix, cor = 0.8)
sigma_matrix

# compute sigma (\eqn{\Sigma}) matrix from geometric distance directly from toy example graph
sigma_matrix <- make_sigma_mat_dist_graph(graph_test, cor = 0.8)
sigma_matrix

# compute sigma (\eqn{\Sigma}) matrix from absolute distance directly from toy example graph
sigma_matrix <- make_sigma_mat_dist_graph(graph_test, cor = 0.8, absolute = TRUE)
sigma_matrix

# compute sigma (\eqn{\Sigma}) matrix from geometric distance with sd = 2
sigma_matrix <- make_sigma_mat_dist_graph(graph_test, cor = 0.8, sd = 2)
sigma_matrix

# construct a synthetic graph network
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)

# compute sigma (\eqn{\Sigma}) matrix from geometric distance directly from synthetic graph network
sigma_matrix_graph_structure <- make_sigma_mat_dist_graph(graph_structure,
                                                          cor = 0.8, absolute = FALSE)
sigma_matrix_graph_structure
# visualise matrix
library("gplots")
heatmap.2(sigma_matrix_graph_structure, scale = "none", trace = "none",
                     col = colorpanel(50, "white", "red"))

# compute sigma (\eqn{\Sigma}) matrix from geometric distance directly from
# synthetic graph network with inhibitions
edge_state <- c(1, 1, -1, 1, 1, 1, 1, -1)
# pass edge state as a parameter
sigma_matrix_graph_structure_inhib <- make_sigma_mat_dist_graph(graph_structure, 
                                                                state = edge_state,
                                                                cor = 0.8,
                                                                absolute = FALSE)
sigma_matrix_graph_structure_inhib
# visualise matrix
library("gplots")
heatmap.2(sigma_matrix_graph_structure_inhib, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))

# compute sigma (\eqn{\Sigma}) matrix from geometric distance directly from 
# synthetic graph network with inhibitions
E(graph_structure)$state <-  c(1, 1, -1, 1, 1, 1, 1, -1)
# pass edge state as a graph attribute
sigma_matrix_graph_structure_inhib <- make_sigma_mat_dist_graph(graph_structure,
                                                                cor = 0.8,
                                                                absolute = FALSE)
sigma_matrix_graph_structure_inhib
# visualise matrix
library("gplots")
heatmap.2(sigma_matrix_graph_structure_inhib, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))

# import graph from package for reactome pathway
# TGF-\eqn{\Beta} receptor signaling activates SMADs (R-HSA-2173789)
TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)

# compute sigma (\eqn{\Sigma}) matrix from geometric distance directly from TGF-\eqn{\Beta} pathway
TFGBeta_Smad_state <- E(TGFBeta_Smad_graph)$state
table(TFGBeta_Smad_state)
# states are edge attributes
 sigma_matrix_TFGBeta_Smad_inhib <- make_sigma_mat_dist_graph(TGFBeta_Smad_graph,
                                                              cor = 0.8,
                                                              absolute = FALSE)
# visualise matrix
library("gplots")
heatmap.2(sigma_matrix_TFGBeta_Smad_inhib, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))

# compute sigma (\eqn{\Sigma}) matrix from geometric distance directly from TGF-\eqn{\Beta} pathway
TGFBeta_Smad_graph <- remove.edge.attribute(TGFBeta_Smad_graph, "state")
# compute with states removed (all negative)
sigma_matrix_TFGBeta_Smad <- make_sigma_mat_dist_graph(TGFBeta_Smad_graph,
                                                       state = -1,
                                                       cor = 0.8,
                                                       absolute = FALSE)
# visualise matrix
library("gplots")
heatmap.2(sigma_matrix_TFGBeta_Smad, scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))
# compute with states removed (all positive)
sigma_matrix_TFGBeta_Smad <- make_sigma_mat_dist_graph(TGFBeta_Smad_graph,
                                                       state = 1,
                                                       cor = 0.8,
                                                       absolute = FALSE)
# visualise matrix
library("gplots")
heatmap.2(sigma_matrix_TFGBeta_Smad, scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))

#restore edge attributes
TGFBeta_Smad_graph <- set_edge_attr(TGFBeta_Smad_graph, "state",
                                    value = TFGBeta_Smad_state)
TFGBeta_Smad_state <- E(TGFBeta_Smad_graph)$state
# states are edge attributes
 sigma_matrix_TFGBeta_Smad_inhib <- make_sigma_mat_dist_graph(TGFBeta_Smad_graph,
                                                              cor = 0.8,
                                                              absolute = FALSE)
# visualise matrix
library("gplots")
heatmap.2(sigma_matrix_TFGBeta_Smad_inhib, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))

}
\seealso{
See also \code{\link[graphsim]{generate_expression}} for computing the simulated data,
\code{\link[graphsim]{make_distance}} for computing distance from a graph object,
and
\code{\link[graphsim]{make_state}} for resolving inhibiting states.

See also \code{\link[graphsim]{plot_directed}} for plotting graphs or 
\code{\link[gplots]{heatmap.2}} for plotting matrices.

See also \code{\link[graphsim]{make_laplacian}}, \code{\link[graphsim]{make_commonlink}}, 
or \code{\link[graphsim]{make_adjmatrix}} for computing input matrices.

See also \code{\link[igraph:aaa-igraph-package]{igraph}} for handling graph objects.

Other graphsim functions: 
\code{\link{generate_expression}()},
\code{\link{make_adjmatrix}},
\code{\link{make_commonlink}},
\code{\link{make_distance}},
\code{\link{make_laplacian}},
\code{\link{make_state}},
\code{\link{plot_directed}()}

Other generate simulated expression functions: 
\code{\link{generate_expression}()},
\code{\link{make_distance}},
\code{\link{make_state}}
}
\author{
Tom Kelly \email{tom.kelly@riken.jp}
}
\concept{generate simulated expression functions}
\concept{graphsim functions}
\keyword{graph}
\keyword{igraph}
\keyword{mvtnorm}
\keyword{network}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/graphsim-package.R
\docType{package}
\name{graphsim-package}
\alias{graphsim-package}
\alias{graphsim}
\title{The graphsim package}
\description{
graphsim is a package to simulate normalised expression data from networks 
for biological pathways using \sQuote{\code{igraph}} objects and multivariate
normal distributions.
}
\details{
This package provides functions to develop simulated continuous data 
(e.g., gene expression) from a Sigma (\eqn{\Sigma}) covariance matrix derived from a 
graph structure in \sQuote{\code{igraph}} objects. Intended to extend 
\sQuote{\code{mvtnorm}} to take 'igraph' structures rather than sigma 
matrices as input. This allows the use of simulated data that correctly
accounts for pathway relationships and correlations. Here we present
a versatile statistical framework to simulate correlated gene expression
data from biological pathways, by sampling from a multivariate normal
distribution derived from a graph structure. This package allows the
simulation of biologicalpathways from a graph structure based on a
statistical model of gene expression, such as simulation of expression
profiles that of log-transformed and normalised data from microarray
and RNA-Seq data.
}
\section{Introduction}{

This package enables the generation of simulated gene expression datasets 
containing pathway relationships from a known underlying network.
These simulated datasets can be used to evaluate various bioinformatics 
methodologies, including statistical and network inference procedures.

These are computed by 1) resolving inhibitory states to derive a consistent
matrix of positive and negative edges, 2) inferring relationships between
nodes from paths in the graph, 3) weighting these in a Sigma (\eqn{\Sigma}) 
covariance matrix and 4) using this to sample a multivariate normal 
distribution.
}

\section{Getting Started}{

The \code{\link[graphsim]{generate_expression}} function is a wrapper 
around all necessary functions to give a final simulated dataset.

Here we set up an example graph object using the 
\code{\link[igraph:aaa-igraph-package]{igraph}} package.

\preformatted{
library("igraph")
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"),c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
}

Then we can call \code{\link[graphsim]{generate_expression}} to return
the simulated data based on the relationships defined in the graph
structure. Various options are available to fine-tune this.

\preformatted{
expr <- generate_expression(100, graph_structure,
                            cor = 0.8,
                            mean = 0,
                            sd = 1,
                            comm = FALSE,
                            dist = TRUE,
                            absolute = FALSE,
                            laplacian = FALSE)
}

Here we can see the final result. The graph
structure defines the covariance matrix used
by \code{\link[mvtnorm:Mvnorm]{rmvnorm}} to
generate a multivariate distribution.

\preformatted{
dim(expr)

library("gplots")
heatmap.2(expr,
          scale = "none",
          trace = "none",
          col = bluered(50),
          colsep = 1:4,
          rowsep = 1:4)
}

This dataset consists of 9 rows (one for each vertex or gene)
in the graph and 100 columns (one for each sample or observation).

Input with an adjacency matrix is available using the
\code{\link[graphsim:generate_expression]{generate_expression_mat}}
function.
}

\section{Creating Input Data}{

Graph structures can be passed directly from the
\code{\link[igraph:aaa-igraph-package]{igraph}} package.
Using this package, you can create an \sQuote{\code{igraph}}
class object.


\preformatted{
> class(graph_structure)
[1] "igraph"

> graph_structure
IGRAPH ba7fa2f DN-- 9 8 -- 
  + attr: name (v/c)
  + edges from ba7fa2f (vertex names):
    [1] A->C B->C C->D D->E D->F F->G F->I H->I
}

This \sQuote{\code{igraph}} object class can be passed
directly to \code{\link[graphsim:generate_expression]{generate_expression}}
shown above and internal functions described below:
\code{\link[graphsim:make_sigma]{make_sigma_mat_graph}},
\code{\link[graphsim:make_sigma]{make_sigma_mat_dist_graph}},
\code{\link[graphsim:make_distance]{make_distance_graph}},
and
\code{\link[graphsim:make_state]{make_state_matrix}}.

The \sQuote{\code{graphsim}} package also supports various
matrix formats and has functions to handle these.
The following functions will compute matrices from an
\sQuote{\code{igraph}} object class:

\itemize{

\item  \code{\link[graphsim:make_adjmatrix]{make_adjmatrix_graph}}
to derive the adjacency matrix for a graph structure.

\item  \code{\link[graphsim:make_commonlink]{make_commonlink_graph}}
to derive the \sQuote{common link} matrix for a graph structure of
mutually shared neighbours.

\item  \code{\link[graphsim:make_laplacian]{make_laplacian_graph}}
to derive the Laplacian matrix for a graph structure.

} 

The following functions will compute matrices from an
\code{\link[graphsim:make_adjmatrix]{adjacency matrix}}:

\itemize{

\item  \code{\link[graphsim:make_commonlink]{make_commonlink_adjmat}}
to derive the \sQuote{common link} matrix for a graph structure of
mutually shared neighbours.

\item  \code{\link[graphsim:make_laplacian]{make_laplacian_adjmat}}
to derive the Laplacian matrix for a graph structure.

} 

We provide some pre-generate pathways from Reactoem database
for testing and demonstrations:

\itemize{

\item  \code{\link[graphsim:RAF_MAP_graph]{RAF_MAP_graph }}
for the interactions in the \dQuote{RAF/MAP kinase} cascade (17 vertices
 and 121 edges).

\item  \code{\link[graphsim:Pi3K_graph]{Pi3K_graph}}
for the phosphoinositide-3-kinase cascade (35 vertices and 251 edges).

\item  \code{\link[graphsim:Pi3K_AKT_graph]{Pi3K_AKT_graph}}
for the phosphoinositide-3-kinase activation of Protein kinase B
pathway \dQuote{PI3K/AKT activation} (275 vertices and 21106 edges).


\item  \code{\link[graphsim:TGFBeta_Smad_graph]{TGFBeta_Smad_graph}}
for the TGF-\eqn{\beta} receptor signaling activates SMADs
pathway (32 vertices and 173 edges).
} 

Please note that demonstrations on larger graph objects. These
can be called directly from the pakage:

\preformatted{
> graphsim::Pi3K_graph
IGRAPH 21437e3 DN-- 35 251 -- 
  + attr: name (v/c)
  + edges from 21437e3 (vertex names):
     [1] AKT1->AKT2  AKT1->AKT3  AKT1->CASP9 AKT1->CASP9
     [5] AKT1->CASP9 AKT1->FOXO1 AKT1->FOXO1 AKT1->FOXO1
     [9] AKT1->FOXO3 AKT1->FOXO3 AKT1->FOXO3 AKT1->FOXO4
     [13] AKT1->FOXO4 AKT1->FOXO4 AKT1->GSK3B AKT1->GSK3B
     [17] AKT1->GSK3B AKT1->NOS1  AKT1->NOS2  AKT1->NOS3 
     [21] AKT1->PDPK1 AKT2->AKT3  AKT2->CASP9 AKT2->CASP9
     [25] AKT2->CASP9 AKT2->FOXO1 AKT2->FOXO1 AKT2->FOXO1
     [29] AKT2->FOXO3 AKT2->FOXO3 AKT2->FOXO3 AKT2->FOXO4
     + ... omitted several edges
     + ... omitted several edges
}

They can also be imported into R:

\preformatted{
data(Pi3K_graph)
}

You can assign them to your local environment
by calling with from the package:

\preformatted{
graph_object <- identity(Pi3K_graph)
}

You can also change the object class directly
from the package:

\preformatted{
library("igraph")
Pi3K_adjmat <- as_adjacency_matrix(Pi3K_graph)
}

 \code{\link[graphsim:Pi3K_AKT_graph]{Pi3K_AKT_graph}} and 
 \code{\link[graphsim:TGFBeta_Smad_graph]{TGFBeta_Smad_graph}}
 contain graph edge attributes for the \sQuote{state} parameter
 described below.
 
 \preformatted{
 > TGFBeta_Smad_graph
 IGRAPH f3eac04 DN-- 32 173 -- 
   + attr: name (v/c), state (e/n)
   + edges from f3eac04 (vertex names):
     [1] BAMBI ->SMAD7  BAMBI ->TGFB1  BAMBI ->TGFBR1 BAMBI ->TGFBR2
     [5] CBL   ->NEDD8  CBL   ->NEDD8  CBL   ->TGFBR2 CBL   ->TGFBR2
     [9] CBL   ->UBE2M  CBL   ->UBE2M  FKBP1A->TGFB1  FKBP1A->TGFBR1
     [13] FKBP1A->TGFBR2 FURIN ->TGFB1  FURIN ->TGFB1  MTMR4 ->SMAD2 
     [17] MTMR4 ->SMAD2  MTMR4 ->SMAD3  MTMR4 ->SMAD3  NEDD4L->RPS27A
     [21] NEDD4L->SMAD7  NEDD4L->SMURF1 NEDD4L->SMURF2 NEDD4L->TGFB1 
     [25] NEDD4L->TGFBR1 NEDD4L->TGFBR2 NEDD4L->UBA52  NEDD4L->UBB   
     [29] NEDD4L->UBC    NEDD8 ->TGFBR2 NEDD8 ->UBE2M  PMEPA1->SMAD2 
     + ... omitted several edges
     
 > E(TGFBeta_Smad_graph)$state
 [1] 2 2 2 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [32] 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
 [63] 2 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [94] 1 1 1 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2 2
 [125] 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1
 [156] 1 1 1 1 1 1 1 1 1 2 2 2 2 2 2 1 1 1
 
 > states <- E(TGFBeta_Smad_graph)$state
 > table(states)
 states
 1   2 
 103  70 
 }
}

\section{Internal Functions}{


The following functions are used by
\code{\link[graphsim:generate_expression]{generate_expression}}
to compute a simulated dataset. They can be called separately
to summarise the steps used to compute the final data matrix
or for troubleshooting.

\itemize{

\item \code{\link[graphsim:make_sigma]{make_sigma_mat_adjmat}},
\code{\link[graphsim:make_sigma]{make_sigma_mat_comm}}, 
\code{\link[graphsim:make_sigma]{make_sigma_mat_laplacian}}, and
\code{\link[graphsim:make_sigma]{make_sigma_mat_graph}} will
compute a Sigma (\eqn{\Sigma}) covariance matrix from an
adjacency matrix, common link matrix, Laplacian matrix,
or an \sQuote{igraph} object. There are computed as above
and passed to \code{\link[mvtnorm:Mvnorm]{rmvnorm}}.

\item \code{\link[graphsim:make_distance]{make_distance_adjmat}},
\code{\link[graphsim:make_distance]{make_distance_comm}}, 
\code{\link[graphsim:make_distance]{make_distance_laplacian}}, and
\code{\link[graphsim:make_distance]{make_distance_graph}} will
compute a distance matrix of relationships from an
adjacency matrix, common link matrix, Laplacian matrix,
or an \sQuote{igraph} object. There are computed as above
and passed to \code{\link[graphsim:make_sigma]{make_sigma}}.

\item \code{\link[graphsim:make_state]{make_state_matrix}}
will compute a \dQuote{state matrix} resolving positive and
negative correlations from a vector of edge properties. This
is called by \code{\link[graphsim:make_sigma]{make_sigma}}
and \code{\link[graphsim:generate_expression]{generate_expression}}
to ensure that the signs of correlations are consistent.
}
}

\section{Examining Step-by-Step}{


These internal functions can be called to compute steps of
the simulation procedure and examine the results.

1. first we create a graph structure and define the input parameters

\preformatted{
library("igraph")
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"),c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
#sample size
data.n <- 100
#data distributions
data.cor <- 0.75
data.mean <- 3
data.sd <- 1.5
#inhibition states
edge_states <- c(1, 1, -1, -1, 1, 1, 1, 1)
}

2. examine the relationships between the genes.

Here we can see which nodes share an edge:

\preformatted{
> adjacency_matrix <- make_adjmatrix_graph(graph_structure)
> adjacency_matrix
  A C B D E F G I H
A 0 1 0 0 0 0 0 0 0
C 1 0 1 1 0 0 0 0 0
B 0 1 0 0 0 0 0 0 0
D 0 1 0 0 1 1 0 0 0
E 0 0 0 1 0 0 0 0 0
F 0 0 0 1 0 0 1 1 0
G 0 0 0 0 0 1 0 0 0
I 0 0 0 0 0 1 0 0 1
H 0 0 0 0 0 0 0 1 0
}

Here we define a geometrically decreasing series of relationships
between genes based on distance by paths in the graph:

\preformatted{
> relationship_matrix <- make_distance_graph(graph_structure, absolute = FALSE)
> relationship_matrix
  A          C          B          D          E          F          G          I          H
A 1.00000000 0.20000000 0.10000000 0.10000000 0.06666667 0.06666667 0.05000000 0.05000000 0.04000000
C 0.20000000 1.00000000 0.20000000 0.20000000 0.10000000 0.10000000 0.06666667 0.06666667 0.05000000
B 0.10000000 0.20000000 1.00000000 0.10000000 0.06666667 0.06666667 0.05000000 0.05000000 0.04000000
D 0.10000000 0.20000000 0.10000000 1.00000000 0.20000000 0.20000000 0.10000000 0.10000000 0.06666667
E 0.06666667 0.10000000 0.06666667 0.20000000 1.00000000 0.10000000 0.06666667 0.06666667 0.05000000
F 0.06666667 0.10000000 0.06666667 0.20000000 0.10000000 1.00000000 0.20000000 0.20000000 0.10000000
G 0.05000000 0.06666667 0.05000000 0.10000000 0.06666667 0.20000000 1.00000000 0.10000000 0.06666667
I 0.05000000 0.06666667 0.05000000 0.10000000 0.06666667 0.20000000 0.10000000 1.00000000 0.20000000
H 0.04000000 0.05000000 0.04000000 0.06666667 0.05000000 0.10000000 0.06666667 0.20000000 1.00000000
}

Here we can see the resolved edge states through paths in the
adjacency matrix:

\preformatted{
> names(edge_states) <- apply(graph_structure_edges, 1, paste, collapse = "-")
> edge_states
A-C B-C C-D D-E D-F F-G F-I H-I 
1   1  -1  -1   1   1   1   1 
> state_matrix <- make_state_matrix(graph_structure, state = edge_states)
> state_matrix
   A  C  B  D  E  F  G  I  H
A  1  1  1 -1  1 -1 -1 -1 -1
C  1  1  1 -1  1 -1 -1 -1 -1
B  1  1  1 -1  1 -1 -1 -1 -1
D -1 -1 -1  1 -1  1  1  1  1
E  1  1  1 -1  1 -1 -1 -1 -1
F -1 -1 -1  1 -1  1  1  1  1
G -1 -1 -1  1 -1  1  1  1  1
I -1 -1 -1  1 -1  1  1  1  1
H -1 -1 -1  1 -1  1  1  1  1
}

3. define a Sigma (\eqn{\Sigma}) covariance matrix

Here we can see that the signs match the \code{state_matrix}
and the covariance is based on the \code{relationship_matrix}
weighted by the correlation (\code{cor}) and standard
deviation (\code{sd}) parameters.

Note that where \code{sd = 1}, the diagonals will be \code{1}
and the off-diagonal terms will be correlations.

\preformatted{
> sigma_matrix <- make_sigma_mat_dist_graph(
+     graph_structure,
+     state = edge_states,
+     cor = data.cor,
+     sd = data.sd,
+     absolute = FALSE
+ )
> sigma_matrix
   A         C         B        D         E        F         G         I         H
A  2.250000  1.687500  0.843750 -0.84375  0.562500 -0.56250 -0.421875 -0.421875 -0.337500
C  1.687500  2.250000  1.687500 -1.68750  0.843750 -0.84375 -0.562500 -0.562500 -0.421875
B  0.843750  1.687500  2.250000 -0.84375  0.562500 -0.56250 -0.421875 -0.421875 -0.337500
D -0.843750 -1.687500 -0.843750  2.25000 -1.687500  1.68750  0.843750  0.843750  0.562500
E  0.562500  0.843750  0.562500 -1.68750  2.250000 -0.84375 -0.562500 -0.562500 -0.421875
F -0.562500 -0.843750 -0.562500  1.68750 -0.843750  2.25000  1.687500  1.687500  0.843750
G -0.421875 -0.562500 -0.421875  0.84375 -0.562500  1.68750  2.250000  0.843750  0.562500
I -0.421875 -0.562500 -0.421875  0.84375 -0.562500  1.68750  0.843750  2.250000  1.687500
H -0.337500 -0.421875 -0.337500  0.56250 -0.421875  0.84375  0.562500  1.687500  2.250000
}

4. generate an expression dataset using this sigma matrix

We use \code{generate_expression} to compute and expression
dataset, simulated using these parameters:

\preformatted{
> expression_data <- generate_expression(
+     n = data.n,
+     graph_structure,
+     state = edge_states,
+     cor = data.cor,
+     mean = data.mean,
+     sd = data.sd,
+     comm = FALSE,
+     dist = FALSE,
+     absolute = FALSE,
+     laplacian = FALSE
+ )
> dim(expression_data)
[1]   9 100
}

Here we also compute the final observed correlations
in the simulated dataset:

\preformatted{
> cor_data <- cor(t(expression_data))
> dim(cor_data)
[1] 9 9
}

These functions are demonstrated in more detail
in the \href{https://CRAN.R-project.org/package=graphsim/vignettes/simulate_expression.html}{main} vignette.
}

\section{Data Visualization}{


Heatmaps can be used from the \code{\link[gplots:heatmap.2]{gplots}}
package to display these simulated datasets.

\preformatted{
library("gplots")
heatmap.2(adjacency_matrix, scale = "none", trace = "none",
          col = colorpanel(50, "white", "black"), key = FALSE)
          
heatmap.2(relationship_matrix, scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))
          
heatmap.2(state_matrix, scale = "none", trace = "none",
          col = colorpanel(50, "royalblue", "palevioletred"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))

heatmap.2(sigma_matrix, scale = "none", trace = "none",
          col = colorpanel(50, "royalblue", "white", "palevioletred"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))
          
heatmap.2(expression_data, scale = "none", trace = "none",
          col = colorpanel(50, "royalblue", "white", "palevioletred"),
          colsep = 1:length(V(graph_structure)),
         rowsep = 1:length(V(graph_structure)))

heatmap.2(cor_data, scale = "none", trace = "none",
          col = colorpanel(50, "royalblue", "white", "palevioletred"),
          colsep = 1:length(V(graph_structure)),
          rowsep = 1:length(V(graph_structure)))
}

In particular we can see here that the expected correlations
show by the \code{sigma_matrix} are similar to the observed
correlations in the \code{cor_data}.
}

\section{Graph Visualization}{


The \sQuote{graphsim} package comes with a built-in plotting
function to display graph objects. 

\preformatted{
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"),c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
plot_directed(graph_structure, layout = layout.kamada.kawai)
}

This supports the \sQuote{state} parameter to display
activating relationships (with positive correlations)
and inhibiting or repressive relationships (with
negative correlations).

\preformatted{
edge_states <- c(1, 1, -1, -1, 1, -1, 1, -1)
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
plot_directed(graph_structure, state = edge_states,
              col.arrow = c("darkgreen", "red")[edge_states / 2 + 1.5]
              layout = layout.kamada.kawai)
}

These states can also be passed from the \sQuote{state} edge
attribute of the graph object.

\preformatted{
graph_pathway <- identity(TGFBeta_Smad_graph)
edge_properties <- E(graph_pathway)$state
plot_directed(graph_pathway,
              col.arrow = c(alpha("navyblue", 0.25),
                            alpha("red", 0.25))[edge_properties],
              fill.node = c("lightblue"),
              layout = layout.kamada.kawai)
}

This plotting function is demonstrated in more detail
in the \href{https://CRAN.R-project.org/package=graphsim/vignettes/plots_directed.html}{plotting} vignette.
}

\section{Further information}{

  The graphsim package is published in the \emph{Journal of Open Source Software}.
  See the paper here for more details:
   \href{https://doi.org/10.21105/joss.02161}{doi: 10.21105/joss.02161}
  
  The graphsim GitHub repository is here:
  \href{https://github.com/TomKellyGenetics/graphsim}{TomKellyGenetics/graphsim}
  You can find the development version and submit an
  \href{https://github.com/TomKellyGenetics/graphsim/issues/new/choose}{issue}
  if you have questions or comments.
}

\section{Citation}{

  
  To cite package 'graphsim' in publications use:
  
  Kelly, S.T. and Black, M.A. (2020). graphsim: An R package for simulating gene
  expression data from graph structures of biological pathways.
  \emph{Journal of Open Source Software}, \bold{5}(51), 2161,
  \url{https://doi.org/10.21105/joss.02161}
  
  A BibTeX entry for LaTeX users is: \preformatted{
  @article{Kelly2020joss02161,
     doi = {10.21105/joss.02161},
     year = {2020},
     publisher = {The Open Journal},
     volume = {5},
     number = {51},
     pages = {2161},
     author = {S. Thomas Kelly and Michael A. Black},
     title = {graphsim: An R package for simulating gene expression data from graph structures of biological pathways},
     journal = {Journal of Open Source Software} 
   }
 }
}

\seealso{
Publication at \emph{Journal of Open Source Software}:
 
 \itemize{
 \item \url{https://doi.org/10.21105/joss.02161}
 }
 
 GitHub repository:
 
 \itemize{
 \item \url{https://github.com/TomKellyGenetics/graphsim/}
 }
 
 Report bugs:
 
 \itemize{
 \item \url{https://github.com/TomKellyGenetics/graphsim/issues}
 }
 
 Contributions:
 
 \itemize{
 \item \url{https://github.com/TomKellyGenetics/graphsim/blob/master/CONTRIBUTING.md}
 }
}
\author{
\bold{Maintainer}:  Tom Kelly  \email{tom.kelly@riken.jp}

Authors:
\itemize{
\item Tom Kelly (RIKEN IMS) \href{https://orcid.org/0000-0003-3904-6690}{ORCID})
\item Mik Black (Otago University) (\href{https://orcid.org/0000-0003-1174-6054}{ORCID})
}

Reviewers:
\itemize{
\item Cory Brunson (UConn) (\href{https://orcid.org/0000-0003-3126-9494}{ORCID})
\item Robrecht Cannoodt (Ghent University) (\href{https://orcid.org/0000-0003-3641-729X}{ORCID})
}

Editor: Mark Jensen (Frederick National Laboratory for Cancer Research)
}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/raf_mapk_graph.R
\docType{data}
\name{RAF_MAP_graph}
\alias{RAF_MAP_graph}
\title{#' RAF/MAP kinase cascade}
\format{
A graph object of 17 vertices and 121 edges:
\describe{
  \item{V}{gene symbol (human)}
  \item{E}{directed relationship for pathway}
}
}
\source{
PathwayCommons \url{https://reactome.org/content/detail/R-HSA-5673001}
}
\usage{
RAF_MAP_graph
}
\description{
Reactome pathway R-HSA-5673001 for the interactions in the RAF/MAP kinase cascade
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tfgb_smad_graph.R
\docType{data}
\name{TGFBeta_Smad_graph}
\alias{TGFBeta_Smad_graph}
\title{TGF-\eqn{\beta} receptor signaling activates SMADs}
\format{
A graph object of 32 vertices and 173 edges:
\describe{
  \item{V}{gene symbol (human)}
  \item{E}{directed relationship for pathway}
  \item{state}{type of relationship (activating or inhibiting) as edge attribute}
}
}
\source{
PathwayCommons \url{https://reactome.org/content/detail/R-HSA-2173789}
}
\usage{
TGFBeta_Smad_graph
}
\description{
Reactome pathway R-HSA-2173789 for the interactions in the TGF-\eqn{\beta} receptor signaling activates SMADs
}
\keyword{datasets}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/state.R
\name{make_state}
\alias{make_state}
\alias{make_state_matrix}
\title{Make State Matrix}
\usage{
make_state_matrix(graph, state = NULL)
}
\arguments{
\item{graph}{An \code{\link[igraph:aaa-igraph-package]{igraph}} object. May be directed or weighted as long as
a shortest path can be computed.}

\item{state}{numeric vector. Vector of length E(graph). Sign used to calculate state matrix,
may be an integer state or inferred directly from expected correlations for each edge. May be
applied a scalar across all edges or as a vector for each edge respectively. May also be
entered as text for "activating" or "inhibiting" or as integers for activating (0,1) or
inhibiting (-1,2). Compatible with inputs for \code{\link[graphsim]{plot_directed}}. Vector 
input is supported either directly calling the function with a value for each edge in 
 \code{E(graph)} or as an edge "attribute" in the igraph object (using 
 \code{E(g)$state <- states}).}
}
\value{
An integer matrix indicating the resolved state
(activating or inhibiting for each edge or path between nodes)
}
\description{
Functions to compute the matrix of states (1 for activating and -1 for inhibiting) 
for link signed correlations, from a vector of edge states to a signed adjacency matrix for use
in \code{\link[graphsim]{generate_expression}}. This resolves edge states to determine the sign
of all correlations between nodes in a network. These are computed interally for sigma matrices
as required.
}
\examples{

# construct a synthetic graph module
library("igraph")
graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)

 # compute state matrix for toy example
state_matrix <- make_state_matrix(graph_test)

# construct a synthetic graph network
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)

# compute state matrix for toy network
graph_structure_state_matrix <- make_state_matrix(graph_structure)
graph_structure_state_matrix

# compute state matrix for toy network with inhibitions
edge_state <- c(1, 1, -1, 1, 1, 1, 1, -1)
# edge states are a variable
graph_structure_state_matrix <- make_state_matrix(graph_structure, state = edge_state)
graph_structure_state_matrix

# compute state matrix for toy network with inhibitions
E(graph_structure)$state <- c(1, 1, -1, 1, 1, 1, 1, -1)
# edge states are a graph attribute
graph_structure_state_matrix <- make_state_matrix(graph_structure)
graph_structure_state_matrix

library("igraph")
graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)
state_matrix <- make_state_matrix(graph_test)

# import graph from package for reactome pathway
# TGF-\eqn{\Beta} receptor signaling activates SMADs (R-HSA-2173789)
TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)

# compute sigma (\eqn{\Sigma}) matrix from geometric distance directly from TGF-\eqn{\Beta} pathway
TFGBeta_Smad_state <- E(TGFBeta_Smad_graph)$state
table(TFGBeta_Smad_state)
# states are edge attributes
state_matrix_TFGBeta_Smad <- make_state_matrix(TGFBeta_Smad_graph)
# visualise matrix
library("gplots")
heatmap.2(state_matrix_TFGBeta_Smad , scale = "none", trace = "none",
          dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = colorpanel(50, "blue", "white", "red"))

# compare the states to the sign of expected correlations in the sigma matrix
sigma_matrix_TFGBeta_Smad_inhib <- make_sigma_mat_dist_graph(TGFBeta_Smad_graph,
                                                             cor = 0.8,
                                                             absolute = FALSE)
# visualise matrix
heatmap.2(sigma_matrix_TFGBeta_Smad_inhib,
          scale = "none", trace = "none",
          dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = colorpanel(50, "blue", "white", "red"))

# compare the states to the sign of final correlations in the simulated matrix
TFGBeta_Smad_data <- generate_expression(100, TGFBeta_Smad_graph, cor = 0.8)
heatmap.2(cor(t(TFGBeta_Smad_data)), scale = "none", trace = "none",
          dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = colorpanel(50, "blue", "white", "red"))


}
\seealso{
See also \code{\link[graphsim]{generate_expression}} for computing the simulated data,
\code{\link[graphsim]{make_sigma}} for computing the Sigma (\eqn{\Sigma}) matrix,
and
\code{\link[graphsim]{make_distance}} for computing distance from a graph object.

See also \code{\link[graphsim]{plot_directed}} for plotting graphs or 
\code{\link[gplots]{heatmap.2}} for plotting matrices.

See also \code{\link[graphsim]{make_laplacian}}, \code{\link[graphsim]{make_commonlink}}, 
or \code{\link[graphsim]{make_adjmatrix}} for computing input matrices.

See also \code{\link[igraph:aaa-igraph-package]{igraph}} for handling graph objects.

Other graphsim functions: 
\code{\link{generate_expression}()},
\code{\link{make_adjmatrix}},
\code{\link{make_commonlink}},
\code{\link{make_distance}},
\code{\link{make_laplacian}},
\code{\link{make_sigma}},
\code{\link{plot_directed}()}

Other generate simulated expression functions: 
\code{\link{generate_expression}()},
\code{\link{make_distance}},
\code{\link{make_sigma}}
}
\author{
Tom Kelly \email{tom.kelly@riken.jp}
}
\concept{generate simulated expression functions}
\concept{graphsim functions}
\keyword{graph}
\keyword{igraph}
\keyword{mvtnorm}
\keyword{network}
\keyword{simulation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/generate.R
\name{generate_expression}
\alias{generate_expression}
\alias{generate_expression_mat}
\title{Generate Simulated Expression}
\usage{
generate_expression(
  n,
  graph,
  state = NULL,
  cor = 0.8,
  mean = 0,
  sd = 1,
  comm = FALSE,
  dist = FALSE,
  absolute = FALSE,
  laplacian = FALSE
)

generate_expression_mat(
  n,
  mat,
  state = NULL,
  cor = 0.8,
  mean = 0,
  sd = 1,
  comm = FALSE,
  dist = FALSE,
  absolute = FALSE,
  laplacian = FALSE
)
}
\arguments{
\item{n}{number of observations (simulated samples).}

\item{graph}{An \code{\link[igraph:aaa-igraph-package]{igraph}} object. May must be 
directed if states are used.}

\item{state}{numeric vector. Vector of length E(graph). Sign used
to calculate state matrix, may be an integer state or inferred directly
from expected correlations for each edge. May be applied a scalar across
all edges or as a vector for each edge respectively. May also be entered
as text for "activating" or "inhibiting" or as integers for activating (0,1)
or inhibiting (-1,2). Compatible with inputs for \code{\link[graphsim]{plot_directed}}.
Also takes a pre-computed state matrix from \code{\link[graphsim]{make_state}}
if applied to the same graph multiple times.}

\item{cor}{numeric. Simulated maximum correlation/covariance of two 
adjacent nodes. Default to 0.8.}

\item{mean}{mean value of each simulated gene. Defaults to 0.
May be entered as a scalar applying to 
all genes or a vector with a separate value for each.}

\item{sd}{standard deviations of each gene. Defaults to 1.
May be entered as a scalar applying to 
all genes or a vector with a separate value for each.}

\item{comm, absolute, laplacian}{logical. Parameters for Sigma matrix
generation. Passed on to \code{\link[graphsim]{make_sigma}} 
or \code{\link[graphsim]{make_sigma}}.}

\item{dist}{logical. Whether a graph distance 
\code{\link[graphsim:make_sigma]{make_sigma_mat_graph}} or derived matrix
\code{\link[graphsim:make_sigma]{make_sigma_mat_dist_graph}} is used to compute the
sigma matrix (using \code{\link[graphsim]{make_distance}}).}

\item{mat}{precomputed adjacency, laplacian, commonlink, or scaled 
distance matrix (generated by \code{\link[graphsim]{make_distance}}).}
}
\value{
numeric matrix of simulated data (log-normalised counts)
}
\description{
Compute simulated continuous expression data from a graph 
network structure. Requires an \code{\link[igraph:aaa-igraph-package]{igraph}} pathway 
structure and a matrix of states (1 for activating and -1 for 
inhibiting) for link signed correlations, from a vector of edge states 
to a signed adjacency matrix for use in 
\code{\link[graphsim]{generate_expression}}. 
Uses graph structure to pass a sigma covariance matrix from 
\code{\link[graphsim:make_sigma]{make_sigma_mat_graph}} or 
\code{\link[graphsim:make_sigma]{make_sigma_mat_dist_graph}} on to 
\code{\link[mvtnorm:Mvnorm]{rmvnorm}}. By default data is generated with a mean of
 0 and standard deviation of 1 for each gene (with correlations between 
 derived from the graph structure).
}
\examples{

# construct a synthetic graph module
library("igraph")
graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)

# compute a simulated dataset for toy example
# n = 100 samples
# cor = 0.8 max correlation between samples
# absolute = FALSE (geometric distance by default)
test_data <- generate_expression(100, graph_test, cor = 0.8)
##' # visualise matrix
library("gplots")
# expression data
heatmap.2(test_data, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))
# correlations
heatmap.2(cor(t(test_data)), scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))
# expected correlations (\eqn{\Sigma})
sigma_matrix <- make_sigma_mat_graph(graph_test, cor = 0.8)
heatmap.2(make_sigma_mat_graph(graph_test, cor = 0.8),
          scale = "none", trace = "none", 
          col = colorpanel(50, "white", "red"))

# compute adjacency matrix for toy example
adjacency_matrix <- make_adjmatrix_graph(graph_test)
# generate simulated data from adjacency matrix input
test_data <- generate_expression_mat(100, adjacency_matrix, cor = 0.8)

# compute a simulated dataset for toy example
# n = 100 samples
# cor = 0.8 max correlation between samples
# absolute = TRUE (arithmetic distance)
test_data <- generate_expression(100, graph_test, cor = 0.8, absolute = TRUE)
##' # visualise matrix
library("gplots")
# expression data
heatmap.2(test_data, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))
# correlations
heatmap.2(cor(t(test_data)),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))
# expected correlations (\eqn{\Sigma})
sigma_matrix <- make_sigma_mat_graph(graph_test, cor = 0.8)
heatmap.2(make_sigma_mat_graph(graph_test, cor = 0.8),
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))

# construct a synthetic graph network
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)

# compute a simulated dataset for toy network
# n = 250 samples
# state = edge_state (properties of each edge)
# cor = 0.95 max correlation between samples
# absolute = FALSE (geometric distance by default)
edge_state <- c(1, 1, -1, 1, 1, 1, 1, -1)
structure_data <- generate_expression(250, graph_structure,
                                      state = edge_state, cor = 0.95)
##' # visualise matrix
library("gplots")
# expression data
heatmap.2(structure_data, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))
# correlations
heatmap.2(cor(t(structure_data)), scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))
# expected correlations (\eqn{\Sigma})
sigma_matrix <- make_sigma_mat_graph(graph_structure,
                                     state = edge_state, cor = 0.8)
heatmap.2(make_sigma_mat_graph(graph_structure,
                               state = edge_state, cor = 0.8),
          scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))

# compute adjacency matrix for toy network
graph_structure_adjacency_matrix <- make_adjmatrix_graph(graph_structure)
# define states for for each edge
edge_state <- c(1, 1, -1, 1, 1, 1, 1, -1)
# generate simulated data from adjacency matrix input
structure_data <- generate_expression_mat(250, graph_structure_adjacency_matrix,
                                          state = edge_state, cor = 0.8)

# compute a simulated dataset for toy network
# n = 1000 samples
# state = TGFBeta_Smad_state (properties of each edge)
# cor = 0.75 max correlation between samples
# absolute = FALSE (geometric distance by default)
 # compute states directly from graph attributes for TGF-\eqn{\Beta} pathway
TGFBeta_Smad_state <- E(TGFBeta_Smad_graph)$state
table(TGFBeta_Smad_state)
# generate simulated data
TGFBeta_Smad_data <- generate_expression(1000, TGFBeta_Smad_graph, cor = 0.75)
##' # visualise matrix
library("gplots")
# expression data
heatmap.2(TGFBeta_Smad_data, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))
# correlations
heatmap.2(cor(t(TGFBeta_Smad_data)), scale = "none", trace = "none",
          dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = colorpanel(50, "blue", "white", "red"))
# expected correlations (\eqn{\Sigma})
sigma_matrix <- make_sigma_mat_dist_graph(TGFBeta_Smad_graph, cor = 0.75)
heatmap.2(make_sigma_mat_dist_graph(TGFBeta_Smad_graph, cor = 0.75),
          scale = "none", trace = "none",
          dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = colorpanel(50, "blue", "white", "red"))


# generate simulated data (absolute distance and shared edges)
TGFBeta_Smad_data <- generate_expression(1000, TGFBeta_Smad_graph,
                                         cor = 0.75, absolute = TRUE, comm = TRUE)
##' # visualise matrix
library("gplots")
# expression data
heatmap.2(TGFBeta_Smad_data, scale = "none", trace = "none",
          col = colorpanel(50, "blue", "white", "red"))
# correlations
heatmap.2(cor(t(TGFBeta_Smad_data)), scale = "none", trace = "none",
          dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = colorpanel(50, "blue", "white", "red"))
# expected correlations (\eqn{\Sigma})
sigma_matrix <- make_sigma_mat_graph(TGFBeta_Smad_graph,
                                     cor = 0.75, comm = TRUE)
heatmap.2(make_sigma_mat_graph(TGFBeta_Smad_graph, cor = 0.75, comm = TRUE),
          scale = "none", trace = "none",
          dendrogram = "none", Rowv = FALSE, Colv = FALSE,
          col = colorpanel(50, "blue", "white", "red"))

}
\seealso{
See also \code{\link[graphsim]{make_sigma}} for computing the Sigma (\eqn{\Sigma}) matrix,
\code{\link[graphsim]{make_distance}} for computing distance from a graph object,
and
\code{\link[graphsim]{make_state}} for resolving inhibiting states.

See also \code{\link[graphsim]{plot_directed}} for plotting graphs or 
\code{\link[gplots]{heatmap.2}} for plotting matrices.

See also \code{\link[graphsim]{make_laplacian}}, \code{\link[graphsim]{make_commonlink}}, 
or \code{\link[graphsim]{make_adjmatrix}} for computing input matrices.

See also \code{\link[igraph:aaa-igraph-package]{igraph}} for handling graph objects.

Other graphsim functions: 
\code{\link{make_adjmatrix}},
\code{\link{make_commonlink}},
\code{\link{make_distance}},
\code{\link{make_laplacian}},
\code{\link{make_sigma}},
\code{\link{make_state}},
\code{\link{plot_directed}()}

Other generate simulated expression functions: 
\code{\link{make_distance}},
\code{\link{make_sigma}},
\code{\link{make_state}}
}
\author{
Tom Kelly \email{tom.kelly@riken.jp}
}
\concept{generate simulated expression functions}
\concept{graphsim functions}
\keyword{graph}
\keyword{igraph}
\keyword{mvtnorm}
\keyword{network}
\keyword{simulation}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/distance.R
\name{make_distance}
\alias{make_distance}
\alias{make_distance_graph}
\alias{make_relationship}
\alias{make_distance_adjmat}
\alias{make_distance_comm}
\alias{make_distance_laplacian}
\title{Generate Distance Matrix}
\usage{
make_distance_graph(graph, directed = FALSE, absolute = FALSE)

make_distance_adjmat(mat, directed = FALSE, absolute = FALSE)

make_distance_comm(mat, directed = FALSE, absolute = FALSE)

make_distance_laplacian(mat, directed = FALSE, absolute = FALSE)
}
\arguments{
\item{graph}{An \code{\link[igraph:aaa-igraph-package]{igraph}} object. May be directed or weighted.}

\item{directed}{logical. Whether directed information is passed to the distance matrix.}

\item{absolute}{logical. Whether distances are scaled as the absolute difference
from the diameter (maximum possible). Defaults to TRUE. The alternative is to
calculate a relative difference from the diameter for a geometric decay in distance.}

\item{mat}{precomputed adjacency or commonlink matrix.}
}
\value{
A numeric matrix of values in the range [0, 1] where higher values are closer in the network
}
\description{
Compute the distance matrix of using shortest paths of a (directed)
\code{\link[igraph:aaa-igraph-package]{igraph}} structure, normalising by the diameter of the network,
preserving node/column/row names (and direction). This is used to compute the
simulatted data for \code{\link[graphsim]{generate_expression}} (when \code{dist = TRUE})
by \code{\link[graphsim:make_sigma]{make_sigma_mat_dist_graph}}.
}
\examples{

# construct a synthetic graph module
library("igraph")
graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)

# compute adjacency matrix for toy example
adjacency_matrix <- make_adjmatrix_graph(graph_test)
# compute nodes with relationships between nodes (geometrically decreasing by default)
distance_matrix_geom <- make_distance_adjmat(adjacency_matrix)
distance_matrix_geom

# compute nodes with relationships between nodes (arithmetically decreasing)
distance_matrix_abs <- make_distance_adjmat(adjacency_matrix, absolute = TRUE)
distance_matrix_abs

# compute Laplacian matrix
laplacian_matrix <- make_laplacian_graph(graph_test)
# compute distances from Laplacian
distance_matrix <- make_distance_laplacian(laplacian_matrix)

# construct a synthetic graph network
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
# compute adjacency matrix for toy network
graph_structure_adjacency_matrix <- make_adjmatrix_graph(graph_structure)
# compute nodes with relationships between nodes (geometrically decreasing by default)
graph_structure_distance_matrix_geom <- make_distance_adjmat(graph_structure_adjacency_matrix)
graph_structure_distance_matrix_geom
# visualise matrix
library("gplots")
heatmap.2(graph_structure_distance_matrix_geom, scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))
# compute nodes with relationships between nodes (arithmetically decreasing)
graph_structure_distance_matrix_abs <- make_distance_adjmat(graph_structure_adjacency_matrix,
                                                            absolute = TRUE)
graph_structure_distance_matrix_abs
# visualise matrix
library("gplots")
heatmap.2(graph_structure_distance_matrix_abs,
          scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))
          
# import graph from package for reactome pathway
# TGF-\eqn{\Beta} receptor signaling activates SMADs (R-HSA-2173789)
TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)
# compute nodes with relationships between nodes (geometrically decreasing by default)
TGFBeta_Smad_adjacency_matrix <- make_adjmatrix_graph(TGFBeta_Smad_graph)
TGFBeta_Smad_distance_matrix_geom <- make_distance_adjmat(TGFBeta_Smad_adjacency_matrix)
# visualise matrix
library("gplots")
heatmap.2(TGFBeta_Smad_distance_matrix_geom, scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))
# compute nodes with relationships between nodes (arithmetically decreasing)
TGFBeta_Smad_distance_matrix_abs <- make_distance_adjmat(TGFBeta_Smad_adjacency_matrix,
                        absolute = TRUE)
# visualise matrix
library("gplots")
heatmap.2(TGFBeta_Smad_distance_matrix_abs, scale = "none", trace = "none",
          col = colorpanel(50, "white", "red"))

}
\seealso{
See also \code{\link[graphsim]{generate_expression}} for computing the simulated data,
\code{\link[graphsim]{make_sigma}} for computing the Sigma (\eqn{\Sigma}) matrix,
\code{\link[graphsim]{make_state}} for resolving inhibiting states.

See also \code{\link[graphsim]{plot_directed}} for plotting graphs or 
\code{\link[gplots]{heatmap.2}} for plotting matrices.

See also \code{\link[graphsim]{make_laplacian}}, \code{\link[graphsim]{make_commonlink}}, 
or \code{\link[graphsim]{make_adjmatrix}} for computing input matrices.

See also \code{\link[igraph:aaa-igraph-package]{igraph}} for handling graph objects.

Other graphsim functions: 
\code{\link{generate_expression}()},
\code{\link{make_adjmatrix}},
\code{\link{make_commonlink}},
\code{\link{make_laplacian}},
\code{\link{make_sigma}},
\code{\link{make_state}},
\code{\link{plot_directed}()}

Other generate simulated expression functions: 
\code{\link{generate_expression}()},
\code{\link{make_sigma}},
\code{\link{make_state}}
}
\author{
Tom Kelly \email{tom.kelly@riken.jp}
}
\concept{generate simulated expression functions}
\concept{graphsim functions}
\keyword{adjacency}
\keyword{graph}
\keyword{igraph}
\keyword{network}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/adjmat.R
\name{make_adjmatrix}
\alias{make_adjmatrix}
\alias{make_adjmatrix_graph}
\title{Generate Adjacency Matrix}
\usage{
make_adjmatrix_graph(graph, directed = FALSE)
}
\arguments{
\item{graph}{An \code{\link[igraph:aaa-igraph-package]{igraph}} object. May be directed or weighted.}

\item{directed}{logical. Whether directed information is passed to the adjacency matrix.}
}
\value{
An adjacency matrix compatible with generating an expression matrix
}
\description{
Compute the adjacency matrix of a (directed) \code{\link[igraph:aaa-igraph-package]{igraph}}
structure, preserving node/column/row names (and direction).
}
\examples{

# construct a synthetic graph module
library("igraph")
graph_test_edges <- rbind(c("A", "B"), c("B", "C"), c("B", "D"))
graph_test <- graph.edgelist(graph_test_edges, directed = TRUE)

# compute adjacency matrix for toy example
adjacency_matrix <- make_adjmatrix_graph(graph_test)
adjacency_matrix

# construct a synthetic graph network
graph_structure_edges <- rbind(c("A", "C"), c("B", "C"), c("C", "D"), c("D", "E"),
                               c("D", "F"), c("F", "G"), c("F", "I"), c("H", "I"))
graph_structure <- graph.edgelist(graph_structure_edges, directed = TRUE)
# compute adjacency matrix for toy network
graph_structure_adjacency_matrix <- make_adjmatrix_graph(graph_structure)
graph_structure_adjacency_matrix

# import graph from package for reactome pathway
# TGF-\eqn{\Beta} receptor signaling activates SMADs (R-HSA-2173789)
TGFBeta_Smad_graph <- identity(TGFBeta_Smad_graph)

# compute adjacency matrix for TGF-\eqn{\Beta} receptor signaling activates SMADs
TGFBeta_Smad_adjacency_matrix <- make_adjmatrix_graph(TGFBeta_Smad_graph)
dim(TGFBeta_Smad_adjacency_matrix)
TGFBeta_Smad_adjacency_matrix[1:12, 1:12]

}
\seealso{
See also \code{\link[graphsim]{generate_expression}} for computing the simulated data,
\code{\link[graphsim]{make_sigma}} for computing the Sigma (\eqn{\Sigma}) matrix,
\code{\link[graphsim]{make_distance}} for computing distance from a graph object,
\code{\link[graphsim]{make_state}} for resolving inhibiting states.

See also \code{\link[graphsim]{plot_directed}} for plotting graphs or 
\code{\link[gplots]{heatmap.2}} for plotting matrices.

See also \code{\link[graphsim]{make_laplacian}}
or  \code{\link[graphsim]{make_commonlink}} for computing input matrices.

See also \code{\link[igraph:aaa-igraph-package]{igraph}} for handling graph objects.

Other graphsim functions: 
\code{\link{generate_expression}()},
\code{\link{make_commonlink}},
\code{\link{make_distance}},
\code{\link{make_laplacian}},
\code{\link{make_sigma}},
\code{\link{make_state}},
\code{\link{plot_directed}()}

Other graph conversion functions: 
\code{\link{make_commonlink}},
\code{\link{make_laplacian}}
}
\author{
Tom Kelly \email{tom.kelly@riken.jp}
}
\concept{graph conversion functions}
\concept{graphsim functions}
\keyword{adjacency}
\keyword{graph}
\keyword{igraph}
\keyword{network}
% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pi3k_akt_graph.R
\docType{data}
\name{Pi3K_AKT_graph}
\alias{Pi3K_AKT_graph}
\title{PI3K/AKT activation}
\format{
A graph object of 275 vertices and 21106 edges:
\describe{
  \item{V}{gene symbol (human)}
  \item{E}{directed relationship for pathway}
  \item{state}{type of relationship (activating or inhibiting) as edge attribute}
}
}
\source{
PathwayCommons \url{https://reactome.org/content/detail/R-HSA-198203}
}
\usage{
Pi3K_AKT_graph
}
\description{
Reactome pathway R-HSA-198203 for the interactions in the phosphoinositide-3-kinase activation of Protein kinase B (PKB), also known as Akt
}
\keyword{datasets}
